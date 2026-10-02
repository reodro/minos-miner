"""The no-container launcher: docker argv parsing, and one real unshare of a tiny rootfs."""
import os
import shutil
import unittest
from unittest import mock
from pathlib import Path

from evaluator.sandbox import DockerSpec, command_for, parse_docker_run, start_in


class ParseTest(unittest.TestCase):
    def test_gatk_command_keeps_mounts_limits_and_java_args(self):
        spec = parse_docker_run([
            "docker", "run", "--rm", "--name", "minos-run-round-abc", "--cpus=1", "--memory=4g",
            "--network", "none", "-v", "/bam:/data/bams:ro", "-v", "/ref:/data/reference:ro",
            "-v", "/out:/data/output", "-e", "RTG_MEM=8g", "broadinstitute/gatk:4.5.0.0",
            "java", "-Xmx3276m", "HaplotypeCaller", "-L", "chr18:1-2",
        ])
        self.assertEqual(spec.name, "minos-run-round-abc")
        self.assertEqual(spec.cpus, 1.0)
        self.assertEqual(spec.memory, "4g")
        self.assertTrue(spec.network_none)
        self.assertEqual(spec.binds, [("/bam", "/data/bams", True), ("/ref", "/data/reference", True),
                                      ("/out", "/data/output", False)])
        self.assertEqual(spec.env, {"RTG_MEM": "8g"})
        self.assertEqual(spec.image, "broadinstitute/gatk:4.5.0.0")
        self.assertEqual(spec.argv[:2], ["java", "-Xmx3276m"])

    def test_entrypoint_is_prepended_and_args_replace_cmd(self):
        spec = DockerSpec(image="hap", argv=["/data/truth.vcf.gz", "/data/query.vcf.gz"])
        cmd = command_for(spec, {"Entrypoint": ["hap.py"], "Cmd": None})
        self.assertEqual(cmd, ["hap.py", "/data/truth.vcf.gz", "/data/query.vcf.gz"])
        bare = command_for(DockerSpec(image="gatk", argv=["java", "HaplotypeCaller"]),
                           {"Entrypoint": None, "Cmd": ["bash"]})
        self.assertEqual(bare, ["java", "HaplotypeCaller"])


class ChrootPathTest(unittest.TestCase):
    def _script(self, sandbox):
        return sandbox._script(Path("/rootfs"), DockerSpec(image="img"), {}, ["echo", "hi"],
                               "/mnt", "/upper", "/work")

    def test_this_machine_uses_a_chroot_that_exists(self):
        from evaluator import sandbox
        path = sandbox._chroot_bin()
        self.assertTrue(os.access(path, os.X_OK))
        self.assertIn(path, self._script(sandbox))

    def test_newer_ubuntu_chroot_in_usr_bin_is_used(self):
        from unittest import mock
        from evaluator import sandbox

        def access(path, mode):
            return path == "/usr/bin/chroot" and mode == os.X_OK

        with mock.patch.object(sandbox.os, "access", access):
            script = self._script(sandbox)
        self.assertIn("/usr/bin/chroot ", script)
        self.assertNotIn("/usr/sbin/chroot", script)

    def test_missing_chroot_is_reported_before_the_script_runs(self):
        from unittest import mock
        from evaluator import sandbox
        with mock.patch.object(sandbox.os, "access", lambda *_a, **_k: False), \
                mock.patch.object(sandbox.shutil, "which", lambda *_a, **_k: None):
            with self.assertRaises(RuntimeError):
                sandbox._chroot_bin()


class CacheHitTest(unittest.TestCase):
    def test_a_cached_image_does_not_call_docker(self):
        import tempfile
        from unittest import mock

        from evaluator import sandbox

        root = Path(tempfile.mkdtemp())
        try:
            ident = "a" * 64
            dest = root / ident
            dest.mkdir()
            (dest / ".ready").write_text(ident + "\n")
            (dest / "config.json").write_text('{"Env": [], "Entrypoint": ["hap.py"], "Cmd": null, "WorkingDir": "/data"}')
            with mock.patch.object(sandbox, "_cache_dir", return_value=root):
                pinned = sandbox.ensure_rootfs(f"genonet/hap-py@sha256:{ident}")
                self.assertEqual(pinned[0], dest)
                sandbox._remember("broadinstitute/gatk:4.5.0.0", ident)

                def boom(*_a, **_k):
                    raise AssertionError("docker was called")

                with mock.patch.object(sandbox.subprocess, "run", boom):
                    tagged = sandbox.ensure_rootfs("broadinstitute/gatk:4.5.0.0")
                self.assertEqual(tagged[0], dest)
        finally:
            shutil.rmtree(root, ignore_errors=True)


def _echo_root(dest: Path) -> None:
    """A rootfs whose /bin/echo is the host binary plus the libraries it needs."""
    for src in ("/bin/echo", "/lib64/ld-linux-x86-64.so.2", "/lib/x86_64-linux-gnu/libc.so.6"):
        target = dest / src[1:]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src, target)


@unittest.skipUnless(os.geteuid() == 0, "unshare --mount needs root")
class RunTest(unittest.TestCase):
    def test_echo_runs_inside_the_root_and_can_be_killed(self):
        import tempfile
        from evaluator.sandbox import kill_name
        root = Path(tempfile.mkdtemp())
        try:
            _echo_root(root)
            shutil.copy("/bin/sleep", root / "bin" / "sleep")
            spec = DockerSpec(image="tiny", argv=["/bin/echo", "hello-sandbox"], name="minos-test-echo",
                              cpus=1, memory="64M", network_none=True)
            proc = start_in(root, {"WorkingDir": "/", "Env": ["PATH=/bin"], "Entrypoint": None, "Cmd": None},
                            spec, stdout=-1)
            out = proc.communicate(timeout=30)[0]
            self.assertEqual(proc.returncode, 0, out)
            self.assertIn(b"hello-sandbox", out)
            sleeper = DockerSpec(image="tiny", argv=["/bin/sleep", "30"], name="minos-test-sleep",
                                 cpus=1, memory="64M", network_none=True)
            running = start_in(root, {"WorkingDir": "/", "Env": ["PATH=/bin"], "Entrypoint": None, "Cmd": None},
                               sleeper, stdout=-1)
            self.assertIsNone(running.poll())
            kill_name("minos-test-sleep")
            running.wait(timeout=10)
            self.assertIsNotNone(running.poll())
        finally:
            shutil.rmtree(root, ignore_errors=True)


class ScoreLimitTest(unittest.TestCase):
    HAPPY = ["docker", "run", "--rm", "-v", "/a:/data/truth", "genonet/hap-py@sha256:abc", "/data/truth/t.vcf.gz"]

    def test_hap_py_gets_the_scoring_cpu_limit(self):
        from evaluator.sandbox import parse_docker_run, with_score_limit
        with mock.patch.dict(os.environ, {"MINOS_SCORE_CPUS": ""}):
            spec = parse_docker_run(with_score_limit(list(self.HAPPY)))
        self.assertEqual(spec.cpus, 1.0)
        self.assertEqual(spec.image, "genonet/hap-py@sha256:abc")

    def test_setting_and_off_switch(self):
        from evaluator.sandbox import parse_docker_run, with_score_limit
        with mock.patch.dict(os.environ, {"MINOS_SCORE_CPUS": "2"}):
            self.assertEqual(parse_docker_run(with_score_limit(list(self.HAPPY))).cpus, 2.0)
        with mock.patch.dict(os.environ, {"MINOS_SCORE_CPUS": "0"}):
            self.assertEqual(with_score_limit(list(self.HAPPY)), self.HAPPY)

    def test_other_images_and_existing_limits_untouched(self):
        from evaluator.sandbox import with_score_limit
        bcf = ["docker", "run", "--rm", "quay.io/biocontainers/bcftools:1.20", "bcftools", "view"]
        self.assertEqual(with_score_limit(list(bcf)), bcf)
        own = ["docker", "run", "--cpus=3", "genonet/hap-py@sha256:abc"]
        self.assertEqual(with_score_limit(list(own)), own)
