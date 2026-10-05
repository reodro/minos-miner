package org.broadinstitute.hellbender;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.hellbender.utils.Utils;

/**
 * minos: one JVM that runs GATK jobs one after another, so each job after the first starts with loaded classes and
 * JIT-compiled code (a cold run spends its first part in slow, not-yet-compiled code: measured 162 s cold vs 119-127 s
 * warm on a full round, identical VCF records).
 *
 *   java ... org.broadinstitute.hellbender.MinosServer <dir> <max jobs> [<idle seconds>, default 600]
 *
 * The runner writes <dir>/job-00001.tsv, job-00002.tsv, ... (atomically): the first line holds the job's JVM
 * properties as tab-separated key=value (may be empty), every further line one GATK argument. The server answers
 * <dir>/job-NNNNN.done: "ok <seconds>" or "error <message>". It writes <dir>/ready once at start, and leaves (writing
 * <dir>/exit with the reason) on <dir>/stop, after <max jobs> jobs, when the heap stays too full after a GC, or after
 * <idle seconds> without a job (a server whose runner is gone - e.g. a finished review's executor - ends by itself).
 *
 * Exactness: before and after every job every static the patch keeps goes back to a new JVM's state
 * (Utils.minosResetAll: both generators at the seed with nothing recorded and no cached gaussian, no overrides, no
 * activity hook, no haplotype-cap report); job properties are removed again. GATK itself is entered exactly as from
 * the command line (Main.instanceMain), which builds every tool, engine and reader anew.
 */
public final class MinosServer {
    private static final double MAX_HEAP_USED_AFTER_GC = 0.6;

    private MinosServer() {}

    public static void main(final String[] argv) throws Exception {
        final Path dir = Paths.get(argv[0]);
        final int maxJobs = Integer.parseInt(argv[1]);
        final long idleNanos = (argv.length > 2 ? Long.parseLong(argv[2]) : 600L) * 1_000_000_000L;
        write(dir.resolve("ready"), "ok\n");
        String reason = "max jobs";
        for (int n = 1; n <= maxJobs; n++) {
            final Path job = dir.resolve(String.format("job-%05d.tsv", n));
            final long waitingSince = System.nanoTime();
            while (!Files.exists(job)) {
                if (Files.exists(dir.resolve("stop"))) {
                    write(dir.resolve("exit"), "stop\n");
                    return;
                }
                if (System.nanoTime() - waitingSince > idleNanos) {
                    write(dir.resolve("exit"), "idle\n");
                    return;
                }
                Thread.sleep(20);
            }
            final List<String> lines = Files.readAllLines(job, StandardCharsets.UTF_8);
            final List<String> keys = new ArrayList<>();
            String answer;
            final long t0 = System.nanoTime();
            try {
                Utils.minosResetAll();
                for (final String kv : lines.get(0).split("\t")) {
                    final int eq = kv.indexOf('=');
                    if (eq > 0) {
                        keys.add(kv.substring(0, eq));
                        System.setProperty(kv.substring(0, eq), kv.substring(eq + 1));
                    }
                }
                final String[] args = lines.subList(1, lines.size()).toArray(new String[0]);
                new Main().instanceMain(args);
                answer = String.format("ok %.3f%n", (System.nanoTime() - t0) / 1e9);
            } catch (final Throwable e) {
                final StringWriter sw = new StringWriter();
                e.printStackTrace(new PrintWriter(sw));
                System.err.println(sw);
                answer = "error " + String.valueOf(e).replace('\n', ' ') + "\n";
            } finally {
                for (final String k : keys) {
                    System.clearProperty(k);
                }
                Utils.minosResetAll();
            }
            write(dir.resolve(String.format("job-%05d.done", n)), answer);
            final Runtime rt = Runtime.getRuntime();
            if (rt.totalMemory() - rt.freeMemory() > MAX_HEAP_USED_AFTER_GC * rt.maxMemory()) {
                System.gc();
                if (rt.totalMemory() - rt.freeMemory() > MAX_HEAP_USED_AFTER_GC * rt.maxMemory()) {
                    reason = "heap";
                    break;
                }
            }
        }
        write(dir.resolve("exit"), reason + "\n");
    }

    private static void write(final Path p, final String text) throws IOException {
        final Path tmp = p.resolveSibling(p.getFileName() + ".tmp");
        Files.write(tmp, text.getBytes(StandardCharsets.UTF_8));
        Files.move(tmp, p, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
    }
}
