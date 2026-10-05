package org.broadinstitute.hellbender.utils;

import java.util.HashMap;
import java.util.Map;

/**
 * minos: what --max-num-haplotypes-in-population did in this run, per assembly (stream and assembly settings).
 *
 * Every k-best haplotype search (GraphBasedKBestHaplotypeFinder) reports the smallest cap that reproduces it and
 * whether the cap stopped it. Per key: the largest such cap over the run, and whether any search was stopped. A run
 * no search of which was stopped makes identical calls under every cap of at least that value (the evaluator reuses
 * its result for those configs). Searches made while no key is set mark the whole run unknown. Counting only:
 * nothing here changes what GATK computes. HaplotypeCaller runs single-threaded, as the runner starts it.
 */
public final class MinosKBest {
    private static String key = null;
    private static boolean unattributed = false;
    private static final Map<String, int[]> STATS = new HashMap<>();   // key -> {cap needed, stopped (0/1)}

    private MinosKBest() {}

    /** minos server: forget everything before the next job in this JVM. */
    public static void reset() {
        key = null;
        unattributed = false;
        STATS.clear();
    }

    /** The assembly whose searches follow ("<stream>|<assembly settings>"); null after it. */
    public static void setKey(final String k) {
        key = k;
    }

    public static void record(final int needed, final boolean stopped) {
        if (key == null) {
            unattributed = true;
            return;
        }
        final int[] s = STATS.computeIfAbsent(key, k -> new int[2]);
        s[0] = Math.max(s[0], needed);
        if (stopped) {
            s[1] = 1;
        }
    }

    /** "needed <n>" (a cap of at least n reproduces this assembly's calls), "stopped" or "unknown". */
    public static String report(final String k) {
        if (unattributed) {
            return "unknown\n";
        }
        final int[] s = STATS.get(k);
        if (s == null) {
            // no search in the whole run: the cap was never used. Searches for other keys only: never guess.
            return STATS.isEmpty() ? "needed 0\n" : "unknown\n";
        }
        return s[1] != 0 ? "stopped\n" : "needed " + s[0] + "\n";
    }
}
