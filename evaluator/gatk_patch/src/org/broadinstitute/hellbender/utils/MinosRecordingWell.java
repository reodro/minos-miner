package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.random.Well19937c;

import java.util.ArrayList;
import java.util.List;

/**
 * minos: the Well19937c behind Utils.getRandomDataGenerator(), with the same record/replay ability as
 * MinosRecordingRandom. Every draw goes through next(bits).
 */
public final class MinosRecordingWell extends Well19937c {
    private static final long serialVersionUID = 1L;
    private transient boolean recording = false;
    private transient List<int[]> draws = new ArrayList<>();

    public MinosRecordingWell(final long seed) { super(seed); }

    @Override
    protected int next(final int bits) {
        final int v = super.next(bits);
        steps++;
        if (recording) {
            draws.add(new int[]{bits, v});
        }
        return v;
    }

    public void setRecording(final boolean on) { recording = on; }

    /** minos server: back to the state of a new JVM's generator (seed position, cached gaussian, nothing recorded). */
    public void minosReset(final long seed) {
        setSeed(seed);                       // AbstractWell.setSeed clears BitsStreamGenerator's cached gaussian
        recording = false;
        draws = new ArrayList<>();
        steps = 0;
    }

    public List<int[]> takeDraws() {
        final List<int[]> out = draws;
        draws = new ArrayList<>();
        return out;
    }

    public boolean replay(final int[] draw) {
        steps++;
        return super.next(draw[0]) == draw[1];
    }

    // minos repair: every next() is one state step of the generator, so a step count is enough to skip a region's draws
    private transient long steps = 0;

    public long steps() { return steps; }

    public void skipSteps(final long n) {
        for (long i = 0; i < n; i++) {
            next(32);
        }
    }

    public MinosRecordingWell copyState() {
        final MinosRecordingWell c = org.apache.commons.lang3.SerializationUtils.clone(this);
        c.recording = false;
        c.draws = new ArrayList<>();
        return c;
    }
}
