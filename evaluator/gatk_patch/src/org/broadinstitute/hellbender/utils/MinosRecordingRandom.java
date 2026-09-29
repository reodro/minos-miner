package org.broadinstitute.hellbender.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * minos: a java.util.Random that behaves exactly like its parent and can record the raw draws (next(bits) calls)
 * made while recording is on, so they can be replayed on other generators to check they would return the same
 * values. Every public Random method funnels through next(bits), so recording and replaying at this level
 * reproduces the generator state exactly.
 */
public final class MinosRecordingRandom extends Random {
    private static final long serialVersionUID = 1L;
    private transient boolean recording = false;
    private transient List<int[]> draws = new ArrayList<>();

    public MinosRecordingRandom(final long seed) { super(seed); }

    @Override
    protected int next(final int bits) {
        final int v = super.next(bits);
        if (recording) {
            draws.add(new int[]{bits, v});
        }
        return v;
    }

    public void setRecording(final boolean on) { recording = on; }

    /** Draws recorded since the last call, then cleared. */
    public List<int[]> takeDraws() {
        final List<int[]> out = draws;
        draws = new ArrayList<>();
        return out;
    }

    /** Replay one recorded draw; true if this generator returns the same value. */
    public boolean replay(final int[] draw) {
        return super.next(draw[0]) == draw[1];
    }

    /** A copy with identical state (seed position and cached Gaussian), not recording. */
    public MinosRecordingRandom copyState() {
        final MinosRecordingRandom c = org.apache.commons.lang3.SerializationUtils.clone(this);
        c.recording = false;
        c.draws = new ArrayList<>();
        return c;
    }
}
