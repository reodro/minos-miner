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
        if (recording) {
            draws.add(new int[]{bits, v});
        }
        return v;
    }

    public void setRecording(final boolean on) { recording = on; }

    public List<int[]> takeDraws() {
        final List<int[]> out = draws;
        draws = new ArrayList<>();
        return out;
    }

    public boolean replay(final int[] draw) {
        return super.next(draw[0]) == draw[1];
    }

    public MinosRecordingWell copyState() {
        final MinosRecordingWell c = org.apache.commons.lang3.SerializationUtils.clone(this);
        c.recording = false;
        c.draws = new ArrayList<>();
        return c;
    }
}
