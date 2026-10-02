"""Model-based answer-set builder (k2/107): candidates from the round BAM, public evidence per candidate, and
gradient-boosted models trained on released rounds. Entry point: K2Provider. Data: see answer_set/k2/util.py."""


def __getattr__(name):
    # lazy, so `python -m answer_set.k2.<module>` does not import the whole builder first
    if name == "K2Provider":
        from answer_set.k2.provider import K2Provider
        return K2Provider
    raise AttributeError(name)
