import os

from common.params import ParamRanges

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


def ranges() -> ParamRanges:
    return ParamRanges.from_file(os.path.join(FIXTURES, "parameter_ranges.json"))
