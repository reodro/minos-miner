"""Config -> callset -> score. Public API: Evaluator, GatkEvaluator, FakeEvaluator."""
from evaluator.base import Evaluator  # noqa: F401
from evaluator.fake import FakeEvaluator  # noqa: F401
from evaluator.gatk import GatkEvaluator  # noqa: F401
