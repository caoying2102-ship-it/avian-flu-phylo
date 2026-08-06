import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "03_trimal" / "scripts" / "trimal.py"

spec = importlib.util.spec_from_file_location("trimal", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class TrimAlValidationTest(unittest.TestCase):
    def test_allows_reduced_sequence_count_when_trimal_drops_gap_only_sequences(self):
        module.validate_output_consistency(
            segment="PB2",
            input_sequence_count=4372,
            output_sequence_count=4364,
            input_alignment_length=1000,
            output_alignment_length=1000,
        )

    def test_raises_when_alignment_length_increases(self):
        with self.assertRaises(RuntimeError):
            module.validate_output_consistency(
                segment="PB2",
                input_sequence_count=4372,
                output_sequence_count=4372,
                input_alignment_length=1000,
                output_alignment_length=1001,
            )


if __name__ == "__main__":
    unittest.main()
