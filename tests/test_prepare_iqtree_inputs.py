import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "04_iqtree_initial" / "scripts" / "prepare_iqtree_inputs.py"

spec = importlib.util.spec_from_file_location("prepare_iqtree_inputs", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PrepareIQTreeInputsTest(unittest.TestCase):
    def test_resolve_segments_supports_all_segments(self):
        self.assertEqual(module.resolve_segments("ALL"), module.SEGMENTS)
        self.assertEqual(module.resolve_segments("mp"), ("MP",))

    def test_deduplicate_fasta_headers_keeps_first_occurrence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.fasta"
            output_path = Path(tmpdir) / "output.fasta"
            input_path.write_text(
                ">seqA\nACGT\n>seqB\nGGGG\n>seqA\nTTTT\n",
                encoding="utf-8",
            )

            module.deduplicate_fasta(input_path, output_path)

            content = output_path.read_text(encoding="utf-8")
            self.assertEqual(content.count(">seqA"), 1)
            self.assertIn(">seqB\nGGGG\n", content)
            self.assertNotIn(">seqA\nTTTT\n", content)


if __name__ == "__main__":
    unittest.main()
