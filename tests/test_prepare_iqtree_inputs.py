import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "04_iqtree_initial" / "scripts" / "prepare_iqtree_inputs.py"
IQTREE_MODULE_PATH = Path(__file__).resolve().parents[1] / "04_iqtree_initial" / "scripts" / "run_iqtree.py"

spec = importlib.util.spec_from_file_location("prepare_iqtree_inputs", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

run_iqtree_spec = importlib.util.spec_from_file_location("run_iqtree", IQTREE_MODULE_PATH)
run_iqtree = importlib.util.module_from_spec(run_iqtree_spec)
run_iqtree_spec.loader.exec_module(run_iqtree)


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

    def test_copy_metadata_to_output_copies_paired_metadata_workbook(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            fasta_dir = tmpdir / "input"
            outdir = tmpdir / "output"
            fasta_dir.mkdir()

            fasta_path = fasta_dir / "HA_reference_dedup_aligned_trimmed.fasta"
            metadata_path = fasta_dir / "HA_metadata_reference_dedup_aligned_trimmed.xlsx"
            fasta_path.write_text(">seq1\nACGT\n", encoding="utf-8")
            metadata_path.write_text("metadata", encoding="utf-8")

            copied = run_iqtree.copy_metadata_to_output(fasta_path, outdir)

            self.assertEqual(copied, outdir / metadata_path.name)
            self.assertTrue((outdir / metadata_path.name).exists())
            self.assertEqual(
                (outdir / metadata_path.name).read_text(encoding="utf-8"),
                "metadata",
            )


if __name__ == "__main__":
    unittest.main()
