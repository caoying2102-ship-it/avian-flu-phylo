import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "01_sample_preparation" / "scripts" / "sample.py"

spec = importlib.util.spec_from_file_location("sample", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SampleRetentionTest(unittest.TestCase):
    def test_main_keeps_required_isolate_for_non_ha_sampling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            script_path = tmpdir_path / "01_sample_preparation" / "scripts" / "sample.py"
            script_path.parent.mkdir(parents=True, exist_ok=True)
            module.__file__ = str(script_path)

            segment = "NA"
            segment_dir = (
                tmpdir_path
                / "01_sample_preparation"
                / "output"
                / "genes"
                / segment
            )
            segment_dir.mkdir(parents=True, exist_ok=True)

            metadata = pd.DataFrame(
                [
                    {
                        "Collection_Date": "2024-01-01",
                        "Subcontinent": "North America",
                        "Host": "Human",
                        "Isolate_Id": "EPI_ISL_1254",
                    },
                    {
                        "Collection_Date": "2024-01-02",
                        "Subcontinent": "North America",
                        "Host": "Human",
                        "Isolate_Id": "EPI_ISL_9999",
                    },
                    {
                        "Collection_Date": "2024-01-03",
                        "Subcontinent": "North America",
                        "Host": "Human",
                        "Isolate_Id": "EPI_ISL_10000",
                    },
                    {
                        "Collection_Date": "2024-01-04",
                        "Subcontinent": "North America",
                        "Host": "Human",
                        "Isolate_Id": "EPI_ISL_10001",
                    },
                ]
            )
            summary = pd.DataFrame(
                [
                    {
                        "Collection_Year": "2024",
                        "Subcontinent": "North America",
                        "Host": "Human",
                        "Count": 4,
                    }
                ]
            )

            metadata_path = segment_dir / f"{segment}_metadata_simplified.xlsx"
            summary_path = segment_dir / f"{segment}_summary.xlsx"
            metadata.to_excel(metadata_path, index=False)
            summary.to_excel(summary_path, index=False)

            original_target = module.TARGET_SAMPLE_SIZE
            module.TARGET_SAMPLE_SIZE = 1
            try:
                with patch.object(sys, "argv", ["sample.py", "--segment", segment]):
                    module.main()
            finally:
                module.TARGET_SAMPLE_SIZE = original_target

            self.assertIn("EPI_ISL_1254", module.REQUIRED_ISOLATE_IDS)

            sampled_metadata = pd.read_excel(
                segment_dir / f"{segment}_metadata_5000sample.xlsx"
            )
            sampled_ids = sampled_metadata["Isolate_Id"].tolist()
            self.assertIn("EPI_ISL_1254", sampled_ids)


if __name__ == "__main__":
    unittest.main()
