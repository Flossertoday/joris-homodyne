import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CURRENT_RAW_DIRECTORY = "DATA/26-07-31-single-sweep"
CURRENT_DERIVED_DIRECTORY = "DATA/26-07-31-derived-single-sweep"
LEGACY_DIRECTORIES = (
    "DATA/26-07-31 single sweep",
    "DATA/26-07-31 single sweep derived",
)


class ProjectDataPathTests(unittest.TestCase):
    def test_maintained_documents_use_current_data_directory_names(self):
        for relative_path in ("README.md", "docs/data.md", "MEMORY.md"):
            with self.subTest(relative_path=relative_path):
                content = (PROJECT_ROOT / relative_path).read_text(encoding="utf-8")
                self.assertIn(CURRENT_RAW_DIRECTORY, content)
                self.assertIn(CURRENT_DERIVED_DIRECTORY, content)
                for legacy_directory in LEGACY_DIRECTORIES:
                    self.assertNotIn(legacy_directory, content)


if __name__ == "__main__":
    unittest.main()
