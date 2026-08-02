import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CURRENT_RAW_DIRECTORY = "DATA/26-07-31-single-sweep"
CURRENT_DERIVED_DIRECTORY = "DATA/26-07-31-derived-single-sweep"
CURRENT_ANALYSIS_DIRECTORY = "analysis/26-07-31-derived-single-sweep"
LEGACY_DIRECTORIES = (
    "DATA/26-07-31 single sweep",
    "DATA/26-07-31 single sweep derived",
)
LEGACY_ANALYSIS_PATHS = (
    "analysis/README.md",
    "analysis/sum_streak_condition.py",
    "analysis/results/26-07-31_central_dot_0p18_sum.png",
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

    def test_maintained_documents_use_experiment_scoped_analysis_paths(self):
        for relative_path in ("README.md", "MEMORY.md"):
            with self.subTest(relative_path=relative_path):
                content = (PROJECT_ROOT / relative_path).read_text(encoding="utf-8")
                self.assertIn(CURRENT_ANALYSIS_DIRECTORY, content)
                for legacy_path in LEGACY_ANALYSIS_PATHS:
                    self.assertNotIn(legacy_path, content)

    def test_analysis_readme_commands_use_experiment_scoped_paths(self):
        readme = PROJECT_ROOT / CURRENT_ANALYSIS_DIRECTORY / "README.md"
        content = readme.read_text(encoding="utf-8")

        self.assertIn(
            f"{CURRENT_ANALYSIS_DIRECTORY}/sum_streak_condition.py",
            content,
        )
        self.assertIn(
            f"{CURRENT_ANALYSIS_DIRECTORY}/results/"
            "26-07-31_central_dot_0p18_sum.png",
            content,
        )
        for legacy_path in LEGACY_ANALYSIS_PATHS[1:]:
            self.assertNotIn(legacy_path, content)

    def test_maintained_environment_docs_use_the_sandbox_initializer(self):
        initializer = "scripts/init_sandbox.ps1"
        for relative_path in (
            "AGENTS.md",
            "README.md",
            "docs/data.md",
            "docs/sandbox-initialization.md",
        ):
            with self.subTest(relative_path=relative_path):
                content = (PROJECT_ROOT / relative_path).read_text(encoding="utf-8")
                self.assertIn(initializer, content)


if __name__ == "__main__":
    unittest.main()
