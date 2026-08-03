import shutil
import subprocess
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "kimi-commit.sh"


class KimiCommitScriptTests(unittest.TestCase):
    def test_script_passes_bash_syntax_check(self):
        executable = shutil.which("bash")
        if executable is None:
            raise AssertionError("bash is required to lint scripts/kimi-commit.sh")

        result = subprocess.run(
            [executable, "-n", str(SCRIPT)],
            text=True,
            capture_output=True,
        )

        self.assertEqual(0, result.returncode, result.stderr)

    def test_script_stamps_kimi_identity_and_trailers(self):
        content = SCRIPT.read_text(encoding="utf-8")

        self.assertIn('user.name="Kimi K3"', content)
        self.assertIn('user.email="kimi-k3@users.noreply.github.com"', content)
        self.assertIn('--trailer "Model:', content)
        self.assertIn('--trailer "Thinking:', content)


if __name__ == "__main__":
    unittest.main()
