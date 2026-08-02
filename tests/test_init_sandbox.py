import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INIT_SCRIPT = PROJECT_ROOT / "scripts" / "init_sandbox.ps1"
REQUIREMENTS = PROJECT_ROOT / "requirements-analysis.txt"


def powershell_quote(value: Path | str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def run_powershell(command: str) -> subprocess.CompletedProcess[str]:
    executable = shutil.which("pwsh")
    if executable is None:
        raise AssertionError("pwsh is required to initialize the Windows sandbox")
    return subprocess.run(
        [
            executable,
            "-NoLogo",
            "-NoProfile",
            "-NonInteractive",
            "-Command",
            command,
        ],
        text=True,
        capture_output=True,
    )


class SandboxInitializationTests(unittest.TestCase):
    def test_analysis_requirements_are_machine_readable(self):
        packages = [
            line.strip()
            for line in REQUIREMENTS.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]

        self.assertEqual(["numpy", "scipy", "h5py", "matplotlib"], packages)

    def test_check_only_reports_python_dependencies_and_local_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            project_root = Path(tmp)
            data_root = project_root / "DATA"
            data_root.mkdir()
            requirements = project_root / "requirements-analysis.txt"
            requirements.write_text(
                "numpy\nscipy\nh5py\nmatplotlib\n",
                encoding="utf-8",
            )
            command = " ".join(
                [
                    "&",
                    powershell_quote(INIT_SCRIPT),
                    "-ProjectRoot",
                    powershell_quote(project_root),
                    "-PythonPath",
                    powershell_quote(sys.executable),
                    "-RequirementsPath",
                    powershell_quote(requirements),
                    "-CheckOnly",
                    "-AsJson",
                ]
            )

            result = run_powershell(command)

            self.assertEqual(0, result.returncode, result.stderr)
            report = json.loads(result.stdout)
            self.assertTrue(report["ready"])
            self.assertFalse(report["environment_created"])
            self.assertEqual([], report["dependencies"]["missing"])
            self.assertEqual(str(project_root.resolve()), report["project_root"])
            self.assertEqual(str(data_root.resolve()), report["data_root"])
            self.assertEqual(
                str(project_root.resolve()).replace("\\", "/"),
                report["git"]["safe_directory"],
            )
            self.assertFalse(report["git"]["is_repository"])

    def test_process_environment_makes_plain_git_and_python_reusable(self):
        command = f"""
& {powershell_quote(INIT_SCRIPT)} `
  -ProjectRoot {powershell_quote(PROJECT_ROOT)} `
  -PythonPath {powershell_quote(sys.executable)} `
  -RequirementsPath {powershell_quote(REQUIREMENTS)} `
  -CheckOnly `
  -AsJson | Out-Null
$gitOutput = @(git status --short --branch 2>&1) -join "`n"
$gitExit = $LASTEXITCODE
[pscustomobject]@{{
  git_exit = $gitExit
  git_output = $gitOutput
  project_python = $env:PROJECT_PYTHON
  data_root = $env:JORIS_DATA_ROOT
  safe_directory = $env:JORIS_SAFE_DIRECTORY
}} | ConvertTo-Json -Compress
"""

        result = run_powershell(command)

        self.assertEqual(0, result.returncode, result.stderr)
        report = json.loads(result.stdout)
        self.assertEqual(0, report["git_exit"], report["git_output"])
        self.assertEqual(str(Path(sys.executable).resolve()), report["project_python"])
        self.assertEqual(
            str(PROJECT_ROOT.resolve()).replace("\\", "/"),
            report["safe_directory"],
        )
        self.assertTrue(Path(report["data_root"]).is_dir())
        self.assertEqual("DATA", Path(report["data_root"]).name)


if __name__ == "__main__":
    unittest.main()
