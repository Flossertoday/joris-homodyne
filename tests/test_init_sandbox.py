import json
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INIT_SCRIPT = PROJECT_ROOT / "scripts" / "init_sandbox.ps1"
REQUIREMENTS = PROJECT_ROOT / "requirements.txt"
LOCAL_ENVIRONMENT = PROJECT_ROOT / ".codex" / "environments" / "environment.toml"


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
    def test_codex_local_environment_runs_initializer_automatically(self):
        content = LOCAL_ENVIRONMENT.read_text(encoding="utf-8")

        self.assertIn("version = 1", content)
        self.assertIn('name = "Joris homodyne"', content)
        self.assertIn("[setup.win32]", content)
        self.assertIn(
            'script = "powershell -NoProfile -ExecutionPolicy Bypass -File '
            '.\\\\scripts\\\\init_sandbox.ps1"',
            content,
        )

    def test_analysis_requirements_are_machine_readable(self):
        packages = [
            line.strip()
            for line in REQUIREMENTS.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]

        self.assertEqual(
            [
                "numpy==1.26.4",
                "scipy==1.11.4",
                "h5py==3.16.0",
                "matplotlib==3.8.0",
                "trackpy==0.7",
                "tomli==2.4.1",
            ],
            packages,
        )
        module_names = [re.split(r"[=<>\[;\s!]", line)[0] for line in packages]
        self.assertEqual(
            ["numpy", "scipy", "h5py", "matplotlib", "trackpy", "tomli"], module_names
        )

    def test_initializer_runs_under_windows_powershell(self):
        executable = shutil.which("powershell")
        if executable is None:
            raise AssertionError("Windows PowerShell is required by environment.toml")

        result = subprocess.run(
            [
                executable,
                "-NoProfile",
                "-ExecutionPolicy",
                "Bypass",
                "-File",
                str(INIT_SCRIPT),
                "-CheckOnly",
                "-AsJson",
            ],
            text=True,
            capture_output=True,
        )

        self.assertEqual(0, result.returncode, result.stderr)
        report = json.loads(result.stdout)
        self.assertEqual(str(PROJECT_ROOT.resolve()), report["project_root"])
        if (PROJECT_ROOT / ".venv" / "Scripts" / "python.exe").is_file():
            self.assertEqual("ready", report["state"])
            self.assertTrue(report["ready"])
        else:
            # Fresh checkout without .venv: readiness must be reported as
            # "not_initialized", not judged from the base interpreter's
            # site-packages — under the isolated-.venv architecture those
            # packages say nothing about the project environment.
            self.assertEqual("not_initialized", report["state"])
            self.assertFalse(report["ready"])
            self.assertEqual([], report["dependencies"]["missing"])

    def test_check_only_without_venv_reports_not_initialized(self):
        with tempfile.TemporaryDirectory() as tmp:
            project_root = Path(tmp)
            data_root = project_root / "DATA"
            data_root.mkdir()
            requirements = project_root / "requirements.txt"
            # trackpy is intentionally included: the base sci interpreter does
            # not provide it, and CheckOnly must still report not_initialized
            # instead of failing on the base environment's package list.
            requirements.write_text(
                "numpy==1.26.4\nscipy==1.11.4\nh5py==3.16.0\nmatplotlib==3.8.0\ntrackpy==0.7\n",
                encoding="utf-8",
            )
            command = " ".join(
                [
                    "&",
                    powershell_quote(INIT_SCRIPT),
                    "-ProjectRoot",
                    powershell_quote(project_root),
                    "-RequirementsPath",
                    powershell_quote(requirements),
                    "-CheckOnly",
                    "-AsJson",
                ]
            )

            result = run_powershell(command)

            self.assertEqual(0, result.returncode, result.stderr)
            report = json.loads(result.stdout)
            self.assertEqual("not_initialized", report["state"])
            self.assertFalse(report["ready"])
            self.assertFalse(report["environment_created"])
            self.assertEqual([], report["dependencies"]["missing"])
            self.assertEqual(str(project_root.resolve()), report["project_root"])

    def test_check_only_reports_python_dependencies_and_local_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            project_root = Path(tmp)
            data_root = project_root / "DATA"
            data_root.mkdir()
            requirements = project_root / "requirements.txt"
            requirements.write_text(
                "numpy==1.26.4\nscipy==1.11.4\nh5py==3.16.0\nmatplotlib==3.8.0\n",
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
            self.assertFalse((project_root / "tmp" / "sandbox-initialization.json").exists())

    def test_setup_persists_report_for_later_agent_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            project_root = Path(tmp)
            data_root = project_root / "DATA"
            data_root.mkdir()
            requirements = project_root / "requirements.txt"
            requirements.write_text(
                "numpy==1.26.4\nscipy==1.11.4\nh5py==3.16.0\nmatplotlib==3.8.0\n",
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
                    "-AsJson",
                ]
            )

            result = run_powershell(command)

            self.assertEqual(0, result.returncode, result.stderr)
            output_report = json.loads(result.stdout)
            state_path = project_root / "tmp" / "sandbox-initialization.json"
            self.assertTrue(state_path.is_file())
            persisted_report = json.loads(state_path.read_text(encoding="utf-8"))
            self.assertEqual(output_report, persisted_report)
            self.assertEqual(str(Path(sys.executable).resolve()), persisted_report["python"])
            self.assertEqual(str(data_root.resolve()), persisted_report["data_root"])

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
