[CmdletBinding()]
param(
    [string]$ProjectRoot,
    [string]$PythonPath,
    [string]$RequirementsPath,
    [string]$DataRoot,
    [switch]$CheckOnly,
    [switch]$AsJson
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if (-not $ProjectRoot) {
    $ProjectRoot = Split-Path -Parent $PSScriptRoot
}

function Add-ProcessGitSafeDirectory {
    param([Parameter(Mandatory)][string]$Directory)

    $safeDirectory = $Directory.Replace("\", "/")
    $count = 0
    if ($env:GIT_CONFIG_COUNT) {
        $parsed = 0
        if (-not [int]::TryParse($env:GIT_CONFIG_COUNT, [ref]$parsed)) {
            throw "GIT_CONFIG_COUNT is not an integer: $env:GIT_CONFIG_COUNT"
        }
        $count = $parsed
    }

    $alreadyConfigured = $false
    for ($index = 0; $index -lt $count; $index++) {
        $key = [Environment]::GetEnvironmentVariable("GIT_CONFIG_KEY_$index", "Process")
        $value = [Environment]::GetEnvironmentVariable("GIT_CONFIG_VALUE_$index", "Process")
        if ($key -eq "safe.directory" -and $value -eq $safeDirectory) {
            $alreadyConfigured = $true
            break
        }
    }

    if (-not $alreadyConfigured) {
        [Environment]::SetEnvironmentVariable(
            "GIT_CONFIG_KEY_$count",
            "safe.directory",
            "Process"
        )
        [Environment]::SetEnvironmentVariable(
            "GIT_CONFIG_VALUE_$count",
            $safeDirectory,
            "Process"
        )
        $env:GIT_CONFIG_COUNT = [string]($count + 1)
    }

    $env:JORIS_SAFE_DIRECTORY = $safeDirectory
    return $safeDirectory
}

function Test-PythonCommand {
    param([Parameter(Mandatory)][string]$Candidate)

    if (-not (Test-Path -LiteralPath $Candidate -PathType Leaf)) {
        return $false
    }
    & $Candidate -c "import sys; raise SystemExit(0)" 2>$null
    return $LASTEXITCODE -eq 0
}

function Test-VenvIsolated {
    param([Parameter(Mandatory)][string]$VenvDirectory)

    $pyvenvCfg = Join-Path $VenvDirectory "pyvenv.cfg"
    if (-not (Test-Path -LiteralPath $pyvenvCfg -PathType Leaf)) {
        return $false
    }
    foreach ($line in Get-Content -LiteralPath $pyvenvCfg) {
        if ($line -match '^\s*include-system-site-packages\s*=\s*true\s*$') {
            return $false
        }
    }
    return $true
}

function Find-BasePython {
    param([Parameter(Mandatory)][string]$Root)

    $candidates = @()
    if ($env:PROJECT_PYTHON) {
        $candidates += $env:PROJECT_PYTHON
    }
    if ($env:USERPROFILE) {
        $candidates += Join-Path $env:USERPROFILE ".conda\envs\sci\python.exe"
        $candidates += Join-Path $env:USERPROFILE "anaconda3\envs\sci\python.exe"
    }
    $pythonCommand = Get-Command python -ErrorAction SilentlyContinue
    if ($pythonCommand) {
        $candidates += $pythonCommand.Source
    }

    foreach ($candidate in $candidates | Select-Object -Unique) {
        if (Test-PythonCommand -Candidate $candidate) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    throw "No usable Python was found. Pass -PythonPath or expose the sci environment."
}

function Get-MissingPythonModules {
    param(
        [Parameter(Mandatory)][string]$Python,
        [Parameter(Mandatory)][string[]]$Modules
    )

    $missing = @()
    foreach ($module in $Modules) {
        & $Python -c (
            "import importlib.util, sys; " +
            "raise SystemExit(0 if importlib.util.find_spec(sys.argv[1]) else 1)"
        ) $module 2>$null
        if ($LASTEXITCODE -ne 0) {
            $missing += $module
        }
    }
    return [object[]]$missing
}

$resolvedProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$safeDirectory = Add-ProcessGitSafeDirectory -Directory $resolvedProjectRoot

if (-not $RequirementsPath) {
    $RequirementsPath = Join-Path $resolvedProjectRoot "requirements.txt"
}
$resolvedRequirements = (Resolve-Path -LiteralPath $RequirementsPath).Path
$requiredModules = @(
    Get-Content -LiteralPath $resolvedRequirements |
        ForEach-Object { $_.Trim() } |
        Where-Object { $_ -and -not $_.StartsWith("#") } |
        ForEach-Object { ($_ -split '[=<>!\[;\s]')[0] } |
        Where-Object { $_ }
)
if ($requiredModules.Count -eq 0) {
    throw "No Python dependencies were declared in $resolvedRequirements"
}

$environmentCreated = $false
$usedBasePythonFallback = $false
$venvDirectory = Join-Path $resolvedProjectRoot ".venv"
$venvPython = Join-Path $venvDirectory "Scripts\python.exe"
if ($PythonPath) {
    if (-not (Test-PythonCommand -Candidate $PythonPath)) {
        throw "The requested Python is not usable: $PythonPath"
    }
    $resolvedPython = (Resolve-Path -LiteralPath $PythonPath).Path
}
elseif ((Test-PythonCommand -Candidate $venvPython) -and (Test-VenvIsolated -VenvDirectory $venvDirectory)) {
    $resolvedPython = (Resolve-Path -LiteralPath $venvPython).Path
}
elseif ($CheckOnly) {
    if (Test-PythonCommand -Candidate $venvPython) {
        # CheckOnly stays read-only: report against the existing .venv even if
        # it still shares the base environment's site-packages.
        $resolvedPython = (Resolve-Path -LiteralPath $venvPython).Path
    }
    else {
        # No .venv yet: the base interpreter only donates the Python runtime a
        # future .venv is created from. Its own site-packages say nothing
        # about project readiness, so the module check below is skipped and
        # the report carries state "not_initialized".
        $resolvedPython = Find-BasePython -Root $resolvedProjectRoot
        $usedBasePythonFallback = $true
    }
}
else {
    # Recreate from scratch: the previous .venv is either unusable or still
    # leaks the base environment through --system-site-packages.
    if (Test-Path -LiteralPath $venvDirectory -PathType Container) {
        Remove-Item -LiteralPath $venvDirectory -Recurse -Force
    }
    $basePython = Find-BasePython -Root $resolvedProjectRoot
    & $basePython -m venv $venvDirectory
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to create the worktree .venv with $basePython"
    }
    $environmentCreated = $true
    $resolvedPython = (Resolve-Path -LiteralPath $venvPython).Path
}

$isWorktreeVenv = $resolvedPython.StartsWith($venvDirectory, [StringComparison]::OrdinalIgnoreCase)
if ($usedBasePythonFallback) {
    # Readiness of the project environment is unknown until a .venv exists;
    # probing the base interpreter's modules would conflate it with the
    # isolated project environment (see docs/sandbox-initialization.md).
    $missingModules = @()
}
else {
    $missingModules = @(Get-MissingPythonModules -Python $resolvedPython -Modules $requiredModules)
}
if (-not $CheckOnly -and $isWorktreeVenv) {
    # The worktree .venv is synced from the pinned requirements on every run;
    # pip no-ops locally when every pin is already satisfied. An explicit
    # external -PythonPath is never installed into.
    $installOutput = @(
        & $resolvedPython -m pip install --disable-pip-version-check -r $resolvedRequirements 2>&1
    )
    if ($LASTEXITCODE -ne 0) {
        throw "Dependency installation failed:`n$($installOutput -join "`n")"
    }
    if (-not $AsJson -and $installOutput.Count -gt 0) {
        $installOutput | ForEach-Object { Write-Host $_ }
    }
    $missingModules = @(
        Get-MissingPythonModules -Python $resolvedPython -Modules $requiredModules
    )
}

$pythonDirectory = Split-Path -Parent $resolvedPython
$pathEntries = @($env:PATH -split ";")
if ($pythonDirectory -notin $pathEntries) {
    $env:PATH = "$pythonDirectory;$env:PATH"
}
$env:PROJECT_PYTHON = $resolvedPython

$gitIsRepository = $false
$gitCommonDirectory = $null
$gitBranch = $null
$gitCommand = Get-Command git -ErrorAction SilentlyContinue
if ($gitCommand) {
    $commonOutput = @(& git -C $resolvedProjectRoot rev-parse --git-common-dir 2>$null)
    if ($LASTEXITCODE -eq 0 -and $commonOutput.Count -gt 0) {
        $gitIsRepository = $true
        $commonPath = $commonOutput[-1]
        if (-not [IO.Path]::IsPathRooted($commonPath)) {
            $commonPath = Join-Path $resolvedProjectRoot $commonPath
        }
        $gitCommonDirectory = (Resolve-Path -LiteralPath $commonPath).Path
        $branchOutput = @(& git -C $resolvedProjectRoot branch --show-current 2>$null)
        if ($LASTEXITCODE -eq 0 -and $branchOutput.Count -gt 0) {
            $gitBranch = $branchOutput[-1]
        }
    }
}

$resolvedDataRoot = $null
if ($DataRoot) {
    $resolvedDataRoot = (Resolve-Path -LiteralPath $DataRoot).Path
}
else {
    $worktreeData = Join-Path $resolvedProjectRoot "DATA"
    if (Test-Path -LiteralPath $worktreeData -PathType Container) {
        $resolvedDataRoot = (Resolve-Path -LiteralPath $worktreeData).Path
    }
    elseif ($gitCommonDirectory) {
        $mainProjectData = Join-Path (Split-Path -Parent $gitCommonDirectory) "DATA"
        if (Test-Path -LiteralPath $mainProjectData -PathType Container) {
            $resolvedDataRoot = (Resolve-Path -LiteralPath $mainProjectData).Path
        }
    }
}

if ($resolvedDataRoot) {
    $env:JORIS_DATA_ROOT = $resolvedDataRoot
}
else {
    Remove-Item Env:JORIS_DATA_ROOT -ErrorAction SilentlyContinue
}

$gitCommonExternal = $false
if ($gitCommonDirectory) {
    $rootWithSeparator = $resolvedProjectRoot.TrimEnd("\") + "\"
    $gitCommonExternal = -not $gitCommonDirectory.StartsWith(
        $rootWithSeparator,
        [StringComparison]::OrdinalIgnoreCase
    )
}

if ($missingModules.Count -gt 0) {
    $readinessState = "missing_dependencies"
}
elseif ($usedBasePythonFallback) {
    $readinessState = "not_initialized"
}
else {
    $readinessState = "ready"
}
$ready = $readinessState -eq "ready"
$report = [ordered]@{
    ready = $ready
    state = $readinessState
    project_root = $resolvedProjectRoot
    environment_created = $environmentCreated
    state_path = $null
    python = $resolvedPython
    dependencies = [ordered]@{
        requirements = $resolvedRequirements
        required = [object[]]$requiredModules
        missing = [object[]]$missingModules
    }
    data_root = $resolvedDataRoot
    data_available = [bool]$resolvedDataRoot
    git = [ordered]@{
        is_repository = $gitIsRepository
        branch = $gitBranch
        safe_directory = $safeDirectory
        common_directory = $gitCommonDirectory
        common_directory_external = $gitCommonExternal
        metadata_write_may_require_approval = $gitCommonExternal
    }
}

if (-not $CheckOnly) {
    $stateDirectory = Join-Path $resolvedProjectRoot "tmp"
    New-Item -ItemType Directory -Path $stateDirectory -Force | Out-Null
    $statePath = Join-Path $stateDirectory "sandbox-initialization.json"
    $partialStatePath = "$statePath.partial"
    $report["state_path"] = $statePath
    $stateJson = $report | ConvertTo-Json -Depth 5 -Compress
    [IO.File]::WriteAllText(
        $partialStatePath,
        $stateJson + [Environment]::NewLine,
        [Text.UTF8Encoding]::new($false)
    )
    Move-Item -LiteralPath $partialStatePath -Destination $statePath -Force
}

if ($AsJson) {
    $report | ConvertTo-Json -Depth 5 -Compress
}
else {
    $state = switch ($readinessState) {
        "ready" { "READY" }
        "not_initialized" { "NOT INITIALIZED (run without -CheckOnly to create .venv)" }
        default { "MISSING DEPENDENCIES" }
    }
    Write-Host "Sandbox initialization: $state"
    Write-Host "  project: $resolvedProjectRoot"
    Write-Host "  branch: $gitBranch"
    Write-Host "  git common dir: $gitCommonDirectory"
    Write-Host "  git metadata external: $gitCommonExternal"
    Write-Host "  python: $resolvedPython"
    Write-Host "  missing modules: $($missingModules -join ', ')"
    Write-Host "  data root: $resolvedDataRoot"
}

if ($readinessState -eq "missing_dependencies") {
    throw "Python dependencies are missing: $($missingModules -join ', ')"
}
