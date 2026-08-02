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
    $RequirementsPath = Join-Path $resolvedProjectRoot "requirements-analysis.txt"
}
$resolvedRequirements = (Resolve-Path -LiteralPath $RequirementsPath).Path
$requiredModules = @(
    Get-Content -LiteralPath $resolvedRequirements |
        ForEach-Object { $_.Trim() } |
        Where-Object { $_ -and -not $_.StartsWith("#") }
)
if ($requiredModules.Count -eq 0) {
    throw "No Python dependencies were declared in $resolvedRequirements"
}

$environmentCreated = $false
$venvPython = Join-Path $resolvedProjectRoot ".venv\Scripts\python.exe"
if ($PythonPath) {
    if (-not (Test-PythonCommand -Candidate $PythonPath)) {
        throw "The requested Python is not usable: $PythonPath"
    }
    $resolvedPython = (Resolve-Path -LiteralPath $PythonPath).Path
}
elseif (Test-PythonCommand -Candidate $venvPython) {
    $resolvedPython = (Resolve-Path -LiteralPath $venvPython).Path
}
elseif ($CheckOnly) {
    $resolvedPython = Find-BasePython -Root $resolvedProjectRoot
}
else {
    $basePython = Find-BasePython -Root $resolvedProjectRoot
    & $basePython -m venv --system-site-packages (Join-Path $resolvedProjectRoot ".venv")
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to create the worktree .venv with $basePython"
    }
    $environmentCreated = $true
    $resolvedPython = (Resolve-Path -LiteralPath $venvPython).Path
}

$missingModules = @(Get-MissingPythonModules -Python $resolvedPython -Modules $requiredModules)
if ($missingModules.Count -gt 0 -and -not $CheckOnly) {
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

$ready = $missingModules.Count -eq 0
$report = [ordered]@{
    ready = $ready
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
    $state = if ($ready) { "READY" } else { "MISSING DEPENDENCIES" }
    Write-Host "Sandbox initialization: $state"
    Write-Host "  project: $resolvedProjectRoot"
    Write-Host "  branch: $gitBranch"
    Write-Host "  git common dir: $gitCommonDirectory"
    Write-Host "  git metadata external: $gitCommonExternal"
    Write-Host "  python: $resolvedPython"
    Write-Host "  missing modules: $($missingModules -join ', ')"
    Write-Host "  data root: $resolvedDataRoot"
}

if (-not $ready) {
    throw "Python dependencies are missing: $($missingModules -join ', ')"
}
