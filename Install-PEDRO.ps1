# Install/update PEDRO from GitHub tags in a private environment. PowerShell 5.1+; Python 3.11+; Git.
[CmdletBinding()]
param(
    [ValidateSet('Install', 'Update', 'Launch')][string]$Action = 'Install',
    [string]$InstallDir = (Join-Path $env:LOCALAPPDATA 'PEDRO'),
    [string]$Python,
    [ValidateSet('beta', 'stable')][string]$Channel,
    [string]$Version,
    [switch]$Unpin,
    [switch]$SkipUpdate,
    [string]$Repository,
    [switch]$NoDeps,
    [switch]$NoBuildIsolation
)
$ErrorActionPreference = 'Stop'
$root = [IO.Path]::GetFullPath($InstallDir)
$configPath = Join-Path $root 'install-config.json'
$statePath = Join-Path $root 'install-state.json'

function Write-JsonAtomic($Path, $Value) {
    $temporary = "$Path.$([Guid]::NewGuid().ToString('N')).tmp"
    $Value | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $temporary -Encoding UTF8
    Move-Item -LiteralPath $temporary -Destination $Path -Force
}

function Invoke-Checked($Executable, [string[]]$Arguments) {
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Command failed ($LASTEXITCODE): $Executable $($Arguments -join ' ')" }
}

function Select-ReleaseTag($Repo, $SelectedChannel, $Pin) {
    $env:GIT_TERMINAL_PROMPT = '0'
    $lines = & git -c http.lowSpeedLimit=1 -c http.lowSpeedTime=8 ls-remote --tags $Repo
    if ($LASTEXITCODE -ne 0) { throw 'Could not read release tags from GitHub.' }
    $tags = @{}
    foreach ($line in $lines) {
        if ($line -match '^([0-9a-f]{40,64})\s+refs/tags/(v?(\d+)\.(\d+)\.(\d+)(?:(a|b|rc)(\d+))?)(\^\{\})?$') {
            $tag = $Matches[2]
            $suffix = $Matches[6]
            $stage = switch ($suffix) { 'a' { 0 }; 'b' { 1 }; 'rc' { 2 }; default { 3 } }
            $pre = if ($Matches[7]) { [int]$Matches[7] } else { 0 }
            $peeled = [bool]$Matches[8]
            if (-not $tags.ContainsKey($tag) -or $peeled) {
                $tags[$tag] = [pscustomobject]@{ Tag=$tag; Version=($tag -replace '^v',''); Revision=$Matches[1];
                    Major=[int]$Matches[3]; Minor=[int]$Matches[4]; Patch=[int]$Matches[5]; Stage=$stage; Pre=$pre }
            }
        }
    }
    $candidates = @($tags.Values)
    if ($Pin) { $candidates = @($candidates | Where-Object { $_.Version -eq ($Pin -replace '^v','') }) }
    elseif ($SelectedChannel -eq 'stable') { $candidates = @($candidates | Where-Object { $_.Stage -eq 3 }) }
    $selected = $candidates | Sort-Object Major,Minor,Patch,Stage,Pre -Descending | Select-Object -First 1
    if (-not $selected) { throw "No matching published version tag exists (channel=$SelectedChannel, pin=$Pin)." }
    return $selected
}

try {
    $config = if (Test-Path -LiteralPath $configPath) { Get-Content -LiteralPath $configPath -Raw | ConvertFrom-Json } else { $null }
    $state = if (Test-Path -LiteralPath $statePath) { Get-Content -LiteralPath $statePath -Raw | ConvertFrom-Json } else { $null }
    if (-not $config -and $Action -ne 'Install') { throw 'Run the Install action first.' }
    if (-not $config) {
        if ((Test-Path -LiteralPath $root) -and @(Get-ChildItem -LiteralPath $root -Force).Count) {
            throw 'The target is not an empty folder or an existing managed PEDRO installation; choose -InstallDir.'
        }
        if (-not $Python) { $Python = 'python' }
        $pythonPath = & $Python -X utf8 -c "import sys; assert sys.version_info >= (3,11), 'Python 3.11+ required'; print(sys.executable)"
        if ($LASTEXITCODE -ne 0) { throw 'Python 3.11 or later is required; use -Python with its full path.' }
        $config = [pscustomobject]@{ Python=([string]$pythonPath).Trim(); Repository='https://github.com/Fisicas/PEDRO.git'; Channel='beta'; Pin='' }
    }
    if ($Repository) {
        if ($state -and $Repository -ne $config.Repository) { throw 'Existing installations keep their configured repository; use a new -InstallDir for another repository.' }
        $config.Repository = $Repository
    }
    if ($Python) {
        $pythonPath = & $Python -X utf8 -c "import sys; assert sys.version_info >= (3,11); print(sys.executable)"
        if ($LASTEXITCODE -ne 0) { throw 'The selected Python is unavailable or older than 3.11.' }
        $config.Python = ([string]$pythonPath).Trim()
    }
    if ($Channel) { $config.Channel = $Channel }
    if ($Version) { $config.Pin = $Version -replace '^v','' }
    if ($Unpin) { $config.Pin = '' }
    New-Item -ItemType Directory -Path $root -Force | Out-Null

    if (-not $SkipUpdate) {
        try {
            if (-not (Get-Command git -ErrorAction SilentlyContinue)) { throw 'Git is required to install or update from GitHub.' }
            $release = Select-ReleaseTag $config.Repository $config.Channel $config.Pin
            if (-not $state -or $state.Revision -ne $release.Revision -or -not (Test-Path -LiteralPath (Join-Path $state.Environment 'Scripts/python.exe'))) {
                $environment = Join-Path $root ("versions/{0}-{1}-{2}" -f $release.Version, $release.Revision.Substring(0,12), [Guid]::NewGuid().ToString('N').Substring(0,8))
                Write-Host "Installing PEDRO $($release.Version) from $($release.Tag) ($($release.Revision))..."
                Invoke-Checked $config.Python @('-X','utf8','-m','venv',$environment)
                $envPython = Join-Path $environment 'Scripts/python.exe'
                $requirement = "pedro-ped[gui] @ git+$($config.Repository)@$($release.Revision)"
                $pipArguments = @('-X','utf8','-m','pip','install','--ignore-installed','--disable-pip-version-check',$requirement)
                if ($NoDeps) { $pipArguments += '--no-deps' }
                if ($NoBuildIsolation) { $pipArguments += '--no-build-isolation' }
                Invoke-Checked $envPython $pipArguments
                Invoke-Checked $envPython @('-X','utf8','-c',"import sys,pathlib,pedro; from PySide6 import QtWidgets; import matplotlib; assert pedro.__version__ == '$($release.Version)', pedro.__version__; assert pathlib.Path(pedro.__file__).is_relative_to(pathlib.Path(sys.prefix))")
                $state = [pscustomobject]@{ Version=$release.Version; Tag=$release.Tag; Revision=$release.Revision;
                    Environment=$environment; InstalledAt=(Get-Date).ToString('o') }
                # Switch the active installation only after installation and import checks succeed.
                Write-JsonAtomic $statePath $state
            } else { Write-Host "PEDRO $($state.Version) is up to date for the selected channel/version." }
            Write-JsonAtomic $configPath $config
        } catch {
            if ($Action -ne 'Launch' -or -not $state) { throw }
            Write-Warning "Update unavailable; starting the existing installation. $($_.Exception.Message)"
        }
    }
    if (-not $state) { throw 'No installed version is available.' }

    # Keep a self-contained updater and launcher beside the private environments, outside any research data.
    $driver = Join-Path $root 'Install-PEDRO.ps1'
    if ([IO.Path]::GetFullPath($PSCommandPath) -ne $driver) { Copy-Item -LiteralPath $PSCommandPath -Destination $driver -Force }
    $launcher = "@echo off`r`npowershell.exe -NoProfile -ExecutionPolicy Bypass -File `"%~dp0Install-PEDRO.ps1`" -Action Launch -InstallDir `"%~dp0.`" %*`r`n"
    Set-Content -LiteralPath (Join-Path $root 'PEDRO.cmd') -Value $launcher -Encoding ASCII
    $updateCommand = $launcher -replace '-Action Launch', '-Action Update'
    Set-Content -LiteralPath (Join-Path $root 'Update-PEDRO.cmd') -Value $updateCommand -Encoding ASCII
    Write-Host "Launch: $(Join-Path $root 'PEDRO.cmd')"
    Write-Host "Installed version: $($state.Version); Git revision: $($state.Revision)"
    if ($Action -eq 'Launch') {
        $pythonw = Join-Path $state.Environment 'Scripts/pythonw.exe'
        if (-not (Test-Path -LiteralPath $pythonw)) { throw 'Installed Python is missing; run Update to repair the installation.' }
        Start-Process -FilePath $pythonw -ArgumentList @('-X','utf8','-m','pedro','gui')
    }
} catch {
    Write-Error $_
    exit 1
}
