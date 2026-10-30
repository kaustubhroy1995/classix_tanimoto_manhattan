param(
    [string]$BuildDirectory
)

$ErrorActionPreference = 'Stop'
if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) {
    throw 'The recorded native kernels currently require Windows x64 and MSVC.'
}
$repositoryRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$sourceDirectory = Join-Path $repositoryRoot 'src\spmv'
$packedSource = Join-Path $sourceDirectory 'spsubmatxvec.c'
foreach ($source in @($packedSource)) {
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
        throw "Required C source is missing: $source"
    }
}

function Import-X64CompilerEnvironment {
    $existingCompiler = Get-Command 'cl.exe' -CommandType Application -ErrorAction SilentlyContinue
    if ($existingCompiler -and $env:VSCMD_ARG_TGT_ARCH -eq 'x64' -and $env:INCLUDE -and $env:LIB) {
        return
    }
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere -PathType Leaf)) {
        throw 'Install Visual Studio C++ Build Tools and a Windows SDK, or run from an x64 Developer PowerShell.'
    }
    $installation = (& $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath | Select-Object -First 1)
    if ($LASTEXITCODE -ne 0 -or -not $installation) {
        throw 'No Visual Studio installation with the x64 C++ toolset was found.'
    }
    $developerScript = Join-Path $installation 'Common7\Tools\VsDevCmd.bat'
    if (-not (Test-Path -LiteralPath $developerScript -PathType Leaf)) {
        throw "Visual Studio developer environment is missing: $developerScript"
    }
    # cmd is used only to read the environment produced by the installed VS script.
    $environmentCommand = 'call "' + $developerScript + '" -no_logo -arch=x64 -host_arch=x64 >nul && set'
    $environmentLines = & $env:ComSpec /d /s /c $environmentCommand
    if ($LASTEXITCODE -ne 0) { throw 'Could not initialize the x64 Visual Studio environment.' }
    foreach ($line in $environmentLines) {
        if ($line -match '^([^=]+)=(.*)$') {
            [Environment]::SetEnvironmentVariable($Matches[1], $Matches[2], 'Process')
        }
    }
}

Import-X64CompilerEnvironment
$compilerCommand = Get-Command 'cl.exe' -CommandType Application -ErrorAction Stop
if ($env:VSCMD_ARG_TGT_ARCH -ne 'x64' -or -not $env:INCLUDE -or -not $env:LIB) {
    throw 'An x64 MSVC environment with Windows SDK INCLUDE and LIB paths is required.'
}
$compiler = $compilerCommand.Source
if (-not $BuildDirectory) {
    $stamp = [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffZ')
    $BuildDirectory = Join-Path $repositoryRoot (Join-Path '_build\native' $stamp)
} elseif (-not [IO.Path]::IsPathRooted($BuildDirectory)) {
    $BuildDirectory = Join-Path $repositoryRoot $BuildDirectory
}
$buildPath = [IO.Path]::GetFullPath($BuildDirectory)
if ($buildPath.TrimEnd('\') -eq $sourceDirectory.TrimEnd('\')) {
    throw 'Choose a separate build directory; generated intermediate files must not replace sources.'
}
foreach ($name in @('spsubmatxvec.dll', 'spsubmatxvec.obj', 'spsubmatxvec.exp', 'spsubmatxvec.lib',
                    'build.log', 'build_metadata.json')) {
    if (Test-Path -LiteralPath (Join-Path $buildPath $name)) {
        throw "Build output already exists; choose a fresh directory: $buildPath"
    }
}
$null = New-Item -ItemType Directory -Path $buildPath -Force
$sourceHashes = @{
    spsubmatxvec = (Get-FileHash -LiteralPath $packedSource -Algorithm SHA256).Hash.ToLowerInvariant()
}
$commonFlags = @('/nologo', '/O2', '/Oi', '/MT', '/LD', '/W4', '/TC')
$packedArguments = $commonFlags + @($packedSource, '/Fospsubmatxvec.obj', '/Fespsubmatxvec.dll', '/link', '/INCREMENTAL:NO')
Push-Location -LiteralPath $buildPath
try {
    & $compiler @packedArguments 2>&1 | Tee-Object -FilePath 'build.log'
    if ($LASTEXITCODE -ne 0) { throw 'Tanimoto kernel compilation failed; inspect build.log.' }
    foreach ($entry in @(@{ Name = 'spsubmatxvec'; Source = $packedSource })) {
        $after = (Get-FileHash -LiteralPath $entry.Source -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($after -ne $sourceHashes[$entry.Name]) { throw "C source changed during compilation: $($entry.Source)" }
    }
    $installed = @()
    foreach ($name in @('spsubmatxvec.dll')) {
        $destination = Join-Path $sourceDirectory $name
        Copy-Item -LiteralPath (Join-Path $buildPath $name) -Destination $destination -Force
        $installed += $destination
    }
    [ordered]@{
        built_at_utc = [DateTime]::UtcNow.ToString('o')
        compiler = $compiler
        compiler_file_version = (Get-Item -LiteralPath $compiler).VersionInfo.FileVersion
        toolset = $env:VCToolsVersion
        windows_sdk = $env:WindowsSDKVersion
        architecture = 'x64'
        common_flags = $commonFlags
        packed_arguments = $packedArguments
        source_sha256 = $sourceHashes.spsubmatxvec
        dll_sha256 = (Get-FileHash -LiteralPath 'spsubmatxvec.dll' -Algorithm SHA256).Hash.ToLowerInvariant()
        installed_files = $installed
        note = 'Built from the published C sources; historical benchmark binaries are not inputs.'
    } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath 'build_metadata.json' -Encoding UTF8
    Write-Output "Built and installed the x64 Tanimoto library. Build record: $(Join-Path $buildPath 'build_metadata.json')"
} finally {
    Pop-Location
}
