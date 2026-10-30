param(
    [string]$BuildDirectory
)

$ErrorActionPreference = 'Stop'
if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) {
    throw 'The recorded native kernels currently require Windows x64 and MSVC.'
}
$repositoryRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$sourceDirectory = Join-Path $repositoryRoot 'src\packed_tanimoto'
$packedSource = Join-Path $sourceDirectory 'packed_kernel.c'
$originalSource = Join-Path $sourceDirectory 'spsubmatxvec.c'
foreach ($source in @($packedSource, $originalSource)) {
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
foreach ($name in @('packed_kernel.dll', 'packed_kernel.obj', 'packed_kernel.exp', 'packed_kernel.lib',
                    'original_spmv.dll', 'original_spmv.obj', 'original_spmv.exp', 'original_spmv.lib',
                    'packed_build.log', 'original_build.log', 'build_metadata.json')) {
    if (Test-Path -LiteralPath (Join-Path $buildPath $name)) {
        throw "Build output already exists; choose a fresh directory: $buildPath"
    }
}
$null = New-Item -ItemType Directory -Path $buildPath -Force
$sourceHashes = @{
    packed_kernel = (Get-FileHash -LiteralPath $packedSource -Algorithm SHA256).Hash.ToLowerInvariant()
    original_spmv = (Get-FileHash -LiteralPath $originalSource -Algorithm SHA256).Hash.ToLowerInvariant()
}
$commonFlags = @('/nologo', '/O2', '/Oi', '/MT', '/LD', '/W4', '/TC')
$packedArguments = $commonFlags + @($packedSource, '/Fopacked_kernel.obj', '/Fepacked_kernel.dll', '/link', '/INCREMENTAL:NO')
$originalArguments = $commonFlags + @($originalSource, '/Fooriginal_spmv.obj', '/Feoriginal_spmv.dll', '/link', '/EXPORT:multiply', '/INCREMENTAL:NO')
Push-Location -LiteralPath $buildPath
try {
    & $compiler @packedArguments 2>&1 | Tee-Object -FilePath 'packed_build.log'
    if ($LASTEXITCODE -ne 0) { throw 'Packed kernel compilation failed; inspect packed_build.log.' }
    & $compiler @originalArguments 2>&1 | Tee-Object -FilePath 'original_build.log'
    if ($LASTEXITCODE -ne 0) { throw 'Original sparse kernel compilation failed; inspect original_build.log.' }
    foreach ($entry in @(@{ Name = 'packed_kernel'; Source = $packedSource },
                          @{ Name = 'original_spmv'; Source = $originalSource })) {
        $after = (Get-FileHash -LiteralPath $entry.Source -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($after -ne $sourceHashes[$entry.Name]) { throw "C source changed during compilation: $($entry.Source)" }
    }
    $installed = @()
    foreach ($name in @('packed_kernel.dll', 'original_spmv.dll')) {
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
        original_arguments = $originalArguments
        packed_source_sha256 = $sourceHashes.packed_kernel
        original_source_sha256 = $sourceHashes.original_spmv
        packed_dll_sha256 = (Get-FileHash -LiteralPath 'packed_kernel.dll' -Algorithm SHA256).Hash.ToLowerInvariant()
        original_dll_sha256 = (Get-FileHash -LiteralPath 'original_spmv.dll' -Algorithm SHA256).Hash.ToLowerInvariant()
        installed_files = $installed
        note = 'Built from the published C sources; historical benchmark binaries are not inputs.'
    } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath 'build_metadata.json' -Encoding UTF8
    Write-Output "Built and installed both x64 DLLs. Build record: $(Join-Path $buildPath 'build_metadata.json')"
} finally {
    Pop-Location
}
