param([Parameter(Mandatory=$true)][string]$OutputDirectory)
$ErrorActionPreference = 'Stop'
$studyDir = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$nativeDir = $OutputDirectory
$source = Join-Path $PSScriptRoot 'silhouette_kernel.c'
$vswhere = 'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
$vsPath = (& $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath)
if (-not $vsPath) { throw 'A local MSVC x64 toolchain was not found.' }
$toolset = Get-ChildItem -LiteralPath (Join-Path $vsPath 'VC\Tools\MSVC') -Directory | Sort-Object Name -Descending | Select-Object -First 1
$compilerDir = Join-Path $toolset.FullName 'bin\Hostx64\x64'
$compiler = Join-Path $compilerDir 'cl.exe'
$kitRoot = 'C:\Program Files (x86)\Windows Kits\10'
$kitVersion = Get-ChildItem -LiteralPath (Join-Path $kitRoot 'Include') -Directory | Sort-Object Name -Descending | Select-Object -First 1
$includes = @((Join-Path $toolset.FullName 'include'), (Join-Path $kitVersion.FullName 'ucrt'), (Join-Path $kitVersion.FullName 'shared'), (Join-Path $kitVersion.FullName 'um'))
$libs = @((Join-Path $toolset.FullName 'lib\x64'), (Join-Path $kitRoot ('Lib\' + $kitVersion.Name + '\ucrt\x64')), (Join-Path $kitRoot ('Lib\' + $kitVersion.Name + '\um\x64')))
New-Item -ItemType Directory -Path $nativeDir -Force | Out-Null
$savedPath = $env:PATH
$savedInclude = $env:INCLUDE
$savedLib = $env:LIB
try {
    $env:PATH = $compilerDir + ';' + $savedPath
    $env:INCLUDE = $includes -join ';'
    $env:LIB = $libs -join ';'
    Push-Location -LiteralPath $nativeDir
    try {
        $flags = @('/nologo', '/O2', '/Oi', '/MT', '/LD', '/W4', '/TC', $source, '/Fosilhouette_kernel.obj', '/Fesilhouette_kernel.dll', '/link', '/INCREMENTAL:NO')
        & $compiler @flags 2>&1 | Tee-Object -FilePath 'build.log'
        if ($LASTEXITCODE -ne 0) { throw 'Silhouette kernel compilation failed.' }
        [ordered]@{
            built_at_utc = [DateTime]::UtcNow.ToString('o')
            compiler = $compiler
            compiler_file_version = (Get-Item -LiteralPath $compiler).VersionInfo.FileVersion
            toolset = $toolset.Name
            windows_sdk = $kitVersion.Name
            architecture = 'x64'
            flags = $flags
            source_sha256 = (Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash.ToLowerInvariant()
            dll_sha256 = (Get-FileHash -LiteralPath 'silhouette_kernel.dll' -Algorithm SHA256).Hash.ToLowerInvariant()
            cpu_policy = 'Runtime POPCNT detection; portable software fallback'
            threading = 'serial; no OpenMP or BLAS'
        } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath 'build_metadata.json' -Encoding UTF8
    } finally { Pop-Location }
} finally {
    $env:PATH = $savedPath
    $env:INCLUDE = $savedInclude
    $env:LIB = $savedLib
}
