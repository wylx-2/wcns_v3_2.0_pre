param(
    [Parameter(Mandatory = $true)]
    [string]$CacheDirectory
)

$ErrorActionPreference = "Stop"
$repository = Split-Path -Parent $PSScriptRoot
$manifestPath = Join-Path $repository "cases/validation/sst_tmr/manifest.json"
$manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
$cache = [System.IO.Path]::GetFullPath($CacheDirectory)

foreach ($asset in $manifest.assets) {
    $target = Join-Path $cache $asset.path
    $parent = Split-Path -Parent $target
    New-Item -ItemType Directory -Force -Path $parent | Out-Null
    if (-not (Test-Path -LiteralPath $target -PathType Leaf)) {
        Write-Host "fetching $($asset.url)"
        curl.exe -L --fail --output $target $asset.url
    }
}

python (Join-Path $repository "tools/verify_sa_tmr_assets.py") `
    --manifest $manifestPath --cache $cache
if ($LASTEXITCODE -ne 0) {
    throw "SST TMR asset verification failed"
}
