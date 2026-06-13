$ErrorActionPreference = "Stop"

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$raw = Join-Path $root "data\raw"
New-Item -ItemType Directory -Force -Path $raw | Out-Null

$files = @(
  @{
    Name = "Supplementary Data 1.xlsx"
    Url = "https://zenodo.org/records/8223844/files/Supplementary%20Data%201.xlsx?download=1"
  },
  @{
    Name = "Supplementary Data 2.xlsx"
    Url = "https://zenodo.org/records/8223844/files/Supplementary%20Data%202.xlsx?download=1"
  },
  @{
    Name = "Supplementary Data 2 subgroup.xlsx"
    Url = "https://zenodo.org/records/8223844/files/Supplementary%20Data%202%20subgroup.xlsx?download=1"
  }
)

foreach ($f in $files) {
  $dest = Join-Path $raw $f.Name
  if (Test-Path $dest) {
    Write-Host "Exists: $dest"
    continue
  }
  Write-Host "Downloading $($f.Name)"
  Invoke-WebRequest -Uri $f.Url -OutFile $dest
}

Write-Host "SMAG metadata download complete."
