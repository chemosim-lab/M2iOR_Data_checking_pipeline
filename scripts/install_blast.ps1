# Installs the blastp.exe binary (NCBI BLAST+) into .venv\Scripts\, so
# scripts/find_protein_mutations.py can shell out to it for identity/mutation
# calculations that match blast.ncbi.nlm.nih.gov exactly. Uses NCBI's
# portable Windows archive (not the win64.exe installer) - just an
# extraction, no installer execution, no admin rights needed. Windows
# counterpart to install_blast.sh.
#
# Re-run this after any `uv sync` that recreates .venv from scratch, since
# .venv isn't version-controlled.

$ErrorActionPreference = "Stop"

$BlastVersion = "2.17.0+"
$RepoRoot = Split-Path -Parent $PSScriptRoot
$VenvScripts = Join-Path $RepoRoot ".venv\Scripts"
$Archive = "ncbi-blast-$BlastVersion-x64-win64.tar.gz"
$Url = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$Archive"

if (-not (Test-Path $VenvScripts)) {
    Write-Error "error: $VenvScripts not found - run 'uv sync' first"
    exit 1
}

$WorkDir = Join-Path $env:TEMP ([System.Guid]::NewGuid())
New-Item -ItemType Directory -Path $WorkDir | Out-Null

try {
    $ArchivePath = Join-Path $WorkDir $Archive

    Write-Host "Downloading $Url ..."
    Invoke-WebRequest -Uri $Url -OutFile $ArchivePath
    # NB: "${Url}.md5" (not "$Url.md5") - a bare dot after an interpolated
    # variable in a double-quoted string is parsed as property access, not
    # literal text; ${} explicitly closes the variable name first.
    Invoke-WebRequest -Uri "${Url}.md5" -OutFile "${ArchivePath}.md5"

    $expectedHash = ((Get-Content "${ArchivePath}.md5") -split '\s+')[0].Trim().ToLower()
    $actualHash = (Get-FileHash -Path $ArchivePath -Algorithm MD5).Hash.ToLower()
    if ($expectedHash -ne $actualHash) {
        throw "MD5 mismatch: expected $expectedHash, got $actualHash"
    }

    # `tar` ships built-in since Windows 10 1803 / Server 2019 and reads
    # .tar.gz directly - same one-liner as the Linux script.
    tar -xzf $ArchivePath -C $WorkDir "ncbi-blast-$BlastVersion/bin/blastp.exe"
    if ($LASTEXITCODE -ne 0) {
        throw "tar extraction failed (exit code $LASTEXITCODE)"
    }

    $BlastpSrc = Join-Path $WorkDir "ncbi-blast-$BlastVersion\bin\blastp.exe"
    $BlastpDst = Join-Path $VenvScripts "blastp.exe"
    Copy-Item -Path $BlastpSrc -Destination $BlastpDst -Force

    Write-Host "Installed: $((& $BlastpDst -version)[0])"
} finally {
    Remove-Item -Recurse -Force $WorkDir
}
