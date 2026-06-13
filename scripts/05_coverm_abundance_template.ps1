$ErrorActionPreference = "Stop"

# This is a template for the final abundance calculation.
# It assumes CoverM, BWA/minimap2, and samtools are available in WSL/Linux/HPC.
# On Windows, run these commands inside WSL or on a cluster.

Write-Host "Expected inputs:"
Write-Host "  data/processed/solibacter_mags/*.fa*"
Write-Host "  data/processed/pmtA_reference_genes.fna"
Write-Host "  data/raw/reads/<sample>_1.fastq.gz and <sample>_2.fastq.gz"
Write-Host ""
Write-Host "Recommended CoverM commands:"
Write-Host @'
# Solibacter MAG abundance
coverm genome \
  --genome-fasta-directory data/processed/solibacter_mags \
  --coupled data/raw/reads/*_1.fastq.gz data/raw/reads/*_2.fastq.gz \
  --methods relative_abundance trimmed_mean covered_fraction \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 90 \
  --output-file results/solibacter_coverm_abundance.tsv

# pmtA gene abundance
coverm contig \
  --reference data/processed/pmtA_reference_genes.fna \
  --coupled data/raw/reads/*_1.fastq.gz data/raw/reads/*_2.fastq.gz \
  --methods rpkm trimmed_mean count \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 90 \
  --output-file results/pmta_coverm_abundance.tsv
'@
