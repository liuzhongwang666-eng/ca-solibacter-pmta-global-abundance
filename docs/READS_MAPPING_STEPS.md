# Reads mapping steps for Candidatus Solibacter abundance

This guide calculates sample-level `Candidatus Solibacter abundance` by mapping
SMAG/SRA reads to Solibacter MAG references. Run these steps on an HPC/server.

## 1. Install tools

```bash
bash scripts/07_install_mapping_env.sh
conda activate solibacter_mapping
```

Required tools: CoverM, BWA/minimap2, samtools, sra-tools, seqkit, pigz.

## 2. Download and extract SMAG MAG archive

```bash
bash scripts/08_download_smag_mag_archive.sh
```

This script queries the Zenodo API, downloads all `mag.tar.gz.*` parts,
combines them into `data/raw/mag.tar.gz`, and extracts MAG FASTA files into
`data/raw/mag/`.

## 3. Extract Solibacter MAG references

```bash
bash scripts/09_extract_solibacter_mags_hpc.sh
```

Input:

- `results/candidate_solibacter_manifest.tsv`
- `data/raw/mag/`

Output:

- `data/processed/solibacter_mags/*.fna`
- `data/processed/solibacter_mag_ids.txt`
- `logs/extract_solibacter_mags.log`

## 4. Download pilot reads

Pilot accessions are listed in `config/pilot_sra_accessions.txt`:

```text
SRR3984960
SRR3985396
SRR11833744
```

Download pilot reads:

```bash
THREADS=8 bash scripts/10_download_sra_reads.sh pilot
```

Output:

- `data/raw/reads/<SRR>_1.fastq.gz`
- `data/raw/reads/<SRR>_2.fastq.gz`

## 5. Run pilot CoverM mapping

```bash
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
```

Outputs:

- `results/coverm_solibacter_pilot/<SRR>.coverm.tsv`
- `results/sample_solibacter_abundance_pilot.tsv`

The sample-level abundance is summarized as:

```text
Candidatus Solibacter abundance per sample
= sum of relative abundance across all Solibacter MAG references
```

## 6. Expand to 28 candidate SRA samples

After the pilot succeeds:

```bash
THREADS=8 bash scripts/10_download_sra_reads.sh all28
THREADS=8 bash scripts/11_run_coverm_solibacter.sh all28
```

Outputs:

- `results/coverm_solibacter_all28/<SRR>.coverm.tsv`
- `results/sample_solibacter_abundance_all28.tsv`

## 7. Quality checks

Check these before using the abundance table in Result 5:

- `relative_abundance` is not all zero.
- `covered_fraction` is not universally low.
- paired reads exist for every SRA sample.
- all expected Solibacter MAG FASTA files were extracted.
- low-quality MAGs can be excluded or used only for sensitivity analysis.

## 8. Notes for manuscript wording

Before reads mapping is complete, use:

```text
Solibacter-related MAG occurrence/distribution
```

After reads mapping is complete, use:

```text
Candidatus Solibacter abundance
```

The current SMAG Supplementary Data 1 and NCBI BioSample records do not provide
reliable soil pH values, so pH correlation still requires additional metadata.
