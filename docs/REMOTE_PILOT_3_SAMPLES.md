# Remote server pilot: 3-sample Candidatus Solibacter mapping

This document is the short PyCharm remote-server version of the workflow. It is
designed for testing three SRA samples first:

```text
SRR3984960
SRR3985396
SRR11833744
```

The final pilot output is:

```text
results/sample_solibacter_abundance_pilot.tsv
```

## 1. Confirm you are on the remote server

Open PyCharm Terminal and run:

```bash
hostname
pwd
df -h
```

Use the largest data disk, not the system disk. Common paths are:

```bash
/root/autodl-tmp
/mnt/data
/data
```

If your largest disk is `/root/autodl-tmp`, use:

```bash
mkdir -p /root/autodl-tmp/solibacter_mapping
cd /root/autodl-tmp/solibacter_mapping
```

If your largest disk is another path, replace `/root/autodl-tmp` with that path.

## 2. Upload this project folder

Upload the local folder:

```text
E:\博士\砂姜黑土项目投稿\03_补充实验与机制证据\Result5_SMAG_global_abundance
```

to the remote server, for example:

```bash
/root/autodl-tmp/solibacter_mapping/Result5_SMAG_global_abundance
```

Then enter the project:

```bash
cd /root/autodl-tmp/solibacter_mapping/Result5_SMAG_global_abundance
ls
ls scripts
```

You should see:

```text
README.md
READS_MAPPING_STEPS.md
REMOTE_PILOT_3_SAMPLES.md
scripts/
results/
config/
data/
```

## 3. Run a preflight check

```bash
bash scripts/14_remote_pilot_preflight.sh
```

This checks:

- current path
- disk space
- pilot accession list
- required metadata tables
- whether conda/mamba and mapping tools are available

Warnings are acceptable before installation, but missing metadata files or tiny
disk space should be fixed before running the pilot.

## 4. Install the mapping environment

If conda is already active:

```bash
bash scripts/07_install_mapping_env.sh
conda activate solibacter_mapping
```

If `conda: command not found`, first run one of these depending on your server:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
```

or:

```bash
source /root/miniconda3/etc/profile.d/conda.sh
```

Then install:

```bash
bash scripts/07_install_mapping_env.sh
conda activate solibacter_mapping
```

Check:

```bash
coverm --version
samtools --version
prefetch --version
fasterq-dump --version
```

## 5. Download SMAG MAGs and extract Solibacter references

```bash
bash scripts/08_download_smag_mag_archive.sh
bash scripts/09_extract_solibacter_mags_hpc.sh
```

Check:

```bash
du -sh data/raw/mag
ls data/processed/solibacter_mags | head
wc -l data/processed/solibacter_mag_ids.txt
cat logs/extract_solibacter_mags.log | head
```

Expected: around 30 Solibacter MAG FASTA files.

## 6. Download the 3 pilot SRA reads

```bash
cat config/pilot_sra_accessions.txt
THREADS=8 bash scripts/10_download_sra_reads.sh pilot
```

Check:

```bash
ls -lh data/raw/reads/*fastq.gz
du -sh data/raw/reads
```

If disk space is tight, run one sample at a time and delete intermediate reads
after mapping.

## 7. Run CoverM mapping

```bash
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
```

Check:

```bash
ls results/coverm_solibacter_pilot
cat results/sample_solibacter_abundance_pilot.tsv
```

The key output columns are:

```text
Candidatus_Solibacter_relative_abundance_sum
Candidatus_Solibacter_trimmed_mean_sum
Candidatus_Solibacter_covered_fraction_mean
Candidatus_Solibacter_mapped_count_sum
```

## 8. Interpret the pilot result

Good signs:

- `mapped_count_sum` is not all zero.
- `relative_abundance_sum` is not all zero.
- `covered_fraction_mean` is not universally near zero.

If all values are zero, check:

```bash
ls data/processed/solibacter_mags
ls data/raw/reads
cat logs/extract_solibacter_mags.log
```

## 9. Expand after pilot success

Only after the 3-sample pilot works:

```bash
THREADS=8 bash scripts/10_download_sra_reads.sh all28
THREADS=8 bash scripts/11_run_coverm_solibacter.sh all28
```

If disk is smaller than 1 TB, process samples in batches and delete fastq files
after each mapping run.

