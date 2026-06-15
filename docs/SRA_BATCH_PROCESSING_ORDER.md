# SRA Batch Processing Order for Candidatus Solibacter and PmtA Marker Abundance

## 1. Purpose

This document records the download and processing order for the global abundance analysis of `Candidatus Solibacter` and functional `PmtA` markers.

The workflow first processes three pilot samples to confirm that SRA download, `CoverM` genome mapping, and PmtA marker `blastx` screening work correctly. After the pilot run, the remaining 25 candidate SRA samples are processed in five batches, with five samples per batch.

## 2. Remote Working Directory

All commands below assume the analysis is run on the remote server under:

```bash
/root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
```

Before running any batch, enter the project directory and activate the mapping environment:

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
export PATH=/root/autodl-tmp/tools/sratoolkit.3.4.1-ubuntu64/bin:$PATH
hash -r
```

The `export PATH` command is important because the conda environment may contain an older `prefetch` version. The workflow should use SRA Toolkit 3.4.1.

## 3. Pilot Samples

The pilot run contains three samples:

```text
SRR3984960
SRR3985396
SRR11833744
```

Pilot commands:

```bash
THREADS=8 bash scripts/10_download_sra_reads.sh pilot
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
THREADS=8 bash scripts/17_run_pmta_marker_blastx_pilot.sh
python3 scripts/18_summarize_pmta_marker_abundance.py
```

The pilot results are used to verify that:

- SRA reads can be downloaded and converted to paired FASTQ files.
- Reads can be mapped to the `Candidatus Solibacter` MAG reference set.
- PmtA marker hits can be detected using translated-read alignment.
- The sample-level abundance summary can be generated.

## 4. Batch Samples

After the pilot samples are completed, the remaining 25 samples are processed in batches of five.

### batch_01

```text
SRR10489815
SRR12756247
SRR3984963
SRR3985190
SRR3985390
```

### batch_02

```text
SRR3985395
SRR3985424
SRR3989121
SRR3989257
SRR5262248
```

### batch_03

```text
SRR5450438
SRR5450747
SRR5664280
SRR5665647
SRR5665811
```

### batch_04

```text
SRR5674476
SRR7097862
SRR7097863
SRR8879130
SRR8879131
```

### batch_05

```text
SRR8879134
SRR8881005
SRR8881006
SRR8881008
SRR8881010
```

## 5. Batch Processing Commands

Generate the batch files:

```bash
python3 scripts/19_make_sra_batches.py
```

Run each batch sequentially:

```bash
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_01.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_02.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_03.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_04.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_05.txt
```

Each batch command performs the following steps for every SRR accession:

1. Download the SRA file with `prefetch`.
2. Convert SRA to paired FASTQ files using `fasterq-dump`.
3. Compress FASTQ files with `pigz` or `gzip`.
4. Map reads to the `Candidatus Solibacter` MAG reference set with `CoverM`.
5. Search reads against the functional PmtA marker database with `DIAMOND blastx`.

## 6. Summary Command

After any batch finishes, update the combined abundance summary:

```bash
python3 scripts/21_summarize_all28_solibacter_pmta.py
```

This script combines:

- `Candidatus Solibacter` relative abundance from `CoverM`.
- `Candidatus Solibacter` mapped reads.
- PmtA marker reads per million.
- PmtA marker RPKM-like values.
- Metadata such as ecosystem and location.

## 7. Important Note for SRR5262248

`SRR5262248` is approximately 60 GB. The default `prefetch` download limit is often 20 GB, so this sample will be skipped unless the maximum size is increased.

The batch script should contain:

```bash
prefetch "${srr}" --output-directory "${READ_DIR}" --max-size 200G
```

If the script still contains the shorter command below, update it before rerunning `batch_02`:

```bash
prefetch "${srr}" --output-directory "${READ_DIR}"
```

## 8. Checking Progress

Check generated `CoverM` files:

```bash
ls -lh results/coverm_solibacter_all28
```

Check generated PmtA marker files:

```bash
ls -lh results/pmta_marker_blastx_all28
```

Check the combined summary:

```bash
head results/all28_solibacter_pmta_summary.csv
```

## 9. Interpretation Boundary

The `Candidatus Solibacter` abundance is estimated by reads mapping to the selected Solibacter MAG reference set.

The PmtA marker abundance is estimated by translated-read alignment against curated functional PmtA-like protein markers. It should be described as `functional PmtA marker abundance`, not as exact nucleotide-level `pmtA gene abundance`.
