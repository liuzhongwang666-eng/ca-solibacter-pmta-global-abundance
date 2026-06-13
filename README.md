# ca-solibacter-pmta-global-abundance

Workflow for estimating global `Candidatus Solibacter` abundance and functional PmtA marker abundance from soil metagenomes.

This repository supports the Result 5 analysis for a study of `Candidatus Solibacter usitatus` Ellin6076, phosphatidylcholine metabolism, and acidic-soil microbial resources. The workflow was designed around the SMAG global soil metagenome catalogue and sample-level metagenomic reads.

## Scientific Aim

The analysis asks whether the experimentally validated Ellin6076 `pmtA`-PC axis has a broader global context in soil metagenomes.

The two main quantities are:

- `Candidatus Solibacter abundance`: estimated by mapping metagenomic reads to curated `Candidatus Solibacter` MAG references.
- `Functional PmtA marker abundance`: estimated by translated read searches against curated protein markers defined by the functional Ellin6076 PmtA protein `gene101949`.

The workflow intentionally avoids treating exact nucleotide `pmtA` mapping as the primary abundance metric, because `pmtA` homologs can be divergent at the DNA level while retaining conserved protein features.

## Conceptual Workflow

```text
SMAG metadata and MAG catalogue
→ screen candidate Candidatus Solibacter MAGs
→ extract Solibacter MAG references
→ map metagenomic reads to Solibacter MAGs with CoverM
→ calculate Candidatus Solibacter relative abundance

Ellin6076 gene101949 functional PmtA protein
→ identify high-confidence PmtA-like protein markers in Solibacter MAGs
→ translated metagenomic read search with DIAMOND blastx
→ calculate functional PmtA marker abundance
→ combine Solibacter abundance and PmtA marker abundance
```

## Repository Structure

```text
config/          Small configuration files and seed sequences
scripts/         Reproducible analysis scripts
docs/            Detailed server, mapping, and PmtA marker notes
results_example/ Small manifest/example files; large outputs are excluded
```

Large files such as SRA reads, FASTQ files, MAG FASTA files, DIAMOND databases, BAM/SAM files, and large intermediate outputs are excluded by `.gitignore`.

## Environment

On the remote server, activate the mapping environment first:

```bash
source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
```

The successful server run used NCBI SRA Toolkit 3.4.1 installed outside conda. Add it before downloading SRA reads:

```bash
export PATH=/root/autodl-tmp/tools/sratoolkit.3.4.1-ubuntu64/bin:$PATH
hash -r
prefetch --version
```

Required tools:

```text
coverm
minimap2
samtools
diamond
prodigal
prefetch
fasterq-dump
pigz or gzip
python3
```

## Main Scripts

### Metadata and Target Preparation

```text
00_download_smag_metadata.ps1
01_prepare_targets.py
02_screen_smag_taxonomy.py
02b_make_solibacter_manifest.py
02c_join_solibacter_metadata.py
02d_make_sra_accession_list.py
02e_fetch_biosample_ph.py
```

These scripts download SMAG metadata, prepare seed sequences, screen candidate `Candidatus Solibacter` MAGs, and build sample manifests.

### Candidatus Solibacter Abundance

```text
08_download_smag_mag_archive.sh
09_extract_solibacter_mags_hpc.sh
11_run_coverm_solibacter.sh
12_summarize_solibacter_abundance.py
```

The Solibacter abundance workflow maps sample reads to extracted Solibacter MAG references using CoverM:

```bash
THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot
```

The key metrics are:

- `relative_abundance_percent`: sum of Solibacter MAG relative abundance after excluding `unmapped`.
- `read_count`: reads mapped to Solibacter MAG references.
- `trimmed_mean`: robust coverage-depth estimate.
- `covered_fraction`: fraction of each MAG covered by mapped reads.

### Functional PmtA Marker Abundance

```text
15_build_functional_pmta_markers.sh
16_filter_pmta_protein_hits.py
17_run_pmta_marker_blastx_pilot.sh
18_summarize_pmta_marker_abundance.py
```

The functional PmtA marker workflow uses `gene101949` as the functional Ellin6076 PmtA seed. `gene102575` is treated as a PmtA-like paralog/control, not as the functional seed.

Candidate PmtA markers are filtered by:

- homology to functional `gene101949`;
- SAM-binding motif `E/DXGXGXG`;
- N-terminal membrane-binding helix proxy;
- high-confidence protein-level classification.

Build functional PmtA markers:

```bash
THREADS=8 bash scripts/15_build_functional_pmta_markers.sh
```

Run translated read searches for pilot samples:

```bash
THREADS=8 bash scripts/17_run_pmta_marker_blastx_pilot.sh pilot
```

The key metrics are:

- `pmta_marker_hit_reads`;
- `pmta_marker_reads_per_million`;
- `pmta_marker_rpkm_like`;
- `n_markers_with_hits`;
- `top_marker`.

Use cautious wording:

```text
functional PmtA marker abundance
PmtA-like functional potential
translated metagenomic read searches against curated functional PmtA protein markers
```

Avoid:

```text
exact nucleotide pmtA abundance
```

## Batch Workflow for Candidate Samples

The 3 pilot samples were:

```text
SRR3984960
SRR3985396
SRR11833744
```

The remaining candidate samples are processed in batches of 5.

Generate batch files:

```bash
python3 scripts/19_make_sra_batches.py
```

Run each batch:

```bash
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_01.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_02.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_03.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_04.txt
THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_05.txt
```

Summarize all completed outputs:

```bash
python3 scripts/21_summarize_all28_solibacter_pmta.py
```

The combined output is:

```text
results/all28_solibacter_pmta_combined_summary.csv
```

## Output Interpretation

The combined summary table includes:

- `solibacter_relative_abundance_percent`: total Solibacter MAG relative abundance for a sample.
- `solibacter_mapped_reads`: reads mapped to Solibacter MAG references.
- `solibacter_trimmed_mean`: summed robust coverage-depth estimate.
- `solibacter_mean_covered_fraction`: mean coverage fraction across Solibacter MAG references.
- `pmta_marker_reads_per_million`: normalized translated-read signal to functional PmtA markers.
- `pmta_marker_rpkm_like`: length-normalized PmtA marker abundance proxy.
- `pmta_RPM_per_Solibacter_percent`: PmtA marker RPM normalized by Solibacter relative abundance.

The ratio `pmta_RPM_per_Solibacter_percent` should be interpreted cautiously. A high value can indicate PmtA-like functional potential in samples where Solibacter MAG abundance is low, but it may also reflect PmtA-like markers from related non-Solibacter microbes.

## Documentation

Additional documentation is available in `docs/`:

- `READS_MAPPING_STEPS.md`
- `REMOTE_PILOT_3_SAMPLES.md`
- `REMOTE_SERVER_AND_GLOBAL_ABUNDANCE_GUIDE.md`
- `PMTA_FUNCTIONAL_MARKER_WORKFLOW.md`

## Data Availability

This repository does not include raw reads, large MAG files, or large intermediate mapping outputs. These files should be regenerated from SMAG/NCBI SRA sources using the scripts in this repository.

Small manifest files are provided in `results_example/` to document the candidate MAGs and samples used by the workflow.
