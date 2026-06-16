# Acidic-soil raw-reads extension workflow

## 1. Rationale

NCBI contains several low-quality `Candidatus Solibacter` MAGs. These MAGs should not be added directly to the mapping reference set because fragmented or incomplete genomes can bias read-mapping abundance estimates.

The corrected strategy is:

```text
Low-quality NCBI Candidatus Solibacter MAGs
-> use only as source-project clues
-> retrieve raw metagenomic SRA reads from soil-confirmed source projects
-> map reads to the existing high-quality Solibacter MAG reference set
-> estimate Candidatus Solibacter abundance and functional PmtA marker abundance
```

## 2. Main rule

The primary analysis is strictly soil-only.

Excluded environments include:

```text
sediment
sludge
activated sludge
groundwater
lake water column
biofilter
bioreactor
rock
hot spring fluid
```

## 3. Reference set

The mapping reference remains the existing downloaded `Candidatus Solibacter` MAG reference directory:

```bash
data/processed/solibacter_mags
```

NCBI MAGs are not added to this directory by default.

The reference manifest is:

```text
config/solibacter_high_quality_reference_manifest.tsv
```

If the remote server contains 29 `.fna` references rather than the 30 local candidates listed in the manifest, the remote `.fna` directory is authoritative.

## 4. NCBI source-project discovery

Soil-confirmed NCBI `Candidatus Solibacter` assemblies are used only to identify source BioProjects and BioSamples.

Configuration file:

```text
config/ncbi_solibacter_soil_source_projects.tsv
```

Generate source-project SRA runs:

```bash
python3 scripts/23_fetch_ncbi_solibacter_source_runs.py
```

Outputs:

```text
results/ncbi_solibacter_soil_project_sra_runs.tsv
results/ncbi_solibacter_soil_project_sra_runs_excluded.tsv
```

## 5. Acidic-soil SRA curation

Curate soil-only acidic SRA candidates:

```bash
python3 scripts/24_curate_acidic_soil_sra_candidates.py
```

Outputs:

```text
results/acidic_soil_sra_extension_candidates.tsv
results/acidic_soil_sra_extension_excluded.tsv
config/acidic_soil_pilot_sra_accessions.txt
```

Current pilot accession list:

```text
SRR10489828
SRR10489818
SRR10489834
SRR10489835
SRR13451248
SRR29932441
SRR29932559
SRR13312976
SRR29932439
SRR29932435
SRR26386574
SRR26386599
SRR26386608
SRR26386572
SRR26386581
```

These represent three soil project groups:

- `PRJNA588342`: permafrost/tundra soil, Chersky, Russia.
- `PRJNA682830`: high-severity burned soil, Wyoming, USA.
- `PRJNA746701`: Antarctic soil, Byers Peninsula, Livingston Island.

This 15-sample set is retained as an NCBI Solibacter source-project pilot, but it is not globally balanced.

## 6. Global balanced acidic-soil pilot

For the first global pilot, use the 40-sample balanced set:

```text
config/global_acidic_soil_balanced_pilot_sra_accessions.txt
results/global_acidic_soil_balanced_candidates.tsv
```

This set covers Asia, Europe, North America, South America, Africa, Oceania, and Antarctica, with five samples per group:

| Group | Continent | Soil type | BioProject | n |
|---|---|---|---|---:|
| Asia red/acid soil | Asia | China red/acid soil | `PRJNA1478117` | 5 |
| Asia acidic paddy soil | Asia | acidic paddy soil | `PRJNA1337480` | 5 |
| Europe acidic forest soil | Europe | acidic forest soil | `PRJNA1431341` | 5 |
| North America peat/permafrost soil | North America | peat/permafrost/tundra soil | `PRJNA1474608` | 5 |
| South America Amazon/tropical soil candidate | South America | Amazon/tropical soil, pH to confirm | `PRJNA1470217` | 5 |
| Africa acidic/Ferralsol/Acrisol candidate | Africa | acidic soil/Ferralsol/Acrisol candidate, pH to confirm | `PRJNA703480` | 5 |
| Oceania acid sulfate soil | Oceania | Australian acid sulfate soil | `PRJNA1016489` | 5 |
| Antarctica polar soil | Antarctica | Antarctic polar soil | `PRJNA746701` | 5 |

The 40-sample balanced set is recommended for the first Result 5 global distribution run.

## 7. Remote processing commands

Run on the remote server:

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance

source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
export PATH=/root/autodl-tmp/tools/sratoolkit.3.4.1-ubuntu64/bin:$PATH
hash -r
```

Process the NCBI Solibacter source-project pilot set:

```bash
THREADS=8 bash scripts/25_process_acidic_soil_extension_batch.sh config/acidic_soil_pilot_sra_accessions.txt
```

Process the global balanced acidic-soil pilot set:

```bash
THREADS=8 bash scripts/25_process_acidic_soil_extension_batch.sh config/global_acidic_soil_balanced_pilot_sra_accessions.txt
```

Summarize the global balanced pilot results:

```bash
python3 scripts/26_summarize_acidic_soil_extension.py \
  --samples config/global_acidic_soil_balanced_pilot_sra_accessions.txt \
  --metadata results/global_acidic_soil_balanced_candidates.tsv \
  --out results/global_acidic_soil_balanced_solibacter_pmta_summary.csv
```

Global balanced output:

```text
results/global_acidic_soil_balanced_solibacter_pmta_summary.csv
```

## 8. T6 paddy-soil MAG clue

The four low-quality T6 paddy-soil MAGs from `PRJNA1004092` remain important China paddy-soil Solibacter clues:

```text
GCA_032665965.1
GCA_032666425.1
GCA_032663185.1
GCA_032662945.1
```

However, current NCBI SRA/E-utilities searches did not return raw SRA runs for `PRJNA1004092` or the corresponding BioSamples. Therefore, these MAGs are not used for reads mapping and are not added to the reference set. The Asia acidic paddy-soil group in the balanced pilot substitutes downloadable paddy-soil metagenome samples for this missing raw-read source.

## 9. Detection thresholds

Continuous abundance values should always be retained.

For map-style presence/absence visualization:

```text
Solibacter detected:
mapped reads >= 100 and mean covered fraction >= 0.01

High-confidence detected:
mapped reads >= 1000 and mean covered fraction >= 0.05
```

PmtA abundance should be described as:

```text
functional PmtA marker abundance
translated-read signal against curated PmtA protein markers
```

Avoid:

```text
exact nucleotide pmtA gene abundance
```

## 10. Leave-one-sample-out mapping

If a future acidic-soil extension sample is also the source of a MAG in the Solibacter reference set, add it to:

```text
config/solibacter_reference_exclusions_by_sample.tsv
```

Format:

```text
sample_id	excluded_mag_ids	reason
SRRxxxxxxx	SRRxxxxxxxbin.1,some_other_mag	self-derived reference removal
```

The extension batch script will create a temporary per-sample reference directory with those MAGs excluded.
