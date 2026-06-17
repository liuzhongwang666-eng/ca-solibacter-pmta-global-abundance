# Global acidic-soil balanced pilot v2

## Purpose

The first global balanced pilot table selected five consecutive SRR runs from each target BioProject. Although the SRR and BioSample identifiers were not exact duplicates, many records represented the same project and nearly the same location. This inflated project-level replication and weakened the intended global comparison.

Version 2 applies BioProject-level and BioSample-level deduplication before mapping:

- each SRR is unique;
- each BioSample is unique;
- each BioProject contributes no more than two selected SRR runs by default;
- all selected records must be soil-derived;
- all selected records must have confirmed latitude and longitude;
- records outside the target continent are excluded.

## Outputs

```text
results/global_acidic_soil_balanced_candidates_v2.tsv
config/global_acidic_soil_balanced_pilot_sra_accessions_v2.txt
results/global_acidic_soil_balanced_excluded_v2.tsv
```

The strict v2 filter currently yields 33 high-quality samples rather than forcing the target to 40. This is intentional: adding more samples under the current metadata would require repeated BioProjects, unclear soil source metadata, missing coordinates, or continent/location mismatches.

## Current v2 summary

```text
Selected samples: 33
Unique SRR: 33
Unique BioSample: 33
BioProjects with more than two selected runs: 0
Samples with confirmed latitude/longitude: 33
```

| Continent | Soil group | n |
|---|---|---:|
| Asia | China red/acid soil | 5 |
| Asia | Acidic paddy soil | 5 |
| Europe | Acidic forest soil | 2 |
| North America | Peat/permafrost/tundra soil | 5 |
| South America | Amazon/tropical acidic soil candidate | 5 |
| Africa | Acidic/Ferralsol/Acrisol soil candidate | 4 |
| Oceania | Acid sulfate soil | 3 |
| Antarctica | Polar soil | 4 |

## Rebuild the v2 table

Run from the project root:

```bash
python3 scripts/27_build_global_balanced_pilot_v2.py --retmax 200 --sleep 0.6
```

The script uses NCBI SRA queries, fetches linked BioSample XML metadata, extracts source evidence and coordinates, and then applies continent-aware soil-only filtering and BioProject-level deduplication.

## Run mapping for v2 samples

On the remote server:

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance

source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
export PATH=/root/autodl-tmp/tools/sratoolkit.3.4.1-ubuntu64/bin:$PATH
hash -r

THREADS=8 bash scripts/25_process_acidic_soil_extension_batch.sh config/global_acidic_soil_balanced_pilot_sra_accessions_v2.txt
```

Summarize:

```bash
python3 scripts/26_summarize_acidic_soil_extension.py \
  --samples config/global_acidic_soil_balanced_pilot_sra_accessions_v2.txt \
  --metadata results/global_acidic_soil_balanced_candidates_v2.tsv \
  --out results/global_acidic_soil_balanced_solibacter_pmta_summary_v2.csv
```

## Interpretation

The v2 set should be treated as the main globally balanced pilot because it prioritizes metadata reliability over sample count. If the final analysis requires 35-42 samples, the next step should be targeted manual curation for the underrepresented continents, especially Europe, Oceania, Africa, and Antarctica, rather than reusing additional consecutive SRR runs from the same BioProject.
