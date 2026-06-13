# Functional PmtA Marker Workflow

This workflow replaces direct nucleotide mapping to full-length `pmtA` genes with a protein-marker strategy.

## Rationale

Direct mapping to all `pmtA` nucleotide sequences can miss divergent homologs. The functional Ellin6076 PmtA signal should therefore be defined at the protein level, using:

- `gene101949` as the functional Ellin6076 PmtA seed.
- `gene102575` as a PmtA-like paralog/control, not as the functional seed.
- Conserved protein features:
  - N-terminal amphipathic/membrane-binding helix proxy.
  - SAM-binding motif: `E/DXGXGXG`.
  - Homology to the functional `gene101949` protein.

The output should be described as `PmtA functional marker abundance` or `PmtA functional potential`, not exact `pmtA gene abundance`.

## Server Commands

Run from the project root:

```bash
cd /root/autodl-tmp/ca_pmtA/Result5_SMAG_global_abundance
source /root/miniconda3/etc/profile.d/conda.sh
conda activate solibacter_mapping
```

If `diamond` or `prodigal` are missing:

```bash
conda install -c bioconda -c conda-forge diamond prodigal -y
```

Build functional PmtA protein markers:

```bash
THREADS=8 bash scripts/15_build_functional_pmta_markers.sh
```

This creates:

```text
data/processed/pmta_markers/functional_pmta_marker_candidates.faa
results/pmta_functional_markers/functional_pmta_marker_candidates.tsv
```

Run translated read searches for the 3 pilot samples:

```bash
THREADS=8 bash scripts/17_run_pmta_marker_blastx_pilot.sh pilot
```

This creates:

```text
results/pmta_marker_blastx_pilot/*.pmta_marker.blastx.tsv
results/sample_pmta_marker_abundance_pilot.csv
```

## Interpretation

Important output columns:

- `pmta_marker_hit_reads`: reads with filtered translated hits to functional PmtA markers.
- `pmta_marker_reads_per_million`: normalized PmtA marker signal per million searched reads.
- `pmta_marker_rpkm_like`: length-normalized protein-marker abundance proxy.
- `n_markers_with_hits`: number of functional PmtA markers detected in the sample.

For manuscript wording, use cautious language:

```text
translated metagenomic read searches against curated functional PmtA protein markers
```

Avoid:

```text
exact pmtA gene abundance from nucleotide mapping
```

## Notes

The N-terminal membrane-binding feature is currently implemented as a heuristic proxy based on N-terminal hydrophobic and basic residue composition. It should be treated as a screening criterion and, where possible, supported by alignment, phylogeny, or structural analysis.
