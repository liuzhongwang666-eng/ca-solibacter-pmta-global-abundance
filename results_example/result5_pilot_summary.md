# Result 5 pilot summary

Date: 2026-06-09

## Completed

- Downloaded SMAG Supplementary Data 1, Supplementary Data 2 and Supplementary Data 2 subgroup from Zenodo.
- Added Ellin6076 16S and `pmtA` seed sequences to `config/targets.fna`.
- Validated the Ellin6076 `pmtA` seed:
  - nucleotide length: 762 bp
  - translated protein length: 253 aa
  - starts with ATG: yes
  - terminal stop codon: yes
  - internal stop codons: 0
- Screened SMAG Supplementary Data 2 with strict Solibacter markers:
  - `g__Solibacter`
  - `s__Solibacter`
  - `f__Solibacteraceae`
  - `Candidatus Solibacter`
  - `Ca. Solibacter`

## Pilot result

- Candidate Solibacter-related SMAG MAG/SGB rows: 30
- Quality groups:
  - high-quality MAGs: 2
  - medium-quality MAGs: 18
  - low/partial-quality MAGs: 10
- Metadata matches in Supplementary Data 1: 30 / 30
- Ecosystems represented:
  - Tundra: 18
  - Wetland: 8
  - Forest: 4
- Continents represented:
  - North America: 28
  - Australia: 1
  - Europe: 1

## Important limitation

SMAG Supplementary Data 1 contains sample ID, project, source, nation,
continent, ecosystem, paired sequencing depth, latitude, longitude,
BioProject and BioSample, but it does not contain an explicit soil pH column.

Therefore:

- current files support a global occurrence/distribution pilot for
  Solibacter-related MAGs;
- true `Candidatus Solibacter abundance` still requires read mapping against
  Solibacter MAG references;
- true `pmtA abundance` still requires pmtA reference construction and read
  mapping;
- soil pH correlation requires pH metadata recovery from BioSample/SRA or the
  original project metadata.

## Next computational step

Download the selected 30 candidate MAG sequences from the SMAG MAG archive and
detect `pmtA` homologs using DIAMOND/BLASTP plus HMMER. In parallel, recover
soil pH for the 30 matched BioSample records, then expand to all SMAG samples if
the metadata recovery is successful.
