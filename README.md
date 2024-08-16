#Metagenome-Assembled Genome of the Obligate Symbiont _Candidatus Stammerula tephritidis_ During Host Plant Shift in the Ecological Divergence of its host, _Tephritis conura_. 

This project involves the quality control, assembly, and analysis of metagenomes, culminating in the identification and classification of Metagenome-Assembled Genomes (MAGs), with a specific focus on identifying the Stammerula MAG.

## Tools and Versions

| Tool             | Version  |
|------------------|----------|
| fastp            | 0.23.4  |
| Kraken2          | 2.1.2   |
| BWA-MEM          | 0.7.17  |
| BEDTools         | 2.31.1  |
| Minimap2         | 2.26    |
| SPAdes           | 3.15.5  |
| metaSPAdes       | 3.15.5  |
| hybridSPAdes     | 3.15.5  |
| Anvio            | 8.0     |
| Prodigal         | 2.6.3   |
| metaWRAP         | 1.3.2   |
| CONCOCT          | 1.0.0   |
| MaxBin2          | 2.2.7   |
| MetaBAT2         | 2.15    |
| CheckM           | 1.0.18  |
| CheckM2          | 1.0.1   |
| BUSCO            | 5.5.0   |
| GTDB-Tk          | 2.4.0   |
| Bakta            | 1.9.3   |
| BLAST+           | 2.15.0  |
| R script dotPlotly | N/A |


## Summary of Procedures

### Quality Control of FASTQ Files
- Raw Illumina reads were trimmed and deduplicated using `fastp`, removing adapter sequences, low-quality reads (Phred score < 20), and reads shorter than 36 bp.
- Host decontamination was performed by aligning reads to the assembled *T. conura* genome using `BWA-MEM` for short reads and `Minimap2` for long reads.
- Unmapped reads were converted back to FASTQ format using the `bamtofastq` module from `BEDTools`.

### Metagenome Co-Assembly
- Co-assemblies were created using `SPAdes` for different population samples, utilizing `metaSPAdes` for short reads and `hybridSPAdes` for mixed read types.

### Metagenome-Assembled Genomes (MAGs)
- MAGs were binned using `metaWRAP`, integrating `CONCOCT`, `MaxBin2`, and `MetaBAT2` to refine bins based on completeness and contamination thresholds.
- The quality of the refined bins was assessed using `CheckM2`, `BUSCO` and `Anvio`.

### Taxonomic Classification and Annotation
- Bins with completeness >50% were taxonomically classified using `GTDB-Tk` and annotated using `Bakta`.
- The Stammerula MAG was identified using a local BLAST search with 16S rRNA sequences against the refined bins.

### Visualization and Validation
- A dot plot was generated using the `R script dotPlotly` to compare the identified MAG with a reference genome.

## References
- (Include the relevant references here as listed in your method section)
