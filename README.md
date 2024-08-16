# Metagenome-Assembled Genome of the Obligate Symbiont _Candidatus Stammerula tephritidis_ During Host Plant Shift in the Ecological Divergence of its host, _Tephritis conura_. 

A metagenomic analysis is exhaustive and complicated as it utilizes a plathera of tools. Through this documentation, the steps taken to produce the _Candidatus Stammerula tephritis_ metagenome-assembled genome (MAG) is described.

There are multiple ways of getting from raw sequences to a MAG. During this research project, a multitude of tools, wrappers, programs and workflow managers were tested.
Thus, the repository is filled with various scripts, and as this research project has ended, my research thesis will start and will be a continuation of this work.
Hence, scripts are abundant and important, as i might require them in the future.

To understand my decision making, fallbacks and wins, please look at the **issues** tab on this repository. It has been used as a log book and might come into hand. As well as my project report.

Finally, before starting this journey, whole-metagenome analysis is quite a hefty mission.
Basic bioinformatical foundamentals must already be understood before getting lost in the jungle of metagenomics.

This guides assumes you have these fundamentals, and therefore, I will speed through the many steps.
I would say there are three major steps

1. `Quality control and Data wrangling`. Many ways to wrangle and structure. Takes time to find optimal.
2. `Assemble metagenome and populate database with contig and alignment stats`. Anvio has so many tools...
3. `Bin curation`: Binning, then populating bins with various statistics, and then finding which one that is best.
   - it is a lot of back and forth, and parsing.
  
Hence, _head my warning_, some scripts are hard coded, some scripts you must give an argument/arguments, some work with `sbatch` (SLURM), some require interactive sessions, some require manual editing script parameters.

From raw reads to assembled metagenome is quite straight forward. 

Setting up Anvio contig and profile databases, and binning this new formatted metagenome is a breeze and well structured.

It is at the moment the bin refinement into a MAG that is clunky and requires some hours of work. It requires a lot of analyzing different results. 

Nonetheless, here we go!


## Tools and Versions
These are the tools and their versions, that was used to perfome this scientific study.
Noteworthy, metaWRAP v.1.3.2 forexample includes the binners `CONCOCT`. `MaxBin2`. and `MetaBAT2`. Depending on the version of `metaWRAP`, the version of the binners might differ.
One can approach this by doing all the steps by themselves. but the wrapper programs such as `Anvio`, `Bakta` and `metaWRAP` makes ones life easier.


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


## Short summary of Procedures

### Quality Control of FASTQ Files
- Raw Illumina reads were trimmed and deduplicated using `fastp`, removing adapter sequences, low-quality reads (Phred score < 20), and reads shorter than 36 bp.
- Host decontamination was performed by aligning reads to the assembled *T. conura* genome using `BWA-MEM` for short reads and `Minimap2` for long reads.
- Unmapped reads were converted back to FASTQ format using the `bamtofastq` module from `BEDTools`.
- Reads where then merged accordingly to "all", "CHST", and "COGE".

### Metagenome Co-Assembly
- Co-assemblies were created using `SPAdes` for different population samples, utilizing `metaSPAdes` for short reads and `hybridSPAdes` for mixed read types.
- k-mer sizes of 21, 33 and 55 were used. (Read https://doi.org/10.1002/cpbi.102, to better understand the importance of these. This set of kmers sizes is standard and is shown to perform well.)
- These co-assembled metagenomes are the starting point of the Anvio workflow.
- These were then, reformatted to get standardized contig names and filtered to only keep contigs > 2500 kb.

### Anvio
- The metagenome is reformatted.
- Contig database is created.
- Reads are aligned to reformatted metagenome.
- Profile database is created.
- Profiles are merged into a merged profile.

### Metagenome-Assembled Genomes (MAGs)
- MAGs were binned using `metaWRAP`, integrating `CONCOCT`, `MaxBin2`, and `MetaBAT2` to refine bins based on completeness and contamination thresholds.
- The quality of the refined bins was assessed using `CheckM2`, `BUSCO` and `Anvio`.

### Taxonomic Classification and Annotation
- Bins with completeness >50% were taxonomically classified using `GTDB-Tk` and annotated using `Bakta`.
- The Stammerula MAG was identified using a local BLAST search with 16S rRNA sequences against the refined bins.


# Trimming
Utilizing metadata files in `doc/` and scripts in `scripts` we will blast through the quality control.
The sample data is split across populations, across samples, and across lanes.
Therefore, run the script `fastp-trim-per-lane.sh`.
It utilizes the `.lst` files that holds paths to corresponding sample lanes.
This script will put `fastp` quality reports in `01-QC` and the trimmed lane `.fastq` files in `02-TRIM/`
These were the settings used for `fastp`
```
    fastp \
        --in1 $R1_IN --in2 $R2_IN \
        --out1 $R1_OUT --out2 $R2_OUT \
        --html "${FASTP_DIR}/${R2_TRIMD}-fastp.html" \
        --json "${FASTP_DIR}/${R2_TRIMD}-fastp.json" \
        --average_qual 20 \
        --length_required 36 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --dedup \
        --thread 16 \
        --verbose 
```
# Decontamination
Run the `bwa` and `minimap` scripts found here:
```
scripts/decontamination/
├── bwa-index-genome.sh
├── bwa-mem-decon.sh
├── kneaddata_run.sh
├── minimap2-decon.sh
└── minimap2-index.sh
```
They will create a corresponding `bwa` and `minimap` index of the _T. conura_ reference genome and place the clean reads in `04-CLEAN-FASTQ`.
The indexes and alignments will be stored in `03-HOST-ALIGNMENT` for all population, as seen below.

```
03-HOST-ALIGNMENT-BAM/
CHES  CHFI  CHSC  CHSK  CHST  COES  COGE  COLI  COSK  CPSC  hifi-pacbio

03-HOST-ALIGNMENT-BAM/
├── CHES
│   ├── P12002_144_S15_L001
│   │   ├── P12002_144_S15_L001.bam
│   │   ├── P12002_144_S15_L001-sorted-unmapped-pairs.bam
│   │   ├── P12002_144_S15_L001-unmapped-pairs.bam
│   │   └── P12002_144_S15_L001-unmapped-reads.bam
│   ├── P12002_144_S60_L003
│   │   ├── P12002_144_S60_L003.bam
│   │   ├── P12002_144_S60_L003-sorted-unmapped-pairs.bam
│   │   ├── P12002_144_S60_L003-unmapped-pairs.bam
│   │   └── P12002_144_S60_L003-unmapped-reads.bam
│   ├── P12002_144_S60_L006
...
```
The clean reads will be stored as following.
```
04-CLEAN-FASTQ/
├── CHES
│   ├── P12002_144_S15_L001_R1-clean.fastq.gz
│   ├── P12002_144_S15_L001_R2-clean.fastq.gz
│   ├── P12002_144_S60_L003_R1-clean.fastq.gz
│   ├── P12002_144_S60_L003_R2-clean.fastq.gz
│   ├── P12002_144_S60_L006_R1-clean.fastq.gz
...
│   └── P12002_153_S63_L006_R2-clean.fastq.gz
├── CHFI
│   ├── P12002_122_S20_L002_R1-clean.fastq.gz
│   ├── P12002_122_S20_L002_R2-clean.fastq.gz
...
```

# Merge clean reads.
Merging is as simple as concatenating, however, depending on how you want to merge they will differ. I merged all reads, population wise, and host plant wise.
Utilize these, and you will have merged `.fastq` files as following.
```
scripts/merge-clean-reads/
├── merge-all-reads.sh
├── merge-clean-reads.sh
└── merge-hostplant-reads.sh

05-CLEAN-MERGED/
├── all_R1-clean.fq.gz
├── all_R2-clean.fq.gz
├── CHES_R1-clean.fq.gz
├── CHES_R2-clean.fq.gz
├── CHFI_R1-clean.fq.gz
├── CHFI_R2-clean.fq.gz
├── CH_R1-clean.fq.gz
├── CH_R2-clean.fq.gz
├── CHSC_R1-clean.fq.gz
├── CHSC_R2-clean.fq.gz
├── CHSK_R1-clean.fq.gz
├── CHSK_R2-clean.fq.gz
├── CHST_R1-clean.fq.gz
├── CHST_R2-clean.fq.gz
├── COES_R1-clean.fq.gz
├── COES_R2-clean.fq.gz
├── COGE_R1-clean.fq.gz
├── COGE_R2-clean.fq.gz
├── COLI_R1-clean.fq.gz
├── COLI_R2-clean.fq.gz
├── CO_R1-clean.fq.gz
├── CO_R2-clean.fq.gz
├── COSK_R1-clean.fq.gz
├── COSK_R2-clean.fq.gz
├── CPSC_R1-clean.fq.gz
└── CPSC_R2-clean.fq.gz
```
# Assemble
Wonderful, done with data wrangling.
Here, you will need to utilize two scripts: `hybridspades-assembly.sh`, and `metaspades-assembly.sh` depending on the data.
`CHST` and `all` utilize long read data in a *hybrid* fashion, thus, these metagenomes were assembled with `hybridspades-assembly.sh`.
If the population only consists of short reads, then use `metaspades-assembly.sh`.

```
scripts/assembly
├── hybridspades-assembly.sh
├── megahit-metagenome-assembly.sh
├── metaspades-assembly.sh
├── spades-lr-metagenome-assembly.sh
├── spades-sr-CHST-metagenome-assembly.sh
└── spades-sr-metagenome-assembly.sh
```
Noteworthy, you will have to specify which population you want to assemple as this script requires an argument.
As our clean and merged reads have standardized names with the prefix of `$POP`. 
Example, you want to assemble the COGE metagenome, then run this: 

```sbatch /scripts/assembly/metaspades-assembly.sh COGE```

The process is straight as an arrow!
The script will store the metagenome as following `06-ASSEMBLY/$POP`

```
# Set variables
POP=$1
SR_DIR="05-CLEAN-MERGED"

# Create directory for trimmed reads if not existing
outdir="/crex/proj/snic2020-6-222/Projects/Tconura/working/Andre/CONURA_WGS/06-ASSEMBLY/${POP}"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# Assembling the metagenome
spades.py \
    --meta \
    --only-assembler \
    -k auto \
    --threads 20 \
    --memory 1000 \
    -1 ${SR_DIR}/${POP}_R1-clean.fq.gz \
    -2 ${SR_DIR}/${POP}_R2-clean.fq.gz \
    -o $outdir
```
# Now the fun starts.
If you need help installing `Anvio`, specifically on a HPC then check my issue here: [Installing Anvio-8, Pavlins approach](https://github.com/ndreey/CONURA_WGS/issues/23#issuecomment-2051736092)

If you want to understand the Anvio steps better, then check [Anvi'o User Tutorial for Metagenomic Workflow](https://merenlab.org/2016/06/22/anvio-tutorial-v2/)

To not clutter the readme, take a peak here at my issue: [REVAMP: Structured standard for MAG curation for all approaches](https://github.com/ndreey/CONURA_WGS/issues/48)

It explains the dynamic scripts i created to make the `Anvio` and `metaWRAP` parts easier and well structured.

# Bin bin bins, lets bin the bins, until we get the optimal bins.
Following previous step you should have this:
```
 01-ALIGNMENT  02-CONTIG-DB  03-PROFILES  04-MERGED-PROFILES  05-metaWRAP  
```
With a bunch of bin results from the different binners in `05-metaWRAP/$POP/`.
In the script folder `mag-curation`, there are a bunch of scripts.
```
scripts/mag-curation/
├── checkm2-check.sh
├── create-contigs2bin.sh
├── gtdb-taxo.sh
├── make-contig-db.sh
├── make-profile-db-FIX.sh
├── make-profile-db.sh
├── metawrap-bin.sh
├── metawrap-refine.sh
├── reformat-fasta-fix.sh
├── reformat-fasta.sh
└── summarize-gtdbtk.sh
```
We want to use `metawrap-refine.sh`.
For my analysis i used the thresholds `COMP=50 CONT=10`. 
I got to this thresholds by first testing a range of different thresholds and i found this to be the optimal one.
Check out the awesome resources on `metaWRAP`'s: [git/metaWRAP](https://github.com/bxlab/metaWRAP),  

And dont forget to read the paper, specifically regarding the `bin_refine` module.
(https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1)

With all the bins, and a cluster of aggregated bins from the `bin_refinement` model from `metaWRAP`. It is time to gather summary statistics of each to determine which bin meets the criteria for a quality MAG.

## CheckM2
`metaWRAP` uses `CheckM1`, and older version, and believe even stopped being updated.
`CheckM2` is great as the _gradient boosted_ algorithm has been shown to perform better assesing novel microbes.
Use the `checkm2-check.sh` script in the `mag-curation` folder.
It will assess all the bins binned by metaWRAP for given population argument and summarize the information into a `.tsv` file.
The output from CheckM2 is also stored in the corresponding population folder.
```
ls 06-CHECKM2/
all  CHST  COGE  results

cat 06-CHECKM2/results/CHST-bin-report.tsv 
Binner  Name    Completeness_General    Contamination   Coding_Density  Contig_N50      Average_Gene_Length     Genome_Size     GC_Content      Total_Coding_Sequences
concoct bin.0   6.14    0.0     0.508   3089    177.1764705882353       88268   0.39    85
concoct bin.1   99.02   0.0     0.733   770501  292.45879556259905      1506883 0.43    1262
concoct bin.10  6.38    0.0     0.96    5441    218.0   5441    0.45    8
```

## Determine taxa for each bin
Here we will use `gtdb-taxo.sh` and `summarize-gtdbtk.sh` to classify all bins and summarize it neatly.
Results will be in their own folder as before with `06-CHECKM2` but with all information from each bin stored in their corresponding folder.

```
07-GTDB-Tk/results/
all-gtdb-report.tsv  CHST-gtdb-report.tsv  COGE-gtdb-report.tsv

ls 07-GTDB-Tk/CHST/metawrap_50_10_bins/
align  classify  gtdbtk.bac120.summary.tsv  gtdbtk.json  gtdbtk.log  gtdbtk.warnings.log  identify
```

## Determining which bin has _Stammerula_ 16S rRNA hits!
Here we will use BLAST+.
First create a database from all the bins you have created using `make-db.sh` and then blast the database with your query using `blastn-localdb.sh`.
Here, many parameters are hard coded. IF you want to change the query from `doc/Erwinia_Stammerula_Wolbachia-16SrRNA.fa`. Then you only need to change the path inside the script.
```
../scripts/blast/
├── align-bin.sh
├── blastn-localdb.sh
├── make-dastool-db.sh
├── make-db.sh
├── make-metawrap-db.sh
└── make-nf-core-db.sh
```
## BUSCO
Lastly, it is BUSCO time!
You will find the script here:
```
../scripts/busco/
└── busco-time.sh
```
Here, you will utilize the path to the bin you want to check, and the taxonomical information you gained from `GTDB-Tk`

## Lastly, lets annotate with Bakta
Here, you first annotate with `bakta-annotate.sh` by setting the `$POP` argument as you run the script.
Notably, the script will avoid any bins that have less than 50% completeness by utilizing the information from earlier steps with `CheckM2`
```
../scripts/bakta/
├── bakta-annotate.sh
└── summarize-bakta.sh
```
Then you can summarize all the bins using the `$POP` argument again.
It will return a `.csv` file holding information of the key genes that are required to achieve **high-quality** by MIMAG's guidelines.
```
Binner,Bin,23S,15S,5S,rRNA,tRNA,ori
concoct,bin.1,1,6,5,12,43,0
concoct,bin.2,1,1,0,2,35,1
maxbin2,bin.0,1,7,5,13,45,5
maxbin2,bin.1,1,1,0,2,35,1
metabat2,bin.1,1,1,0,2,35,1
metabat2,bin.2,1,7,5,13,45,0
metawrap,bin.1,1,6,5,12,43,0
metawrap,bin.2,1,1,0,2,35,1
```


# Now, it is time to summarize all this data.
I used R, i loaded in the files and did some data wrangling.

```
{r}
library(tidyverse)
library(flexplot)


# Load in data
all_bin <- read_delim("results/CHST-bin-report.tsv", delim = "\t")
all_tax <- read_delim("results/CHST-gtdb-report.tsv", delim = "\t")
all_mag <- read_delim("results/CHST-MAG-assembly-stats.csv", delim = ",")
all_anno <- read_delim("results/CHST-bakta-report.csv", delim = ",")

# Rename Name column in all_bin to Bin
all_bin <- all_bin %>% rename(Bin = Name)
all_tax <- all_tax %>% rename(Bin = user_genome)
all_mag <- all_mag %>% rename(Bin = bin, Binner = binner)

# Discard unwanted columns


# Merge the dataframes
merged_data <- all_bin %>%
  inner_join(all_tax, by = c("Binner", "Bin")) %>%
  inner_join(all_mag, by = c("Binner", "Bin")) %>%
  inner_join(all_anno, by = c("Binner","Bin"))

# Change column names
merged_data <- merged_data %>%
  rename(Comp. = Completeness_General,
         Cont. = Contamination,
         CDS_Dens. = Coding_Density,
         N50 = ctg_N50,
         L50 = ctg_L50,
         Avg.Gene.Len. = Average_Gene_Length,
         GC = GC_Content,
         Tot.CDS = Total_Coding_Sequences,
         Taxa_tree = classification,
         Closest_ref = closest_genome_reference,
         ANI = closest_genome_ani,
         AF = closest_genome_af,
         Contigs = n_contigs,
         Contig_max = ctg_max,
         Size = Genome_Size)

merged_data <- merged_data %>% 
  select(-closest_genome_taxonomy, -Contig_N50) %>% 
  select(Binner, Bin, Comp., Cont., Contigs, N50, L50, Contig_max, GC, Size, CDS_Dens.,
         Tot.CDS, Avg.Gene.Len., `23S`, `15S`, `5S`, rRNA, tRNA, ori, ANI, 
         AF, Closest_ref, everything())

# Separate Taxa_tree into individual columns
merged_data <- merged_data %>%
  separate(Taxa_tree, into = c("Domain", "Phyla", "Clade", "Order", "Family", 
                               "Genus", "Species"), sep = ";", fill = "right")

# Remove the prefixes (d__, p__, c__, o__, f__, g__, s__) from the new columns
merged_data <- merged_data %>%
  mutate(across(c(Domain, Phyla, Clade, Order, Family, Genus, Species), ~ sub("^[a-z]__", "", .)))

# Replace missing Genus values with "unspecified"
merged_data <- merged_data %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", "unspecified", Genus))


# Metawrap dataframe
metawrap_data <- merged_data %>%
  select(Binner, Bin, Comp., Cont., Contigs, N50, L50, Contig_max, GC, Size,
         CDS_Dens., Tot.CDS, Avg.Gene.Len., `23S`, `15S`, `5S`, rRNA, tRNA, ori,
         ANI, AF, Closest_ref, Domain, Phyla, Clade, Order, Family, Genus, 
         Species) %>% 
  filter(Binner == "metawrap") %>% 
  select(-Binner)


# Convert all columns to character type
metawrap_data <- metawrap_data %>% mutate(across(everything(), as.character))

# Step 1: Convert to longer format, excluding Bin
wrap_long_df <- metawrap_data %>%
  pivot_longer(cols = -Bin, names_to = "Metric", values_to = "value")

# Step 2: Pivot to wider format using Bin as the new headers
metawrap_data <- wrap_long_df %>%
  pivot_wider(names_from = Bin, values_from = value)

# Save the metawrap_data as a csv file
write_csv(metawrap_data, "results/COGE-metawrap-MAG-stats.csv")

```

