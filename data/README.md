# Data generation

The following is a description of some data analysis steps used for the paper
and not included in the supplied R code. It aims at providing all necessary
information to go from the raw FASTQ files output by the Illumina sequencing
machines to the processed data used to build all figures in the paper.

Some tasks were performed on the NYU HPC platform (running the
 [slurm](https://slurm.schedmd.com/) workload manager).

### List of included analyses

* [Downloading raw data](#downloading-raw-data)
* [Preparation of hybrid SK1-S288C yeast genome](README.md#Preparation-of-hybrid-SK1-S288C-yeast-genome)
* [List of SNPs between SK1 and S288C](README.md#List-of-SNPs-between-SK1-and-S288C)
* [Proteome sequence identity between SK1 and S288C](README.md#Proteome-sequence-identity-between-SK1-and-S288C)
* [Count number of aligned reads per chromosome](README.md#Count-number-of-aligned-reads-per-chromosome)
* [FASTQ subsampling and counting aligned reads](README.md#FASTQ-subsampling-and-counting-aligned-reads)
* [Summary of aligned read coverage per genomic position](README.md#Summary-of-aligned-read-coverage-per-genomic-position)
* [Red1 peaks on SK1 and S288C genomes](README.md#Red1-peaks-on-SK1-and-S288C-genomes)
* [SK1 dataset with simulated S288C spike-in](README.md#SK1-dataset-with-simulated-S288C-spike-in)
* [Influence of read length on signal gaps](README.md#Influence-of-read-length-on-signal-gaps)
* [Influence of strain divergence](README.md#Influence-of-strain-divergence)

## Downloading raw data

The raw sequencing data can be obtained from NCBI's Gene Expression Omnibus
through GEO Series accession number
[GSE115092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115092)
(private until publication). The following code uses NCBI's Entrez Direct and
SRA Toolkit (see [Download Tools](https://www.ncbi.nlm.nih.gov/home/tools/)) to
download all SRR files in the project and dump them to FASTQ files (will work
  after publication of the paper, once the submission is public).

> _Please note that the following code will download many relatively large
files. Alternatively, you can download selected files using their individual
SRA code. For example, download "WT Input rep1" by running:_
`$ fasterq-dump SRX4140900`

```bash
# Use GSE115092's BioProject name in query
esearch -db sra -query PRJNA473777 | \
efetch --format runinfo | cut -d ',' -f 1 | \
grep SRR | xargs fasterq-dump
```

## Preparation of hybrid SK1-S288C yeast genome

```bash
mkdir S288C_SK1_Yue_hybrid_genome && cd S288C_SK1_Yue_hybrid_genome

# Download nuclear genomes
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/SK1.all_feature.gff.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/S288C.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/S288C.all_feature.gff.gz
gunzip *

# SK1 FASTA has 32 lines, S288c has 221656;
# Unlikely to make any difference, but delete excess new lines
# in S288C.genome.fa for consistency
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' S288C.genome.fa \
>> new_S288C.genome.fa
mv new_S288C.genome.fa S288C.genome.fa

# Add strain name to chr names (to make them distinguishable)
sed -i -E 's/(>chr[IVX]+)/\1_S288C/' S288C.genome.fa
sed -i -E 's/(>chr[IVX]+)/\1_SK1/' SK1.genome.fa

# Concatenate genomes (cat alone would not introduce new line between files)
cat S288C.genome.fa <(echo) SK1.genome.fa > S288c_SK1_Yue.fa
```


## List of SNPs between SK1 and S288C

Align SK1 and S288C genome assemblies
([Yue _et al._, Nat Genet 2017](https://www.ncbi.nlm.nih.gov/pubmed/28416820))
and get SNPs using [Mauve](http://darlinglab.org/mauve/mauve.html)
([Darling _et al._, Genome Res 2014](https://www.ncbi.nlm.nih.gov/pubmed/15231754))
as described below.

#### Download genomes

```bash
mkdir Yue_SK1_v_S288c_SNPs
cd Yue_SK1_v_S288c_SNPs

wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/SK1.all_feature.gff.gz

wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/S288C.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/S288C.all_feature.gff.gz

gunzip *
```

#### Generate list of SNPs

* Select **`Align with progressiveMauve...`** in the Mauve GUI
    * Select genome files and define output file name:
        - `S288C.genome.fa`
        - `SK1.genome.fa`
        - Output file name: `S288c_v_SK1.aln`
* Select **`Tools > Export > Export SNPs`** in the Mauve GUI
    * Define output file name:
        - Output file names: `S288c_v_SK1.snp`

#### Format SNP table for use

Run the following R code.

```r
# Load SNPs
S288cvSK1_snp <- readr::read_tsv(here('data/S288c_v_SK1.snp'))
head(S288cvSK1_snp, 10)

message('Number of identified SNPs: ', nrow(S288cvSK1_snp))

# SNPs on different contigs:
S288cvSK1_snp_diff_contig <- S288cvSK1_snp[S288cvSK1_snp$sequence_1_Contig != S288cvSK1_snp$sequence_2_Contig, ]

head(S288cvSK1_snp_diff_contig)

message('Number of SNPs on different chromosomes: ', nrow(S288cvSK1_snp_diff_contig))
message('(Remaining SNPs on matching chromosomes: ',
        nrow(S288cvSK1_snp[S288cvSK1_snp$sequence_1_Contig == S288cvSK1_snp$sequence_2_Contig, ]))

# Drop the 236 non-matched SNPs
S288cvSK1_snp <- S288cvSK1_snp[S288cvSK1_snp$sequence_1_Contig == S288cvSK1_snp$sequence_2_Contig, ]

# Make tidy table
SNPs <- S288cvSK1_snp[, c('sequence_1_Contig', 'sequence_1_PosInContg',
                          'sequence_2_Contig', 'sequence_2_PosInContg')]

SNPs$sequence_1_Contig <- paste0(SNPs$sequence_1_Contig, '_S288C')
SNPs$sequence_2_Contig <- paste0(SNPs$sequence_2_Contig, '_SK1')

colnames(SNPs) <- c('chr_S288C', 'position_S288C', 'chr_SK1', 'position_SK1')

S288C <- subset(SNPs, select=c(chr_S288C, position_S288C))
SK1 <- subset(SNPs, select=c(chr_SK1, position_SK1))
colnames(S288C) <- c('chr', 'position')
colnames(SK1) <- c('chr', 'position')
SNPs <- rbind(S288C, SK1)

head(SNPs)

# Save new table to file
readr::write_tsv(SNPs, here('data/S288c_v_SK1.snp'))
```

## Proteome sequence identity between SK1 and S288C

Align sequences for several example genes in order to highlight differences in
homology at the DNA and amino acid level. The task was performed using the
following steps:

* Extract gene sequences from whole genome FASTA files of the S288C and SK1
strains;
* Download proteome sequences;
* Align both DNA and amino acid sequences and calculate average homology at the
DNA and protein level.

#### Gene sequences

Download genome files as above. Keep only genes in GFF file and use those to get
gene sequences from whole genome FASTA files.

```bash
# Check available feature types
cut -f 3 S288C.all_feature.gff | sort | uniq -c
cut -f 3 SK1.all_feature.gff | sort | uniq -c

# Keep only 'gene' features
grep 'gene' S288C.all_feature.gff > S288c.genes.gff
grep 'gene' SK1.all_feature.gff > SK1.genes.gff

# Convert GFF to BED format in order to save the gene ID in output FASTA
cut -f1,4-5 S288c.genes.gff > S288c.genes.temp.bed
paste S288c.genes.temp.bed <(cut -f9 S288c.genes.gff | cut -d "=" -f3) > S288c.genes.bed

cut -f1,4-5 SK1.genes.gff > SK1.genes.temp.bed
paste SK1.genes.temp.bed <(cut -f9 SK1.genes.gff | cut -d "=" -f3) > SK1.genes.bed

rm *.temp.bed

# Extract DNA sequences
bedtools getfasta -name -fi S288C.genome.fa -bed S288c.genes.bed -fo S288c.genes.fa
bedtools getfasta -name -fi SK1.genome.fa -bed SK1.genes.bed -fo SK1.genes.fa
```

#### Download proteomes

```bash
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/S288c.pep.fa.gz

wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_PEP/SK1.pep.fa.gz

gunzip *
```

#### Align example genes

```bash
# Get sequences
GENES=(YAL058W YAL056W YAL055W YAL054C)

for Gene in ${GENES[@]}
do
  grep -A1 $Gene S288C.genes.fa SK1.genes.fa S288c.pep.fa SK1.pep.fa
done

for FILE in SK1.genes.fa S288c.genes.fa
do
  grep -A1 "YBR030W" $FILE | tail -n1
done
```

## Count number of aligned reads per chromosome

The count of aligned reads per chromosome, used to calculate spike-in
normalization factors, was obtained using `samtools idxstats`. The following is
example code that can be run in a directory containing aligned read maps in BAM
format (sorted by position) to generate read count text files.

```bash
for ALN in *_sorted.bam
do
  echo ">>>>> $ALN"
  samtools idxstats $ALN | cut -f 1,3 > stats_${ALN%_sorted.bam}.txt
done
```

## FASTQ subsampling and counting aligned reads

In order to test the robustness of the method to sequencing depth, as well as
unbalanced sequencing depths between strains and between IP and input samples,
we generated subsamples of each FASTQ. This will simulate a range of sequencing
depths. We used the following reference and test strains:

| Strain code | Genotype                |
|-------------|-------------------------|
| AH119       | Wild type               |
| AH7011      | _red1<sub>ycs4S</sub>_  |
| AH9048      | _hphMX4::p(G162A)-Red1_ |

The task was done in the following steps:

* Concatenate replicate FASTQ files from two experiments (in order to have at
least 10 million reads per sample);
* Generate random samples of 1 to 10 million reads without replacement for each
sample
* Calculate the spike-in normalization factor for all combinations of read count
subsamples

Use raw, non-trimmed read sequencing output FASTQ files. The following is the
`bash` code used.

```console
# Concatenate FASTQ files
$ mkdir Mapped_subsamples
$ cd Mapped_subsamples/

$ cat *ah119spikea*gz > AH119_input.fastq.gz
$ cat *ah119spikeb*gz > AH119_ip.fastq.gz
$ cat *ah9048spikea*gz > AH9048_input.fastq.gz
$ cat *ah9048spikeb*gz > AH9048_ip.fastq.gz
$ cat *ah7011spike*inp*gz > AH7011_input.fastq.gz
$ cat *ah7011spike*chip*gz > AH7011_ip.fastq.gz

# Check line counts
$ for FASTQ in *.fastq.gz
> do
>   echo "${FASTQ}:"
>   zcat ${FASTQ} | wc -l
> done

AH119_input.fastq.gz:
55116728
AH119_ip.fastq.gz:
302938588
AH7011_input.fastq.gz:
64043348
AH7011_ip.fastq.gz:
70163716
AH9048_input.fastq.gz:
72933436
AH9048_ip.fastq.gz:
83094632
HLYHHAFXX_n01_ah119spikea-062817.fastq.gz:
14477064
HLYHHAFXX_n01_ah119spikeb-062817.fastq.gz:
129054988
HLYHHAFXX_n01_ah9048spikea-062817.fastq.gz:
47212668
HLYHHAFXX_n01_ah9048spikeb-062817.fastq.gz:
60823480
HNKY2AFXX_n01_ah119spikea-100917.fastq.gz:
40639664
HNKY2AFXX_n01_ah119spikeb-100917.fastq.gz:
173883600
HNKY2AFXX_n01_ah9048spikea-100917.fastq.gz:
25720768
HNKY2AFXX_n01_ah9048spikeb-100917.fastq.gz:
22271152
HW373AFXX_n01_ah7011spike-chip_122017.fastq.gz:
44903352
HW373AFXX_n01_ah7011spike-inp_122017.fastq.gz:
37517104
HWKGGAFXX_n01_ah7011spike-red1-chip_0118.fastq.gz:
25260364
HWKGGAFXX_n01_ah7011spike-red1-inp_0118.fastq.gz:
26526244
```

```bash
# Remove original files
rm H*
```

Get the pipeline and run in `for loop` across subsample read number with a 
nested `for loop` across FASTQ files:

```bash
# Unzip files
gunzip *fastq.gz

### Run pipeline for all read subsample sizes and FASTQ files
### (using Slurm on the NYU HPC)

# Download pipeline
git clone https://github.com/hochwagenlab/SNP-ChIP.git
mv shell_scripts/Spike-in_subset_mapping.slurm .
rm -rf shell_scripts

# Alternatively, just run:
wget https://raw.githubusercontent.com/hochwagenlab/SNP-ChIP/master/shell_scripts/Spike-in_subset_mapping.slurm

# Need to avoid having more than one instance of the pipeline access the same
# file at the same time: submit each read count in "for loop" across files;
# wait for all 4 to finish before launching the next one

# Example for 1 million reads:
for FASTQ in *.fastq
do
  F_NAME=${FASTQ%.fastq}
  echo ">>> ${F_NAME}..."
  sbatch --export EXPID="${F_NAME}",\
    RUNDIR="Mapped_subsamples",\
    FQ=${FASTQ},NREADS=1,\
    GZIP=F Spike-in_subset_mapping.slurm
done
```

The resulting `stats*.txt` files can be found in this directory alongside this
`README.md` file.

## Summary of aligned read coverage per genomic position

Read coverage at each genomic position (or read pileup) was obtained using
`bedtools genomecov`. It was used, for example, for an alternative way of
computing spike-in normalization factors, to compare with the standard way
(using total number of mapped reads).

The following is example code that can be run in a directory containing aligned
read maps in BAM format (sorted by position) to generate coverage files in
BedGraph format.

```bash
for FILE in *_sorted.bam
do
  bedtools genomecov -ibam ${FILE} -bg > Pileup_${FILE%_sorted.bam}.bdg
done
```

The generated bedtools pileups were used to calculate spike-in normalization
factors using different input data. The R code can be found in file
`helper_spikein_normalization_factor_alternative_calculations.R`. The final
results table was saved to file
 `data/spikein_normalization_factors_using-different_methods.csv`.

## Red1 peaks on SK1 and S288C genomes

Generate reference lists of Red1 peak annotations by aligning reads from three
replicate Red1 ChIP-seq experiments of wild type SK1 strain to both the SK1 and
S288C genome assemblies
([Yue _et al._, Nat Genet 2017](https://www.ncbi.nlm.nih.gov/pubmed/28416820))
and then calling peaks using [MACS2](https://github.com/taoliu/MACS)
([Zhang _et al._, Genome Biol 2008](https://www.ncbi.nlm.nih.gov/pubmed/18798982)).

The following is example code illustrating the procedure for the S288C genome
assembly.

```bash
# Get S288C Yue genome
mkdir S288C && cd S288C_genome

wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/S288C.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/S288C.all_feature.gff.gz

# Build bowtie index
bowtie-build -f S288C.genome.fa S288c_Yue

# Align all fastq files (3 replicates of WT Red1 ChIP-seq)
cd ..
mkdir Unaligned
mkdir Bowtie

for FILE in *.fastq
do
  echo "--------------"
  echo ">>>>> Align ${FILE}"
  echo "--------------"
  bowtie -q -m 1 -v 0 -p 8 -S \
    --un Unaligned/${FILE%.fastq}_S288c_Unaligned.fastq \
    --max Unaligned/${FILE%.fastq}_S288c_Max.fastq \
    S288C_genome/S288c_Yue ${FILE} Bowtie/${FILE%.fastq}_S288c.sam
done

# Call peaks with MACS2
EXPID=Red1-wildtype-Reps-S288C

mkdir $EXPID
cd $EXPID

CHIP="../Bowtie/WT_A_ip_S288c.sam ../Bowtie/WT_B_ip_S288c.sam ../Bowtie/WT_C_ip_S288c.sam"
INPUT="../Bowtie/WT_A_input_S288c.sam ../Bowtie/WT_B_input_S288c.sam ../Bowtie/WT_C_input_S288c.sam"

macs2 callpeak -t $CHIP -c $INPUT \
  --keep-dup="auto" -B --nomodel \
  --extsize 200 --SPMR -g 1.2e7 \
  -n $EXPID &> ${EXPID}.out # save "stdout" and "stderr" to file
```

## SK1 dataset with simulated S288C spike-in

In order to compare data with and without spike-in directly, we generated a
simulated spike-in experiment by taking Red1 ChiP-seq experiment data and adding
a controlled amount of S288C strain reads _in silico_. We used data from a
cohesin subunit Scc1 (6HA-tagged) ChIP-seq experiment using metaphase-arrested
wild type S288c cells
([Hinshaw _et al._, eLife 2015](https://elifesciences.org/content/4/e06057/)).

#### Get data

Downloaded the IP and input data, using the
[`SRA Toolkit`](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc).
The datasets were one of our lab's previously published Red1 ChIP-seq
experiments in wild type SK1
(GEO record [GSM1695718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1695718))
and a Scc1-6HA sample in S288C from Hinshaw _et al._
(GEO record [GSE68573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68573)):

* GSM1695716	AH119 input at 3h (19'838'251 reads)
* GSM1695718	AH119 anti-Red1 at 3h (22'907'545 reads)
* GSM1675773	AM1145_Scc1-6HA_2hrNocodazole_Input (17'812'053 reads)
* GSM1675776	AM1145_Scc1-6HA_2hrNocodazole_IP (20'734'852 reads)

```bash
mkdir Spike-in_simulation
cd Spike-in_simulation

# Check first 5 spots of dataset
fastq-dump -X 5 -Z SRR2040146
fastq-dump -X 5 -Z SRR2040148
fastq-dump -X 5 -Z SRR2009025
fastq-dump -X 5 -Z SRR2009028

# Download
fastq-dump SRR2040146
fastq-dump SRR2040148
fastq-dump SRR2009025
fastq-dump SRR2009028
```

#### Combine SK1 and S288C reads

In order to simulate a spike-in experiment, we combined reads from the SK1 and
the S288C experiments (at a proportion of 80%:20%: 8 million plus 2 million
reads).


```bash
# Take 8M / 2M reads from each file:
seqtk sample -s100 SRR2040146.fastq 8000000 > SRR2040146_sub.fastq
seqtk sample -s100 SRR2040148.fastq 8000000 > SRR2040148_sub.fastq
seqtk sample -s100 SRR2009025.fastq 2000000 > SRR2009025_sub.fastq
seqtk sample -s100 SRR2009028.fastq 2000000 > SRR2009028_sub.fastq

# Combine reads
cat SRR2040146_sub.fastq SRR2009025_sub.fastq > input.fastq
cat SRR2040148_sub.fastq SRR2009028_sub.fastq > chip.fastq
```

#### Align reads and call peaks

For direct comparison, run the read alignment and MACS2 peak calling pipelines
as follows:

* Align reads of the simulated spike-in sample to the hybrid SK1-S288C genome;
* Align reads of the subsetted Red1 ChiP-seq samples to the SK1 genome.

Finally, compare called peaks between the two cases.

Get the read alignment and MACS2 pipelines
(`shell_scripts/ChIPseq_Pipeline_hybrid_genome/ChIPseq_Pipeline_hybrid_genome.sbatch`)
and run FASTQ files.


```bash
### Run pipelines (using Slurm on the NYU HPC)
sbatch --export EXPID="Simulated_spike-in_YueSK1_S288C_PM_SPMR",\
  RUNDIR="Spike-in_simulation",\
  CHIP="chip.fastq",\
  INPUT="input.fastq",\
  GENNAME="S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue" \
  Spike-in_simulation/ChIPseq_Pipeline_hybrid_genome.sbatch

sbatch --export \
  INPUT=SRR2040146_sub.fastq,\
  CHIP=SRR2040148_sub.fastq,\
  GEN="SK1Yue-PM",\
  TAGI="AH119_input",\
  TAGC="AH119_Red1_chip",\
  PEAK="BOTH" \
  ~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

## Influence of read length on signal gaps

In order to test the impact of increasing read length on the signal gaps due to
lack of sequence polymorphisms, we used data from a separate ChIP-seq study
sequencing hybrid SK1/S288C strains using 150-nt long reads. These produce a mix
of about 50-50% reads from each genome. To obtain an unbiased comparison we
performed the following steps:

* Produce read pileups from the original 150-nt long reads;
* Trim the reads down to 100 nt and produce signal pileups;
* Trim the reads down to 50 nt and produce signal pileups;
* Compare percentage of positions with zero signal between the three cases and
with regular, non-hybrid data.

#### Raw data

Download data from study
[GSE114731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114731),
namely sample GSM3148970 (SK1 WT / S288c WT Input rep 1). Also download a
regular ChIP-seq sample (no spike-in) for comparison. Use sample GSM2320246
(Wild type input rep3) from study
[GSE87060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87060),

Get random subsets of 15 million reads.

```bash
seqtk sample -s42 BS00251A_S1_R1_001.fastq.gz 15000000 > Reads.fastq
seqtk sample -s42 SRR4254718.fastq 15000000 > Control_reads.fastq
```

#### Trim reads to shorter lengths

```bash
fastx_trimmer -l 50 -v -i Reads.fastq -o Reads_50bp.fastq
fastx_trimmer -l 100 -v -i Reads.fastq -o Reads_100bp.fastq
```

#### Obtain read pileups

```console
$ srun -t1:00:00 --mem=32000 --pty /bin/bash

# Align reads to hybrid genome
$ module load samtools/intel/1.3.1
$ module load bowtie/gnu/1.2.0

$ for FILE in *.fastq
> do
>   bowtie -q -m 1 -v 0 -p 8 \
>   -S /home/lv38/Library/S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue \
>   ${FILE} | samtools view -bS | samtools sort -o ${FILE%.fastq}.bam
>done

# Alignment stats
# reads processed: 15000000
# reads with at least one reported alignment: 1897743 (12.65%)
# reads that failed to align: 11534512 (76.90%)
# reads with alignments suppressed due to -m: 1567745 (10.45%)
Reported 1897743 alignments to 1 output stream(s)

# reads processed: 15000000
# reads with at least one reported alignment: 3421738 (22.81%)
# reads that failed to align: 7194285 (47.96%)
# reads with alignments suppressed due to -m: 4383977 (29.23%)
Reported 3421738 alignments to 1 output stream(s)

# reads processed: 15000000
# reads with at least one reported alignment: 2885917 (19.24%)
# reads that failed to align: 4065307 (27.10%)
# reads with alignments suppressed due to -m: 8048776 (53.66%)
Reported 2885917 alignments to 1 output stream(s)

# Align control reads (allow two mapping positions)
$ bowtie -q -m 2 -v 0 -p 8 \
> -S /home/lv38/Library/S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue \
> Control_reads.fastq | samtools view -bS | samtools sort -o Control_reads.bam

# reads processed: 15000000
# reads with at least one reported alignment: 11786362 (78.58%)
# reads that failed to align: 2448865 (16.33%)
# reads with alignments suppressed due to -m: 764773 (5.10%)
Reported 11786362 alignments to 1 output stream(s)

# Compute read pileups
$ for FILE in *.bam
> do
>   bedtools genomecov -ibam ${FILE} -bg > Pileup_${FILE%.bam}.bdg
>   bedtools genomecov -ibam Control_reads_no_m.bam -bg > Pileup_Control_reads_no_m.bdg
> done
```

#### Compare percentage of 0-signal positions

Start by expanding the bedgraph to single-nt signal. Then compute percentage of
zeros in the genome.

Run the following R code.

```r
library(here)
library(stringr)
library(GenomicRanges)
source(here('helper_io.R'))

# Load pileups
pileups <- list(
  '50bp'=import_bedGraph('~/Desktop/Pileups/Pileup_Reads_50bp.bdg'),
  '100bp'=import_bedGraph('~/Desktop/Pileups/Pileup_Reads_100bp.bdg'),
  '150bp'=import_bedGraph('~/Desktop/Pileups/Pileup_Reads.bdg'),
  'Control'=import_bedGraph('~/Desktop/Pileups/Pileup_Control_reads.bdg')
)

# # Rename chromosomes in control sample to match other samples
# pileups$Control <- renameSeqlevels(
#   pileups$Control, paste0(seqlevels(pileups$Control), '_SK1'))

### Expand to single-bp genomic positions (for direct comparison)
# Prepare genome coordinates
SK1_gff <- rtracklayer::import.gff(here('data/GFF/SK1.all_feature.gff'))
S288C_gff <- rtracklayer::import.gff(here('data/GFF/S288C.all_feature.gff'))

SK1_start <- SK1_gff[SK1_gff$type == 'centromere']@ranges@start
SK1_end <- SK1_start + SK1_gff[SK1_gff$type == 'centromere']@ranges@width - 1

S288C_start <- S288C_gff[S288C_gff$type == 'centromere']@ranges@start
S288C_end <-
  S288C_start + S288C_gff[S288C_gff$type == 'centromere']@ranges@width - 1

coord_table <- data.frame(
  "Chromosome" = c(paste0('chr', as.roman(1:16), '_SK1'),
                   paste0('chr', as.roman(1:16), '_S288C')),
  "Start" = c(SK1_start, S288C_start),
  "End" = c(SK1_end, S288C_end),
  "LenChr" = c(
    c(228861, 829469, 340914, 1486921, 589812, 299318,
      1080440, 542723, 449612, 753937, 690901, 1054145,
      923535, 791982, 1053869, 946846),
    c(219929, 813597, 341580, 1566853, 583092, 271539,
      1091538, 581049, 440036, 751611, 666862, 1075542,
      930506, 777615, 1091343, 954457)))

genome_info <- with(
  coord_table, GRanges(Chromosome, IRanges(Start + 1, End),
                       seqlengths=setNames(LenChr, Chromosome)))

# Work with SK1 chromosomes only
SK1_chrs <- paste0('chr', as.roman(1:16), '_SK1')

for (i in seq_along(pileups)) {
  pileups[[i]] <- keepSeqlevels(pileups[[i]], SK1_chrs, pruning.mode="coarse")
}

genome_info <- keepSeqlevels(genome_info, SK1_chrs, pruning.mode="coarse")

# Sort sequences and levels to make sure they match
genome_info <- sortSeqlevels(genome_info)
pileups <- lapply(pileups, sortSeqlevels)
pileups <- lapply(pileups, sort)

# Add info to signal objects
for (i in seq_along(pileups)) {
  seqlengths(pileups[[i]]) <- seqlengths(genome_info)
}

# Compute 1-bp tiling windows
bins <- tileGenome(
  seqlengths(pileups[[1]]), tilewidth=1, cut.last.tile.in.chrom=TRUE)

# Get signal as "RleList"; the signal is stored in the "score" metadata column
scores <- lapply(pileups, coverage, weight="score")

# Get signal per bp
binned_scores <- list()
for (i in seq_along(scores)) {
  binned_scores[[i]] <- binnedAverage(bins, scores[[i]], "binned_score")
}

# Percentage of positions with signal 0 in each sample (increasing read length)
results <- list()
results <- lapply(
  binned_scores, function(x) sum(x$binned_score == 0) / length(x) * 100)

names(results) <- names(pileups)

for (i in seq_along(results)) message(
  names(results)[i], ': ', round(results[[i]], 2), '%')
```

## Influence of strain divergence

Test the minimal required level of genetic divergence between the test and the
spike-in strains, using an _in-silico_ experiment. This uses the SK1 strain to
simulate pairs of strains with increasing divergence. The analysis includes the
following steps:

* Introduce increasing numbers of random mutations in the SK1 genome, to
simulate spike-in genomes with increasing genetic distance;
* Generate sequencing reads for the original and mutated sequences and mix them
at two different proportions: 80:20% as the “wild type” condition and 50:50% for
a simulated low target protein condition;
* Map read mixes to hybrid genomes made of the original and the mutated
sequences
* Compute spike-in normalization factors and check stability with decreasing
genetic divergence.

#### Introduce random mutations

Use a small application written for this task:
[SNPr](https://github.com/luisvalesilva/snpr).

```bash
mkdir Genetic_divergence && cd Genetic_divergence

# Run interactive session on HPC
srun -t5:00:00 --mem=32000 --pty /bin/bash

# Download nuclear SK1 genome
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz
gunzip SK1.genome.fa.gz

# Get SNPr
git clone https://github.com/luisvalesilva/snpr.git
cd snpr

# Install dependencies and activate virtual environment
module load python3/intel/3.6.3
pip install --user pipenv
pipenv install
pipenv shell

# Generate mutated SK1 genomes:
# introduce progressively larger numbers of mutations, starting from 0.001% of
# the positions

for FREQ in 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05
do
  echo ">>>>> SNP frequency: $FREQ..."
  python snpr.py -f $FREQ -r 42 ../SK1.genome.fa > ../SK1_${FREQ}mut.fa
done

# Exit virtual environment
exit
```

#### Generate sequencing reads

Use the NIH program
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm).

```bash
# Download and unpack ART binaries
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
tar -xvzf artbinmountrainier2016.06.05linux64.tgz
rm artbinmountrainier2016.06.05linux64.tgz

# Single-end read simulation
# Generate 100-bp reads from HiSeq 2000 Illumina sequencing system
# (in order to be able to generate 100-bp reads)
# Total number of 5'000'000 reads (target + spike-in)

# Reads for original SK1 genome
for NREADS in 2500000 4000000 4500000
do
  echo
  echo ">>>>> Read #: $NREADS"
  art_bin_MountRainier/art_illumina -ss HS20 -i SK1.genome.fa \
  -l 100 -c `expr $NREADS / 16` -o SK1_${NREADS}reads
done

# Reads for simulated (spike-in) mutated SK1 genome
for REF in SK1_0.0*
do
  echo
  echo ">>>>> Running $REF"
  for NREADS in 500000 1000000 2500000
  do
    echo ">>>>>      Read #: $NREADS"
    art_bin_MountRainier/art_illumina -ss HS20 -i $REF \
  -l 100 -c `expr $NREADS / 16` -o ${REF%.fa}_${NREADS}reads
  done
done
```

#### Mix simulated reads

Mixing simulated reads for the target reference genome and the simulated
(mutated) spike-in genomes at the proportions of 80:20%, 50:50%, and 90:10%
(target:spike-in%) will simulate an experiment in which we added 20% spike-in
cells and got the same ratio of sequencing reads in the “wild type” condition
and changes to 50:50% for a condition with low target protein and and 90:10% for
a condition with high target protein.

Here I will simply concatenate pairs of FASTQ files. The input samples are
simply mixes at a 80:20% proportion in both cases (use the "Spiked_WT" FASTQ
file for both).

```bash
# "Wild-type" condition
# Concatenate FASTQ files (make sure files end in newline "\n" character)

for FQ in SK1_0*1000000reads.fq
do
  echo ">>>>> Running $FQ"
  cat SK1_4000000reads.fq $FQ > Spiked_WT_${FQ%_1000000reads.fq}.fq
done

# "Low-target" condition
for FQ in SK1_0*2500000reads.fq
do
  echo ">>>>> Running $FQ"
  cat SK1_2500000reads.fq $FQ > Spiked_low_target_${FQ%_2500000reads.fq}.fq
done

# "High-target" condition
for FQ in SK1_0*_500000reads.fq
do
  echo ">>>>> Running $FQ"
  cat SK1_4500000reads.fq $FQ > Spiked_high_target_${FQ%_500000reads.fq}.fq
done
```

#### Prepare hybrid genomes

Rename chromosomes, in order to make them distinguishable and then concatenate
genome FASTA files (make sure files end in newline "\n" character).

```bash
# Add "mut" to mutated genome chr names (to make them distinguishable)
for FA in SK1*mut.fa
do
  echo ">>>>> Running $FA"
  sed -i -E 's/(>chr[IVX]+)/\1_mut/' $FA
done

# Concatenate genomes
for FA in SK1*mut.fa
do
  echo ">>>>> Running $FA"
  cat SK1.genome.fa $FA > Hybrid_${FA}
done
```

#### Map read mixes to hybrid genomes

Map the generated read mixes at the different proportions to the hybrid genomes
prepared by concatenating the respective reference sequences.

```bash
# Build bowtie index for each genome
mkdir Hyb_ref_gen
mv Hybrid_SK1_0.0* Hyb_ref_gen/
cd Hyb_ref_gen

# Run interactive session on HPC
srun -t5:00:00 --mem=32000 --pty /bin/bash
module load bowtie/gnu/1.2.0

for FA in Hybrid*.fa
do
  echo ">>>>> Running $FA"
  bowtie-build -f $FA ${FA%.fa}
done

# Align mixed FASTQ files
cd ..
mkdir Bowtie

for FREQ in 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.05
do
  for FA in Hybrid_SK1_${FREQ}mut.fa
  do
    for FQ in Spiked*${FREQ}mut.fq
    do
      echo "--------------"
      echo ">>>>> Align ${FQ} to ${FA}"
      echo "--------------"
      bowtie -q -m 1 -v 0 -p 8 -S \
      Hyb_ref_gen/${FA%.fa} ${FQ} Bowtie/${FQ%.fq}.sam
    done
  done
done
```

And here's the console output.

```console
--------------
>>>>> Align Spiked_high_target_SK1_0.00001mut.fq to Hybrid_SK1_0.00001mut.fa
--------------
# reads processed: 4994818
# reads with at least one reported alignment: 2387 (0.05%)
# reads that failed to align: 2636579 (52.79%)
# reads with alignments suppressed due to -m: 2355852 (47.17%)
Reported 2387 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.00001mut.fq to Hybrid_SK1_0.00001mut.fa
--------------
# reads processed: 4995005
# reads with at least one reported alignment: 2414 (0.05%)
# reads that failed to align: 2640716 (52.87%)
# reads with alignments suppressed due to -m: 2351875 (47.08%)
Reported 2414 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.00001mut.fq to Hybrid_SK1_0.00001mut.fa
--------------
# reads processed: 4994695
# reads with at least one reported alignment: 2462 (0.05%)
# reads that failed to align: 2637559 (52.81%)
# reads with alignments suppressed due to -m: 2354674 (47.14%)
Reported 2462 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.00005mut.fq to Hybrid_SK1_0.00005mut.fa
--------------
# reads processed: 4994874
# reads with at least one reported alignment: 10981 (0.22%)
# reads that failed to align: 2636543 (52.78%)
# reads with alignments suppressed due to -m: 2347350 (47.00%)
Reported 10981 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.00005mut.fq to Hybrid_SK1_0.00005mut.fa
--------------
# reads processed: 4994866
# reads with at least one reported alignment: 11358 (0.23%)
# reads that failed to align: 2639704 (52.85%)
# reads with alignments suppressed due to -m: 2343804 (46.92%)
Reported 11358 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.00005mut.fq to Hybrid_SK1_0.00005mut.fa
--------------
# reads processed: 4994773
# reads with at least one reported alignment: 11135 (0.22%)
# reads that failed to align: 2637390 (52.80%)
# reads with alignments suppressed due to -m: 2346248 (46.97%)
Reported 11135 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.0001mut.fq to Hybrid_SK1_0.0001mut.fa
--------------
# reads processed: 4994832
# reads with at least one reported alignment: 22624 (0.45%)
# reads that failed to align: 2636081 (52.78%)
# reads with alignments suppressed due to -m: 2336127 (46.77%)
Reported 22624 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.0001mut.fq to Hybrid_SK1_0.0001mut.fa
--------------
# reads processed: 4994872
# reads with at least one reported alignment: 22775 (0.46%)
# reads that failed to align: 2637113 (52.80%)
# reads with alignments suppressed due to -m: 2334984 (46.75%)
Reported 22775 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.0001mut.fq to Hybrid_SK1_0.0001mut.fa
--------------
# reads processed: 4994731
# reads with at least one reported alignment: 22514 (0.45%)
# reads that failed to align: 2637034 (52.80%)
# reads with alignments suppressed due to -m: 2335183 (46.75%)
Reported 22514 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.0005mut.fq to Hybrid_SK1_0.0005mut.fa
--------------
# reads processed: 4994847
# reads with at least one reported alignment: 108900 (2.18%)
# reads that failed to align: 2636575 (52.79%)
# reads with alignments suppressed due to -m: 2249372 (45.03%)
Reported 108900 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.0005mut.fq to Hybrid_SK1_0.0005mut.fa
--------------
# reads processed: 4994856
# reads with at least one reported alignment: 110993 (2.22%)
# reads that failed to align: 2639761 (52.85%)
# reads with alignments suppressed due to -m: 2244102 (44.93%)
Reported 110993 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.0005mut.fq to Hybrid_SK1_0.0005mut.fa
--------------
# reads processed: 4994755
# reads with at least one reported alignment: 109471 (2.19%)
# reads that failed to align: 2636996 (52.80%)
# reads with alignments suppressed due to -m: 2248288 (45.01%)
Reported 109471 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.001mut.fq to Hybrid_SK1_0.001mut.fa
--------------
# reads processed: 4994816
# reads with at least one reported alignment: 215447 (4.31%)
# reads that failed to align: 2635989 (52.77%)
# reads with alignments suppressed due to -m: 2143380 (42.91%)
Reported 215447 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.001mut.fq to Hybrid_SK1_0.001mut.fa
--------------
# reads processed: 4994881
# reads with at least one reported alignment: 218503 (4.37%)
# reads that failed to align: 2639082 (52.84%)
# reads with alignments suppressed due to -m: 2137296 (42.79%)
Reported 218503 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.001mut.fq to Hybrid_SK1_0.001mut.fa
--------------
# reads processed: 4994750
# reads with at least one reported alignment: 216693 (4.34%)
# reads that failed to align: 2636355 (52.78%)
# reads with alignments suppressed due to -m: 2141702 (42.88%)
Reported 216693 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.005mut.fq to Hybrid_SK1_0.005mut.fa
--------------
# reads processed: 4994810
# reads with at least one reported alignment: 889903 (17.82%)
# reads that failed to align: 2634059 (52.74%)
# reads with alignments suppressed due to -m: 1470848 (29.45%)
Reported 889903 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.005mut.fq to Hybrid_SK1_0.005mut.fa
--------------
# reads processed: 4994795
# reads with at least one reported alignment: 908599 (18.19%)
# reads that failed to align: 2636734 (52.79%)
# reads with alignments suppressed due to -m: 1449462 (29.02%)
Reported 908599 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.005mut.fq to Hybrid_SK1_0.005mut.fa
--------------
# reads processed: 4994750
# reads with at least one reported alignment: 893300 (17.88%)
# reads that failed to align: 2636146 (52.78%)
# reads with alignments suppressed due to -m: 1465304 (29.34%)
Reported 893300 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.01mut.fq to Hybrid_SK1_0.01mut.fa
--------------
# reads processed: 4994835
# reads with at least one reported alignment: 1429870 (28.63%)
# reads that failed to align: 2633564 (52.73%)
# reads with alignments suppressed due to -m: 931401 (18.65%)
Reported 1429870 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.01mut.fq to Hybrid_SK1_0.01mut.fa
--------------
# reads processed: 4994877
# reads with at least one reported alignment: 1457878 (29.19%)
# reads that failed to align: 2635716 (52.77%)
# reads with alignments suppressed due to -m: 901283 (18.04%)
Reported 1457878 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.01mut.fq to Hybrid_SK1_0.01mut.fa
--------------
# reads processed: 4994763
# reads with at least one reported alignment: 1435302 (28.74%)
# reads that failed to align: 2633808 (52.73%)
# reads with alignments suppressed due to -m: 925653 (18.53%)
Reported 1435302 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_high_target_SK1_0.05mut.fq to Hybrid_SK1_0.05mut.fa
--------------
# reads processed: 4994848
# reads with at least one reported alignment: 2233008 (44.71%)
# reads that failed to align: 2636258 (52.78%)
# reads with alignments suppressed due to -m: 125582 (2.51%)
Reported 2233008 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_low_target_SK1_0.05mut.fq to Hybrid_SK1_0.05mut.fa
--------------
# reads processed: 4994872
# reads with at least one reported alignment: 2279059 (45.63%)
# reads that failed to align: 2639713 (52.85%)
# reads with alignments suppressed due to -m: 76100 (1.52%)
Reported 2279059 alignments to 1 output stream(s)
--------------
>>>>> Align Spiked_WT_SK1_0.05mut.fq to Hybrid_SK1_0.05mut.fa
--------------
# reads processed: 4994712
# reads with at least one reported alignment: 2245253 (44.95%)
# reads that failed to align: 2636740 (52.79%)
# reads with alignments suppressed due to -m: 112719 (2.26%)
Reported 2245253 alignments to 1 output stream(s)
```

#### Compute spike-in normalization factors

Start by getting the count of aligned reads per chromosome, using
`samtools idxstats`. The alignment maps need to be sorted first.

```bash
function sam_to_sorted_and_indexed_bam() # Convert to BAM file, sort, and index
 {
     local BASE=${1%.sam}
     local BAM=${BASE}.bam
     local S_BAM=${BASE}_sorted.bam

     samtools view -bS $1 > $BAM
     samtools sort -o $S_BAM $BAM
     samtools index $S_BAM

     # Clean up
     rm $1
     rm $BAM

     echo "$S_BAM"
}

# Convert all to indexed and sorted BAM files
cd Bowtie

# Run interactive session on HPC
srun -t6:00:00 --mem=32000 --pty /bin/bash
module load samtools/intel/1.3.1

for ALN in Spiked*.sam
do
    echo ">>>>> $ALN"
    sam_to_sorted_and_indexed_bam $ALN
done

# Count aligned reads per chromosome
for ALN in *_sorted.bam
do
    echo ">>>>> $ALN"
    samtools idxstats $ALN | cut -f 1,3 > stats_${ALN%_sorted.bam}.txt
done
```

Calculate spike-in normalization factors using the following R code.

```r
library(stringr)
library(readr)

#' Compute spike-in normalization factor from total read counts
#'
#' Computes spike-in normalization factor between two spiked-in samples using
#' total counts of aligned reads. Inputs paths to text files containing counts
#' of aligned reads per chromosome of a hybrid SK1:S288C genome.
#' @param ref_chip_counts Either a single or a list of paths to reference ChIP
#' samples' read counts file. No default.
#' @param ref_input_counts Either a single or a list of paths to reference input
#' samples' read counts file. No default.
#' @param test_chip_counts Either a single or a list of paths to test ChIP
#' samples' read counts file. No default.
#' @param test_input_counts Either a single or a list of paths to test input
#' samples' read counts file. No default.
#' @param return_counts Logical indicating whether to return the computed read
#' counts instead of the normalization factor. Defaults to \code{FALSE}.
#' @return Numeric normalization factor.
#' @examples
#' \dontrun{
#' spikein_normalization_factor_from_counts(
#'      ref_chip_counts='Counts_AH119_chip.txt',
#'      ref_input_counts='Counts_AH119_input.txt',
#'      test_chip_counts='Counts_AH8104_chip.txt',
#'      test_input_counts='Counts_AH8104_input.txt')
#'
#' spikein_normalization_factor_from_counts(
#'     ref_chip_counts=list('Counts_AH119_chip_1.txt',
#'                          'Counts_AH119_chip_2.txt',
#'                          'Counts_AH119_chip_3.txt'),
#'     ref_input_counts=list('Counts_AH119_inp_1.txt',
#'                           'Counts_AH119_inp_2.txt',
#'                           'Counts_AH119_inp_3.txt'),
#'     test_chip_counts='Counts_AH8104_chip.txt',
#'     test_input_counts='Counts_AH8104_input.txt')
#' }
#' @export
spikein_normalization_factor_from_counts <- function(
  ref_chip_counts, ref_input_counts, test_chip_counts, test_input_counts,
  return_counts=FALSE) {

  # Put paths in list
  files <- list(ref_chip=ref_chip_counts, ref_input=ref_input_counts,
                test_chip=test_chip_counts, test_input=test_input_counts)

  # Convert each element into list, if not one already
  for (i in seq_along(files)) {
    if (!is.list(files[[i]])) files[[i]] <- list(files[[i]])
  }

  # Print files to read to console
  message('>>> Read alignment count files:')
  for (i in seq_along(files)) {
    for (file in files[[i]]) {
      message('   ', basename(file))
    }
  }    

  message()
  # Read files into tibble in list
  tables <- list()
  for (i in seq_along(files)) {
    tables[[i]] <- sapply(files[[i]], FUN=read_tsv, col_names=F,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(tables) <- names(files)

  message()
  # Get read counts per chromosome
  message('>>> Count reads per genome:')
  counts <- list()
  for (i in seq_along(tables)) {
    counts[[i]] <- sapply(tables[[i]], FUN=sum_per_genome,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(counts) <- names(tables)

  # Add-up counts for replicates (results in nested lists)
  for (i in seq_along(counts)) {
    if (length(counts[[i]]) > 1) {
      total <- counts[[i]][[1]]
      for (j in 2:length(counts[[i]])) {
        total <- total + counts[[i]][[j]]
      }
      counts[[i]] <- total
    } else counts[[i]] <- unlist(counts[[i]])
  }

  if (return_counts) {
    message('---')
    message('Done!')
    return(counts)
  }

  # Compute normalization factor
  result <- normalization_factor(ctrl_input=counts$ref_input,
                                 ctrl_chip=counts$ref_chip,
                                 test_input=counts$test_input,
                                 test_chip=counts$test_chip)

  message('---')
  message('Done!')

  return(result)
}

# Helper functions
sum_per_genome <- function(df) {
  # Compute sum of reads aligned to each genome
  Mut <- sum(
    df[apply(df, 1, function(x) str_detect(x[1],'_mut')), 2])
  SK1 <- sum(
    df[apply(df, 1, function(x) str_detect(x[1], 'chr[XVI]$')), 2])

  # Print result to console
  message('  Mut: ', formatC(Mut, big.mark=",",
                               drop0trailing=TRUE, format="f"))
  message('  SK1: ', formatC(SK1, big.mark=",",
                             drop0trailing=TRUE, format="f"))
  message('      ', round(Mut * 100 / (SK1 + Mut), 1), '% spike-in reads')

  # Return result as named vector
  c('Mut'=Mut, 'SK1'=SK1)
}


normalization_factor <- function(ctrl_input, ctrl_chip,
                                 test_input, test_chip) {
  # Compute Q values
  Q_ctrl_input <- ctrl_input['Mut'] / ctrl_input['SK1']
  Q_ctrl_chip <- ctrl_chip['Mut'] / ctrl_chip['SK1']

  Q_test_input <- test_input['Mut'] / test_input['SK1']
  Q_test_chip <- test_chip['Mut'] / test_chip['SK1']

  # Compute normalization factors
  a_ctrl <- Q_ctrl_input / Q_ctrl_chip
  a_test <- Q_test_input / Q_test_chip

  # Return reference strain-centric normalization factor
  a_test/ a_ctrl
}
```

Run the function on the obtained read counts.

```r
#read_counts <- data.frame(
#  SNP_freq=c('0.001', '0.0025', '0.005', '0.0075', '0.01', '0.025', '0.05'),
#  NF=c(
#    spikein_normalization_factor_from_counts(
#    ref_chip_counts='stats_Spiked_WT_SK1_0.001mut.txt',
#    ref_input_counts='stats_Spiked_WT_SK1_0.001mut.txt',
#    test_chip_counts='stats_Spiked_low_target_SK1_0.001mut.txt',
#    test_input_counts='stats_Spiked_WT_SK1_0.001mut.txt'),
#    spikein_normalization_factor_from_counts(
#      ref_chip_counts='stats_Spiked_WT_SK1_0.0025mut.txt',
#      ref_input_counts='stats_Spiked_WT_SK1_0.0025mut.txt',
#      test_chip_counts='stats_Spiked_low_target_SK1_0.0025mut.txt',
#      test_input_counts='stats_Spiked_WT_SK1_0.0025mut.txt'),
#    spikein_normalization_factor_from_counts(
#      ref_chip_counts='stats_Spiked_WT_SK1_0.005mut.txt',
#      ref_input_counts='stats_Spiked_WT_SK1_0.005mut.txt',
#      test_chip_counts='stats_Spiked_low_target_SK1_0.005mut.txt',
#      test_input_counts='stats_Spiked_WT_SK1_0.005mut.txt'),
#    spikein_normalization_factor_from_counts(
#      ref_chip_counts='stats_Spiked_WT_SK1_0.0075mut.txt',
#      ref_input_counts='stats_Spiked_WT_SK1_0.0075mut.txt',
#      test_chip_counts='stats_Spiked_low_target_SK1_0.0075mut.txt',
#      test_input_counts='stats_Spiked_WT_SK1_0.0075mut.txt'),
#    spikein_normalization_factor_from_counts(
#      ref_chip_counts='stats_Spiked_WT_SK1_0.01mut.txt',
#      ref_input_counts='stats_Spiked_WT_SK1_0.01mut.txt',
#      test_chip_counts='stats_Spiked_low_target_SK1_0.01mut.txt',
#      test_input_counts='stats_Spiked_WT_SK1_0.01mut.txt'),
#    spikein_normalization_factor_from_counts(
#      ref_chip_counts='stats_Spiked_WT_SK1_0.025mut.txt',
#      ref_input_counts='stats_Spiked_WT_SK1_0.025mut.txt',
#      test_chip_counts='stats_Spiked_low_target_SK1_0.0025mut.txt',
#      test_input_counts='stats_Spiked_WT_SK1_0.025mut.txt'),
#    spikein_normalization_factor_from_counts(
#      ref_chip_counts='stats_Spiked_WT_SK1_0.05mut.txt',
#      ref_input_counts='stats_Spiked_WT_SK1_0.05mut.txt',
#      test_chip_counts='stats_Spiked_low_target_SK1_0.05mut.txt',
#      test_input_counts='stats_Spiked_WT_SK1_0.05mut.txt')
#  )
#)

read_counts <- data.frame(
  Condition=rep(c('Low target', 'High target'), 8),
  SNP_freq=c(rep('0.00001', 2), rep('0.00005', 2), rep('0.0001', 2),
             rep('0.0005', 2), rep('0.001', 2), rep('0.005', 2), rep('0.01', 2),
             rep('0.05', 2)),
  NF=c(
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.00001mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.00001mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.00001mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.00001mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.00001mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.00001mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.00001mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.00001mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.00005mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.00005mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.00005mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.00005mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.00005mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.00005mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.00005mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.00005mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.0001mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.0001mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.0001mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.0001mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.0001mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.0001mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.0001mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.0001mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.0005mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.0005mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.0005mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.0005mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.0005mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.0005mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.0005mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.0005mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.001mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.001mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.001mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.001mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.001mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.001mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.001mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.001mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.005mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.005mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.005mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.005mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.005mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.005mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.005mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.005mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.01mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.01mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.01mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.01mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.01mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.01mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.01mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.01mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.05mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.05mut.txt',
      test_chip_counts='stats_Spiked_low_target_SK1_0.05mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.05mut.txt'),
    spikein_normalization_factor_from_counts(
      ref_chip_counts='stats_Spiked_WT_SK1_0.05mut.txt',
      ref_input_counts='stats_Spiked_WT_SK1_0.05mut.txt',
      test_chip_counts='stats_Spiked_high_target_SK1_0.05mut.txt',
      test_input_counts='stats_Spiked_WT_SK1_0.05mut.txt')
  )
)
```

And here's the output.

```console
> read_counts
     Condition SNP_freq        NF
1   Low target  0.00001 0.2599659
2  High target  0.00001 2.2833024
3   Low target  0.00005 0.2978224
4  High target  0.00005 2.5157366
5   Low target   0.0001 0.2639502
6  High target   0.0001 2.1984121
7   Low target   0.0005 0.2535488
8  High target   0.0005 2.2272501
9   Low target    0.001 0.2513364
10 High target    0.001 2.2514906
11  Low target    0.005 0.2510153
12 High target    0.005 2.2133372
13  Low target     0.01 0.2510616
14 High target     0.01 2.2204894
15  Low target     0.05 0.2497198
16 High target     0.05 2.2437482
```

Plot results.

```r
ggplot(read_counts, aes(SNP_freq, NF * 100, colour=Condition, fill=Condition)) +
  geom_hline(aes(yintercept = 100), linetype = 3) +
  geom_point(size=1.5) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title = '', x = '', y = 'Red1 amount\n(% of WT)')
```
