
```markdown
# Bioinformatics Workflow for E. coli X Assembly

## 1. Downloading and Checking Files

### Libraries to Download

Download three libraries from the TY2482 sample with the following insert sizes and orientation:

| Library ID   | Type         | Insert Size | Size       |
|--------------|--------------|-------------|------------|
| SRR292678    | Paired End   | 470 bp      | 400 Mb each |
| SRR292862    | Mate Pair    | 2 kb        | 200 Mb each |
| SRR292770    | Mate Pair    | 6 kb        | 200 Mb each |

### Download Commands

Use the following commands to download the files:

```bash
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R2_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz 
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R2_001.fastq.gz
```

## Quality Control

### FastQC Installation

Download FastQC software via `apt-get` and apply the following commands for quality control:

```bash
fastqc -o ~/bio/project3/ ~/bio/project3/SRR292678sub_S1_L001_R1_001.fastq.gz
fastqc -o ~/bio/project3/ ~/bio/project3/SRR292678sub_S1_L001_R2_001.fastq.gz
fastqc -o ~/bio/project3/ ~/bio/project3/SRR292770_S1_L001_R1_001.fastq.gz
fastqc -o ~/bio/project3/ ~/bio/project3/SRR292770_S1_L001_R2_001.fastq.gz
fastqc -o ~/bio/project3/ ~/bio/project3/SRR292862_S2_L001_R1_001.fastq.gz
fastqc -o ~/bio/project3/ ~/bio/project3/SRR292862_S2_L001_R2_001.fastq.gz
```

### Generated HTML Files

These commands will generate the following HTML files for quality check:

| HTML File                                           |
|-----------------------------------------------------|
| SRR292678sub_S1_L001_R1_001_fastqc.html            |
| SRR292678sub_S1_L001_R2_001_fastqc.html            |
| SRR292770_S1_L001_R1_001_fastqc.html               |
| SRR292770_S1_L001_R2_001_fastqc.html               |
| SRR292862_S2_L001_R1_001_fastqc.html               |
| SRR292862_S2_L001_R2_001_fastqc.html               |

#### Quality Checking 

Check the quality of the sequences using the generated FastQC reports.

- SRR292678sub_S1_L001_R1_001_fastqc
- SRR292678sub_S1_L001_R2_001_fastqc
- SRR292770_S1_L001_R1_001_fastqc
- SRR292770_S1_L001_R2_001_fastqc
- SRR292862_S2_L001_R1_001_fastqc
- SRR292862_S2_L001_R2_001_fastqc

#### Common Feature

The common feature of all three runs is the CG-content.

## 2. Assemble Genomes Using SPAdes

### Initial SPAdes Test

Assemble a single library of sequencing (paired end) reads from E. coli X from the library SRR292678 (forward and reverse).

First, try SPAdes:

```bash
spades.py --test
```

**Output**: Thank you for using SPAdes! If you use it in your research, please cite:
> Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A. and Korobeynikov, A., 2020. Using SPAdes de novo assembler. Current protocols in bioinformatics, 70(1), p.e102. doi.org/10.1002/cpbi.102

### Target Files Assembly

Now, assemble our target files. Refer to [SPAdes documentation](https://github.com/ablab/spades/blob/main/docs/running.md) for more details. The command is as follows:

```bash
spades.py --pe1-1 SRR292678sub_S1_L001_R1_001.fastq.gz --pe1-2 SRR292678sub_S1_L001_R2_001.fastq.gz -o spades_output -t 4
```

### Memory Issue

At this point, I ran out of memory and gave up. Instead, I downloaded the pre-made data.

### Download QUAST

Follow these steps to download QUAST:

```bash
mamba create -n quast_env quast -c bioconda -c conda-forge
conda activate quast_env
quast.py contigs.fasta -o quast_results
```

### QUAST Result Folder Contents

The QUAST result folder contains the following files:

```plaintext
drwxr-xr-x 4 lena lena   4096 Dec  3 00:28 ./
drwxrwxr-x 3 lena lena   4096 Dec  3 00:28 ../
drwxr-xr-x 2 lena lena   4096 Dec  3 00:28 basic_stats/
-rw-r--r-- 1 lena lena  53346 Dec  3 00:28 icarus.html
drwxr-xr-x 2 lena lena   4096 Dec  3 00:28 icarus_viewers/
-rw-r--r-- 1 lena lena   2943 Dec  3 00:28 quast.log
-rw-r--r-- 1 lena lena 379072 Dec  3 00:28 report.html
-rw-r--r-- 1 lena lena   32001 Dec  3 00:28 report.pdf
-rw-r--r-- 1 lena lena   1212 Dec  3 00:28 report.tex
-rw-r--r-- 1 lena lena    520 Dec  3 00:28 report.tsv
-rw-r--r-- 1 lena lena   1010 Dec  3 00:28 report.txt
-rw-r--r-- 1 lena lena   1045 Dec  3 00:28 transposed_report.tex
-rw-r--r-- 1 lena lena    520 Dec  3 00:28 transposed_report.tsv
-rw-r--r-- 1 lena lena   995 Dec  3 00:28 transposed_report.txt
```

### Opening Report

Open `report.html`. All statistics are based on contigs of size >= 500 bp, unless otherwise noted.

#### Statistics Without Reference

| Statistic                          | Count  |
|------------------------------------|--------|
| # contigs                          | 210    |
| # contigs (>= 0 bp)               | 386    |
| # contigs (>= 1000 bp)            | 159    |
| # contigs (>= 5000 bp)            | 81     |
| # contigs (>= 10000 bp)           | 67     |
| # contigs (>= 25000 bp)           | 50     |
| # contigs (>= 50000 bp)           | 32     |
| Largest contig                     | 300763 |
| Total length                       | 5295721|
| Total length (>= 0 bp)            | 5334575|
| Total length (>= 1000 bp)         | 5259101|
| Total length (>= 5000 bp)         | 5076685|
| Total length (>= 10000 bp)        | 4977737|
| Total length (>= 25000 bp)        | 4714504|
| Total length (>= 50000 bp)        | 4035821|
| N50                                | 111860 |
| N90                                | 18506  |
| auN                                | 131921 |
| L50                                | 14     |
| L90                                | 53     |
| GC (%)                             | 50.56  |
| # N's per 100 kbp                 | 0      |
| # N's                              | 0      |

### Running SPAdes with Multiple Libraries

As mentioned, there is an advantage in using multiple libraries with different insert sizes, since the library with a small insert size can resolve short repeats, whereas the library with a larger insert size can resolve longer repeats.

This time, we will run SPAdes again by consolidating three libraries:

```bash
spades.py --pe1-1 SRR292678sub_S1_L001_R1_001.fastq.gz --pe1-2 SRR292678sub_S1_L001_R2_001.fastq.gz --mp1-1 SRR292862_S2_L001_R1_001.fastq.gz --mp1-2 SRR292862_S2_L001_R2_001.fastq.gz --mp2-1 SRR292770_S1_L001_R1_001.fastq.gz --mp2-2 SRR292770_S1_L001_R2_001.fastq.gz -o spades_output_new -t 4
```

### Download Ready-to-Use Data

Run QUAST again:

```bash
conda activate quast_env
quast.py contigs.fasta -o quast_results
```

#### Statistics Without Reference

| Statistic                          | Count  |
|------------------------------------|--------|
| # contigs                          | 105    |
| # contigs (>= 0 bp)               | 369    |
| # contigs (>= 1000 bp)            | 79     |
| # contigs (>= 5000 bp)            | 33     |
| # contigs (>= 10000 bp)           | 30     |
| # contigs (>= 25000 bp)           | 26     |
| # contigs (>= 50000 bp)           | 22     |
| Largest contig                     | 698474 |
| Total length                       | 5350156|
| Total length (>= 0 bp)            | 5403327|
| Total length (>= 1000 bp)         | 5331230|
| Total length (>= 5000 bp)         | 5202939|
| Total length (>= 10000 bp)        | 5183802|
| Total length (>= 25000 bp)        | 5133691|
| Total length (>= 50000 bp)        | 4975501|
| N50                                | 335515 |
| N90                                | 79998  |
| auN                                | 319603 |
| L50                                | 6      |
| L90                                | 20     |
| GC (%)                             | 50.59  |
| # N's per 100 kbp                 | 0      |
| # N's                              | 0      |

### Using Barrnap to Find 16S RNA Sequence

We then use Barrnap software to find the 16S rRNA sequence of the annotated strain. Barrnap predicts the location of ribosomal RNA genes in genomes.

#### Barrnap Installation

Install Barrnap via:

```bash
mamba install -c bioconda barrnap
```

#### Running Barrnap

Run Barrnap as follows:

```bash
barrnap --kingdom bac --outseq 16S_rRNA_sequences.fasta scaffolds.fasta
```

#### Resulting 16S rRNA Record Sample

```plaintext
>16S_rRNA::NODE_1_length_2815616_cov_74.3819_ID_564387:326358-327896(-)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
```

### Using BLAST to Search for the Genome

To search for the genome in the RefSeq database with 16S rRNA that is most similar to this 16S rRNA, follow these steps:

1. Open [BLAST homepage](http://blast.ncbi.nlm.nih.gov).
2. Select “Nucleotide blast”.
3. To perform the search against complete genomes in the RefSeq database, select the “Reference Genome Database (refseq_genomes)” in the “Database” field.
4. Select Escherichia coli in the “Organism” field.
5. Use `PDAT` in the "Entrez Query" field: `1900/01/01:2011/01/01[PDAT]`.

### Resulting Top Score Strain

The resulting top score strain is **Escherichia coli 55989**  
**NCBI Reference Sequence**: NC_011748.1  
[Link to Article](https://link.springer.com/article/10.1007/s00203-011-0725-6)

### Aligning with Mauve

Install Mauve:

```bash
conda install -c bioconda mauve
```

Run Mauve with:

```bash
progressiveMauve --output=alignment_output.xmfa ecoli55989.fasta scaffolds.fasta
```
```

