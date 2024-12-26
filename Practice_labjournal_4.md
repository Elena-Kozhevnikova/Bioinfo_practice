
# Tardigrades: The Resilient Water Bears

## Introduction

Tardigrades, also known as water bears, pudgy wudgies, or moss piglets, are microscopic animals capable of surviving some of the harshest environmental conditions. These eight-legged creatures can be found across diverse habitats, from the Himalayas (at elevations above 20,000 ft) to the deep sea (below 13,000 ft). 

### Extremophiles

Tardigrades are classified as "extremophiles" due to their ability to withstand:

- **Freezing temperatures** (down to 1°K)
- **Total dehydration**
- **Pressure** (over 1,200 atmospheres)
- **Radiation** (1,000 times more than other animals)

Such resilience indicates an efficient DNA repair mechanism, which was a mystery until 2016 when the genome of *Ramazzottius varieornatus*, a highly stress-tolerant tardigrade species, was sequenced.

## Project Objectives

This project aims to analyze the genome of the YOKOZUNA-1 strain of *Ramazzottius varieornatus* sequenced at the University of Tokyo, specifically focusing on understanding its unique DNA repair mechanisms.

### Genome Download

The genome can be downloaded from:
- [NCBI Genome FTP](ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz)

## Gene Prediction Using AUGUSTUS

We will employ **AUGUSTUS**, a tool for gene prediction in eukaryotic genomes, to analyze the genome. AUGUSTUS uses a probabilistic model to predict gene locations and outputs the predictions in formats like GFF (General Feature Format) and FASTA for protein sequences.

### Files Overview

| **File Name**                                     | **Description**                              |
|---------------------------------------------------|----------------------------------------------|
| GCA_001949185.1_Rvar_4.0_genomic.fna.gz          | Original genome sequence                     |
| augustus.whole.aa                                 | Predicted protein sequences                   |
| augustus.whole.gff                                | Gene structure annotations                    |

## Extracting Protein Sequences

We will extract protein sequences from the AUGUSTUS output using the following command:

```bash
perl getAnnoFasta.pl --seqfile=GCA_001949185.1_Rvar_4.0_genomic.fna.gz augustus.whole.gff
```

### Number of Obtained Proteins

To count the number of proteins:

```bash
grep ">" augustus.whole.aa | wc -l
```

**Result**: 16,435 proteins

## Hypothesis on DNA Repair

Given that DNA is a primary target of radiation damage, we hypothesize that tardigrades possess unique proteins to protect and repair their DNA. We will combine genomic and proteomic data to investigate this further.

## BLAST Analysis

### Install BLAST

```bash
mamba install -c bioconda blast
```

### Create Database

```bash
makeblastdb -in augustus.whole.aa -dbtype prot -out my_protein_db
```

### Align Mass Spectrometry Peptides

```bash
blastp -db my_protein_db -query peptides.fa -outfmt 6 -out blast_results.txt
```

### Format Results

Using the following command, we format the headers of the results:

```bash
sed -i '1i Query ID\tSubject ID\t% Identity\tAlignment Length\tMismatches\tGap Openings\tQuery Start\tQuery End\tSubject Start\tSubject End\tE-value\tBit Score' blast_results.txt
```

### Example of Blast Results

| **Query ID** | **Subject ID** | **% Identity** | **Alignment Length** | **Mismatches** | **Gap Openings** | **Query Start** | **Query End** | **Subject Start** | **Subject End** | **E-value** | **Bit Score** |
|--------------|-----------------|----------------|----------------------|----------------|-------------------|------------------|---------------|-------------------|-----------------|--------------|----------------|
| 1            | g5641.t1        | 100.000        | 9                    | 0              | 0                 | 1                | 9             | 25                | 33              | 2.0          | 21.6           |

## Protein IDs of Interest

We extracted the unique protein IDs of interest using the command:

```bash
awk 'NR > 1 {print $2}' blast_results.txt | sort -u > unique_query_ids.txt
```

### Resulting Protein IDs

| **Protein ID** |
|-----------------|
| g10513.t1       |
| g10514.t1       |
| g11320.t1       |
| ...             |
| g16368.t1       |

## Subcellular Localization Prediction

Using **WoLF PSORT** (https://wolfpsort.hgc.jp/), we will predict the subcellular localization of proteins. Warning messages indicate proteins that do not start with Methionine (M).

### Predicted Localizations

| **Gene ID** | **Details**                                         |
|--------------|-----------------------------------------------------|
| g2203.t1     | plas: 29, nucl: 2, golg: 1                          |
| g3428.t1     | mito: 18, cyto: 11, extr: 2, nucl: 1               |
| g5616.t1     | extr: 31, mito: 1                                   |
| g5927.t1     | nucl: 30.5, cyto_nucl: 16.5, cyto: 1.5             |
| ...          | ...                                                 |

## Selection of Nuclear Localization Proteins

We will select proteins with nuclear localization and create a file `wolf_preselected_ids.txt`:

### Selected Protein IDs

| **Selected Protein ID** |
|--------------------------|
| g2203.t1                 |
| g3428.t1                 |
| ...                      |

## Amino Acid Sequence Extraction

To recover amino acid sequences for these IDs:

```bash
seqtk subseq augustus.whole.aa wolf_preselected_ids.txt > wolf_proteins.fasta
```

## Domain Prediction

We will submit to the domain prediction service available at [DTU Health Tech](https://services.healthtech.dtu.dk/).

### TargetP Predictions

| **ID**      | **Prediction** | **SP**    | **mTP**   | **CS Position**             |
|-------------|----------------|-----------|-----------|-----------------------------|
| g702.t1     | SP             | 0.000347  | 0.000001  | CS pos: 16-17, ALA-AN       |
| g1285.t1    | SP             | 0.003029  | 0.000173  | CS pos: 16-17, ASA-TS       |
| g2203.t1    | OTHER          | 0.999869  | 0.000100  |                             |
| g3428.t1    | OTHER          | 0.999903  | 0.000064  |                             |
| ...         | ...            | ...       | ...       |                             |

## Pfam Domain Predictions

We will download the Pfam database and perform domain predictions.

### Unique Gene Predictions from HMMER

| **Target Name**      | **Accession**   | **Query Name** | **E-value** | **Score** | **Description**                                |
|----------------------|------------------|-----------------|-------------|-----------|------------------------------------------------|
| CBM_14               | PF01607.23       | g702.t1         | 3.6e-10     | 39.7      | Chitin binding Peritrophin-A domain            |
| Glyco_hydro_31       | PF01055.25       | g2203.t1        | 3.5e-42     | 145.0     | Glycosyl hydrolases family 31                  |
| Roc                  | PF08477.12       | g2203.t1        | 0.098       | 12.8      | Ras of Complex, Roc, domain of DAPkinase      |
| EF-hand_7            | PF13499.5        | g3428.t1        | 1.1e-14     | 54.5      | EF-hand domain pair                            |
| EF-hand_1            | PF00036.31       | g3428.t1        | 3.1e-14     | 51.3      | EF hand                                         |
| ...                  | ...              | ...             | ...         | ...       | ...                                            |

## Conclusion

This project aims to explore the unique adaptations of tardigrades through genomic and proteomic analyses. Their exceptional resilience to extreme conditions not only provides insights into their biology but may also have implications for biotechnology and medicine.
```

