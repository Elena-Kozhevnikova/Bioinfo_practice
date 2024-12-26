# Tardigrades: The Resilient Water Bears

## Introduction

Tardigrades, also known as water bears, pudgy wudgies, or moss piglets, are microscopic animals capable of withstanding some of the most severe environmental conditions. These water-dwelling, eight-legged creatures can be found throughout the world, from the Himalayas (above 20,000 ft) to the deep sea (below 13,000 ft). Classified as "extremophiles", they can survive freezing (up to 1°K), total dehydration, pressure (more than 1,200 atmospheres), and radiation (1,000 times more radiation than other animals).

### DNA Repair Mechanism

It means that they can efficiently repair damage to their DNA. But how?! It was a mystery until 2016, when we were able to sequence the genome of *Ramazzottius varieornatus*, one of the most stress-tolerant species of tardigrades. In this project, we will analyze this genome and understand this secret.

## Project Overview

For this project, we will use a sequence of the *Ramazzottius varieornatus*, the YOKOZUNA-1 strain (sequenced in the University of Tokyo and named after the highest rank in professional sumo).

### Genome Download

The genome was downloaded from the following link:
- [NCBI Genome FTP](ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz)

## Gene Prediction Using AUGUSTUS

To look for homologous proteins and conserved domains to assign potential functions to the genes and proteins, we can use **AUGUSTUS**. AUGUSTUS is a widely used software program for gene prediction in eukaryotic genomes. It employs a probabilistic model to predict the locations of genes and the structures of their transcripts, including coding regions and untranslated regions. 

### AUGUSTUS Features

- **Gene Annotations**: Can take advantage of known gene annotations.
- **Training**: Can be trained on specific organisms to improve its accuracy.
- **Output Formats**: Outputs predictions in formats like GFF (General Feature Format) and FASTA format for the predicted protein sequences.

## Summary of Files Obtained

After running AUGUSTUS, we ended up with the following four files:

| **File Name**                                     | **Description**                           |
|---------------------------------------------------|-------------------------------------------|
| GCA_001949185.1_Rvar_4.0_genomic.fna.gz          | Genome sequence                           |
| augustus.whole.aa                                 | Predicted protein sequences                |
| augustus.whole.gff                                | Gene annotations                          |
| augustus.whole.aa (processed)                     | Extracted protein sequences from GFF     |

### Extracting Protein Sequences

To extract protein sequences, we use the `getAnnoFasta.pl` script:

```bash
perl getAnnoFasta.pl --seqfile=GCA_001949185.1_Rvar_4.0_genomic.fna.gz augustus.whole.gff
```

#### Resulting Protein File
The resulting protein file is `augustus.whole.aa`. The head of the file looks like this:

```
>g1.t1
AATYSKHPRSYRENIVVWRICMGHCDGRYLALYNLNGWSERDMIDRWALMVTLLVPPAFLLPARDIDNLLTAGHETSPQTGGSDVMSEYVTPFCYVCDSL
VDADCRESQDWMPRRRNGLRKKVMWTGQGEDDHGRLAFDWKEVYIRLEKFKQRCPSNNGCTTHVFLNSDRVIRRCFDPKEPPCSGHEWMVNRRNPKKRSK
VKVVSRLTCADWPFCNVDSAYYWPNGTVEALLSPHDDPQLKVRSAKLSCYECDSILQPDCTGQNFVKAVRKPNDAVAHSVDWVGLRKGMEKFRTRCPVES
GCFFQGDHMTNRVIRGCRNMTTDLPRHLLVKDRRKESSVFRTEEVCRESECNSRVLFFSHSGAFLDSDWHFNNSSAGRDNGQGEAKKEVYQGDAPATEET
GFLPHVKQSGQGKRKEGRLQKRTIGDVEREERRADEARQRRGKPKAEGRGGAGTRVLSFGHRSRKTSAGGKLDSCWLLLLLPFFLL
>g2.t1
MTFAVTVATTSQSCAIAFDVFWLTTIKLFFIAKYFGSEFTPFDKVQSLTCPHCAQQGIEYKDLLPHLETNHADTSDFLVVSDNITIIWDVAHCI
>g3.t1
MLDKAKLVLLNKGGRPSIRSAEDGWKVFSDSLLCRASSDGRVFIGKGGGNGGAVAKAEVAAIGVVLVHHMLQQVVRHGEVVIVSANTNWTPLPY
```

### Counting the Number of Proteins

To count the number of obtained proteins:

```bash
grep ">" augustus.whole.aa | wc -l
```

**Result**: 16,435

### Hypothesis on Tardigrades

Considering that DNA is a major target of radiation damage, we hypothesize that tardigrades might have unique proteins associated with their DNA to protect and/or effectively repair it. To explore this possibility, we will combine genomic and proteomic data.

## BLAST Analysis

### Installing BLAST

```bash
mamba install -c bioconda blast
```

### Creating a BLAST Database

Create the database for BLAST:

```bash
makeblastdb -in augustus.whole.aa -dbtype prot -out my_protein_db
```

### Aligning Mass Spectrometry Peptides

To align mass spectrometry peptides to the newly created database:

```bash
blastp -db my_protein_db -query peptides.fa -outfmt 6 -out blast_results.txt
```

### Formatting the Results

We’ll format the output to add headers:

```bash
sed -i '1i Query ID\tSubject ID\t% Identity\tAlignment Length\tMismatches\tGap Openings\tQuery Start\tQuery End\tSubject Start\tSubject End\tE-value\tBit Score' blast_results.txt
```

### Example of BLAST Results

| **Query ID** | **Subject ID** | **% Identity** | **Alignment Length** | **Mismatches** | **Gap Openings** | **Query Start** | **Query End** | **Subject Start** | **Subject End** | **E-value** | **Bit Score** |
|--------------|-----------------|----------------|----------------------|----------------|-------------------|------------------|---------------|-------------------|-----------------|--------------|----------------|
| 1            | g5641.t1        | 100.000        | 9                    | 0              | 0                 | 1                | 9             | 25                | 33              | 2.0          | 21.6           |
| 1            | g15153.t1       | 100.000        | 9                    | 0              | 0                 | 1                | 9             | 25                | 33              | 2.1          | 21.6           |
| 2            | g5641.t1        | 100.000        | 9                    | 0              | 0                 | 1                | 9             | 36                | 44              | 1.5          | 21.9           |
| 2            | g12562.t1       | 88.889         | 9                    | 1              | 0                 | 1                | 9             | 36                | 44              | 3.3          | 20.8           |
| 2            | g5616.t1        | 88.889         | 9                    | 1              | 0                 | 1                | 9             | 36                | 44              | 6.7          | 20.0           |
| 2            | g15153.t1       | 77.778         | 9                    | 2              | 0                 | 1                | 9             | 36                | 44              | 7.5          | 20.0           |
| 3            | g5641.t1        | 100.000        | 13                   | 0              | 0                 | 1                | 13            | 155               | 167             | 0.26         | 24.6           |
| 4            | g5641.t1        | 100.000        | 11                   | 0              | 0                 | 1                | 11            | 144               | 154             | 0.053        | 26.2           |
| 4            | g12562.t1       | 72.727         | 11                   | 3              | 0                 | 1                | 11            | 144               | 154             | 1.6          | 21.9           |

### Explanation of Results

- **E-value**: Indicates the number of hits one can "expect" to see by chance when searching a database of a particular size. Lower E-values indicate more significant matches.
- **Bit Score**: Represents the score of the alignment based on a specific scoring system. Higher scores generally indicate better matches.

## Extracting Unique Protein IDs

### Installation of seqtk

```bash
mamba install -c bioconda seqtk
```

Make sure that you have the bioconda channel added to your conda configuration, as seqtk is hosted there:

```bash
conda config --add channels bioconda
```

### Get Unique Protein IDs

Get the IDs of the proteins of interest:

```bash
awk 'NR > 1 {print $2}' blast_results.txt | sort -u > unique_query_ids.txt
```

### Resulting Protein IDs

The resulting protein IDs are:

| **Protein ID** |
|-----------------|
| g10513.t1       |
| g10514.t1       |
| g11320.t1       |
| g11513.t1       |
| g11806.t1       |
| g11960.t1       |
| g12388.t1       |
| g12510.t1       |
| g12562.t1       |
| g1285.t1        |
| g13530.t1      |
| g14472.t1      |
| g15153.t1      |
| g15484.t1      |
| g16318.t1      |
| g16368.t1      |
| g2203.t1       |
| g3428.t1       |
| g3679.t1       |
| g4106.t1       |
| g4970.t1       |
| g5237.t1       |
| g5443.t1       |
| g5467.t1       |
| g5502.t1       |
| g5503.t1       |
| g5510.t1       |
| g5616.t1       |
| g5641.t1       |
| g5927.t1       |
| g702.t1        |
| g7861.t1       |
| g8100.t1       |
| g8312.t1       |

### Extracting Amino Acid Sequences

To extract their amino acid sequences with seqtk:

```bash
seqtk subseq augustus.whole.aa unique_query_ids.txt > extracted_unique_query_proteins.fasta
```

## Predicting Subcellular Localization

Using **WoLF PSORT** (https://wolfpsort.hgc.jp/), we can predict the subcellular localization of proteins.

### Predictions and Warnings

| **Gene ID** | **Localization Details**                                            |
|-------------|---------------------------------------------------------------------|
| g16318.t1   | Warning: sequence does not start with M                            |
| g1285.t1    | extr: 25, plas: 5, mito: 1, lyso: 1                                 |
| g2203.t1    | plas: 29, nucl: 2, golg: 1                                         |
| g3428.t1    | mito: 18, cyto: 11, extr: 2, nucl: 1                               |
| g3679.t1    | extr: 26, mito: 2, lyso: 2, plas: 1, E.R.: 1                       |
| g4106.t1    | E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1 |
| g4970.t1    | plas: 32                                                           |
| g5237.t1    | plas: 24, mito: 8                                                  |
| g5443.t1    | extr: 28, nucl: 3, cyto: 1                                         |
| g5467.t1    | extr: 27, plas: 4, mito: 1                                         |
| g5502.t1    | extr: 31, lyso: 1                                                  |
| g5503.t1    | extr: 29, plas: 1, mito: 1, lyso: 1                                |
| g5510.t1    | plas: 23, mito: 7, E.R.: 1, golg: 1                                |
| g5616.t1    | extr: 31, mito: 1                                                 |
| g5641.t1    | extr: 31, lyso: 1                                                 |
| g5927.t1    | nucl: 30.5, cyto_nucl: 16.5, cyto: 1.5                            |
| g7861.t1    | nucl: 16, cyto_nucl: 14, cyto: 8, plas: 5, pero: 1, cysk: 1, golg: 1 |
| g8100.t1    | nucl: 16.5, cyto_nucl: 12.5, cyto: 7.5, plas: 5, extr: 2, E.R.: 1 |
| g8312.t1    | nucl: 15.5, cyto_nucl: 15.5, cyto: 12.5, mito: 2, plas: 1, golg: 1 |
| g10513.t1   | nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1 |
| g10514.t1   | nucl: 19, cyto_nucl: 15, cyto: 9, extr: 3, mito: 1               |
| g11320.t1   | plas: 24.5, extr_plas: 16, extr: 6.5, lyso: 1                     |
| g11513.t1   | cyto: 17, cyto_nucl: 12.8333, cyto_mito: 9.83333, nucl: 7.5, E.R.: 3, mito: 1.5, plas: 1, pero: 1, golg: 1 |
| g11806.t1   | nucl: 18, cyto_nucl: 11.8333, mito: 5, extr: 4, cyto: 3.5, cyto_pero: 2.66667, cysk_plas: 1 |
| g11960.t1   | nucl: 32                                                          |
| g12388.t1   | extr: 25, plas: 4, mito: 2, lyso: 1                                |
| g12510.t1   | plas: 29, cyto: 3                                                 |
| g12562.t1   | extr: 30, lyso: 2                                                 |
| g13530.t1   | extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5 |
| g14472.t1   | nucl: 28, plas: 2, cyto: 1, cysk: 1                               |
| g15153.t1   | extr: 32                                                          |
| g15484.t1   | nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1 |
| g16318.t1   | nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1 |
| g16368.t1   | nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1 |

### Selecting Nuclear Localization Proteins

We shall select those with nuclear localization and create a `wolf_preselected_ids.txt` file with the protein IDs:

### Selected Protein IDs

| **Protein ID** |
|-----------------|
| g2203.t1        |
| g3428.t1        |
| g5443.t1        |
| g5927.t1        |
| g7861.t1        |
| g8100.t1        |
| g8312.t1        |
| g10513.t1       |
| g10514.t1       |
| g11513.t1       |
| g11806.t1       |
| g11960.t1       |
| g13530.t1       |
| g14472.t1       |
| g15484.t1       |
| g16318.t1       |
| g16368.t1       |

### Extracting Amino Acid Sequences for Selected IDs

To recover amino acid sequences for these IDs:

```bash
seqtk subseq augustus.whole.aa wolf_preselected_ids.txt > wolf_proteins.fasta
```

### Domain Prediction

Next, we submit the data to the domain prediction service available at DTU Health Tech.

## TargetP Predictions

### Results of Domain Prediction
SP	- is a numerical score that indicates the probability of the protein containing a signal peptide. Values close to 1 suggest strong evidence for a signal peptide, implying it may be directed towards the secretory pathway.

| **ID**      | **Prediction** | **OTHER** | **SP**    | **mTP**   | **CS Position**                     |
|-------------|----------------|-----------|-----------|-----------|-------------------------------------|
| g702.t1     | SP             | 0.000347  | 0.999652  | 0.000001  | CS pos: 16-17. ALA-AN               |
| g1285.t1    | SP             | 0.003029  | 0.996798  | 0.000173  | CS pos: 16-17. ASA-TS               |
| g2203.t1    | OTHER          | 0.999869  | 0.000031  | 0.000100  |                                     |
| g3428.t1    | OTHER          | 0.999903  | 0.000033  | 0.000064  |                                     |
| g3679.t1    | SP             | 0.001755  | 0.998023  | 0.000222  | CS pos: 18-19. TFA-AR               |
| g4106.t1    | OTHER          | 0.729658  | 0.266917  | 0.003425  |                                     |
| g4970.t1    | OTHER          | 0.999996  | 0.000003  | 0.000001  |                                     |
| g5237.t1    | OTHER          | 0.999545  | 0.000345  | 0.000111  |                                     |
| g5443.t1    | OTHER          | 0.952853  | 0.043784  | 0.003363  |                                     |
| g5467.t1    | SP             | 0.000096  | 0.999845  | 0.000059  | CS pos: 16-17. ASA-GS               |
| g5502.t1    | SP             | 0.001134  | 0.998823  | 0.000043  | CS pos: 16-17. ASA-GS               |
| g5503.t1    | SP             | 0.001222  | 0.998720  | 0.000058  | CS pos: 16-17. ASA-GS               |
| g5510.t1    | OTHER          | 0.999108  | 0.000016  | 0.000876  |                                     |
| g5616.t1    | SP             | 0.000067  | 0.999933  | 0.000000  | CS pos: 16-17. ACA-AN               |
| g5641.t1    | SP             | 0.000130  | 0.999869  | 0.000001  | CS pos: 16-17. ACA-AS               |
| g5927.t1    | OTHER          | 0.999995  | 0.000001  | 0.000004  |                                     |
| g7861.t1    | OTHER          | 0.999975  | 0.000004  | 0.000022  |                                     |
| g8100.t1    | OTHER          | 0.999955  | 0.000024  | 0.000021  |                                     |
| g8312.t1    | OTHER          | 0.999930  | 0.000065  | 0.000004  |                                     |
| g10513.t1   | OTHER          | 0.999999  | 0.000001  | 0.000000  |                                     |
| g10514.t1   | OTHER          | 0.999543  | 0.000349  | 0.000107  |                                     |
| g11320.t1   | SP             | 0.000184  | 0.999816  | 0.000000  | CS pos: 20-21. AYS-AG               |
| g11513.t1   | OTHER          | 0.999434  | 0.000401  | 0.000164  |                                     |
| g11806.t1   | OTHER          | 0.998977  | 0.000887  | 0.000136  |                                     |
| g11960.t1   | OTHER          | 0.999996  | 0.000002  | 0.000002  |                                     |
| g12388.t1   | SP             | 0.000490  | 0.999481  | 0.000029  | CS pos: 16-17. ASA-SS               |
| g12510.t1   | OTHER          | 0.999738  | 0.000099  | 0.000163  |                                     |
| g12562.t1   | SP             | 0.000076  | 0.999923  | 0.000001  | CS pos: 16-17. SYA-AN               |
| g13530.t1   | SP             | 0.116007  | 0.883840  | 0.000153  | CS pos: 19-20. TIP-FT               |
| g14472.t1   | OTHER          | 0.999999  | 0.000001  | 0.000000  |                                     |
| g15153.t1   | SP             | 0.000014  | 0.999986  | 0.000000  | CS pos: 16-17. AYA-AN               |
| g15484.t1   | OTHER          | 0.999980  | 0.000010  | 0.000010  |                                     |
| g16318.t1   | OTHER          | 0.997047  | 0.002953  | 0.000000  |                                     |
| g16368.t1   | OTHER          | 0.996693  | 0.003307  | 0.000000  |                                     |

### Downloading and Scanning Pfam

To perform additional domain prediction, download the Pfam database:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm
hmmpress Pfam-A.hmm
hmmscan --domtblout output.domtbl Pfam-A.hmm ~/bio/project4/extracted_unique_query_proteins.fasta
```

### Unique Genes from HMMER

| **Target Name** | **Accession**   | **Query Name** | **E-value** | **Score** | **Description**                              |
|------------------|-----------------|-----------------|-------------|-----------|----------------------------------------------|
| CBM_14           | PF01607.23      | g702.t1         | 3.6e-10     | 39.7      | Chitin binding Peritrophin-A domain         |
| Glyco_hydro_31   | PF01055.25      | g2203.t1        | 3.5e-42     | 145.0     | Glycosyl hydrolases family 31                |
| Roc               | PF08477.12      | g2203.t1        | 0.098       | 12.8      | Ras of Complex, Roc, domain of DAPkinase    |
| EF-hand_7        | PF13499.5       | g3428.t1        | 1.1e-14     | 54.5      | EF-hand domain pair                          |
| EF-hand_1        | PF00036.31      | g3428.t1        | 3.1e-14     | 51.3      | EF hand                                      |
| UPF0154          | PF03672.12      | g3428.t1        | 0.00016     | 21.5      | Uncharacterised protein family (UPF0154)    |
| SurA_N_2         | PF13623.5       | g3428.t1        | 0.002       | 17.9      | SurA N-terminal domain                       |
| Toxin_TOLIP      | PF00087.20      | g5237.t1        | 0.22        | 11.7      | Snake toxin and toxin-like protein           |
| CBM_14           | PF01607.23      | g5467.t1        | 1e-11       | 44.7      | Chitin binding Peritrophin-A domain         |
| CBM_14           | PF01607.23      | g5502.t1        | 3.1e-11     | 43.2      | Chitin binding Peritrophin-A domain         |
| Ldl_recept_a     | PF00057.17      | g4970.t1        | 1.5e-15     | 57.0      | Low-density lipoprotein receptor domain class A |
| zf-CCHC          | PF00098.22      | g4970.t1        | 5.6e-14     | 51.4      | Zinc knuckle                                 |
| CBM_14           | PF01607.23      | g5616.t1        | 1.5e-11     | 44.2      | Chitin binding Peritrophin-A domain         |
| TF_AP-2         | PF03299.13      | g15484.t1       | 0.54        | 10.1      | Transcription factor AP-2                   |
| Mucin            | PF01456.16      | g15484.t1       | 5.2         | 7.0       | Mucin-like glycoprotein                     |


### After bringing the data data from all three tables together, the final result is this:

| **Gene ID**    | **Query ID / E-value / Bit Score**                                       | **Nuclear Localization** | **Domain / E-value / Description**                                     |
|----------------|-------------------------------------------------------------------------|-------------------------|-----------------------------------------------------------------------|
| g2203.t1       | Matches: 3 (g2203.t1, E-value: 3.5e-42, Bit Score: 145.0)               | Nuclear                 | Glyco_hydro_31 / PF01055.25 / Glycosyl hydrolases family 31          |
|                |                                                                         |                         | Roc / PF08477.12 / Ras of Complex, Roc, domain of DAPkinase          |
| g3428.t1       | Matches: 5 (g3428.t1, E-value: 1.1e-14, Bit Score: 54.5)               | Nuclear                 | EF-hand_7 / PF13499.5 / EF-hand domain pair                           |
|                |                                                                         |                         | EF-hand_1 / PF00036.31 / EF hand                                     |
|                |                                                                         |                         | UPF0154 / PF03672.12 / Uncharacterised protein family (UPF0154)      |
|                |                                                                         |                         | SurA_N_2 / PF13623.5 / SurA N-terminal domain                         |
| g5237.t1       | Matches: 1 (g5237.t1, E-value: 0.004, Bit Score: 29.6)                  | Nuclear                 | Toxin_TOLIP / PF00087.20 / Snake toxin and toxin-like protein         |
| g5467.t1       | Matches: 2 (g5467.t1, E-value: 1e-11, Bit Score: 44.7)                  | Nuclear                 | CBM_14 / PF01607.23 / Chitin binding Peritrophin-A domain            |
| g5502.t1       | Matches: 3 (g5502.t1, E-value: 3.1e-11, Bit Score: 43.2)                | Nuclear                 | CBM_14 / PF01607.23 / Chitin binding Peritrophin-A domain            |
|                |                                                                         |                         | Ldl_recept_a / PF00057.17 / Low-density lipoprotein receptor domain class A |
| g5616.t1       | Matches: 5 (g5616.t1, E-value: 1.5e-11, Bit Score: 44.2)                | Nuclear                 | CBM_14 / PF01607.23 / Chitin binding Peritrophin-A domain            |
| g10513.t1      | Matches: 3 (g10513.t1, E-value: 0.003, Bit Score: 30.0)                 | Nuclear                 | -                                                                     |
| g10514.t1      | Matches: 3 (g10514.t1, E-value: 0.97, Bit Score: 22.7)                  | Nuclear                 | -                                                                     |
| g11513.t1      | Matches: 2 (g11513.t1, E-value: 8.3, Bit Score: 20.0)                   | Nuclear                 | -                                                                     |
| g11806.t1      | Matches: 1 (g11806.t1, E-value: 3.9, Bit Score: 21.2)                   | Nuclear                 | -                                                                     |
| g11960.t1      | Matches: 1 (g11960.t1, E-value: 3.5, Bit Score: 21.6)                   | Nuclear                 | -                                                                     |
| g13530.t1      | Matches: 2 (g13530.t1, E-value: 0.002, Bit Score: 30.8)                 | Nuclear                 | -                                                                     |
| g14472.t1      | Matches: 1 (g14472.t1, E-value: 0.002, Bit Score: 30.8)                 | Nuclear                 | -                                                                     |
| g15484.t1      | Matches: 3 (g15484.t1, E-value: 0.54, Bit Score: 10.1)                  | Nuclear                 | TF_AP-2 / PF03299.13 / Transcription factor AP-2                      |
|                |                                                                         |                         | Mucin / PF01456.16 / Mucin-like glycoprotein                         |
| g16318.t1      | Matches: 1 (g16318.t1, E-value: 4.9, Bit Score: 21.2)                   | Nuclear                 | -                                                                     |
| g16368.t1      | Matches: 1 (g16368.t1, E-value: 5.2, Bit Score: 20.8)                   | Nuclear                 | -                                                                     |

### From these we only select proteins with a potential signal peptide:

| **Gene ID**    | **Query ID / E-value / Bit Score**                                    | **Nuclear Localization** | **Domain / E-value / Description**                                     |
|----------------|----------------------------------------------------------------------|-------------------------|-----------------------------------------------------------------------|
| g5467.t1       | Matches: 2 (g5467.t1, E-value: 1e-11, Bit Score: 44.7)               | Nuclear                 | CBM_14 / PF01607.23 / Chitin binding Peritrophin-A domain            |
| g5502.t1       | Matches: 3 (g5502.t1, E-value: 3.1e-11, Bit Score: 43.2)             | Nuclear                 | CBM_14 / PF01607.23 / Chitin binding Peritrophin-A domain            |
|                |                                                                      |                         | Ldl_recept_a / PF00057.17 / Low-density lipoprotein receptor domain class A |
| g5616.t1       | Matches: 5 (g5616.t1, E-value: 1.5e-11, Bit Score: 44.2)              | Nuclear                 | CBM_14 / PF01607.23 / Chitin binding Peritrophin-A domain            |
| g13530.t1      | Matches: 2 (g13530.t1, E-value: 0.002, Bit Score: 30.8)              | Nuclear                 | -                                                                     |

### The interesting candidates are g13530.t1 and g3428.t1

The first one is a completely unknown protein. Let's extract its amino acid sequence:

```bash
awk '/^>g13530\.t1/{flag=1; next} /^>/{flag=0} flag' augustus.whole.aa
```
We then align it in blast, and retrieve no putative domains:
![image](https://github.com/user-attachments/assets/77cf162b-6ef1-4ff8-9daf-c55b45cf9765)

It partially aligng to another Tardigrade Hypsibius exemplaris, whic is either a philogenetic relation, or possibly a common adaptive feature.

This protein might be the one with the unknown function and does not align with known proteins or their domains. So that it might serve a candidate to search for a completely novel adaptive mechanism observed in Tardigrades. 

## Conclusion

This project aims to explore the unique adaptations of tardigrades through genomic and proteomic analyses. Their exceptional resilience to extreme conditions not only provides insights into their biology but may also have implications for biotechnology and medicine.
```

