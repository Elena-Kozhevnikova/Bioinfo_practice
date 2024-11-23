```markdown
# Viral Data Analysis

## Data Source
Data was published there and labeled **SRR1705851**, so you can download it from the SRA FTP server and unpack it on your machine:

[Download Data](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/)

Download the reference gene from here:

[Download Reference Gene](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta)

## Unzipping the FASTQ File
Unzip the fastq file:
```bash
gunzip -c SRR1705851.fastq.gz > roommate.fastq
```

## Indexing and Aligning
Index the reference file and align your roommate’s viral data to the reference sequence. 

### Generate Index Files for `HA_gene.fasta`
```bash
bwa index HA_gene.fasta
```

This results in the following files:
- HA_gene.fasta.bwt
- HA_gene.fasta.pac
- HA_gene.fasta.ann
- HA_gene.fasta.sa
- HA_gene.fasta.amb
- HA_gene.fasta.bai

### Aligning the Data
```bash
bwa mem HA_gene.fasta SRR1705851.fastq > output.sam
```

**Convert SAM to BAM**
```bash
samtools view -S -b output.sam > output.bam
```

**Sort the BAM file**
```bash
samtools sort -o HA_result.bam output.bam
```

**Index the BAM file**
```bash
samtools index HA_result.bam
```
This results in creating `HA_result.bam.bai`.

## Create an MPileup
Using the `samtools mpileup` command, set the depth limit to 20,000:
```bash
samtools mpileup -d 20000 -o output.mpileup -f HA_gene.fasta HA_result.bam
```
The resulting file is `output.mpileup`.

## Count Unmapped Reads
```bash
samtools view -f 4 -c HA_result.bam > unmapped_read_count.txt
```
**Result:** 
```
233
```

## SNP Mapping with VarScan
Use VarScan to map the SNPs:
```bash
java -jar VarScan.v2.4.6.jar mpileup2snp output.mpileup --min-var-freq 0.95 --output-var 1 --min-coverage 1 --p-value 0.01 > variants.snp 
```

### Extract High-Quality Variants
```bash
awk '$1 !~ /^#/ && $6 > 30' variants.snp > filtered_variants.txt
```

### Output High-Quality Variants
Results from `awk` for filtering:
```bash
cat filename.vcf | awk 'NR>24 {print $1, $2, $4, $5}'
```
**Results:**
```
KF848938.1 72 A G
KF848938.1 117 C T
KF848938.1 774 T C
KF848938.1 999 C T
KF848938.1 1260 A C
```

## Mapping Mutations to the Gene
### Silent Mutations:
- **Position 72:** ACA -> ACG
- **Position 117:** GCC -> GCT
- **Position 774:** TTT -> TTC
- **Position 999:** GGC -> GGT
- **Position 1260:** CTA -> GCT

*All five mutations are silent, so they cannot explain how the virus escaped the immune system of a patient. Thus we are to search for rare variants.*

## Rare Variants Analysis
```bash
java -jar VarScan.v2.4.6.jar mpileup2snp output.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_rare.vcf
```

**Output Results:**
- 1665 bases in pileup file
- 18 variant positions (16 SNP, 2 indel)
- 0 were failed by the strand-filter
- 16 variant positions reported (16 SNP, 0 indel)

### Frequencies Extraction
Use the following command to get the frequencies:
```bash
cat VarScan_results_rare.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}'
```

**Extracted Columns (without frequencies):**
```
KF848938.1 72 A G
KF848938.1 117 C T
KF848938.1 165 T C
KF848938.1 307 C T
KF848938.1 389 T C
KF848938.1 722 A G
KF848938.1 774 T C
KF848938.1 802 A G
KF848938.1 913 T C
KF848938.1 915 T C
KF848938.1 999 C T
KF848938.1 1086 A G
KF848938.1 1213 A G
KF848938.1 1260 A C
KF848938.1 1280 T C
KF848938.1 1458 T C
```

### Checking Results with SnapGene Viewer
Here’s a summary of some variants:
| Position | Original Codon | Changed Codon | Mutation Type          |
|----------|----------------|----------------|-------------------------|
| 307      | CCG            | TCG            | aa substitution, 0.91%  |
| 389      | GTC            | GCC            | aa substitution, 0.22%  |
| 722      | GAC            | GGC            | aa substitution, 0.22%  |
| 802      | ATG            | GTG            | aa substitution, 0.27%  |
| 913      | TGT            | CGT            | aa substitution, 0.21%  |
| 1213     | AGA            | GGA            | aa substitution, 0.27%  |
| 1280     | CTT            | CCT            | aa substitution, 0.2%   |

## Comparison with Control Sequences
Download the Illumina sequences from these plasmids:
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz
```

### Count Read Numbers
```bash
echo "SRR1705858 has $(( $(zcat SRR1705858.fastq.gz | wc -l) / 4 )) reads"
echo "SRR1705859 has $(( $(zcat SRR1705859.fastq.gz | wc -l) / 4 )) reads"
echo "SRR1705860 has $(( $(zcat SRR1705860.fastq.gz | wc -l) / 4 )) reads"
```

**Results:**
```
SRR1705858 has 256586 reads
SRR1705859 has 233327 reads
SRR1705860 has 249964 reads
```

### Coverage Calculation
```bash
samtools faidx HA_gene.fasta
cat HA_gene.fasta.fai
```
**Result:**
```
KF848938.1      1665    104     70      71
```

**Estimate Coverage:**
Coverage can be estimated as:
\[
\text{Coverage} = \frac{\text{# of reads} \times 70}{1665}
\]

| Sample        | Coverage Estimate |
|---------------|-------------------|
| SRR1705858    | 10,800            |
| SRR1705859    | 9,800             |
| SRR1705860    | 10,500            |

### Align Control FASTQ Files
Align each control fastq file to the reference:
```bash
bwa mem HA_gene.fasta SRR1705858.fastq > control_sample1_aln.sam
bwa mem HA_gene.fasta SRR1705859.fastq > control_sample2_aln.sam
bwa mem HA_gene.fasta SRR1705860.fastq > control_sample3_aln.sam
```

### Convert and Sort
```bash
samtools view -Sb control_sample1_aln.sam | samtools sort -o control_sample1_sorted.bam
samtools view -Sb control_sample2_aln.sam | samtools sort -o control_sample2_sorted.bam
samtools view -Sb control_sample3_aln.sam | samtools sort -o control_sample3_sorted.bam
```

### Run VarScan on Each Alignment
1. Generate a pileup file for each sorted BAM file:
```bash
samtools mpileup -f HA_gene.fasta control_sample1_sorted.bam > control_sample1.pileup
samtools mpileup -f HA_gene.fasta control_sample2_sorted.bam > control_sample2.pileup
samtools mpileup -f HA_gene.fasta control_sample3_sorted.bam > control_sample3.pileup
```

2. Run VarScan on each pileup file:
```bash
java -jar VarScan.v2.4.6.jar mpileup2cns control_sample1.pileup --min-var-freq 0.001 --output-vcf 1 > control_sample1_variants.vcf
java -jar VarScan.v2.4.6.jar mpileup2cns control_sample2.pileup --min-var-freq 0.001 --output-vcf 1 > control_sample2_variants.vcf
java -jar VarScan.v2.4.6.jar mpileup2cns control_sample3.pileup --min-var-freq 0.001 --output-vcf 1 > control_sample3_variants.vcf
```

### Output Results
**Control Sample 1:**
```
1665 bases in pileup file
28 variant positions (28 SNP, 0 indel)
1 were failed by the strand-filter
28 variant positions reported (28 SNP, 0 indel)
```

**Control Sample 2:**
```
1665 bases in pileup file
24 variant positions (24 SNP, 0 indel)
1 were failed by the strand-filter
24 variant positions reported (24 SNP, 0 indel)
```

**Control Sample 3:**
```
1665 bases in pileup file
27 variant positions (27 SNP, 0 indel)
1 were failed by the strand-filter
27 variant positions reported (27 SNP, 0 indel)
```

### Parse Each VCF File
Create three commands to extract the required fields and save them to separate files.
```bash
awk '$1 !~ /^#/ && $5 != "." { print $2, $3, $4, $5, $10 }' control_sample1_variants.vcf > filtered_variants_control1.csv
```
<table>
    <tr>
        <td>

### Control File 1 Summary
| Position | Original | Alternative | Frequency |
|----------|----------|-------------|-----------|
| 125      | T        | C           | 0.43%     |
| 129      | T        | C           | 0.41%     |
| 139      | T        | C           | 0.31%     |
| 151      | A        | G           | 0.28%     |
| 183      | A        | G           | 0.33%     |
| 276      | A        | G           | 0.30%     |
| 370      | A        | G           | 0.27%     |
| 389      | T        | C           | 0.26%     |
| 409      | T        | C           | 0.30%     |
| 566      | A        | G           | 0.32%     |
| 595      | G        | T           | 0.29%     |
| 722      | A        | G           | 0.30%     |
| 744      | A        | G           | 0.28%     |
| 774      | T        | C           | 0.33%     |
| 859      | A        | G           | 0.26%     |
| 898      | A        | G           | 0.26%     |
| 915      | T        | C           | 0.41%     |
| 1086     | A        | G           | 0.28%     |
| 1209     | A        | G           | 0.34%     |
| 1213     | A        | G           | 0.31%     |
| 1260     | A        | C           | 0.30%     |
| 1264     | T        | C           | 0.26%     |
| 1280     | T        | C           | 0.26%     |
| 1339     | T        | C           | 0.47%     |
| 1358     | A        | G           | 0.41%     |
| 1366     | A        | G           | 0.29%     |
| 1460     | A        | G           | 0.32%     |
| 1523     | A        | G           | 0.27%     |

        </td>

        <td>

### Control File 2 Summary
| Position | Original | Alternative | Frequency |
|----------|----------|-------------|-----------|
| 125      | T        | C           | 0.54%     |
| 158      | A        | G           | 0.29%     |
| 165      | T        | C           | 0.29%     |
| 222      | T        | C           | 0.26%     |
| 235      | T        | C           | 0.26%     |
| 291      | T        | C           | 0.30%     |
| 319      | T        | C           | 0.26%     |
| 370      | A        | G           | 0.27%     |
| 499      | A        | G           | 0.27%     |
| 566      | A        | G           | 0.28%     |
| 609      | A        | G           | 0.27%     |
| 722      | A        | G           | 0.28%     |
| 859      | A        | G           | 0.37%     |
| 913      | T        | C           | 0.27%     |
| 1031     | A        | G           | 0.36%     |
| 1209     | A        | G           | 0.33%     |
| 1213     | A        | G           | 0.37%     |
| 1358     | A        | G           | 0.31%     |
| 1366     | A        | G           | 0.30%     |
| 1460     | A        | G           | 0.41%     |
| 1482     | A        | G           | 0.26%     |
| 1517     | A        | G           | 0.30%     |
| 1520     | T        | C           | 0.27%     |
| 1600     | T        | C           | 0.37%     |

        </td>

        <td>

### Control File 3 Summary
| Position | Original | Alternative | Frequency |
|----------|----------|-------------|-----------|
| 105      | A        | G           | 0.94%     |
| 139      | T        | C           | 0.36%     |
| 158      | A        | G           | 0.36%     |
| 165      | T        | C           | 0.34%     |
| 235      | T        | C           | 0.37%     |
| 276      | A        | G           | 0.38%     |
| 297      | T        | C           | 0.40%     |
| 370      | A        | G           | 0.27%     |
| 414      | T        | C           | 0.32%     |
| 421      | A        | G           | 0.26%     |
| 566      | A        | G           | 0.40%     |
| 660      | A        | G           | 0.32%     |
| 722      | A        | G           | 0.31%     |
| 759      | T        | C           | 0.26%     |
| 859      | A        | G           | 0.26%     |
| 898      | A        | G           | 0.26%     |
| 915      | T        | C           | 0.33%     |
| 1031     | A        | G           | 0.31%     |
| 1086     | A        | G           | 0.28%     |
| 1209     | A        | G           | 0.43%     |
| 1213     | A        | G           | 0.33%     |
| 1358     | A        | G           | 0.37%     |
| 1366     | A        | G           | 0.33%     |
| 1421     | A        | G           | 0.35%     |
| 1460     | A        | G           | 0.34%     |
| 1482     | A        | G           | 0.26%     |

        </td>
    </tr>
</table>


