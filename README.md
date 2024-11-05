```markdown
# Lab Journal: Practical Training in Bioinformatics
**Author:** Elena Kozhevnikova  
**Date:** 02.11.2024  

## Project Overview
The aim of this project is to find putative antibiotic resistance genes in a laboratory *E. coli* strain compared to a natural non-resistant strain. We are provided with sequence reads of the strain under study and a reference wild-type genome.

### File Download
Download gz files using the `wget` command from the following URLs:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
wget https://figshare.com/articles/dataset/amp_res_2_fastq_zip/10006541/3?file=23769692
```

After downloading, we obtain the following files:
- `GCF_000005845.2_ASM584v2_genomic.fna.gz`
- `GCF_000005845.2_ASM584v2_genomic.gff.gz`
- `amp_res_1.fastq.gz` (forward read)
- `amp_res_2.fastq.gz` (reverse read)

### File Extraction
We gunzip the amp archives:
```bash
gunzip amp_res_1.fastq.gz
gunzip amp_res_2.fastq.gz
```

### Quality Control
Download FastQC software via `apt-get` and apply the following command for quality control:
```bash
fastqc -o ~/bio/project1/ ~/bio/project1/amp_res_1.fastq
fastqc -o ~/bio/project1/ ~/bio/project1/amp_res_2.fastq
```

These commands produce HTML files with QC reports:
- `amp_res_1_fastqc.html`
- `amp_res_2_fastqc.html`

#### QC Reports Indicate:
1. **Per Base Quality**: Warnings due to the lower quartile for bases being less than acceptable thresholds.
2. **Tile Quality Scores**: Issues might arise from transient or permanent problems on flowcell.
3. **Per Base Sequence Content**: Warnings occur if base proportions are out of expected ranges.

### Sequence Trimming
Sequence trimming is performed using Trimmomatic with the following parameters:
```bash
trimmomatic PE amp_res_1.fastq.gz amp_res_2.fastq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
LEADING:3 TRAILING:3 MINLEN:36
```

#### Summary of Reads:
- Input Read Pairs: 455,876
- Both Surviving: 445,940 (97.82%)
- Forward Only Surviving: 7,096 (1.56%)
- Reverse Only Surviving: 2,397 (0.53%)
- Dropped: 443 (0.10%)

### Indexing and Alignment
Build the index file:
```bash
bwa index GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Gunzip the paired files:
```bash
gunzip output_forward_paired.fq.gz
gunzip output_reverse_paired.fq.gz
```

Align to reference genome:
```bash
bwa mem GCF_000005845.2_ASM584v2_genomic.fna.gz output_forward_paired.fq output_reverse_paired.fq > alignment.sam
```

### Convert to BAM and QC
Compress the SAM file:
```bash
samtools view -S -b alignment.sam > alignment.bam
```

Get stats:
```bash
samtools flagstat alignment.bam
```
- Total: 892,165
- Mapped: 892,017 (99.98%)

### Sort and Index BAM File
Sort the BAM file:
```bash
samtools sort alignment.bam -o alignment_sorted.bam
```

Index the BAM file:
```bash
samtools index alignment_sorted.bam
```

### Variant Calling
Perform variant calling:
```bash
samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna alignment_sorted.bam > my.mpileup
```

Use VarScan to find putative mutations:
```bash
java -jar ~/bio/project1/varscan-2.4.6/VarScan.v2.4.0.jar mpileup2snp my.mpileup --min-var-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf
```

### SNP Annotation
Automatically annotate SNP using SnpEff:
1. Download GenBank format:
   [Download GenBank](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz)

2. Create an empty config file `snpEff.config` and add:
   ```
   k12.genome : ecoli_K12
   ```

3. Build the database:
```bash
mkdir -p data/k12
gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
mv GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk
snpEff build -dataDir ~/bio/project1/data -genbank -v k12
```

4. Annotate using SnpEff:
```bash
snpEff ann -dataDir ~/bio/project1/data k12 VarScan_results.vcf > VarScan_results_annotated.vcf
```

### Final Results
Annotated SNP results indicate:
**6 Variant Positions Reported**
| **Chromosome**  | **Position** | **Ref** | **Alt** | **Filter** | **ADP** | **WT** | **HET** | **HOM** | **NC** | **Annotation**                                                  |
|------------------|--------------|---------|---------|------------|---------|--------|---------|---------|--------|------------------------------------------------------------------|
| NC_000913.3      | 93043        | C       | G       | PASS       | 17      | 0      | 0       | 1       | 0      | G|missense_variant|MODERATE|ftsI|b0084|transcript         |
| NC_000913.3      | 482698       | T       | A       | PASS       | 17      | 0      | 0       | 1       | 0      | A|missense_variant|MODERATE|acrB|b0462|transcript        |
| NC_000913.3      | 852762       | A       | G       | PASS       | 15      | 0      | 0       | 1       | 0      | G|upstream_gene_variant|MODIFIER|glnH|b0811|transcript     |
| NC_000913.3      | 1905761      | G       | A       | PASS       | 15      | 0      | 0       | 1       | 0      | A|missense_variant|MODERATE|mntP|b1821|transcript      |
| NC_000913.3      | 3535147      | A       | C       | PASS       | 17      | 0      | 0       | 1       | 0      | C|missense_variant|MODERATE|envZ|b3404|transcript      |
| NC_000913.3      | 4390754      | G       | T       | PASS       | 16      | 0      | 0       | 1       | 0      | T|synonymous_variant|LOW|rsgA|b4161|transcript 
---

This concludes my lab journal entry for this practical training in bioinformatics.
```

