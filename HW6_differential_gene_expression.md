
# Yeast RNA Expression Changes During Fermentation

## Overview
We will explore how RNA expression levels change as yeast undergo fermentation to make bread rise. There are two replicates of RNA-seq data from yeast before and during fermentation, and our goal is to find out if the yeast express different genes during fermentation than they do under normal growth.

## Data Download
### RNA-seq Data
- **SRR941816**: Fermentation 0 minutes replicate 1  
  [Download](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941816/SRR941816.fastq.gz) (413 Mb)  
- **SRR941817**: Fermentation 0 minutes replicate 2  
  [Download](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941817/SRR941817.fastq.gz) (455 Mb)  
- **SRR941818**: Fermentation 30 minutes replicate 1  
  [Download](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941818/SRR941818.fastq.gz) (79.3 Mb)  
- **SRR941819**: Fermentation 30 minutes replicate 2  
  [Download](ftp.sra.ebi.ac.uk/vol1/fastq/SRR941/SRR941819/SRR941819.fastq.gz) (282 Mb)  

### Reference Genome
As a reference genome, we will use *Saccharomyces cerevisiae*, in the genome database at NCBI. Make sure you have strain S288c and assembly R64. Download the reference genome in FASTA format and annotation in GFF format.

- **Reference Genome File**:  
  [GCF_000146045.2_R64_genomic.fna.gz](ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz)

- **Annotation File**:  
  [GCF_000146045.2_R64_genomic.gff.gz](ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz)

## Index the Genome
```bash
conda install -c bioconda hisat2
gunzip GCF_000146045.2_R64_genomic.fna.gz
hisat2-build GCF_000146045.2_R64_genomic.fna yeast_index
```

## Align Each Sample
```bash
hisat2 -x yeast_index -U SRR941816.fastq.gz -S SRR941816.sam
hisat2 -x yeast_index -U SRR941817.fastq.gz -S SRR941817.sam
hisat2 -x yeast_index -U SRR941818.fastq.gz -S SRR941818.sam
hisat2 -x yeast_index -U SRR941819.fastq.gz -S SRR941819.sam
```

## Convert SAM to BAM
```bash
samtools view -Sb SRR941816.sam > SRR941816.bam
samtools view -Sb SRR941817.sam > SRR941817.bam
samtools view -Sb SRR941818.sam > SRR941818.bam
samtools view -Sb SRR941819.sam > SRR941819.bam
```

## Count Reads per Gene
```bash
gunzip GCF_000146045.2_R64_genomic.gff.gz
conda install -c bioconda subread
featureCounts -a GCF_000146045.2_R64_genomic.gff -o counts.txt -g "locus_tag" -t "CDS" *.bam
```

### Analyzing Counts in R
We then use `counts.txt` to analyze in R using edgeR package. It allows for the Differential expression analysis of RNA-seq expression profiles with biological replication. Implements a range of statistical methodology based on the negative binomial distributions, including empirical Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests. As well as RNA-seq, it be applied to differential signal analysis of other types of genomic data that produce read counts, including ChIP-seq, ATAC-seq, Bisulfite-seq, SAGE and CAGE.

```r
# Load necessary library
library(edgeR)

# Set working directory
setwd("C:/Users/Елена/Downloads/Command_line")

# Load count data
count_data <- read.table("counts.txt", header=TRUE, row.names=1, comment.char="#")
counts <- count_data[, 6:ncol(count_data)]

# Define experimental groups
group <- factor(c("0min", "0min", "30min", "30min"))

# Create DGEList object
dge <- DGEList(counts=counts, group=group)

# Filter low-expressed genes
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize the data
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~group)
dge <- estimateDisp(dge, design)

# Fit the model and test for differential expression
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit)

# View significant genes
print(topTags(qlf))

# Collect results into a dataframe
all_results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
all_results$comparison <- "30min_vs_0min"

# Output results
write.csv(all_results, "edgeR_all_results_with_comparison.csv")

sig_genes <- all_results[all_results$FDR < 0.05 & abs(all_results$logFC) > 1, ]
write.csv(sig_genes, "edgeR_significant_genes.csv", row.names = TRUE)

clean_results <- all_results[, c("logFC", "PValue", "FDR")]
write.csv(clean_results, "edgeR_clean_results.csv", row.names = TRUE)

# Visualization
plotMD(qlf)
summary(decideTests(qlf))
```

## Explanation of Main Code
### **1. Creating DGEList**
- **Purpose**: Creates a DGEList object (edgeR's core data structure).
- **Inputs**: Count data and group information.

### **2. Filtering Low-Expressed Genes**
```r
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes=FALSE]
```
- Converts counts to CPM and identifies genes with CPM > 1 in each sample.  
- Retains genes with CPM > 1 in at least 2 samples.

### **3. Normalization** 
```r
dge <- calcNormFactors(dge)
```
- Adjusts for library composition bias using TMM normalization.

### **4. Design Matrix** 
```r
design <- model.matrix(~group)
```
- Defines the experimental design for linear modeling.

### **5. Dispersion Estimation** 
```r
dge <- estimateDisp(dge, design)
```
- Estimates dispersions for the count data.

### **6. Model Fitting & Testing** 
```r
fit <- glmQLFit(dge, design)  # Fits quasi-likelihood model
qlf <- glmQLFTest(fit)        # Tests for differential expression
```

### **Key Notes**
1. **CPM Filtering**: Removes noise from lowly expressed genes.
2. **TMM Normalization**: Accounts for sample-specific biases.
3. **Quasi-Likelihood**: Robust to overdispersion in RNA-seq counts.
4. **Output**: Use `topTags(qlf)` to extract significant DEGs.

## Visualization of Results
```r
# Heatmap Generation
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

setwd("C:/Users/Елена/Downloads/Command_line")

# Log2 CPM normalization
norm_counts <- cpm(dge, log = TRUE)
sig_genes <- rownames(all_results[all_results$FDR < 0.05 & abs(all_results$logFC) > 1.2, ])
sig_norm_counts <- norm_counts[sig_genes, ]
scaled_data <- t(scale(t(sig_norm_counts)))

# Prepare heatmap data
heatmap_data <- as.data.frame(sig_norm_counts) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression")

# Sample annotation
sample_annot <- data.frame(
  Sample = colnames(sig_norm_counts),
  Group = c("0min", "0min", "30min", "30min")
)

heatmap_data <- heatmap_data %>% left_join(sample_annot, by = "Sample")

# Generating PDF
pdf("DEG_heatmap.pdf", width = 5, height = 7.5)
heatmap(
  scaled_data,
  col = colorRampPalette(c("royalblue3", "white", "firebrick3"))(100),
  scale = "none",
  labRow = NA,
  margins = c(8, 2),
  main = paste("Differentially Expressed Genes (n =", nrow(scaled_data), ")")
)

# Legend
legend("topright",
       legend = c("0 min", "30 min"),
       fill = c("#1F78B4", "#E31A1C"),
       border = NA,
       bty = "n")

dev.off()

# Volcano Plot Generation
ggplot(plot_data, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c("Down" = "royalblue3", "NS" = "gray80", "Up" = "firebrick3"),
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  labs(
    x = expression(log[2]~("Fold Change")),
    y = expression(-log[10]~("FDR")),
    title = "Differential Expression: 30min vs 0min"
  ) +
  theme_minimal() +
  ggsave("volcano_clean.pdf", width = 6, height = 5)
```

![image](https://github.com/user-attachments/assets/3521962f-7158-48fd-81f7-76a43a69d9a6)

![image](https://github.com/user-attachments/assets/96e7c3d0-818d-40d9-a893-9cb9a43e38c9)

## Gene Ontology Analysis
Then we used gene names to analyze them for gene ontology using the online service ShinyGO ([link](https://bioinformatics.sdstate.edu/go/)). The gene names were extracted as follows:

```python
import pandas as pd
from goatools import obo_parser

datafile = "/home/lena/bio/project6/edgeR_clean_results.csv"
df = pd.read_csv(datafile)
df.columns = ["Gene"] + list(df.columns[1:])
significant_genes = df[df['FDR'] < 0.05].sort_values('FDR')
gene_list = significant_genes['Gene'].tolist()

print("Significant genes sorted by FDR:")
print(significant_genes)
print("\nGene list:")
print(gene_list)
```

![image](https://github.com/user-attachments/assets/57c81b3e-43bc-4458-8b28-4446aff552a1)

### Results of GO Analysis
The result of GO analysis revealed the following pathways:
| Enrichment FDR | nGenes | Pathway                                                 | Fold Enrichment |
|----------------|--------|---------------------------------------------------------|-----------------|
| 4.7E-05        | 28     | Citrate cycle (TCA cycle)                               | 1.8             |
| 2.4E-05        | 32     | Proteasome                                             | 1.8             |
| 9.5E-07        | 60     | Ribosome biogenesis in eukaryotes                      | 1.6             |
| 3.5E-02        | 16     | Fatty acid degradation                                   | 1.6             |
| 2.3E-12        | 121    | Ribosome                                               | 1.6             |
| 1.7E-02        | 24     | RNA polymerase                                         | 1.5             |
| 2.5E-02        | 22     | Glyoxylate and dicarboxylate metabolism                | 1.5             |
| 1.6E-03        | 42     | Purine metabolism                                       | 1.5             |
| 3.8E-03        | 38     | Pyruvate metabolism                                     | 1.5             |
| 2.2E-02        | 28     | Longevity regulating pathway-multiple species          | 1.5             |
| 2.3E-02        | 29     | Starch and sucrose metabolism                           | 1.4             |
| 2.3E-02        | 29     | Glycerophospholipid metabolism                          | 1.4             |
| 2.3E-02        | 29     | Peroxisome                                             | 1.4             |
| 2.5E-02        | 30     | Nucleotide metabolism                                   | 1.4             |
| 1.8E-02        | 38     | Glycolysis/Gluconeogenesis                              | 1.4             |
| 4.8E-04        | 77     | Carbon metabolism                                       | 1.4             |
| 5.6E-10        | 234    | Biosynthesis of secondary metabolites                   | 1.3             |
| 1.6E-03        | 83     | Biosynthesis of amino acids                             | 1.3             |
| 1.2E-02        | 60     | Autophagy-yeast                                        | 1.3             |
| 4.9E-15        | 505    | Metabolic pathways                                      | 1.3             |

## Conclusion
The citrate cycle, also known as the tricarboxylic acid (TCA) cycle or Krebs cycle, is a crucial metabolic pathway that plays a central role in cellular respiration and energy production in aerobic organisms, including yeast. It is responsible for the oxidative metabolism of carbohydrates, fats, and proteins into carbon dioxide and water, while generating energy-rich molecules such as ATP, NADH, and FADH₂.

### Key Roles of the Citrate Cycle
1. **Energy Generation**: The TCA cycle processes acetyl-CoA, releasing energy stored in acetyl-CoA as it is oxidized to CO₂.
2. **Biosynthetic Precursors**: TCA cycle intermediates serve as precursors for various biosynthetic pathways.
3. **Feedback Regulation**: The cycle is regulated by substrate availability and the energy needs of the cell.

### Impact of Fermentation on the TCA Cycle in Yeast
- **Aerobic vs. Anaerobic Conditions**: Yeast relies on glycolysis and alcoholic fermentation under anaerobic conditions, with reduced reliance on the TCA cycle.
- **Nutrient Availability**: Substrate availability affects TCA cycle activity based on the conversion of sugars into acetyl-CoA.
- **Stress Response and Adaptation**: Fermentation induces stress responses that can alter the expression of TCA cycle genes.

**Conclusion**: The increased activity of the citrate cycle in your analysis suggests significant adjustments in the metabolic state of yeast during or after fermentation, reflecting changes in energy demand, substrate availability, and metabolic regulation. The findings are critical for understanding metabolic regulation and cell growth in yeast, particularly in fermentation processes.
```
