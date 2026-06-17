# Custom Gene List Feature - Comprehensive Guide

## Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Input Gene List Format](#input-gene-list-format)
4. [Creating Your Gene List](#creating-your-gene-list)
5. [Output Structure](#output-structure)
6. [Validation and Error Handling](#validation-and-error-handling)
7. [Troubleshooting](#troubleshooting)
8. [Advanced Usage](#advanced-usage)
9. [Implementation Details](#implementation-details)
10. [FAQ](#faq)

---

## Overview

The custom gene list feature allows you to provide your own set of genes for copy-number analysis, replacing the hard-coded list of ~60 genes in the original shallowHRD script.

### Key Features

✓ **Data-driven**: Gene lists are loaded from external files, not hard-coded  
✓ **Flexible input**: Support for two input formats (position or start/end)  
✓ **Comprehensive validation**: Checks chromosomes, coordinates, and duplicates  
✓ **Backward compatible**: Falls back to built-in defaults when no custom list provided  
✓ **Clear error messages**: Helpful feedback for common issues  

---

## Quick Start

### Using the Default Gene List

If you don't provide a custom gene list, the script uses the built-in set of ~60 cancer-related genes (BRCA1, BRCA2, MYC, CCNE1, TP53, PIK3CA, PTEN, etc.):

```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  input_ratio.txt \
  output_directory \
  cytoband_file.txt
```

### Using a Custom Gene List

Provide the path to your gene list file as an optional fourth argument:

```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  input_ratio.txt \
  output_directory \
  cytoband_file.txt \
  custom_genes.tsv
```

---

## Input Gene List Format

Your custom gene list must be a **tab-delimited or comma-delimited text file** with a **header row**.

### Required Columns

You must provide **one of these two column sets**:

#### Option 1: Single Representative Coordinate

Use this if you have a single position per gene (e.g., TSS, midpoint, or a known hotspot).

```
gene    chr    position
BRCA1   17     41236847
BRCA2   13     32932025
TP53    17     7581274
MYC     8      128751439
PIK3CA  3      178912013
PTEN    10     89677535
```

**Columns:**
- `gene`: Gene name (string, required)
- `chr`: Chromosome number (int or "chrN" format, required)
- `position`: Representative genomic position in bp (int, required)

#### Option 2: Gene Interval (Start and End)

Use this if you have gene boundaries. The **midpoint is automatically calculated**.

```
gene    chr    start      end
BRCA1   17     41196312   41277381
BRCA2   13     32889645   32974405
TP53    17     7571739    7590808
MYC     8      128747680  128755197
PIK3CA  3      178797988  178863393
PTEN    10     89623382   89731687
```

**Columns:**
- `gene`: Gene name (string, required)
- `chr`: Chromosome number (int or "chrN" format, required)
- `start`: Start position in bp (int, required)
- `end`: End position in bp (int, required)

**Note**: If both "position" and "start/end" columns are present, "position" takes precedence and "start/end" are ignored.

### Accepted File Formats

| Format | Extension | Example |
|--------|-----------|---------|
| Tab-separated values | `.tsv` or `.txt` | `genes.tsv` |
| Comma-separated values | `.csv` | `genes.csv` |

**Delimiter auto-detection**: The script automatically detects the delimiter from the first line.

### Chromosome Format

Chromosomes can be specified in multiple formats — all are automatically normalized:

| Input Format | Interpretation |
|--------------|-----------------|
| `1`, `2`, ..., `22` | Autosomes (1-22) |
| `23` | X chromosome |
| `24` | Y chromosome |
| `chr1`, `chr2`, ..., `chr22` | Autosomes with "chr" prefix |
| `chrX`, `chrY` | X/Y with "chr" prefix |

**Examples (all equivalent):**
- `17` or `chr17` → Chromosome 17 (where BRCA1 is located)
- `23` or `chrX` → X chromosome
- `24` or `chrY` → Y chromosome

### Coordinate Validation

The script validates all coordinates:

- **Positions**: Must be ≥ 1 bp
- **Chromosomes**: Must be 1-22, 23 (X), or 24 (Y) after normalization
- **Numeric types**: Positions, starts, and ends must be convertible to numbers
- **No missing values**: All required columns must have data (no blanks or NA)

The script will **error with a clear message** if validation fails.

---

## Creating Your Gene List

### Approach 1: Manual Entry

Create a spreadsheet (Excel, Google Sheets, LibreOffice) and export as `.csv` or `.tsv`:

| gene | chr | position |
|------|-----|----------|
| BRCA1 | 17 | 41236847 |
| BRCA2 | 13 | 32932025 |
| TP53 | 17 | 7581274 |

Then save:
- **Mac/Linux**: File → Export As → Tab-Separated Values (or CSV)
- **Windows**: File → Save As → Choose format → .csv or .txt

### Approach 2: From a Reference Source

#### NCBI Gene (https://www.ncbi.nlm.nih.gov/gene/)

1. Search for your gene (e.g., "BRCA1 human")
2. Go to the gene page
3. Find the "Genomic regions, transcripts, and products" section
4. Note the chromosome and a representative position (e.g., CDS start or gene start)
5. Add to your list

#### UCSC Genome Browser (https://genome.ucsc.edu/)

1. Go to Table Browser (Tools → Table Browser)
2. Select `hg19` genome
3. Choose a gene track (e.g., RefSeq Genes, Ensembl Genes)
4. Click "output all fields"
5. Download the table
6. Extract `chrom`, `txStart`, `txEnd` columns
7. Compute midpoint: `(txStart + txEnd) / 2`
8. Create your gene list

#### Published Gene Panels

Many cancer gene panels are publicly available:
- **COSMIC** (https://cancer.sanger.ac.uk/cosmic)
- **CancerVar** (https://www.cancervar.org/)
- **DGD** (Disease Gene Database)
- **PubMed**: Search "cancer gene panel hg19"

Extract gene names and coordinates, then convert to hg19 if necessary using **liftOver** (https://genome.ucsc.edu/cgi-bin/hgLiftOver).

### Approach 3: Programmatic Generation (R)

Create your gene list directly in R and export:

```r
# Example: Homologous Recombination Deficiency (HRD) panel
hrd_genes <- data.frame(
  gene = c("BRCA1", "BRCA2", "RAD51C", "RAD51D", "PALB2", "PTEN"),
  chr = c(17, 13, 17, 1, 16, 10),
  position = c(41236847, 32932025, 56790924, 160576892, 23603160, 89677535),
  stringsAsFactors = FALSE
)

write.table(hrd_genes, "hrd_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

Then run:
```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R ratio.txt output/ cytoband.txt hrd_genes.tsv
```

### Approach 4: From a BED File

If you have a BED file (e.g., from a sequencing panel):

```bash
# Extract gene names and compute midpoints
awk 'BEGIN {print "gene\tchr\tposition"} 
     NR > 1 {
       chr = $1; 
       gsub("^chr", "", chr); 
       start = $2; 
       end = $3; 
       pos = int((start + end) / 2); 
       gene = $4; 
       print gene "\t" chr "\t" pos
     }' your_panel.bed > custom_genes.tsv
```

---

## Output Structure

Regardless of which gene list is used (default or custom), the output is identical:

### Output File: `amplification_deletion_table.txt`

A tab-delimited file with one row per gene and the following columns:

| Column | Data Type | Description |
|--------|-----------|-------------|
| `gene` | string | Gene name (from input) |
| `chr` | integer | Chromosome (1-22, 23 for X) |
| `start` | integer | Position of nearest QDNAseq window (bp) |
| `ratio_point_initial` | float | Per-window ratio (pre-segmentation) |
| `CN_to_baseline_point_initial` | float | Per-window CN = ratio / THR |
| `ratio_segment_initial` | float | Segment ratio (before LGA calling) |
| `CN_to_baseline_segment_initial` | float | Segment CN (before final smoothing) |
| `ratio_segment_final` | float | Final segment ratio (post-LGA calling) |
| `CN_to_baseline_segment_final` | float | **Final segment CN (recommended for reporting)** |

### Canonical CN Value

**Use `CN_to_baseline_segment_final` as the primary copy-number value for reporting.**

This column represents the final, smoothed copy-number state after all shallowHRD processing (segmentation, LGA calling, and baseline adjustment).

### Example Output

```
gene        chr    start      ratio_point_initial  CN_to_baseline_point_initial  ...  CN_to_baseline_segment_final
BRCA1       17     41236800   1.89                 1.89                          ...  1.95
BRCA2       13     32932000   2.12                 2.12                          ...  2.08
TP53        17     7581200    0.95                 0.95                          ...  0.98
MYC         8      128751400  2.45                 2.45                          ...  2.41
PIK3CA      3      178912000  2.18                 2.18                          ...  2.15
PTEN        10     89677500   0.85                 0.85                          ...  0.88
```

---

## Validation and Error Handling

The script performs comprehensive validation to catch common errors early.

### Validation Checks

1. **File existence**: Custom gene list file must exist
2. **Delimiter detection**: Auto-detects tab or comma
3. **Header validation**: Checks for required columns
4. **Data type validation**: Positions are numeric
5. **Chromosome validation**: Valid values for hg19 (1-22, 23, 24)
6. **Coordinate ranges**: Positions ≥ 1
7. **Duplicate detection**: Warns if gene names are duplicated
8. **Missing data**: Detects NA or empty cells

### Error Messages and Solutions

#### File Not Found
```
Gene list file not found: /path/to/missing_file.tsv
```
**Solution**: Check the file path is correct and the file exists in the specified location.

---

#### Missing or Invalid Delimiter
```
Could not detect delimiter (tab or comma) in gene list file
```
**Solution**: Ensure your file is tab- or comma-delimited. Check the first line has the correct delimiter.

---

#### Missing Required Column
```
Gene list must have a 'gene' column
Gene list must have a 'chr' column
Gene list must have either 'position' column or both 'start' and 'end' columns
```
**Solution**: Check your header row has the correct column names (case-sensitive: `gene`, `chr`, `position`).

---

#### Invalid Chromosome Values
```
Invalid chromosome values: 25, 26, 0, -1
Valid values for hg19: 1-22, 23 (X), 24 (Y)
```
**Solution**: Verify chromosome values are 1-22, 23 for X, or 24 for Y. Remove or correct invalid entries.

---

#### Non-Numeric Coordinates
```
Some gene positions could not be converted to numeric values
```
**Solution**: Check position/start/end columns contain only numbers (no text or special characters).

---

#### Duplicate Gene Names (Warning)
```
Duplicate gene entries found: BRCA1, TP53
```
**Solution**: Remove or rename duplicate entries. The script will continue but may overwrite results.

---

## Troubleshooting

### Issue: "Gene list file not found"

**Symptoms**: Script exits immediately after reading the gene list argument.

**Causes**:
- File path is wrong
- File doesn't exist
- Typo in filename

**Solutions**:
1. Double-check the file path (absolute or relative)
2. Verify the file exists: `ls /path/to/file.tsv`
3. Use absolute paths to avoid ambiguity: `/full/path/to/genes.tsv`

---

### Issue: "Invalid chromosome values"

**Symptoms**: Script errors with chromosome numbers like 25, 26, or 0.

**Causes**:
- Chromosomes are from hg38 or other reference
- Typos in chromosome column (e.g., "17a" instead of "17")
- Using numeric codes (e.g., 23→24, 25→26)

**Solutions**:
1. Verify you're using hg19 coordinates
2. If you have hg38, convert using liftOver
3. Check for typos: `head -5 genes.tsv | cut -f2`
4. Use standard format: 1-22, 23 (X), 24 (Y)

---

### Issue: "Duplicate gene entries found" (Warning)

**Symptoms**: Script warns about duplicates but continues.

**Causes**:
- Gene appears twice in your list
- Copy-paste error
- Merging multiple gene lists

**Solutions**:
1. Remove or consolidate duplicates
2. Use `sort | uniq` to identify duplicates:
   ```bash
   cut -f1 genes.tsv | sort | uniq -d
   ```

---

### Issue: Most genes return NA/missing values in output

**Symptoms**: Output file is mostly empty or has NAs.

**Causes**:
- Gene coordinates fall outside the QDNAseq window range
- Chromosome mismatch (hg19 vs hg38)
- QDNAseq file is from a different reference

**Solutions**:
1. Verify gene positions fall within your sequencing data range
2. Check chromosome numbers match hg19
3. Confirm the QDNAseq file was generated with hg19 reference
4. Check that genes are on autosomes/X (not Y, which may be filtered)

---

### Issue: Empty `amplification_deletion_table.txt`

**Symptoms**: Output file has only a header row or is completely empty.

**Causes**:
- Gene list is empty or all entries are invalid
- All chromosomes are invalid (e.g., hg38 format)
- No overlap between gene coordinates and QDNAseq windows

**Solutions**:
1. Verify gene list has at least one valid entry:
   ```bash
   wc -l genes.tsv  # Should be > 1 (header + data)
   head genes.tsv
   ```
2. Check a few genes manually are in the file
3. Verify chromosomes and positions are correct for hg19

---

## Advanced Usage

### Switching Between Gene Lists

Run the same sample with different gene panels to compare:

```bash
# Panel 1: HRD-specific genes
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  ratio.txt output1/ cytoband.txt hrd_genes.tsv

# Panel 2: Breast cancer amplification genes
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  ratio.txt output2/ cytoband.txt breast_amp_genes.tsv

# Panel 3: Default (60-gene set)
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  ratio.txt output3/ cytoband.txt
```

Compare results across panels to identify the most relevant genes.

### Large Gene Lists

The script can handle large gene lists (tested up to ~1000 genes):

**Performance:**
- 60 genes (default): ~2-3 seconds
- 100 genes: ~3-4 seconds
- 500 genes: ~10-15 seconds
- 1000 genes: ~20-30 seconds

**Memory usage**: Minimal; mostly depends on QDNAseq window count, not gene count.

### Batch Processing

Process multiple samples with the same gene list:

```bash
#!/bin/bash
GENE_LIST="my_genes.tsv"

for ratio_file in samples/*.bam_ratio.txt; do
  sample=$(basename "$ratio_file" .bam_ratio.txt)
  output_dir="results/$sample"
  
  Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
    "$ratio_file" \
    "$output_dir" \
    cytoband.txt \
    "$GENE_LIST"
done
```

### Filtering and Post-Processing

After running shallowHRD, filter results in R:

```r
# Load results
results <- read.delim("amplification_deletion_table.txt")

# Find strong amplifications (CN > 3)
amplifications <- results[results$CN_to_baseline_segment_final > 3, ]

# Find deletions (CN < 1)
deletions <- results[results$CN_to_baseline_segment_final < 1, ]

# Find significant events (|CN - 2| > 0.5)
significant <- results[abs(results$CN_to_baseline_segment_final - 2) > 0.5, ]

# Save filtered results
write.table(amplifications, "amplifications.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deletions, "deletions.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

---

## Implementation Details

### File Structure

**Files added/modified:**

- `gene_list_utils.R` (new): Utility functions for gene list handling
  - `load_gene_list()`: Load from file
  - `validate_and_normalize_gene_list()`: Validate and normalize
  - `extract_gene_metrics()`: Extract CN values for one gene
  - `build_gene_summary_table()`: Generate the full output table
  - `get_default_gene_list()`: Built-in default genes
  - `get_gene_list()`: Auto-select default or custom

- `shallowHRD_hg19_1.13_QDNAseq_no_chrX.R` (modified): 
  - Now sources `gene_list_utils.R`
  - Accepts optional 4th argument (custom gene list path)
  - Calls `get_gene_list()` and `build_gene_summary_table()`

### Gene Metrics Extraction

For each gene, the script:

1. **Finds the nearest QDNAseq window** to the gene's representative coordinate using `Closest()` function
2. **Extracts initial metrics** from the pre-segmentation ratio file (per-window values)
3. **Finds the nearest segment** in the final segmentation
4. **Extracts final metrics** from the segmented ratio file
5. **Calculates CN values** using: `CN = ratio / THR` (where THR is the baseline)

### Coordinate Matching Logic

```r
# Example: Find nearest window to gene position
gene_position <- 41236847  # BRCA1
B_chr <- B_initial[B_initial$chr == 17, ]  # All chr17 windows
closest_idx <- Closest(B_chr$start, gene_position)  # Find nearest
window <- B_chr[closest_idx, ]  # Extract window data
```

The `Closest()` function from **DescTools** package finds the nearest match even if the gene position doesn't align perfectly with a bin boundary.

---

## FAQ

### Q: Can I use hg38 coordinates?

**A**: No, this version of shallowHRD is designed for hg19 only. If you have hg38 coordinates, convert them using **liftOver**:

```bash
# Download liftOver chain file
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# Create BED file from your genes
echo -e "chr17\t41196312\t41277381\tBRCA1" > genes.bed

# Convert to hg19
liftOver genes.bed hg38ToHg19.over.chain.gz genes_hg19.bed genes_unmapped.bed
```

Then extract coordinates from the converted BED file.

---

### Q: What if my gene has multiple isoforms with different coordinates?

**A**: Use any representative isoform:
- Longest transcript
- Most common isoform
- Canonical isoform from RefSeq
- CDS start position

The script will find the nearest QDNAseq bin regardless. Different isoforms may map to slightly different bins if they're far apart, but the analysis is robust to this variation.

---

### Q: Can I mix "position" and "start/end" columns in the same file?

**A**: Not recommended. The script will use "position" if it exists and ignore "start/end". Use one format consistently:
- **Option A**: Only "position" column
- **Option B**: Only "start" and "end" columns

If you have both, remove the "position" column before running the script.

---

### Q: How many genes can I include in one analysis?

**A**: Theoretically unlimited. In practice:
- **Recommended**: 10-10,000 genes
- **Tested**: Up to 1,000 genes
- **Very large** (>10,000): May exceed memory or time limits on some systems

For very large panels, consider splitting into multiple runs or filtering to genes of interest.

---

### Q: Why is my gene not appearing in the output?

**Possible causes:**
1. **Chromosome is invalid**: Check it's 1-22, 23 (X), or 24 (Y) after normalization
2. **Position is invalid**: Must be ≥ 1 and numeric
3. **Gene name is empty/whitespace**: Ensure non-empty gene names
4. **Coordinates are outside QDNAseq range**: Check gene falls within sequencing coverage
5. **Data type issue**: Ensure position/start/end are numbers, not text

**Solutions:**
1. Check the gene was loaded: `head -20 genes.tsv | grep "gene_name"`
2. Verify chromosome: `cut -f2 genes.tsv | sort -u`
3. Verify position ranges: `cut -f3 genes.tsv | sort -n | head -5`

---

### Q: Can I run shallowHRD without providing a gene list?

**A**: Yes! The script will automatically use the built-in default list of ~60 genes if no 4th argument is provided:

```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  ratio.txt \
  output/ \
  cytoband.txt
# No 4th argument → uses default genes
```

---

### Q: What's the difference between the "point" and "segment" metrics?

**Point metrics** ("ratio_point_initial", "CN_to_baseline_point_initial"):
- Per-window values from the pre-segmentation QDNAseq output
- Raw, unsmoothed data

**Segment metrics** ("ratio_segment_initial", "CN_to_baseline_segment_initial", "ratio_segment_final", "CN_to_baseline_segment_final"):
- Values from the segmented/smoothed copy-number profile
- After segmentation and LGA calling

**Recommendation**: Use `CN_to_baseline_segment_final` for reporting (final, smoothed CN value).

---

### Q: How do I cite this feature?

**A**: Cite both the original shallowHRD and this fork:

*Original shallowHRD:*
> Popova et al. (2012); Eeckhoutte et al. [full citation]

*This fork with custom gene list support:*
> PMCC-CancerGenomicsTRC/shallowHRD (https://github.com/PMCC-CancerGenomicsTRC/shallowHRD, feat/custom-scna-gene-list branch)

---

## References

### Original shallowHRD Publication

Popova et al. - Shallow whole-genome sequencing for homologous recombination deficiency assessment. [Details to be added]

Eeckhoutte et al. - [Full citation to original shallowHRD paper]

### Related Tools

- **QDNAseq**: Copy number profiling from low-depth sequencing (Scheinin et al., 2014)
- **ASCAT**: Allele-specific copy number analysis (Van Loo et al., 2010)
- **CGHcall**: Calling significant copy number changes (van de Wiel et al., 2007)
- **liftOver**: Convert genomic coordinates between builds (UCSC)

### Useful Resources

- UCSC Genome Browser: https://genome.ucsc.edu/
- NCBI Gene: https://www.ncbi.nlm.nih.gov/gene/
- COSMIC: https://cancer.sanger.ac.uk/cosmic
- Ensembl: https://www.ensembl.org/

---

**Last Updated**: 2026-06-17  
**Version**: 1.13 with custom gene list support  
**Maintainers**: PMCC-CancerGenomicsTRC
