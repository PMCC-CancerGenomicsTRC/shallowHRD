# shallowHRD

> **Note:** This is a forked repository with added support for custom gene lists. See "Custom Gene List Support" section below. For the original shallowHRD implementation, see the [primary repository](https://github.com/aeeckhou/shallowHRD).

## Fork Attribution

This fork extends the original **shallowHRD** software with data-driven gene list support. Original credit:

- **Original Authors & Publication:** Popova et al. (2012); Eeckhoutte et al.
- **Original Repository:** [aeeckhou/shallowHRD](https://github.com/aeeckhou/shallowHRD)
- **This Fork:** PMCC-CancerGenomicsTRC/shallowHRD (feat/custom-scna-gene-list branch)

### Modifications in This Fork

- **Custom gene list support**: Replaces hard-coded gene annotations with external, user-provided gene lists
- **Data-driven workflow**: Gene-specific metrics are now generated programmatically
- **Enhanced validation**: Comprehensive checks for input gene list format and chromosome values
- **Backward compatibility**: Falls back to built-in 60-gene default when no custom list is provided

---

## Custom Gene List Support

### Quick Start

**Using the default gene list (original behavior):**
```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  input_ratio.txt \
  output_directory \
  cytoband_file.txt
```

**Using a custom gene list:**
```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R \
  input_ratio.txt \
  output_directory \
  cytoband_file.txt \
  custom_genes.tsv
```

### Input Gene List Format

Provide a tab- or comma-delimited file with header row. Choose one of:

**Option A: Single representative coordinate per gene**
```
gene    chr    position
BRCA1   17     41236847
BRCA2   13     32932025
TP53    17     7581274
```

**Option B: Gene interval (start/end) — midpoint is derived automatically**
```
gene    chr    start      end
BRCA1   17     41196312   41277381
BRCA2   13     32889645   32974405
TP53    17     7571739    7590808
```

### Input Validation

- ✓ Chromosomes: numeric (1-22, 23 for X, 24 for Y) or "chr1"-"chr22", "chrX", "chrY"
- ✓ Positions: numeric, ≥ 1 bp
- ✓ Delimiter: auto-detected (tab or comma)
- ✓ Duplicates: detected and warned
- ✓ Missing values: detected and reported

### Output

Same output structure as original:
- **File**: `amplification_deletion_table.txt`
- **Columns**: gene, chr, start, ratio_point_initial, CN_to_baseline_point_initial, ratio_segment_initial, CN_to_baseline_segment_initial, ratio_segment_final, **CN_to_baseline_segment_final** (final CN value)
- **One row per gene**

### Implementation

- **`gene_list_utils.R`**: Utility functions for loading, validating, and processing gene lists
- **`shallowHRD_hg19_1.13_QDNAseq_no_chrX.R`**: Updated main script (now sources gene_list_utils.R)

### Example: Creating a Custom Gene List

```r
# Minimal HRD panel
hrd_genes <- data.frame(
  gene = c("BRCA1", "BRCA2", "RAD51C", "PTEN"),
  chr = c(17, 13, 17, 10),
  position = c(41236847, 32932025, 56790924, 89677535)
)

write.table(hrd_genes, "hrd_genes.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

Then run:
```bash
Rscript shallowHRD_hg19_1.13_QDNAseq_no_chrX.R ratio.txt output/ cytoband.txt hrd_genes.tsv
```

### Full Documentation

For comprehensive details on the custom gene list feature, validation rules, troubleshooting, and examples, see [CUSTOM_GENE_LIST.md](CUSTOM_GENE_LIST.md).

---

# shallowHRD (Original Documentation)

This method uses shallow Whole Genome Sequencing (sWGS > 0.3x) and the segmentation of a tumor genomic profile to infer the Homologous Recombination status of a breast and ovarian tumor based on th[...]

## Introduction

*shallowHRD* is a R script that can be launched from the command line. It relies on a ratio file characterizing the normalized read counts of a shallow Whole Genome Sequencing (>0.3x) in sliding wi[...]

Softwares such as QDNAseq count reads in sliding windows, normalize read count and then segment the genomic profile. *shallowHRD*, based on a inferred CNA cut-off representing a one copy difference[...]

**IMPORTANT : This GitHub contains the *first* version of *shallowHRD* (v1.13). Since its publication, the software has been under continuous developpement and the *shallowHRDv2* has been publishe[...]
The version 2.0 is now available under licence. Please contact Marc-Henri Stern : marc-henri.stern@curie.fr** 

## Requirements

* R installed (tested with v.4.1.0)
* The following packages installed : 
  * ggpubr (tested with v.0.4.0)
  * gridExtra (tested with v.2.3)
  * DescTools (tested with v.0.99.42)
  * GenomicRanges (tested with v.1.44.0)
  * ks (tested with v.1.13.2)
  * ggrepel (tested with v.0.9.1)

Tested on Linux, Mac and Windows.

## Prerequisities

First, FASTQ files should be aligned to a reference genome (hg19 or hg38) (using [BWA-MEM](https://github.com/lh3/bwa) for instance) and supplementary & duplicate reads removed from the BAM files,[...]

**IMPORTANT: Please only use chromosomes 1 to 22 (plus the Chromosome X if you want to) for the alignment step. Additionnal chromosomes (contigs) might introduce errors.**

Then, the BAM file should then be processed by a software such as ControlFREEC. The recommended options for controlFREEC are indicated in a config file example in the repository (*controlfreec_con[...]

Finally, the file *cytoBand_adapted_hg19.csv* or *cytoBand_adapted_hg38.csv* (available in the repository) has to be downloaded. 

The R packages needed can be installed with the script *install_packages.R* (in repository) and the command line :

```
/path/to/Rscript /path/to/install.packages.R
```

## Run *shallowHRD*

To run *shallowHRD* only one ratio file is needed (formated in ControlFREEC's output).

The name of the file should be in this format : *SAMPLE_NAME.bam_ratio.txt*. <br/>

*shallowHRD* will rely on the *first four columns* of the input file (tabulated and with column Chromosome *in number*) : <br/>
Chromosome &nbsp; Start &nbsp; Ratio &nbsp; RatioMedian <br/>
1 &nbsp;&nbsp; 1 &nbsp;&nbsp;&nbsp; -1 &nbsp;&nbsp;&nbsp; -1 <br/>
1 &nbsp;&nbsp; 20001 &nbsp;&nbsp; -1 &nbsp;&nbsp; -1 <br/>
. &nbsp;&nbsp; . &nbsp;&nbsp; . &nbsp;&nbsp; . <br/>
. &nbsp;&nbsp; . &nbsp;&nbsp; . &nbsp;&nbsp; . <br/>

The command line to launch *shallowHRD* is (absolute or relative paths) :

```
/path/to/Rscript /path/to/shallowHRD_hg19.R /path/to/SAMPLE_NAME.bam_ratio.txt /path/to/output_directory /path/to/cytoBand_adapted_hg19.csv
```
For Windows, it will be with /path/to/Rscript.exe.

Two examples in hg19 and one example in hg38 are downloadable in the repository to try *shallowHRD*.

## Outputs

All the figures and files created by the script will be available in the output directory.

The summary plot figure recapitulating all the information will look like this :

![alt text](https://github.com/aeeckhou/shallowHRD/blob/master/example_1_QDNAseq_final_summary_plot_hg19.jpeg)

A : Genomic profile with LGAs in green (the entire processed segmentation is represented in red if there are no LGA) <br/>
B : Density representing pairwise comparison between large segments used to fix the difference for a copy level <br/>
C : Graphe representing the value of each final segment (small blue circles) - <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If the segmentation is good, the different copy number should appear clearly with disctinct steps <br/>
D : Table recapitulating different data, including the case quality and the final diagnostic for the HR status

## Nota Bene

1. The scripts for *QDNAseq* and *controlfreec* have been updated to the 1.13 version. They harbor a more robust CNA cut-off detection and overall optimization of the profiles   

2. The 1.13 version of *shallowHRD* is more robust and reliable but takes a longer time to run compared to older version (~1 hour by sample)

2. Different scripts for *QDNAseq* and *controlfreec* are available depending on whether the chromosome X is included in the ratio file 

3. *shallowHRD* can be adapted to other softwares with slight modification of outputs to match *shallowHRD* intput format <br/> 

4. The overall pipeline works also on WGS with a higher coverage   

## Contact

Don't hesitate to contact us for any questions regarding the method !

eeckhoutte.alexandre@gmail.com <br/>
tatiana.popova@curie.fr <br/>
marc-henri.stern@curie.fr <br/>

**Reach out to Marc-Henri Stern to ask about obtaining version 2.0.**

## Publications
*shallowHRD* publication :

Alexandre Eeckhoutte, Alexandre Houy, Elodie Manié, Manon Reverdy, Ivan Bièche, Elisabetta Marangoni, Oumou Goundiam, Anne Vincent-Salomon, Dominique Stoppa-Lyonnet, François-Clément Bidard, [...]

*shallowHRDv2* publication :

Celine Callens, Manuel Rodrigues, Adrien Briaux, Eleonore Frouin, Alexandre Eeckhoutte, Eric Pujade-Lauraine, Victor Renault, Dominique Stoppa-Lyonnet, Ivan Bieche, Guillaume Bataillon, Lucie Kar[...]
