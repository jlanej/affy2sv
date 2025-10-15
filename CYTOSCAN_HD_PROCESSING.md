# CytoScan HD Array Processing Guide

This document describes how **Affymetrix CytoScan HD arrays** (and CytoScan 750K arrays) are processed by the `affy2sv` package to extract usable output including **genotypes**, **BAF (B Allele Frequency)**, **LRR (Log R Ratio)**, and **X/Y probe intensities (A/B signals)**.

## Overview

The processing workflow consists of three main steps:

1. **Raw CEL file processing** using Affymetrix Power Tools (APT)
2. **Data parsing and extraction** from APT output files
3. **Output generation** in various formats (MAD/R-GADA, PennCNV, SnpMatrix, PLINK)

## Prerequisites

### Required Files

Before processing CytoScan HD arrays, you need the following files from Affymetrix:

1. **Raw CEL files** - Raw intensity data from the array scanner
2. **Library files** (typically in a directory like `lib_cytoHD/`):
   - `CytoScanHD_Array.cdf` - Chip Definition File
   - `CytoScanHD_Array.chrXprobes` - X chromosome probe locations
   - `CytoScanHD_Array.chrYprobes` - Y chromosome probe locations
   - `CytoScanHD_Array.r1.qca` - QC analysis file
   - `CytoScanHD_Array.r1.qcc` - QC control file
   - `CytoScanHD_Array.snplist.txt` - SNP probe list
   - `CytoScanHD_Array.na32.3.annot.db` - Annotation database
   - `CytoScanHD_Array.na32.3.v1.REF_MODEL` - Reference model
   - `CytoScanHD_Array.na32.3.annot.csv` - Annotation in CSV format (for downstream processing)

### Software Requirements

- R (>= 3.0.0)
- Python 2.7 with pandas and numpy
- Affymetrix Power Tools (APT) binaries (included in the package)

## Step 1: Raw CEL File Processing with APT

The first step uses the `Cyto2APT` function to process raw CEL files using Affymetrix Power Tools.

### Function: `Cyto2APT()`

This function runs the `apt-copynumber-cyto` command-line tool to process CEL files and generate `.cychp.txt` files.

**Input:** Raw `.CEL` files from the scanner  
**Output:** `.cychp.txt` files containing probe-level data

### Example Usage

```r
library(affy2sv)

# Create APT parameter object
aptParam <- APTparam(
  type = "cytoscan", 
  level = "standard",
  cel.list = "/path/to/cel/files/", 
  output.path = "/path/to/apt/output/", 
  analysis.path = "/path/to/lib_cytoHD/",
  cdf = "CytoScanHD_Array.cdf", 
  chrX = "CytoScanHD_Array.chrXprobes", 
  chrY = "CytoScanHD_Array.chrYprobes", 
  qca = "CytoScanHD_Array.r1.qca", 
  qcc = "CytoScanHD_Array.r1.qcc", 
  snp = "CytoScanHD_Array.snplist.txt", 
  annot.db = "CytoScanHD_Array.na32.3.annot.db", 
  refmodel = "CytoScanHD_Array.na32.3.v1.REF_MODEL"
)

# Run APT processing
Cyto2APT(aptParam, verbose = TRUE)
```

### What APT Does

The APT tool performs:
- **Normalization** of raw probe intensities
- **Copy number analysis** to calculate Log2 ratios
- **Genotype calling** (AA, AB, BB, or NoCall)
- **Allele-specific signal extraction** (A and B probe intensities)
- **Quality control** metrics calculation

### APT Output: `.cychp.txt` Files

The `.cychp.txt` files contain multiple data tables (sections), including:

- **cnds** (Copy Number Data Set): Contains Log2Ratio, WeightedLog2Ratio, SmoothSignal for each probe
- **genotype**: Contains genotype calls (AA/AB/BB/NC), confidence scores, and **A and B signal intensities**
- **apds** (Allele Peaks Data Set): Contains AllelePeaks information for BAF calculation

## Step 2: Data Parsing and Extraction

The `.cychp.txt` files are parsed using a Python script (`parser.py`) to extract relevant data tables.

### Python Parser

The parser extracts the following tables from each `.cychp.txt` file:

1. **genotype** table:
   - ProbeSetName (probe ID)
   - Call (genotype: AA, AB, BB, NC)
   - Confidence
   - ForcedCall
   - **ASignal** - A allele probe intensity
   - **BSignal** - B allele probe intensity
   - SignalStrength
   - Contrast

2. **cnds** table:
   - ProbeSetName
   - Chromosome
   - Position
   - **Log2Ratio** (LRR equivalent)
   - WeightedLog2Ratio
   - SmoothSignal

3. **apds** table:
   - ProbeSetName
   - Chromosome
   - Position
   - AllelePeaks0
   - AllelePeaks1

The parser filters and aligns these tables to ensure all three contain data for the same set of probes, creating filtered output files (`.genotype.f.txt`, `.cnds.f.txt`, `.apds.f.txt`).

## Step 3: Output Generation

The parsed data can be converted into multiple output formats depending on your downstream analysis needs.

### Option A: MAD/R-GADA Format (Includes Genotypes, BAF, and LRR)

Function: `Cyto2Mad()`

Creates files compatible with the MAD and R-GADA R packages for CNV detection.

**Output format:** Tab-delimited files with columns:
- `Name` - Probe Set ID
- `Chr` - Chromosome
- `Position` - Physical position
- `Log.R.Ratio` - Copy number log ratio (LRR)
- `GType` - Genotype call (AA=0, AB=1, BB=2, NC=NA)
- `B.Allele.Freq` - B Allele Frequency (BAF)

#### Example Usage

```r
# Create MAD format files (one file per sample)
Cyto2Mad(
  cychp.files = "/path/to/apt/output/",
  output.name = "/path/to/mad/output/",
  annotation.file = "/path/to/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
  output.type = "mad",
  cychp.attime = 4,  # Process 4 files in parallel
  verbose = TRUE
)
```

#### BAF Calculation Details

The B Allele Frequency (BAF) is calculated from the `AllelePeaks0` values in the `apds` table:

1. **Quantile normalization** is applied to AllelePeaks0:
   - Bottom quantile (default 0.25) and top quantile (default 0.75) are calculated
   - Values below bottom quantile are set to bottom quantile
   - Values above top quantile are set to top quantile

2. **Normalization to [0, 1] range:**
   ```
   normalized = (AllelePeaks0 - bottom_quantile) / (top_quantile - bottom_quantile)
   ```

3. **Inversion** to get BAF:
   ```
   BAF = 1 - normalized
   ```

This results in BAF values where:
- BAF ≈ 0 indicates AA genotype (mostly A allele)
- BAF ≈ 0.5 indicates AB genotype (heterozygous)
- BAF ≈ 1 indicates BB genotype (mostly B allele)

### Option B: PennCNV Format (Includes Genotypes, BAF, and LRR)

Function: `Cyto2Mad()` with `output.type = "penncnv"`

Creates a single file compatible with PennCNV for CNV detection.

**Output format:** Tab-delimited file with columns for each sample:
- `Name` - Probe Set ID
- `Chr` - Chromosome
- `Position` - Physical position
- `SampleID.Log R Ratio` - LRR for sample
- `SampleID.GType` - Genotype for sample (0/1/2/NA)
- `SampleID.B Allele Freq` - BAF for sample

#### Example Usage

```r
# Create PennCNV format file (all samples in one file)
Cyto2Mad(
  cychp.files = "/path/to/apt/output/",
  output.name = "/path/to/output/dataset",  # Will create dataset.penncnv
  annotation.file = "/path/to/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
  output.type = "penncnv",
  cychp.attime = 4,
  verbose = TRUE
)
```

### Option C: SnpMatrix Format (Genotypes Only)

Function: `Cyto2SnpMatrix()`

Creates a SnpMatrix object (from the `snpStats` package) containing genotypes.

**Output:** An R list with two components:
- `map` - Data frame with probe annotation
- `genotype` - SnpMatrix object with genotype calls

#### Example Usage

```r
# Create SnpMatrix container
smc <- Cyto2SnpMatrix(
  cychp.files = "/path/to/apt/output/", 
  annotation.file = "/path/to/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
  output.type = "snpmatrix",
  cychp.attime = 4,
  verbose = TRUE
)

# Access genotypes
genotypes <- smc$genotype  # SnpMatrix object
map_info <- smc$map        # Probe annotation
```

### Option D: PLINK Format (Genotypes Only)

Function: `Cyto2SnpMatrix()` with `output.type = "plink"`

Creates a `.tped` file compatible with PLINK.

**Output:** A transposed pedigree file (`.tped`) with alleles in forward strand orientation.

#### Example Usage

```r
# Create PLINK .tped file
Cyto2SnpMatrix(
  cychp.files = "/path/to/apt/output/", 
  annotation.file = "/path/to/lib_cytoHD/CytoScanHD_Array.na32.3.annot.csv", 
  output.type = "plink",
  output.name = "/path/to/output/dataset",  # Will create dataset.tped
  cychp.attime = 4,
  verbose = TRUE
)

# Note: You need to create the .tfam file separately
```

## Extracting A/B Probe Intensities (X/Y Signals)

The **A and B signal intensities** for each probe are available in the **genotype** table of the `.cychp.txt` files. These raw signal intensities represent the hybridization signal for the A and B alleles at each probe location.

### Accessing A/B Signal Intensities

To extract the raw A and B signal intensities, you need to:

1. **Run APT processing** to generate `.cychp.txt` files
2. **Use the Python parser** to extract the genotype table
3. **Read the genotype.f.txt files** which contain the ASignal and BSignal columns

### Custom Extraction Example

```r
library(data.table)

# Parse a cychp.txt file
system(paste0("python ", 
              system.file("exec/parser.py", package="affy2sv"),
              " -i /path/to/sample.cyhd.cychp.txt",
              " -o /path/to/output/"))

# Read the genotype data with A/B signals
genotype_data <- fread("/path/to/output/sample.cyhd.cychp.genotype.f.txt")

# The genotype_data table contains:
# - ProbeSetName: Probe identifier
# - Call: Genotype call (AA/AB/BB/NC)
# - Confidence: Call confidence score
# - ForcedCall: Alternative call
# - ASignal: A allele intensity (X signal)
# - BSignal: B allele intensity (Y signal)
# - SignalStrength: Overall signal
# - Contrast: Contrast metric

# Extract A and B signals
a_intensities <- genotype_data$ASignal
b_intensities <- genotype_data$BSignal
probe_names <- genotype_data$ProbeSetName

# Create a data frame with signals
signal_data <- data.frame(
  Probe = probe_names,
  A_Intensity = as.numeric(a_intensities),
  B_Intensity = as.numeric(b_intensities),
  Genotype = genotype_data$Call
)
```

### A/B Signal Characteristics

- **ASignal** and **BSignal** are raw, normalized intensity values
- For AA genotypes: ASignal is high, BSignal is low
- For AB genotypes: Both ASignal and BSignal are intermediate
- For BB genotypes: ASignal is low, BSignal is high
- These signals can be used for:
  - Custom BAF calculation
  - Quality control
  - Allele-specific copy number analysis
  - Custom genotype calling algorithms

## Data Flow Summary

```
Raw CEL Files
      ↓
[Cyto2APT - APT Processing]
      ↓
.cychp.txt files (contains all data)
      ↓
[Python Parser - Table Extraction]
      ↓
Filtered tables:
├── genotype.f.txt (Genotypes + A/B Signals)
├── cnds.f.txt (Log2Ratio = LRR)
└── apds.f.txt (AllelePeaks for BAF)
      ↓
[Output Generation Functions]
      ↓
Final Output:
├── MAD/R-GADA: Individual files (LRR + BAF + Genotypes)
├── PennCNV: Single file (LRR + BAF + Genotypes)
├── SnpMatrix: R object (Genotypes)
└── PLINK: .tped file (Genotypes)
```

## Complete Workflow Example

Here's a complete example processing CytoScan HD arrays from raw CEL files to MAD output:

```r
library(affy2sv)

# Step 1: Set up paths
cel_path <- "/data/cytoscan/cel_files/"
lib_path <- "/data/cytoscan/lib_cytoHD/"
apt_output <- "/data/cytoscan/apt_output/"
mad_output <- "/data/cytoscan/mad_output/"

# Step 2: Configure APT parameters
aptParam <- APTparam(
  type = "cytoscan", 
  level = "standard",
  cel.list = cel_path, 
  output.path = apt_output, 
  analysis.path = lib_path,
  cdf = "CytoScanHD_Array.cdf", 
  chrX = "CytoScanHD_Array.chrXprobes", 
  chrY = "CytoScanHD_Array.chrYprobes", 
  qca = "CytoScanHD_Array.r1.qca", 
  qcc = "CytoScanHD_Array.r1.qcc", 
  snp = "CytoScanHD_Array.snplist.txt", 
  annot.db = "CytoScanHD_Array.na32.3.annot.db", 
  refmodel = "CytoScanHD_Array.na32.3.v1.REF_MODEL"
)

# Step 3: Run APT processing
Cyto2APT(aptParam, verbose = TRUE)

# Step 4: Generate MAD output (includes LRR, BAF, Genotypes)
Cyto2Mad(
  cychp.files = apt_output,
  output.name = mad_output,
  annotation.file = file.path(lib_path, "CytoScanHD_Array.na32.3.annot.csv"),
  output.type = "mad",
  cychp.attime = 4,
  bottom.quantile = 0.25,
  top.quantile = 0.75,
  verbose = TRUE
)

# Step 5: (Optional) Also create PennCNV format
Cyto2Mad(
  cychp.files = apt_output,
  output.name = "/data/cytoscan/dataset",
  annotation.file = file.path(lib_path, "CytoScanHD_Array.na32.3.annot.csv"),
  output.type = "penncnv",
  cychp.attime = 4,
  verbose = TRUE
)

# Step 6: (Optional) Create SnpMatrix for genotype analysis
smc <- Cyto2SnpMatrix(
  cychp.files = apt_output, 
  annotation.file = file.path(lib_path, "CytoScanHD_Array.na32.3.annot.csv"), 
  output.type = "snpmatrix",
  cychp.attime = 4,
  verbose = TRUE
)
```

## Output File Descriptions

### MAD Format Output

Each sample gets an individual file with the following structure:

```
Name            Chr  Position    Log.R.Ratio  GType  B.Allele.Freq
CN_001231       1    742429      -0.0523      AB     0.4892
CN_001232       1    744210      0.1234       BB     0.9876
CN_001233       1    745123      -0.0891      AA     0.0234
...
```

### PennCNV Format Output

Single file with all samples:

```
Name        Chr  Position  Sample1.Log R Ratio  Sample1.GType  Sample1.B Allele Freq  Sample2.Log R Ratio  ...
CN_001231   1    742429    -0.0523              1              0.4892                 0.0234               ...
CN_001232   1    744210    0.1234               2              0.9876                 0.8765               ...
...
```

## Key Metrics Explained

### Log R Ratio (LRR)
- **Source:** Log2Ratio column from cnds table
- **Meaning:** Relative copy number in log2 scale
- **Interpretation:**
  - LRR ≈ 0: Normal copy number (2 copies)
  - LRR > 0: Copy number gain
  - LRR < 0: Copy number loss
- **Use:** CNV detection, copy number analysis

### B Allele Frequency (BAF)
- **Source:** Calculated from AllelePeaks0 in apds table
- **Meaning:** Proportion of B allele signal
- **Interpretation:**
  - BAF ≈ 0: Homozygous AA
  - BAF ≈ 0.5: Heterozygous AB
  - BAF ≈ 1: Homozygous BB
- **Use:** Allelic imbalance detection, LOH analysis, CNV genotyping

### Genotype (GType)
- **Source:** Call column from genotype table
- **Values:** AA (or 0), AB (or 1), BB (or 2), NC/NA (No Call)
- **Use:** SNP genotyping, association studies, population genetics

### A/B Signal Intensities
- **Source:** ASignal and BSignal columns from genotype table
- **Meaning:** Raw probe hybridization intensities for each allele
- **Use:** Custom normalization, quality control, allele-specific analysis

## Quality Control

The APT processing includes quality metrics that can be reviewed:

1. Check the QC output files generated by APT
2. Examine call rates (proportion of non-NC genotypes)
3. Look for outlier samples in LRR and BAF distributions
4. Use `CytoQCView()` function for visualization (if available)

## Troubleshooting

### Common Issues

1. **Missing library files:** Ensure all required files from Affymetrix are present
2. **Python errors during parsing:** Verify Python 2.7 and pandas are installed
3. **Memory issues:** Reduce `cychp.attime` parameter to process fewer samples in parallel
4. **Path issues:** Use absolute paths or ensure working directory is correct

## References

- [affy2sv paper](http://www.biomedcentral.com/1471-2105/16/167): Hernandez-Ferrer et al., BMC Bioinformatics 2015, 16:167
- Affymetrix CytoScan HD Array documentation
- [Affymetrix Power Tools (APT) documentation](https://www.thermofisher.com/affymetrix)

## Citation

If you use this workflow, please cite:

> Carles Hernandez-Ferrer, Ines Quintela Garcia, Katharina Danielski, Ángel Carracedo, Luis A. Pérez-Jurado and Juan R. González. "affy2sv: an R package to pre-process Affymetrix CytoScan HD and 750K arrays for SNP, CNV, inversion and mosaicism calling." BMC Bioinformatics 2015, 16:167. doi:10.1186/s12859-015-0608-y
