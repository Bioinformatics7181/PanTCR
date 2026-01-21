# PanTCR: Robust Inference of Novel TCR Alleles from Fragmented RNA-seq Data via Pan-clonotype Graph Priors

This repository contains the implementation and reproduction pipeline for **PanTCR**, a tool designed to infer novel TCR alleles from fragmented RNA-seq data using pan-clonotype graph priors.

---

## Prerequisites

Before running the pipeline, ensure your system meets the following requirements:

- **OS**: Linux or macOS  
- **Python**: 3.8 or higher  
- **Conda**: Recommended for environment management  
- **Java**: Required for MiXCR (v4.0 or higher recommended)

---

## Installation & Setup

### 1. Create and Activate Environment

We recommend using Conda to create an isolated environment to avoid dependency conflicts.

```bash
# Create a new conda environment named 'pantcr'
conda create -n pantcr python=3.9 -y

# Activate the environment
conda activate pantcr
```

---

### 2. Install Python Dependencies

Install the required Python libraries (mainly `numpy` and `pandas`).

```bash
pip install -r requirements.txt
```

If you do not yet have a `requirements.txt`, create one with the following content:

```txt
numpy>=1.21.0
pandas>=1.3.0
openpyxl>=3.0.0
```

---

### 3. Install External Bioinformatics Tools

This pipeline relies on **MiXCR** (for alignment and allele calling) and **ART** (for NGS simulation).

#### ART (Illumina Simulator)

```bash
conda install -c bioconda art
```

#### MiXCR (Version 4.0+)

MiXCR must be installed and available in your system `PATH`.

- **Option A (Official Zip)**:  
  Download from the MiXCR website, unzip, and add it to your `PATH`.

- **Option B (Conda – verify version)**:  
  ```bash
  conda install -c bioconda mixcr
  ```
  Ensure the installed version is **v4.0+**. Otherwise, manual installation is required.

#### Verification

```bash
art_illumina --help
mixcr --version   # Must be >= 4.0
```

---

## Data Preparation

Before running the simulation, unzip the reference CDR3 dictionary.

```bash
unzip simu/ref/cdr3_dict.zip -d simu/ref/
```

---

## Usage

The repository includes a comprehensive script to reproduce the experiments described in the paper.

### Reproducing 5-Fold Cross-Validation (Default)

```bash
# Enter the simulation directory
cd simu

# Grant execution permissions (if necessary)
chmod +x simu_validation_pipeline.sh

# Run the pipeline
./simu_validation_pipeline.sh
```

### Customizing Experiments

You can modify experimental parameters (e.g., read length, depth, population, degradation levels) by:

- Editing the script `simu_validation_pipeline.sh`
- Passing command-line arguments if supported, e.g.:
  ```bash
  ./simu_validation_pipeline.sh --expr-id expr_1 --p-deg 0.9
  ```

Refer to the script comments for supported options.

---

## Output

Results are generated in the `simu/samples` directory (or a configured output path), organized by experiment ID and population.

Log files are available in:

```text
simu/logs
```

---

## License

MIT License: only for acdamic use.

---

## Contact

Xinyang Qian: qianxy@stu.xjtu.edu.cn
