# Raman Spectral Quantification Toolkit

[![Python
3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![NumPy](https://img.shields.io/badge/NumPy-1.19+-orange.svg)](https://numpy.org/)
[![SciPy](https://img.shields.io/badge/SciPy-1.5+-green.svg)](https://scipy.org/)

A Python toolkit for molecular quantification analysis using Raman spectra and Non-Negative Least Squares (NNLS) fitting.


## Overview

This toolkit consists of two main components:

1.  **`prepare_reference.py`** - Constructs a reference spectral matrix from individual compound spectra
2.  **`spectral_fitting.py`** - Performs NNLS fitting to quantify molecular components in complex spectra


## Installation

### Requirements

``` bash
pip install numpy pandas scipy matplotlib
```

### Clone Repository

``` bash
git clone https://github.com/MicrobeLab/raman-molecular-quantification.git
cd raman-molecular-quantification
```


## Step 1: Reference Matrix Construction

### `prepare_reference.py`

Builds a reference spectral matrix from individual compound files. `reference_library` contains ready-to-use spectra (converted to integrated area) for each compound. `compound_list.txt` is a text file with compound names (one per line). `ref_spec.txt` is a tab-separated table where each column represents one compound (in the order specified in `compound_list.txt`).

#### Usage

``` bash
python prepare_reference.py --input-dir reference_library  --ref-list compound_list.txt --output ref_spec.txt
```

## Step 2: Spectral Fitting and Quantification

### `spectral_fitting.py`

Performs NNLS fitting. `cell_spectrum.txt` is the path to spectra (converted to integrated area) from cells or other mixtures. `ref_spec.txt` is the path to reference matrix from Step 1. `compound_list.txt` is a text file with compound names (same as Step 1). `quantification.txt` is the output file for quantification coefficients (fitting weights). `residuals.txt` is the output file for residual analysis.

#### Usage

``` bash
python spectral_fitting.py --input cell_spectrum.txt --ref ref_spec.txt --ref-list compound_list.txt \
--quant_out quantification.txt --residual_out residuals.txt
```

## Notes
Wavenumber ranges should be consistent across all reference and input spectra.

## Example

An example dataset is provided in the `example` directory, demonstrating quantitative analysis of Raman spectra for *Bacillus licheniformis* cells under three conditions: control, citrate-treated, and spores.

## Other Resources
The workflow for Raman spectral pre-processing is available [here](https://github.com/MicrobeLab/scCulturePrec-data).
Bugs and difficulties in using the toolkit are welcome on [the issue tracker](https://github.com/MicrobeLab/raman-molecular-quantification/issues).

