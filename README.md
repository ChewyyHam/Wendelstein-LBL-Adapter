# LBL-Wendelstein

A lightweight adapter and analysis toolkit for running the Line-by-Line (LBL) radial-velocity pipeline on Wendelstein/FOCES spectra.

---

## Overview

This repository provides two lightweight tools for working with reduced Wendelstein/FOCES spectra within the LBL framework.

- **`wendelstein_lbl_adapter.py`** converts reduced Wendelstein `*_ods.fits` spectra into FITS files compatible with the LBL Generic Instrument.
- **`lbl_analysis.py`** provides interactive visualization and quick-look analysis of the resulting LBL radial-velocity (`.rdb`) products.

LBL itself is **not included** in this repository and should be installed separately following the official documentation.

---

## Workflow

```text
Reduced Wendelstein spectra
            |
            v
wendelstein_lbl_adapter.py
            |
            v
LBL Generic Instrument
            |
            v
lbl_*.rdb / lbl2_*.rdb
            |
            v
lbl_analysis.py
```

---

## Repository Contents

### `wendelstein_lbl_adapter.py`

- Interactive configuration
- Inclusive physical-order selection
- Automatic alignment of missing orders
- SNR estimation using `DER_SNR`
- Interactive target-coordinate handling
- Science and blaze FITS generation

### `lbl_analysis.py`

- Relative RV time series
- Lomb–Scargle periodograms
- Phase-folded RV curves
- Batch processing of `.rdb` files

---

## Requirements

- Python 3.10+
- NumPy
- Astropy
- Pandas
- Matplotlib
- LBL

```bash
pip install numpy astropy pandas matplotlib
```

LBL should be installed separately following the official documentation.

---

# Quick Start

Run

```bash
python wendelstein_lbl_adapter.py
```

and

```bash
python lbl_analysis.py
```

See **MANUAL.md** for detailed configuration, troubleshooting and implementation notes.

---

## Citation

If this software contributes to your work, please cite the original LBL publication together with the relevant Wendelstein/FOCES references.
