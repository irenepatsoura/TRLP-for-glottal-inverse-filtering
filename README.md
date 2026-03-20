# TRLP for Glottal Inverse Filtering

This repository collects MATLAB and Python code for glottal inverse filtering experiments and downstream classification of pathological speech.

The core goal is to compare three inverse-filtering approaches:

- IAIF
- QCP
- TRLP

and evaluate whether glottal-flow-derived features help distinguish healthy control and pathological speech recordings.

## Overview

The project has two main parts.

1. MATLAB analysis and feature extraction
   - Frame-based glottal inverse filtering on speech signals
   - Comparison against ground-truth glottal signals on synthetic data
   - Extraction of engineered glottal and acoustic descriptors into CSV files

2. Python classification
   - Nested speaker-grouped cross-validation
   - Comparison across SVM, Random Forest, and XGBoost
   - Export of summary tables and diagnostic plots

## End-to-End Workflow

1. Run MATLAB inverse-filtering experiments on reference signals in `data/Glottal_signals_db/`.
2. Use the MATLAB feature extraction pipeline to generate `features_*.csv` tables.
3. Run `classify_tasks.py` on those CSV files.
4. Review per-task and global classification summaries under `classification_results/`.

## Main Components

### MATLAB side

- `run_iaif.m`: frame-by-frame IAIF analysis on example speech signals
- `run_qcp.m`: frame-by-frame QCP analysis
- `run_trlp.m`: frame-by-frame TRLP analysis
- `run_metrics_analysis.m`: error analysis for NAQ, QOQ, HRF, and H1-H2 against ground truth
- `run_h1h2_analysis.m`: focused H1-H2 estimation error analysis
- `extract_features_all_tasks.m`: feature extraction from task-organized recordings into `features_*.csv`
- `mytrlp.m`: TRLP coefficient estimation routine

### Python side

- `classify_tasks.py`: classification entry point for all `features_*.csv` files in the repository root

The classifier script:

- expects precomputed CSV feature tables
- does not require MATLAB or PC-GITA once the CSVs already exist
- evaluates multiple feature subsets
- writes fold-level and summary metrics to timestamped output folders

## Requirements

### MATLAB requirements

- MATLAB with Signal Processing functionality used by the scripts
- access to the external datasets used in the experiments
- task-organized speech recordings for feature extraction

The feature extraction script currently expects task data under:

- `data/PC-GITA_per_task_44100Hz`

### Python requirements

To run `classify_tasks.py`, install:

- `numpy`
- `pandas`
- `matplotlib`
- `scikit-learn`
- `xgboost`

Example setup:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install numpy pandas matplotlib scikit-learn xgboost
python classify_tasks.py
```

## How to Run Classification

If you already have the exported `features_*.csv` files, the Python classification step can be run without MATLAB.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python classify_tasks.py
```

This will scan the repository root for `features_*.csv` files and write outputs under `classification_results/`.

## Expected Inputs and Outputs

### Inputs for classification

The Python classification step expects one or more files named like:

- `features_Vowels_IAIF.csv`
- `features_Vowels_QCP.csv`
- `features_Vowels_TRLP.csv`

Each CSV should include metadata columns such as:

- `speaker`
- `label`
- optionally `file_name` and `task`

plus numeric engineered feature columns.

### Outputs from classification

Running `classify_tasks.py` creates:

- `classification_results/run_<timestamp>/...`
- per-task `summary_metrics.csv`
- per-task `fold_metrics.csv`
- per-task plots such as `summary_barplots.png` and `summary_heatmap.png`
- run-level `global_summary_metrics.csv`
- run-level `global_fold_metrics.csv`

## Project Structure

```text
TRLP-for-glottal-inverse-filtering/
├── classify_tasks.py
├── extract_features_all_tasks.m
├── mytrlp.m
├── run_h1h2_analysis.m
├── run_iaif.m
├── run_metrics_analysis.m
├── run_qcp.m
├── run_trlp.m
├── presentation_slides.html
├── speech_classification_regression_guide.ipynb
├── SVM_vowels.ipynb
├── core/
│   ├── dap.m
│   ├── find_f0.m
│   ├── find_f0_yin.m
│   ├── gci.m
│   ├── iaif.m
│   ├── qcp.m
│   ├── VTfix.m
│   └── wlp.m
├── eval/
│   ├── compare_with_ground_truth.m
│   ├── compute_glottal_metrics.m
│   ├── compute_h1h2.m
│   ├── plot_filter_evolution.m
│   ├── visualize_frame_results.m
│   └── visualize_iaif_results.m
├── framework/
│   ├── at.m
│   ├── lpc_signal.m
│   ├── mergestruct.m
│   ├── signal.m
│   ├── time.m
│   ├── valid.m
│   └── win.m
├── pipeline/
│   ├── create_adaptive_frames.m
│   ├── create_fixed_frames.m
│   ├── estimate_pitch.m
│   └── reconstruct_signal.m
├── data/
│   └── Glottal_signals_db/
├── HOWTO_ENVELOPE/
└── HOWTO_IF/
```

## Reading List

This section is intended to replace local storage of publisher PDFs in `HOWTO_ENVELOPE/` and `HOWTO_IF/`. Prefer DOI landing pages, official conference pages, or other legal public links over redistributing PDFs.

1. Glottal wave analysis with pitch synchronous iterative adaptive inverse filtering
   - Authors: Paavo Alku
   - Year: 1991
   - Venue: 2nd European Conference on Speech Communication and Technology (Eurospeech 1991)
   - DOI/URL: https://doi.org/10.21437/Eurospeech.1991-257

2. Quasi closed phase analysis for glottal inverse filtering
   - Authors: Manu Airaksinen, Brad Story, Paavo Alku
   - Year: 2013
   - Venue: Interspeech 2013
   - DOI/URL: https://doi.org/10.21437/Interspeech.2013-55

3. Quasi Closed Phase Glottal Inverse Filtering Analysis With Weighted Linear Prediction
   - Authors: Manu Airaksinen, Tuomo Raitio, Brad Story, Paavo Alku
   - Year: 2014
   - Venue: IEEE/ACM Transactions on Audio, Speech, and Language Processing
   - DOI/URL: https://doi.org/10.1109/TASLP.2013.2294585

4. Glottal inverse filtering using stabilised weighted linear prediction
   - Authors: George P. Kafentzis, Yannis Stylianou, Paavo Alku
   - Year: 2011
   - Venue: 2011 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)
   - DOI/URL: https://doi.org/10.1109/ICASSP.2011.5947581

5. Glottal Inverse Filtering Using Probabilistic Weighted Linear Prediction
   - Authors: Achuth Rao M. V., Prasanta Kumar Ghosh
   - Year: 2019
   - Venue: IEEE/ACM Transactions on Audio, Speech, and Language Processing
   - DOI/URL: https://doi.org/10.1109/TASLP.2018.2873897

6. Regularized Linear Prediction of Speech
   - Authors: L. Anders Ekman, W. Bastiaan Kleijn, Manohar N. Murthi
   - Year: 2008
   - Venue: IEEE Transactions on Audio, Speech, and Language Processing
   - DOI/URL: https://doi.org/10.1109/TASL.2007.909448

7. DC-constrained linear prediction for glottal inverse filtering
   - Authors: Paavo Alku, Carlo Magi, Tom Backstrom
   - Year: 2008
   - Venue: Interspeech 2008
   - DOI/URL: https://doi.org/10.21437/Interspeech.2008-698

## License and Data Use

This repository includes an `All rights reserved` license in the `LICENSE` file.

Code, scripts, figures, and documentation in this repository may not be reused, modified, or redistributed without prior written permission from the copyright holder.

Third-party datasets, papers, and other externally sourced materials remain subject to their own licenses, copyright terms, and access restrictions.