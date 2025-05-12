# Figure 2 Analysis

This repository contains code for generating the analysis and plots used in **Figure 2** of Houghton et al 2024 on BioRxIV. The notebook processes protein sequences, extracts regions of interest (e.g. transmembrane domains) and analyzes associations with observed subcellular localization.

## Contents

- `Figure_2_analysis.ipynb`: Jupyter notebook.
- `Figure_2_analysis.py`: Python script version of the notebook.
- `input/`: Folder containing input FASTA and CSV files.
- `output/`: Folder where results will be saved (auto-created if not present).

## Requirements

This project uses Python and requires the following packages:

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `biopython`
- `peptides`

You can install dependencies via:

```bash
pip install -r requirements.txt
```

## Usage

Run the notebook or script in a Python environment with the necessary dependencies. The script will automatically create the `output/` directory if it does not exist.

```bash
# Run the notebook interactively
jupyter notebook Figure_2_analysis.ipynb

# Or execute the script
python Figure_2_analysis.py
```

## License

This code is shared for academic purposes. Please cite the Houghton et al. bioRxiv 2024.08.28.610198 if used.

