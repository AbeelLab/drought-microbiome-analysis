# A pipeline to process amplicon data for the root and soil microbiome

This repository is associated with the following manuscript:

... (put citation here)

A schematic of the pipeline is shown below. We used QIIME2 to download and process amplicon sequencing data up to taxonomic profiles at genus level, and performed downstream analysis in Python (and R). Processed datasets used in our analysis are available in the ```datasets/``` folder.

[to do: add overview figure]

## Installation

[to do: conda environment]

## QIIME2 pipeline example usage

We use QIIME2 with the q2-fondue plugin to process amplicon datasets from NCBI. You can pass an accession or an existing study directory as input arguments.
If a study directory with an existing `accession.tsv` file already exists:
```bash
cd bash_scripts
./pipeline.sh --study_id study_name
```

If a study directory does not exist, you can provide a study accession and a study name, which will create a study fodler with an accession file:
```bash
./pipeline.sh --study_id study_name --accession EXAMPLE123456,EXAMPLE234567
```

To run only selected steps of the pipeline:
```bash
./pipeline.sh --study_id study_name --accession EXAMPLE123456,EXAMPLE234567 --run_download --run_cutadapt
```

For the full list of pipeline steps, options and flags, run:
```bash
./pipeline.sh --help
```

## Signature analysis and figure reproduction

You can download processed DADA2 feature tables here: [to do]

To process feature tables and obtain the merged dataset, you can run:
```bash
cd python_scripts && python3 process_datasets.py
```

To reproduce the analysis and associated figures in the paper, you can run:
```bash
python3 plot_figure_[figure_number].py
```
