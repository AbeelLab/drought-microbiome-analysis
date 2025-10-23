# A pipeline to process amplicon data for the root and soil microbiome

This repository is associated with the following manuscript:

... (put citation here)

## QIIME2 pipeline example usage

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

To do.
