#!/bin/bash

cat <<EOF
Usage: pipeline.sh [options]

Required arguments:
  --study_id STUDY_ID
	      Unique identifier for the study. A sub-directory should exist under data: ${data_path}/study_id

Optional arguments:
  --accession ACCESSION
	      NCBI accession(s) for the study, compatible with q2-fondue.
  	      Provide this if ${data_path}/study_id/accession.tsv does not exist.
  	      If multiple accessions are provided, separate with commas (i.e. PRJNA435634,PRJNA435643). 
  	      This should be a .tsv file with contents of the form:
	      id
	      accession_1
	      accession_2

  --data_path PATH
	      Path to the data directory (default: "../data").

  --classifier_path PATH
	      Path to the QIIME2 artifact classifier file (default: "../silva-138-99-nb-classifier.qza").
  	      The pre-trained SILVA QIIME2 classifier for full-length sequences is available via:
	      https://resources.qiime2.org/

You can run the following pipeline steps. If no steps are specified, the entire pipeline is run.
Otherwise, only the specified steps are run.
  --run_download true|false
  	      Download data from NCBI in batches with QIIME q2-fondue plugin.
	      Batch size: 75. (TO DO: make this tweakable)

  --run_cutadapt true|false
  	      Remove primers with QIIME Cutadapt plugin.
	      We use default Cutadapt parameters.

  --run_figaro true|false
	      Extract trimming parameters using FIGARO.
	      We use default FIGARO parameters.

  --run_dada2 true|false
  	      Perform denoising and ASV detection with QIIME DADA2 plugin.
	      We use DADA2 defaults, with truncation lengths determine dby FIGARO.
	      These are stored in a file: ${data_path}/study_id/figaro_trim_params.txt with contents of the form:
	      forward_length
	      reverse_length

	      (TO DO: this could be done better if the user can specify truncation parameters as well.)

  --run_taxonomic_profiling
  	      Run taxonomic classification using QIIME2 pre-trained classifier.
	      This is a QIIME2 artifact specified by --classifier_path.

Other:
  --help      Show this help message and exit.

Example usage:
  If a study directory with accession.tsv already exists:
  ./pipeline.sh --study_id xu2018drought

  If a study directory does not exist:
  ./pipeline.sh --study_id xu2018drought --accession PRJNA435634,PRJNA435643

  To run only some steps of the pipeline:
  ./pipeline.sh --study_id xu2018drought --accession PRJNA435634,PRJNA435643 --run_download true --run_cutadapt true
EOF
