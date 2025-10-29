# Datasets

This collection of drought datasets was compiled from 13 public amplicon 16S datasets, as cited below.

## Studies included in the meta-analysis

| Study                         | Study ID                      | DOI                                                                                  | Accession(s)             |
|-------------------------------|-------------------------------|--------------------------------------------------------------------------------------|--------------------------|
| **Drought datasets**          |                               |                                                                                      |                          |
| Xu et al. (2018)              | xu2018drought                 | [10.1073/pnas.1717308115](https://doi.org/10.1073/pnas.1717308115)                   | PRJNA435634, PRJNA435643 |
| Simmons et al. (2020)         | simmons2020drought            | [10.3389/fpls.2020.00599](https://doi.org/10.3389/fpls.2020.00599)                   | PRJNA607579              |
| Fitzpatrick et al. (2018)     | fitzpatrick2018assembly       | [10.1073/pnas.1717617115](https://doi.org/10.1073/pnas.1717617115)                   | PRJNA428446              |
| Azarbad et al. (2020)         | azarbad2020four               | [10.1093/femsec/fiaa098](https://pubmed.ncbi.nlm.nih.gov/32440671/)                  | PRJNA526458              |
| Santos-Medellin et al. (2021) | santos-medellin2021prolonged  | [10.1038/s41477-021-00967-1](https://www.nature.com/articles/s41477-021-00967-1)     | PRJNA551661              |
| Bandopadhyay et al. (2024)    | bandopadhyay2024disentangling | [10.1038/s41467-024-50463-1](https://www.nature.com/articles/s41467-024-50463-1)     | PRJNA862978              |
| Naylor et al. (2017)          | naylor2017drought             | [10.1038/ismej.2017.118](https://pubmed.ncbi.nlm.nih.gov/28753209/)                  | PRJNA369551              |
| Santos-Medellin et al. (2017) | santos-medellin2017drought    | [10.1128/mbio.00764-17](https://journals.asm.org/doi/10.1128/mbio.00764-17)          | PRJNA386367              |
| Azarbad et al. (2022)         | azarbad2022response           | [10.1038/s43705-022-00151-2](https://www.nature.com/articles/s43705-022-00151-2)     | PRJNA736197              |
| **Inoculation datasets**      |                               |                                                                                      |                          |
| Moore et al. (2023)           | moore2023microbial            | [10.1128/spectrum.01476-22](https://journals.asm.org/doi/10.1128/spectrum.01476-22)  | PRJNA780613              |
| Swift et al. (2025)           | swift2025drought              | [10.1007/s11104-024-06853-x](https://journals.asm.org/doi/10.1128/spectrum.01476-22) | PRJNA913622              |
| Zhang et al. (2022)           | zhang2022cross                | [10.1111/tpj.15775](https://onlinelibrary.wiley.com/doi/abs/10.1111/tpj.15775)       | PRJNA839620              |
| Munoz-Ucros et al. (2022)     | munoz-ucros2022drought        | [10.​1007/​s11104-​021-​05227-x](https://doi.org/10.1007/s11104-021-05227-x)             | PRJEB41348               |


## File and folder documentation
| File/folder                                                     | Description                                                                                                                                                                   |
|-----------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| counts_[taxonomic_resolution].tsv                               | Count table where the samples (NCBI accessions) are rows and the columns are features at a given taxonomic resolution, from phylum to genus. Obtained by merging all studies. |
| counts_[taxonomic_resolution]_for_signature.tsv                 | Count table obtained by merging all drought studies, and retaining only features occuring in 30% of studies (rounded up).                                                     |
| counts_[taxonomic_resolution]_for_signature_batch_corrected.tsv | Count table (as before) with batch correction for samples originating from different studies.                                                                                 |
| metadata.tsv                                                    | Metadata for all samples (NCBI accessions). See below for a more detailed description.                                                                                        |
| [study_id]/                                                     | Folder with all processed counts and metadata .tsv files for a given study.                                                                                                   |
| [compartment_name]/                                             | Folder with all processed datasets for a given compartment (endosphere, rhizosphere, bulks soil), obtained by merging samples across studies.                                 |


### Metadata columns
| Column          | Description                                                                                                                                          |
|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| #SampleID       | NCBI sample accession.                                                                                                                               |
| StudyID         | Study ID as given above.                                                                                                                             |
| StudyName       | Study name as given above.                                                                                                                           |
| Treatment       | Drought or Control.                                                                                                                                  |
| Host            | Plant host clade. Soil for bulk soil samples.                                                                                                        |
| HostSpecific    | Plant host. Soil for bulk soil samples.                                                                                                              |
| RootCompartment | Rhizosphere, Endosphere or Bulk soil.                                                                                                                |
| Inoculum        | Dry, Reference or Sterile. Only for inoculation studies. Check Supplementary Table S1 in the publication for information on how these were selected. |
| Primers         | Primers used for amplicon sequencing in the corresponding study, manually extracted from each study.                                                 |
	
