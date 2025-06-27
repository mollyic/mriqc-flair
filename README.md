# FLAIR support for MRIQC

This repository extends the standard [MRIQC](https://github.com/poldracklab/mriqc) framework to support quality control and conversion workflows for 3D **FLAIR** MRI sequences.

## Issues, comments queries
Please log an issue on github or contact me at `molly.ireland@florey.edu.au`

## Version note:
The current main repository is a working version of the ongoing merge with the most recent MRIQC pipeline version 24. To access the original code used in the OHBM poster, please checkout on the branch `feature/mriqc-flair_v23.1.0`:

```
git checkout feature/mriqc-flair_v23.1.0
```

Please note that the current workflow with the most recent MRQIC pipeline has not yet been tested as extensively as the original pipeline. Results from this testing will be made available as this arises and documented here.

### Change log: mriqc-flair_v24
- **Likely change:** Reversion to MRIQC SpatialNormalization function given new ANTs settings files are outputting better segmentations


## Overview
As a common scan in clinical and research settings, this project was undertaken to extend the quantitative assessments of scan quality using the MRIQC framework for FLAIR. This project introduces:

- **FLAIR-specific image handling** added to core workflows
- Inclusion of FLAIR images in outputted MRQIC **quality control (QC) reports**

## OHBM Poster 2025

This work was presented at the Organization for Human Brain Mapping (OHBM) Annual Meeting.  
A pdf of this poster can be found in the [here](docs/mriqc-flair/OHBM_eposter_mriqcflair.pdf).

![Poster preview](docs/mriqc-flair/OHBM_eposter_mriqcflair.png)

## How to Use

To run MRIQC with FLAIR support:

```bash
python3 -m mriqc.cli.run /path/to/bids_dataset /path/to/output participant \
  --modalities FLAIR \
  --no-sub
```

*Note:* Files must be BIDS compliant e.g. `*_FLAIR.nii.gz`


## Acknowledgements

This extension builds on the incredible work of the [MRIQC team](https://mriqc.readthedocs.io). FLAIR support was inspired by QC needs for the growing multisite cohort in the Australian Epilepsy Project.

## License

MIT License. See `LICENSE` for details. Updates made to `NOTICE` file.
