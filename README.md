# Chemical Curation

This repository provides a Python script for parsing and curation of chemical data. The curation process follows the procedure recommended in the following paper:

https://pubs.acs.org/doi/10.1021/ci100176x

Eventually, this repository will include the procedures recommended in the following:

https://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00129

## Curation process

### 1. Data aggregation
  - sdf and csv files are supported, some basic searching for desired fields is implemented
### 2. Removal of mixtures, inorganics, and salts
  - currently using MolVS Standardizer
### 3. Structure standardization
  - currently using MolVS Standardizer
### 4. Duplicate averaging/removal
  - InChI keys are used to determine if structures are the same
  - Users can provide the threshold for determining if activities for the same molecule are close enough to average, default threshold is one log unit
### 5. Manual review
  - Any issues throughout the process are stored in a text file for manual review


