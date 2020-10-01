# LSPpred
Random Forest prediction model for leaderless secretory proteins, based on features from curated lists of positive training data.

## Overview

`LSPpred` has two sub-modules: 

- The main module LSPpred which is based on a curated list of likely unconventionally secreted proteins in Arabidopsis 

- A secondary module SPLpred which is SecretomeP-like, in that it is based on classically secreted proteins with their signal peptide removed - notably from Arabidopsis, as opposed to the human data in the eukaryotic version of SecretomeP

Result of both modules are included in the .CSV output file

Each module has a built in cutoff to call a true-positive LSP with an estimated false-positive rate of 0.05.

The output file contains the following columns:

- Sequence - sequence ID  
- LSPpred_probability - model predicated probability of LSP status
- LSPpred - TRUE if LSPpred_probability exceeds builtin cutoff
- SPLpred_probability - - model predicated probability of LSP status
- SPLpred - TRUE if SPLpred_probability exceeds builtin cutoff
- Either - TRUE if either of LSPpred and SPLpred are TRUE
- Consensus - TRUE if both  of LSPpred and SPLpred are TRUE

Optionally the  `--low` flag can be used to add low confidence results for each module, where `LSPpred_probability` or `SPLpred_probability` exceeds 0.5

## Usage

```
usage: python lsppred.py [-h] [--version] [--low] [--output OUT_FILE]
                  [--log LOG_FILE]
                  FASTA_FILE

Predict LSPs in a FASTA file

positional arguments:
  FASTA_FILE         Input FASTA protein file

optional arguments:
  -h, --help         show this help message and exit
  --version          show program's version number and exit
  --low              Add low confidence predictions
  --output OUT_FILE  save output in CSV format to OUT_FILE
  --log LOG_FILE     record program progress in LOG_FILE

```


## Citation

```
Definitely Maybe: Predicting leaderless secretory proteins in plants without a gold standard 
Andrew Lonsdale, Melissa J. Davis, Monika S. Doblin and Tony Bacic. 
In Preparation.
```


