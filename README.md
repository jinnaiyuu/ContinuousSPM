# KDD2017 codes
Codes to reproduce results of our KDD2017 paper:
**Significant Pattern Mining on Continuous Variables** by Mahito Sugiyama (Osaka Univ.) and Karsten Borgwadt (ETH Zürich).

The proposed method with the Bonferroni correction version is implemented in C++ and its source code is `src/cc/main.cc`.  
The binarization-based methods are composed of two parts: Binarization is implemented in Python 3 and the code is `src/python/binarize.py`, and significant itemset mining is implemented in C and codes are in the directory `src/c`, which is obtained as a part of [Westfall-Young light](https://www.bsse.ethz.ch/mlcb/research/machine-learning/wylight.html).

You can reproduce results in our paper by the R script `src/R/script_for_paper.R`.

## Usege:

First, compile C and C++ codes:
```
$ cd src/cc
$ make
$ cd ../c/lcm_lamp_fisher
$ make
$ cd ../lcm_comp_pvalues_fisher
$ make
$ cd ../../../
```
The Boost libary is required to compile the C++ code.
Then go to the R directory and start R:
```
$ cd src/R
$ R
```
In the R environment, load the script file `script_for_paper.R` and run the function `run.script` with appropriate arguments:
```
> source("script_for_paper.R")
> res <- run.script()
  > Generate synthetic data
  > Perform the proposed mehtod
    Compute precision and recall
  > Perform the Bonferroni correction mehtod
  > Perform median discretization + significant itemset mining
    Binarize data
    Perform significant itemset mining to count the number of testable patterns
    Enumerate significant patterns
    Compute precision and recall
> res
                       Proposal   Bonferroni Median-bin
#testable patterns 3661.0000000 1.048580e+06 1.5547e+04
Runtime               0.0280180 2.232519e+00 8.6458e-02
Final freq. thr.      0.0649643 0.000000e+00 1.9000e+01
Precision             0.1013514           NA 1.3930e-03
Recall                0.9375000           NA 5.0000e-01
```
### Arguments for `run.script()`
* `N`: The number of samples for synthetic data (default: `N = 200`)
* `n`: The number of features for synthetic data (default: `n = 20`)
* `r0`: Minor class ratio (default: `r0 = 0.5`)
* `data.name`: Name of a dataset. If `data.name = "synth"`, a synthetic dataset is generated, and if it is one of the following names, the corresponding real data is used (default: `data.name = "synth"`)
  * `("ctg", "faults", "ionosphere", "segment", "waveform", "wdbc")`
* `base.dir`: A directory for data and results (default: `base.dir = "./"`)
* `alpha`: The FWER level (default: `alpha = 0.05`)
* `bin.type`: A type of binarization, which should be one of the following list (default: `bin.type = 1`)
  * `bin.type = 1`: median-based binarization
  * `bin.type = 2`: interordinal scaling
  * `bin.type = 3`: interval scaling
* `run.prop`: Run the proposed method if `TRUE` (default: `run.prop = TRUE`)
* `run.bonf`: Run the bonferroni correction method if `TRUE` (default: `run.bonf = TRUE`)
* `run.bin`: Run the binarization-based method if `TRUE` (default: `run.bin = TRUE`)

## Contact
Author: Mahito Sugiyama  
Affiliation: ISIR, Osaka University, Japan  
E-mail: mahito@ar.sanken.osaka-u.ac.jp
