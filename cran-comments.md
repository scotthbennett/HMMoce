## RESUBMISSION - 27 OCT 2017
This is a resubmission. This package differs significantly from package 'moveHMM' in that 'moveHMM' relies on animal tag data with negligible measurement error. The purpose of 'HMMoce' is specifically to address animal tag data with significant error. It does so (by necessity) in a significantly different way than how 'moveHMM' functions.

The supporting peer-reviewed publication for this package is not yet in press and, therefore, cannot be cited in the 'Description' field of the DESCRIPTION file. One requirement for publication, as requested by the journal editor, is for the package to be available on CRAN. Once the package is available and the manuscript is published, the DESCRIPTION file can be updated with the appropriate citation of the publication.

In this version, I have:
* Added 2 contributors to the Authors@R field of the DESCRIPTION file
* Modified several examples to run short "toy" examples; however, many calculations performed in this package require large data downloads and manipulations of these files. In these cases, examples are provided but not run.

## INITIAL SUBMISSION - 27 OCT 2017
## Test environments
* local OS X install, R 3.3.2
* ubuntu 14.04 (on travis-ci), R 3.4.2
* win-builder, R 3.4.2 and R-devel

## R CMD check results
There were no ERRORs or WARNINGs.

There is 1 NOTE when "checking CRAN incoming feasibility" that indicates this is a new submission.

## Downstream dependencies
There are currently no downstream dependencies for this package.