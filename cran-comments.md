## Comments

Addresses current notes on CRAN:

- adds imports for functions in recommended packages
- registered C routines
- avoid warning re recycling length 1 array
- avoid BibTeX errors by switching to jss.bst

Plus various minor bug fixes.

## Test environments

* Local
 - Windows 8, R 3.5.0
 
* Via R-hub 
 - Mac OS 10.11 El Capitan, R-release (experimental)
 - Ubuntu Linux 16.04 LTS, R-release, GCC
 - Debian Linux, R-devel, GCC
 
## Results

On Ubuntu 16.04, 2 notes:

 - vcdExtra unavailable for testing
 - Author field differs from that derived from Authors@R
     - difference due to orcid, yet these are specified as on https://cran.r-project.org/web/packages/submission_checklist.html