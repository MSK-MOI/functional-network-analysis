## Test environments
* local install Arch Linux 9.3.0-1, R 4.0.0
* win-builder


## R CMD check results
Checked with R-release 4.0.0 and with R-devel 2020-05-04 r78358.

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘James Mathews <mathewj2@mskcc.org>’

New submission

* checking top-level files ... NOTE
Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.
Non-standard file/directory found at top level:
  ‘cran-comments.md’

The second NOTE only appeared when checking with R-devel. R-devel indicated that 'pandoc' cannot be installed for R-devel, so I am hoping this part is not a problem.
As for the "new submission" indication, this is correct.
Also, all of the documentation seems to say that 'cran-comments.md' is a standard file (it is this file, in fact), and in several examples it seems to appear 'at the top-level'. So I left it where it is.


## Downstream dependencies
There are no downstream dependencies.
