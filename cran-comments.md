## Test environments
* Linux
  - apt packages: build-essential, libmpfr-dev, libblas-dev, liblapack-dev
* Windows
  - [precompiled libraries](http://www.stats.ox.ac.uk/pub/Rtools/)
* OS X
  - brew packages: mpfr
* Travis-CI
  - Linux: libmpfr-dev
  - OS X: mpfr

## R CMD check results
There were no ERRORs, WARNING or NOTEs.
Status: OK

## Update
This is an update fixing warnings generated from using sprintf.
