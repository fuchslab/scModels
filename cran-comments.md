## Test environments
* Linux
  - R 3.2.3
  - apt packages: build-essential, libmpfr-dev, libblas-dev, liblapack-dev
* Windows
  - R 3.4.1, R 3.5.2
  - [precompiled libraries](http://www.stats.ox.ac.uk/pub/Rtools/)
* OS X
  - R 3.5.3
  - brew packages: mpfr
* Travis-CI
  - R 3.5.2 (Linux)
  - R 3.5.3 (OS X)

## R CMD check results
There were no ERRORs, WARNING or NOTEs.
Status: OK

## Update
This is an update fixing warnings generated for deprecated functions in mpfr. 
