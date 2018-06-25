[![Build Status](https://travis-ci.com/fuchslab/mpb2.svg?token=oRmho23ZpPzFsxXkFqih&branch=master)](https://travis-ci.com/fuchslab/mpb2)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)


mpb2
================

mpb2 is an R package for the poisson-beta distribution. At its current stage, the package has implementations for the density, distribution function, quantile function and random number generation. The package depends on the MPFR/GMP libraries.



### Installation
--------------

#### MPFR on Linux
On Ubuntu, the mpfr libraries can be installed with the following command:

```bash
sudo apt-get install -qy libmpfr-dev
```

For RPM based distributions, the corresponding package may be found at [RPMfind](https://fr2.rpmfind.net/linux/rpm2html/search.php?query=mpfr-devel).

#### MPFR on OS X
MPFR libraries can be installed on OS X with the [Homebrew](https://brew.sh/) package manager.

```bash
brew install mpfr
```

Installing Homebrew is easy. The following command does it.

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

#### MPFR on Windows
Installation on Windows requires building MPFR with the recommended toolchain for R.

Pre-compiled libraries have been put together by Prof. Brian Ripley and are available on his [webpage](http://www.stats.ox.ac.uk/pub/Rtools/). The files are under [goodies/multilib](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/) directory. This includes a 'local' tree, named as [local323.zip](http://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/local323.zip) at the time this manual was prepared, that contains the headers and static libraries for MPFR amongst others.

The following instructions are vaild for Windows 8 and above.

1. The local323.zip has to be extracted at a convenient place. The C: drive is the location for the screenshots used in this manual.
2. The location of the library needs to be added to the environment variables of the system. There are two ways to get to the list of environment variables:

    * One of the ways is to select Computer -> System properties in the Explorer. In the next window, select Advanced system settings.
    <img src="images/windows-1.png" width=400/>
    <img src="images/windows-2.png" width=400/>
    
    * Another option is to use the Windows search bar as shown in the picture below:  
    <img src="images/windows-3.png" width=300/>
3. The location needs to be added against the variable LIB_MPFR.  
    <img src="images/windows-4.png" width=300/>
    <img src="images/windows-5.png" width=300/>
    <img src="images/windows-6.png" width=300/>
    <img src="images/windows-7.png" width=300/>
4. Rstudio (or the current R session) has to be restarted so that these environment variables can be read.


### Running tests for the package
---------------

The package comes with a comprehensive test suite. In order to run the tests, the tarball must be extracted. Once that is done, the following command will run all tests:

```r
devtools::test("/path/to/mpb2/")
```
