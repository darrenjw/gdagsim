GDAGsim
=======

### C Library for conditional simulation of Gaussian DAG models

#### Please note that this repo is an archive of the source of a very old software library that is no longer supported or maintained. It is largely of historical interest.

[gdagsim-03.tgz](gdagsim-03.tgz) is the current release of `GDAGsim`.
This is a C library for analysis and block-sampling of Gaussian DAG
models. Full documentation and installation instructions are included
with the package - see the README in the top-level directory. To use the
library you will need to be reasonably familiar with programming in C,
and it helps (but is not vital) if you\'ve used the [GNU Scientific
Library](http://sourceware.cygnus.com/gsl/) (GSL) before. Note that you
must have the GSL (\>=1.0) installed and working correctly before
attempting to install GDAGsim. Since there are (currently) no sparse
matrix algorithms in the GSL, this version of GDAGsim also depends on
the [Meschach](ftp://ftpmaths.anu.edu.au/pub/meschach/meschach.html)
matrix library. However, no familiarity with Meschach is required.

The main change from 0.2 is a change of syntax for the stochastic
simulation functions. This syntax change is annoying, but the decision
to introduce it has not been taken lightly! The library is now
thread-safe. This will ease the development of multi-threaded and
parallel codes which use the library, and also has some other benefits.

If you use Linux, then pre-compiled binary packages are available for
both the GSL and Meschach - just install the library packages and the
corresponding \"-dev\" packages before installing GDAGsim, and
everything should work fine. I hope to put together binary packages for
GDAGsim for use with Debian GNU/Linux sometime in the next few months,
but for now, building from source should be very straightforward.

The latest [documentation](gdag.pdf) is available in PDF format for
on-line browsing. The theory behind GDAGsim is explained in the
following paper:

Wilkinson, DJ & Yeung, SKH (2004).
:   [A sparse matrix approach to Bayesian computation in large linear
    models.](http://dx.doi.org/10.1016/S0167-9473(02)00252-9)
    *Computational Statistics and Data Analysis*, **44**(3):493-516.
