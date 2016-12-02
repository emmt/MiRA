# Installation of MiRA

There are 2 possible ways to use MiRA:

* with Yorick interpreter or in your Yorick code with `#include "mira.i"` (or
  similar);

* from the command line with the `mira` command.

Provided you have installed required software (see
[*Prerequisites*](#prerequisites) below), MiRA software is usable right after
unpacking the archive (see
[*Using MiRA without installing*](#using-mira-without-installing) below) or can
be installed for easier access (see
[*Installation of MiRA software*](#installation-of-mira-software) below).


## Prerequisites

You must have [Yorick](http://dhmunro.github.io/yorick-doc/) (at least
version 2.1), [Yeti](https://github.com/emmt/Yeti) (at least version 6.2.0) and
[OptimPack](https://cral.univ-lyon1.fr/labo/perso/eric.thiebaut/?Software/OptimPack)
(version 1.2) installed.  For faster operations, you may want to use
[YNFFT](https://github.com/emmt/ynfft), a Yorick plugin for the *Nonequispaced
Fast Fourier Transform*.


## Installation of MiRA software

To install MiRA, three parameters are required:

* `INCDIR` is the directory where to copy MiRA code files;

* `BINDIR` is the directory where to copy MiRA executable;

* `YORICK` is the path to Yorick executable.

These parameters can be directly specified from the command line:

    make install INCDIR=... BINDIR=... YORICK=...

Another possibility is to use the `configure` script before installing:

    ./configure ...
    make install ...

The configuration script is able to automatically find Yorick executable and,
by default, set `INCDIR` to be `Y_SITE/i` where `Y_SITE` is the platform
independent installation directory of Yorick.  With these defaults, it is
sufficient to do:

    #include "mira.i"

to use MiRA in your Yorick code.  To have a description of available options:

    ./configure --help


*Remarks:*

* Setting installation parameters in the `make install ...` command line
  override the value set by the configuration script.

* If `BINDIR` is empty then `YORICK` may be empty and the MiRA executable is
  not installed.

* `make clean` removes the file `install.cfg` where the configuration
  script stores the installation parameters.


## Using MiRA without installing

### Usage with Yorick interpreter

MiRA Yorick code takes care of locating its own directory and expects that all
Yorick files distributed with MiRA are in this directory.

To use MiRA from Yorick, you just have to do:

    include, "INCDIR/mira.i"

where `INCDIR` is the location of the MiRA source files.  If you just unpack
the archive, `INCDIR` is the `src` directory of the archive.  You may also copy
all the `src/*.i` files in some other directory at your convenience.


### Usage from the command line

MiRA can can be run from the command line with the `mira` script.  This script
can be installed anywhere but must know the directory where MiRA source files
are copied and the path to Yorick executable.

If Yorick executable is in your shell command path, you can directly run the
`mira` command from the `bin` directory of the distribution.

To install the `mira` command in a given directory, say `BINDIR`, you just have
to copy the script `src/mira` in `BINDIR` and edit it so that variables
`INCDIR` and `YORICK` are correctly set, the former to the directory where MiRA
sources have been copied, the latter to the path of the Yorick interpreter
executable.  For instance:

    cp src/mira "$BINDIR"
    edit "$BINDIR/mira"
    chmod 755 "$BINDIR/mira"

