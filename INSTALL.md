# Installation of MiRA

There are 2 possible ways to use MiRA:

* with Yorick interpreter or in your Yorick code with `#include "mira.i"` (or
  similar);

* from the command line with the `ymira` command.

* via the Docker images `ferreol/mira`

Provided you have installed required software (see
[*Prerequisites*](#prerequisites) below), MiRA software is usable right after
unpacking the archive or cloning the repository and its sub-modules (see
[*Using MiRA without installing*](#using-mira-without-installing) below) or can
be installed for easier access (see
[*Installation of MiRA software*](#installation-of-mira-software) below).


## Getting MiRA software

MiRA software can be retrieved in two ways: as a Git repository, or as an
archive file.


### Cloning Git repository

To clone MiRA Git repository, do:

```sh
git clone https://github.com/emmt/MiRA.git
```

or (depending what work best for you):

```sh
git clone git@github.com:emmt/MiRA.git
```

In order to pull the last versions of the sources, you can do:

```sh
git pull
```

### Getting a source archive

Select a version in https://github.com/emmt/MiRA/releases, download it and
unpack it.  For instance:

```sh
wget https://github.com/emmt/MiRA/releases/mira-${VERSION}.tar.bz2
tar jxvf mira-${VERSION}.tar.bz2
cd mira-${VERSION}
```

## Prerequisites

You must have installed the following software:

- [Yorick](http://dhmunro.github.io/yorick-doc/) (version 2.2.04 or superior);
- [Yeti](https://github.com/emmt/Yeti) (version 6.3.2 or superior);
- [OptimPackLegacy](https://github.com/emmt/OptimPackLegacy) (version 1.4.0 or
   superior);
- [YOIFITS](https://github.com/emmt/YOIFITS) for OI-FITS files;
- [YLIB](https://github.com/emmt/ylib) for various general purpose utilities;
- [IPY](https://github.com/emmt/IPY) for tools useful to solve inverse
  problems.

For faster operations, it is also recommended to install:

- [YNFFT](https://github.com/emmt/ynfft) (version 1.0.3 or superior), a Yorick
  plugin for the *Nonequispaced Fast Fourier Transform*.


## Installation of MiRA software

To install MiRA, three parameters are required:

* `INCDIR` is the directory where to copy MiRA code files;

* `BINDIR` is the directory where to copy MiRA executable;

* `YORICK` is the path to Yorick executable.

These parameters can be directly specified from the command line:

```sh
make install INCDIR=... BINDIR=... YORICK=...
```

Another possibility is to use the `configure` script before installing:

```sh
./configure ...
make install ...
```

The configuration script is able to automatically find Yorick executable and,
by default, set `INCDIR` to be `Y_SITE/i` where `Y_SITE` is the platform
independent installation directory of Yorick.  With these defaults, it is
sufficient to do:

```c
#include "mira.i"
```

to use MiRA in your Yorick code.  To have a description of available options:

```sh
./configure --help
```

*Remarks:*

* Setting installation parameters in the `make install ...` command line
  override the value set by the configuration script.

* If `BINDIR` is empty then `YORICK` may be empty and the MiRA executable is
  not installed.

* `make distclean` removes the file `install.cfg` where the configuration
  script stores the installation parameters.


## Using MiRA without installing

### Usage with Yorick interpreter

MiRA Yorick code takes care of locating its own directory and expects that all
Yorick files distributed with MiRA are in this directory.

To use MiRA from Yorick, you just have to do:

```c
include, "INCDIR/mira.i"
```

where `INCDIR` is the location of the MiRA source files.  If you just unpack
the archive, `INCDIR` is the `src` directory of the archive.  You may also copy
all the `src/*.i` files in some other directory at your convenience.


### Usage from the command line

MiRA can can be run from the command line with the `mira` script.  This script
can be installed anywhere but must know the directory where MiRA source files
are copied and the path to Yorick executable.

If Yorick executable is in your shell command path, you can directly run the
`ymira` command from the `bin` directory of the distribution.

To install the `ymira` command in a given directory, say `BINDIR`, you just
have to copy the script `bin/ymira` in `BINDIR` and edit it so that variables
`INCDIR` and `YORICK` are correctly set, the former to the directory where MiRA
sources have been copied, the latter to the path of the Yorick interpreter
executable.  For instance:

```sh
cp bin/ymira "$BINDIR"
edit "$BINDIR/ymira"
chmod 755 "$BINDIR/ymira"
```


###  MiRA  via Docker

MiRA is available as a Docker image that can be run on any system without any further installation.
