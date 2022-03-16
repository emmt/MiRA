# Installation of MiRA

There are different possible ways to use MiRA:

* with Yorick interpreter or in your Yorick code with `#include "mira2.i"` (or
  similar);

* from the command line with the `ymira` command.

* via the Docker image `ferreol/mira`.

Below we explain how to install MiRA with
[EasyYorick](https://github.com/emmt/EasyYorick) and how to manually install
MiRA.


## Installation with EasyYorick

Installation of MiRA by [EasyYorick](https://github.com/emmt/EasyYorick) is
fully supported.  This is the easiest and recommended way to install MiRA.

To install EasyYorick, you can do:

```sh
prefix=choose_some_directory
mkdir -p "$prefix/src"
cd "$prefix/src"
git clone https://github.com/emmt/EasyYorick.git ypkg
cd ypkg
./configure --prefix="$prefix"
make install
```

where `$prefix` is some writable directory where to install the software:
executables will go into `$prefix/bin`, sources will go to `$prefix/src`, etc.
For example, you may choose something like `prefix=$HOME/apps`.  Once
`EasyYorick` is installed, the command `ypk` is available in`$prefix/bin/ypkg`.
You may add `$prefix/bin` to your environment variable `PATH`:

```
export PATH="$prefix/bin:$PATH"
```

Assuming you have installed EasyYorick, you just have to execute:

```sh
ypkg install yorick yeti optimpacklegacy yoifits ylib ipy mira
```

which should install MiRA and its dependencies
([Yorick](http://dhmunro.github.io/yorick-doc/),
[Yeti](https://github.com/emmt/Yeti),
[OptimPackLegacy](https://github.com/emmt/OptimPackLegacy),
[YOIFITS](https://github.com/emmt/YOIFITS),
[YLib](https://github.com/emmt/ylib) and
[IPY](https://github.com/emmt/IPY)) if not yet installed.  It is recommended
to also install [YNFFT](https://github.com/emmt/ynfft), a Yorick plugin for the
*Nonequispaced Fast Fourier Transform*.  This is done by adding `ynfft` to the
above list of dependencies, or by executing:

```sh
ypkg install ynfft
```

You may need to install the NFFT library (and header files) before the `ynfft`
Yorick extension; on Ubuntu-like systems, this is done by:

```sh
apt install libnfft3-dev
```

MiRA command is available at `$prefix/bin/ymira`.  You may type:

```sh
ymira --help
```

for a short help and to check that the command `ymira` is available.  To read
the manual page, type:

```sh
man "$prefix/src/mira/doc/ymira.1
```

To upgrade to the lastest MiRA version:

```sh
ypkg upgrade mira
```


## Manual installation

You must first retrieve the code source (see [*Getting MiRA
software*](#getting-mira-software) below).  Then, provided you have installed
required software (see [*Prerequisites*](#prerequisites) below), MiRA software
is usable directly from the sources (see [*Using MiRA without
installing*](#using-mira-without-installing) below) or can be installed for
easier access (see [*Installation of MiRA
software*](#installation-of-mira-software) below).


### Getting MiRA software

MiRA software can be retrieved in two ways: as a Git repository, or as an
archive file.


* To clone MiRA Git repository, do:

  ```sh
  git clone https://github.com/emmt/MiRA.git
  ```

  or (depending what works best for you):

  ```sh
  git clone git@github.com:emmt/MiRA.git
  ```

  In order to pull the last versions of the code, you can do:

  ```sh
  git pull
  ```

* An alternative is to just download and unpack an archive with the code
  source.  First select a version in https://github.com/emmt/MiRA/releases,
  then download it and unpack it.  For instance:

  ```sh
  wget https://github.com/emmt/MiRA/releases/mira-${VERSION}.tar.bz2
  tar jxvf mira-${VERSION}.tar.bz2
  cd mira-${VERSION}
  ```

### Prerequisites

You must have installed the following software:

- [Yorick](http://dhmunro.github.io/yorick-doc/) (version 2.2.04 or superior);
- [Yeti](https://github.com/emmt/Yeti) (version 6.3.2 or superior);
- [OptimPackLegacy](https://github.com/emmt/OptimPackLegacy) (version 1.4.0 or
   superior);
- [YOIFITS](https://github.com/emmt/YOIFITS) for OI-FITS files;
- [YLib](https://github.com/emmt/ylib) for various general purpose utilities;
- [IPY](https://github.com/emmt/IPY) for tools useful to solve inverse
  problems.

For faster operations, it is also recommended to install:

- [YNFFT](https://github.com/emmt/ynfft) (version 1.0.3 or superior), a Yorick
  plugin for the *Nonequispaced Fast Fourier Transform*.


### Installation of MiRA software

To install MiRA, three parameters are required:

* `INCDIR` is the directory where to copy MiRA code files;

* `BINDIR` is the directory where to copy MiRA executable;

* `YORICK` is the path to Yorick executable.

These parameters can be directly specified from the command line:

```sh
make install INCDIR=... BINDIR=... YORICK=...
```

A more simple possibility is to use the `configure` script before installing:

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
  overrides the value set by the configuration script.

* If `BINDIR` is empty then `YORICK` may be empty and the MiRA executable is
  not installed.

* `make distclean` removes the file `install.cfg` where the configuration
  script stores the installation parameters.


### Using MiRA without installing

#### Usage with Yorick interpreter

MiRA Yorick code takes care of locating its own directory and expects that all
Yorick files distributed with MiRA are in this directory.

To use MiRA from Yorick, you just have to do:

```c
include, "INCDIR/mira.i"
```

where `INCDIR` is the location of the MiRA source files.  If you just unpack
the archive, `INCDIR` is the `src` directory of the archive.  You may also copy
all the `src/*.i` files in some other directory at your convenience.


#### Usage from the command line

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

MiRA is available as a Docker image that can be run on any system without any
further installation.
