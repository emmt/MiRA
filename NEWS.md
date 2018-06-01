# History of MiRA

## Version 2.0.0 (2018-06-01)
* Plugins can be loaded by MiRA to modify the model of the complex
  visibilities.  See file `mira2_plugin_central_star.i` for a concrete example
  or have a look at [`SPARCO`](https://github.com/kluskaj/mira-sparco) plugin.

## Version 2.0.0b (2018-05-04)
* Can save model complex visibilities (options `--save_visibilities`).
* Add soft-thresholding (skip thresholding if no pixels above soft-threshold
  level).
* Add multi-thread support.
* Use environment variables `MIRA_SRCDIR` and `MIRA_YORICK`.

## Version 2.0.0a (2018-02-12)
* Accounting of bandwidth smearing (with customizable shape and importance
  factor).
* Improved accounting of missing/partial data.
* Multiple possible choices for the objective function used for specific kind
  of (partial) data.
* Speedup building of non-separable linear model (with or without smearing).
* Global setting of the debug and quiet modes.
* Colored messages.

## Version 1.1.2 (2018-05-04)

* Fix accounting of bad data in OI-FITS file.
* Ability to resample images.
* Improve graphics.
* Fix multiple wavelenghts case.

## Version 1.1.1 (2017-01-26)
* Better algorithm to compute the dirty map and the dirty beam.
* Improve command line help, provide a manual page.
* Speed up plots.

## Version 1.1.0 (2017-01-10)
* Provide a command line interface (options in command line have units).
* Provide configuration and installation scripts.
* Now available on GitHub (https://github.com/emmt/MiRA).
* Use Git submodules (see
  https://chrisjean.com/git-submodules-adding-using-removing-and-updating/ for
  a tutorial) for YLib, YOIFITS, IPY to share files with other repositories.

## Old releases
2015-05-04:
* Version 1.0.1 released.
* Add mask in total variation regularization.

2010-12-xx:
* Now use `fits.i` distributed with (CVS) Yorick.

2010-07-05:
* MiRA 0.9.10 released.
* Changed `T_CHAR`, `T_SHORT`, etc. to `Y_CHAR`, `Y_SHORT`, etc. to work with
  Yorick >= 2.2 and Yeti >= 6.3.1

2009-05-14:
* MiRA 0.9.9 released.

2009-04-23:
* MIRA 0.9.8 released.
* Changes to use OptimPack-1.3 and to display Fdata and Fprior.
* Complex visibilities can be fitted with (or without) Goodman
  approximation.

2008-12-12:
* MIRA 0.9.7 released.
* Fixed oifits.i for empty ARRNAME or INSNAME.

2008-12-09:
* MIRA 0.9.6 released.
* Fixed a bug in `oifits.i` with `oifits_new_*` when master is provided
  (thanks to Thibaut Paumard).
* Fixed a bug in `mira_new_fft_xform` with spatial frequencies exactly
  equal to zero (thanks to Stéphanie Renard).

2008-10-03:
* MIRA 0.9.5 released.
* Fix bug when no target is specified in `mira_add_oidata`.
* In `mira_add_oidata/mira_new`, keyword `cleanup_bad_data` can be set
  to 0, 1, 2 to achieve different levels of filtering of invalid data.
* Various fixes in `oifits.i` for creating/saving OI-FITS files (thanks
  to Sylvestre Lacour).

2008-09-26:
* MIRA 0.9.4 released.
* Added the possibility to select the target in `mira_new`.
* Hack for FLAG column to deal with AMBER data.

2008-09-23:
* MIRA 0.9.3 released.
* Fixed a bug in `__mira_build_coordinate_list__` which prevent to use
  data files with central frequency (0,0) measured.
* New function mira_dirac.

2008-09-07:
* Massive rewrite of `oifits.i` to optimize the code, make it easier to read
  and let the users create OI-FITS data on the fly and save it to a file.

2008-09-04:
* MIRA 0.9.2 released.  This is the version used for the MIRA demonstration
  (see `mira-demo.i`) at SPIE 2008 Conference in Marseille (France).
* Some change in plots.

2008-07-12:
* MIRA 0.9.1 released.  This is the version used for the VLTI 2008 Summer
  School at Keszthely (Hungary).
* Missing `fft_utils.i` is now part of the distribution.

2008-06-06:
* MIRA 0.9 released.

2008-05-02:
* MIRA 0.8 released.

2008-01-31:
* MIRA 0.7 released.
* Orientation of u-v and image coordinates fixed to match astronmical
  conventions.
* Use OptimPack1.
* Noisy works again.
* New Monte-Carlo optimizer.
* Fix a bug in `mira_add_oidata` with loading of complex visibilities.


2007:
* Penalty with respect to data is now computed for cartesian representation
  of a complex (not a polar).
* Huge rewrite of code: all data of given type are collected in a single
  table (more compact representation and faster computation/plotting when a
  lot of interferometric data-blocks are fitted).
* Wavelength selection is more flexible (you can specify the effective
  wavelength and bandwidth to consider).
* Very fast approximation of the Fourier transform by means of FFT (or FFTW)
  and spectral interpolation.
