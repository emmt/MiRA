/*
 * mira.i -
 *
 * Implement MiRA (Multi-aperture Image Reconstruction Algorithm) in
 * Yeti/Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2001-2010, Éric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 * This file is part of MiRA: a Multi-aperture Image Reconstruction
 * Algorithm.
 *
 * MiRA is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License version 2 as published by the Free
 * Software Foundation.
 *
 * MiRA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 *-----------------------------------------------------------------------------
 */

MIRA_VERSION = "1.0.0";

local MIRA_VERSION;
local mira;
/* DOCUMENT MiRA: a Multi-aperture Image Reconstruction Algorithm.
  
     MiRA (Multi-aperture Image Reconstruction Algorithm) is a software tool
     for image reconstruction from interferometric data.
  
     Global variable MIRA_VERSION is a stting with the current version of
     MiRA.
  
   SEE ALSO: mira_new, mira_config,
 */

if (is_void(MIRA_HOME) && strcase(0, get_env("USER")) == "eric") {
  write, format="*** FIXME: %s\n",
    ["(auto)recentering to solve translation degeneracy",
     "use OIFITS FLAGS for extracting relevant data",
     "peek computation of diagonal of Hessian in old mira.i"];
}

/*---------------------------------------------------------------------------*/
/* INITIALIZATION OF MIRA */

/* First of all MiRA requires hash-tables (and some other functions)
   provided by Yeti: */
if (is_func(h_new) != 2) {
  include, "yeti.i", 3;
  if (is_func(h_new) != 2) {
    error, "Yeti is mandatory to run MiRA.";
  }
}

local MIRA_HOME;
/* DOCUMENT MIRA_HOME
 *   Global variable used to store the full path to the directory where MiRA
 *   software suite is installed.
 *
 * SEE ALSO: setup_package.
 */
MIRA_HOME = setup_package();

func mira_include(sym, src)
/* DOCUMENT mira_include, sym, src;
 *   Include source file SRC if symbol SYM is not a function.  This is
 *   a shortcut to:
 *      if (! is_func(SYM)) include, SRC, 1;
 *
 * SEE ALSO: include, is_func, require.
 */
{
  code = is_func(sym); 
  if (code != 1 && code != 2) {
    include, src, 1;
  }
}

/* Loading of utility functions: */
mira_include, glob, MIRA_HOME + "utils.i";

/* MiRA requires OI-FITS support: */
mira_include, fits_open, "fits.i";
mira_include, oifits_load, MIRA_HOME + "oifits.i";

/* MiRA requires linear operator class: */
mira_include, linop_new, MIRA_HOME + "linop.i";

/* MiRA requires regularization operators: */
mira_include, rgl_new, MIRA_HOME + "rgl.i";

/* MiRA requires FFT_UTILS: */
mira_include, fft_indgen, MIRA_HOME + "fft_utils.i";

/* MiRA requires OptimPack1: */
mira_include, op_vmlmb_next, "OptimPack1.i";
mira_include, fmin, MIRA_HOME + "fmin.i";

/* MiRA requires additional plot functions.  Unfortunately many people have
   already included their own (old) "plot.i" file so we force loading the
   version delivered with MiRA: */
include, MIRA_HOME + "plot.i", 1;

/* Load some files from the standard Yorick installation. */
mira_include, bessj0, Y_SITE + "i/bessel.i";
mira_include, random_n, Y_SITE + "i/random.i";

/* Some constants. */
local MIRA_PI, MIRA_MICRON;
local MIRA_DEGREE, MIRA_ARCSECOND, MIRA_MILLIARCSECOND;
/* DOCUMENT MIRA_PI             = 3.1415.....
 *     -or- MIRA_MICRON         = micron to meter conversion factor
 *     -or- MIRA_DEGREE         = degree to radian conversion factor
 *     -or- MIRA_ARCSECOND      = arcsecond to radian conversion factor
 *     -or- MIRA_MILLIARCSECOND = milliarcsecond to radian conversion factor
 *
 * SEE ALSO: mira.
 */
MIRA_PI = 3.141592653589793238462643383279503;
MIRA_TWO_PI = 2.0*MIRA_PI;
MIRA_DEGREE = MIRA_PI/180.0;
MIRA_ARCSECOND = MIRA_DEGREE/3600.0;
MIRA_MILLIARCSECOND = 1e-3*MIRA_ARCSECOND;
MIRA_MICRON = 1e-6;

/* Various global options. */
MIRA_SPARSE = 1n; /* use Yeti sparse matrix */
MIRA_DEBUG = 0n; /* print out some debug messages */
MIRA_FLAGS = 0n;

MIRA_USE_NORMALIZE      = 1;
MIRA_USE_PHASE          = 2;
MIRA_USE_AMPLITUDE      = 4;
MIRA_USE_POWER_SPECTRUM = 8;

MIRA_COMPLEX_VISIBILITY           = 1; /* real and imaginary parts of complex visibility */
MIRA_COMPLEX_VISIBILITY_POLAR     = 2; /* amplitude and phase of complex visibility */
MIRA_COMPLEX_VISIBILITY_AMPLITUDE = 3; /* amplitude of complex visibilities */
MIRA_COMPLEX_VISIBILITY_PHASE     = 4; /* phase of complex visibilities */
MIRA_POWER_SPECTRUM               = 5; /* power-spectrum */
MIRA_BISPECTRUM                   = 6; /* real and imaginary parts of bispectrum */
MIRA_BISPECTRUM_POLAR             = 7; /* amplitude and phase of bispectrum */
MIRA_BISPECTRUM_AMPLITUDE         = 8; /* amplitude of bispectrum */
MIRA_BISPECTRUM_PHASE             = 9; /* phase of bispectrum (phase closure) */

/*---------------------------------------------------------------------------*/

/*
 * Functions
 * ~~~~~~~~~
 *   mira_new - create a new MiRA opaque object
 *   mira_add_oidata - append interferometric data from OI-FITS file or handle
 *   mira_config - setup options for image reconstruction
 *
 *   mira_data_penalty - compute misfit value and gradient
 *   mira_dirty_beam - estimate the "dirty beam" of actual u-v coverage
 *   mira_update
 *   mira_recenter
 *   mira_solve - driver for image reconstruction
 *
 *   mira_get_fov - get width of the model image field of view (in radians)
 *   mira_get_dim - get number of pixels per side of model image
 *   mira_get_pixelsize - get the pixel size (in radians) for the model image
 *   mira_get_w - returns wavelength(s)
 *   mira_get_x - returns sky X-coordinates
 *   mira_get_y - returns sky Y-coordinates
 *   mira_get_ndata - get number of measurements used in last fit
 *
 *   mira_dirac - make an image of a point-like object
 *
 *   mira_get_one_integer - get an integer scalar
 *   mira_get_one_real - get a real scalar
 *
 *   mira_new_exact_xform - build exact Fourier transform operator
 *   mira_new_fft_xform - build FFT-based Fourier transform operator
 *
 *   mira_classify - classify values
 *   mira_digitize - digitize values
 *
 *
 * Obsolete or unused?
 * ~~~~~~~~~~~~~~~~~~~
 *   mira_bicubic_ft
 *   mira_bilinear_ft
 *   mira_cast_complex_as_real
 *   mira_cast_real_as_complex
 *   mira_color_bar
 *   mira_fit_profile
 *   mira_gauss_ft
 *   mira_glob
 *   mira_include
 *   mira_plot_baselines
 *   mira_plot_frequencies
 *   mira_plot_image
 *   mira_plot
 *   mira_regul0mask
 *   mira_regul1
 *   mira_regul2
 *   mira_regul3
 *   mira_regul4
 *   mira_regul_fov_l2
 *   mira_regul_mem2
 *   mira_regul_mem3
 *   mira_regul_mem
 *   mira_regul_roughness
 *   mira_relative_absolute_difference
 *   mira_rescale
 *   mira_stdev_to_weight
 *   mira_trash
 *   mira_udisk_ft
 *   mira_weight_to_stdev
 *
 *
 * Master Opaque Object
 * ~~~~~~~~~~~~~~~~~~~~
 *   MASTER.vis - complex visibility data (see below for layout)
 *   MASTER.vis2 - powerspectrum data (see below for layout)
 *   MASTER.vis3 - bispectrum data (see below for layout)
 *   MASTER.monochromatic - monochromatic or gray case?
 *   MASTER.monochromatic_option - monochromatic option as set by user
 *   MASTER.u - list of measured spatial frequencies
 *   MASTER.v - list of measured spatial frequencies
 *   MASTER.w - list of measured wavelenghts
 *   MASTER.update_pending - mira_update must be called
 *
 *
 * Interferometric Coordinates
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   In principle, there are up to four different coordinates: u, v,
 *   lambda, and t.  In MiRA description, the object brightness
 *   distribution is assumed to be "static" hence the time t is not
 *   considered.  Besides, in the polychromatic case, spatial (u,v)
 *   and spectral (lambda) coordinates are assumed to be separable
 *   (because the "image" is a 3D array and because the Fourier
 *   transform is approximated by FFT's in that case).
 *
 *     - specify orientation for (u,v)
 *
 *     - specify resolution for (u,v,lambda,t)
 *
 *     - coordinates for all measurements are collected and a minimal
 *       set of "unique" coordinates is build (restricted to the half
 *       u-v plane), all measurement coordinate are related to the set
 *       of unique coordinates by an index and a sign for (u,v);
 *
 *
 * Interferometric Data
 * ~~~~~~~~~~~~~~~~~~~~
 *   Complex visibility:
 *     THIS.type - MIRA_COMPLEX_VISIBILITY
 *     THIS.t - exposure time
 *     THIS.u - spatial frequency
 *     THIS.v - spatial frequency
 *     THIS.w - effective wavelength [m]
 *     THIS.re - real part of (calibrated) complex visibility
 *     THIS.im - imaginary part of (calibrated) complex visibility
 *     THIS.wrr - weight for real part
 *     THIS.wri - weight for cross term real-imaginary part
 *     THIS.wii - weight for imaginary part
 *     THIS.idx - index of coordinates in global list
 *     THIS.sgn - sign of (u,v) coordinates w.r.t. to global list
 *
 *   Note that other members are allowed but not considered by the
 *   image reconstruction.
 *
 *   Bispectrum:
 *     THIS.type - MIRA_BISPECTRUM
 *     THIS.t - exposure time
 *     THIS.u1 - spatial frequency
 *     THIS.v1 - spatial frequency
 *     THIS.u2 - spatial frequency
 *     THIS.v2 - spatial frequency
 *     THIS.w - effective wavelength [m]
 *     THIS.re - real part of bispectrum
 *     THIS.im - imaginary part of bispectrum
 *     THIS.wrr - weight for real part
 *     THIS.wri - weight for cross term real-imaginary parts
 *     THIS.wii - weight for imaginary part
 *     THIS.idx1 - index of coordinates in global list for (u1,v1)
 *     THIS.sgn1 - sign of (u1,v1) coordinates w.r.t. to global list
 *     THIS.idx2 - index of coordinates in global list for (u2,v2)
 *     THIS.sgn2 - sign of (u2,v2) coordinates w.r.t. to global list
 *     THIS.idx3 - index of coordinates in global list for (u3,v3)
 *     THIS.sgn3 - sign of (u3,v3) coordinates w.r.t. to global list
 *
 *     where by definition (u1,v1) + (u2,v2) + (u3,v3) = (0,0),
 *     hence: u3 = - u1 - u2 and v3 = - v1 - v2.
 *
 *   Powerspectrum (squared visibility):
 *     THIS.type - MIRA_POWERSPECTRUM
 *     THIS.t - exposure time
 *     THIS.u - spatial frequency
 *     THIS.v - spatial frequency
 *     THIS.w - effective wavelength [m]
 *     THIS.vis2data - powerspectrum data
 *     THIS.vis2err - standard deviation of powerspectrum data
 *     THIS.vis2wght - weight of powerspectrum data
 *     THIS.idx - index of coordinates in global list
 *     THIS.sgn - sign of (u,v) coordinates w.r.t. to global list
 *
 *   The index (idx) and sign (sgn) members work as follows:
 *
 *     DB.u ~ MASTER.u(DB.idx)*DB.sgn
 *     DB.v ~ MASTER.v(DB.idx)*DB.sgn
 *     DB.w ~ MASTER.w(DB.idx)            // not in monochromatic mode
 *
 *   where MASTER is the parent of datablock DB
 *
 */

func mira_new(.., eff_wave=, eff_band=, wave_tol=,
              quiet=, base_tol=, monochromatic=,
              noise_method=, noise_level=,
              cleanup_bad_data=, target=, goodman=)
/* DOCUMENT mira_new(filename1, filename2, ...)
  
     Return a new instance of MiRA data opaque structure filled with data read
     from files FILENAME1, FILENAME2, ...  Input data files must follow
     OI-FITS standard.  The following keywords are allowed:
  
       EFF_WAVE = Effective wavelength (units: meters), data with wavelengths
           in the range EFF_WAVE +/- 0.5*EFF_BAND will be selected for image
           reconstruction.  If EFF_WAVE is not specified, the average
           wavelength from the first block of data is used.
  
       EFF_BAND = Effective spectral bandwidth (units: meters), default value
           is 1e-7 (0.1 micron).
  
       WAVE_TOL = Tolerance for wavelength grouping (units: meters).  Default
           value is 1e-10 (1 Ångström).  This tolerance is used to decide
           whether different wavelengths correspond to the same one.
  
       BASE_TOL = Tolerance for baseline grouping (units: meters).  Default
           value is 1e-3 (1 millimeter).  This tolerance is used to decide
           whether different positions correspond to the same baseline.
  
       MONOCHROMATIC = True if a monochromatic (gray) model of the object
           brightness distribution is to be reconstructed.
  
       QUIET = Turn off informational messages?  Set 2nd bit to one to also
           turn off the printing of the filename which is loaded and summary
           of u-v coverage and maximum pixel size.
  
       TARGET = Target name as a glob-style pattern (see strglob), the
           comparison is case-insensitive and leading/trailing spaces are
           ignored.  Using TARGET="*" will select all available targets (only
           recommended if you are sure that all the data really come from the
           same object).
  
       CLEANUP_BAD_DATA = Delete invalid data?  If true, data with invalid
           error bars (ERR <= 0) get removed.  If CLEANUP_BAD_DATA > 1, then
           data with out of range amplitudes (AMP < 0 or AMP > 1) get also
           removed.
  
       GOODMAN = True to force Goodman approximation for the penalty with
           respect to complex visibility or bispectrum data (default is
           false).
  
     It is possible to add noise to the data by specifying the keyword
     NOISE_METHOD and, possibly, the keyword NOISE_LEVEL:
  
       NOISE_METHOD = 1 or "generate" to add noise to noiseless data; the
           standard deviation of the noise is taken from the contents of
           MASTER; in this case, the NOISE_LEVEL argument must be nil or
           omitted.
  
       NOISE_METHOD = 2 or "snr" to add noise to noiseless data; the standard
           deviation of noise is computed to achieve a signal-to-noise ratio
           equal to the value of NOISE_LEVEL.
  
       NOISE_METHOD = 3 or "amplify" to add noise to noisy data so that the
           standard deviation of total noise (existing one plus added one) is
           multiplied by the value of NOISE_LEVEL which must be greater or
           equal one.  The standard deviation of the noise prior to the
           amplification is taken from the contents of MASTER.
  
  
   SEE ALSO: mira_add_oidata, mira_config, mira_solve.
 */
{
  if (is_void(quiet)) quiet = 0;

  /* Get spectral bandwidth parameters (in meters). */
  if (! is_void(eff_wave) &&
      (mira_get_one_real(eff_wave) || eff_wave <= 0.0)) {
    error, "bad value for effective wavelength (EFF_WAVE in meters)";
  }
  if (mira_get_one_real(eff_band, 1e-7) ||
      eff_band <= 0.0 || eff_band >= 1e2) {
    error, "bad value for effective spectral bandwidth (EFF_BAND in meters)";
  }

  /* Absolute tolerance for wavelengths (in meters). */
  if (mira_get_one_real(wave_tol, 1e-10) ||
      wave_tol < 0.0 || wave_tol >= 1e-2) {
    error, "bad value for wavelength tolerance (WAVE_TOL in meters)";
  }

  /* Absolute tolerance for baselines (in meters). */
  if (mira_get_one_real(base_tol, 1e-3) ||
      base_tol < 0.0 || base_tol >= 1.0) {
    error, "bad value for absolute baseline tolerance (BASE_TOL in meters)";
  }

  if (! is_void(target)) {
    if (is_string(target) && is_scalar(target)) {
      target = strtrim(target);
    } else {
      error, "invalid value for keyword TARGET";
    }
  }

  master = h_new(eff_wave = eff_wave,
                 eff_band = eff_band,
                 wave_tol = wave_tol,
                 base_tol = base_tol,
                 monochromatic_option = monochromatic,
                 flags = 0,
                 flux_weight = 0.0,
                 flux_mean = 1.0,
                 update_pending = 1n,
                 target = target,
                 plot_model_options = h_new(color="green",
                                            width=1,
                                            type="solid",
                                            symbol=6,
                                            size=0.3,
                                            fill=0n,
                                            ticks=0n),
                 plot_data_options = h_new(color="red",
                                           width=1,
                                           type="solid",
                                           symbol=2,
                                           size=0.3,
                                           fill=0n,
                                           ticks=1n,
                                           residuals=1n));

  /* Load OI-FITS data file(s). */
  local arg;
  while (more_args()) {
    eq_nocopy, arg, next_arg();
    if (is_string(arg)) {
      n = numberof(arg);
      for (i = 1; i <= n; ++i) {
        if (! (quiet & 2)) write, format="Loading file \"%s\"...\n", arg(i);
        mira_add_oidata, master, arg(i), quiet=quiet,
          noise_method=noise_method, noise_level=noise_level,
          cleanup_bad_data=cleanup_bad_data, goodman=goodman;
      }
    } else {
      mira_add_oidata, master, arg, quiet=quiet,
          noise_method=noise_method, noise_level=noise_level,
          cleanup_bad_data=cleanup_bad_data, goodman=goodman;
    }
  }

  /* Build list of coordinates. */
  local u, v;
  _mira_build_coordinate_list, master;
  eq_nocopy, u, master.u;
  eq_nocopy, v, master.v;
  if (is_void(u)) {
      write, format="WARNING %s\n", "no data";
  } else {
    freq_max = max(max(u), -min(u), max(v), -min(v));
    freq_len = abs(u, v);
    if (! (quiet & 2)) {
      shannon = 0.5/freq_max;
      write, format="CUTOFF FREQUENCY = %g per radians\n", freq_max;
      write, format="==> MAX PIXEL SIZE = %g radians\n", shannon;
      write, format="                   = %g milliarcseconds\n",
        shannon/MIRA_MILLIARCSECOND;
    }
  }
  return h_set(master, freq_max=freq_max, freq_len=freq_len);
}

func mira_add_oidata(this, .., quiet=, noise_method=, noise_level=,
                     cleanup_bad_data=, goodman=)
/* DOCUMENT mira_add_oidata, this, data, ...;
  
     Append interferometric OI-FITS data to MiRA handle THIS (as created by
     mira_new, which see).  DATA and subsequent arguments are either OI-FITS
     file names or opaque OI-FITS handles as returned by oifits_load (which
     see).
  
     If keyword QUIET is true, the operation is performed silently.
  
     Keywords NOISE_METHOD and NOISE_LEVEL can be used to add some noise to
     the data.
  
     If keyword CLEANUP_BAD_DATA is true, then data with invalid error bars
     (ERR <= 0) get removed.  If CLEANUP_BAD_DATA > 1, then data with out of
     range amplitudes (AMP < 0 or AMP > 1) get also removed.
  
     Keyword GOODMAN can be used to force Goodman approximation for the
     penalty with respect to complex visibility or bispectrum data (default is
     false).
  
   SEE ALSO: mira_new, mira_noise, oifits_load.
 */
{
  /* For the following tests, CLEANUP_BAD_DATA must be a numerical scalar. */
  if (is_void(cleanup_bad_data)) {
    cleanup_bad_data = 0;
  } else if (! is_scalar(cleanup_bad_data) || ! is_numerical(cleanup_bad_data)
             || is_complex(cleanup_bad_data)) {
    error, "CLEANUP_BAD_DATA must be nil or an integer/real scalar";
  }

  /* Process all inputs. */
  local data;
  while (more_args()) {
    eq_nocopy, data, next_arg();

    /* Load OIFITS file. */
    if (is_string(data)) {
      data = oifits_load(data, quiet=quiet, errmode=1n);
    }
    if (! is_void(noise_method)) {
      oifits_add_noise, data, noise_method, noise_level;
    }

    /* Search the numerical identifier of the target.  FIXME: targets could be
       identified by their coordinates RA/DEC? */
    flag = 0n;
    for (db = oifits_first(data); db; db = oifits_next(data, db)) {
      if (oifits_get_type(db) == OIFITS_TYPE_TARGET) {
        flag = 1n;
        break;
      }
    }
    if (! flag) error, "missing OI-FITS TARGET HDU";
    data_target_id = oifits_get_target_id(data, db);
    data_target = strtrim(oifits_get_target(data, db));
    if (is_void(this.target)) {
      if (anyof(data_target != (target = data_target(1)))) {
        error, "specify a target (too many different targets in data file)";
      }
      h_set, this, target = target;
      if (! (quiet & 4)) {
        write, format="Selecting data for target \"%s\".\n", target;
      }
    }
    j = where(strglob(this.target, data_target, case=0, path=3, esc=0));
    if (! is_array(j)) {
       if (! (quiet & 4)) {
         write, format="WARNING - No data for target \"%s\" in this file.",
           this.target;
       }
       continue;
    }
    target_id = data_target_id(j);

    /* Explicitely declare as local the variables which are shared with
       _mira_grow_freqlist: */ local u_list, v_list, freq_tol;

    /* Get parameters from MiRA handle. */
    eq_nocopy, u_list, this.u;
    eq_nocopy, v_list, this.v;
    global_flags = this.flags; /* overall flags */
    flux_weight = this.flux_weight;

    /* Load all datablocks. */
    db_count = 0;
    for (db = oifits_first(data) ; db ; db = oifits_next(data, db)) {

      ++db_count;
      if (! oifits_is_data(db)) continue;

      /* Target selection -- the test is done in a very pedestrian way, but
         this is the only sure one since there is no rules for the min/max
         target ID. */
      temp = oifits_get_target_id(data, db);
      target_select = array(int, dimsof(temp));
      for (j = numberof(target_id); j >= 1; --j) {
        target_select |= (temp == target_id(j));
      }
      target_select = where(target_select);
      if (! is_array(target_select)) continue;

      /* Wavelength selection. */
      data_eff_wave = oifits_get_eff_wave(data, db);
      if (is_void(this.eff_wave)) {
        avg_wave = avg(data_eff_wave);
        _mira_warn, ("taking average data wavelength (" +
                     swrite(format="%.3f microns", avg_wave/MIRA_MICRON) +
                     ") as effective wavelength for reconstruction");
        h_set, this, eff_wave=avg_wave;
      }
      wave_select = where(abs(data_eff_wave - this.eff_wave)
                          <= 0.5*this.eff_band);
      if (numberof(wave_select) == 1) {
        wave_select = wave_select(1);
      } else if (is_array(wave_select)) {
        wave_select = wave_select(-,);
      } else {
        if (! quiet) {
          _mira_warn, swrite(format=("skipping datablock number %d: "+
                                     "no matching wavelengths"), db_count);
        }
        continue;
      }

      /* Selection index in 2D array. */
      stride = numberof(oifits_get_time(data, db));
      select = (target_select + stride*(wave_select - 1))(*);

      /* Create MiRA datablock child object. */
      wavelength = data_eff_wave(wave_select);
      w = wavelength(-:1:numberof(target_select),)(*);
      one_over_wavelength = 1.0/wavelength;
      freq_tol = this.base_tol/this.eff_wave;
      type = oifits_get_type(db);
      if (type == OIFITS_TYPE_VIS) {
        local u, v, w;
        u = (one_over_wavelength*oifits_get_ucoord(data, db)(target_select))(*);
        v = (one_over_wavelength*oifits_get_vcoord(data, db)(target_select))(*);
        amp    = oifits_get_visamp(data, db)(select);
        amperr = oifits_get_visamperr(data, db)(select);
        phi    = oifits_get_visphi(data, db)(select)*MIRA_DEGREE;
        phierr = oifits_get_visphierr(data, db)(select)*MIRA_DEGREE;
        if (cleanup_bad_data) {
          good = ((amperr > 0.0)&(phierr > 0.0));
          if (cleanup_bad_data > 1) good &= ((amp >= 0.0)&(amp <= 1.0));
          good = where(good);
          if (! is_array(good)) continue;
          u      = u(good);
          v      = v(good);
          w      = w(good);
          amp    = amp(good);
          amperr = amperr(good);
          phi    = phi(good);
          phierr = phierr(good);
        }
        temp = mira_polar_to_cartesian(amp, amperr, phi, phierr,
                                       "visibility", goodman=goodman);
        child = _mira_get_datablock(this, "vis", MIRA_COMPLEX_VISIBILITY);
        h_grow, child, flatten = 1,
          "u",      u,
          "v",      v,
          "w",      w,
          "amp" ,   amp,
          "amperr", amperr,
          "phi",    phi,
          "phierr", phierr,
          "re",     temp.re,
          "im",     temp.im,
          "wrr",    temp.wrr,
          "wri",    temp.wri,
          "wii",    temp.wii;
        // FIXME: flux_weight += sum(child.amp_weight*(child.amp_data)^2);
      } else if (type == OIFITS_TYPE_VIS2) {
        local u, v, w;
        u = (one_over_wavelength*oifits_get_ucoord(data, db)(target_select))(*);
        v = (one_over_wavelength*oifits_get_vcoord(data, db)(target_select))(*);
        vis2data = oifits_get_vis2data(data, db)(select);
        vis2err  = oifits_get_vis2err(data, db)(select);
        if (cleanup_bad_data) {
          good = (vis2err > 0.0);
          if (cleanup_bad_data > 1) good &= ((vis2data >= 0.0)&(vis2data <= 1.0));
          good = where(good);
          if (! is_array(good)) continue;
          u        = u(good);
          v        = v(good);
          w        = w(good);
          vis2data = vis2data(good);
          vis2err  = vis2err(good);
        }
        child = _mira_get_datablock(this, "vis2", MIRA_POWER_SPECTRUM);
        h_grow, child, flatten = 1,
          "u",        u,
          "v",        v,
          "w",        w,
          "vis2data", vis2data,
          "vis2err",  vis2err,
          "vis2wght", mira_stdev_to_weight(vis2err);
        // FIXME: flux_weight += 4.0*sum(child.vis2wght*(child.vis2data)^2);
      } else if (type == OIFITS_TYPE_T3) {
        local u1, u2, v1, v2, w;
        u1 = (one_over_wavelength*oifits_get_u1coord(data, db)(target_select))(*);
        v1 = (one_over_wavelength*oifits_get_v1coord(data, db)(target_select))(*);
        u2 = (one_over_wavelength*oifits_get_u2coord(data, db)(target_select))(*);
        v2 = (one_over_wavelength*oifits_get_v2coord(data, db)(target_select))(*);
        amp    = oifits_get_t3amp(data, db)(select);
        amperr = oifits_get_t3amperr(data, db)(select);
        phi    = oifits_get_t3phi(data, db)(select)*MIRA_DEGREE;
        phierr = oifits_get_t3phierr(data, db)(select)*MIRA_DEGREE;
        if (cleanup_bad_data) {
          good = ((amperr > 0.0)&(phierr > 0.0));
          if (cleanup_bad_data > 1) good &= ((amp >= 0.0)&(amp <= 1.0));
          good = where(good);
          if (! is_array(good)) continue;
          u1     = u1(good);
          v1     = v1(good);
          u2     = u2(good);
          v2     = v2(good);
          w      = w(good);
          amp    = amp(good);
          amperr = amperr(good);
          phi    = phi(good);
          phierr = phierr(good);
        }
        temp = mira_polar_to_cartesian(amp, amperr, phi, phierr,
                                       "bispectrum", goodman=goodman);
        child = _mira_get_datablock(this, "vis3", MIRA_BISPECTRUM);
        h_grow, child, flatten = 1,
          "u1",       u1,
          "v1",       v1,
          "u2",       u2,
          "v2",       v2,
          "w",        w,
          "t3amp",    amp,
          "t3amperr", amperr,
          "t3phi",    phi,
          "t3phierr", phierr,
          "t3re",     temp.re,
          "t3im",     temp.im,
          "wrr",      temp.wrr,
          "wri",      temp.wri,
          "wii",      temp.wii;
        h_pop, child, "cl_matrix";
        h_pop, child, "cl_weight";
        // FIXME: flux_weight += 9.0*sum(child.t3amp_weight*(child.t3amp_data)^2);
      } else {
        _mira_warn, "unsupported OIFITS data type";
        continue;
      }
    }
  }

  /* Update MiRA handle. */
  return h_set(this,
               u = u_list, v = v_list,
               flux_weight = flux_weight,
               flags = global_flags,
               update_pending = 1n);
}

func _mira_get_datablock(this, key, type)
{
  if (! h_has(this, key)) {
    h_set, this, key, h_new(type = type);
  }
  h_set, this, update_pending = 2;
  return h_get(this, key);
}

func _mira_build_coordinate_list(master)
/* DOCUMENT _mira_build_coordinate_list, master;
  
     Build or rebuild global list of "unique" sampled coordinates in MiRA
     opaque object MASTER.  This must be done after any addition/removal of
     data and prior to any attempt of image reconstruction.  Normally this
     operation is automatically triggered by mira_update (which see).
  
   SEE ALSO: mira_new, mira_config, mira_add_oidata, mira_update.
 */
{
  local u_list, v_list, w_list;

  /*
  ** Destroy previous information.
  */

  h_set, master, update_pending = 2;
  h_delete, master, "u", "v", "w";


  /*
  ** Collect all coordinates for all different types of data.
  */

  collect = _mira_build_coordinate_list_pass1;

  /* Complex visibilities. */
  key = "vis";
  db = h_get(master, key);
  if (db) {
    collect, master, key, "idx", "sgn", db.u, db.v, db.w;
  }

  /* Powerspectrum data. */
  key = "vis2";
  db = h_get(master, key);
  if (db) {
    collect, master, key, "idx", "sgn", db.u, db.v, db.w;
  }

  /* Bispectrum data. */
  key = "vis3";
  db = h_get(master, key);
  if (db) {
    collect, master, key, "idx1", "sgn1", db.u1, db.v1, db.w;
    collect, master, key, "idx2", "sgn2", db.u2, db.v2, db.w;
    collect, master, key, "idx3", "sgn3",
      -(db.u1 + db.u2), -(db.v1 + db.v2), db.w;
  }

  /* Total number of sampled coordinates before reduction. */
  number = numberof(master.u);
  if (numberof(master.v) != number || numberof(master.w) != number) {
    error, "incompatible number of coordinates (BUG)";
  }
  if (number < 1) {
    return;
  }

  /* Figure out whether or not we are in "monochromatic" mode. */
  w_digit = mira_digitize(master.w, master.wave_tol);
  number_of_wavelengths = numberof(w_digit.value);
  if (number_of_wavelengths > 1) {
    monochromatic = (master.monochromatic_option ? 1n : 0n);
    w = w_digit.value;
  } else {
    monochromatic = 1n;
    w = w_digit.value(1);
  }
  h_set, master, w = w, monochromatic = monochromatic;


  /*
  ** Make a list of "unique" coordinates using a *slow* O(N^2) algorithm.
  */

  if (monochromatic) {

    local u_inp, v_inp;
    eq_nocopy, u_inp, master.u;
    eq_nocopy, v_inp, master.v;
    number = numberof(u_inp);
    u_out = array(double, number); /* maximum size */
    v_out = array(double, number); /* maximum size */
    n_out = array(long, number); /* maximum size */
    idx = array(long, number);
    sgn = array(long, number);
    mid_wavelength = 0.5*(max(master.w) + min(master.w));
    freq_tol = master.base_tol/mid_wavelength;

    j = k = 1;
    u_out(k) = u_inp(j);
    v_out(k) = v_inp(j);
    n_out(k) = 1;
    idx(j) = k;
    sgn(j) = 1;
    while (++j <= number) {

      /* Get j-th position. */
      u = u_inp(j);
      v = v_inp(j);

      /* Search +/-position among list of positions. */
      u_tmp = u_out(1:k);
      v_tmp = v_out(1:k);
      rp = (temp = u - u_tmp)*temp + (temp = v - v_tmp)*temp;
      rn = (temp = u + u_tmp)*temp + (temp = v + v_tmp)*temp;
      rp_min = min(rp);
      rn_min = min(rn);
      if (min(rp_min, rn_min) > freq_tol) {
        /* Got a new position. */
        idx(j) = ++k;
        sgn(j) = 1;
        u_out(k) = u;
        v_out(k) = v;
        n_out(k) = 1;
      } else if (rp_min <= rn_min) {
        idx(j) = (kp = rp(mnx));
        sgn(j) = 1;
        np1 = (n = n_out(kp)) + 1;
        u_out(kp) = (n*u_out(kp) + u)/np1;
        v_out(kp) = (n*v_out(kp) + v)/np1;
        n_out(kp) = np1;
      } else {
        idx(j) = (kp = rn(mnx));
        sgn(j) = 1;
        np1 = (n = n_out(kp)) + 1;
        u_out(kp) = (n*u_out(kp) - u)/np1;
        v_out(kp) = (n*v_out(kp) - v)/np1;
        n_out(kp) = np1;
      }
    }
    if (k < number) {
      u_out = u_out(1:k);
      v_out = v_out(1:k);
      n_out = n_out(1:k);
    }
    write, format="There are %d sampled frequencies out of %d measurements.\n",
        k, number;

  } else {
    error, "only monochromatic mode is implemented by MiRA";
  }


  /*
  ** Store the global list of coordinates and indirection tables.
  */

  db = master.vis;
  if (db) {
    /* Complex visibilities. */
    h_set, db, idx = idx(db.idx), sgn = sgn(db.idx)*db.sgn;
  }
  db = master.vis2;
  if (db) {
    /* Powerspectrum data. */
    h_set, db, idx = idx(db.idx), sgn = sgn(db.idx)*db.sgn;
  }
  db = master.vis3;
  if (db) {
    /* Bispectrum data. */
    h_set, db,
      idx1 = idx(db.idx1), sgn1 = sgn(db.idx1)*db.sgn1,
      idx2 = idx(db.idx2), sgn2 = sgn(db.idx2)*db.sgn2,
      idx3 = idx(db.idx3), sgn3 = sgn(db.idx3)*db.sgn3;
  }

  return h_set(master, update_pending = 1, u = u_out, v = v_out);
}

func _mira_build_coordinate_list_pass1(this, key, idx, sgn, u, v, w)
{
  local u_list, v_list, w_list;

  eq_nocopy, u_list, this.u;
  eq_nocopy, v_list, this.v;
  eq_nocopy, w_list, this.w;

  offset = numberof(u_list);
  dims = dimsof(u, v, w);
  if (is_void(dims)) {
    error, "non-conformable coordinates";
  }
  for (number = 1, k = numberof(dims); k >= 2; --k) {
    number *= dims(k);
  }

  /* Spatial frequency sign to keep only 1/2 (u,v) plane -- same
     choice as FFTW. */
  s = double(1 - 2*((u < 0.0) | ((u == 0.0)&(v < 0.0))));
  // FIXME: check geometry
  h_set, this(key), idx, offset + 1 : offset + number, sgn, s;

  /* Make sure S and W have the correct number of elements and
     dimension lists. */
  if (numberof(s) != number) {
    s += array(double, dims);
  }
  if (numberof(w) != number) {
    w += array(double, dims);
  }

  /* Append coordinates to external lists and fix the sign of the
     spatial frequencies. */
  grow, u_list, (s*u)(*);
  grow, v_list, (s*v)(*);
  grow, w_list, w(*);
  h_set, this, u = u_list, v = v_list, w = w_list;
}

/*---------------------------------------------------------------------------*/
/* CONFIGURATION */

func mira_config(this, pixelsize=, dim=, xform=)
/* DOCUMENT mira_config, this, keyword=value, ...;
 *
 *   Configure or change parameters of MiRA opaque handle THIS.  All
 *   configurable parameters are specified by keywords (see below).  Return
 *   THIS when called as a function.
 *
 * KEYWORDS:
 *   PIXELSIZE = size of pixel in radians;
 *   DIM       = number of pixels accross the field of view;
 *   XFORM     = method to compute the Fourier transform:
 *     "exact"    to use exact (but slow) transform;
 *     "nfft"     to use nonequispaced FFT;
 *     "fft"      to use FFT(W) followed by linear interpolation.
 *
 *  SEE ALSO: mira_new, mira_solve.
 */
{
  update = 0n;
  if (is_void(dim)) {
    dim = this.dim;
    if (is_void(dim)) {
      error, "DIM must be specified";
    }
  } else {
    if (mira_get_one_integer(dim) || dim < 2) {
      error, "DIM must be an integer greater or equal 2";
    }
    if (this.dim != dim) {
      update = 1n;
    }
  }
  if (is_void(pixelsize)) {
    pixelsize = this.pixelsize;
    if (is_void(pixelsize)) {
      error, "PIXELSIZE must be specified";
    }
  } else {
    if (mira_get_one_real(pixelsize) || pixelsize <= 0.0) {
      error, "PIXELSIZE must be a strictly positive real";
    }
    if (this.pixelsize != pixelsize) {
      update = 1n;
    }
  }
  fov = pixelsize*(dim - 1.0);

  if (is_void(xform)) {
    xform = this.xform_name;
    if (is_void(xform)) {
      xform = "exact";
      update = 1n;
    }
  } else {
    if (! is_scalar(xform) || ! is_string(xform)) {
      error, "XFORM must be a scalar string";
    }
    if (this.xform_name != xform) {
      update = 1n;
    }
  }
  if (update) {
    h_set, this, update_pending=1n, dim=dim, fov=fov,
      pixelsize=pixelsize, xform_name=xform;
  }
  return this;
}

func mira_update(this)
/* DOCUMENT mira_update, this;
     Updates internals of MiRA instance THIS.  This function is normally
     automatically called whenever any parameters of THIS have changed which
     require to recompute some cached values (most importantly the
     coefficients of the Fourier transform).
  
   SEE ALSO: mira_config.
 */
{
  if (this.update_pending > 1) {
    _mira_build_coordinate_list, this;
  }

  /* Get sky coordinates (in radians). */
  if (! (dim = this.dim) ||
      ! (fov = this.fov) ||
      ! (pixelsize = this.pixelsize)) {
    error, "FOV, DIM and PIXELSIZE must be set before (see mira_config)";
  }

  /* Compute transform. */
  if (this.update_pending || is_void(this.xform)) {
    if (is_void(this.xform_name)) {
      h_set, this, xform_name="exact";
    }
    if (this.xform_name == "exact") {
      xform =  mira_new_exact_xform(this.u, this.v,
                                    pixelsize, dim, dim);
    } else if (this.xform_name == "nfft") {
      xform =  mira_new_nfft_xform(this.u, this.v,
                                   pixelsize, dim, dim);
    } else if (this.xform_name == "fft") {
      xform =  mira_new_fft_xform(this.u, this.v,
                                  pixelsize, dim, dim);
    } else {
      error, "unknown value for XFORM";
    }
    return h_set(this, xform=xform, update_pending=0n);
  }
}

/*---------------------------------------------------------------------------*/
/* EXACT FOURIER TRANSFORM */

local mira_new_exact_xform, _mira_apply_exact_xform;
/* DOCUMENT xform = mira_new_exact_xform(u, v, pixelsize, nx, ny)
 *
 *   Creates a linear operator to compute the "exact" linear transform
 *   between an image and the measured complex visibilities.  Arguments U
 *   and V give the coordinates of the measured spatial frequencies,
 *   PIXELSIZE is the pixel size in the image plane, NX and NY are the
 *   number of pixels along the two first image dimensions.  The
 *   coordinates U and V must be vectors of same length and
 *   NFREQS=numberof(U)=numberof(V) is the number of measured complex
 *   visibilities.
 *
 *   The returned operator can be used as follows:
 *
 *      vis = XFORM(img)
 *
 *   to compute the model visibilities VIS, such that VIS(1,..)  and
 *   VIS(2,..) are respectively the real and imaginary parts of the complex
 *   visibilities.  The transpose of the operator can also be applied (for
 *   instance to compute the gradient of the likelihood):
 *
 *      XFORM(inp, 1)
 *
 *   where input array INP is a 2-by-NFREQS array, with NFREQS the number
 *   of complex visibilities.
 *
 *
 * SEE ALSO: mira_update, mira_new_fft_xform.
 */

func mira_new_exact_xform(u, v, pixelsize, nx, ny)
{
  /* The argument of the complex exponent in the Fourier transform is:
   *
   *     Q = -2*PI*(RA*U + DEC*V)
   *
   * FIXME:
   * where RA is the relative right ascension and DEC is the relative
   * declination.  The relationships with image coordinates (X,Y) are:
   *
   *     RA  = X
   *     DEC = Y
   *
   * hence:
   *
   *   Q = -2*PI*(U*X + V*Y)
   */
  if (! is_vector(u) || ! is_vector(v) || numberof(u) != numberof(v)) {
    error, "arguments U and V must be vectors of same length";
  }
  x = pixelsize*(indgen(nx) - (nx + 1)/2.0);
  y = pixelsize*(indgen(ny) - (ny + 1)/2.0);
  q = (-2.0*MIRA_PI)*(u*x(-,..) + v*y(-,-,..));
  a = array(double, 2, dimsof(q));
  a(1,..) = cos(q);
  a(2,..) = sin(unref(q));
  obj = h_new(a=a);
  h_evaluator, obj, "_mira_apply_exact_xform";
  return obj;
}

func _mira_apply_exact_xform(this, x, job)
{
  return mvmult(this.a, x, job);
}

/*---------------------------------------------------------------------------*/
/* NONEQUISPACED FAST FOURIER TRANSFORM */

local mira_new_nfft_xform, _mira_apply_nfft_xform;
/* DOCUMENT xform = mira_new_nfft_xform(u, v, pixelsize, nx, ny)
 *
 *   Creates a linear operator to compute the "exact" linear transform
 *   between an image and the measured complex visibilities.  Arguments U
 *   and V give the coordinates of the measured spatial frequencies,
 *   PIXELSIZE is the pixel size in the image plane, NX and NY are the
 *   number of pixels along the two first image dimensions.  The
 *   coordinates U and V must be vectors of same length and
 *   NFREQS=numberof(U)=numberof(V) is the number of measured complex
 *   visibilities.
 *
 *   The returned operator can be used as follows:
 *
 *      vis = XFORM(img)
 *
 *   to compute the model visibilities VIS, such that VIS(1,..)  and
 *   VIS(2,..) are respectively the real and imaginary parts of the complex
 *   visibilities.  The transpose of the operator can also be applied (for
 *   instance to compute the gradient of the likelihood):
 *
 *      XFORM(inp, 1)
 *
 *   where input array INP is a 2-by-NFREQS array, with NFREQS the number
 *   of complex visibilities.
 *
 *
 * SEE ALSO: mira_update, mira_new_exact_xform, mira_new_fft_xform.
 */

func mira_new_nfft_xform(u, v, pixelsize, nx, ny)
{
  if (! is_func(nfft_new)) {
    include, "nfft.i", 1;
  }
  if (! is_vector(u) || ! is_vector(v) || numberof(u) != numberof(v)) {
    error, "arguments U and V must be vectors of same length";
  }

  /* The coordinate system is slighlty different between NFFT and MiRA.
   *
   * In MiRA, (0,0) is at the geometrical center of the field of view (FOV):
   *   x = pixelsize*(indgen(nx) - (nx + 1)/2.0);
   *   y = pixelsize*(indgen(ny) - (ny + 1)/2.0);
   *
   * In NFFT:
   *   x = [-nx/2, 1-nx/2, ..., nx/2-1]*pixelsize
   *   y = [-ny/2, 1-ny/2, ..., ny/2-1]*pixelsize
   * plus NX and NY must be even.
   */
  local r1, r2;
  if (nx % 2 == 1) {
    n1 = nx + 1;
    r1 = 2 : nx;
  } else {
    n1 = nx;
  }
  if (ny % 2 == 1) {
    n2 = ny + 1;
    r2 = 2 : ny;
  } else {
    n2 = ny;
  }
  nodes = [u*pixelsize, v*pixelsize];
  dims = [n1, n2];
  obj = h_new(nfft = nfft_new(dims, nodes),
              n1 = n1, r1 = r1,
              n2 = n2, r2 = r2,
              sub = (n1 != nx || n2 != ny));
  h_evaluator, obj, "_mira_apply_nfft_xform";
  return obj;
}

func _mira_apply_nfft_xform(this, x, job)
{
  local z;
  if (job) {
    /* adjoint operator */
    z = this.nfft(mira_cast_real_as_complex(x), 1n);
    if (this.sub) {
      return double(z)(this.r1, this.r2);
    } else {
      return double(z);
    }
  } else {
    /* direct operator */
    if (this.sub) {
      tmp = array(complex, this.n1, this.n2);
      tmp(this.r1, this.r2) = x;
      eq_nocopy, tmp, x;
    }
    reshape, z, &this.nfft(x), double, 2, this.nfft.num_nodes;
    return z;
  }
}

/*---------------------------------------------------------------------------*/
/* FFT BASED FOURIER TRANSFORM */

local mira_new_fft_xform, _mira_apply_fft_xform, _mira_fft_xform_builder;
/* DOCUMENT xform = mira_new_fft_xform(u, v, pixelsize, nx, ny)

     Creates a linear operator to approximate the linear transform between an
     image and the measured complex visibilities using FFT and interpolation
     of the frequencies.  Otherwise, these functions behave the same as
     mira_new_exact_xform (which see).


   SEE ALSO: mira_update, mira_new_exact_xform.
 */

func mira_new_fft_xform(u, v, pixelsize, nx, ny)
{
  /* Build FFT operator.  The complex visibility array is seen as a
   * 2-by-NFREQS array of reals:
   *    VIS(1,) = real part
   *    VIS(2,) = imaginary part
   *
   * Output of FFT is a STRIDE -by- NY array where:
   *    STRIDE = NX         with FFT
   *    STRIDE = (NX/2 + 1) with FFTW
   */
  if (nx % 2 != 0 || ny % 2 != 0) {
    /* Even dimensions are required to allow for the trick of rolling the
       image by simply multiplying the FFT spectrum by +/-1 (see remarks in
       _mira_fft_xform_builder). */
    error, "dimensions must be even numbers";
  }
  dims = [2, nx, ny];
  if (is_func(fftw)) {
    stride = nx/2 + 1; /* number of "positive" frequencies along X */
    f = linop_new_fftw(dims=dims, real=1);
    half = 1n;
  } else {
    stride = nx;
    f = linop_new_fft(dims, real=1);
    half = 0n;
  }

  /* Convert (U,V) into frequel units -- FFT frequency sampling is
     1/(N*PIXELSIZE). */
  if (! is_vector(u) || ! is_vector(v)
      || (nfreqs = numberof(u)) != numberof(v)) {
    error, "arguments U and V must be vectors of same length";
  }
  u *= (pixelsize*nx); /* RA  = X */
  v *= (pixelsize*ny); /* DEC = Y */

  /* Compute integer bounding box of frequencies such that:
   *   U0 <= U < U0 + 1
   *   V0 <= V < V0 + 1
   */
  u0 = floor(u);
  v0 = floor(v);

  /* Check that there is no aliasing. */
  umax = (nx - 1)/2; /* Nyquist frequency along U */
  vmax = (ny - 1)/2; /* Nyquist frequency along V */
  if (min(u0) < -umax || max(u0) >= +umax ||
      min(v0) < -vmax || max(v0) >= +vmax) {
    error, "pixel size must be reduced to avoid aliasing";
  }

  /* Compute the weights (and indices) of the bilinear interpolation.
   *   W(c,1,k) = weight of c-th corner for real part of k-th
   *              measured frequency
   *   W(c,2,k) = ditto for imaginary part.
   *
   * Depending on the position of the interpolated frequency:
   *   W(c,2,k) = +W(c,1,k) if conjugate not taken;
   *   W(c,2,k) = -W(c,1,k) if conjugate taken.
   */
  i1 = (i0 = long(u0)) + 1;
  j1 = (j0 = long(v0)) + 1;
  w00 = 1.0 - (w01 = u - u0);
  w10 = 1.0 - (w11 = v - v0);
  w = array(double, 4, 2, nfreqs); /* array of weights */
  j = array(long, 4, 2, nfreqs); /* array of indices in input space */
  builder = _mira_fft_xform_builder; /* shortcut */
  builder, w, j, 1, w00*w10, nx, i0, ny, j0, stride;
  builder, w, j, 2, w01*w10, nx, i1, ny, j0, stride;
  builder, w, j, 3, w00*w11, nx, i0, ny, j1, stride;
  builder, w, j, 4, w01*w11, nx, i1, ny, j1, stride;
  w00 = w01 = w10 = w11 = i0 = j0 = i1 = j1 = 0; /* free some memory */
  i = array(long, 4, 2, nfreqs); /* array of indices in output space */
  i(*) = indgen(1:2*nfreqs)(-:1:4,)(*);

  /* Get rid of zero-weights (if any) and make sparse interpolation matrix. */
  k = where(w);
  if (numberof(k) < numberof(w)) {
    w = w(k);
    i = i(k);
    j = j(k);
  }
  s = sparse_matrix(w, [2, 2, nfreqs], i, [3, 2, stride, ny], j);

  /* Build indirection tables for U ~ (U1,U2) and V ~ (-U1,-U2) for making the
     gradient hermitian --- in the notation U ~ (U1,U2) means: index U
     coresponding to spatial frequency (U1,U2). */
  n1 = nx;
  n2 = ny;
  local idx0, idx1, idx2;
  if (half) {
    /* FFTW case.  Only zero-th and, maybe, Nyquist frequencies along the
       first dimension has possibly a negative counterpart in the same
       array. */
    stride = n1/2 + 1;

    idx0 = 0; /* (0,0) must be real */
    if (n1 >= 2 && n1 % 2 == 0) {
      grow, idx0, (n1/2); /* (n1/2,0) must be real */
    }
    if (n2 >= 2 && n2 % 2 == 0) {
      grow, idx0, stride*(n2/2); /* (0,n2/2) must be real */
      if (n1 >= 2 && n1 % 2 == 0) {
        grow, idx0, (n1/2) + stride*(n2/2); /* (n1/2,n2/2) must be real */
      }
    }
    ++idx0;
    if (n2 >= 3) {
      idx1 = stride*indgen(n2 - (n2/2) - 1);
      idx2 = stride*indgen(n2 - 1 : (n2/2) + 1 : -1);
      if (n1 >= 2 && n1 % 2 == 0) {
        grow, idx1, idx1 + (n1/2);
        grow, idx2, idx2 + (n1/2);
      }
      ++idx1;
      ++idx2;
    }
  } else {
    /* FFT case. */
    stride = n1;
    u1 = indgen(0:n1-1);
    v1 = (n1 - u1)%n1;
    u2 = indgen(0:n2-1);
    v2 = (n2 - u2)%n2;
    u = 1 + u1 + (stride*u2)(-,); /* + 1 because Yorick index */
    v = 1 + v1 + (stride*v2)(-,); /* + 1 because Yorick index */
    j = where(u == v);
    idx0 = (is_array(j) ? u(j) : []);
    j = where(u < v);
    idx1 = (is_array(j) ? u(j) : []);
    idx2 = (is_array(j) ? v(j) : []);
  }

  /* Build object. */
  obj = h_new(f=f, s=s, idx0=idx0, idx1=idx1, idx2=idx2);
  h_evaluator, obj, "_mira_apply_fft_xform";
  return obj;
}

func _mira_fft_xform_builder(wa, ia, c, w, nx, kx, ny, ky, stride)
{
  /* Change sign of frequency (KX,KY) where conjugate complex visibility is
     taken; that is, where KX < 0 so that interpolation works for FFT and for
     FFTW (real-complex transform). */
  if (sign(0) != 1) error, "assertion failed";
  sx = sign(kx); /* -1 where conjugate frequency is taken, +1 elsewhere */
  kx = sx*kx;
  ky = sx*ky;

  /* Convert (KX,KY) into zero-based indices in FFT array (already done above
     for KX) and compute offsets in array of pairs of reals (re,im) --- hence
     the factor of 2. */
  ky = ky + ny*(ky < 0);
  offset = 2*(kx + stride*ky);

  /* Multiply the weights by +/-1 to account for a shift by (NX/2,NY/2),
     i.e. assume image coordinate origin is the center of the field of
     view --- this trick only works for arrays with even dimensions. */
  w *= (1 - 2*(kx%2))*(1 - 2*(ky%2));

  /* Set weights and input indices of the linear transform. */
  wa(c, 1, ) =    w; /* weights for real part */
  wa(c, 2, ) = sx*w; /* weights for imaginary part */
  ia(c, 1, ) = 1 + offset; /* indices for real part */
  ia(c, 2, ) = 2 + offset; /* indices for imaginary part */
}

MIRA_MAKE_HERMITIAN = 0n;

func _mira_apply_fft_xform(this, x, job)
{
  if (! job) {
    /* direct transform */
    return this.s(mira_cast_complex_as_real(this.f(x)));
  } else if (job == 1) {
    /* transpose transform for the gradient */
    temp = this.s(x, 1);
    if (MIRA_MAKE_HERMITIAN) {
      /* FIXME: only work for 2-D images */
      if (is_array(this.idx0)) {
        temp(2,this.idx0) = 0.0;
      }
      if (is_array(this.idx1)) {
        wre = (temp(1,this.idx1) + temp(1,this.idx2));
        wim = (temp(2,this.idx1) + temp(2,this.idx2));
        temp(1,this.idx1) =  wre;
        temp(1,this.idx2) =  wre;
        temp(2,this.idx1) =  wim;
        temp(2,this.idx2) = -wim;
      }
    }
    return this.f(mira_cast_real_as_complex(temp), 1);
  } else {
    error, "unsupported value for JOB";
  }
}

/*---------------------------------------------------------------------------*/

func mira_plot(master, model, what)
{
  // FIXME: optimize

  /* Apply linear transform to compute the model of the complex visibility
     of image MODEL for _all_ the measured frequencies. */
  vis = master.xform(model);
  vis_re = vis(1,..);
  vis_im = vis(2,..);
  vis_amp = abs(vis_re, vis_im);
  vis_phi = atan(vis_im, vis_re);
  nil = []; /* for cleanup */


  /* Complex visibilities. */
  key = "vis";
  if (is_void(what) || what == key) {
    db = h_get(master, key);
    if (db) {
      f = master.freq_len(db.idx);
      plp, db.amp, f, dy=db.amperr, ticks=1, color="magenta", symbol=1, size=0.3, fill=1;
      plp, vis_amp(db.idx), f, color="blue", symbol=3, size=0.4, fill=1;
      xytitles, "spatial frequency", "amplitude";
    }
  }

  /* Powerspectrum data. */
  if (is_void(what) || strpart(what, 1:4) == "vis2") {
    local idx, amp, amperr, relerr;
    if (what == "vis2-2D") {
      /* Collect data from powerspectrum. */
      db = h_get(master, "vis2");
      if (db) {
        temp = vis_amp(db.idx);
        grow, relerr, (db.vis2data - temp*temp)/db.vis2err;
        grow, idx, db.idx;
      }

      /* Collect data from complex visibilities. */
      db = h_get(master, "vis");
      if (db) {
        eq_nocopy, amp, db.amp;
        eq_nocopy, amperr,  db.amperr;
        if (is_void(amp)) {
          amp = abs(db.re, db.im);
          amperr = sqrt(0.5/db.wrr + 0.5/db.wii);
          h_set, db, amp=amp, amperr=amperr;
        }
        sel = where(amperr > 0.0);
        if (is_array(sel)) {
          idx_sel = db.idx(sel);
          grow, relerr, (amp(sel) - vis_amp(idx_sel))/amperr(sel);
          grow, idx, idx_sel;
        }
      }

      /* Plot 2-D amplitude. */
      z = abs2(fft(model)); // FIXME: save time: use FFTW, the FFT of the model may already be computed
      if (z(1) > 0.0) z = log(z + 1e-6*z(1));
      fft_pli, z, scale=[-1.0,1.0]/master.fov; // FIXME: this is weird
      mira_fix_image_axis;
      if (! is_void(idx)) {
        u = master.u(idx);
        v = master.v(idx);
        color_map = [char(indgen(0:255)),
                     char(indgen(255:0:-1)),
                     array(char, 256)];
        color_index = min(long((255/3.0)*abs(relerr) + 1.5), 256);
        for (i = numberof(color_index); i >= 1; --i) {
          plp, v(i)*[1,-1], u(i)*[-1,1], color=color_map(color_index(i),),
            symbol=1, size=0.1;
        }
      }
    } else if (what == "vis2") {
      /* Collect data from complex visibilities. */
      db = h_get(master, "vis");
      if (db) {
        eq_nocopy, idx, db.idx;
        eq_nocopy, amp, db.amp;
        eq_nocopy, amperr,  db.amperr;
        if (is_void(amp)) {
          // FIXME:
          amp = abs(db.re, db.im);
          amperr = sqrt(0.5/db.wrr + 0.5/db.wii);
          h_set, db, amp=amp, amperr=amperr;
        }
      }

      /* Collect data from powerspectrum. */
      db = h_get(master, "vis2");
      if (db) {
        sel = where(db.vis2data > 0.0);
        if (is_array(sel)) {
          temp = sqrt(db.vis2data(sel));
          grow, idx, db.idx(sel);
          grow, amp, temp;
          grow, amperr, 0.5*db.vis2err(sel)/temp;
        }
      }

      /* Plot 1-D amplitude. */
      if (! is_void(idx)) {
        f = master.freq_len(idx);
        plp, amp, f, dy=amperr, ticks=1, color="red", symbol=1,
          size=0.3, fill=1;
        plp, vis_amp(idx), f, color="green", symbol=3, size=0.3, fill=1;
        xytitles, "spatial frequency", "amplitude";
      }
    }

  }

  /* Bispectrum phase. */
  key = "vis3";
  if (is_void(what) || what == key) {
    db = h_get(master, key);
    if (db && ! master.zap_phase) {

      /* normalized residuals */
      res = arc(db.t3phi - (db.sgn1*vis_phi(db.idx1) +
                            db.sgn2*vis_phi(db.idx2) +
                            db.sgn3*vis_phi(db.idx3)))/sqrt(db.t3phierr);

      binsize = 0.2;
      index = long(floor((1.0/binsize)*res + 0.5));
      inf = min(index) - 1;
      sup = max(index) + 1;
      hy = histogram(index + (1 - inf), top=(sup - inf + 1));
      hx = binsize*indgen(inf:sup);

      plh, hy, hx, marks=0;
      //require, "modulo.i";
      //vis_CP = modulo(vis_CP+pi,2*pi)-pi;
#if 0
      vis_CP = arc(vis_CP);

      plp, CP*180./pi, dy = CPerr*180./pi, ticks=1, color="red", symbol=1, size=0.3, fill=1;
      plp, vis_CP*180./pi, color="green", symbol=3, size=0.3, fill=1;
      flag2 = 1;
#endif
    }
  }



  return; // FIXME: cleanup

  if (! is_void(img)) {
    vis = master.xform(img);
    re = vis(1,..);
    im = vis(2,..);
    old_ps = h_pop(master, "ps");
    h_set, master, ps=(re*re + im*im);
  }
  for (db = master.data ; db ; db = db.next) {
    op = db.plot; /* work around Yorick parser bug */
    if (is_symlink(op)) {
      op = value_of_symlink(op);
    }
    if (is_func(op)) {
      op, master, db;
    }
  }
  if (! is_void(master.ps)) {
    opt = master.plot_model_options;
    plp, master.ps, master.freq_len,
      color=opt.color, symbol=opt.symbol, size=opt.size,
      fill=opt.fill, ticks=opt.ticks, type=opt.type,
      width=opt.width;
  }
  if (! is_void(img)) {
    /* restore PowerSpectrum */
    h_set, master, ps=old_ps;
  }
}

/* To compute penalty w.r.t. complex data, there are different cases:
 *   1. use quadratic penalty in (RE,IM) coordinates
 *   2. compute penalty in (AMP,PHI) coordinates
 *      a) use only amplitude data
 *      b) use only phase data
 *      c) use both
 *
 * In addition there are different kinds of data:
 *   1. complex visibilities (in polar or cartesian coordinates)
 *   2. power spectrum
 *   3. bispectrum (in polar or cartesian coordinates)
 */

func mira_data_penalty(master, model, &grd)
/* DOCUMENT mira_data_penalty(data, model, grd)
 *
 *   Compute misfit penalty w.r.t. to interferometric data.  DATA is MiRA
 *   opaque handle which stores the interferometric data.  MODEL is a 2-D
 *   or 3-D model of the brightness distribution.  GRD is an optional
 *   output variable to store the gradient of the penalty w.r.t. the model
 *   parameters.
 *
 * SEE ALSO: mira_new.
 */
{
  if (master.update_pending) mira_update, master;


  /* Apply linear transform to compute the model of the complex visibility
     of image MODEL for _all_ the measured frequencies. */
  vis = master.xform(model);
  grd = array(double, dimsof(vis));
  vis_re = vis(1,..);
  vis_im = vis(2,..);
  nil = []; /* for cleanup */


  /*
  ** Compute fitting error and gradient for every different type of
  ** data.
  */

  /* Complex visibilities. */
  key = "vis";
  err1 = 0.0;
  ndata1 = 0;
  db = h_get(master, key);
  if (db) {
    /* RE = real part of residuals
     * IM = imaginary part of residuals
     * ERR = sum(WRR*RE^2 + 2*WRI*RE*IM + WII*IM*IM)
     *     = sum((WRR*RE + WRI*IM)*RE + (WRI*RE + WII*IM)*IM
     *     = 0.5*sum((dERR/dRE)*RE) + 0.5*sum((dERR/dIM)*IM)
     * dERR/dRE = 2*(WRR*RE + WRI*IM)
     * dERR/dIM = 2*(WRI*RE + WII*IM)
     */
    local idx; eq_nocopy, idx, db.idx;
    local sgn; eq_nocopy, sgn, db.sgn;

    re = vis_re(idx) -     db.re; /* real part of residuals */
    im = vis_im(idx) - sgn*db.im; /* imaginary part of residuals */
    temp_re = db.wrr*re + db.wri*im;
    temp_im = db.wri*re + db.wii*im;
    err1 = sum(temp_re*re) + sum(temp_im*im);
    grd(1, idx) += (temp_re + temp_re);
    grd(2, idx) += (temp_im + temp_im);
    ndata1 = 2*numberof(idx);
    im = re = temp_re = temp_im = nil; /* cleanup */
  }

  /* Powerspectrum data. */
  key = "vis2";
  err2 = 0.0;
  ndata2 = 0;
  db = h_get(master, key);
  if (db && ! master.zap_amplitude) {
    /* Error and gradient w.r.t. powerspectrum data:
     *        ERR = sum(W*E^2)
     *   dERR/dPS = 2*W*E
     * where E = PS - VIS2DATA. The convertion of the gradient
     * writes:
     *   dERR/dRE = (dPS/dRE) * (dERR/dPS) = (2*RE) * (2*W*E)
     *   dERR/dIM = (dPS/dIM) * (dERR/dPS) = (2*IM) * (2*W*E)
     * finally:
     *   dERR/dRE = 4*W*E*RE
     *   dERR/dIM = 4*W*E*IM
     */
    local idx; eq_nocopy, idx, db.idx;
    re = vis_re(idx);
    im = vis_im(idx);
    e = (re*re + im*im) - db.vis2data;
    q = db.vis2wght*e;
    err2 = sum(q*e);
    ndata2 = numberof(e);
    e = nil; /* cleanup */
    q *= 4.0;
    grd(1, idx) += q*re;
    grd(2, idx) += q*im;
    im = re = nil; /* cleanup */
  }

  /* Bispectrum data. */
  key = "vis3";
  err3 = 0.0;
  ndata3 = 0;
  db = h_get(master, key);
  if (db) {
    use_phasor = 1;
    if (use_phasor && ! master.zap_phase) {

      /* Compute power-spectrum and skip all computations if there is not
         at least one non-zero model amplitude. */
      ps = vis_re*vis_re + vis_im*vis_im;
      if (max(ps) > 0.0) {

        /* Get all indirections. */
        local idx1; eq_nocopy, idx1, db.idx1;
        local idx2; eq_nocopy, idx2, db.idx2;
        local idx3; eq_nocopy, idx3, db.idx3;
        local sgn1; eq_nocopy, sgn1, db.sgn1;
        local sgn2; eq_nocopy, sgn2, db.sgn2;
        local sgn3; eq_nocopy, sgn3, db.sgn3;

        /* Compute the model phase.  Safe mode is when the model amplitude
           is non-zero everywhere. */
        safe = (min(ps) > 0.0);
        if (safe) {
          phi = atan(vis_im, vis_re);
          one_over_ps = 1.0/ps;
        } else {
          not_zero = (ps > 0.0);
          select = where(not_zero);
          phi = array(double, dimsof(ps));
          phi(select) = atan(vis_im(select), vis_re(select));
          one_over_ps = array(double, dimsof(ps));
          one_over_ps(select) = 1.0/ps(select);
        }

        /* Get/compute the weights W. */
        local w;
        if (is_void(db.cl_weight)) {
          w = mira_stdev_to_weight(db.t3phierr);
          h_set, db, cl_weight=w;
        } else {
          eq_nocopy, w, db.cl_weight;
        }
        if (! safe) {
          w *= (not_zero(idx1) & not_zero(idx2) & not_zero(idx3));
        }

        /* Compute weighted residuals E. */
        ndata3 = numberof(idx1);
        if (MIRA_SPARSE) {
          c = db.cl_matrix;
          if (is_void(c)) {
            if (MIRA_DEBUG) _mira_debug, "computing sparse phase closure matrix...";
            i = indgen(ndata3);
            c = sparse_matrix([sgn1,sgn2,sgn3],
                              dimsof(idx1), [i,i,i],
                              [1,numberof(phi)], [idx1,idx2,idx3]);
            h_set, db, cl_matrix=c;
          }
          e = c(phi) - db.t3phi;
        } else {
          e = (sgn1*phi(idx1) + sgn2*phi(idx2) + sgn3*phi(idx3)) - db.t3phi;
        }

        if (master.haniff) {
          /* Use length of arc to compute penalty.
           *   ERR = sum(W*arc(A)^2)
           *   GRD ~ 2*W*arc(A)
           * where E is the phase closure error and W is the weight.
           */
          e = arc(e);
          w = (w + w)*e;
          err3 = 0.5*sum(w*e);
        } else if (! (MIRA_FLAGS & 1)) {
          /* Use length of cord to compute penalty (this also corresponds
           * to the Von Mises distribution a.k.a. circular normal distribution
           * -- http://en.wikipedia.org/wiki/Von_Mises_distribution):
           *
           *   ERR = sum(W*(2*sin(E/2))^2)
           *       = 4*sum(W*sin(E/2)^2)
           *       = 2*sum(W*(1 - cos(E)))
           *   GRD = 4*W*sin(E/2)*cos(E/2)
           *       = 2*W*sin(E)
           *
           * where E is the phase closure error and W the weight.  Maybe
           * 1-cos(E) is more subject to rounding errors than sin(E/2)^2
           * especially for small residual errors.
           */
          if (1) {
            /* avoids rounding errors for small residuals */
            e *= 0.5;
            s = sin(e);
            w *= 4.0*s;
            err3 = sum(w*s);
            s = [];
            w *= cos(e);
          } else {
            err3 = 2.0*sum(w - w*cos(e));
            w = (w + w)*sin(e);
          }
        } else {
          /* Use convexified criterion:
           *   ERR = sum(W*sin(E)^2)
           *   GRD = 2*W*sin(E)*cos(E)
           * where E is the phase closure error and W the weight.  Maybe
           * 1-cos(E) is more subject to rounding errors than sin(E/2)^2.
           */
          s = sin(e);
          w *= s;
          err3 = sum(w*s);
          s = [];
          w = (w + w)*cos(e);
        }
        e = [];

        /* Compute gradient with respect to phase by applying transpose of
           phase closure operator. */
        if (MIRA_SPARSE) {
          w = c(w, 1);
        } else {
          top = numberof(ph);
          w = (histogram(idx1, sgn1*w, top=top) +
               histogram(idx2, sgn2*w, top=top) +
               histogram(idx3, sgn3*w, top=top));
        }

        /* Convert gradient w.r.t. phase into gradient w.r.t. real and
         * imaginary parts:
         *   dERR/dRE = (dPH/dRE) * (dERR/dPH) = (-IM/PS) * GRD_PH
         *   dERR/dIM = (dPH/dIM) * (dERR/dPH) = (+RE/PS) * GRD_PH
         */
        w *= one_over_ps;
        grd(1,..) -= w*vis_im;
        grd(2,..) += w*vis_re;

      }

    }

    if (! use_phasor) {
      /* Get all indirections. */
      local idx1; eq_nocopy, idx1, db.idx1;
      local idx2; eq_nocopy, idx2, db.idx2;
      local idx3; eq_nocopy, idx3, db.idx3;
      local sgn1; eq_nocopy, sgn1, db.sgn1;
      local sgn2; eq_nocopy, sgn2, db.sgn2;
      local sgn3; eq_nocopy, sgn3, db.sgn3;
      zmult_re = _mira_zmult_re;
      zmult_im = _mira_zmult_im;

      re1 = vis_re(idx1);
      im1 = vis_im(idx1)*sgn1;
      re2 = vis_re(idx2);
      im2 = vis_im(idx2)*sgn2;
      re3 = vis_re(idx3);
      im3 = vis_im(idx3)*sgn3;

      re23 = zmult_re(re2, im2, re3, im3);
      im23 = zmult_im(re2, im2, re3, im3);

      re = zmult_re(re1, im1, re23, im23) - db.t3re;
      im = zmult_im(re1, im1, re23, im23) - db.t3im;

      temp_re = db.wrr*re + db.wri*im;
      temp_im = db.wri*re + db.wii*im;
      err3 = sum(temp_re*re) + sum(temp_im*im);
      ndata3 = 2*numberof(re);
      re = im = nil;

      temp_re = temp_re + temp_re; /* dERR3 / dRe(Z) */
      temp_im = temp_im + temp_im; /* dERR3 / dIm(Z) */

      grd(1, idx1) +=  temp_re*re23 + temp_im*im23;
      grd(2, idx1) += (temp_im*re23 - temp_re*im23)*sgn1;
      re23 = im23 = nil;

      re13 = zmult_re(re1, im1, re3, im3);
      im13 = zmult_im(re1, im1, re3, im3);
      grd(1, idx2) +=  temp_re*re13 + temp_im*im13;
      grd(2, idx2) += (temp_im*re13 - temp_re*im13)*sgn2;
      re13 = im13 = nil;

      re12 = zmult_re(re1, im1, re2, im2);
      im12 = zmult_im(re1, im1, re2, im2);
      grd(1, idx3) +=  temp_re*re12 + temp_im*im12;
      grd(2, idx3) += (temp_im*re12 - temp_re*im12)*sgn3;
      re12 = im12 = nil;
    }

  }

  /* Convert gradient with to real and imaginary parts into gradient with
     respect to pixel values. */
  grd = master.xform(grd, 1);
  h_set, master, ndata = (ndata1 + ndata2 + ndata3);
  return (err1 + err2 + err3);
}

func _mira_zmult_re(re1, im1, re2, im2) { return re1*re2 - im1*im2; }
func _mira_zmult_im(re1, im1, re2, im2) { return re1*im2 + im1*re2; }

/*---------------------------------------------------------------------------*/
/* IMAGE RECONSTRUCTION */

local _mira_window;
func mira_window_new
{
  extern _mira_window;
  _mira_window = -1;
}

func _mira_solve_viewer(x, extra)
{
  extern _mira_window;
  flags = extra.view;
  if ((flags & (1|2|4|8)) == 0) {
    /* Nothing to display. */
    return;
  }

  /* Setup for graphics. */
  if (is_void(_mira_window) || ! window_exists(_mira_window)) {
    vp = [[5,45,55,95],[55,95,55,95],[5,65,5,45],[75,95,5,25]]*1e-2;
    xwindow, dpi=100, viewport=vp, units=1, size=6, xopt=, yopt=,
      frame = [0,0,1,1];
    _mira_window = current_window();
    palette, "earth.gp";
  }
  if (is_void(x)) {
    /* Nothing to display. */
    return;
  }

  /* Re-normalization. */
  normalization = extra.normalization;
  if (normalization && (xsum = sum(x)) > 0) {
    xscl = normalization/double(xsum);
    if (xscl != 1) {
      x *= xscl;
    }
  }

  /* Make graphics. */
  window, _mira_window;
  fma;
  if ((flags & 2) != 0) {
    plsys, 3;
    mira_plot, extra.master, x, "vis2";
  }
  if ((flags & 4) != 0) {
    plsys, 2;
    mira_plot, extra.master, x, "vis2-2D";
 }
  if ((flags & 8) != 0) {
    plsys, 4;
    mira_plot, extra.master, x, "vis3";
  }
  
  /* Display current solution. */
  if ((flags & 1) != 0) {
    plsys, 1;
    // FIXME: CMIN/CMAX
    cmin = 0.0;
    cmax = max(x);
    //if (min(x) != cmax) {
    //  cmax = max(x(where(x != cmax)));
    //}
    pli, x, cmin=cmin, cmax=cmax;
    //mira_plot_image, /*log(1e-6*max(x) + x)*/ x, extra.master,
    //  cmin=cmin, cmax=cmax, clear=1, zformat="%+.2e";
    //title = extra.title;
    //if (structof(title) == string) pltitle, title;
  }

  pause, 1;
}

func _mira_solve_cost2(self, x, &g)
{
  return _mira_solve_cost(x, g, self);
}

func _mira_solve_cost(x, &grd, extra)
{
  /* Evaluate the function and the gradient at X. */

  /* Re-normalization. */
  normalization = extra.normalization;
  if (normalization && (xsum = sum(x)) > 0) {
    xscl = normalization/double(xsum);
    if (xscl != 1) {
      x *= xscl;
    }
  } else {
    normalization = 0n;
  }
  if (extra.zap_data) {
    data_err = 0.0;
  } else {
    data_err = mira_data_penalty(extra.master, x, grd);
  }
  if (extra.mu) {
    rgl_update, extra.regul, x;
    regul_err = rgl_get_penalty(extra.regul, x);
    regul_grd = rgl_get_gradient(extra.regul, x);
    if (extra.zap_data) {
      eq_nocopy, grd, regul_grd;
    } else {
      grd += regul_grd;
    }
  } else {
    regul_err = 0.0;
  }
  if (normalization) {
    /* Fix the gradient. */
    grd = xscl*grd - sum(grd*x)/xsum;
  }
  total_err = data_err + regul_err;
  if (! h_has(extra, total_err=) || extra.total_err > total_err) {
    h_set, extra, total_err = total_err, data_err = data_err,
      regul_err = regul_err, grd = grd;
  }
  return total_err;
}

func mira_select(this, select)
/* DOCUMENT other = mira_select(this, cutoff);
 *     -or- other = mira_select(this, select);
 *   Returns a clone of MiRA instance THIS with only data below frequency
 *   CUTOFF (a real scalar) or with only data at spatial frequencies for
 *   which SELECT (a vector of same size as THIS.u or THIS.v) is true.  If
 *   all data are selected, THIS is returned rather than a clone of it.
 *
 * SEE ALSO: mira_new, mira_solve.
 */
{
  if (is_void(select)) return this;
  s = structof(select);
  if (is_scalar(select) && (s == double || s == float) && select > 0.0) {
    /* Get rid of data outside cutoff frequency. */
    select = (this.freq_len <= select);
  } else if (is_vector(select) && numberof(select) == numberof(this.u)) {
    select = !(!select);
  } else {
    error, "expecting a scalar cutoff frequency or a boolean vector";
  }
  other = h_new();
  for (key = h_first(this); key; key = h_next(this, key)) {
    if (key == "vis") {
      db = h_get(this, key);
      j = where(select(db.idx));
      if (is_array(j)) {
        h_set, other, key, h_new(idx=db.idx(j),
                                 sgn=db.sgn(j),
                                 u=db.u(j),
                                 v=db.v(j),
                                 w=db.w(j),
                                 amp=db.amp(j),
                                 amperr=db.amperr(j),
                                 phi=db.phi(j),
                                 phierr=db.phierr(j),
                                 re=db.re(j),
                                 im=db.im(j),
                                 wii=db.wii(j),
                                 wri=db.wri(j),
                                 wrr=db.wrr(j));
      }
    } else if (key == "vis2") {
      db = h_get(this, key);
      j = where(select(db.idx));
      if (is_array(j)) {
        h_set, other, key, h_new(idx=db.idx(j),
                                 sgn=db.sgn(j),
                                 u=db.u(j),
                                 v=db.v(j),
                                 w=db.w(j),
                                 vis2data=db.vis2data(j),
                                 vis2err=db.vis2err(j),
                                 vis2wght=db.vis2wght(j));
      }
    } else if (key == "vis3") {
      db = h_get(this, key);
      j = where(select(db.idx1)&select(db.idx2)&select(db.idx3));
      if (is_array(j)) {
        h_set, other, key, h_new(idx1=db.idx1(j),
                                 idx2=db.idx2(j),
                                 idx3=db.idx3(j),
                                 sgn1=db.sgn1(j),
                                 sgn2=db.sgn2(j),
                                 sgn3=db.sgn3(j),
                                 u1=db.u1(j),
                                 u2=db.u2(j),
                                 v1=db.v1(j),
                                 v2=db.v2(j),
                                 w=db.w(j),
                                 /* FIXME: cl_weight=db.cl_weight(j),*/
                                 t3amp=db.t3amp(j),
                                 t3amperr=db.t3amperr(j),
                                 t3im=db.t3im(j),
                                 t3phi=db.t3phi(j),
                                 t3phierr=db.t3phierr(j),
                                 t3re=db.t3re(j),
                                 wii=db.wii(j),
                                 wri=db.wri(j),
                                 wrr=db.wrr(j));
      }
    } else {
      h_set, other, key, h_get(this, key);
    }
  }
  return other;
}

func mira_new_functor(args)
/* DOCUMENT obj = mira_new_functor(fn, key1=val1, keu2=val2, ...);

      The function mira_new_functor() creates a functor which calls
      function FN with itself prepended to its argument list:

         obj(arg1, arg2, ...)

      is the same as:

         fn(obj, arg1, arg2, ...)
      
      Argument FN can be a name or any object callable as a function
      (including another functor).  Any given keywords will be stored into the
      returned object (which can be used as a hash-table) and are accessible
      by OBJ.KEY1, OBJ.KEY2, etc.

      The object is built as:

         obj = h_new(key1=val1, keu2=val2, ...);
         h_evaluator, obj, fn;
      

   SEE ALSO: h_new, h_evaluator, wrap_args.
 */
{
  local f;
  nargs = args(*);
  if (nargs != 1) {
    error, "expecting 1 positional argument";
  }
  eq_nocopy, fn, args(1);
  if (is_void(fn)) {
    error, "invalid function argument";
  }
  obj = h_new();
  keys = args(-);
  nkeys = numberof(keys);
  for (k = 1; k <= nkeys; ++k) {
    key = keys(k);
    //if (h_has(obj, key)) {
    //  error, ("keyword \"" + key + "\" is reserved");
    //}
    h_set, obj, key, args(key);
  }
  h_evaluator, obj, fn;
  return obj;
}
wrap_args, mira_new_functor;

local mira_solve; /* only to provide the documentation */
/* DOCUMENT img = mira_solve(data, key1=value1, key2=value2, ...);
 *      or: img = mira_solve(data, img_init, key1=value1, key2=value2, ...);
 *      or: img = mira_solve(data, img_init, penalty, key1=value1, ...);
 *
 *   Builds an image from the interferometric data stored into instance DATA
 *   (see mira_new) which must have been properly initialized for
 *   reconstruction (see mira_config).  The reconstruction is done by
 *   minimizing an objective function which is the sum of the likelihood
 *   penalty (which enforces agreement of the model image with the data) and a
 *   regularization penalty (which enforces some priors about the image).
 *
 *   Optional argument IMG_INIT is the initial image.  Most arguments are
 *   provided by keywords (see below). At least, you want to specify the
 *   regularization method and level by the keywords REGUL and MU.  Usually,
 *   you also want to impose a positivity constraint with XMIN=0.0 (use a
 *   small but strictly positive value if regularization such as entropy with
 *   a logarithm is used) and a normalization constraint with
 *   NORMALIZATION=1.0 (which is suitable for OI-FITS data, otherwise the
 *   actual value should be equal to the total flux in the image).
 *
 *   Optional argument PENALTY is a simple symbol name to store the values of
 *   the penalty terms at the final solution as a vector of 5 values: [NDATA,
 *   DATA_COST, REGUL_WGHT, REGUL_COST, TOTAL_COST, GPNORM] where NDATA is the
 *   number of measurements, DATA_COST is the data penalty per datum,
 *   REGUL_WGHT and REGUL_COST are the regularization weight and penalty,
 *   TOTAL_COST = NDATA*DATA_COST + REGUL_WGHT*REGUL_COST and GPNORM is the
 *   Euclidean norm of the (projected) gradient.
 *
 *
 * KEYWORDS
 *   XMIN - minimum allowed value in the image; can be a scalar or a
 *          pixel-wise value; there is no limit if this option is not set.
 *   XMAX - maximum allowed value in the image; can be a scalar or a
 *          pixel-wise value; there is no limit if this option is not set.
 *   NORMALIZATION - value of the sum of the image pixels to impose a
 *          normalization constraint.
 *   REGUL - regularization method (see rgl_new).
 *   MU   - regularization level; the higher is MU the more the solution is
 *          influenced by the priors set by the regularization.
 *   MAXITER - maximum number of iterations, unlimited if not set.
 *   MAXEVAL - maximum number of evaluations of the objective function,
 *          unlimited if not set.
 *   OUTPUT - output stream/file for iteration informations; default is
 *          standard text output (see keyword VERB).
 *   VERB - verbose level: informations get printed out every VERB iterations
 *          and at convergence.
 *   VIEW - Bitwise mask to specify which graphics to show.  If unset, all
 *          graphics are shown.  Otherwise, if bit 0x01 is set, the current
 *          image is drawn into 1st sub-window; if bit 0x02 is set, a radial
 *          plot of the visibility amplitudes is drawn into 2nd sub-window; if
 *          bit 0x04 is set, a 2-D plot of the complex visibilities is drawn
 *          into 4th sub-window; if bit 0x03 is set, an histogram of the phase
 *          closure resiudals is drawn into 3rd sub-window.
 *   SELECT - if set, select data to fit; the value can be the cutoff
 *          frequency or a boolean vector set true for frequencies to keep
 *          (see mira_select).
 *   ZAP_AMPLITUDE - if true, the image reconstruction is performed without
 *          any visibility amplitude data (FIXME: not yet tested).
 *   ZAP_PHASE - if true, the image reconstruction is performed without any
 *          phase data.
 *   ZAP_DATA - if true, the image reconstruction is performed without any
 *          data, that is with only the regularization.  Useful to get the
 *          default solution set by the priors.
 *   HANIFF - use Haniff's method to account for phase modulo 2-PI
 *          uncertainty; the default is to use phasors which yields better
 *          convergence properties.
 *   MEM  - control the memory usage by the optimizer; the value is the
 *          number of corrections and gradient differences memorized by the
 *          variable metric algorithm; by default, MEM=7 (see op_mnb).
 *   FTOL - relative function tolerance for the stopping criterion of
 *          the optimizer; default value is: FTOL = 1e-12 (see op_mnb).
 *   GTOL - gradient tolerance for the stopping criterion of the
 *          optimizer; default value is: GTOL = 0.0 (see op_mnb).
 *   SFTOL, SGTOL, SXTOL - control the stopping criterion of the
 *          line-search method in the optimizer (see op_mnb).
 *
 * SEE ALSO:
 *   mira_new, mira_config.
 */
func mira_solve(master, x, &penalty, reset=, fix=,
                normalization=,
                haniff=, zap_phase=, zap_data=, zap_amplitude=,
                cubic=,
                view=, title=,
                cmin=, cmax=,
                select=,
                regul=, mu=,
                data_cost=, data_hyper=,
                png_format=,
                colortable=, movie_file=, movie_fps=,
                /* options for OptimPack */
                xmin=, xmax=,
                method=, mem=, verb=, factor=, wait=,
                maxiter=, maxeval=, output=,
                ftol=, gtol=, sftol=, sgtol=, sxtol=,
                gpnormconv=, get_cost=)
{
  if (is_void(view)) view = -1; /* default is to show every graphics */

  /* Set default values for optimizer. */
  if (is_void(ftol)) ftol =  1e-12;
  if (is_void(gtol)) gtol = 0.0;
  if (is_void(mem)) mem = 7;

  /* Update internal cache. */
  if (master.update_pending) {
    mira_update, master;
  }

  /* Data selection. */
  h_set, master, haniff=(haniff ? 1n : 0n),
    zap_amplitude=(zap_amplitude ? 1n : 0n),
    zap_phase=(zap_phase ? 1n : 0n);
  if (! is_void(select)) {
    master = mira_select(master, select);
  }

  /* Default solution. */
  dim = master.dim;
  if (is_void(x)) {
    //x = array(1.0/(dim^2), dim, dim);
    x = random(dim,dim);
  } else if (numberof(x) != dim*dim) {
    x = mira_rescale(x, dim, dim, cubic=0);
  }
  dims = dimsof(x);

  /* Get cost function for the data. */
  if (is_void(data_cost)) {
    data_cost = cost_l2;
  } else if (is_scalar(data_cost) && is_string(data_cost)) {
    data_cost = strcase(1, data_cost);
    if (data_cost == "L2") {
      data_cost = cost_l2;
    } else if (data_cost == "L2-L0") {
      data_cost = cost_l2l0;
    } else if (data_cost == "L2-L1") {
      data_cost = cost_l2l0;
    } else {
      error, "bad cost function name for data";
    }
  } else if (! is_func(data_cost)) {
    error, "data cost must be a name or a function";
  }

  /* Get hyper-parameters for data. */
  if (is_void(data_hyper)) {
    data_hyper = 1.0;
  } else {
    if (data_hyper(1) != 1) {
      error, "first element of DATA_HYPER must be equal to 1";
    }
  }

  /* Setup for regularization. */
  if (is_void(regul)) {
    mu = 0.0;
  } else {
    temp = rgl_get_global_weight(regul);
    if (is_void(mu)) {
      mu = temp;
    } else if (is_scalar(mu) && (is_real(mu) || is_integer(mu)) && mu >= 0.0) {
      mu = double(mu);
      if (mu != temp) {
        rgl_set_global_weight, regul, mu;
      }
    } else {
      error, "global regularization weight MU must be a non negative scalar";
    }
  }

  extra = h_new(normalization=normalization,
                zap_data=zap_data,
                mu=mu,
                regul=regul,
                master=master,
                view=view,
                wait=wait,
                title=title);
  if (get_cost) {
    h_evaluator, extra, "_mira_solve_cost2";
    return extra;
  }
  cost = _mira_solve_cost;
  
  viewer = _mira_solve_viewer;
  printer = _mira_solve_printer;
  if (verb) viewer, , extra;

  // FIXME: fix bug in op_mnb:
  if (structof(x) != double || (is_void(xmin) && is_void(xmax))) {
    x = double(x); // force copy and conversion
  }

  x = mira_mnb(cost, x, penalty, extra=extra,
               xmin=xmin, xmax=xmax, method=method, mem=mem,
               verb=verb, viewer=viewer, printer=printer,
               maxiter=maxiter, maxeval=maxeval, output=,
               frtol=ftol, fatol=0.0, sftol=sftol, sgtol=sftol, sxtol=sftol,
               gpnormconv=gpnormconv);

  /* Re-normalization. */
  if (normalization && (xsum = sum(x)) > 0) {
    xscl = normalization/double(xsum);
    if (xscl != 1) {
      x *= xscl;
    }
  }
  gpnorm =  mira_projected_gradient_norm(x, extra.grd, xmin=xmin, xmax=xmax);
  penalty = [(master.ndata ? master.ndata : 0),
             (master.ndata ? extra.data_err/master.ndata : 0.0),
             (extra.mu ? extra.mu : 0.0),
             (extra.mu ? extra.regul_err/extra.mu : 0.0),
             extra.total_err,
             gpnorm];
  return x;
}

func _mira_solve_printer(output, iter, eval, cpu, fx, gnorm, steplen, x, extra)
{
  if (eval == 1) {
    write, output, format="# %s\n# %s\n",
      "ITER  EVAL   CPU (ms)        FUNC               <FDATA>     FPRIOR    GNORM   STEPLEN",
      "-------------------------------------------------------------------------------------";
  }
  avg_data_err = (extra.master.ndata > 0 ?
                  extra.data_err/extra.master.ndata : 0.0);
  write, output,
    format=" %5d %5d %10.3f  %+-24.15e%-11.3e%-11.3e%-9.1e%-9.1e\n",
    iter, eval, cpu, fx, avg_data_err,
    extra.regul_err, gnorm, step;
  if (extra.wait) {
    pause, extra.wait;
  }
}

func mira_mnb(f, x, &fx, &gx, fmin=,
              xmin=, xmax=, method=, mem=, verb=, quiet=,
              viewer=, printer=, extra=,
              maxiter=, maxeval=, output=,
              frtol=, fatol=, sftol=, sgtol=, sxtol=,
              gpnormconv=)
/* DOCUMENT mira_mnb(f, x)
 *     -or- mira_mnb(f, x, fout, gout)
 *
 *   Returns a minimum of a multivariate function by an iterative
 *   minimization algorithm (conjugate gradient or limited memory variable
 *   metric) possibly with simple bound constraints on the parameters.
 *   Arguments are:
 *
 *     F - User defined function to optimize.
 *         The prototype of F is:
 *           func F(x, &gx) {
 *             fx = ....; // compute function value at X
 *             gx = ....; // store gradient of F in GX
 *             return fx; // return F(X)
 *           }
 *
 *     X - Starting solution (a floating point array).
 *
 *     FOUT - Optional output variable to store the value of F at the
 *         minimum.
 *
 *     GOUT - optional output variable to store the value of the gradient
 *         of F at the minimum.
 *
 *   If the multivariate function has more than one minimum, which minimum
 *   is returned is undefined (although it depends on the starting
 *   parameters X).
 *
 *   In case of early termination, the best solution found so far is
 *   returned.
 *
 *
 * KEYWORDS
 *
 *   EXTRA - Supplemental argument for F; if non-nil, F is called as
 *       F(X,GX,EXTRA) so its prototype must be: func F(x, &gx, extra).
 *
 *   XMIN, XMAX  - Lower/upper bounds for  X.  Must be  conformable with X.
 *       For instance with XMIN=0, the non-negative solution will be
 *       returned.
 *
 *   METHOD - Scalar integer which  defines the optimization method to use.
 *       Conjugate  gradient   algorithm  is  used  if  one   of  the  bits
 *       OP_FLAG_POLAK_RIBIERE,         OP_FLAG_FLETCHER_REEVES,         or
 *       OP_FLAG_HESTENES_STIEFEL  is  set;  otherwise,  a  limited  memory
 *       variable  metric algorithm  (VMLM-B) is  used.  If  METHOD  is not
 *       specified and  if MEM=0, a conjugate gradient  search is attempted
 *       with flags: (OP_FLAG_UPDATE_WITH_GP |
 *                    OP_FLAG_SHANNO_PHUA    |
 *                    OP_FLAG_MORE_THUENTE   |
 *                    OP_FLAG_POLAK_RIBIERE  |
 *                    OP_FLAG_POWELL_RESTART)
 *       otherwise VMLM-B is used with flags: (OP_FLAG_UPDATE_WITH_GP |
 *                                             OP_FLAG_SHANNO_PHUA    |
 *                                             OP_FLAG_MORE_THUENTE).
 *       See documentation  of op_get_flags to  figure out the  allowed bit
 *       flags and their meaning.
 *
 *   MEM - Number of previous directions used in variable metric limited
 *       memory method (default min(7, numberof(X))).
 *
 *   MAXITER - Maximum number of iterations (default: no limits).
 *
 *   MAXEVAL - Maximum number of function evaluations (default: no limits).
 *
 *   FTOL - Relative function change tolerance for convergence (default:
 *       1.5e-8).
 *
 *   GTOL - Gradient tolerance for convergence (default: 3.7e-11).
 *
 *   VERB - Verbose mode?  If non-nil and non-zero, print out information
 *       every VERB iterations and for the final one.
 *
 *   QUIET - If true and not in verbose mode, do not print warning nor
 *       convergence error messages.
 *
 *   OUPTPUT - Output for verbose mode.  For instance, text file stream
 *       opened for writing.
 *
 *   VIEWER - User defined subroutine to call every VERB iterations (see
 *       keyword VERB above)to display the solution X.  The subroutine will
 *       be called as:
 *          viewer, x, extra;
 *       where X is the current solutio and EXTRA is the value of keyword
 *       EXTRA (which to see).  If the viewer uses Yorick graphics
 *       window(s) it may call "pause, 1;" before returning to make sure
 *       that graphics get correctly updated.
 *
 *   PRINTER - User defined subroutine to call every VERB iterations (see
 *       keyword VERB above) to printout iteration information.
 *       The subroutine will be called as:
 *          printer, output, iter, eval, cpu, fx, gnorm, steplen, x, extra;
 *       where OUTPUT is the value of keyword OUTPUT (which to see), ITER
 *       is the number of iterations, EVAL is the number of function
 *       evaluations, CPU is the elapsed CPU time in seconds, FX is the
 *       function value at X, GNORM is the Euclidean norm of the gradient
 *       at X, STEPLEN is the length of the step along the search
 *       direction, X is the current solution and EXTRA is the value of
 *       keyword EXTRA (which to see).
 *
 *   SFTOL, SGTOL, SXTOL, SXBIG - Line   search   tolerance  and  safeguard
 *      parameters (see op_csrch).
 *
 * SEE ALSO: op_get_flags, op_csrch,
 *           op_cgmnb_setup, op_cgmnb_next,
 *           op_vmlmb_setup, op_vmlmb_next.
 */
{
  local result, gx;

  /* Get function. */
  if (! is_func(f)) {
    error, "expecting a function for argument F";
  }
  use_extra = (! is_void(extra));

  /* Starting parameters. */
  if ((s = structof(x)) != double && s != float && s != long &&
      s != int && s != short && s != char) {
    error, "expecting a numerical array for initial parameters X";
  }
  n = numberof(x);
  dims = dimsof(x);

  /* Bounds on parameters. */
  bounds = 0;
  if (! is_void(xmin)) {
    if (is_void((t = dimsof(x, xmin))) || t(1) != dims(1)
        || anyof(t != dims)) {
      error, "bad dimensions for lower bound XMIN";
    }
    if ((convert = (s = structof(xmin)) != double) && s != float &&
        s != long && s != int && s != short && s != char) {
      error, "bad data type for lower bound XMIN";
    }
    if (convert || (t = dimsof(xmin))(1) != dims(1) || anyof(t != dims)) {
      xmin += array(double, dims);
    }
    bounds |= 1;
  }
  if (! is_void(xmax)) {
    if (is_void((t = dimsof(x, xmax))) || t(1) != dims(1)
        || anyof(t != dims)) {
      error, "bad dimensions for lower bound XMAX";
    }
    if ((convert = (s = structof(xmax)) != double) && s != float &&
        s != long && s != int && s != short && s != char) {
      error, "bad data type for lower bound XMAX";
    }
    if (convert || (t = dimsof(xmax))(1) != dims(1) || anyof(t != dims)) {
      xmax += array(double, dims);
    }
    bounds |= 2;
  }

  /* Output stream. */
  if (! is_void(output)) {
    if (structof(output) == string) {
      output = open(output, "a");
    } else if (typeof(output) != "text_stream") {
      error, "bad value for keyword OUTPUT";
    }
  }

  /* Maximum number of iterations and function evaluations. */
  check_iter = (! is_void(maxiter));
  check_eval = (! is_void(maxeval));

  /* Viewer and printer subroutines. */
  if (is_void(printer)) {
    use_printer = 0n;
  } else if (mira_is_func(printer)) {
    use_printer = 1n;
  } else {
    error, "bad value for keyword PRINTER";
  }
  if (is_void(viewer)) {
    use_viewer = 0n;
  } else if (mira_is_func(viewer)) {
    use_viewer = 1n;
  } else {
    error, "bad value for keyword VIEWER";
  }


  /* Choose minimization method. */
  //if (is_void(frtol)) frtol = 1e-10;
  //if (is_void(fatol)) fatol = 1e-10;
  if (! method) {
    /* Variable metric. */
    if (is_void(mem)) mem = min(n, 7);
    if (is_void(fmin)) fmin = 0.0;
    method = 0;
    method_name = swrite(format="Limited Memory BFGS (VMLM with MEM=%d)",
                         mem);
    ws = op_vmlmb_setup(n, mem, /*fmin=fmin,*/
                        fatol=fatol, frtol=frtol,
                        sftol=sftol, sgtol=sgtol, sxtol=sxtol);
  } else if (method < 0) {
    if (is_void(mem)) mem = min(n, 7);
    method_name = swrite(format="Limited Memory BFGS (LBFGS with MEM=%d)",
                         mem);
    ws = op_lbfgs_setup(n, mem);
  } else if (method >= 1 && method <= 15) {
    /* Conjugate gradient. */
    mem = 2;
    error, "conjugate-gradient method not yet implemented";
    method_name = swrite(format="Conjugate Gradient (%s)",
                         ["Fletcher-Reeves", "Polak-Ribiere",
                          "Polak-Ribiere with non-negative BETA"](method&3));
    ws = optim_cgmn_setup(method, fmin=fmin, fatol=fatol, frtol=frtol);
  } else {
    error, "bad METHOD";
  }
  step = 0.0;
  task = 1;
  eval = iter = 0;
  stop = 0n;
  if (verb) {
    elapsed = array(double, 3);
    timer, elapsed;
    cpu_start = elapsed(1);
  }
  for (;;) {
    local gx; /* to store the gradient */
    if (task == 1) {
      /* Evaluate function and gradient. */
      if (bounds) {
        if (bounds & 1) {
          x = max(x, xmin);
        }
        if (bounds & 2) {
          x = min(x, xmax);
        }
      }
      fx = (use_extra ? f(x, gx, extra) : f(x, gx));
      ++eval;
      if (bounds) {
        /* Figure out the set of free parameters:
         *   ACTIVE(i) = 0 if X(i) has a lower bound XMIN(i)
         *                 and X(i) = XMIN(i) and GX(i) >= 0
         *               0 if X(i) value has an upper bound XMAX(i)
         *                 and X(i) = XMAX(i) and GX(i) <= 0
         *               1 (or any non-zero value) otherwise
         */
        if (bounds == 1) {
          active = ((x > xmin) | (gx < 0.0));
        } else if (bounds == 2) {
          active = ((x < xmax) | (gx > 0.0));
        } else {
          active = (((x > xmin) | (gx < 0.0)) & ((x < xmax) | (gx > 0.0)));
        }
      }
    }

    /* Check for convergence. */
    if (task != 1 || eval == 1) {
      gnorm =  mira_projected_gradient_norm(x, gx, xmin=xmin, xmax=xmax);
      if (task > 2) {
        stop = 1n;
        msg = op_vmlmb_msg(ws);
      } else if (! is_void(gpnormconv) && gnorm <= gpnormconv) {
        stop = 1n;
        msg = "GP norm <= pre-set level";
      } else if (check_iter && iter > maxiter) {
        stop = 1n;
        msg = swrite(format="warning: too many iterations (%d)\n", iter);
      } else if (check_eval && eval > maxeval) {
        stop = 1n;
        msg = swrite(format="warning: too many function evaluations (%d)\n",
                     eval);
      }
      if (verb) {
        if (eval == 1 && ! use_printer) {
          write, output, format="# Method %d (MEM=%d): %s\n#\n",
            method, mem, method_name;
          write, output, format="# %s\n# %s\n",
            "ITER  EVAL   CPU (ms)        FUNC               GNORM   STEPLEN",
            "---------------------------------------------------------------";
        }
        if (stop || ! (iter % verb)) {
          timer, elapsed;
          cpu = 1e3*(elapsed(1) - cpu_start);
          // FIXME: this is a hack...
          if (use_printer) {
            printer, output, iter, eval, cpu, fx, gnorm, steplen, x, extra;
          } else {
            write, output, format=" %5d %5d %10.3f  %+-24.15e%-9.1e%-9.1e\n",
              iter, eval, cpu, fx, gnorm, step;
          }
          if (use_viewer) {
            viewer, x, extra;
          }
        }
      }
      if (stop) {
        if (msg && (verb || (task != 3 && ! quiet))) {
          write, output, format="# %s\n", strtrim(msg, 2, blank=" \t\v\n\r");
        }
        return x;
      }
    }

    /* Call optimizer. */
    if (! method) {
      task = op_vmlmb_next(x, fx, gx, ws, active);
      iter = (*ws(2))(7);
      step = (*ws(3))(22);
    } else if (method < 0) {
      task = op_lbfgs_next(x, fx, gx, ws);
      if (task == 2 || task == 3) ++iter;
      step = -1.0;
    }
  }
}

func mira_projected_gradient_norm(x, gx, xmin=, xmax=)
{
  // FIXME: normalization not take into account
  local gp;
  if (is_void(xmin)) {
    eq_nocopy, gp, gx;
  } else {
    gp = gx*((gx < 0.0)|(x > xmin));
  }
  if (! is_void(xmax)) {
    gp *= ((gx > 0.0)|(x < xmax));
  }
  return sqrt(sum(gp*gp));
}


/*---------------------------------------------------------------------------*/
/* PSEUDO-OBJECT MANAGEMENT */

local mira_get;
local mira_get_w, mira_get_x, mira_get_y;
/* DOCUMENT mira_get_w(this) - returns wavelength(s)
 *     -or- mira_get_x(this) - returns sky X-coordinates
 *     -or- mira_get_y(this) - returns sky Y-coordinates
 *
 *   These functions can be used to query internals of MiRA master
 *   object THIS.
 *
 * SEE ALSO:
 */
func mira_get_w(this) { return this.w; }
func mira_get_x(this)
{
  if (this.update_pending) mira_update, this;
  return this.pixelsize*(indgen(this.dim) - 0.5*(this.dim + 1));
}
func mira_get_y(this)
{
  if (this.update_pending) mira_update, this;
  return this.pixelsize*(indgen(this.dim) - 0.5*(this.dim + 1));
}

func mira_get_ndata(this)
/* DOCUMENT mira_get_ndata(this);
 *   Get the number of valid measurements which have been used during last
 *   fit or image reconstruction (e.g. by the mira_solve routine) involving
 *   MiRA opaque handle THIS.
 *
 * SEE ALSO: mira_solve.
 */
{
  return ((ndata = this.ndata) ? ndata : 0);
}

func mira_get_dim(this) { return this.dim; }
/* DOCUMENT mira_get_dim(this);
 *   Get the number of pixels per side for the model image assumed by MiRA
 *   opaque handle THIS.
 *
 * SEE ALSO: mira_config, mira_get_fov, mira_get_pixelsize.
 */

func mira_get_fov(this) { return this.fov; }
/* DOCUMENT mira_get_fov(this);
 *   Get the width of the field of view (in radians) for the model image
 *   assumed by MiRA opaque handle THIS.
 *
 * SEE ALSO: mira_config, mira_get_dim, mira_get_pixelsize.
 */

func mira_get_pixelsize(this) { return this.pixelsize; }
/* DOCUMENT mira_get_pixelsize(this);
 *   Get the pixel size (in radians) for the model image assumed by MiRA
 *   opaque handle THIS.
 *
 * SEE ALSO: mira_config, mira_get_dim, mira_get_pixelsize.
 */

/*---------------------------------------------------------------------------*/
/* PRIVATE UTILITIES */

func _mira_warn(s) { write, format="WARNING: %s\n", s; }
func _mira_info(s) { write, format="INFO: %s\n", s; }
func _mira_debug(s) { write, format="DEBUG: %s\n", s; }

func _mira_add(master, key, value, index)
{
  tmp = h_pop(master, key);
  if (is_void(index)) {
    h_set, master, key, (is_void(tmp) ? value : tmp + value);
  } else {
    tmp(index) += value;
    h_set, master, key, tmp;
  }
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func mira_plot_image(img, this, clear=, cmin=, cmax=, zformat=, keeplimits=,
                     normalize=, pixelsize=, pixelunits=)
/* DOCUMENT mira_plot_image, img;
 *     -or- mira_plot_image, img, this;
 *
 *   Plot image IMG in current graphics window.  The pixelsize is taken
 *   from MiRA instance THIS if it is provided.
 *
 *   Unless keyword KEEPLIMITS is true, the axis orientations are set to
 *   match conversions in astronomy (see mira_fix_image_axis).
 *
 *   Keyword CLEAR can be used to call fma command (which see): if CLEAR >
 *   0, fma is called prior to drawing the image; if CLEAR < 0, fma is
 *   called after drawing the image (useful in animate mode); if CLEAR = 0
 *   or undefined, then fma is not called.
 *
 *   Keyword ZFORMAT can be used to specify the format for the color bar
 *   labels (see mira_color_bar).
 *
 *   Keyword PIXELSIZE and PIXELUNITS can be used to specify the size of
 *   the pixel in given angular units and the name of these units.  Both
 *   must be specified to be taken into account.  If MiRA instance THIS is
 *   provided, these keywords are ignored.
 *
 *   Keywords CMIN and CMAX can be used to specify cut levels for the
 *   display (see pli).  If keyword NORMALIZE (see below) is true, CMIN and
 *   CMAX are in units of normalized intensity.

 *   If keyword NORMALIZE is true, the flux is normalized (divided by
 *   PIXELSIZE^2).
 *
 *
 * SEE ALSO:
 *   pli, fma, animate, mira_color_bar, mira_fix_image_axis, xytitles.
 */
{
  dims = dimsof(img);
  if (is_void(dims) || dims(1) != 2) {
    error, "expecting a 2-D image";
  }
  width = dims(2);
  height = dims(3);

  if (is_hash(this)) {
    pixelsize = this.pixelsize/MIRA_MILLIARCSECOND;
    pixelunits = "milliarcseconds";
  } else {
    if (is_void(pixelsize) || is_void(pixelunits)) {
      pixelsize = 1.0;
      pixelunits = "pixels";
    }
  }
  if (is_void(normalize)) {
    scl = 1.0;
  } else {
    scl = 1.0/pixelsize^2;
    if (scl != 1) {
      img *= scl;
    }
  }

  if (is_void(cmin)) cmin = min(img);
  if (is_void(cmax)) cmax = max(img);

  x0 = -(x1 = 0.5*pixelsize*width);
  y0 = -(y1 = 0.5*pixelsize*height);
  if (clear && clear > 0) {
    fma;
  }
  pli, img, x0, y0, x1, y1, cmin=cmin, cmax=cmax;
  if (! keeplimits) mira_fix_image_axis;
  local red, green, blue;
  palette, red, green, blue, query=1;
  ncolors = numberof(red);
  levs = span(cmin, cmax, ncolors + 1);
  colors = indgen(ncolors);
  mira_color_bar, cmin=cmin, cmax=cmax, vert=1, nlabs=11, format=zformat;
  xytitles, "relative !a ("+pixelunits+")", "relative !d ("+pixelunits+")";
  if (clear && clear < 0) {
    fma;
  }
}

func mira_fix_image_axis
/* DOCUMENT mira_fix_image_axis;
     Fix orientation of horizontal axis (right ascension, RA) and vertical
     axis (declination, DEC) in current window to match the conventions in
     astronomy to display RA positive toward East (that is left of graphics)
     and DEC positive toward North (that is top of graphics).  In other words,
     the horizontal axis is reversed with respect to mathematical convention.
  
   SEE ALSO: limits, mira_plot_image, mira_plot_baselines.
 */
{
  lm = limits();
  x0 = lm(1);
  x1 = lm(2);
  y0 = lm(3);
  y1 = lm(4);
  if (x0 < x1 || y0 > y1) {
    limits, max(x0, x1), min(x0, x1), min(y0, y1), max(y0, y1);
  }
}

func mira_color_bar(z, cmin=, cmax=, vert=, nlabs=, adjust=,
                    color=, font=, height=, opaque=, orient=,
                    width=, ticklen=, thickness=, vport=, format=)
/* DOCUMENT mira_color_bar, z;
 *     -or- mira_color_bar, cmin=CMIN, cmax=CMAX;
 *
 *   Draw a color bar below the current coordinate system the colors and
 *   the associated label values are from min(Z) to max(Z) -- alternatively
 *   keywords CMIN and CMAX can be specified.  With the VERT=1 keyword the
 *   color bar appears to the left of the current coordinate system (vert=0
 *   is the default).
 *
 *   Keyword NLABS can be used to choose the number of displayed labels; by
 *   default, NLABS=11 which correspond to a label every 10% of the
 *   dynamic; use NLABS=0 to suppress all labels.  The format of the labels
 *   can be specified with keyword FORMAT; by default FORMAT= "%.3g".  The
 *   font type, font height and text orientation for the labels can be set
 *   with keywords FONT (default "helvetica"), HEIGHT (default 14 points)
 *   and ORIENT respectively.
 *
 *   By default the colorbar is drawn next to the current viewport; other
 *   viewport coordinates can be given by VPORT=[xmin,xmax,ymin,ymax].
 *   Keyword ADJUST can be used to move the bar closer to (adjust<0) or
 *   further from (adjust>0) the viewport.
 *
 *   Keyword COLOR can be used to specify the color of the labels, the
 *   ticks and the frame of the colorbar.  Default is foreground color.
 *
 *   Keyword WIDTH can be used to set the width of the lines used to draw
 *   the frame and the ticks of the colorbar.
 *
 *   Keyword TICKLEN can be used to set the lenght (in NDC units) of the
 *   ticks.  Default is 0.007 NDC.
 *
 *   Keyword THICKNESS can be used to set the thickness of the colorbar
 *   (in NDC units).  Default is 0.020 NDC.
 *
 *
 *  SEE ALSO: pli, plt, pldj, plg, viewport.
 */
{
  if (is_void(cmin)) cmin = min(z);
  if (is_void(cmax)) cmax = max(z);
  if (is_void(vport)) vport = viewport();
  if (is_void(adjust)) adjust = 0.0;
  if (is_void(ticklen)) ticklen = 0.007;
  if (is_void(thickness)) thickness = 0.020;
  if (is_void(nlabs)) nlabs = 11;

  local red, green, blue;
  palette, red, green, blue, query=1;
  ncolors = numberof(red);
  if (ncolors < 2) {
    ncolors = 240;
  }
  levs = span(cmin, cmax, ncolors + 1);
  cells = char(indgen(0 : ncolors - 1));

  linetype = 1; /* "solid" */

  if (vert) {
    x0 = vport(2) + adjust + 0.022;
    x1 = x0 + thickness;
    y0 = vport(3);
    y1 = vport(4);
    cells = cells(-,);
  } else {
    x0 = vport(1);
    x1 = vport(2);
    y0 = vport(3) - adjust - 0.045;
    y1 = y0 - thickness;
    cells = cells(,-);
  }
  sys = plsys(0);
  pli, cells, x0, y0, x1, y1;
  if (is_void(width) || width != 0) {
    plg, [y0,y0,y1,y1], [x0,x1,x1,x0], closed=1,
      color=color, width=width, type=linetype, marks=0;
  }

  if (nlabs) {
    if (is_void(format)) format= "%.3g";
    text = swrite(format=format, span(cmin, cmax, nlabs));

    local lx0, lx1, lx2, ly0, ly1, ly2;
    if (vert) {
      lx0 = array(x1, nlabs);
      lx1 = array(x1 + ticklen, nlabs);
      lx2 = array(x1 + 1.67*ticklen, nlabs);
      ly0 = span(y0, y1, nlabs);
      eq_nocopy, ly1, ly0;
      eq_nocopy, ly2, ly0;
      justify = "LH";
    } else {
      ly0 = array(y1, nlabs);
      ly1 = array(y1 - ticklen, nlabs);
      ly2 = array(y1 - 1.67*ticklen, nlabs);
      lx0 = span(x0, x1, nlabs);
      eq_nocopy, lx1, lx0;
      eq_nocopy, lx2, lx0;
      justify = "CT";
    }
    if (ticklen && (is_void(width) || width != 0)) {
      pldj, lx0, ly0, lx1, ly1,
        color=color, width=width, type=linetype;
    }
    for (i = 1; i <= nlabs; ++i) {
      plt, text(i), lx2(i), ly2(i), tosys=0, color=color, font=font,
        height=height, opaque=opaque, orient=orient, justify=justify;
    }
  }
  plsys, sys;
}

/*---------------------------------------------------------------------------*/
/* ESTIMATION OF DIRTY BEAM */

func mira_dirty_beam(this)
/* DOCUMENT mira_dirty_beam(this);
 *   Computes dirty beam for MiRA instance THIS.  The result has the same
 *   geometry as an image reconstructed from THIS, i.e. it depends on the
 *   u-v coverage and on the synthetic field of view parameters as set by
 *   mira_config.  There is however a small (1/2 pixel) offset in the X and
 *   Y (RA and dec) directions for for even dimensions.
 *
 * EXAMPLE:
 *   mira_config, this, dim=128, pixelsize=0.3*MIRA_MILLIARCSECOND;
 *   dirty = mira_dirty_beam(this);
 *   mira_plot_image, dirty, this;
 *
 * SEE ALSO: mira_plot_image, mira_config, fft, roll.
 */
{
  /* Get sky coordinates (in radians). */
  if (! (dim = this.dim) ||
      ! (pixelsize = this.pixelsize)) {
    error, "DIM and PIXELSIZE must be set before (see mira_config)";
  }

  /* convert (u,v) into frequels */
  s = dim*pixelsize;
  u = long(floor(s*this.u + 0.5));
  v = long(floor(s*this.v + 0.5));

  /* convert (u,v) into FFT indices */
  fmax = dim/2;
  fmin = fmax + 1 - dim;
  j = where((u >= fmin)&(u <= fmax)&(v >= fmin)&(v <= fmax));
  if (numberof(j) != numberof(u)) {
    _mira_warn, "Nyquist frequency too small";
    u = u(j);
    v = v(j);
  }
  u += (u < 0)*dim; // zero-based FFT U-index
  v += (v < 0)*dim; // zero-based FFT V-index
  w = array(double, dim, dim);
  w(1) = 1.0;
  w(1 + u + v*dim) = 1.0;
  w(1 + (dim - u)%dim + ((dim - v)%dim)*dim) = 1.0;
  return roll((1.0/(dim*dim))*double(fft(w, -1)), [dim/2, dim/2]);
}

/*---------------------------------------------------------------------------*/
/* PLOTTING OF U-V COVERAGE */

local mira_plot_baselines, mira_plot_frequencies;
/* DOCUMENT mira_plot_baselines, this;
 *     -or- mira_plot_frequencies, this;
 *
 *   Plot all observed baselines or spatial frequencies in MiRA opaque
 *   handle THIS.  The (U,V) coordinates are baselines projected onto the
 *   sky in meters or spatial frequencies projected onto the sky in cycles
 *   per radian.
 *
 *
 * KEYWORDS
 *
 *   Unless keyword KEEPLIMITS is true, the axis orientations are set to
 *   match conversions in astronomy (see mira_fix_image_axis).
 *
 *   Keywords COLOR, SYMBOL and SIZE are passed to plp (which see).
 *
 *   Keyword NOTITLE can be set true to disable axis titles.
 *
 *
 * SEE ALSO: plp, mira_new, mira_plot_baselines.
 */

func mira_plot_frequencies(this, color=, symbol=, size=, fill=,
                           notitle=, keeplimits=)
{
  if (is_void(symbol)) symbol = 4;
  if (is_void(size)) size = 0.33;
  local u, v;
  eq_nocopy, u, this.u;
  eq_nocopy, v, this.v;
  w = this.eff_wave;
  i = where((u != 0.0)|(v != 0.0));
  if (is_array(i)) {
    grow, u, -u(i);
    grow, v, -v(i);
  }
  limits, square=1;
  plp, v*1e-6, u*1e-6, color=color, symbol=symbol, size=size, fill=fill;
  if (! keeplimits) mira_fix_image_axis;
  if (! notitle) {
    xytitles,
      "u-spatial frequency (Mega-cycles/radian)",
      "v-spatial frequency (Mega-cycles/radian)";
  }
}

func mira_plot_baselines(this, color=, symbol=, size=, fill=,
                         keeplimits=)
{
  if (is_void(symbol)) symbol = 4;
  if (is_void(size)) size = 0.33;
  _mira_warn, "fix wavelength";
  w = avg(this.eff_wave);
  u = w*this.u;
  v = w*this.v;
  i = where((u != 0.0)|(v != 0.0));
  if (is_array(i)) {
    grow, u, -u(i);
    grow, v, -v(i);
  }
  limits, square=1;
  plp, v, u, color=color, symbol=symbol, size=size, fill=fill;
  if (! keeplimits) mira_fix_image_axis;
  if (! notitle) {
    xytitles, "u-baseline (meters)", "v-baseline (meters)";
  }
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func mira_polar_to_cartesian(amp, amperr, phi, phierr, what, goodman=)
/* DOCUMENT data = mira_polar_to_cartesian(amp, amperr, phi, phierr, what)
 *
 *   Convert complex data given in polar coordinates (AMP,PHI) with
 *   their standard deviations (AMPERR,PHIERR) into cartesian
 *   coordinates (RE,IM) and associated noise model.  The result is
 *   a hash table:
 *
 *     DATA.re = real part of complex data
 *     DATA.im = imaginary part of complex data
 *     DATA.crr = variance of real part of complex data
 *     DATA.cii = variance of imaginary part of complex data
 *     DATA.cri = covariance of real and imaginary parts of complex data
 *     DATA.wrr = statistical weight for real part of residuals
 *     DATA.wii = statistical weight for imaginary part of residuals
 *     DATA.wri = statistical weight for real times imaginary parts of residuals
 *
 *   The quadratic penalty writes:
 *
 *     ERR =     DATA.wrr*(DATA.re - re)^2
 *           +   DATA.wii*(DATA.im - im)^2
 *           + 2*DATA.wri*(DATA.re - re)*(DATA.im - im);
 *
 * SEE ALSO: mira_data_penalty.
 */
{
  if (min(amp) < 0.0) {
    _mira_warn, swrite(format="There are negative %s amplitudes!!! [FIXED]", what);
    k = where(amp < 0.0);
    phi(k) += MIRA_PI;
    amp(k) *= -1.0;
  }
  cos_phi = cos(phi);
  sin_phi = sin(phi);
  re = amp*cos_phi;
  im = amp*sin_phi;
  case = (amperr > 0.0) + 2*(phierr > 0.0);
  n0 = numberof((j0 = where(case == 0)));
  n1 = numberof((j1 = where(case == 1)));
  n2 = numberof((j2 = where(case == 2)));
  n3 = numberof((j3 = where(case == 3)));
  if (n0 || n1 || n2) {
    _mira_warn, swrite(format="there are %d out of %d invalid complex data",
                       (n0 + n1 + n2), numberof(case));
  }
  if (goodman) {
    /* Use Goodman approximation with SIGMA such that the area of ellipsoids
     * at one standard deviation are the same:
     *    PI*SIGMA^2 = PI*(RHO*SIGMA_THETA)*SIGMA_RHO
     * hence:
     *    SIGMA = sqrt(RHO*SIGMA_THETA*SIGMA_RHO)
     * If SIGMA_THETA or SIGMA_RHO is invalid, take:
     *    SIGMA = SIGMA_RHO
     *    SIGMA = RHO*SIGMA_THETA
     * accordingly.
     */
    wrr = array(double, dimsof(case));
    wri = array(double, dimsof(case));
    if (n3) {
      /* Valid amplitude and phase data. */
      wrr(j3) = 1.0/(amp(j3)*phierr(j3)*amperr(j3));
    }
    if (n2) {
      /* Only phase data. */
      if (0) {
        /* FIXME: requires amplitude. */
        err = amp(j2)*phierr(j2);
        j = where(err > 0.0);
        if (is_array(j)) {
          j2 = j2(j);
          err = err(j);
          wrr(j2) = 1.0/(err*err);
        }
      } else {
        err = phierr(j2);
        j = where(err > 0.0);
        if (is_array(j)) {
          j2 = j2(j);
          err = err(j);
          wrr(j2) = 1.0/(err*err);
        }
      }
    }
    if (n1) {
      /* Only amplitude data (FIXME: requires phase). */
      err = amperr(j1);
      wrr(j1) = 1.0/(err*err);
    }
    wii = wrr;
  } else {
    /* Use convex quadratic local approximation. */
    wrr = wri = wii = array(double, dimsof(case));
    if (n3) {
      /* Valid amplitude and phase data. */
      err1 = amperr(j3);
      var1 = err1*err1;
      err2 = amp(j3)*phierr(j3);
      var2 = err2*err2;
      cs = cos_phi(j3);
      sn = sin_phi(j3);
      crr = cs*cs*var1 + sn*sn*var2;
      cri = cs*sn*(var1 - var2);
      cii = sn*sn*var1 + cs*cs*var2;
      a = 1.0/(var1*var2);
      wrr(j3) =  a*cii;
      wri(j3) = -a*cri;
      wii(j3) =  a*crr;
    }
    if (n2) {
      /* Only phase data (FIXME: requires amplitude). */
      err2 = amp(j2)*phierr(j2);
      j = where(err2 > 0.0);
      if (is_array(j)) {
        j2 = j2(j);
        err2 = err2(j);
        var2 = err2*err2;
        cs = cos_phi(j2);
        sn = sin_phi(j2);
        a = 1.0/var2;
        wrr(j2) =  a*sn*sn;
        wri(j2) = -a*sn*cs;
        wii(j2) =  a*cs*cs;
      }
    }
    if (n1) {
      /* Only amplitude data (FIXME: requires phase). */
      err1 = amperr(j1);
      var1 = err1*err1;
      cs = cos_phi(j1);
      sn = sin_phi(j1);
      a = 1.0/var1;
      wrr(j1) =  a*cs*cs;
      wri(j1) =  a*cs*sn;
      wii(j1) =  a*sn*sn;
    }
  }

  return h_new(re = re, im = im,
               wrr = wrr, wri = wri, wii = wii);
}

func mira_stdev_to_weight(s)
/* DOCUMENT mira_stdev_to_weight(s)
     Compute weighting array from standard deviation S.
   SEE ALSO: mira_weight_to_stdev. */
{
  if (is_array((i = where(s > (w = array(double, dimsof(s))))))) {
    s = s(i);
    w(i) = 1.0/(s*s);
  }
  return w;
}

func mira_weight_to_stdev(w)
/* DOCUMENT mira_weight_to_stdev(w)
     Compute standard deviation from weighting array W.
   SEE ALSO: mira_stdev_to_weight. */
{
  if (is_array((i = where(w > (s = array(double, dimsof(w))))))) {
    s(i) = 1.0/sqrt(w(i));
  }
  return s;
}

func mira_relative_absolute_difference(a, b)
/* DOCUMENT mira_relative_absolute_difference(a, b)
 *   Returns elementwise relative absolute difference between A and B defined
 *   as: 0                                             if A(i) = B(i)
 *       2*abs(A(i) - B(i))/(abs(A(i)) + abs(B(i))     otherwise
 *
 * SEE ALSO: mira.
 */
{
  diff = a - b;
  rdif = array(double, dimsof(diff));
  if (! is_array((i = where(diff)))) return rdif;
  rdif(i) = diff(i)/(abs(a) + abs(b))(i);
  return rdif + rdif;
}

func mira_rescale(a, .., scale=, rgb=, cubic=)
/* DOCUMENT mira_rescale(a, dimlist)
       -or- mira_rescale(a, scale=FACT)
     Return an array obtained by interpolation of A with new dimension list
     as given by DIMLIST or with all its dimension multiplied by a scaling
     factor FACT if keyword SCALE is specified.  If keyword RGB is true the
     first dimension of A and of the interpolated array must be 3 and the
     interpolated array is converted into char.  If keyword CUBIC is true
     cubic spline interpolation is used.

   SEE ALSO: interp, spline, transpose. */
{
  /* explicit */ extern spline;

  /* Get dimension lists. */
  if (! is_array(a)) error, "unexpected non-array argument";
  old_dimlist = dimsof(a);
  ndims = old_dimlist(1);
#if 0
  if (is_void(rgb)) {
    /* can be used to automatically guess RGB image */
    rgb = (structof(a) == char && ndims == 3 && old_dimlist(2) == 3);
  }
#endif
  if (rgb && (ndims != 3 || old_dimlist(2) != 3)) {
    error, "bad dimension list for RGB image";
  }
  if (! is_void(scale)) {
    if (more_args()) error, "two many arguments with scale option";
    if (! is_scalar(scale) || ! (is_real(scale) || is_integer(scale)) ||
        scale <= 0) {
      error, "SCALE must be a strictly positive scalar";
    }
    new_dimlist = array(long, ndims + 1);
    new_dimlist(1) = ndims; // FIXME: not needed?
    if (rgb) {
      new_dimlist(2) = 3;
      k = 3;
    } else {
      k = 2;
    }
    new_dimlist(k:) = max(long(scale*old_dimlist(k:) + 0.5), 1);
  } else {
    /* Build dimension list. */
    new_dimlist = [0];
    while (more_args()) {
      local arg;
      eq_nocopy, arg, next_arg();
      if (is_integer(arg)) {
        if (is_scalar(arg)) {
          grow, new_dimlist, arg;
        } else if (is_vector(arg) && (n = numberof(arg)) == arg(1) + 1) {
          /* got a vector which is a valid dimension list */
          if (n >= 2) {
            arg = arg(2:);
            if (min(arg) <= 0) {
              error, "negative value in dimension list";
            }
            grow, new_dimlist, arg;
          }
        } else {
          error, "bad dimension list";
        }
      } else if (! is_void(arg)) {
        error, "unexpected data type in dimension list";
      }
    }
    new_dimlist(1) = numberof(new_dimlist) - 1;
    if (rgb) {
      if (new_dimlist(1) != 3 || old_dimlist(2) != 3) {
        error, "bad new dimension list for RGB image";
      }
    } else {
      if (new_dimlist(1) != ndims) {
        error, "bad number of new dimensions";
      }
    }
  }

  if (cubic) {
    /* use cubic interpolation */
    if (! is_func(spline)) {
      require, "spline.i";
    }
    if (rgb) {
      a = transpose(a, 0);
      k = 1;
    } else {
      k = 0;
    }
    while (++k <= ndims) {
      old_dimlist = dimsof(a);
      n0 = old_dimlist(2);
      n1 = new_dimlist(k + 1);
      if (n1 == 1) {
        a = a(avg,..)(..,-);
      } else {
        if (n1 != n0) {
          old_dimlist(2) = n1;
          b = array(double, old_dimlist);
          x0 = (indgen(n0) - (n0 + 1)/2.0)/n0;
          x1 = (indgen(n1) - (n1 + 1)/2.0)/n1;
          n = numberof(a)/n0;
          for (i=1 ; i<=n ; ++i) {
            b(..,i) = spline(a(..,i), x0, x1);
          }
        }
        eq_nocopy, a, b;
      }
      if (ndims > 1) a = transpose(a, 0);
    }
  } else {
    for (j = 1 ; j <= ndims ; ++j) {
      k = j + 1;
      n0 = old_dimlist(k);
      n1 = new_dimlist(k);
      if (n1 == n0) continue;
      if (n1 == 1) {
        /* fix output dimensions equal to 1 (FIXME: use transpose)*/
        /**/ if (j ==  1) a = a(avg, ..);
        else if (j ==  2) a = a(,avg, ..);
        else if (j ==  3) a = a(,,avg, ..);
        else if (j ==  4) a = a(,,,avg, ..);
        else if (j ==  5) a = a(,,,,avg, ..);
        else if (j ==  6) a = a(,,,,,avg, ..);
        else if (j ==  7) a = a(,,,,,,avg, ..);
        else if (j ==  8) a = a(,,,,,,,avg, ..);
        else if (j ==  9) a = a(,,,,,,,,avg, ..);
        else if (j == 10) a = a(,,,,,,,,,avg, ..);
        else if (j == 11) a = a(,,,,,,,,,,avg, ..);
        else if (j == 12) a = a(,,,,,,,,,,,avg, ..);
        else if (j == 13) a = a(,,,,,,,,,,,,avg, ..);
        else if (j == 14) a = a(,,,,,,,,,,,,,avg, ..);
        else if (j == 15) a = a(,,,,,,,,,,,,,,avg, ..);
        else if (j == 16) a = a(,,,,,,,,,,,,,,,avg, ..);
        else if (j == 17) a = a(,,,,,,,,,,,,,,,,avg, ..);
        else if (j == 18) a = a(,,,,,,,,,,,,,,,,,avg, ..);
        else if (j == 19) a = a(,,,,,,,,,,,,,,,,,,avg, ..);
        else if (j == 20) a = a(,,,,,,,,,,,,,,,,,,,avg, ..);
        else error, "too many dimensions";
      } else {
        /* use linear interpolation */
        x0 = (1.0/n0)*(indgen(n0) - (n0 + 1)/2.0);
        x1 = (1.0/n1)*(indgen(n1) - (n1 + 1)/2.0);
        a = interp(a, x0, x1, j);
      }
    }
  }
  if (rgb) return char(max(min(floor(a + 0.5), 255.0), 0.0));
  return a;
}

func mira_recenter(x, quiet=)
/* DOCUMENT mira_recenter(x);

     Recenter model image X at its photo-centre rounded to nearest pixel.
     Argument X must have at least 2 dimensions, the first 2 dimensions of X
     are considered to be the angular direction.  Extra dimensions are ignored
     for the recentering (they can represent other coordinates for instance
     the wavelength or the time).  Along a dimensions of lenght N, the center
     is at (N - N/2)-th pixel -- with integer division -- which corresponds to
     the model of the Fourier transform assumed by MiRA.

     Unless keyword QUIET is true, the coordinates of the center get printed
     out.

   SEE ALSO: mira_solve.
 */
{
  sx = sum(x);
  if (sx <= 0.0) {
    return x;
  }

  dimlist = dimsof(x);
  n1 = dimlist(2);
  n2 = dimlist(3);
  o1 = n1 - (n1/2);
  o2 = n2 - (n2/2);
  c1 = sum(double(indgen(n1) - o1)     * x)/sx;
  c2 = sum(double(indgen(n2) - o2)(-,) * x)/sx;
  if (! quiet) {
    write, format="Offsets of photo-center: (%+.1f, %+.1f) pixels.\n", c1, c2;
  }
  i1 = lround(c1);
  i2 = lround(c2);
  if (i1 == 0 && i2 == 0) {
    return x;
  }
  if (i1 > 0) {
    dst1 = 1:n1-i1;
    src1 = 1+i1:n1;
  } else {
    dst1 = 1-i1:n1;
    src1 = 1:n1+i1;
  }
  if (i2 > 0) {
    dst2 = 1:n2-i2;
    src2 = 1+i2:n2;
  } else {
    dst2 = 1-i2:n2;
    src2 = 1:n2+i2;
  }
  xp = array(structof(x), dimsof(x));
  xp(dst1, dst2, ..) = x(src1, src2, ..);
  return xp;
}

func mira_azimuthal_average(img, x0=, y0=, scale=, profile=)
{
  if (! is_func(histo2)) {
    require, "histo.i";
  }
  dims = dimsof(img);
  xdim = dims(2);
  ydim = dims(3);
  if (is_void(x0)) x0 = 0.5*(xdim + 1);
  if (is_void(y0)) y0 = 0.5*(ydim + 1);
  if (is_void(scale)) scale = 1.0;
  x = (1.0/scale)*(indgen(xdim) - x0);
  y = (1.0/scale)*(indgen(ydim) - y0);
  r = abs(x, y(-,));
  local px;
  py = histo2(r, px, weight=img, average=1, interp=1);
  if (profile) {
    return [px, py];
  }
  return interp(py, px, r);
}

func mira_dirac(dim)
/* DOCUMENT img = mira_dirac(dim);
     Returns a 2-D DIM-by-DIM image with a point-like object approximately
     centered suitable as an initial solution for mira_solve (which see).
   SEE ALSO: mira_solve.
 */
{
  cen = dim - (dim/2); // position of center
  (img = array(double, dim, dim))(cen, cen) = 1.0;
  return img;
}

local mira_cast_real_as_complex, mira_cast_complex_as_real;
/* DOCUMENT z = mira_cast_real_as_complex(x);
 *     -or- x = mira_cast_complex_as_real(z);
 *
 *   The first function converts a 2-by-any real array X into a complex
 *   array Z such that:
 *
 *      Z.re = X(1,..)
 *      Z.im = X(2,..)
 *
 *   the second function does the inverse operation.
 *
 * SEE ALSO: reshape.
 */
func mira_cast_real_as_complex(x)
{
  local z;
  if (structof(x) != double) {
    if (! is_real(x) && ! is_integer(x)) {
      error, "bad data type (expecting real or integer)";
    }
    x = double(unref(x));
  }
  if ((ndims = (dimlist = dimsof(x))(1)) < 1 || dimlist(2) != 2) {
    error, "incompatible dimension list";
  }
  reshape, z, &x, complex, (ndims == 1 ? [0] : grow(ndims - 1, dimlist(3:0)));
  return z;
}
func mira_cast_complex_as_real(z)
{
  local x;
  if (structof(z) != complex) {
    error, "bad data type (expecting complex)";
  }
  reshape, x, &z, double, make_dimlist(2, dimsof(z));
  return x;
}

func mira_glob(pat)
/* DOCUMENT mira_glob(pat)
 *   Returns a list of files matching glob-style pattern PAT.  Only the
 *   'file' part of PAT can have wild characters.
 *
 * SEE ALSO: lsdir, strglob.
 */
{
  i = strfind("/", pat, back=1);
  if ((i = i(2)) >= 1) {
    dir = strpart(pat, 1:i);
    pat = strpart(pat, i+1:0);
  } else {
    dir = "./";
  }
  list = lsdir(dir);
  list = list(where(strglob(pat, list)));
  if (! is_void(list)) return dir + list;
}

local mira_get_one_integer, mira_get_one_real;
/* DOCUMENT mira_get_one_integer(symbol)
 *     -or- mira_get_one_integer(symbol, defval)
 *     -or- mira_get_one_real(symbol)
 *     -or- mira_get_one_real(symbol, defval)
 *
 *   Make sure SYMBOL (a simple variable reference) is a real or integer
 *   scalar.  If SYMBOL is undefined on entry, DEFVAL is used instead.
 *   These functions returns true (non-zero) if the assertion failed (for
 *   SYMBOL or, if SYMBOL is undefined, for DEFVAL); otherwise, on return,
 *   the value stored by SYMBOL is a scalar of type long (integer) or
 *   double (real).
 *
 * SEE ALSO: is_integer, is_real, is_scalar.
 */
func mira_get_one_real(&symbol, defval)
{
  if (is_void(symbol)) symbol = defval;
  if (is_scalar(symbol) && (is_real(symbol) || is_integer(symbol))) {
    symbol = double(symbol);
    return 0n;
  }
  return 1n;
}
func mira_get_one_integer(&symbol, defval)
{
  if (is_void(symbol)) symbol = defval;
  if (is_scalar(symbol) && is_integer(symbol)) {
    symbol = long(symbol);
    return 0n;
  }
  return 1n;
}

func mira_is_func(f){ return (is_func(f) || (is_hash(f) && h_evaluator(f))); }
/* DOCUMENT mira_is_func(f)
     Check whether F is an object callable as a function.
   SEE ALSO: is_func, is_hash, h_evaluator.
 */

/*---------------------------------------------------------------------------*/
/* SAVE/LOAD IMAGES */
func mira_save_image(filename, img, master, single=, overwrite=,
                     comment=, history=)
{
  fix_axis = 1n;
  if (is_void(comment)) {
    comment = "Image created by MiRA.";
  }
  if (single) {
    bitpix = -32;
    type = float;
  } else {
    bitpix = -64;
    type = double;
  }
  dims = dimsof(img);
  naxis = dims(1);
  fh = fits_create(filename, bitpix=bitpix, dimlist=dims,
                   comment=comment, history=history, overwrite=overwrite);
  cr1 = mira_get_x(master);
  cr2 = mira_get_y(master);
  pixelsize = mira_get_pixelsize(master);
  cdelt1 = (cr1(1) < cr1(0) ? pixelsize : -pixelsize);
  cdelt2 = (cr2(1) < cr2(0) ? pixelsize : -pixelsize);
  if (fix_axis) {
    /* Change axis coordinates according to astronomical conventions
       (CDELT1 < 0 and CDELT2 > 0). */
    flip = ((cdelt1 > 0.0 ? 1n : 0n) | (cdelt2 < 0.0 ? 2n : 0n));
  } else {
    flip = 0n;
  }
  if (flip == 1n) {
    cr1 = cr1(::-1);
    cdelt1 = -cdelt1;
    img = type(img)(::-1,..);
  } else if (flip == 2n) {
    cr2 = cr2(::-1);
    cdelt2 = -cdelt2;
    img = type(img)(,::-1,..);
  } else if (flip == 3n) {
    cr1 = cr1(::-1);
    cdelt1 = -cdelt1;
    cr2 = cr2(::-1);
    cdelt2 = -cdelt2;
    img = type(img)(::-1,::-1,..);
  } else if (structof(img) != type) {
    img = type(img);
  }
  i1 = abs(cr1)(mnx);
  fits_set, fh, "CRPIX1", i1, "1-based index of reference pixel along 1st axis";
  fits_set, fh, "CRVAL1", cr1(i1), "coordinate of reference pixel along 1st axis";
  fits_set, fh, "CDELT1", cdelt1, "increment along 1st axis";
  fits_set, fh, "CTYPE1", "radian", "units for 1st axis";
  i2 = abs(cr2)(mnx);
  fits_set, fh, "CRPIX2", i2, "1-based index of reference pixel along 2nd axis";
  fits_set, fh, "CRVAL2", cr2(i2), "coordinate of reference pixel along 2nd axis";
  fits_set, fh, "CDELT2", cdelt2, "increment along 2nd axis";
  fits_set, fh, "CTYPE2", "radian", "units for 2nd axis";
  fits_write_header, fh;
  fits_write_array, fh, img;
}

/*---------------------------------------------------------------------------*/
/* DIGITIZATION AND CLASSIFICATION */

func mira_digitize(data, precision)
/* DOCUMENT bin = mira_digitize(data);
         or bin = mira_digitize(data, precision);
  
     This function digitizes values in DATA and returns a hash table object
     BIN with the following members:

        BIN.index - has same dimension list as DATA and is the index
                    associated with each value (running from 1 to N, where N
                    is the number of significantly different values);
  
        BIN.count - is a N-element vector set with the number of elements in
                    each set of data values;
  
        BIN.value - is a N-element vector set with the central value in each
                    set of data values;

     The result is such that:
     
        abs(DATA - BIN.value(BIN.index)) <= 0.5*PRECISION

     where the optional absolute precision is 0 by default.

  
   SEE ALSO: mira_classify, heapsort, histogram.
 */
{
  local rounded_data;

  if (is_void(precision)) {
    precision = 0.0;
  } else if (! is_scalar(precision)
             || identof(precision) > Y_DOUBLE
             || precision < 0.0) {
    error, "bad value for PRECISION";
  }
  if (precision > 0.0) {
    /* Round DATA to PRECISION. */
    data = precision*round((1.0/precision)*data);
  }

  index = array(long, dimsof(data));
  order = heapsort(data);
  sorted_data = data(order);
  data = [];
  test = (sorted_data(dif) > 0);
  index(order) = 1 + long(test)(cum);
  order = [];
  if (anyof(test)) {
    value = sorted_data(grow(0, where(test)) + 1);
  } else {
    value = [sorted_data(1)];
  }
  test = [];
  count = histogram(index);
  return h_new(index=index, count=count, value=value);
}

func mira_classify(data, threshold)
/* DOCUMENT obj = mira_classify(data);
         or obj = mira_classify(data, threshold);
  
     Classify DATA in different cliques.  Optional THRESHOLD (default is 0) is
     the absolute minimum distance between the values taken in different
     cliques.  In words, if
  
        abs(DATA(INDEX(j)) - DATA(INDEX(j+1))) > THRESHOLD
  
     where INDEX = sort(DATA(*)), then DATA(INDEX(j)) and DATA(INDEX(j+1))
     belongs to a different clique; otherwise they belong to the same clique.
  
     The result is a hash table object such that:
  
        OBJ.region - has same dimension list as DATA and is the clique
                     index associated with each datum (running from 1 to N,
                     where N is the number of different cliques);
  
        OBJ.count  - is a N-element vector set with the number of
                     elements in each clique;
  
        OBJ.mean   - is a N-element vector set with the mean data value in
                     each clique;
  
        OBJ.stdv   - is a N-element vector set with the standard deviation
                     of data values in each clique;
  
     Beware that the classification only works for cliques well separated.
     It would fail if DATA is continuously varying.  In that case, the only
     solution is to use THRESHOLD = 0.
  
  
   SEE ALSO: mira_digitize, heapsort.
 */
{
  if (is_void(threshold)) {
    threshold = 0;
  } else if (! is_scalar(threshold) ||
             identof(threshold) > Y_DOUBLE ||
             threshold < 0) {
    error, "bad value for THRESHOLD";
  }
  order = heapsort(data);
  region = array(long, dimsof(data));
  region(order) = 1 + (data(order)(dif) > threshold)(cum);
  order = [];
  count = histogram(region);
  mean = histogram(region, data)/count;
  diff = data - mean(region);
  stdv = histogram(region, diff*diff)/count; // FIXME
  return h_new(region=region, count=count, mean=mean, stdv=stdv);
}

/*---------------------------------------------------------------------------*/
/* EMULATION OF MISSING FOR YETI FUNCTIONS */

func _mira_symlink_to_name(s) { return link(s + ""); }
if (is_func(symlink_to_name) == 2 && is_func(is_symlink) == 2 &&
    is_func(name_of_symlink) == 2 && is_func(value_of_symlink) == 2) {
  /* Yeti 6.2.1 and newer */
  _mira_symlink_to_name = [];
 } else if (is_func(link) == 2 && is_func(is_link) == 2 &&
           is_func(link_name) == 2 && is_func(solve_link) == 2) {
  /* Yeti 6.2.0 */
  _mira_warn, "old symbolic link functions (consider upgrading Yeti).";
  symlink_to_name = _mira_symlink_to_name;
  symlink_to_variable = link;
  name_of_symlink = link_name;
  value_of_symlink = solve_link;
  is_symlink = is_link;
} else {
  error, "symbolic link not implemented (upgrade Yeti)";
}

func _mira_grow_members(obj, .., flatten=)
{
  local key, value;
  if (flatten) {
    while (more_args()) {
      eq_nocopy, key, next_arg();
      h_set, obj, key, grow(h_get(obj, key), next_arg()(*));
    }
  } else {
    while (more_args()) {
      eq_nocopy, key, next_arg();
      h_set, obj, key, grow(h_get(obj, key), next_arg());
    }
  }
}
if (is_func(h_grow)) {
  _mira_grow_members = [];
} else {
  _mira_warn, "h_grow provided by MiRA (consider upgrading Yeti).";
  h_grow = _mira_grow_members;
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 78
 * coding: utf-8
 * End:
 */
