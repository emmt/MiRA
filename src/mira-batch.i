/*
 * mira-batch.i -
 *
 * Run MiRA in batch mode all arguments are parsed from the command line.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2001-2016, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

/* The following initialization function is executed once. */
func _mira_batch_init {
  if (! is_func(mira_require) || is_void(MIRA_HOME)) {
    /* Locate installation directory and include main MiRA code. */
    path = current_include();
    dir = (i = strfind("/", path, back=1)(2)) > 0 ? strpart(path, 1:i) : "./";
    include, dir + "mira.i", 1;
  }
}
_mira_batch_init;
_mira_batch_init = [];

mira_require, "opt_init", MIRA_HOME + ["", "../lib/ylib/"] + "options.i";

func mira_central_index(dim) { return dim/2; }
/* DOCUMENT mira_central_index(dim);
     Get the index of the central element along a dimension of length DIM.
     The same conventions as for fftshift are used.
 */

func mira_extract_region(arr, dim)
/* DOCUMENT mira_extract_region(arr, dim);

     Extract the central DIM-by-DIM region of 2D or 3D array ARR.
*/
{
  if (is_array(arr)) {
    dims = dimsof(arr);
    rank = dims(1);
  } else {
    dims = [];
    rank = -1;
  }
  if (rank != 2 && rank != 3) {
    error, "expecting a 2D or 3D array";
  }
  if (! is_scalar(dim) || ! is_integer(dim) || dim < 1) {
    error, "invalid dimension";
  }

  size = dims(1:2);
  off = (size/2) - (dim/2);
  i1 = max(   1 + off, 1);
  i2 = min(size + off, size);
  o1 = max(   1 - off, 1);
  o2 = min( dim - off, dim);
  if (max(o1) == 1 && min(o2) == dim) {
    /* No needs for zero-padding. */
    return arr(i1(1):i2(1),i1(2):i2(2),..);
  } else {
    dims(2:3) = dim;
    dst = array(structof(arr), dims);
    if (allof(o1 <= o2)) {
      dst(o1(1):o2(1),o1(2):o2(2),..) = arr(i1(1):i2(1),i1(2):i2(2),..);
    }
    return dst;
  }
}

func mira_save_result(master, img, filename, overwrite=, bitpix=)
{
  polychromatic = 0n;
  dim = mira_get_dim(master);
  pixelsize = mira_get_pixelsize(master);
  dims = dimsof(img);
  if (polychromatic) {
    if (dims(1) != 3) {
      error, "expecting a 3D array";
    }
  } else {
    if (dims(1) != 2) {
      error, "expecting a 2D array";
    }
  }
  width = dims(2);
  height = dims(3);
  depth = (polychromatic ? dims(4) : 1);
  if (width != dim || height != dim) {
    error, "incompatible dimensions";
  }

  if (is_void(bitpix)) bitpix = -32;
  fh = fits_create(filename, overwrite=overwrite, bitpix=bitpix, extend=1);

  /* see http://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html */
  crpix1 = crpix2 = 0.0;
  ctype1 = ctype2 = "milliarcsecond";
  delta =  pixelsize/MIRA_MILLIARCSECOND;

  fits_set, fh, "CRPIX1", mira_central_index(dim), "coordinate system reference pixel";
  fits_set, fh, "CRVAL1", 0.0, "coordinate system value at reference pixel";
  fits_set, fh, "CDELT1", -delta, "coordinate increment along axis";
  fits_set, fh, "CTYPE1", "milliarcsecond", "name of the coordinate axis";

  fits_set, fh, "CRPIX2", mira_central_index(dim), "coordinate system reference pixel";
  fits_set, fh, "CRVAL2", 0.0, "coordinate system value at reference pixel";
  fits_set, fh, "CDELT2", +delta, "coordinate increment along axis";
  fits_set, fh, "CTYPE2", "milliarcsecond", "name of the coordinate axis";

  if (polychromatic) {
    eff_wave = mira_get_eff_wave(master);
    if (numberof(eff_wave) == 1) {
      stp_wave = 0.0;
    } else {
      stp_wave = (eff_wave(0) - eff_wave(1))/(numberof(eff_wave) - 1.0);
    }
    fits_set, fh, "CRPIX3", 1.0, "coordinate system reference pixel";
    fits_set, fh, "CRVAL3", eff_wave(1), "coordinate system value at reference pixel";
    fits_set, fh, "CDELT3", stp_wave, "coordinate increment along axis";
    fits_set, fh, "CTYPE3", "wavelength", "name of the coordinate axis";
  }

  fits_set, fh, "HISTORY", "This image has been reconstructed by MiRA";
  fits_set, fh, "HISTORY", "freely available at <https://github.com/emmt/MiRA>";
  fits_write_header, fh;
  fits_write_array, fh, img;
  fits_pad_hdu, fh;
  fits_close, fh;
}

/* regularization
 * clique (mu, region)
 * entropy (mu, type, normalized, prior, epsilon)
 * l2l1_smoothness (mu, threshold)
 * lpnorm (mu, power, epsilon)
 * qsmooth (mu, prior)
 * quadratic (mu, A, b, W)
 * roughness (mu, threshold, cost, periodic)
 * simple (mu, threshold, cost)
 * smoothness (mu)
 * totvar (mu, epsilon, isotropic)
 * xsmooth (mu, threshold, cost, dimlist, isotropic)
 */
/*
 *
 *    regul
 *
 *    "totvar"     "mu"         real >= 0.0      regul_mu
 *                 "epsilon"    real > 0.0       regul_epsilon
 *                 "isotropic"  true/false       regul_isotropic
 *
 *    "lpnorm"     "mu"         real >= 0.0      regul_mu
 *                 "power"      real > 0.0       regul_power
 *                 "epsilon"    real > 0.0       regul_epsilon
 *
 *    "roughness"  "mu"         real >= 0.0      regul_mu
 *                 "threshold"  real > 0.0       regul_threshold
 *                 "cost"       name             regul_cost
 *                 "periodic"   true/false       regul_periodic
 *
 *    "entropy"    "mu"         real >= 0.0      regul_mu
 *                 "type"       "sqrt" or "log"  regul_type
 *                 "prior"      "none"           regul_cost
 *                 "normalized" true/false       regul_normalized
 *                 "epsilon"    real > 0.0       regul_epsilon
 *
 *    possible costs are:  "l1", "l2", "l2l1", "l2l0", or "cauchy";
 *                          default is "l2" (i.e. quadratic)
 */

NULL = []; /* FIXME: make a private function to hide local variables */
_MIRA_OPTIONS = opt_init\
  ("Usage: mira [OPTIONS] INPUT [...] OUTPUT",
   "Image reconstruction.  INPUT and [...] are the OI-FITS data file and OUTPUT" +
   "is the result saved into a FITS file.",
   _lst(
        /* DATA SELECTION */
        _lst("target",           NULL,  "NAME",   OPT_STRING,   "Name of the astrophysical object"),
        _lst("effwave",          NULL,  "VALUE",  OPT_REAL,     "Effective wavelength (microns)"),
        _lst("effband",          NULL,  "VALUE",  OPT_REAL,     "Effective bandwidth (microns)"),
        _lst("minwave",          NULL,  "VALUE",  OPT_REAL,     "Minimum wavelength (microns)"),
        _lst("maxwave",          NULL,  "VALUE",  OPT_REAL,     "Maximum wavelength (microns)"),
        /* IMAGE PARAMETERS */
        _lst("pixelsize",        NULL,  "VALUE",  OPT_REAL,     "Size of the pixel (milliarcseconds)"),
        _lst("dim",              NULL,  "NUMBER", OPT_INTEGER,  "Number of pixels per side of the image"),
        _lst("xform",           "nfft", "NAME",   OPT_STRING,   "Method to compute the Fourier transform"),
        _lst("normalization",    NULL,  "VALUE",  OPT_REAL,     "Flux normalization"),
        _lst("min",              NULL,  "LOWER",  OPT_REAL,     "Lower bound for the flux"),
        _lst("max",              NULL,  "UPPER",  OPT_REAL,     "Upper bound for the flux"),
        _lst("overwrite",        NULL,  NULL,     OPT_FLAG,     "Overwrite output if it exists"),
        /* REGULARIZATION SETTINGS */
        _lst("regul",            NULL,  "NAME",  OPT_STRING,    "Name of regularization method"),
        _lst("regul_mu",         0.0,   "VALUE", OPT_REAL,      "Global regularization weight"),
        _lst("regul_epsilon",    1E-6,  "VALUE", OPT_REAL,      "Edge preserving threshold"),
        _lst("regul_threshold",  1E-3,  "VALUE", OPT_REAL,      "Threshold for the cost function"),
        _lst("regul_power",      2.0,   "VALUE", OPT_REAL,      "Power for the Lp-norm"),
        _lst("regul_normalized", NULL,  NULL,    OPT_FLAG,      "Image is normalized"), /* FIXME: */
        _lst("regul_cost",       "l2",  "NAME",  OPT_STRING,    "Cost function for the regularization"),
        _lst("regul_type",       "log", "NAME",  OPT_STRING,    "Subtype for the regularization"),
        _lst("regul_periodic",   NULL,  NULL,    OPT_FLAG,      "Use periodic conditions for the regularization"),
        _lst("regul_isotropic",  NULL,  NULL,    OPT_FLAG,      "Use isotropic version of the regularization"),
        /* ALGORITHM PARAMETERS */
        _lst("initial",          "random", "NAME", OPT_STRING, "FITS file or method for initial image"),
        _lst("seed",             NULL, "VALUE",    OPT_REAL,   "Seed for the random generator"),
        /* can also be "random", "Dirac", "Gauss", or "Cauchy-Lorentz" */
        _lst("bootstrap", 20,   "COUNT", OPT_INTEGER, "Number of bootstrapping iterations"),
        _lst("maxiter",   NULL, "COUNT", OPT_INTEGER, "Maximum number of iterations"),
        _lst("maxeval",   NULL, "COUNT", OPT_INTEGER, "Maximum number of evaluations of the objective function"),
        _lst("verb",      NULL, "COUNT", OPT_INTEGER, "Verbose level"),
        _lst("view",      0,    "MASK",  OPT_INTEGER, "Bitwise mask to specify which graphics to show"),
        _lst("mem",       NULL, "COUNT", OPT_INTEGER, "Number of previous steps to memorize in VMLMB"),
        _lst("ftol",      NULL, "REAL",  OPT_REAL,    "Function tolerance for the global convergence"),
        _lst("gtol",      NULL, "REAL",  OPT_REAL,    "Gradient tolerance for the global convergence"),
        _lst("sftol",     NULL, "REAL",  OPT_REAL,    "Function tolerance for the line search"),
        _lst("sgtol",     NULL, "REAL",  OPT_REAL,    "Gradient tolerance for the line search"),
        _lst("sxtol",     NULL, "REAL",  OPT_REAL,    "Step tolerance for the line search"),
        /* MISCELLANEOUS */
        _lst("help",     NULL, NULL, OPT_HELP, "Print out this help"),
        _lst("version",  MIRA_VERSION, NULL, OPT_VERSION, "Print out version number")));
/* FIXME: merge regul_cost and regul_type */

func mira_scale_option(value, factor)
{
  if (is_void(value)) return;
  return factor*value;
}

func mira_main(argv0, argv)
{
  FALSE = 0n;
  TRUE = 1n;

  opt = opt_parse(_MIRA_OPTIONS, argv);
  if (is_void(opt)) {
    /* Options "--help", or "--usage", or "--version" have been set. */
    return;
  }
  argc = numberof(argv);
  if (argc < 2) {
    opt_usage, _MIRA_OPTIONS;
    return;
  }
  final_filename = argv(0);

  /* Setup the regularization. */
  regul_name = opt.regul;
  if (! is_void(regul_name)) {
    if (opt.regul_mu < 0.0) {
      opt_error, "value of `--regul_mu` must be >= 0.0";
    }
    if (opt.regul_threshold <= 0.0) {
      opt_error, "value of `--regul_threshold` must be > 0.0";
    }
    if (opt.regul_epsilon <= 0.0) {
      opt_error, "value of `--regul_epsilon` must be > 0.0";
    }
    if (opt.regul_power <= 0.0) {
      opt_error, "value of `--regul_power` must be > 0.0";
    }
    regul = rgl_new(regul_name);
    regul_keywords = rgl_info(regul_name);
    for (k = numberof(regul_keywords); k >= 1; --k) {
      name = regul_keywords(k);
      option = "regul_"+name;
      value = h_get(opt, option);
      if (is_void(value)) {
        if (! h_has(opt, option)) {
          if (regul_name == "entropy" && name == "prior") {
            /* skip this one */
            continue;
          }
          opt_error, "unsupported regularization \""+regul_name+"\"";
        }
        value = FALSE;
      }
      rgl_config, regul, name, value;
    }
  }

  initial_random = FALSE;
  initial_dirac = FALSE;
  initial_gauss = FALSE;
  initial_cauchy = FALSE;
  initial_filename = [];
  initial_name = opt.initial;
  if (initial_name == "random") {
    initial_random = TRUE;
  } else if  (initial_name == "Gauss") {
    initial_gauss = TRUE;
  } else if  (initial_name == "Dirac") {
    initial_dirac = TRUE;
  } else if  (initial_name == "Cauchy-Lorentz") {
    initial_cauchy = TRUE;
  } else {
    initial_filename = initial_name;
  }

  local dim, pixelsize;
  if (is_void(opt.pixelsize)) {
    if (! initial_filename) {
      opt_error, "option -pixelsize=VALUE or -inital=FILENAME must be specified";
    }
  } else {
    if (opt.pixelsize <= 0.0) {
      opt_error, "bad value for PIXELSIZE";
    }
    pixelsize = opt.pixelsize*MIRA_MILLIARCSECOND;
  }
  if (is_void(opt.dim)) {
    if (! initial_filename) {
      opt_error, "option -dim=NUMBER or -inital=FILENAME must be specified";
    }
  } else {
    if (opt.dim <= 0) {
      opt_error, "bad value for DIM";
    }
    dim = opt.dim;
  }

  if (initial_filename) {
    fh = fits_open(initial_filename, 'r');
    naxis = fits_get(fh, "NAXIS");
    if (naxis != 2) {
      opt_error, "expecting a 2-D initial image";
    }
    naxis1 = fits_get(fh, "NAXIS1");
    naxis2 = fits_get(fh, "NAXIS2");
    if (! is_void(dim) && dim != max(naxis1, naxis2)) {
      opt_error, "option --dim cannot be specified with --initial=FILENAME";
    }
    dim = max(naxis1, naxis2);
    if (dim % 2) {
      ++dim;
    }

    crpix1 = fits_get(fh, "CRPIX1");
    if (is_void(crpix1)) crpix1 = mira_central_index(dim);
    crval1 = fits_get(fh, "CRVAL1");
    if (is_void(crval1)) crval1 = 0.0;
    cdelt1 = fits_get(fh, "CDELT1");

    crpix2 = fits_get(fh, "CRPIX2");
    if (is_void(crpix2)) crpix2 = mira_central_index(dim);
    crval2 = fits_get(fh, "CRVAL2");
    if (is_void(crval2)) crval2 = 0.0;
    cdelt2 = fits_get(fh, "CDELT2");

    if (is_void(pixelsize)) {
      if (is_void(cdelt1) || is_void(cdelt2)) {
        opt_error, "PIXELSIZE must be specified (CDELT# missing in FITS file)";
      }
      pixelsize = avg(abs(cdelt1), abs(cdelt2));
      if (max(abs(abs(cdelt1) - pixelsize),
              abs(abs(cdelt2) - pixelsize)) > 1E-6*pixelsize) {
        opt_error, "|CDELT1| and |CDELT2| must be equal in FITS file";
      }
    } else {
      /* Overwrite the pixel size written in the FITS file. */
      cdelt1 = ((is_void(cdelt1) || cdelt1 < 0.0) ? -pixelsize :  pixelsize);
      cdelt2 = ((is_void(cdelt2) || cdelt2 > 0.0) ?  pixelsize : -pixelsize);
    }

    initial = fits_read_array(fh);
    fits_close, fh;
    fh = [];

    /* Fix the orientation of the image and its size. */
    if ((cdelt1 > 0.0) || (cdelt2 < 0.0)) {
      local flip1, flip2;
      if (cdelt1 > 0.0) flip1 = ::-1;
      if (cdelt2 < 0.0) flip2 = ::-1;
      initial = initial(flip1, flip2);
    }
    if (dim > naxis1 || dim >= naxis2) {
      require, MIRA_HOME+"img.i";
      initial = img_pad(initial, dim, dim, just=1);
    }

  } else if (initial_random) {
    if (! is_void(opt.seed)) {
      if (opt.seed <= 0.0 || opt.seed >= 1.0) {
        opt_error, "seed value must be between 0.0 and 1.0 non-inclusive";
      }
      random_seed, opt.seed;
    }
    initial = random(dim, dim);
  } else if (initial_dirac) {
    initial = array(double, dim, dim);
    initial(dim/2 +1, dim/2 + 1) = 1.0;
  } else {
    opt_error, "option --initial=\""+initial_name+"\" not yet implemented";
  }

  /* FIXME: take the default pixelsize from the data file */

  master = mira_new(argv(1:-1),
                    target = opt.target,
                    eff_wave = mira_scale_option(opt.effwave, 1E-6),
                    eff_band = mira_scale_option(opt.effband, 1E-6),
                    minwave = mira_scale_option(opt.minwave, 1E-6),
                    maxwave = mira_scale_option(opt.maxwave, 1E-6),
                    monochromatic = 1n);

  mira_config, master, dim=dim, pixelsize=pixelsize, xform=opt.xform;

  local maxeval, maxiter;
  if (opt.bootstrap < 0) {
    opt_error, "bad value for option --bootstrap";
  }
  maxiter0 = opt.bootstrap;
  if (! is_void(opt.maxiter)) {
    if (opt.maxiter < 0) {
      opt_error, "bad value for option --maxiter";
    }
    maxiter = opt.maxiter;
    if (maxiter > maxiter0) {
      maxiter0 = maxiter;
      maxiter = 0;
    }
  }
  maxeval0 = 2*maxiter0;
  if (! is_void(opt.maxeval)) {
    if (opt.maxeval < 0) {
      opt_error, "bad value for option --maxeval";
    }
    maxeval = opt.maxeval;
    if (maxeval > maxeval0) {
      maxeval0 = maxeval;
      maxeval = 0;
    }
  }

  /* Bootstrap. */
  local final;
  if (maxeval0 >= 1 && maxiter0 >= 1) {
    final = mira_recenter(mira_solve(master, initial,
                                     maxeval = maxeval0,
                                     maxiter = maxiter0,
                                     verb = opt.verb,
                                     view = opt.view,
                                     xmin = opt("min"),
                                     xmax = opt("max"),
                                     normalization = opt.normalization,
                                     regul = regul,
                                     mem = opt.mem,
                                     ftol = 0.0,
                                     gtol = 0.0,
                                     sftol = opt.sftol,
                                     sgtol = opt.sgtol,
                                     sxtol = opt.sxtol));
  } else {
    eq_nocopy, final, initial;
  }

  /* Final reconstruction. */
  final = mira_solve(master, final,
                     maxeval = opt.maxeval,
                     maxiter = opt.maxiter,
                     verb = opt.verb,
                     view = opt.view,
                     xmin = opt("min"),
                     xmax = opt("max"),
                     normalization = opt.normalization,
                     regul = regul,
                     mem = opt.mem,
                     ftol = opt.ftol,
                     gtol = opt.gtol,
                     sftol = opt.sftol,
                     sgtol = opt.sgtol,
                     sxtol = opt.sxtol);


  /* Save the result. */

  polychromatic = 0;
  if (is_void(opt.bitpix)) h_set, opt, bitpix = -32;  // FIXME:
  dim = mira_get_dim(master);
  pixelsize = mira_get_pixelsize(master);
  naxis = 2;
  naxis1 = naxis2 = dim;
  dimlist = [naxis,naxis1,naxis2];
  fh = fits_create(final_filename, dimlist=dimlist,
                   overwrite=opt.overwrite,
                   bitpix=opt.bitpix, extend=1);

  /* see http://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html */
  cdelt1 = -pixelsize/MIRA_MILLIARCSECOND;
  cdelt2 =  pixelsize/MIRA_MILLIARCSECOND;

  for (hdu = 1; hdu <= 2; ++hdu) {
    if (hdu > 1) {
      fits_new_hdu, fh, "IMAGE";
      fits_set, fh, "BITPIX", opt.bitpix, "bits per pixel";
      fits_set_dims, fh, dimlist;

    }
    fits_set, fh, "CRPIX1", 0.5*naxis1, "coordinate system reference pixel";
    fits_set, fh, "CRVAL1", 0.0, "coordinate system value at reference pixel";
    fits_set, fh, "CDELT1", cdelt1, "[rd] coordinate increment along axis";
    fits_set, fh, "CTYPE1", "RA", "right ascension";

    fits_set, fh, "CRPIX2", 0.5*naxis2, "coordinate system reference pixel";
    fits_set, fh, "CRVAL2", 0.0, "coordinate system value at reference pixel";
    fits_set, fh, "CDELT2", cdelt2, "coordinate increment along axis";
    fits_set, fh, "CTYPE2", "DEC", "declination";

    if (hdu == 1) {
      fits_set, fh, "COMMENT", "Image reconstructed by MiRA.";
      // FIXME: add regularization parameters, etc.
    } else {
      fits_set, fh, "COMMENT", "Initial image used by MiRA.";
      // FIXME: add regularization parameters, etc.
    }

    fits_write_header, fh;
    fits_write_array, fh, (hdu == 1 ? final : initial);
    fits_pad_hdu, fh;
  }
  fits_close, fh;


  if (opt.view && batch()) {
    write, format="%s\n", "Hit Ctrl-C twice to finish.";
    while (1) {
      pause, 5000;
    }
  }
}

if (batch()) {
  pldefault, style="boxed.gs", marks=0, width=1, palette="gray.gp",
    maxcolors=240, legends=0;
  argv = get_argv();
  argv = (numberof(argv) >= 2 ? argv(2:) : []); /* strip executable path */
  mira_main, current_include(), argv;
  quit;
}
