/*
 * mira2_batch.i -
 *
 * Run MiRA in batch mode all arguments are parsed from the command line.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of MiRA, a "Multi-aperture Image Reconstruction
 * Algorithm", <https://github.com/emmt/MiRA>.
 *
 * Copyright (C) 2001-2018, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 * MiRA is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License version 2 as published by the Free
 * Software Foundation.
 *
 * MiRA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

/* The following initialization function is executed once. */
local MIRA_BATCH_HOME;
func _mira_batch_init(path)
{
  path = filepath(path);
  i = strfind("/", path, back=1)(2);
  dir = strpart(path, 1:i);
  path = dir+"mira2.i";
  if (! open(path, "r", 1)) {
    path = dir+"../src/mira2.i";
    if (! open(path, "r", 1)) {
      throw, "\"mira2.i\" not found";
    }
  }
  include, path, 1;
  return dir;
}
MIRA_BATCH_HOME = _mira_batch_init(current_include());
_mira_batch_init = [];


/* Make sure "options" package is loaded. */
mira_require, "opt_init", "ylib", "options.i";

if (! is_func(nfft_new)) {
  /* Attempt to pre-load YNFFT. */
  include, "nfft.i", 3;
}
_mira_batch_xform = (is_func(nfft_new) ? "nfft" : "separable");

_MIRA_OPTIONS = opt_init\
  ("Usage: mira2 [OPTIONS] INPUT [...] OUTPUT\n",
   "Image reconstruction.  INPUT and [...] are the OI-FITS data file and\n"+
   "OUTPUT is the result saved into a FITS file.\n\n"+
   "Environment variables MIRA_SRCDIR and MIRA_YORICK may be set to specify\n"+
   "the directory where are installed the sources and the path to the\n"+
   "Yorick interpreter.",
   _lst\
   (
    "\nData selection:",
    _lst("target", [], "NAME", OPT_STRING,
         "Name of the astrophysical object"),
    _lst("effwave", [], "LENGTH", OPT_STRING,
         "Effective wavelength, e.g. 1.6micron"),
    _lst("effband", [], "LENGTH", OPT_STRING,
         "Effective bandwidth, e.g. 200nm"),
    _lst("wavemin", [], "LENGTH", OPT_STRING,
         "Minimum wavelength, e.g. 1.5µm"),
    _lst("wavemax", [], "LENGTH", OPT_STRING,
         "Maximum wavelength, e.g. 1.7µm"),
    "\nWhich data to fit and how:",
    _lst("visamp", "yes", "yes|no", OPT_STRING,
         "Fit complex visibility amplitudes?"),
    _lst("visphi", "yes", "yes|no", OPT_STRING,
         "Fit complex visibility phases?"),
    _lst("vis2", "yes", "yes|no", OPT_STRING,
         "Fit powerspectrum data?"),
    _lst("t3amp", "no", "yes|no", OPT_STRING,
         "Fit bispectrum amplitudes?"),
    _lst("t3phi", "yes", "yes|no", OPT_STRING,
         "Fit bispectrum phases?"),
    _lst("convexcost", "yes", "yes|no", OPT_STRING,
         "Use convex approximation for fitting complex data?"),
    _lst("phasecost", "vonmises", "vonmises|haniff|convexlimit", OPT_STRING,
         "Co-log-likelihood approximation to use for phase data"),
    "\nImage parameters and direct model:",
<<<<<<< HEAD
    _lst("pixelsize", [], "ANGLE", OPT_STRING,
         "Angular size of pixels, e.g. 0.1mas"),
    _lst("fov", [], "ANGLE", OPT_STRING,
         "Angular size of the field of view, e.g. 20mas"),
    _lst("imagesize", [], "DIM|DIM1xDIM2", OPT_STRING,
         "Number of pixels per side of the image"),
    _lst("normalization", [], "VALUE", OPT_REAL,
         "Flux normalization (sum of pixels = VALUE)"),
    _lst("min", [], "LOWER", OPT_REAL,
         "Lower bound for the pixel values"),
    _lst("max", [], "UPPER", OPT_REAL,
         "Upper bound for the pixel values"),
    _lst("xform", _mira_batch_xform, "nfft|separable|nonseparable", OPT_STRING,
         "Method to compute the Fourier transform"),
    _lst("smearingfunction", "none", "none|sinc|gauss", OPT_STRING,
         "Method to model the effects of the spectral bandwidth smearing"),
    _lst("smearingfactor", 1.0, "VALUE", OPT_REAL,
         "Factor to scale the effects of the spectral bandwidth smearing"),
    "\nRegularization settings:",
    _lst("regul", [], "NAME", OPT_STRING,
         "Name of regularization method (-regul=help for more information)"),
    _lst("mu", 0.0, "VALUE(S)", OPT_REAL_LIST,
         "Global regularization weight"),
    _lst("tau", 1E-6, "VALUE", OPT_REAL,
         "Edge preserving threshold"),
    _lst("eta", 1.0, "VALUE(S)", OPT_REAL_LIST,
         "Gradient scales along dimensions"),
    _lst("gamma", [], "FWHM", OPT_STRING,
         "A priori full half width at half maximum, e.g. 15mas"),
    "\nInitial image:",
    _lst("initial", "random", "NAME", OPT_STRING,
         "FITS file or method for initial image"),
    _lst("seed", [], "VALUE", OPT_REAL,
         "Seed for the random generator"),
    "\nOutput image:",
    _lst("overwrite", [], [], OPT_FLAG,
         "Overwrite output if it exists"),
    _lst("bitpix", -32, "BITPIX", OPT_INTEGER,
         "Bits per pixel"),
    _lst("save_initial", [], [], OPT_FLAG,
         "Save initial image as a secondary HDU in the result"),
    "\nReconstruction strategy:",
    _lst("bootstrap", [], "COUNT", OPT_INTEGER,
         "Number of bootstrapping iterations"),
    _lst("recenter", [], [], OPT_FLAG,
         "Recenter result of bootstrapping iterations"),
    _lst("threshold", [], "FRACTION", OPT_REAL,
         "Level for soft-thresholding input image(s)"),
    "\nMessages:",
    _lst("quiet", [], [], OPT_FLAG,
         "Suppress most messages"),
    _lst("verb", [], "COUNT", OPT_INTEGER,
         "Verbose level"),
    "\nOptimizer Settings:",
    _lst("mem", [], "COUNT", OPT_INTEGER,
         "Number of previous steps to memorize in VMLMB"),
    _lst("ftol", [], "REAL", OPT_REAL,
         "Function tolerance for the global convergence"),
    _lst("gtol", [], "REAL", OPT_REAL,
         "Gradient tolerance for the global convergence"),
    _lst("maxiter", [], "COUNT", OPT_INTEGER,
         "Maximum number of iterations"),
    _lst("maxeval", [], "COUNT", OPT_INTEGER,
         "Maximum number of evaluations of the objective function"),
    "\nLine search:",
    _lst("sftol", [], "REAL", OPT_REAL,
         "Function tolerance for the line search"),
    _lst("sgtol", [], "REAL", OPT_REAL,
         "Gradient tolerance for the line search"),
    _lst("sxtol", [], "REAL", OPT_REAL,
         "Step tolerance for the line search"),
    "\nMiscellaneous:",
    _lst("debug", [], [], OPT_FLAG,
         "Debug mode"),
    _lst("help", [], [], OPT_HELP,
         "Print out this help"),
    _lst("version", MIRA_VERSION, [], OPT_VERSION,
         "Print out version number")));

func _mira_is_nonnegative(x) { return x >= 0; }
func _mira_is_strictly_positive(x) { return x > 0; }

func mira_get_angle(opt, key, check)
/* DOCUMENT val = mira_get_angle(opt, key);
         or val = mira_get_angle(opt, key, check);

     yields the value of an angular parameter from the field KEY in option
     table OPT.  If the value is not a valid angle, or if CHECK is provided
     and CHECK(VAL) is not true, an error is thrown.

   SEE ALSO: opt_error, mira_get_length, mira_get_yesno.
 */
{
  local val;
  eq_nocopy, val, h_get(opt, key);
  if (! is_void(val)) {
    if (! mira_angle(val)) {
      opt_error, "Expecting value and angular units for `-"+key+"=...`";
    }
    if (is_func(check) && ! check(val)) {
      opt_error, "Invalid value for `-"+key+"=...`";
    }
    return val;
  }
}
errs2caller, mira_get_angle;

func mira_get_length(opt, key, check)
/* DOCUMENT val = mira_get_length(opt, key);
         or val = mira_get_length(opt, key, check);

     yields the value of a length parameter from the field KEY in option table
     OPT.  If the value is not a valid length, or if CHECK is provided and
     CHECK(VAL) is not true, an error is thrown.

   SEE ALSO: opt_error, mira_get_angle, mira_get_yesno.
 */
{
  local val;
  eq_nocopy, val, h_get(opt, key);
  if (! is_void(val)) {
    if (! mira_length(val)) {
      opt_error, "Expecting value and length units for `-"+key+"=...`";
    }
    if (is_func(check) && ! check(val)) {
      opt_error, "Invalid value for `-"+key+"=...`";
    }
    return val;
  }
}
errs2caller, mira_get_length;

func mira_get_yesno(opt, key)
/* DOCUMENT bool = mira_get_yesno(opt, key);
     yields true/false depending whether field KEY is "yes" or "no" (or not
     set) in option table OPT; otherwise, throw an error.

   SEE ALSO: opt_error, mira_get_angle, mira_get_length.
 */
{
  local str;
  eq_nocopy, str, h_get(opt, key);
  if (is_scalar(str) && is_string(str)) {
    if (str == "yes") {
      return 1n;
    } else if (str == "no") {
      return 0n;
    }
  } else if (is_void(str)) {
    return 0n;
  }
  opt_error, "Invalid value for `-"+key+"=...`, expecting `yes` or `no`";
}
errs2caller, mira_get_yesno;

func mira_main(argv0, argv)
{
  /* Constants and shortcuts. */
  FALSE = 0n;
  TRUE = 1n;
  get_angle = mira_get_angle;
  get_length = mira_get_length;
  get_yesno = mira_get_yesno;
  nonnegative = _mira_is_nonnegative;
  strictly_positive = _mira_is_strictly_positive;

  /* Pre-parse options. */
  arguments = mira_concatenate_arguments(argv);
  opt = opt_parse(_MIRA_OPTIONS, argv);
  if (is_void(opt)) {
    /* Options "-help", or "-usage", or "-version" have been set. */
    return;
  }
  if (opt.regul == "help") {
    write, format="\n%s\n", "Available regularizations:";
    write, format="\n  -regul=%s -tau=...\n  %s\n  %s\n",
      "hyperbolic",
      "Edge-preserving smoothness with hyperbolic norm and threshold TAU,",
      "approximates total variation (TV) with TAU very small.";
    write, format="\n  -regul=%s -gamma=...\n  %s\n",
      "compactness",
      "Quadratic compactness with full-width at half maximum GAMMA.";
    write, format="%s\n", "";
    if (batch()) {
      quit;
    }
  }
  argc = numberof(argv);
  if (argc < 2) {
    opt_usage, _MIRA_OPTIONS;
    return;
  }
  final_filename = argv(0);
  if (! opt.debug) {
    error = opt_error;
  }
  if (! opt.overwrite && open(final_filename, "r", 1)) {
    opt_error, ("output file \""+final_filename+"\" already exists");
  }

  /* Data fidelity functions and data selection. */
  flags = 0;
  if (get_yesno(opt, "visamp")) {
    flags |= MIRA_FIT_VISAMP;
  }
  if (get_yesno(opt, "visphi")) {
    flags |= MIRA_FIT_VISPHI;
  }
  if (get_yesno(opt, "vis2")) {
    flags |= MIRA_FIT_VIS2;
  }
  if (get_yesno(opt, "t3amp")) {
    flags |= MIRA_FIT_T3AMP;
  }
  if (get_yesno(opt, "t3phi")) {
    flags |= MIRA_FIT_T3PHI;
  }
  if (get_yesno(opt, "convexcost")) {
    flags |= MIRA_CONVEX_APPROX;
  }
  if (opt.phasecost == "vonmises") {
    flags |= MIRA_VON_MISES_APPROX;
  } else if (opt.phasecost == "haniff") {
    flags |= MIRA_HANIFF_APPROX;
  } else if (opt.phasecost == "convexlimit") {
    flags |= MIRA_CONVEX_LIMIT;
  } else {
    opt_error, "Bad value for `-phasecost`, expecting `vonmises`, `haniff` or `convexlimit`";
  }

  /* Check bitpix. */
  if (opt.bitpix != 8 && opt.bitpix != 16 && opt.bitpix != 32 &&
      opt.bitpix != -32 && opt.bitpix != -64) {
    opt_error, "Invalid value for `-bitpix`";
  }

  /* Initial comment. */
  comment = ("Image reconstructed by MiRA algorithm" +
             " <https://github.com/emmt/MiRA>");
  grow, comment, "Arguments: "+arguments;

  /* Setup the regularization. */
  regul_post = FALSE;
  regul_name = opt.regul;
  if (! is_void(regul_name)) {
    if (min(opt.mu) < 0.0) {
      opt_error, "Value(s) of `-mu` must be >= 0.0";
    }

#if 0
    if (strpart(regul_name, 1:11) == "hyperbolic_") {
      regul_loss = "hyperbolic";
      regul_type = strpart(regul_name, 12:0);
    } else if (strpart(regul_name, 1:6) == "Huber_") {
      regul_loss = "Huber";
      regul_type = strpart(regul_name, 7:0);
    } else if (strpart(regul_name, 1:10) == "quadratic_") {
      regul_loss = "quadratic";
      regul_type = strpart(regul_name, 11:0);
    } else {
      regul_loss = "quadratic";
      regul_type = regul_name;
    }
    if (regul_loss != "quadratic" && opt.tau <= 0.0) {
      opt_error, "Value of `-tau` must be > 0.0";
    }
    if (regul_type == "smoothness") {
      /* Quadratic or edge-preserving smoothness prior. */
      regul = rgl_new("newsmoothness");
      if (regul_loss == "quadratic") {
        grow, comment, swrite(format="Regularization: \"%s\" with MU=%s",
                              regul_name);
      } else {
        if (min(opt.eta) <= 0.0) {
          opt_error, "Values of `-eta` must be > 0.0";
        }
        rgl_config, regul, tau=opt.tau, eta=opt.eta, loss=regul_loss;
        grow, comment,
          swrite(format="Regularization: \"%s\" with MU=%s, TAU=%g, ETA=%s",
                 regul_name, mira_format(opt.mu), opt.tau,
                 mira_format(opt.eta));
      }
    } else {
      ...;
    }
#endif

    if (regul_name == "hyperbolic") {
      /* Quadratic or edge-preserving smoothness prior. */
      regul = rgl_new("hyperbolic");
      if (opt.tau <= 0.0) {
        opt_error, "Value of `-tau` must be > 0.0";
      }
      if (min(opt.eta) <= 0.0) {
        opt_error, "Values of `-eta` must be > 0.0";
      }
      rgl_config, regul, tau=opt.tau, eta=opt.eta; //, loss=regul_loss;
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s, TAU=%g, ETA=%s",
               regul_name, mira_format(opt.mu), opt.tau, mira_format(opt.eta));
    } else if (regul_name == "compactness") {
      /* Compactness prior. */
      if (is_void(opt.gamma)) {
        opt_error, "Option `-gamma=...` must be specified with \"compactness\" regularization";
      }
      value = get_angle(opt, "gamma", strictly_positive);
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s and GAMMA=%s",
               regul_name, mira_format(opt.mu), opt.gamma);
      h_set, opt, gamma=value;
      regul_post = TRUE;

#if 0
    } else {
      if (opt.regul_threshold <= 0.0) {
        opt_error, "Value of `-regul_threshold` must be > 0.0";
      }
      if (opt.regul_epsilon <= 0.0) {
        opt_error, "Value of `-regul_epsilon` must be > 0.0";
      }
      if (opt.regul_power <= 0.0) {
        opt_error, "Value of `-regul_power` must be > 0.0";
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
            opt_error, "Unsupported regularization \""+regul_name+"\"";
          }
          value = FALSE;
        }
        rgl_config, regul, name, value;
      }
#endif
    } else {
      opt_error, "Unsupported regularization \""+regul_name+"\"";
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

  /* Get image dimensions and pixel size, if specified. */
  local dim1, dim2, dims, pixelsize, imagesize, fov;
  pixelsize = get_angle(opt, "pixelsize", strictly_positive);
  fov = get_angle(opt, "fov", strictly_positive);
  eq_nocopy, imagesize, opt.imagesize;
  choice = ((is_void(fov)       ? 0 : 1) |
            (is_void(pixelsize) ? 0 : 2) |
            (is_void(imagesize) ? 0 : 4));
  if ((choice & 4) != 0) {
    dim1 = 0;
    dim2 = 0;
    dummy = "";
    if (strmatch(imagesize, "x")) {
      n = sread(imagesize, format="%d x %d %s", dim1, dim2, dummy);
    } else if (strmatch(imagesize, "×")) {
      n = sread(imagesize, format="%d × %d %s", dim1, dim2, dummy);
    } else {
      if (sread(imagesize, format="%d %s", dim1, dummy) == 1) {
        dim2 = dim1;
        n = 2;
      }
    }
    if (n != 2 || dim1 <= 0 || dim2 <= 0) {
      opt_error, "Bad value for `-imagesize=...`";
    }
    dims = [2, dim1, dim2]
  }
  if (choice == 3) {
    dim1 = dim2 = lround(fov/pixelsize);
    dims = [2, dim1, dim2];
  } else if (choice == 5) {
    pixelsize = fov/max(dim1, dim2);
  } else if (choice == 7) {
    opt_error, ("At most two of `-fov=...`, `-pixelsize=...` or " +
                "`-imagesize=...` can be specified");
  } else if (! initial_filename) {
    opt_error, ("When the initial image is not from a file, image " +
                "dimensions must be specified with exactly two of " +
                "`-fov=...`, `-pixelsize=...` or `-imagesize=...`");
  }

  /* Check some options. */
  if (! is_void(opt.threshold) && (opt.threshold < 0 || opt.threshold >= 1)) {
    opt_error, "Invalid value for `-threshold=...`";
  }

  /* Initial image. */
  local initial;
  if (initial_filename) {
    /* Read initial image. */
    img = mira_read_image(initial_filename);
    naxis = img.naxis;
    if (naxis < 2 || naxis > 3) {
      opt_error, "Expecting a 2D/3D initial image";
    }
    naxis1 = img.naxis1;
    naxis2 = img.naxis2;
    naxis3 = (naxis >= 3 ? img.naxis3 : 1);
    if (naxis == 3) {
      if (naxis3 != 1) {
        warn, "Converting 3D image into a 2D grayscaled image";
      }
      h_set, img, arr = img.arr(,,avg), naxis=2;
      h_pop, img, "naxis3";
      h_pop, img, "crpix3";
      h_pop, img, "crval3";
      h_pop, img, "cdelt3";
      h_pop, img, "cunit3";
      h_pop, img, "ctype3";
      naxis = 2;
    }

    /* Maybe resample initial image. */
    siunit = "radian"; // pixelsize and FOV are in SI units
    cdelt1 = mira_convert_units(img.cunit1, siunit)*img.cdelt1;
    cdelt2 = mira_convert_units(img.cunit2, siunit)*img.cdelt2;
    if (is_void(pixelsize)) {
      pixelsize = min(cdelt1, cdelt2);
    }
    if (is_void(dims)) {
      if (is_void(fov)) {
        fov1 = cdelt1*naxis1;
        fov2 = cdelt2*naxis2;
        fov = max(fov1, fov2);
        dim1 = lround(fov1/pixelsize);
        dim2 = lround(fov2/pixelsize);
      } else {
        dim1 = lround(fov/pixelsize);
        dim2 = dim1;
      }
      dims = [2, dim1, dim2];
    }
    cunit = "mas";
    cdelt = mira_convert_units(siunit, cunit)*pixelsize;
    crpix1 = (dim1/2) + 1;
    crpix2 = (dim2/2) + 1;
    crval = 0.0;
    img = mira_resample_image(img, pad=0, norm=1n,
                              naxis1=dim1,   naxis2=dim2,
                              crpix1=crpix1, crpix2=crpix2,
                              crval1=crval,  crval2=crval,
                              cdelt1=cdelt,  cdelt2=cdelt,
                              cunit1=cunit,  cunit2=cunit);
    eq_nocopy, initial, img.arr;

  } else if (initial_random) {
    if (! is_void(opt.seed)) {
      if (opt.seed <= 0.0 || opt.seed >= 1.0) {
        opt_error, "Seed value must be between 0.0 and 1.0 non-inclusive";
      }
      random_seed, opt.seed;
    }
    initial = random(dims);
  } else if (initial_dirac) {
    initial = array(double, dim1, dim2);
    initial((dim1/2) + 1, (dim2/2) + 1) = 1.0;
  } else {
    opt_error, "Option -initial=\""+initial_name+"\" not yet implemented";
  }

  /* Parse wavelenght settings. */
  wavemin = get_length(opt, "wavemin", strictly_positive);
  wavemax = get_length(opt, "wavemax", strictly_positive);
  effwave = get_length(opt, "effwave", strictly_positive);
  effband = get_length(opt, "effband", nonnegative);
  choice = ((is_void(effwave) ? 0 : 1) |
            (is_void(effband) ? 0 : 2) |
            (is_void(wavemin) ? 0 : 4) |
            (is_void(wavemax) ? 0 : 8));
  if (choice == 3) {
    if (2*effband > effwave) {
      opt_error, "Too large value for `-effband=...`";
    }
    wavemin = effwave - 0.5*effband;
    wavemax = effwave + 0.5*effband;
  } else if (choice == 12) {
    if (wavemin > wavemax) {
      opt_error, "Incompatible values for `-wavemin=...` and `-wavemax=...` ";
    }
  } else if (choice != 0 && choice != 4 && choice != 8) {
    opt_error, "Specify `-effwave` and `-effband`, or `-wavemin` anr/or `-wavemax`";
  }

  /* Read input data. */
  master = mira_new(argv(1:-1),
                    target = opt.target,
                    wavemin = wavemin,
                    wavemax = wavemax,
                    quiet = opt.quiet,
                    pixelsize = pixelsize,
                    dims = dims,
                    flags = flags,
                    xform = opt.xform,
                    smearingfunction = opt.smearingfunction,
                    smearingfactor = opt.smearingfactor);
  inform, "flags = "+mira_format_flags(master.flags);

  /* Post operations for the regularization. */
  if (regul_post) {
    if (regul_name == "compactness") {
      regul = mira_new_compactness_prior(master, opt.gamma);
    }
  }

  local mu;
  if (is_void(opt.bootstrap)) {
    eq_nocopy, mu, opt.mu;
    n = numberof(mu);
  } else {
    if (opt.bootstrap < 0) {
      opt_error, "Bad value for option `-bootstrap`";
    }
    n = opt.bootstrap + 1;
    m = numberof(opt.mu);
    if (m == n) {
      eq_nocopy, mu, opt.mu;
    } else if (m == 1) {
      mu = array(opt.mu, n);
    } else if (m == 2 && n >= 2) {
      mu = spanl(opt.mu(1), opt.mu(2), n);
    } else if (m == 3 && n >= 3) {
      mu = grow(spanl(opt.mu(1), opt.mu(2), n - 1), opt.mu(3));
    } else if (m == 3 && n == 2) {
      mu = [opt.mu(1), opt.mu(3)];
    } else {
      opt_error, "Bad number of values for option `-mu`";
    }
  }
  if (! is_void(opt.maxiter) && opt.maxiter < 0) {
      opt_error, "Bad value for option `-maxiter`";
  }
  if (! is_void(opt.maxeval) && opt.maxeval < 0) {
    opt_error, "Bad value for option `-maxeval`";
  }

  /* Reconstructions. */
  local final;
  eq_nocopy, img, initial;
  for (k = 1; k <= n; ++k) {
    if (! is_void(opt.threshold) && opt.threshold > 0) {
      img = mira_soft_threshold(img, opt.threshold, opt.normalization);
    }
    img = mira_solve(master, img,
                     maxeval = opt.maxeval,
                     maxiter = opt.maxiter,
                     verb = opt.verb,
                     xmin = opt("min"),
                     xmax = opt("max"),
                     normalization = opt.normalization,
                     regul = regul,
                     mu = mu(k),
                     mem = opt.mem,
                     ftol = opt.ftol,
                     gtol = opt.gtol,
                     sftol = opt.sftol,
                     sgtol = opt.sgtol,
                     sxtol = opt.sxtol);
    if (opt.recenter && k < n) {
      img = mira_recenter(img, quiet=opt.quiet);
    }
  }

  /* Save the result. */
  fh = mira_save_image(mira_wrap_image(img, master), final_filename,
                       overwrite=opt.overwrite, bitpix=opt.bitpix,
                       comment=comment);
  if (opt.save_initial) {
    mira_save_image, initial, fh, bitpix=opt.bitpix,
      extname="INITIAL_IMAGE",
      comment="Initial image used by MiRA";
  }
  fits_close, fh;


  if (opt.view && batch()) {
    /* Just make a substancial pause? */
    write, format="\n%s\n", "Close the graphic window to finish.";
    while (window_exists(_mira_window)) {
      pause, 100;
    }
  }
}

func mira_concatenate_arguments(argv)
/* DOCUMENT mira_concatenate_arguments(argv);

     yields command line arguments ARGV as a single string with common style
     (i.e., double dashes for options) and strings with spaces protected by
     double quotes.

   SEE ALSO:
 */
{
  spaces = "*[ \t\n\r\v\f]*";
  str = "";
  opt = 1n;
  for (i = 1; i <= numberof(argv); ++i) {
    arg = argv(i);
    if (opt) {
      if (strpart(arg, 1:1) != "-") {
        opt = 0n;
      } else if (arg == "--") {
        /* last option */
        opt = 0n;
      } else if (arg != "-" && strpart(arg, 1:2) != "--") {
        arg = "-" + arg;
      }
    }
    if (strglob(spaces, arg)) {
      /* Simple attempt to protect spaces (not double quotes though). */
      if (opt) {
        i = strfind("=", arg)(1);
        if (i >= 0) {
          arg = strpart(arg, 1:i+1)+"\""+strpart(arg, i+2:0)+"\"";
        }
      } else {
        arg = "\""+arg+"\"";
      }
    }
    if (strlen(str) > 0) {
      str += " " + arg;
    } else {
      eq_nocopy, str, arg;
    }
  }
  return str;
}

if (batch()) {
  pldefault, style="boxed.gs", marks=0, width=1, palette="gray.gp",
    maxcolors=240, legends=0;
  argv = get_argv();
  argv = (numberof(argv) >= 2 ? argv(2:) : []); /* strip executable path */
  mira_main, current_include(), argv;
  quit;
}
