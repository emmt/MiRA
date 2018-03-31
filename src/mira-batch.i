/*
FE modified 20170922
*/
/*
 * mira-batch.i -
 *
 * Run MiRA in batch mode all arguments are parsed from the command line.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of MiRA, a "Multi-aperture Image Reconstruction
 * Algorithm", <https://github.com/emmt/MiRA>.
 *
 * Copyright (C) 2001-2016, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

if (! is_func(nfft_new)) {
  /* Attempt to pre-load YNFFT. */
  include, "nfft.i", 3;
}
_mira_batch_xform = (is_func(nfft_new) ? "nfft" : "exact");

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
 *    "hyperbolic" "mu"   real >= 0.0      -mu
 *                 "tau"  real > 0.0       -tau
 *                 "eta"  values           -eta // FIXME:
 *                 "mask" image
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
  ("Usage: ymira [OPTIONS] INPUT [...] OUTPUT",
   "Image reconstruction.  INPUT and [...] are the OI-FITS data file and OUTPUT " +
   "is the result saved into a FITS file.",
   _lst\
   (
    "\nData selection:",
    _lst("target", NULL, "NAME", OPT_STRING, "Name of the astrophysical object"),
    _lst("effwave", NULL, "LENGTH", OPT_STRING, "Effective wavelength, e.g. 1.6micron"),
    _lst("effband", NULL, "LENGTH", OPT_STRING, "Effective bandwidth, e.g. 200nm"),
    _lst("wavemin", NULL, "LENGTH", OPT_STRING, "Minimum wavelength, e.g. 1.5µm"),
    _lst("wavemax", NULL, "LENGTH", OPT_STRING, "Maximum wavelength, e.g. 1.7µm"),
    "\nImage parameters and direct model:",
    _lst("pixelsize", NULL, "ANGLE", OPT_STRING, "Angular size of pixels, e.g. 0.1mas"),
    _lst("fov", NULL, "ANGLE", OPT_STRING, "Angular size of the field of view, e.g. 20mas"),
    _lst("dim", NULL, "NUMBER", OPT_INTEGER, "Number of pixels per side of the image"),
    _lst("normalization", NULL, "VALUE", OPT_REAL, "Flux normalization (sum of pixels = VALUE)"),
    _lst("min", NULL, "LOWER", OPT_REAL, "Lower bound for the pixel values"),
    _lst("max", NULL, "UPPER", OPT_REAL, "Upper bound for the pixel values"),
    _lst("xform", _mira_batch_xform, "NAME", OPT_STRING, "Method to compute the Fourier transform"),
    "\nRegularization settings:",
    _lst("regul", NULL, "NAME", OPT_STRING, "Name of regularization method (use -regul=help for more information)"),
    _lst("mu", 0.0, "VALUE(S)", OPT_REAL_LIST, "Global regularization weight"),
    _lst("tau", 1E-6, "VALUE", OPT_REAL, "Edge preserving threshold"),
    _lst("eta", 1.0, "VALUE(S)", OPT_REAL_LIST, "Gradient scales along dimensions"),
    _lst("gamma", NULL, "FWHM", OPT_STRING, "A priori full half width at half maximum, e.g. 15mas"),
 /* FE 20170922 */
    _lst("regul_threshold", 1E-3, "VALUE", OPT_REAL, "Threshold for the cost function"),
    _lst("regul_epsilon", 1E-6, "VALUE", OPT_REAL, "Edge preserving threshold"),
    _lst("regul_isotropic", NULL, NULL, OPT_FLAG, "Use isotropic version of the regularization"),
    _lst("no_vis", NULL, NULL, OPT_FLAG, "Not not use complex visibilities"),
    _lst("no_vis2", NULL, NULL, OPT_FLAG, "Not not use square visibilities"),
    _lst("no_t3", NULL, NULL, OPT_FLAG, "Not not use closure quantities"),
/* FE end */
     /*
    _lst("regul_power", 2.0, "VALUE", OPT_REAL, "Power for the Lp-norm"),
    _lst("regul_normalized", NULL, NULL, OPT_FLAG, "Image is normalized"),
    _lst("regul_cost", "l2", "NAME", OPT_STRING, "Cost function for the regularization"),
    _lst("regul_type", "log", "NAME", OPT_STRING, "Subtype for the regularization"),
    _lst("regul_periodic", NULL, NULL, OPT_FLAG, "Use periodic conditions for the regularization"),
    */
    "\nInitial image:",
    _lst("initial", "random", "NAME", OPT_STRING, "FITS file or method for initial image"),
    _lst("seed", NULL, "VALUE", OPT_REAL, "Seed for the random generator"),
    "\nOutput image:",
    _lst("overwrite", NULL, NULL, OPT_FLAG, "Overwrite output if it exists"),
    _lst("bitpix", -32, "BITPIX", OPT_INTEGER, "Bits per pixel"),
    _lst("save_initial", NULL, NULL, OPT_FLAG, "Save initial image as a secondary HDU in the result"),
    "\nReconstruction strategy:",
    _lst("bootstrap", NULL, "COUNT", OPT_INTEGER, "Number of bootstrapping iterations"),
    _lst("recenter", NULL, NULL, OPT_FLAG, "Recenter result of bootstrapping iterations"),
    _lst("threshold", NULL, "FRACTION", OPT_REAL, "Level for soft-thresholding input image(s)"),
    "\nMessages and graphics:",
    _lst("quiet", NULL, NULL, OPT_FLAG, "Suppress most messages"),
    _lst("verb", NULL, "COUNT", OPT_INTEGER, "Verbose level"),
    _lst("view", 0, "MASK", OPT_INTEGER, "Bitwise mask to specify which graphics to show"),
    "\nOptimizer Settings:",
    _lst("mem", NULL, "COUNT", OPT_INTEGER, "Number of previous steps to memorize in VMLMB"),
    _lst("ftol", NULL, "REAL", OPT_REAL, "Function tolerance for the global convergence"),
    _lst("gtol", NULL, "REAL", OPT_REAL, "Gradient tolerance for the global convergence"),
    _lst("maxiter", NULL, "COUNT", OPT_INTEGER, "Maximum number of iterations"),
    _lst("maxeval", NULL, "COUNT", OPT_INTEGER, "Maximum number of evaluations of the objective function"),
    "\nLine search:",
    _lst("sftol", NULL, "REAL", OPT_REAL, "Function tolerance for the line search"),
    _lst("sgtol", NULL, "REAL", OPT_REAL, "Gradient tolerance for the line search"),
    _lst("sxtol", NULL, "REAL", OPT_REAL, "Step tolerance for the line search"),
    "\nMiscellaneous:",
    _lst("debug", NULL, NULL, OPT_FLAG, "Debug mode"),
    _lst("help", NULL, NULL, OPT_HELP, "Print out this help"),
    _lst("version", MIRA_VERSION, NULL, OPT_VERSION, "Print out version number")));

func _mira_cli_parse_angle(str, opt)
{
  if (is_void(opt)) {
    opt = "angle";
  }
  dummy = units = string();
  value = 0.0;
  if (sread(str, value, units, dummy) != 2) {
    opt_error, "Expecting value and units for " + opt;
  }
  fact = mira_parse_angular_units(units, 0);
  if (fact == 0) {
    opt_error, "Invalid units for " + opt;
  }
  return value*fact;
}

func _mira_cli_parse_length(str, opt)
{
  if (is_void(opt)) {
    opt = "length";
  }
  dummy = units = string();
  value = 0.0;
  if (sread(str, value, units, dummy) != 2) {
    opt_error, "Expecting value and units for " + opt;
  }
  fact = mira_parse_length_units(units, 0);
  if (fact == 0) {
    opt_error, "Invalid units for " + opt;
  }
  return value*fact;
}

func mira_format(val) { return sum(print(val)); }
/* DOCUMENT mira_format(arg);
     Yield a string representation of ARG.

     SEE ALSO; print.
*/

func mira_main(argv0, argv)
{
  /* Constants and shortcuts. */
  FALSE = 0n;
  TRUE = 1n;
  format = mira_format;

  opt = opt_parse(_MIRA_OPTIONS, argv);
  if (is_void(opt)) {
    /* Options "-help", or "-usage", or "-version" have been set. */
    return;
  }
  if (opt.regul == "help") {
    write, format="\n%s\n", "Available regularizations:";
    write, format="\n  -regul=%s -tau=...\n  %s\n  %s\n",
      "hyperbolic_smoothness",
      "Edge-preserving smoothness with hyperbolic norm and threshold TAU,",
      "approximates total variation (TV) with TAU very small.";
    write, format="\n  -regul=%s -tau=...\n  %s\n",
      "Huber_smoothness",
      "Edge-preserving smoothness with Huber semi-norm and threshold TAU.";
    write, format="\n  -regul=%s\n  %s\n",
      "quadratic_smoothness",
      "Quadratic smoothness.";
    write, format="\n  -regul=%s -gamma=...\n  %s\n",
      "quadratic_compactness",
      "Quadratic compactness with full-width at half maximum GAMMA.";
    write, format="\n  -regul=%s -gamma=...\n  %s\n",
      "Huber_compactness",
      "L1-L2 compactness with Huber norm and full-width at half maximum GAMMA.";
    write, format="\n  -regul=%s -gamma=...\n  %s\n",
      "hyperbolic_compactness",
      "L1-L2 compactness with Huber norm and full-width at half maximum GAMMA.";
    write, format="%s\n", "";
    quit;
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

  /* Check bitpix. */
  if (opt.bitpix != 8 && opt.bitpix != 16 && opt.bitpix != 32 &&
      opt.bitpix != -32 && opt.bitpix != -64) {
    opt_error, "Invalid value for `-bitpix`";
  }

  /* Initial comment. */
  comment = ("Image reconstructed by MiRA algorithm" +
             " <https://github.com/emmt/MiRA>");
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
                 regul_name, format(opt.mu), opt.tau, format(opt.eta));
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
               regul_name, format(opt.mu), opt.tau, format(opt.eta));
    } else if (regul_name == "compactness") {
      /* Compactness prior. */
      if (is_void(opt.gamma)) {
        opt_error, "Option `-gamma=...` must be specified with \"compactness\" regularization";
      }
      value = _mira_cli_parse_angle(opt.gamma, "`-gamma=...`");
      if (value <= 0) {
        opt_error, "Invalid value for `-gamma=...`";
      }
      h_set, opt, gamma=value;
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s and GAMMA=%g",
               regul_name, format(opt.mu), opt.gamma);
      regul_post = TRUE;

 /* begin FE 20170922 */
   } else if (regul_name == "l2l1_smoothness") {
      /* l2l1_smoothness option added by FE 20170822. */
      if (opt.mu <= 0.0) {
        opt_error, "Value of `-mu` must be > 0.0";
      }
      if (opt.regul_threshold <= 0.0) {
        opt_error, "Value of `-regul_threshold` must be > 0.0";
      }
      regul = rgl_new("l2l1_smoothness","threshold",opt.regul_threshold);
      h_set, opt, mu=opt.mu;
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s, REGUL_THRESHOLD=%s",
               regul_name, format(opt.mu), format(opt.regul_threshold));

   } else if (regul_name == "totvar") {
      /* total variance option added by FE 20170922 */
      if (opt.regul_epsilon <= 0.0) {
        opt_error, "Value of `-regul_epsilon` must be > 0.0";
      }
      regul = rgl_new("totvar","epsilon",opt.regul_epsilon,"isotropic",opt.regul_isotropic);
      h_set, opt, mu=opt.mu;
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s, REGUL_EPSILON=%s, REGUL_ISOTROPIC=%s",
               regul_name, format(opt.mu), format(opt.regul_epsilon), format(opt.regul_isotropic));
 /* end FE 20170922 */

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
  local dim, pixelsize;
  if (! is_void(opt.pixelsize)) {
    pixelsize = _mira_cli_parse_angle(opt.pixelsize, "`-pixelsize=...`");
    if (pixelsize <= 0) {
      opt_error, "Invalid value for `-pixelsize=...`";
    }
  }
  choice = ((is_void(opt.fov) ? 0 : 1) |
            (is_void(opt.dim) ? 0 : 2));
  if (choice == 1) {
    fov = _mira_cli_parse_angle(opt.fov, "`-fov=...`");
    if (fov <= 0) {
      opt_error, "Invalid value for `-fov=...`";
    }
    if (! is_void(opt.pixelsize)) {
      dim = lround(fov/pixelsize);
    }
  } else if (choice == 2) {
    if (opt.dim <= 0) {
      opt_error, "Bad value for `-dim=...`";
    }
    dim = opt.dim;
  } else if (choice == 3) {
    opt_error, "Only one of `-fov=...` or `-dim=...` can be specified";
  } else if (! initial_filename) {
    opt_error, ("Image dimensions must be specified with `-dim=...` "
                + "or `-fov=..` when the initial image is not from a file");
  }

  /* Check some options. */
  if (! is_void(opt.threshold) && (opt.threshold < 0 || opt.threshold >= 1)) {
    opt_error, "Invalid value for `-threshold=...`";
  }

  /* Initial image. */
  local initial;
  if (initial_filename) {
    /* Read initial image and fix orientation. */
    img = mira_read_image(initial_filename);
    naxis = img.naxis;
    for (i = 1; i <= naxis; ++i) {
      cdelt_i = swrite(format="cdelt%d", i);
      if (h_get(img, cdelt_i) < 0) {
        /**/ if (i ==  1) h_set, img, arr = img.arr(::-1,..);
        else if (i ==  2) h_set, img, arr = img.arr(,::-1,..);
        else if (i ==  3) h_set, img, arr = img.arr(,,::-1,..);
        else if (i ==  4) h_set, img, arr = img.arr(,,,::-1,..);
        else if (i ==  5) h_set, img, arr = img.arr(,,,,::-1,..);
        else if (i ==  6) h_set, img, arr = img.arr(,,,,,::-1,..);
        else if (i ==  7) h_set, img, arr = img.arr(,,,,,,::-1,..);
        else if (i ==  8) h_set, img, arr = img.arr(,,,,,,,::-1,..);
        else if (i ==  9) h_set, img, arr = img.arr(,,,,,,,,::-1,..);
        else if (i == 10) h_set, img, arr = img.arr(,,,,,,,,,::-1,..);
        crpix_i = swrite(format="crpix%d", i);
        naxis_i = swrite(format="naxis%d", i);
        h_set, img,
          cdelt_i, -h_get(img, cdelt_i),
          crpix_i, h_get(img, naxis_i) + 1 - h_get(img, crpix_i);
      }
    }

    /* Check number of axis. */
    if (naxis < 2 || naxis > 3) {
      opt_error, "Expecting a 2D/3D initial image";
    }
    naxis1 = img.naxis1;
    naxis2 = img.naxis2;
    naxis3 = (naxis >= 3 ? img.naxis3 : 1);
    if (naxis == 3) {
      if (naxis3 != 1) {
        write, format="WARNING - %s\n",
          "Converting 3D image into a 2D grayscaled image";
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
    if (is_void(dim)) {
      if (is_void(fov)) {
        fov1 = cdelt1*naxis1;
        fov2 = cdelt2*naxis2;
        fov = max(fov1, fov2);
      }
      dim = lround(fov/pixelsize);
    }
    cunit = "mas";
    cdelt = mira_convert_units(siunit, cunit)*pixelsize;
    crpix = (dim/2) + 1;
    crval = 0.0;
    img = mira_resample_image(img, pad=0, norm=1n,
                              naxis1=dim,   naxis2=dim,
                              crpix1=crpix, crpix2=crpix,
                              crval1=crval, crval2=crval,
                              cdelt1=cdelt, cdelt2=cdelt,
                              cunit1=cunit, cunit2=cunit);
    eq_nocopy, initial, img.arr;

  } else if (initial_random) {
    if (! is_void(opt.seed)) {
      if (opt.seed <= 0.0 || opt.seed >= 1.0) {
        opt_error, "Seed value must be between 0.0 and 1.0 non-inclusive";
      }
      random_seed, opt.seed;
    }
    initial = random(dim, dim);
  } else if (initial_dirac) {
    initial = array(double, dim, dim);
    initial((dim/2) + 1, (dim/2) + 1) = 1.0;
  } else {
    opt_error, "Option -initial=\""+initial_name+"\" not yet implemented";
  }

  /* Parse wavelenght settings. */
  keys = ["effband", "effwave", "wavemin", "wavemax"];
  for (k = 1; k <= numberof(keys); ++k) {
    local str;
    key = keys(k);
    eq_nocopy, str, h_get(opt, key);
    if (! is_void(str)) {
      value = _mira_cli_parse_length(str, "`-"+key+"=...`");
      if (value <= 0) {
        opt_error, "Invalid value for `-"+key+"=...`";
      }
      h_set, opt, key, value;
    }
  }

  /* Read input data. */

  master = mira_new(argv(1:-1),
                    target = opt.target,
                    eff_wave = opt.effwave,
                    eff_band = opt.effband,
                    wavemin = opt.wavemin,
                    wavemax = opt.wavemax,
/* FE 20170922 Start */
                    no_vis = opt.no_vis,
                    no_vis2 = opt.no_vis2,
                    no_t3 = opt.no_t3,
/* FE 20170922 End */
                    quiet = opt.quiet);

  mira_config, master, dim=dim, pixelsize=pixelsize, xform=opt.xform;


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
                     view = opt.view,
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

  /* OP20180213 get normalized dirty beam*/
  dirty=mira_dirty_beam(master);
  dirty = dirty / max(dirty);
  fh = mira_save_image(mira_wrap_image(dirty, master), "dirty_"+final_filename,
                       overwrite=opt.overwrite, bitpix=opt.bitpix,
                       comment=comment);

  /* Define a clean beam from this dirty beam */
  print,opt.pixelsize;
  print,opt.dim;
  pixelsize = _mira_cli_parse_angle(opt.pixelsize, "`-pixelsize=...`");

  x = span(-pixelsize*opt.dim/2,pixelsize*opt.dim/2,opt.dim)(,-:1:opt.dim);
  y = transpose(x);
  /* restrict clean beam to <4mas */
  res_op=4.848e-9 * 4;
  clean = dirty * (dirty>0.1) * (abs(x,y)<res_op);

  imgS = convoln(img, clean);
  imgS = imgS / max(imgS);

  fh = mira_save_image(mira_wrap_image(imgS, master), "Convol_"+final_filename,
                       overwrite=opt.overwrite, bitpix=opt.bitpix,
                       comment=comment);
  /* OP20180213 End of modification*/

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

if (batch()) {
  pldefault, style="boxed.gs", marks=0, width=1, palette="gray.gp",
    maxcolors=240, legends=0;
  argv = get_argv();
  argv = (numberof(argv) >= 2 ? argv(2:) : []); /* strip executable path */
  mira_main, current_include(), argv;
  quit;
}
