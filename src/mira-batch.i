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
 *    "hyperbolic" "mu"   real >= 0.0      --mu
 *                 "tau"  real > 0.0       --tau
 *                 "eta"  values           --eta // FIXME:
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
  ("Usage: mira [OPTIONS] INPUT [...] OUTPUT",
   "Image reconstruction.  INPUT and [...] are the OI-FITS data file and OUTPUT " +
   "is the result saved into a FITS file.",
   _lst\
   (/* DATA SELECTION */
    _lst("target", NULL, "NAME", OPT_STRING, "Name of the astrophysical object"),
    _lst("effwave", NULL, "VALUE", OPT_STRING, "Effective wavelength, e.g. 1.6micron"),
    _lst("effband", NULL, "VALUE", OPT_STRING, "Effective bandwidth, e.g. 200nm"),
    _lst("wavemin", NULL, "VALUE", OPT_STRING, "Minimum wavelength, e.g. 1.5µm"),
    _lst("wavemax", NULL, "VALUE", OPT_STRING, "Maximum wavelength, e.g. 1.7µm"),
    /* IMAGE PARAMETERS */
    _lst("pixelsize", NULL, "ANGLE", OPT_STRING, "Angular size of pixels, e.g. 0.1mas"),
    _lst("fov", NULL, "ANGLE", OPT_STRING, "Angular size of the field of view, e.g. 20mas"),
    _lst("dim", NULL, "NUMBER", OPT_INTEGER, "Number of pixels per side of the image"),
    _lst("xform", "nfft", "NAME", OPT_STRING, "Method to compute the Fourier transform"),
    _lst("normalization", NULL, "VALUE", OPT_REAL, "Flux normalization (sum of pixels = VALUE)"),
    _lst("min", NULL, "LOWER", OPT_REAL, "Lower bound for the pixel values"),
    _lst("max", NULL, "UPPER", OPT_REAL, "Upper bound for the pixel values"),
    _lst("overwrite", NULL, NULL, OPT_FLAG, "Overwrite output if it exists"),
    /* REGULARIZATION SETTINGS */
    _lst("regul", NULL, "NAME", OPT_STRING, "Name of regularization method"),
    _lst("mu", 0.0, "VALUE(S)", OPT_REAL_LIST, "Global regularization weight"),
    _lst("tau", 1E-6, "VALUE", OPT_REAL, "Edge preserving threshold"),
    _lst("eta", 1.0, "VALUE(S)", OPT_REAL_LIST, "Gradient scales along dimensions"),
    _lst("gamma", NULL, "FWHM", OPT_STRING, "A priori full half width at half maximum, e.g. 15mas"),
    /*
    _lst("regul_epsilon", 1E-6, "VALUE", OPT_REAL, "Edge preserving threshold"),
    _lst("regul_threshold", 1E-3, "VALUE", OPT_REAL, "Threshold for the cost function"),
    _lst("regul_power", 2.0, "VALUE", OPT_REAL, "Power for the Lp-norm"),
    _lst("regul_normalized", NULL, NULL, OPT_FLAG, "Image is normalized"),
    _lst("regul_cost", "l2", "NAME", OPT_STRING, "Cost function for the regularization"),
    _lst("regul_type", "log", "NAME", OPT_STRING, "Subtype for the regularization"),
    _lst("regul_periodic", NULL, NULL, OPT_FLAG, "Use periodic conditions for the regularization"),
    _lst("regul_isotropic", NULL, NULL, OPT_FLAG, "Use isotropic version of the regularization"),
    */
    /* ALGORITHM PARAMETERS */
    _lst("initial", "random", "NAME", OPT_STRING, "FITS file or method for initial image"),
    _lst("seed", NULL, "VALUE", OPT_REAL, "Seed for the random generator"),
    _lst("save_initial", NULL, NULL, OPT_FLAG, "Save initial image as a secondary HDU in result"),
    /* can also be "random", "Dirac", "Gauss", or "Cauchy-Lorentz" */
    _lst("bootstrap", NULL, "COUNT", OPT_INTEGER, "Number of bootstrapping iterations"),
    _lst("recenter", NULL, NULL, OPT_FLAG, "Recenter result of bootstrapping iterations"),
    _lst("maxiter", NULL, "COUNT", OPT_INTEGER, "Maximum number of iterations"),
    _lst("maxeval", NULL, "COUNT", OPT_INTEGER, "Maximum number of evaluations of the objective function"),
    _lst("quiet", NULL, NULL, OPT_FLAG, "Suppress most messages"),
    _lst("verb", NULL, "COUNT", OPT_INTEGER, "Verbose level"),
    _lst("view", 0, "MASK", OPT_INTEGER, "Bitwise mask to specify which graphics to show"),
    _lst("mem", NULL, "COUNT", OPT_INTEGER, "Number of previous steps to memorize in VMLMB"),
    _lst("ftol", NULL, "REAL", OPT_REAL, "Function tolerance for the global convergence"),
    _lst("gtol", NULL, "REAL", OPT_REAL, "Gradient tolerance for the global convergence"),
    _lst("sftol", NULL, "REAL", OPT_REAL, "Function tolerance for the line search"),
    _lst("sgtol", NULL, "REAL", OPT_REAL, "Gradient tolerance for the line search"),
    _lst("sxtol", NULL, "REAL", OPT_REAL, "Step tolerance for the line search"),
    /* MISCELLANEOUS */
    _lst("bitpix", -32, "BITPIX", OPT_INTEGER, "Bits per pixel"),
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
    opt_error, "expecting value and units for " + opt;
  }
  fact = mira_parse_angular_units(units, 0);
  if (fact == 0) {
    opt_error, "invalid units for " + opt;
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
    opt_error, "expecting value and units for " + opt;
  }
  fact = mira_parse_length_units(units, 0);
  if (fact == 0) {
    opt_error, "invalid units for " + opt;
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
  monochromatic = TRUE;
  format = mira_format;

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

  /* Check bitpix. */
  if (opt.bitpix != 8 && opt.bitpix != 16 && opt.bitpix != 32 &&
      opt.bitpix != -32 && opt.bitpix != -64) {
    opt_error, "invalid value for `--bitpix`";
  }

  /* Initial comment. */
  comment = ("Image reconstructed by MiRA algorithm" +
             " <https://github.com/emmt/MiRA>");
  /* Setup the regularization. */
  regul_post = FALSE;
  regul_name = opt.regul;
  if (! is_void(regul_name)) {
    if (min(opt.mu) < 0.0) {
      opt_error, "value(s) of `--mu` must be >= 0.0";
    }
    if (regul_name == "hyperbolic") {
      /* Edge-preserving smoothness prior. */
      if (opt.tau <= 0.0) {
        opt_error, "value of `--tau` must be > 0.0";
      }
      if (min(opt.eta) <= 0.0) {
        opt_error, "values of `--eta` must be > 0.0";
      }
      regul = rgl_new("hyperbolic");
      rgl_config, regul, tau=opt.tau, eta=opt.eta;
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s, TAU=%g, ETA=%s",
               regul_name, format(opt.mu), opt.tau, format(opt.eta));
    } else if (regul_name == "compactness") {
      /* Compactness prior. */
      if (is_void(opt.gamma)) {
        opt_error, "option `--gamma=...` must be specified with \"compactness\" regularization";
      }
      value = _mira_cli_parse_angle(opt.gamma, "`--gamma=...`");
      if (value <= 0) {
        opt_error, "invalid value for `--gamma=...`";
      }
      h_set, opt, gamma=value;
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s and GAMMA=%g",
               regul_name, format(opt.mu), opt.gamma);
      regul_post = TRUE;

#if 0
    } else {
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
#endif
    } else {
      opt_error, "unsupported regularization \""+regul_name+"\"";
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
    pixelsize = _mira_cli_parse_angle(opt.pixelsize, "`--pixelsize=...`");
    if (pixelsize <= 0) {
      opt_error, "invalid value for `--pixelsize=...`";
    }
  }
  if (! is_void(opt.fov)) {
    if (! is_void(opt.dim)) {
      opt_error, "only one of `--fov=...` or `--dim=...` can be specified";
    }
    if (is_void(opt.pixelsize)) {
      opt_error, "option `--pixelsize=...` must be specified with `--fov=...`";
    }
    fov = _mira_cli_parse_angle(opt.fov, "`--fov=...`");
    if (fov <= 0) {
      opt_error, "invalid value for `--fov=...`";
    }
    dim = lround(fov/pixelsize);
  } else if (! is_void(opt.dim)) {
    if (opt.dim <= 0) {
      opt_error, "bad value for `--dim=...`";
    }
    dim = opt.dim;
  } else if (! initial_filename) {
    opt_error, ("image dimension must be specified with `--dim=...` "
                + "or `--fov=..` when no initial image is given");
  }

  /* Initial image. */
  local initial;
  if (initial_filename) {
    /* Read initial image. */
    img = mira_read_image(initial_filename);
    naxis = img.naxis;
    naxis1 = img.naxis1;
    naxis2 = img.naxis2;
    naxis3 = (naxis >= 3 ? img.naxis3 : 1);
    eq_nocopy, initial, img.arr;
    if (monochromatic) {
      if (naxis == 3) {
        if (naxis3 != 1) {
          opt_error, "expecting a 2D initial image";
        }
        initial = initial(,,avg);
      }
    }
    // FIXME: only the pixel size is considered...

    /* Figure out pixel size if not overriden by command line argument. */
    if (is_void(pixelsize)) {
      /* Get pixel size from FITS header and fix the orientation of the
         image. */
      cdelt1 = img.cdelt1;
      cdelt2 = img.cdelt2;
      cunit1 = mira_parse_angular_unit(img.cunit1, 0);
      cunit2 = mira_parse_angular_unit(img.cunit2, 0);
      if (cunit1 == 0 || cunit2 == 0) {
        opt_error, ("Unknown units in CUNIT1 or CUNIT2, "
                    + "can be overriden by `--pixelsize`");
      }
      pixsiz1 = cdelt1*cunit1;
      pixsiz2 = cdelt2*cunit2;
      if (abs(pixsiz1 - pixsiz2) > 1e-7*max(abs(pixsiz1), abs(pixsiz2))) {
        opt_error, "non-square pixels not supported";
      }
      pixelsize = (abs(pixsiz1) + abs(pixsiz2))/2;
      if (pixsiz1 < 0) initial = unref(initial)(::-1,..);
      if (pixsiz2 < 0) initial = unref(initial)(,::-1,..);
    }

    /* Fix dimensions of initial image. */
    if (is_void(dim)) {
      dim = max(img.naxis1, img.naxis2);
    }
    if (dim != img.naxis1 || dim != img.naxis2) {
      initial = mira_extract_region(initial, dim);
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

  /* Parse wavelenght settings. */
  keys = ["effband", "effwave", "wavemin", "wavemax"];
  for (k = 1; k <= numberof(keys); ++k) {
    local str;
    key = keys(k);
    eq_nocopy, str, h_get(opt, key);
    if (! is_void(str)) {
      value = _mira_cli_parse_length(str, "`--"+key+"=...`");
      if (value <= 0) {
        opt_error, "invalid value for `--"+key+"=...`";
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
                    monochromatic = monochromatic,
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
      opt_error, "bad value for option `--bootstrap`";
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
      opt_error, "bad number of values for option `--mu`";
    }
  }
  if (! is_void(opt.maxiter) && opt.maxiter < 0) {
      opt_error, "bad value for option `--maxiter`";
  }
  if (! is_void(opt.maxeval) && opt.maxeval < 0) {
    opt_error, "bad value for option `--maxeval`";
  }

  /* Reconstructions. */
  local final;
  eq_nocopy, final, initial;
  for (k = 1; k <= n; ++k) {
    final = mira_solve(master, final,
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
    if (opt.recenter && k <= opt.bootstrap) {
      final = mira_recenter(final, quiet=opt.quiet);
    }
  }

  /* Save the result. */
  fh = mira_save_image(mira_wrap_image(final, master), final_filename,
                       overwrite=opt.overwrite, bitpix=opt.bitpix,
                       comment=comment);
  if (opt.save_initial) {
    mira_save_image, initial, fh, bitpix=opt.bitpix,
      extname="INITIAL_IMAGE",
      comment="Initial image used by MiRA";
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
