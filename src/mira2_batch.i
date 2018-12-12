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


if (! is_func(nfft_new)) {
  /* Attempt to pre-load YNFFT. */
  include, "nfft.i", 3;
}
_mira_batch_xform = (is_func(nfft_new) ? "nfft" : "separable");
_MIRA_CL_USAGE = "Usage: ymira [-help] [OPTIONS] INPUT [...] OUTPUT\n";
_MIRA_CL_BRIEF = \
  "Run the MiRA image reconstruction algorithm.  INPUT and [...] are the\n"+
  "OI-FITS data file and OUTPUT is the result saved into a FITS file.\n\n"+
  "Environment variables MIRA_SRCDIR and MIRA_YORICK may be set to specify\n"+
  "the directory where are installed the sources and the path to the\n"+
  "Yorick interpreter.";
_MIRA_CL_OPTS = _lst\
  ("\nData selection:",
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
   _lst("use_vis", [], "all|none|amp|phi", OPT_STRING,
        "Complex visibility data to consider"),
   _lst("use_vis2", [], "all|none", OPT_STRING,
        "Powerspectrum data to consider"),
   _lst("use_t3", [], "all|none|amp|phi", OPT_STRING,
        "Bispectrum data to consider"),
   _lst("convexcost", [], "yes|no", OPT_STRING,
        "Use convex approximation for fitting complex data?"),
   _lst("phasecost", "vonmises", "vonmises|haniff|convexlimit", OPT_STRING,
        "Co-log-likelihood approximation to use for phase data"),
   "\nImage parameters and direct model:",
   _lst("pixelsize", [], "ANGLE", OPT_STRING,
        "Angular size of pixels, e.g. 0.1mas"),
   _lst("fov", [], "ANGLE", OPT_STRING,
        "Angular size of the field of view, e.g. 20mas"),
   _lst("imagesize", [], "DIM|DIM1xDIM2", OPT_STRING,
        "Number of pixels per side of the image"),
   _lst("flux", [], "VALUE", OPT_REAL,
        "Total flux (sum of pixels = VALUE)"),
   _lst("fluxerr", [], "VALUE", OPT_REAL,
        "Assumed standard deviation for the total flux"),
   _lst("min", [], "LOWER", OPT_REAL,
        "Lower bound for the pixel values"),
   _lst("max", [], "UPPER", OPT_REAL,
        "Upper bound for the pixel values"),
   _lst("xform", [], "nfft|separable|nonseparable", OPT_STRING,
        "Method to compute the Fourier transform"),
   _lst("smearingfunction", "none", "none|sinc|gauss", OPT_STRING,
        "Method to model the effects of the spectral bandwidth smearing"),
   _lst("smearingfactor", 1.0, "VALUE", OPT_REAL,
        "Factor to scale the effects of the spectral bandwidth smearing"),
   _lst("nthreads", 1, "NUMBER", OPT_INTEGER,
        "Number of threads for the fast Fourier transform"),
   "\nRegularization settings:",
   _lst("regul", [], "NAME", OPT_STRING,
        "Name of regularization method (-regul=help for more information)"),
   _lst("mu", [], "VALUE(S)", OPT_REAL_LIST,
        "Global regularization weight"),
   _lst("tau", [], "VALUE", OPT_REAL,
        "Edge preserving threshold"),
   _lst("eta", [], "VALUE(S)", OPT_REAL_LIST,
        "Gradient scales along dimensions"),
   _lst("gamma", [], "FWHM", OPT_STRING,
        "A priori full half width at half maximum, e.g. 15mas"),
   "\nInitial image:",
   _lst("initial", [], "NAME", OPT_STRING,
        "FITS file or method for initial image"),
   _lst("initialhdu", [], "HDUNAME", OPT_STRING,
        "Name of FITS extension with initial image"),
   _lst("seed", [], "VALUE", OPT_REAL,
        "Seed for the random generator"),
   "\nOutput file:",
   _lst("save_visibilities", [], [], OPT_FLAG,
        "Save model complex visibilities"),
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
   _lst("oi-imaging", [], [], OPT_FLAG,
        "OI-Imaging mode"),
   _lst("plugin", [], "NAME", OPT_STRING,
        "Name of plugin"),
   _lst("debug", [], [], OPT_FLAG,
        "Debug mode"),
   _lst("help", [], [], OPT_HELP,
        "Print out this help"),
   _lst("version", MIRA_VERSION, [], OPT_VERSION,
        "Print out version number"));

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

func mira_get_yesno(opt, key, def)
/* DOCUMENT bool = mira_get_yesno(opt, key);
         or bool = mira_get_yesno(opt, key, def);

     yields true/false depending whether field KEY is "yes" or "no" (or not
     set) in option table OPT; otherwise, throw an error.  If KEY is not set
     DEF is returned.

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
    return !(!def);
  }
  opt_error, "Invalid value for `-"+key+"=...`, expecting `yes` or `no`";
}
errs2caller, mira_get_yesno;

func _mira_get_use_polar(tab, key, use_amp, use_phi, def)
{
  value = h_get(tab, key);
  if (is_void(value))  return def;
  if (value == "none") return 0n;
  if (value == "amp")  return use_amp;
  if (value == "phi")  return use_phi;
  if (value == "all")  return (use_amp|use_phi);
  throw, "invalid value for option `--" + key + "=all|none|amp|phi`";
}
errs2caller, _mira_get_use_polar;

func _mira_get_use_allnone(tab, key, bits, def)
{
  value = h_get(tab, key);
  if (is_void(value))  return def;
  if (value == "none") return 0n;
  if (value == "all")  return bits;
  throw, "invalid value for option `--" + key + "=all|none`";
}
errs2caller, _mira_get_use_allnone;

func _mira_fetch_plugin(&argv, &options)
{
  n = numberof(argv);
  if (n < 1) {
    return;
  }
  test = strglob("--", argv);
  imax = anyof(test) ? where(test)(1) - 1 : n;
  if (imax < 1) {
    return;
  }
  sel = strgrep("^--?plugin=(.*)", argv, sub=1, n=1);
  test = (sel(2,) > 0);
  if (noneof(test)) {
    return;
  }
  i = where(test)(1);
  if (i > imax) {
    return;
  }
  name = strpart(argv(i), sel(,i));
  if (strlen(name) < 1) {
    throw, "Missing plugin name";
  }
  if (strpart(name, -1:0) == ".i") {
    file = name;
    sel = strgrep("^(|.*/)mira2_plugin_(.+)\.i$", file, sub=2, n=1);
    name = strpart(file, sel);
    if (strlen(name) < 1) {
      throw, "Plugin file name must be \"SOMEDIR/mira2_plugin_SOMENAME.i\"";
    }
    if (strglob("*/*", file) && ! strglob("[~/]*", file)) {
      /* Assume file path is relative. */
      file = cd(".") + file;
    }
  } else {
    file = "mira2_plugin_" + name + ".i";
    tmp = MIRA_HOME + file;
    write, tmp;
    if (open(tmp, "r", 1n)) {
      eq_nocopy, file, tmp;
    }
  }
  include, file, 3;
  init = "mira_plugin_" + name + "_init";
  if (! symbol_exists(init)) {
    throw, ("File \"" + file + "\" not readable/found or function \"" +
            init + "\" not defined in this file.");
  }
  plugin = symbol_def(init)();
  if (! is_hash(plugin) || ! h_has(plugin, "__vops__")) {
    throw, ("Invalid initialization of plugin \"" + name + "\"");
  }
  options = _cat(plugin.__vops__.options, options);
  return plugin;
}

func mira_get_fits_use_polar(fh, kwd, use_amp, use_phi)
{
  value = mira_get_fits_string(fh, kwd);
  id = identof(value);
  if (id == Y_VOID) return [];
  if (is_scalar(id)) {
    if (id == Y_STRING) {
      if (value == "NONE") return "none";
      if (value == "AMP")  return "amp";
      if (value == "PHI")  return "phi";
      if (value == "ALL")  return "all";
    } else if (id == Y_CHAR) {
      if (value == 'T') return "all";
      if (value == 'F') return "none";
    }
  }
  throw, "invalid value for FITS keyword " + kwd;
}
errs2caller, mira_get_fits_use_polar;

func mira_batch_read_initial_image(filename, hdu)
{
  fh = fits_open(filename);
  number = mira_parse(long, hdu);
  if (is_void(number)) {
    while (! fits_eof(fh)) {
      value = fits_get(fh, "HDUNAME");
      if (! is_void(value) && strcase(1, strtrim(value, 2)) == hdu) {
        break;
      }
    }
  } else {
    fits_goto_hdu, fh, number;
  }
  if (fits_eof(fh)) {
    error, "initial hdu not found in \""+filename+"\"";
  }
  if (fits_get_xtension(fh) != "IMAGE") {
    error, "initial hdu is not an IMAGE in \""+filename+"\"";
  }
  img = mira_read_image(fh);
  fits_close, fh;
  return img;
}

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

  /* Pre-parse options.  A plugin may insert its own options *before* the list
     of standard options so that a new list can be created without perturbating
     the standard one. */
  arguments = mira_concatenate_arguments(argv);
  options = _MIRA_CL_OPTS; /* to start a new list of options */
  plugin = _mira_fetch_plugin(argv, options);
  options = opt_init(_MIRA_CL_USAGE, _MIRA_CL_BRIEF, options);
  opt = opt_parse(options, argv);
  if (is_void(opt)) {
    /* Options "-help", or "-usage", or "-version" have been set. */
    return;
  }
  h_set, opt, "oi_imaging", h_pop(opt, "oi-imaging");
  h_set, opt, flags = 0n;
  if (is_hash(plugin)) {
    subroutine = plugin.__vops__.parse_options;
    subroutine, plugin, opt;
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
    opt_usage, options;
    return;
  }
  final_filename = argv(0);
  if (! opt.debug) {
    error = opt_error;
  }
  if (! opt.overwrite && open(final_filename, "r", 1)) {
    opt_error, ("output file \""+final_filename+"\" already exists");
  }

  /* In OI-Imaging mode, the initial settings and the data are given by the
     input FITS file.  Command line options have precedence. */
  if (opt.oi_imaging) {
    if (argc != 2) {
      opt_error, "option `-oi-imaging` requires exactly 2 positional arguments";
      return;
    }
    mira_set_defaults, opt, mira_read_input_params(argv(1));
  }

  /* Default settings. */
  mira_set_defaults, opt,
    h_new("min", 0.0, flux=1.0, fluxerr=0.0, xform=_mira_batch_xform,
          mu=1.0, tau=1e-6, eta=1.0, gamma=20*MIRA_MILLIARCSECOND);

  /* Data fidelity functions and data selection. */
  flags = opt.flags;
  flags |= _mira_get_use_polar(opt, "use_vis",
                               MIRA_FIT_VISAMP, MIRA_FIT_VISPHI,
                               (MIRA_FIT_VISAMP|MIRA_FIT_VISPHI));
  flags |= _mira_get_use_allnone(opt, "vis2", MIRA_FIT_VIS2, MIRA_FIT_VIS2);
  flags |= _mira_get_use_polar(opt, "use_t3",
                               MIRA_FIT_T3AMP, MIRA_FIT_T3PHI,
                               (MIRA_FIT_T3AMP|MIRA_FIT_T3PHI));
  if (get_yesno(opt, "convexcost", 1n)) {
    flags |= MIRA_CONVEX_APPROX;
  }
  if (is_void(opt.phasecost) || opt.phasecost == "vonmises") {
    flags |= MIRA_VON_MISES_APPROX;
  } else if (opt.phasecost == "haniff") {
    flags |= MIRA_HANIFF_APPROX;
  } else if (opt.phasecost == "convexlimit") {
    flags |= MIRA_CONVEX_LIMIT;
  } else {
    opt_error, ("Bad value for `--phasecost`, expecting `vonmises`, " +
                "`haniff` or `convexlimit`");
  }

  /* Check bitpix. */
  if (opt.bitpix != 8 && opt.bitpix != 16 && opt.bitpix != 32 &&
      opt.bitpix != -32 && opt.bitpix != -64) {
    opt_error, "Invalid value for `--bitpix`";
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
        opt_error, ("Option `-gamma=...` must be specified with " +
                    "\"compactness\" regularization");
      }
      value = get_angle(opt, "gamma", strictly_positive);
      grow, comment,
        swrite(format="Regularization: \"%s\" with MU=%s and GAMMA=%s",
               regul_name, mira_format(opt.mu), opt.gamma);
      h_set, opt, gamma=value;
      regul_post = TRUE;
    } else {
      opt_error, "Unsupported regularization \""+regul_name+"\"";
    }
  }

  /* Initial image. */
  if (is_void(opt.initial)) {
    if (! opt.oi-imaging || is_void(opt.initialhdu)) {
      opt_error, ("An initial image must be specified, e.g. with " +
                  "`--initial=Dirac|random|FILENAME`");
    }
    h_set, opt, initial = mira_batch_read_initial_image(argv(1),
                                                        opt.initialhdu);
  }
  if (is_string(opt.initial)) {
    if (opt.initial != "random" && opt.initial != "Dirac") {
      /* Read initial image. */
      h_set, opt, initial = mira_read_image(opt.initial);
    }
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
  } else if (is_string(opt.initial)) {
    opt_error, ("When the initial image is not from a file, image " +
                "dimensions must be specified with exactly two of " +
                "`-fov=...`, `-pixelsize=...` or `-imagesize=...`");
  }

  /* Check some options. */
  if (! is_void(opt.threshold) && (opt.threshold < 0 || opt.threshold >= 1)) {
    opt_error, "Invalid value for `-threshold=...`";
  }

  /* Fix/generate initial image. */
  if (is_string(opt.initial)) {
    image = h_new(naxis = 2,
                  naxis1 = dim1,
                  crpix1 = (dim1/2) + 1,
                  crval1 = 0.0,
                  ctype1 = "RA---TAN",
                  cdelt1 = pixelsize/MIRA_MILLIARCSECOND,
                  cunit1 = "mas",
                  naxis2 = dim2,
                  crpix2 = (dim2/2) + 1,
                  crval2 = 0.0,
                  ctype2 = "DEC--TAN",
                  cdelt2 = pixelsize/MIRA_MILLIARCSECOND,
                  cunit2 = "mas");
    if (opt.initial == "random") {
      if (! is_void(opt.seed)) {
        if (opt.seed <= 0.0 || opt.seed >= 1.0) {
          opt_error, "Seed value must be between 0.0 and 1.0 non-inclusive";
        }
        random_seed, opt.seed;
      }
      h_set, image, arr = random(dims);
    } else if (opt.initial == "Dirac") {
      h_set, image, arr = array(double, dim1, dim2);
      image.arr(image.crpix1, image.crpix2) = 1.0;
    } else {
      opt_error, "Option -initial=\""+opt.initial+"\" not yet implemented";
    }
    h_set, opt, initial = image;
  }

  /* Fix dimensions, sampling and normalization of initial image. */
  image = opt.initial;
  naxis = image.naxis;
  if (naxis < 2 || naxis > 3) {
    opt_error, "Expecting a 2D/3D initial image";
  }
  naxis1 = image.naxis1;
  naxis2 = image.naxis2;
  naxis3 = (naxis >= 3 ? image.naxis3 : 1);
  if (naxis == 3) {
    if (naxis3 != 1) {
      warn, "Converting 3D image into a 2D grayscaled image";
    }
    h_set, image, arr = image.arr(,,avg), naxis=2;
    h_pop, image, "naxis3";
    h_pop, image, "crpix3";
    h_pop, image, "crval3";
    h_pop, image, "cdelt3";
    h_pop, image, "cunit3";
    h_pop, image, "ctype3";
    naxis = 2;
  }
  siunit = "radian"; // pixelsize and FOV are in SI units
  cdelt1 = mira_convert_units(image.cunit1, siunit)*image.cdelt1;
  cdelt2 = mira_convert_units(image.cunit2, siunit)*image.cdelt2;
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
  image = mira_resample_image(image, pad=0, norm=1n,
                              naxis1=dim1,   naxis2=dim2,
                              crpix1=crpix1, crpix2=crpix2,
                              crval1=crval,  crval2=crval,
                              cdelt1=cdelt,  cdelt2=cdelt,
                              cunit1=cunit,  cunit2=cunit);
  h_set, opt, initial = image;

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
    wavemin = effwave - effband/2.0;
    wavemax = effwave + effband/2.0;
  } else if (choice == 12) {
    if (wavemin > wavemax) {
      opt_error, "Incompatible values for `-wavemin=...` and `-wavemax=...` ";
    }
  } else if (choice != 0 && choice != 4 && choice != 8) {
    opt_error, ("Specify `-effwave` and `-effband`, "+
                "or `-wavemin` and/or `-wavemax`");
  }
  h_pop, opt, "effwave";
  h_pop, opt, "effband";
  h_set, opt, wavemin=wavemin, wavemax=wavemax;

  /* Read input data. */
  master = mira_new(argv(1:-1),
                    plugin = plugin,
                    target = opt.target,
                    wavemin = opt.wavemin,
                    wavemax = opt.wavemax,
                    quiet = opt.quiet,
                    pixelsize = pixelsize,
                    dims = dims,
                    flags = flags,
                    xform = opt.xform,
                    nthreads = opt.nthreads,
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

  /* Run image reconstruction stages. */
  local initial_arr, final_arr;
  eq_nocopy, initial_arr, image.arr;
  current_arr = image.arr; /* copy */
  for (k = 1; k <= n; ++k) {
    if (! is_void(opt.threshold) && opt.threshold > 0) {
      current_arr = mira_soft_threshold(current_arr, opt.threshold, opt.flux);
    }
    current_arr = mira_solve(master, current_arr,
                             maxeval = opt.maxeval,
                             maxiter = opt.maxiter,
                             verb = opt.verb,
                             xmin = opt("min"),
                             xmax = opt("max"),
                             flux = opt.flux,
                             fluxerr = opt.fluxerr,
                             regul = regul,
                             mu = mu(k),
                             mem = opt.mem,
                             ftol = opt.ftol,
                             gtol = opt.gtol,
                             sftol = opt.sftol,
                             sgtol = opt.sgtol,
                             sxtol = opt.sxtol);
    if (opt.recenter && k < n) {
      current_arr = mira_recenter(current_arr, quiet=opt.quiet);
    }
  }
  eq_nocopy, final_arr, current_arr;

  /* Save the result. */
  fh = mira_save_image(h_set(image, arr = final_arr), final_filename,
                       overwrite=opt.overwrite, bitpix=opt.bitpix,
                       hduname="IMAGE-OI FINAL IMAGE", comment=comment);
  h_set, opt, init_img = (opt.save_initial ? "IMAGE-OI INITIAL IMAGE" : []);
  mira_write_input_params, fh, master, opt;
  if (opt.save_initial) {
    mira_save_image, h_set(image, arr = initial_arr), fh,
      bitpix=opt.bitpix, hduname="IMAGE-OI INITIAL IMAGE",
      comment="Initial image used by MiRA";
  }
  if (opt.save_visibilities) {
    inform, "Saving model complex visibilities...";
    mira_update, master, final_arr;
    mira_save_visibilities, master, fh;
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
