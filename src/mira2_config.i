/*
 * mira2_config.i -
 *
 * Management of configuration and options in MiRA.
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

if (! is_scalar(MIRA_HOME) || ! is_string(MIRA_HOME)) {
  error, "include \"mira2.i\" first";
}

/*
 * Configuration
 * =============
 *
 * Configuration of MiRA is done in several steps:
 *
 *  - Stage 0: merging of OI-FITS data for a given target.
 *
 *  - Stage 1: selection of individual measurements.  This depends on
 *    parameters such as `wavemin`, `wavemax`, `freqmin`, `freqmax`, etc.  The
 *    list of observed baselines (or frequencies) is built during this stage.
 *    Changing anything here also result in rebuilding the linear operator (see
 *    below).
 *
 *  - Stage 2: reduction of the number of coordinates depending on whether the
 *    model of the linear transform between the image and the observed
 *    frequencies depends on (u/λ,v/λ), on (u,v,λ), or on (u,v,λ,Δλ).
 *
 *  - Stage 3: building the `xform` operator which computes the nonuniform
 *    Fourier transform of the image.  This depends on parameters such as the
 *    image dimensions, the pixel size and the specific model for the transform
 *    to use.
 */
func mira_new(.., target=, wavemin=, wavemax=, flags=, pixelsize=, dims=,
              xform=, nthreads=, smearingfunction=, smearingfactor=,
              atol=, rtol=, quiet=, plugin=,
              noise_method=, noise_level=, baseline_precision=)
/* DOCUMENT obj = mira_new(oidata, ..., target=...);

     Create a new MiRA instance for fitting the data in the `oidata`,
     ... arguments.

     When a new instance is created, its default parameters are:

       - minimal wavelength: wavemin = 0
       - maximal wavelength: wavemax = +Inf
       - angular pixel size: pixelsize = 5 milliarcseconds
       - image dimensions: dims = [2,128,128]
       - model of the nonuniform Fourier transform: xform = "nfft"

     these can be changed by the corresponding keywords (see `mira_config`).

     Keywords ATOL and RTOL can be used to specify the absolute and relative
     tolerances when comparing numerical values which should be identical
     (i.e., when merging different OI-FITS files).


   SEE ALSO mira_config.
*/
{
  /* Get precision (in meters) for rounding baselines. */
  if (! scalar_double(baseline_precision, 1e-5) || baseline_precision <= 0) {
    error, "absolute baseline precision must be a strictly positive length";
  }

  /* Create an empty MiRA instance. */
  if (! scalar_string(target, "*")) {
    error, "target must be a string";
  }
  master = h_new(oidata = oifits_new(), target = target,
                 stage = 0, xform="nfft", nthreads = 1,
                 smearingfactor = 1.0, smearingfunction = "none",
                 dims=[2,128,128], pixelsize = 5*MIRA_MILLIARCSECOND,
                 wavemin = 0.0, wavemax = MIRA_HUGE,
                 first = [],
                 model = h_new(),
                 baseline_precision = baseline_precision,
                 flags = (MIRA_FIT_VISAMP | MIRA_FIT_VISPHI | MIRA_FIT_VIS2 |
                          MIRA_FIT_T3AMP | MIRA_FIT_T3PHI |
                          MIRA_CONVEX_APPROX | MIRA_VON_MISES_APPROX));

  /* Apply configuration options. */
  mira_config, master, wavemin=wavemin, wavemax=wavemax, flags=flags,
    pixelsize=pixelsize, dims=dims, xform=xform, nthreads=nthreads,
    smearingfunction=smearingfunction, smearingfactor=smearingfactor,
    plugin=plugin;

  /* Load OI-FITS data file(s). */
  while (more_args()) {
    mira_add_oidata, master, next_arg(), quiet=quiet,
      noise_method=noise_method, noise_level=noise_level,
      atol=atol, rtol=rtol;
  }

  return master;
}

/*---------------------------------------------------------------------------*/
/* CONFIGURING AN INSTANCE */

func mira_config(master, wavemin=, wavemax=, dims=, pixelsize=, xform=,
                 nthreads=, flags=, smearingfunction=, smearingfactor=,
                 plugin=)
/* DOCUMENT mira_config, master, key=val, ...;

     Configure options in MiRA instance `master`.  All options are passed as
     keywords.

     Keywords `wavemin` and `wavemax` specify the wavelength range of the data
     to take into account.  Unless units are specified, these values are in SI
     units (i.e., meters).

     Keyword `pixelsize` specifies the angular pixel size of the restored
     image.  Unless units are specified, this value is in SI units (i.e.,
     radians).

     Keyword `xform` is the name of the model to use for computing the
     nonuniform Fourier transform.  Possibilities are: "simple" or "nfft".

     If keyword `xform` is set to "simple", keywords `smearingfactor` and
     `smearingfunction` can be set to tune the accounting of bandwidth
     smearing.  By default, `smearingfactor=1` and `smearingfunction=sinc`.
     Setting `smearingfactor=0` or `xform="nfft"`, disables the accounting of
     bandwidth smearing.

     Keyword `nthreads` can be specified with the number of threads to
     use for computing the fast Fourier transform.

     Keyword `dims` specifies the dimension(s) of the restored image.  If it is
     a scalar, it specifies the width and height of the image; otherwise,
     `dims` can be `[width,height]` or `[2,width,height]`.

     Keyword `plugin` specifies an optional plugin.  The value is a hash table
     initialized by `mira_new_plugin`.


   SEE ALSO mira_new, mira_new_plugin.
*/
{
  /* Check all options. */
  if (! mira_length(wavemin, master.wavemin) || wavemin < 0.0) {
    throw, "minimum wavelength must be a nonnegative length";
  }
  if (! mira_length(wavemax, master.wavemax) || wavemax < 0.0) {
    throw, "maximum wavelength must be a nonnegative length";
  }
  if (wavemin > wavemax) {
    throw, "minimum wavelength > maximum wavelength";
  }
  if (! mira_angle(pixelsize, master.pixelsize) || pixelsize <= 0.0) {
    throw, "pixel size must be a strictly positive angle";
  }
  master_xform = mira_xform_name(master);
  if (! scalar_string(xform, master_xform)) {
    throw, "`xform` must be a string";
  }
  if (xform != "nfft" && xform != "nonseparable" && xform != "separable") {
    throw, ("unknown `xform` name (must be one of \"nfft\", \"separable\" " +
            "or \"nonseparable\"");
  }
  if (! scalar_long(nthreads, master.nthreads) || nthreads < 1) {
    throw, "`nthreads` must be a strictly nonnegative integer";
  }
  if (! scalar_double(smearingfactor, master.smearingfactor) ||
      smearingfactor < 0) {
    throw, "`smearingfactor` must have a nonnegative value";
  }
  if (! scalar_string(smearingfunction, master.smearingfunction)) {
    throw, "`smearingfunction` must be a string";
  }
  if (smearingfunction != "sinc" && smearingfunction != "gauss" &&
      smearingfunction != "none") {
    throw, ("unknown `smearingfunction` name (must be one of \"none\", " +
            "\"sinc\" or \"gauss\"");
  }
  bad_dims = 1n;
  if (is_void(dims)) {
    dims = mira_image_size(master);
  }
  if (is_integer(dims)) {
    if (is_scalar(dims) && dims > 0) {
      dims = [2, long(dims), long(dims)];
      bad_dims = 0n;
    } else if (is_vector(dims) && min(dims) > 0) {
      if (numberof(dims) == 2) {
        dims = grow(2, long(dims));
        bad_dims = 0n;
      } else if (numberof(dims) == 3 && dims(1) == 2) {
        dims = long(dims);
        bad_dims = 0n;
      }
    }
  }
  if (bad_dims) {
    throw, "invalid image dimension(s)";
  }
  if (is_void(flags)) {
    flags = master.flags;
  }
  flags = mira_fix_flags(flags);

  /* Apply changes if any. */
  changes = 0;
  if (is_hash(plugin)) {
    // FIXME: Find means to unload plugin.
    h_set, master, plugin = plugin;
  }
  if (master.wavemin != wavemin || master.wavemax != wavemax) {
    h_set, master, model = h_new(), stage = min(master.stage, 0),
      wavemin = wavemin, wavemax = wavemax;
    changes |= 2;
  }
  if (master.flags != flags) {
    h_set, master, model = h_new(), stage = min(master.stage, 0),
      flags = flags;
  }
  if (master.pixelsize != pixelsize) {
    h_set, master, model = h_new(), stage = min(master.stage, 2),
      pixelsize = pixelsize;
    changes |= 1;
  }
  if (xform != master_xform || nthreads != master.nthreads ||
      smearingfunction != master.smearingfunction ||
      smearingfactor != master.smearingfactor) {
    h_set, master, model = h_new(), stage = min(master.stage, 2),
      xform = xform, nthreads = nthreads, smearingfunction = smearingfunction,
      smearingfactor = smearingfactor;
  }
  if (! mira_same_dimensions(master.dims, dims)) {
    h_set, master, model = h_new(), stage = min(master.stage, 2),
      dims = dims;
    changes |= 1;
  }

  if ((changes&1) != 0 ||
      ! h_has(master, "img_x") ||  ! h_has(master, "img_y")) {
    pixelsize = mira_pixel_size(master);
    nx = mira_image_size(master, 1);
    ny = mira_image_size(master, 2);
    h_set, master,
      img_x = mira_sky_coordinates(nx, pixelsize),
      img_y = mira_sky_coordinates(ny, pixelsize);
  }
  if ((changes&2) != 0 || ! h_has(master, "img_wave")) {
    // FIXME: should use mean/median data wavelength.
    h_set, master, img_wave = (master.wavemin + master.wavemax)/2.0;
  }

  return master;
}

func mira_xform_name(master)
/* DOCUMENT mira_xform_name(master);

     yields the name of the nonuniform Fourier transform model used in MiRA
     instance `master`.  The nonuniform Fourier transform model can be changed
     with `mira_config`.

   SEE ALSO: mira_config.
 */
{
  return (is_string(master.xform) ? master.xform : master.xform.name);
}

local mira_pixel_size, mira_maximum_pixel_size;
/* DOCUMENT mira_pixel_size(master);
         or mira_maximum_pixel_size(master);
         or mira_maximum_pixel_size, master;

     The first function yields the angular size (in radians) of the pixels for
     the image restored with MiRA instance `master`.  The image pixel size can
     be changed with `mira_config`.

     The second function yields the maximum pixel size allowed to avoid
     aliasing.  The returned value depends on the current set of selected data
     and image dimensions.  When called as a subroutine, the maximum pixel size
     is nicely printed.

   SEE ALSO: mira_config.
 */

func mira_pixel_size(master) /* DOCUMENTED */
{
  return master.pixelsize;
}


local mira_image_x, mira_image_x, mira_image_wave;
/* DOCUMENT mira_image_x(master);
         or mira_image_x(master);
         or mira_image_wave(master);

     The two first functions yield the relative right ascension (RA) and declination (DEC)
     (in radians) of the pixels for the image restored with MiRA instance
     `master`.  The image pixel size and dimensions can be changed with
     `mira_config`.

     The third function yields the wavelength(s) of the restored image along
     its 3rd axis.

   SEE ALSO: mira_config, mira_image_size, mira_pixel_size.
 */

func mira_image_x(master) /* DOCUMENTED */
{
  return master.img_x;
}

func mira_image_y(master) /* DOCUMENTED */
{
  return master.img_y;
}

func mira_image_wave(master) /* DOCUMENTED */
{
  return master.img_wave;
}

func mira_maximum_pixel_size(master) /* DOCUMENTED */
{
  /* Make sure to apply data selection. */
  if (master.stage < 1) {
    _mira_select_data, master;
  }

  local ufreq, vfreq;
  coords = master.coords;
  if (h_has(coords, "ufreq")) {
    eq_nocopy, ufreq, coords.ufreq;
    eq_nocopy, vfreq, coords.vfreq;
  } else {
    ufreq = coords.u/coords.wave;
    vfreq = coords.v/coords.wave;
  }
  theta = 0.5/max(max(abs(ufreq)), max(abs(vfreq)));
  if (am_subroutine()) {
    inform, "Maximum pixel size = %g rad = %g mas\n",
      theta, theta/MIRA_MILLIARCSECOND;
  }
  return theta;
}

func mira_image_size(master, i)
/* DOCUMENT mira_image_size(master);
         or mira_image_size(master, i);

     yields the dimensions of the model or of its `i`-th dimension.  The image
     size can be changed with `mira_config`.

   SEE ALSO: mira_config, dimsof.
 */
{
  return is_void(i) ? master.dims : master.dims(i+1);
}

local mira_smearing_factor, mira_smearing_function;
/* DOCUMENT mira_smearing_factor(master, 0/1);
         or mira_smearing_function(master, 0/1);

     yield the value of the smearing factor or the name of the smearing
     function in MiRA instance MASTER.  If 2nd argument is true, the returned
     value is that of the image to complex visibilities transform.

   SEE ALSO: mira_update, mira_config.
 */

func mira_smearing_factor(master, xform) /* DOCUMENTED */
{
  if (xform) {
    mira_update, master;
    return master.xform.smearingfactor;
  } else {
    return master.smearingfactor;
  }
}

func mira_smearing_function(master, xform) /* DOCUMENTED */
{
  if (xform) {
    mira_update, master;
    return master.xform.smearingfunction;
  } else {
    return master.smearingfunction;
  }
}

local mira_minimal_wavelength, mira_maximal_wavelength;
/* DOCUMENT mira_minimal_wavelength(master);
         or mira_maximal_wavelength(master);

     yield the minimal and maximal wavelengths (in meters) for selected data in
     MiRA instance `master`.  The range of selected wavelengths can be changed
     with `mira_config`.

   SEE ALSO: mira_config, mira_image_size, mira_xform_name, mira_pixel_size.
 */
func mira_minimal_wavelength(master) { return master.wavemin; }
func mira_maximal_wavelength(master) { return master.wavemax; }

local MIRA_FIT_VISAMP, MIRA_FIT_VISPHI, MIRA_FIT_VIS2;
local MIRA_FIT_T3AMP, MIRA_FIT_T3PHI;
local MIRA_CONVEX_APPROX;
local MIRA_HANIFF_APPROX, MIRA_CONVEX_LIMIT, MIRA_VON_MISES_APPROX;
local MIRA_KEEP_WAVELENGTH, MIRA_KEEP_BANDWIDTH;
local mira_fix_flags, mira_format_flags;
func mira_flags(master)
/* DOCUMENT mira_flags(master);
         or mira_fix_flags(flags);
         or mira_format_flags(flags);

     The first function yields the current flags used by MiRA instance
     `master`.  The flags may be changed with `mira_config` and the following
     bits can be combined into the flags:

         MIRA_FIT_VISAMP ......... fit amplitude of complex visibilities;
         MIRA_FIT_VISPHI ......... fit phase of complex visibilities;
         MIRA_FIT_VIS2 ........... fit powerspectrum data;
         MIRA_FIT_T3AMP .......... fit amplitude of bispectrum data;
         MIRA_FIT_T3PHI .......... fit phase of bispectrum data;
         MIRA_CONVEX_APPROX ...... use convex approximation for fitting polar
                                   data;
         MIRA_HANIFF_APPROX ...... use Haniff approximation for phase data;
         MIRA_CONVEX_LIMIT ....... use limit of the convex approximation for
                                   fitting phase data;
         MIRA_VON_MISES_APPROX ... use von Mises approximation for phase data;
         MIRA_KEEP_WAVELENGTH .... keep wavelength coordinate;
         MIRA_KEEP_BANDWIDTH ..... keep bandwidth coordinate;

      The second function returns FLAGS after checking its type, bits and
      setting defaults.

      The third function yields a textual representation of FLAGS which is
      assumed to be a bitwise combination of options in MiRA.


   SEE ALSO: mira_config.
 */
{
  return master.flags;
}

MIRA_FIT_VISAMP       = (1n << 0);
MIRA_FIT_VISPHI       = (1n << 1);
MIRA_FIT_VIS2         = (1n << 2);
MIRA_FIT_T3AMP        = (1n << 3);
MIRA_FIT_T3PHI        = (1n << 4);
MIRA_CONVEX_APPROX    = (1n << 5);

/* bits for phase-only approximation of the objective function */
MIRA_HANIFF_APPROX    = (1n << 6);
MIRA_CONVEX_LIMIT     = (2n << 6);
MIRA_VON_MISES_APPROX = (3n << 6);

MIRA_KEEP_WAVELENGTH  = (1n << 8);
MIRA_KEEP_BANDWIDTH   = (1n << 9);

MIRA_FIT_VIS = (MIRA_FIT_VISAMP | MIRA_FIT_VISPHI);
MIRA_FIT_T3  = (MIRA_FIT_T3AMP | MIRA_FIT_T3PHI);

_MIRA_PHASE_ONLY_BITS = (MIRA_HANIFF_APPROX |
                         MIRA_CONVEX_LIMIT |
                         MIRA_VON_MISES_APPROX);

func mira_fix_flags(flags) /* DOCUMENTED */
{
  if (is_void(flags)) {
    flags = 0n;
  }
  if (! is_scalar(flags) || ! is_integer(flags) || int(flags) != flags) {
    throw, "invalid flags";
  }

  /* Clear unused bits. */
  flags &= (MIRA_FIT_VISAMP | MIRA_FIT_VISPHI | MIRA_FIT_VIS2 |
            MIRA_FIT_T3AMP | MIRA_FIT_T3PHI | MIRA_CONVEX_APPROX |
            _MIRA_PHASE_ONLY_BITS |
            MIRA_KEEP_WAVELENGTH | MIRA_KEEP_WAVELENGTH);

  /* Check phase-only data bits. */
  bits = (flags & _MIRA_PHASE_ONLY_BITS);
  if (bits == 0) {
    /* Von Mises is the default approximation for phase-only data. */
    flags |= MIRA_VON_MISES_APPROX;
  } else if (bits != MIRA_HANIFF_APPROX &&
             bits != MIRA_CONVEX_LIMIT &&
             bits != MIRA_VON_MISES_APPROX) {
    throw, "only a single approximation can be used for phase-only data";
  }

  return int(flags);
}
errs2caller, mira_fix_flags;

func mira_format_flags(flags) /* DOCUMENTED */
{
  str = [];
  _mira_format_flags_helper1, MIRA_FIT_VISAMP,    "MIRA_FIT_VISAMP";
  _mira_format_flags_helper1, MIRA_FIT_VISPHI,    "MIRA_FIT_VISPHI";
  _mira_format_flags_helper1, MIRA_FIT_VIS2,      "MIRA_FIT_VIS2";
  _mira_format_flags_helper1, MIRA_FIT_T3AMP,     "MIRA_FIT_T3AMP";
  _mira_format_flags_helper1, MIRA_FIT_T3PHI,     "MIRA_FIT_T3PHI";
  _mira_format_flags_helper1, MIRA_CONVEX_APPROX, "MIRA_CONVEX_APPROX";
  bits = (flags & _MIRA_PHASE_ONLY_BITS);
  if (bits == MIRA_VON_MISES_APPROX) {
    _mira_format_flags_helper2, "MIRA_VON_MISES_APPROX";
  } else if (bits == MIRA_HANIFF_APPROX) {
    _mira_format_flags_helper2, "MIRA_HANIFF_APPROX";
  } else if (bits == MIRA_CONVEX_LIMIT) {
    _mira_format_flags_helper2, "MIRA_CONVEX_LIMIT";
  }
  _mira_format_flags_helper1, MIRA_KEEP_WAVELENGTH, "MIRA_KEEP_WAVELENGTH";
  _mira_format_flags_helper1, MIRA_KEEP_BANDWIDTH, "MIRA_KEEP_BANDWIDTH";
  return is_void(str) ? "0" : str;
}

func _mira_format_flags_helper1(bits, literal)
{
  extern flags;
  if ((flags & bits) == bits) {
    _mira_format_flags_helper2, literal;
  }
}

func _mira_format_flags_helper2(literal)
{
  extern str;
  if (is_void(str)) {
    str = literal;
  } else {
    str += "|" + literal;
  }
}

/*---------------------------------------------------------------------------*/
/* UPDATING RESSOURCES AND MODEL */

func mira_update(master, img, force=)
/* DOCUMENT mira_update, master;
         or mira_update, master, img;

     Update the internal ressources of MiRA instance `master` and, if image
     `img` is provided, update the current model of the complex visibilities.

     When called as a function, `master` is returned.

     Set keyword `force` true to force updating form the beginning (data
     selection and building of the model operator).


   SEE ALSO: mira_cost, mira_cost_and_gradient, mira_model_vis, mira_model_vis2,
             mira_model_re, mira_model_im, mira_model_amp, mira_model_phi.
 */
{
  if (master.tweaking_visibilities) {
    throw, "mira_update cannot be called while tweaking visibilities";
  }
  if (force) {
    h_set, master, stage = 0;
  }
  if (master.stage < 1) {
    _mira_select_data, master;
  }
  if (master.stage < 3) {
    _mira_define_xform, master;
  }
  if (! is_void(img)) {
    /* Compute complex visibilities (derived quantities will be computed on
       demand). */
    vis = master.xform(img);
    plugin = mira_plugin(master);
    if (is_hash(plugin)) {
      h_set, master, tweaking_visibilities=1n;
      vis = plugin.__vops__.tweak_visibilities(master, vis);
      h_set, master, tweaking_visibilities=0n;
    }
    h_set, master.model, img=img, vis=vis, re=[], im=[],
      vis2=[], amp=[], phi=[], reciprocal_amp=[], reciprocal_vis2=[];
  }
  return master;
}

local mira_model_vis, mira_model_vis2;
local mira_model_re, mira_model_im;
local mira_model_amp, mira_model_phi;
/* DOCUMENT mira_model_vis(master [, idx]);
         or mira_model_vis2(master [, idx]);
         or mira_model_re(master [, idx]);
         or mira_model_im(master [, idx]);
         or mira_model_amp(master [, idx]);
         or mira_model_phi(master [, idx]);

     These functions yield the current model of (respectively) the complex
     visibilities, the powerspectrum, the real and imaginary parts of the
     complex visibilities, the amplitude and phase of the complex visibilities.
     If optional argument `idx` is provided, only a subset of these quantities
     is returned.

   SEE ALSO: mira_update, mira_cost_and_gradient, mira_cost.
 */
func mira_model_vis(master, idx)
{
  model = master.model;
  if (is_void(model.vis)) {
    error, "no complex visibilities have been computed";
  }
  return (is_void(idx) ? model.vis : model.vis(,idx));
}

func mira_model_re(master, idx)
{
  model = master.model;
  if (is_void(model.re)) {
    h_set, model, re = mira_model_vis(master)(1,..);
  }
  return model.re(idx);
}

func mira_model_im(master, idx)
{
  model = master.model;
  if (is_void(model.im)) {
    h_set, model, im = mira_model_vis(master)(2,..);
  }
  return model.im(idx);
}

func mira_model_vis2(master, idx)
{
  model = master.model;
  if (is_void(model.vis2)) {
    local re, im;
    eq_nocopy, re, mira_model_re(master);
    eq_nocopy, im, mira_model_im(master);
    h_set, model, vis2 = re*re + im*im;
  }
  return model.vis2(idx);
}

func mira_model_amp(master, idx)
{
  model = master.model;
  if (is_void(model.amp)) {
    h_set, model, amp = sqrt(mira_model_vis2(master));
  }
  return model.amp(idx);
}

func mira_model_phi(master, idx)
{
  model = master.model;
  if (is_void(model.phi)) {
    h_set, model, phi = mira_atan(mira_model_im(master),
                                  mira_model_re(master));
  }
  return model.phi(idx);
}

func mira_model_reciprocal_amp(master, idx)
{
  model = master.model;
  if (is_void(model.reciprocal_amp)) {
    h_set, model, reciprocal_amp = mira_reciprocal(mira_model_amp(master));
  }
  return model.reciprocal_amp(idx);
}

func mira_model_reciprocal_vis2(master, idx)
{
  model = master.model;
  if (is_void(model.reciprocal_vis2)) {
    h_set, model, reciprocal_vis2 = mira_reciprocal(mira_model_vis2(master));
  }
  return model.reciprocal_vis2(idx);
}

local mira_model_ufreq, mira_model_vfreq;
local mira_model_u, mira_model_v;
/* DOCUMENT mira_model_ufreq(master);
         or mira_model_vfreq(master);
         or mira_model_u(master);
         or mira_model_v(master);
         or mira_model_wave(master);
         or mira_model_band(master);

     These functions yield the coordinates where the model complex visibilities
     are computed.  UFREQ and VFREQ are the spatial frequencies, U and V are
     the baseline coordinates (in meters), WAVE and BAND are the wavelength and
     spectral bandwidth (in meters).

     Depending on the coordinate selection mode, not all coordinates may be
     available: UFREQ and VFREQ are always available, U, V and WAVE are only
     available if `MIRA_KEEP_WAVELENGTH` is set in the "flags" of master (which
     can be set by `mira_config` or `mira_new`) or if spectral bandwidth
     smearing is taken into account, BAND is only available if
     `MIRA_KEEP_BANDWIDTH` is set in the "flags" of master or if spectral
     bandwidth smearing is taken into account.  If the requested coordinate is
     not available, void is returned.

   SEE ALSO: mira_config, mira_update, mira_cost_and_gradient, mira_cost.
 */
func mira_model_ufreq(master) { return _mira_model_coords(master).ufreq; }
func mira_model_vfreq(master) { return _mira_model_coords(master).vfreq; }
func mira_model_u(master)     { return _mira_model_coords(master).u;     }
func mira_model_v(master)     { return _mira_model_coords(master).v;     }
func mira_model_wave(master)  { return _mira_model_coords(master).wave;  }
func mira_model_band(master)  { return _mira_model_coords(master).band;  }

func _mira_model_coords(master)
{
  if (master.stage < 3) {
    mira_update, master;
  }
  return master.coords;
}

/*---------------------------------------------------------------------------*/

func mira_new_plugin(nil, options=, parse_options=,
                     tweak_visibilities=, tweak_gradient=,
                     add_keywords=, add_extensions=)
/* DOCUMENT plugin = mira_new_plugin(...)

     Create a new MiRA plugin.  All settings are passed by keywords.

   KEYWORDS

     OPTIONS: Additional command line options for the plugin.  This should be a
         list of lists.  See the documentation of `opt_init` for the syntax.

     PARSE_OPTIONS: Function called by the command line version of MiRA to let
         the pluging examines the command line options and tunes its behavior
         accordingly.  The function is called with 2 arguments, the plugin
         instance and a hash table with parsed options, and shall return
         nothing.  If not specified, the default is do do nothing.

     TWEAK_VISIBILITIES: Function called to tweak the model complex
         visibilities computed from the pixelized image.  The function is
         called with 2 arguments `master`, the MiRA instance, and `vis`, the
         complex visibilities computed from the image, and shall return the
         modified complex visbilities.  If not specified, the default is do do
         nothing.

     TWEAK_GRADIENT: Function called to tweak the gradient of the objective
         function with respect to the complex visibilities computed from the
         model complex visibilities.  The function is called with 2 arguments
         `master`, the MiRA instance, and `grd`, the gradient of the objective
         function, and shall return the modified complex gradient.  If not
         specified, the default is do do nothing.

     ADD_KEYWORDS: Function called to add keywords in the primary HDU of the
         output FITS file.  The function is called with 2 arguments `master`,
         the MiRA instance, and `fh` a FITS handle.  If not specified, the
         default is do do nothing.

     ADD_EXTENSIONS: Function called to add extension(s) to the output FITS
         file.  The function is called with 2 arguments `master`, the MiRA
         instance, and `fh` a FITS handle.  If not specified, the default is do
         do nothing.

   SEE ALSO: opt_init, mira_config.
 */
{
  if (! is_void(nil)) {
    throw, "All settings are passed by keywords";
  }
  vops = h_new();
  if (is_list(options)) {
    h_set, vops, options=options;
  } else if (! is_void(options)) {
    throw, "keyword `options` must be set with a list of lists";
  }
  h_set, vops,
    parse_options = (is_void(parse_options) ?
                       symlink_to_name("_mira_plugin_parse_options") :
                       parse_options),
    tweak_visibilities = (is_void(tweak_visibilities) ?
                          symlink_to_name("_mira_plugin_tweak_visibilities") :
                          tweak_visibilities),
    tweak_gradient = (is_void(tweak_gradient) ?
                      symlink_to_name("_mira_plugin_tweak_gradient") :
                      tweak_gradient),
    add_keywords = (is_void(add_keywords) ?
                    symlink_to_name("_mira_plugin_add_keywords") :
                    add_keywords),
    add_extensions = (is_void(add_extensions) ?
                      symlink_to_name("_mira_plugin_add_extensions") :
                      add_extensions);
  return h_new(__vops__=vops);
}

func _mira_plugin_parse_options(plugin, opt) { /* do nothing */ }
func _mira_plugin_tweak_visibilities(master, vis) { return vis; }
func _mira_plugin_tweak_gradient(master, grd) { return grd; }
func _mira_plugin_add_extensions(master, fh) { /* do nothing */ }
func _mira_plugin_add_keywords(master, fh) { /* do nothing */ }

func mira_plugin(master)
{
  return master.plugin;
}
