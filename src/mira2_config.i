/*
 * mira2_config.i -
 *
 * Management of configuration and options in MiRA.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2001-2018, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 * This file is part of MiRA: a Multi-aperture Image Reconstruction
 * Algorithm (version 2).
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

if (! is_scalar(MIRA_SRCDIR) || ! is_string(MIRA_SRCDIR)) {
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
 *  - Stage 1: selecting individual measurements.  This depends on parameters
 *    such as `wavemin`, `wavemax`, `freqmin`, `freqmax`, etc.  The list of
 *    observed baselines (or frequencies) is built during this stage.  Changing
 *    anything here also result in rebuilding the linear operator (see below).
 *
 *  - Stage 2: building the `xform` operator which computes the nonuniform
 *    Fourier transform of the image.  This depends on parameters such as the
 *    image dimensions, the pixel size and the specific model for the transform
 *    to use.
 */
func mira_new(.., target=, wavemin=, wavemax=, pixelsize=, dims=, xform=,
              smearingfunction=, smearingfactor=,
              quiet=, noise_method=, noise_level=, baseline_precision=)
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

   SEE ALSO mira_config.
*/
{
  /* Get precision (in meters) for rounding baselines. */
  if (! scalar_double(baseline_precision, 1e-3) || baseline_precision <= 0) {
    error, "absolute baseline precision must be a strictly positive length";
  }

  /* Create an empty MiRA instance. */
  if (! scalar_string(target, "*")) {
    error, "target must be a string";
  }
  master = h_new(oidata = oifits_new(), target = target,
                 stage = 0, xform="nfft",
                 smearingfactor = 1.0, smearingfunction = sinc,
                 dims=[2,128,128], pixelsize = 5*MIRA_MILLIARCSECOND,
                 wavemin = 0.0, wavemax = MIRA_INF,
                 first = [],
                 model = h_new(),
                 baseline_precision = baseline_precision,
                 flags = (MIRA_FIT_VISAMP | MIRA_FIT_VISPHI | MIRA_FIT_VIS2 |
                          MIRA_FIT_T3AMP | MIRA_FIT_T3PHI |
                          MIRA_CONVEX_APPROX | MIRA_VON_MISES_APPROX));

  /* Apply configuration options. */
  mira_config, master, wavemin=wavemin, wavemax=wavemax, flags=flags,
    pixelsize=pixelsize, dims=dims, xform=xform,
    smearingfunction=smearingfunction, smearingfactor=smearingfactor;

  /* Load OI-FITS data file(s). */
  while (more_args()) {
    mira_add_oidata, master, next_arg(), quiet=quiet,
      noise_method=noise_method, noise_level=noise_level;
  }

  return master;
}

/*---------------------------------------------------------------------------*/
/* CONFIGURING AN INSTANCE */

func mira_config(master, wavemin=, wavemax=, dims=, pixelsize=, xform=,
                 flags=, smearingfunction=, smearingfactor=)
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

     Keyword `dims` specifies the dimension(s) of the restored image.  If it is
     a scalar, it specifies the width and height of the image; otherwise,
     `dims` can be `[width,height]` or `[2,width,height]`.


   SEE ALSO mira_new.
*/
{
  /* Check all options. */
  if (! mira_length(wavemin, master.wavemin) || wavemin < 0.0) {
    error, "minimum wavelength must be a nonnegative length";
  }
  if (! mira_length(wavemax, master.wavemax) || wavemax < 0.0) {
    error, "maximum wavelength must be a nonnegative length";
  }
  if (wavemin > wavemax) {
    error, "minimum wavelength > maximum wavelength";
  }
  if (! mira_angle(pixelsize, master.pixelsize) || pixelsize <= 0.0) {
    error, "pixel size must be a strictly positive angle";
  }
  master_xform = mira_xform_name(master);
  if (! scalar_string(xform, master_xform)) {
    error, "`xform` must be a string";
  }
  if (xform != "simple" && xform != "nfft") {
    error, "unknown `xform` name";
  }
  if (xform != "simple" && (! is_void(smearingfactor) &&
                            smearingfactor != 0)) {
    warn, "use xform=\"simple\" to account for bandwidth smearing";
  }
  if (! scalar_double(smearingfactor, master.smearingfactor) ||
      smearingfactor < 0) {
    error, "keyword SMEARINGFACTOR must have a nonnegative value";
  }
  if (is_void(smearingfunction)) {
    smearingfunction = master.smearingfunction;
  }
  if (! is_func(smearingfunction) && ! is_hash(smearingfunction)) {
    error, "keyword SMEARINGFUNCTION must be set with a function";
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
    error, "invalid image dimension(s)";
  }
  master_flags = master.flags;
  if (! is_scalar(master_flags) || structof(master_flags) != int) {
    if (! is_scalar(master_flags) || int(master_flags) != master_flags) {
      error, "invalid current flags";
    }
    master_flags = int(master.flags);
    h_set, master, flags = master_flags;
  }
  if (is_void(flags)) {
    flags = master_flags;
  }
  if (! is_scalar(flags) || ! is_integer(flags) || int(flags) != flags) {
    error("invalid flags");
  }
  flags = int(flags);

  /* Apply changes if any. */
  if (master.wavemin != wavemin || master.wavemax != wavemax) {
    h_set, master, model = h_new(), stage = max(master.stage, 0),
      wavemin = wavemin, wavemax = wavemax;
  }
  if (master_flags != flags) {
    h_set, master, model = h_new(), stage = max(master.stage, 0),
      flags = flags;
  }
  if (master.pixelsize != pixelsize) {
    h_set, master, model = h_new(), stage = max(master.stage, 1),
      pixelsize = pixelsize;
  }
  if (master_xform != xform || smearingfunction != master.smearingfunction ||
      smearingfactor != master.smearingfactor) {
    h_set, master, model = h_new(), stage = max(master.stage, 1),
      xform = xform, smearingfunction = smearingfunction,
      smearingfactor = smearingfactor;
  }
  if (! mira_same_dimensions(master.dims, dims)) {
    h_set, master, model = h_new(), stage = max(master.stage, 1),
      dims = dims;
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

func mira_pixel_size(master)
/* DOCUMENT mira_pixel_size(master);

     yields the angular size (in radians) of the pixels for the image restored
     with MiRA instance `master`.  The image pixel size can be changed with
     `mira_config`.

   SEE ALSO: mira_config.
 */
{
  return master.pixelsize;
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
func mira_flags(master)
/* DOCUMENT mira_flags(master);

     yields the current flags used by MiRA instance `master`.  The flags may be
     changed with `mira_config` and the following bits can be combined into the
     flags:

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
MIRA_HANIFF_APPROX    = (1n << 6);
MIRA_CONVEX_LIMIT     = (2n << 6);
MIRA_VON_MISES_APPROX = (3n << 6);

MIRA_FIT_VIS = (MIRA_FIT_VISAMP | MIRA_FIT_VISPHI);
MIRA_FIT_T3  = (MIRA_FIT_T3AMP | MIRA_FIT_T3PHI);

/*---------------------------------------------------------------------------*/
/* UPDATING RESSOURCES AND MODEL */

func mira_update(master, img)
/* DOCUMENT mira_update, master;
         or mira_update, master, img;

     Update the internal ressources of MiRA instance `master` and, if image
     `img` is provided, update the current model of the complex visibilities.

     When called as a function, `master` is returned.

   SEE ALSO: mira_cost, mira_cost_and_gradient, mira_model_vis, mira_model_vis2,
             mira_model_re, mira_model_im, mira_model_amp, mira_model_phi.
 */
{
  if (master.stage < 1) {
    _mira_select_data, master;
  }
  if (master.stage < 2) {
    _mira_define_xform, master;
  }
  if (! is_void(img)) {
    /* Compute complex visibilities (derived quantities will be computed on
       demand). */
    vis = mira.xform(img);
    h_set, master.model, img=img, vis=vis, re=[], im=[],
      vis2=[], amp=[], phi=[];
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

/*---------------------------------------------------------------------------*/
