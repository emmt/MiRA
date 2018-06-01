/*
 * mira2_xform.i -
 *
 * Models of the image to complex visibilities transform in MiRA (nonuniform
 * Fourier transform and, optionally, spectral bandwidth smearing).
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


/*---------------------------------------------------------------------------*/
/* MODELS OF THE IMAGE TO COMPLEX VISIBILITIES TRANSFORM */

local _mira_apply_separable_xform;
local _mira_apply_nonseparable_xform
local _mira_apply_nfft_xform;
func _mira_define_xform(master)
/* DOCUMENT _mira_define_xform, master;

   Private routine to define the operator which implements the linear
   transform of the image into the complex visibilities in MiRA instance
   `master`.  This routine must be called after data selection (by
   `_mira_select_data`) and after changing any options related to the model.
   Attribute `master.stage` is 1 on entry and is set to 3 on exit.  When
   called as a fucntion, `master` is returned.

   The operator, say `xform`, is an object (actually a hash table with its own
   evaluator) which can be used as follows:

       vis = xform(img)

   to compute the model visibilities `vis`, such that `vis(1,..)` and
   `vis(2,..)` are respectively the real and imaginary parts of the complex
   visibilities.  The transpose of the operator can also be applied (for
   instance to compute the gradient of the likelihood):

       xform(inp, 1)

   where input array `inp` is a 2-by-`nfreqs` array, with `nfreqs` the number
   of complex visibilities.

   SEE ALSO: mira_update, _mira_select_data.
*/
{
  /* Check assumptions and get the parameters needed to define the
     transform. */
  if (master.stage < 1 || master.stage > 2) {
    error, "private routine called at wrong stage";
  }
  coords = master.coords;
  prec = master.baseline_precision;
  xform = mira_xform_name(master);
  pixelsize = mira_pixel_size(master);
  smearingfactor = master.smearingfactor;
  smearingfunction = master.smearingfunction;


  /* Figure out whether or not to account for spectral bandwidth smearing. */
  smearing = (smearingfunction != "none" && smearingfactor > 0);
  if (smearing && (xform == "nfft" || xform == "separable")) {
    if (! MIRA_QUIET) {
      warn, ("Spectral bandwidth smearing is ignored with xform=\"%s\" "+
           "(set smearingfunction=\"none\" or smearingfactor=0 to avoid "+
           "this message)."), xform;
    }
    smearing = 0n;
  }
  if (! MIRA_QUIET && ! smearing && xform != "nfft") {
    warn, ("When spectral bandwidth smearing is ignored, xform=\"nfft\" "+
           "is faster that xform=\""+xform+"\".");
  }

  /* Select reduced list of coordinates. */
  flags = master.flags & (MIRA_KEEP_WAVELENGTH|MIRA_KEEP_BANDWIDTH);
  if (smearing) {
    flags |= (MIRA_KEEP_WAVELENGTH|MIRA_KEEP_BANDWIDTH);
  }
  if (flags == (MIRA_KEEP_WAVELENGTH|MIRA_KEEP_BANDWIDTH)) {
    mode = 3;
  } else if (flags == MIRA_KEEP_WAVELENGTH) {
    mode = 2;
  } else {
    mode = 1;
  }
  if (master.stage != 1 && master.coords.mode != mode) {
    /* Must redo data selection. */
    h_set, master, stage=0;
    _mira_select_data, master;
  }
  if (master.stage < 2) {
    _mira_reduce_coordinates, master, mode;
  }

  /* Check pixel size. */
  maxpixelsize = mira_maximum_pixel_size(master);
  if (pixelsize > maxpixelsize) {
    warn, ("Current pixelsize (%g mas) is too large which will yield " +
           "aliasing.  Maximum pixelsize to avoid aliasing is %g mas."),
      pixelsize/MIRA_MILLIARCSECOND, maxpixelsize/MIRA_MILLIARCSECOND;
  }

  /* Get the coordinates (and check them). */
  local u, v, wave, band, ufreq, vfreq, x, y;
  eq_nocopy, x, mira_image_x(master);
  eq_nocopy, y, mira_image_y(master);
  nx = numberof(x);
  ny = numberof(y);
  if (! vector_double(x)) error, "X must be a vector of reals";
  if (! vector_double(y)) error, "Y must be a vector of reals";
  if (mode == 1) {
    eq_nocopy, ufreq, master.coords.ufreq;
    eq_nocopy, vfreq, master.coords.vfreq;
    if (! vector_double(ufreq)) error, "UFREQ must be a vector of reals";
    if (! vector_double(vfreq)) error, "VFREQ must be a vector of reals";
    if (numberof(vfreq) != numberof(ufreq)) {
      error, "UFREQ and VFREQ must have the same size";
    }
  } else if (mode == 2) {
    eq_nocopy, u,    master.coords.u;
    eq_nocopy, v,    master.coords.v;
    eq_nocopy, wave, master.coords.wave;
    if (! vector_double(u)) error, "U must be a vector of reals";
    if (! vector_double(v)) error, "V must be a vector of reals";
    if (! vector_double(wave) || min(wave) <= 0) {
      error, "WAVE must be a vector of strictly positive reals";
    }
    m = numberof(u);
    if (numberof(v) != m || numberof(wave) != m) {
      error, "U, V and WAVE must have the same size";
    }
  } else {
    eq_nocopy, u,    master.coords.u;
    eq_nocopy, v,    master.coords.v;
    eq_nocopy, wave, master.coords.wave;
    eq_nocopy, band, master.coords.band;
    if (! vector_double(u)) error, "U must be a vector of reals";
    if (! vector_double(v)) error, "V must be a vector of reals";
    if (! vector_double(wave) || min(wave) <= 0) {
      error, "WAVE must be a vector of strictly positive reals";
    }
    if (! vector_double(band) || min(band) < 0) {
      error, "BAND must be a vector of nonnegative reals";
    }
    m = numberof(u);
    if (numberof(v) != m || numberof(wave) != m || numberof(band) != m) {
      error, "U, V, WAVE and BAND must have the same size";
    }
  }
  if (mode > 1) {
    ufreq = u/wave;
    vfreq = v/wave;
    h_set, master.coords, ufreq=ufreq, vfreq=vfreq;
  }

  /* Convert smearing function name into a function. */
  if (smearingfunction == "none") {
    smearingfunction = mira_no_smearing();
  } else if (smearingfunction == "gauss") {
    smearingfunction = mira_gaussian_smearing;
  } else if (smearingfunction == "sinc") {
    smearingfunction = sinc;
  } else {
    throw, "unknown `smearingfunction` name";
  }

  /*
   * The transform model is linear and can be implemented by a simple matrix
   * multiplication whose coefficients are given by:
   *
   *     exp(-i⋅2⋅π⋅(x⋅u + y⋅v)/λ)⋅s(γ⋅Δλ⋅(x⋅u + y⋅v)/λ²)
   *
   * where `x` is the relative right ascension (RA), `y` is the relative
   * declination (DEC), `(u,v)` are the baselines divided, `λ` is the
   * wavelength, `Δλ` is the effective spectral bandwidth, `γ≈1` is a
   * correction factor and `s` is the smearing function (specified by the
   * `smearingfactor` and the `smearingfunction` options).  The complex
   * exponential is given by Eq. (1) in Pauls, Young, Cotton & Monnier, "A Data
   * Exchange Standard for Optical (Visible/IR) Interferometry", Publications
   * of the Astronomical Society of the Pacific (PASP), vol. 117, pp. 1255-1262
   * (2005).
   *
   * The nonuniform Fourier transform is separable in RA/DEC.  Exploiting this
   * can lead to a much smaller number of coefficients if spectral bandwidth
   * smearing is ignored.
   */

  if (xform == "separable") {

    local A1re, A1im, A2re, A2im;
    q = -2*MIRA_PI;
    _mira_cos_sin, A1re, A1im, q*ufreq*x(-,);
    _mira_cos_sin, A2re, A2im, q*vfreq*y(-,);
    xform = h_new(name=xform, A1re=A1re, A1im=A1im, A2re=A2re, A2im=A2im);
    h_evaluator, xform, "_mira_apply_separable_xform";

  } else if (xform == "nonseparable") {

    local A1re, A1im, A2re, A2im;
    if (smearing) {
      q = -2*MIRA_PI/wave;
      _mira_cos_sin, A1re, A1im, q*u*x(-,);
      _mira_cos_sin, A2re, A2im, q*v*y(-,);
      if (smearingfactor > 0 && max(band) > 0) {
        r = smearingfactor*band/(wave*wave);
        fct = smearingfunction(r*u*x(-,) + r*v*y(-,-,));
      } else {
        smearing = 0n;
        fct = [];
      }
    } else {
      q = -2*MIRA_PI;
      _mira_cos_sin, A1re, A1im, q*ufreq*x(-,);
      _mira_cos_sin, A2re, A2im, q*vfreq*y(-,);
    }
    Are = A1re*A2re(,-,) - A1im*A2im(,-,);
    Aim = A1re*A2im(,-,) + A1im*A2re(,-,);
    if (smearing) {
      Are = unref(Are)*fct;
      Aim = unref(Aim)*unref(fct);
    }
    A = _mira_fake_complex(unref(Are), unref(Aim));
    xform = h_new(name=xform, A=A);
    h_evaluator, xform, "_mira_apply_nonseparable_xform";

  } else if (xform == "nfft") {

    /* This model of the nonuniform Fourier transform is similar to the
       "simple" model above but is much faster as it uses `nfft` to compute the
       transform. */

    /* Make sure NFFT plugin is installed and loaded. */
    if (! is_func(nfft_new)) {
      include, "nfft.i", 3;
      if (! is_func(nfft_new)) {
        inform, ("YNFFT plugin is not installed, "+
                 "see \"https://github.com/emmt/ynfft\".");
        inform, ("Depending on your system, you may just have to "+
                 "install package \"yorick-ynfft\".");
        throw, "plugin YNFFT is not installed";
      }
    }

    /*
     * In NFFT:
     *   x = [-nx/2, 1-nx/2, ..., nx/2-1]*step
     *   y = [-ny/2, 1-ny/2, ..., ny/2-1]*step
     * plus NX and NY must be even.
     */
    local r1, r2;
    if ((nx & 1) == 1) {
      n1 = nx + 1;
      r1 = 2 : n1;
    } else {
      n1 = nx;
    }
    if ((ny & 1) == 1) {
      n2 = ny + 1;
      r2 = 2 : n2;
    } else {
      n2 = ny;
    }
    flags = NFFT_SORT_NODES;
    nodes = [pixelsize*ufreq, pixelsize*vfreq];
    dims = [n1, n2];
    xform = h_new(name = xform, ufreq = ufreq, vfreq = vfreq,
                  nfft = nfft_new(dims, nodes, flags=flags,
                                  nthreads=master.nthreads),
                  n1 = n1, r1 = r1,
                  n2 = n2, r2 = r2,
                  sub = (n1 != nx || n2 != ny));
    h_evaluator, xform, "_mira_apply_nfft_xform";

  } else {
    throw, "invalid `xform` name";
  }


  /* Install the evaluator. */
  if (smearing) {
    h_set, xform, smearingfactor=smearingfactor,
      smearingfunction=master.smearingfunction;
  } else {
    h_set, xform, smearingfactor=1.0, smearingfunction="none";
  }
  return h_set(master, xform = xform, stage = 3);
}

/*---------------------------------------------------------------------------*/
/* APPLY THE MODELS */

func _mira_apply_separable_xform(this, x, job)
{
  if (! job) {
    /* Direct transform. */
    re = this.A1re(,+)*x(+,);
    im = this.A1im(,+)*x(+,);
    return _mira_fake_complex((this.A2re*re)(,sum) - (this.A2im*im)(,sum),
                              (this.A2re*im)(,sum) + (this.A2im*re)(,sum));
  } else if (job == 1) {
    /* Apply adjoint operator */
    if (numberof((dims = dimsof(x))) != 3 || dims(2) != 2) {
      throw, "expecting a 2-by-m array";
    }
    re = x(1,);
    im = x(2,);
    return (this.A1re(+,)*(this.A2re*re + this.A2im*im)(+,) +
            this.A1im(+,)*(this.A2re*im - this.A2im*re)(+,));
  } else {
    error, "unsupported value for JOB";
  }
}

func _mira_fake_complex(re, im)
{
  z = array(double, 2, dimsof(re, im));
  z(1,..) = re;
  z(2,..) = im;
  return z;
}

func _mira_apply_nonseparable_xform(this, x, job)
{
  return mvmult(this.A, x, job);
}

func _mira_apply_nfft_xform(this, x, job)
{
  local z;
  if (! job) {
    /* direct operator */
    if (this.sub) {
      tmp = array(complex, this.n1, this.n2);
      tmp(this.r1, this.r2) = x;
      eq_nocopy, x, tmp;
    }
    reshape, z, &this.nfft(x), double, 2, this.nfft.num_nodes;
    return z;
  } else if (job == 1) {
    /* adjoint operator */
    z = this.nfft(mira_cast_real_as_complex(x), 1n);
    if (this.sub) {
      return double(z)(this.r1, this.r2);
    } else {
      return double(z);
    }
  } else {
    error, "unsupported value for JOB";
  }
}

/*---------------------------------------------------------------------------*/
/* SMEARING FUNCTIONS */

local MIRA_GAUSS_FWHM, MIRA_SINC_FWHM;
func mira_gaussian_smearing(t)
/* DOCUMENT mira_gaussian_smearing(t);

     This function approximates the cardinal sine function by a Gaussian which
     has the same full width at half-maximum (FHWM) as the sinc function.  The
     approximation error in smaller than 2% (relative to the peak) on the
     interval t ∈ [-0.682,+0.682].

     Constants MIRA_GAUSS_FWHM and MIRA_SINC_FWHM are set with the respective
     FHWM of a normal distribution and of the sinc fucntion.

   SEE ALSO: sinc.
 */
{
  /* The following constant is to have the same FWHM as the sinc function. */
  K = MIRA_GAUSS_FWHM/MIRA_SINC_FWHM;
  return exp((-K*K/2.0)*t*t);
}

MIRA_SINC_FWHM = 1.2067091288032283;
MIRA_GAUSS_FWHM = sqrt(log(256));

func mira_no_smearing(t)
{
  return array(1.0, dimsof(t));
}

/*---------------------------------------------------------------------------*/
/* UTILITIES FOR COMPUTING XFORM MODELS */

func _mira_cos_sin(&re, &im, phi)
/* DOCUMENT _mira_cos_sin, re, im, phi;
     overwrites caller's variables RE and IM with:

         re = cos(PHI)
         im = sin(PHI)

     otherwise.
*/
{
  /* Note: use `unref` to reduce the memory requirements (only works if
     arguments are temporary expressions). */
  re = cos(phi);
  im = sin(unref(phi));
}

func mira_save_visibilities(master, fh)
{
  local idx;
  coords = master.coords;
  if (h_has(coords, "unique")) {
    coords = coords.unique;
  }

  fits_new_bintable, fh;
  fits_set, fh, "EXTNAME", "MODEL-VISIBILITIES";
  ptr = array(pointer, 6);

  fits_set, fh, "TTYPE1", "model_visamp", "Model of the visibility amplitudes";
  fits_set, fh, "TFORM1", "1D";
  ptr(1) = &mira_model_amp(master, idx);

  fits_set, fh, "TTYPE2", "model_visphi", "Model of the visibility phases";
  fits_set, fh, "TFORM2", "1D";
  fits_set, fh, "TUNIT2", "deg";
  ptr(2) = &(mira_model_phi(master, idx)*(180/pi));

  fits_set, fh, "TTYPE3", "ucoord", "U coordinates of baselines";
  fits_set, fh, "TFORM3", "1D";
  fits_set, fh, "TUNIT3", "m";
  ptr(3) = &(coords.u);

  fits_set, fh, "TTYPE4", "vcoord", "V coordinates of baselines";
  fits_set, fh, "TFORM4", "1D";
  fits_set, fh, "TUNIT4", "m";
  ptr(4) = &(coords.v);

  fits_set, fh, "TTYPE5", "eff_wave", "Effective wavelengths";
  fits_set, fh, "TFORM5", "1D";
  fits_set, fh, "TUNIT5", "nm";
  ptr(5) = &(coords.wave*1e9);

  fits_set, fh, "TTYPE6", "eff_band", "Effective spectral bandwidths";
  fits_set, fh, "TFORM6", "1D";
  fits_set, fh, "TUNIT6", "nm";
  ptr(6) = &(coords.band*1e9);

  fits_write_bintable, fh, ptr;
  fits_pad_hdu, fh;
}

/*---------------------------------------------------------------------------*/

func _mira_reduce_coordinates(master, mode, debug=)
/* DOCUMENT _mira_reduce_coordinates, master;
         or _mira_reduce_coordinates, master, mode;

     Private subroutine to reduce the list of coordinates stored in MiRA
     instance MASTER to a minimum.  MODE=1 (the default) to only consider
     unique values of (u/wave,v/wave), MODE=2 to consider unique values of
     (u,v,wave) and MODE=3 to consider unique values of (u,v,wave,band).  Here
     u, v, wave and band correspond to the contents of the fields of
     MASTER.coords.

*/
{
  /* Check stage and mode. */
  if (master.stage != 1) {
    throw, "private routine called at wrong stage";
  }
  if (! scalar_long(mode, 1) || mode < 1 || mode > 3) {
    error, "invalid MODE";
  }

  /* Get coordinates (should be flat vectors of doubles or void). */
  coords = master.coords;
  local u, v, wave, band;
  eq_nocopy, u, coords.u;
  eq_nocopy, v, coords.v;
  eq_nocopy, wave, coords.wave;
  eq_nocopy, band, coords.band;

  /* Uniquely sort all coordinates. */
  local sel, rev;
  _mira_unique, sel, rev, wave, band, u, v;
  wave = wave(sel);
  band = band(sel);
  u = u(sel);
  v = v(sel);
  coords = h_new(u=u, v=v, wave=wave, band=band);
  if (mode < 3) {
    /* Uniquely sort a restricted list of coordinates (only the ones which have
       an influence on the computation of the model complex visibilities). */
    local subrev;
    unique = coords;
    idx = indgen(numberof(sel));
    if (mode == 1) {
      ufreq = u/wave;
      vfreq = v/wave;
      _mira_unique, sel, subrev, ufreq, vfreq;
      coords = h_new(ufreq = ufreq(sel), vfreq = vfreq(sel));
    } else {
      _mira_unique, sel, subrev, wave, u, v;
      coords = h_new(wave = wave(sel), u = u(sel), v = v(sel));
    }
    idx = subrev(idx);
    rev = subrev(rev);
    h_set, coords, unique = h_set(unique, idx = idx);
  }

  /* Overwrite coordinates and indirections for every data block. */
  h_set, master, stage=-1; // in case of interrupts
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    if (h_has(db, "idx")) {
      h_set, db, idx = rev(db.idx);
    }
  }

  /* Overwrite coordinates in MASTER and update stage. */
  return h_set(master, coords=h_set(coords, mode=mode), stage=2);
}

func _mira_unique(&sel, &rev, x1, x2, x3, x4)
/* DOCUMENT local sel, rev; _mira_unique, sel, rev, x1, x2, ...;

     Yields arrays of indices sel and rev such that (x1(sel), x2(sel), ...)
     form a list of unique multiplets and such that x1(sel)(rev) == x1
     everywhere and similarly for x2 and the subsequent arguments.  The values
     x1(sel), x2(sel), ... are sorted with x1 the most significant key, x2 the
     next most significant key, and so on.  Arguments sel and rev must be
     simple variables, not expression nor values.

   SEE ALSO: msort.
*/
{
  /* Check inputs. */
  if (! is_void(x4)) {
    mode = 4;
  } else if (! is_void(x3)) {
    mode = 3;
  } else if (! is_void(x2)) {
    mode = 2;
  } else if (! is_void(x1)) {
    mode = 1;
  } else {
    error, "expecting at least one non-void argument";
  }
  if (! is_vector(x1)) {
    error, "X1 must be a vector";
  }
  n = numberof(x1);
  badlength = 0n;
  if (mode >= 2) {
    if (! is_vector(x2)) {
      error, "X2 must be a vector";
    }
    badlength = (badlength || (numberof(x2) != n));
  }
  if (mode >= 3) {
    if (! is_vector(x3)) {
      error, "X3 must be a vector";
    }
    badlength = (badlength || (numberof(x3) != n));
  }
  if (mode >= 4) {
    if (! is_vector(x4)) {
      error, "X4 must be a vector";
    }
    badlength = (badlength || (numberof(x4) != n));
  }
  if (badlength) {
    error, "all arguments must be vectors of same length";
  }

  /* Multi-key sort of coordinates. */
  if (mode == 1) {
    i = msort(x1);
  } else if (mode == 2) {
    i = msort(x1, x2);
  }  else if (mode == 3) {
    i = msort(x1, x2, x3);
  } else if (mode == 4) {
    i = msort(x1, x2, x3, x4);
  }

  /* Select unique coordinates. */
  i1 = i(1:-1);
  i2 = i(2:0);
  uniq = (x1(i2) != x1(i1));
  if (mode > 1) {
    uniq |= (x2(i2) != x2(i1));
    if (mode > 2) {
      uniq |= (x3(i2) != x3(i1));
      if (mode > 3) {
        uniq |= (x4(i2) != x4(i1));
      }
    }
  }
  uniq = grow(1n, uniq);

  /* Build indirection indices and extract unique coordinates. */
  (rev = array(long, numberof(i)))(i) = long(uniq)(psum);
  sel = i(where(uniq));
}

/*---------------------------------------------------------------------------*/
