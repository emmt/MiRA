/*
 * mira2_xform.i -
 *
 * Models of the image to complex visibilities transform in MiRA (nonuniform
 * Fourier transform and, optionally, spectral bandwidth smearing).
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


/*---------------------------------------------------------------------------*/
/* MODELS OF THE IMAGE TO COMPLEX VISIBILITIES TRANSFORM */

local _mira_apply_simple_xform, _mira_apply_nfft_xform;
func _mira_define_xform(master)
/* DOCUMENT _mira_define_xform, master;

   Private routine to define the operator which implements the linear
   transform of the image into the complex visibilities in MiRA instance
   `master`.  This routine must be called after data selection (by
   `_mira_select_data`) and after changing any options related to the model.
   Attribute `master.stage` is 1 on entry and is set to 2 on exit.  When
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
  if (master.stage != 1) {
    error, "private routine called at wrong stage";
  }
  coords = master.coords;
  prec = master.baseline_precision;
  xform = mira_xform_name(master);
  pixelsize = mira_pixel_size(master);
  smearingfactor = master.smearingfactor;
  smearingfunction = master.smearingfunction;

  /* Figure out whether or not account for spectral bandwidth smearing. */
  smearing = (smearingfunction != "none" && smearingfactor > 0);
  if (smearing && xform == "nfft") {
    warn, ("Spectral bandwidth smearing is ignored with xform=\"nfft\" "+
           "(set smearingfunction=\"none\" or smearingfactor=0 to avoid "+
           "this message).");
    smearing = 0n;
  }
  if (! smearing && xform != "nfft") {
    warn, ("when spectral bandwidth smearing is ignored, xform=\"nfft\" "+
           "is faster that xform=\""+xform+"\".");
  }
  if (smearing && xform == "separable" && smearingfunction != "gauss") {
    warn, ("Non-Gaussian smearing function with separable model is not "+
           "recommended (set smearingfunction=\"gauss\" or "+
           "xform=\"nonseparable\" to avoid this message).");
  }

  /* Select reduced list of coordinates. */
  if (smearing) {
    mode = 3;
  } else {
    mode = 1;
  }
  _mira_reduce_coordinates, master, mode;

  /* Check pixel size. */
  maxpixelsize = mira_maximum_pixel_size(master);
  if (pixelsize > maxpixelsize) {
    warn, ("Current pixelsize (%g mas) is too large which will yield " +
           "aliasing.  Maximum pixelsize to avoid aliasing is %g mas."),
      pixelsize/MIRA_MILLIARCSECOND, maxpixelsize/MIRA_MILLIARCSECOND;
  }

  /* Get the coordinates (and check them). */
  local u, v, wave, band, ufreq, vfreq;
  vector_double = mira_vector_double;
  nx = mira_image_size(master, 1);
  ny = mira_image_size(master, 2);
  x = mira_sky_coordinates(nx, pixelsize);
  y = mira_sky_coordinates(ny, pixelsize);
  if (! vector_double(x)) error, "X must be a vector of reals";
  if (! vector_double(y)) error, "Y must be a vector of reals";
  if (mode == 1) {
    eq_nocopy, ufreq, master.coords.ufreq;
    eq_nocopy, vfreq, master.coords.vfreq;
    if (! vector_double(ufreq)) error, "UFREQ must be a vector of reals";
    if (! vector_double(vfreq)) error, "VFREQ must be a vector of reals";
    if (numberof(vfreq) != numberof(ufreq)) {
      error, "UFREQ and VFREQ must have the same number of elements";
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
      error, "U, V, WAVE and BAND must have the same number of elements";
    }
  }

  /* Convert smearing function name into a function. */
  if (smearingfunction == "none") {
    smearingfunction = _mira_one();
  } else if (smearingfunction == "gauss") {
    smearingfunction = mira_gaussian_smearing;
  } else if (smearingfunction == "sinc") {
    smearingfunction = sinc;
  } else {
    throw, "unknown `smearingfunction` name";
  }

  /* This model of the nonuniform Fourier transform is implemented by a
   * simple matrix multiplication whose coefficients are given by:
   *
   *     exp(-i*2*pi*(x*u + y*v))
   *
   * where `x` is the relative right ascension (RA), `y` is the relative
   * declination (DEC) and `(u,v)` are the spatial frequencies (baselines
   * divided by wavelength).  See Eq. (1) in Pauls, Young, Cotton & Monnier,
   * "A Data Exchange Standard for Optical (Visible/IR) Interferometry",
   * Publications of the Astronomical Society of the Pacific (PASP),
   * vol. 117, pp. 1255-1262 (2005).
   *
   * FIXME: The nonuniform Fourier transform is separable in RA/DEC, the
   *        smearing is also separable provided a Gaussian is used (instead
   *        of a sinc).  Exploiting this can lead to a much smaller number of
   *        coefficients.
   */

  if (xform == "separable") {

    if (smearing) {
      q = -2*MIRA_PI/wave;
      if (smearingfactor > 0 && max(band) > 0) {
        r = smearingfactor*band/(wave*wave);
        A1 = _mira_xform_coefs(smearingfunction(r*u*x(-,)), q*u*x(-,));
        A2 = _mira_xform_coefs(smearingfunction(r*v*y(-,)), q*v*y(-,));
      } else {
        smearing = 0n;
        A1 = _mira_xform_coefs([], q*u*x(-,));
        A2 = _mira_xform_coefs([], q*v*y(-,));
      }
    } else {
      q = -2*MIRA_PI;
      A1 = _mira_xform_coefs([], q*ufreq*x(-,));
      A2 = _mira_xform_coefs([], q*vfreq*y(-,));
    }
    xform = h_new(name=xform, A1=A1, A2=A2);
    h_evaluator, xform, "_mira_apply_separable_xform";

  } else if (xform == "nonseparable") {

    if (smearing) {
      q = -2*MIRA_PI/wave;
      if (smearingfactor > 0 && max(band) > 0) {
        r = smearingfactor*band/(wave*wave);
        fct = smearingfunction(r*u*x(-,) + r*v*y(-,-,));
      } else {
        smearing = 0n;
        fct = [];
      }
      A = _mira_xform_coefs(unref(fct), q*u*x(-,) + q*v*y(-,-,));
    } else {
      q = -2*MIRA_PI;
      A = _mira_xform_coefs([], q*ufreq*x(-,) + q*vfreq*y(-,-,));
    }
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
                  nfft = nfft_new(dims, nodes, flags=flags),
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
  return h_set(master, xform = xform, stage = 2);
}

/*---------------------------------------------------------------------------*/
/* APPLY THE MODELS */

func _mira_apply_separable_xform(this, x, job)
{
  if (! job) {
    /* Direct transform. */
    return ((this.A1(,,+)*x(+,))*this.A2)(,,sum);
  } else if (job == 1) {
    /* Apply adjoint operator */
    return ((this.A1*x)(*,))(+,)*((this.A2)(*,))(+,);
  } else {
    error, "unsupported value for JOB";
  }
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

func _mira_one(t)
{
  return array(1.0, dimsof(t));
}

/*---------------------------------------------------------------------------*/
/* UTILITIES FOR COMPUTING XFORM MODELS */

func _mira_xform_coefs(rho, phi)
/* DOCUMENT z = _mira_xform_coefs(rho, phi);
     yields a 2-by-dimsof(PHI) array Z such that:

         Z(1,..) = RHO*cos(PHI)
         Z(2,..) = RHO*sin(PHI)

     if RHO is not void and:

         Z(1,..) = cos(PHI)
         Z(2,..) = sin(PHI)

     otherwise. */
{
  /* Note: use `unref` to reduce the memory requirements (only works if
     arguments are temporary expressions). */
  if (is_void(rho)) {
    z = array(double, 2, dimsof(phi));
    z(1,..) = cos(phi);
    z(2,..) = sin(unref(phi));
  } else {
    z = array(double, 2, dimsof(rho, phi));
    z(1,..) = rho*cos(phi);
    z(2,..) = unref(rho)*sin(unref(phi));
  }
  return z;
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
  /* Check stage. */
  if (master.stage != 1) {
    throw, "private routine called at wrong stage";
  }

  /* Get coordinates (should be flat vectors of doubles or void). */
  coords = master.coords;
  local u, v, wave, band;
  eq_nocopy, u, coords.u;
  eq_nocopy, v, coords.v;
  eq_nocopy, wave, coords.wave;
  eq_nocopy, band, coords.band;

  /* Check mode. */
  if (! scalar_long(mode, 1) || mode < 1 || mode > 3) {
    error, "invalid MODE";
  }

  /* Sort coordinates. */
  if (mode == 1) {
    ufreq = u/wave;
    vfreq = v/wave;
    i = msort(ufreq, vfreq);
  } else if (mode == 2) {
    i = msort(wave, u, v);
  } else {
    i = msort(wave, band, u, v);
  }

  /* Select unique coordinates. */
  i1 = i(1:-1);
  i2 = i(2:0);
  if (mode == 1) {
    uniq = (ufreq(i2) != ufreq(i1)) | (vfreq(i2) != vfreq(i1));
  } else {
    uniq = (u(i2) != u(i1)) | (v(i2) != v(i1));
  }
  if (mode > 1) uniq |= (wave(i2) != wave(i1));
  if (mode > 2) uniq |= (band(i2) != band(i1));
  uniq = grow(1n, uniq);

  /* Build indirection indices and extract unique coordinates. */
  (idx = array(long, numberof(i)))(i) = long(uniq)(psum);
  i = i(where(uniq));
  if (mode == 1) {
    coords = h_new(ufreq = ufreq(i), vfreq = vfreq(i));
    if (debug && (max(abs(ufreq - coords.ufreq(idx))) != 0 ||
                  max(abs(vfreq - coords.vfreq(idx))) != 0)) {
      throw, "bug: something wrong with `(u/λ,v/λ)`!";
    }
  } else {
    coords = h_new(u = u(i), v = v(i));
    if (debug && (max(abs(u - coords.u(idx))) != 0 ||
                  max(abs(v - coords.v(idx))) != 0)) {
      throw, "bug: something wrong with `(u,v)`!";
    }
    if (mode > 1) {
      h_set, coords, wave = wave(i);
      if (debug && max(abs(wave - coords.wave(idx))) != 0) {
        throw, "bug: something wrong with `wave`!";
      }
    }
    if (mode > 2) {
      h_set, coords, band = band(i);
      if (debug && max(abs(band - coords.band(idx))) != 0) {
        throw, "bug: something wrong with `band`!";
      }
    }
  }

  /* Overwrite coordinates and indirections for every data block. */
  h_set, master, stage=-1; // in case of interrupts
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    if (h_has(db, "idx")) {
      h_set, db, idx = idx(db.idx);
    }
  }

  /* Overwrite coordinates in MASTER and update stage. */
  return h_set(master, coords=h_set(coords, mode=mode), stage=2);
}

/*---------------------------------------------------------------------------*/
