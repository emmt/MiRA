/*
 * mira2_dirty.i -
 *
 * Dirty beam, dirty map and residual map computations in MiRA.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of MiRA, a "Multi-aperture Image Reconstruction
 * Algorithm", <https://github.com/emmt/MiRA>.
 *
 * Copyright (C) 2001-2020, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

func mira_compute_dirty_beam(master, method)
/* DOCUMENT psf = mira_compute_dirty_beam(master);
         or psf = mira_compute_dirty_beam(master, method);

     yields the "dirty beam", that the equivalent point spread function (PSF),
     given by the (u,v) coverage in the interferometric data stored by MASTER.

     If optional argument METHOD is "xform", the actual pixel to complex
     visibilities transform is used to compute the dirty beam.  If optional
     argument METHOD is not specified or is "exact", the equations are applied
     with no other approximations.

*/
{
  /*
   * The complex visibilities are given by:
   *
   *     vis(u,v) = sum_{x,y} H(x,y,u,v,λ,Δλ) × img(x,y)                    (1)
   *
   * where `x` is the relative right ascension (RA), `y` is the relative
   * declination (DEC), `(u,v)` are the baseline coordinates, `λ` is the
   * wavelength, `Δλ` is the effective spectral bandwidth and `H(x,y,u,v,λ,Δλ)`
   * is given by:
   *
   *     H(x,y,u,v,λ,Δλ) = exp(-i⋅2⋅π⋅(x⋅u + y⋅v)/λ) × s(γ⋅Δλ⋅(x⋅u + y⋅v)/λ²)
   *
   * with `γ≈1` a correction factor and `s` the smearing function (specified by
   * the `smearingfactor` and the `smearingfunction` options).
   *
   * Neglecting bandwidth smearing:
   *
   *     H(x,y,u,v,λ,Δλ) = exp(-i⋅2⋅π⋅(x⋅u + y⋅v)/λ)                        (2)
   *
   * is the complex exponential is given by Eq. (1) in Pauls, Young, Cotton &
   * Monnier, "A Data Exchange Standard for Optical (Visible/IR)
   * Interferometry", Publications of the Astronomical Society of the Pacific
   * (PASP), vol. 117, pp. 1255-1262 (2005).  Then Eq. (1) is similar to the
   * Fourier transform of the convolution of the image by a point spread
   * function (PSF) called the "dirty beam".
   *
   * The pseudo-inverse is given by:
   *
   *     inv(x,y) = sum_{u,v} Re(conj(H(x,y,u,v,λ,Δλ)) × vis(u,v))          (3)
   *
   * Assuming `vis(u,v) = 1` for all `(u,v)` yields the "dirty beam":
   *
   *     psf(x,y) = sum_{u,v} Re(conj(H(x,y,u,v,λ,Δλ)))
   *              = sum_{u,v} cos(2⋅π⋅(x⋅u + y⋅v)/λ)                        (4)
   *
   * Exploit the fact that:
   *
   *     cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
   *
   * to compute the PSF with separable expressions.  This saves computations
   * and memory.  Then, using Yorick's notation, computing the PSF is pretty
   * straightforward:
   */
  if (is_void(method) || method == "exact") {
    /* Extract necessary parameters. */
    local u, v, w, x, y, wave, band;
    mira_update, master;
    eq_nocopy, x,    mira_image_x(master);
    eq_nocopy, y,    mira_image_y(master);
    eq_nocopy, u,    master.coords.unique.u;
    eq_nocopy, v,    master.coords.unique.v;
    eq_nocopy, wave, master.coords.unique.wave;
    eq_nocopy, band, master.coords.unique.band;
    nx = numberof(x);
    ny = numberof(y);
    if (! vector_double(x)) error, "X must be a vector of reals";
    if (! vector_double(y)) error, "Y must be a vector of reals";
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

    /* Apply equations. FIXME: Add zero-th frequency. */
    q = MIRA_TWO_PI/wave;
    a = q*u*x(-,);
    b = q*v*y(-,);
    psf = (cos(a))(+,)*(cos(b))(+,) - (sin(a))(+,)*(sin(b))(+,);

  } else if (method == "xform") {
    /* Apply the adjoint of the pixel to visibilities transform. FIXME: Add
       zero-th frequency. */
    mira_update, master;
    nfreqs = numberof(master.coords.unique.u);
    (vis = array(double, 2, nfreqs))(1,) = 1.0;
    psf = master.xform(vis, 1n);
  } else {
    error, "invalid method \""+method+"\" for the dirty beam";
  }

  return psf;
}
