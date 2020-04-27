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

/*
 * The complex visibilities are given by:
 *
 *     vis(u,v) = sum_{x,y} H(x,y,u,v,λ,Δλ) × img(x,y)                      (1)
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
 *     H(x,y,u,v,λ,Δλ) = exp(-i⋅2⋅π⋅(x⋅u + y⋅v)/λ)                          (2)
 *
 * is the complex exponential is given by Eq. (1) in Pauls, Young, Cotton &
 * Monnier, "A Data Exchange Standard for Optical (Visible/IR) Interferometry",
 * Publications of the Astronomical Society of the Pacific (PASP), vol. 117,
 * pp. 1255-1262 (2005).  Then Eq. (1) is similar to the Fourier transform of
 * the convolution of the image by a point spread function (PSF) called the
 * "dirty beam".
 *
 * The pseudo-inverse is given by:
 *
 *     inv(x,y) = sum_{u,v} Re(conj(H(x,y,u,v,λ,Δλ)) × vis(u,v))            (3)
 *
 * Assuming `vis(u,v) = 1` for all `(u,v)` yields the "dirty beam":
 *
 *     psf(x,y) = sum_{u,v} Re(conj(H(x,y,u,v,λ,Δλ)))
 *              = sum_{u,v} cos(2⋅π⋅(x⋅u + y⋅v)/λ)                          (4)
 *
 * Exploit the fact that:
 *
 *     cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
 *
 * to compute the PSF with separable expressions.  This saves computations and
 * memory.  Then, using Yorick's notation, computing the PSF is pretty
 * straightforward:
 */

func mira_apply_adjoint_of_exact_xform(master, vis)
/* DOCUMENT psf = mira_apply_adjoint_of_exact_xform(master);
         or img = mira_apply_adjoint_of_exact_xform(master, vis);

      yields the result of applying the adjoint of the exact pixels to complex
      visibilities transform in MASTER to the complex visibilities VIS (an
      array of dimensions 2-by-NFREQS with NFREQS the number of measured
      frequencies).  If VIS is not specified, it is assumed that VIS = 1 for
      every frequencies, that is the "dirty beam" is returned.

   SEE ALSO: mira_config, mira_compute_dirty_beam, mira_compute_dirty_map.
 */
{
  /* Make sure everything up to date. */
  mira_update, master;

  /* Extract necessary parameters. */
  local x, y, ufreq, vfreq;
  eq_nocopy, x,     mira_image_x(master);
  eq_nocopy, y,     mira_image_y(master);
  eq_nocopy, ufreq, mira_model_ufreq(master);
  eq_nocopy, vfreq, mira_model_vfreq(master);
  nx = numberof(x);
  ny = numberof(y);

  /* Compute trigonometric factors. */
  ux = MIRA_TWO_PI*ufreq*x(-,);
  cos_ux = cos(ux);
  sin_ux = sin(ux);
  ux = [];
  vy = MIRA_TWO_PI*vfreq*y(-,);
  cos_vy = cos(vy);
  sin_vy = sin(vy);
  vy = [];

  /*
   * Compute for all (x,y):
   *
   *     img(x,y) = sum_k Re(conj(H(k,x,y))*vis(k))
   *              = sum_k Re(exp(+2*i*pi*(u(k)*x + v(k)*y))*vis(k))
   *              =   sum_k cos(2*pi*(u(k)*x + v(k)*y))*vis_re(k)
   *                - sum_k sin(2*pi*(u(k)*x + v(k)*y))*vis_im(k)
   *
   * using the trigonometric identities:
   *
   *     cos(ux + vy) = cos(ux)*cos(vy) - sin(ux)*sin(vy)
   *     sin(ux + vy) = cos(ux)*sin(vy) + sin(ux)*cos(vy)
   *
   * to exploit separability and reduce the number of operation and memory
   * footprint.
   */
  if (is_void(vis)) {
    /* Compute "dirty beam", i.e. assume vis(k) = 1 for all k . */
    return (cos_ux)(+,)*(cos_vy)(+,) - (sin_ux)(+,)*(sin_vy)(+,);
  } else {
    re = vis(1,);
    im = vis(2,);
    return (((cos_ux)(+,)*(cos_vy*re)(+,) - (sin_ux)(+,)*(sin_vy*re)(+,)) -
            ((cos_ux)(+,)*(sin_vy*im)(+,) + (sin_ux)(+,)*(cos_vy*im)(+,)));
  }
}

func mira_compute_dirty_beam(master, method=)
/* DOCUMENT psf = mira_compute_dirty_beam(master);

     yields the "dirty beam", that the equivalent point spread function (PSF),
     given by the (u,v) coverage in the interferometric data stored by MASTER.

     If keyword METHOD is not specified or is "xform", the actual pixels to
     complex visibilities transform is used to compute the dirty beam.  If
     keyword METHOD is "exact", the equations are applied with no other
     approximations.

   SEE ALSO: mira_config, mira_compute_dirty_map, mira_compute_residual_map,
     mira_apply_adjoint_of_exact_xform.
*/
{
  /* Apply the adjoint of the pixels to complex visibilities transform to
     complex visibilities all equal to 1.  FIXME: Add zero-th frequency. */
  if (is_void(method) || method == "xform") {
    /* Apply the adjoint of the pixel to visibilities transform. */
    nfreqs = numberof(mira_model_ufreq(master));
    (vis = array(double, 2, nfreqs))(1,) = 1.0;
    psf = master.xform(vis, 1n);
  } else if (method == "exact") {
    /* Apply adjoint of exact transform. */
    psf = mira_apply_adjoint_of_exact_xform(master);
  } else {
    error, "invalid method \""+method+"\" for the dirty beam";
  }

  return psf;
}

func mira_compute_dirty_map(master, method=, tol=)
/* DOCUMENT img = mira_compute_dirty_map(master);

     yields the "dirty map", that the inverse Fourier transform of the complex
     visibility data stored by MASTER.  If there are no measured complex
     visibilities, the returned image is an array of zeros.

     If keyword METHOD is not specified or is "xform", the actual pixel to
     complex visibilities transform is used to compute the dirty beam.  If
     keyword METHOD is "exact", the equations are applied with no other
     approximations.

   SEE ALSO: mira_config, mira_compute_dirty_beam, mira_compute_residual_map,
     mira_apply_adjoint_of_exact_xform.
*/
{
  /* Compute the weighted average of the measured complex visibilities. */
  vis = mira_compute_mean_visibilities(master, tol=tol);
  if (noneof(vis)) {
    /* No valid data. */
    return array(double, mira_image_size(master));
  }

  /* Apply the adjoint of the pixels to complex visibilities transform to the
     measured complex visibilities.  FIXME: Add zero-th frequency. */
  if (is_void(method) || method == "xform") {
    /* Apply the adjoint of the pixel to visibilities transform. */
    img = master.xform(vis, 1n);
  } else if (method == "exact") {
    /* Apply adjoint of exact transform. */
    img = mira_apply_adjoint_of_exact_xform(master, vis);
  } else {
    error, "invalid method \""+method+"\" for the dirty map";
  }

  return img;
}

func mira_compute_residual_map(master, img, method=, tol=)
/* DOCUMENT img = mira_compute_residual_map(master);
         or img = mira_compute_residual_map(master, img);

     yields the "residual map", that the inverse Fourier transform of the
     complex visibility data stored by MASTER minus the model complex
     visibilities resulting from image IMG.  If image IMG is not supplied, the
     last image model set by `mira_update` is assumed.

     If keyword METHOD is not specified or is "xform", the actual pixel to
     complex visibilities transform is used to compute the dirty beam.  If
     keyword METHOD is "exact", the equations are applied with no other
     approximations.

   SEE ALSO: mira_config, mira_update, mira_compute_dirty_beam,
     mira_compute_dirty_map, mira_apply_adjoint_of_exact_xform.
*/
{
  /* Compute the weighted average of the measured complex visibilities. */
  vis = mira_compute_mean_visibilities(master, tol=tol);

  /* Subtract the model complex visibilities.  FIXME: it does not make sense to
     use another transform for that. */
  vis -= (is_void(img) ? mira_model_vis(master) : master.xform(img));

  /* Apply the adjoint of the pixels to complex visibilities transform to the
     measured complex visibilities.  FIXME: Add zero-th frequency. */
  if (is_void(method) || method == "xform") {
    /* Apply the adjoint of the pixel to visibilities transform. */
    img = master.xform(vis, 1n);
  } else if (method == "exact") {
    /* Apply adjoint of exact transform. */
    img = mira_apply_adjoint_of_exact_xform(master, vis);
  } else {
    error, "invalid method \""+method+"\" for the residual map";
  }

  return img;
}

func mira_compute_mean_visibilities(master, tol=)
/* DOCUMENT vis = mira_compute_mean_visibilities(master);

     yields the weighted mean of the complex visibilities measured in MASTER.
     If there are no measured complex visibilities, the result is an array of
     zeros.
*/
{
  /* Prepare the result. */
  nfreqs = numberof(mira_model_ufreq(master));
  vis = array(double, 2, nfreqs);

  /* Collect measured complex visibilities. */
  db = mira_vis_data(master);
  if (is_void(db)) {
    /* No complex visibilities data. */
    return vis;
  }

  /* Extract weights and coordinates. */
  local re, im, wrr, wri, wrr, idx;
  eq_nocopy, re,  db.re;
  eq_nocopy, im,  db.im;
  eq_nocopy, wrr, db.wrr;
  eq_nocopy, wri, db.wri;
  eq_nocopy, wii, db.wii;
  eq_nocopy, idx, db.idx;

  /* Compute the weighted sum of the measured complex visibilities and the sum
     of weights. */
  wre = histogram(idx, wrr*re + wri*im, top=nfreqs);
  wim = histogram(idx, wri*re + wii*im, top=nfreqs);
  wrr = histogram(idx, wrr, top=nfreqs);
  wri = histogram(idx, wri, top=nfreqs);
  wii = histogram(idx, wii, top=nfreqs);

  /* Apply the inverse of the inverse block-diagonal weighting matrix with 2x2
     blocks to the weighted sum of the measured complex visibilities. */
  det = wrr*wii - wri*wri;
  if (is_void(tol)) tol = 1e-20;
  if (tol > 0) {
    nrm = max(abs(wrr), abs(wii));
    i = where(det > tol*nrm*nrm);
  } else {
    i = where(det > 0.0);
  }
  if (is_array(i)) {
    (q = array(double, nfreqs))(i) = 1.0/det(i);
    vis(1,) = q*(wii*wre - wri*wim);
    vis(2,) = q*(wrr*wim - wri*wre);
  }
  return vis;
}
