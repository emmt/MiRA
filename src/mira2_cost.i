/*
 * mira2_cost.i -
 *
 * Objectve functions in MiRA.
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

/*--------------------------------------------------------------------------*/
/* TOPLEVEL ROUTINES */

func mira_ndata(master)
/* DOCUMENT mira_ndata(master);

      yields the number of valid data measurements in MiRA instance `master`.

   SEE ALSO: mira_new.
 */
{
 // FIXME: check stage?
  n = 0;
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    n += db.ops.ndata(master, db);
  }
  return n;
}

func mira_cost(master, x)
{
  /* Update model and integrate cost for each datablock. */
  mira_update, master, x;
  cost = 0.0;
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    cost += db.ops.cost(master, db);
  }
  return cost;
}

func mira_cost_and_gradient(master, x, &grd)
{
  /* Update model and integrate cost and gradient with respect to the complex
     visibilities. */
  mira_update, master, x;
  cost = 0.0;
  grd = array(double, dimsof(mira_model_vis(master)));
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    cost += db.ops.cost(master, db, grd);
  }

  /* Convert gradient with respect to the complex visibilities into gradient
     with respect to pixel values. */
  grd = master.xform(grd, 1);
  return cost;
}

func mira_polar_to_cartesian(amp, amperr, phi, phierr, what,
                             goodman=, use_covar=)
/* DOCUMENT tbl = mira_polar_to_cartesian(amp, amperr, phi, phierr, what)

     Convert complex data given in polar coordinates `(amp,phi)` with their
     standard deviations `(amperr,phierr)` into Cartesian coordinates `(re,im)`
     and associated noise model.  The result is a hash table:

       tbl.re  = real part of complex data
       tbl.im  = imaginary part of complex data
       tbl.wrr = statistical weights for real part of residuals
       tbl.wii = statistical weights for imaginary part of residuals
       tbl.wri = statistical weights for real times imaginary parts of
                 residuals

     The quadratic penalty writes:

       err =     tbl.wrr*(tbl.re - re)^2
             +   tbl.wii*(tbl.im - im)^2
             + 2*tbl.wri*(tbl.re - re)*(tbl.im - im);

     With Goodman approximation, wrr = wii = wgt and wri = 0.  The returned
     table has then the following fields:

       tbl.re  = real part of complex data
       tbl.im  = imaginary part of complex data
       tbl.wgt = statistical weights

   SEE ALSO: mira_cost.
 */
{
  if (min(amp) < 0.0) {
    k = where(amp < 0.0);
    warn, "There are %d negative %s amplitude(s)!!! [FIXED]",
      numberof(k), what;
    phi = phi; // force copy
    phi(k) += MIRA_PI;
    amp = abs(amp);
  }
  cos_phi = cos(phi);
  sin_phi = sin(phi);
  re = amp*cos_phi;
  im = amp*sin_phi;
  tbl = h_new(re=re, im=im);

  /* Deal with valid nonzero complex visibilities. */
  k = where((amperr > 0.0) & (phierr > 0.0) & (amp > 0.0));
  if (goodman) {
    /*
     * Use Goodman approximation with SIGMA such that the area of ellipsoids at
     * one standard deviation are the same (this also amount to having the same
     * determinant of the covariance matrices):
     *
     *     PI*SIGMA^2 = PI*(RHO*SIGMA_THETA)*SIGMA_RHO
     *
     * hence:
     *
     *     SIGMA = sqrt(RHO*SIGMA_THETA*SIGMA_RHO)
     *
     * If SIGMA_THETA or SIGMA_RHO is invalid, take:
     *
     *     SIGMA = SIGMA_RHO
     *     SIGMA = RHO*SIGMA_THETA
     *
     * accordingly.
     */
    wgt = array(double, dimsof(re));
    h_set, tbl, goodman=1n, wgt=wgt;
    if (is_array(k)) {
      wgt(k) = 1.0/(amp(k)*phierr(k)*amperr(k));
    }
  } else {
    /*
     * Use convex quadratic local approximation.
     *
     * We consider complex data in polar form: ρ⋅exp(i⋅ϕ)
     * and a complex model: z.
     *
     * The amplitude and phase standard deviations directly yield the
     * covariance anti-residuals in a frame whose first axis is parallel to
     * the complex data vector:
     *
     *     err = err_para⋅u_para + err_perp⋅u_perp
     *
     * with u_para = [cosϕ, sinϕ] and v_para = [-sinϕ, cosϕ] the unit vectors
     * parallel and perpendicular to the data vector.  Assuming independent
     * errors along these 2 axes, the covariance writes:
     *
     *     Cov(err) = var_para⋅u_para⋅u_para' + var_perp⋅u_perp⋅u_perp'
     *              = [cov_rr, cov_ri; cov_ri, cov_im]
     *
     * with:
     *
     *     cov_rr = cos²ϕ⋅var_para + sin²ϕ⋅var_perp,
     *     cov_ri = cosϕ⋅sinϕ⋅(var_para - var_perp),
     *     cov_ii = cos²ϕ⋅var_perp + sin²ϕ⋅var_para,
     *
     * Remarks: (i) det(Cov(err)) = var_para⋅var_perp; (ii) as expected, for a
     * positive angle ϕ ≤ 90°, Re(err) and Im(err) are correlated if var_para >
     * var_perp, uncorrelated if var_para = var_perp and anti-correlated if
     * var_para < var_perp.
     */
    wrr = wri = wii = array(double, dimsof(re));
    h_set, tbl, goodman=0n, wrr=wrr, wri=wri, wii=wii;
    if (is_array(k)) {
      std_para = amperr(k);
      std_perp = amp(k)*phierr(k);
      cs = cos_phi(k);
      sn = sin_phi(k);
      cs_cs = cs*cs;
      sn_sn = sn*sn;
      if (use_covar) {
        /* Using the covariances, it is possible to have 1 division instead of
           2 but 2 more multiplications and a negation. */
        var_para = std_para*std_para;
        var_perp = std_perp*std_perp;
        cov_rr = cs_cs*var_para + sn_sn*var_perp;
        cov_ri = cs*sn*(var_para - var_perp);
        cov_ii = cs_cs*var_perp + sn_sn*var_para;
        one_over_det = 1.0/(var_para*var_perp);
        wrr(k) =  one_over_det*cov_ii;
        wri(k) = -one_over_det*cov_ri;
        wii(k) =  one_over_det*cov_rr;
      } else {
        wgt_para = 1.0/(std_para*std_para);
        wgt_perp = 1.0/(std_perp*std_perp);
        wrr(k) = cs_cs*wgt_para + sn_sn*wgt_perp;
        wri(k) = cs*sn*(wgt_para - wgt_perp);
        wii(k) = cs_cs*wgt_perp + sn_sn*wgt_para;
      }
    }
  }

  /* Deal with valid zero complex visibilities. */
  k = where((amperr > 0.0) & (amp == 0.0));
  if (is_array(k)) {
    std_para = amperr(k);
    wgt_para = 1.0/(std_para*std_para);
    if (goodman) {
      tbl.wgt(k) = wgt_para;
    } else {
      tbl.wrr(k) = wgt_para;
      tbl.wii(k) = wgt_para;
    }
  }

  return tbl;
}

/*---------------------------------------------------------------------------*/
/* POWERSPECTRUM DATA */

func _mira_vis2_ndata(master, db)
{
  return numberof(db.idx);
}

func _mira_vis2_cost(master, db, grd)
{
  /* Error and gradient w.r.t. powerspectrum data:
   *           COST = sum(WGT*ERR^2)
   *     ∂COST/∂MDL = 2*WGT*ERR
   * where ERR = MDL - DAT. The convertion of the gradient
   * writes:
   *     ∂COST/∂RE = (∂MDL/∂RE)⋅(∂COST/∂MDL) = (2⋅RE)⋅(2⋅WGT⋅ERR)
   *     ∂COST/∂IM = (∂MDL/∂IM)⋅(∂COST/∂MDL) = (2⋅IM)⋅(2⋅WGT⋅ERR)
   * finally:
   *     ∂COST/∂RE = 4⋅WGT⋅ERR⋅RE
   *     ∂COST/∂IM = 4⋅WGT⋅ERR⋅IM
   */
  local idx; eq_nocopy, idx, db.idx;
  re = mira_model_re(master, idx);
  im = mira_model_im(master, idx);
  err = (re*re + im*im) - db.dat;
  wgt_err = db.wgt*err;
  cost = sum(wgt_err*err);
  if (! is_void(grd)) {
    fct = 4.0*unref(wgt_err);
    _mira_update_gradient, grd, idx, fct*re, fct*im;
  }
  return cost;
}

/*---------------------------------------------------------------------------*/
/* COMPLEX VISIBILITY DATA */

func _mira_vis_ndata(master, db)
{
  return 2*numberof(db.idx);
}

func _mira_visamp_ndata(master, db)
{
  return numberof(db.idx);
}

func _mira_visphi_ndata(master, db)
{
  return numberof(db.idx);
}

func _mira_vis_cost(master, db, grd)
{
  local idx; eq_nocopy, idx, db.idx;
  return _mira_visibility_cost(master.flags,
                               mira_model_re(master, idx),
                               mira_model_im(master, idx),
                               db.wrr, db.wri, db.wii,
                               db.re, db.im,
                               idx, grd);
}

func _mira_visamp_cost(master, db, grd)
{
  local idx; eq_nocopy, idx, db.idx;
  return _mira_amplitude_cost(master.flags,
                              mira_model_re(master, idx),
                              mira_model_im(master, idx),
                              mira_model_amp(master, idx),
                              db.wgt, db.dat, idx, grd);
}

func _mira_visphi_cost(master, db, grd)
{
  local idx; eq_nocopy, idx, db.idx;
  return _mira_phase_cost(master.flags,
                          mira_model_re(master, idx),
                          mira_model_im(master, idx),
                          mira_model_vis2(master, idx),
                          mira_model_phi(master, idx),
                          db.wgt, db.dat, idx, grd);
}

/*---------------------------------------------------------------------------*/
/* BISPECTRUM DATA */

func _mira_t3_ndata(master, db)
{
  return 2*numberof(db.idx)/3;
}

func _mira_t3amp_ndata(master, db)
{
  return numberof(db.idx)/3;
}

func _mira_t3phi_ndata(master, db)
{
  return numberof(db.idx)/3;
}

func _mira_t3_cost(master, db, grd)
{
  /* Declare variables shared by pre/post helper subroutines. */
  local idx, sgn, re, im;
  local t3re, t3im, t3idx, t3grd;
  local re1, im1, re2, im2, re3, im3;
  local re12, im12;

  /* Initialize shared variables. */
  _mira_t3_pre_cost, master, db;

  /* Compute the cost (and its gradient) with respect to the bispectrum. */
  cost = _mira_visibility_cost(master.flags, t3re, t3im,
                               db.wrr, db.wri, db.wii,
                               db.re, db.im,
                               t3idx, t3grd);

  /* If needed, update gradient w.r.t. the model complex visibilities. */
  if (! is_void(grd)) {
    _mira_t3_post_cost, grd;
  }
  return cost;
}

func _mira_t3amp_cost(master, db, grd)
{
  /* Declare variables shared by pre/post helper subroutines. */
  local idx, sgn, re, im;
  local t3re, t3im, t3idx, t3grd;
  local re1, im1, re2, im2, re3, im3;
  local re12, im12;

  /* Initialize shared variables. */
  _mira_t3_pre_cost, master, db;

  /* Compute the model of the bispectrum amplitude. */
  t3amp = abs(t3re, t3im);

  /* Compute the cost (and its gradient) with respect to the bispectrum. */
  cost = _mira_amplitude_cost(master.flags, t3re, t3im, t3amp,
                              db.wgt, db.dat, t3idx, t3grd);

  /* If needed, update gradient w.r.t. the model complex visibilities. */
  if (! is_void(grd)) {
    _mira_t3_post_cost, grd;
  }
  return cost;
}

func _mira_t3phi_cost(master, db, grd)
{
  /* Declare variables shared by pre/post helper subroutines. */
  local idx, sgn, re, im;
  local t3re, t3im, t3idx, t3grd;
  local re1, im1, re2, im2, re3, im3;
  local re12, im12;

  /* Initialize shared variables. */
  _mira_t3_pre_cost, master, db;

  /* Compute the model of the phase closures. */
  t3phi = sgn(+,)*mira_model_phi(master, idx)(+,);

  /* Compute the cost (and its gradient) with respect to the bispectrum. */
  cost = _mira_phase_cost(master.flags, t3re, t3im,
                          t3re*t3re + t3im*t3im, t3phi,
                          db.wgt, db.dat, t3idx, t3grd);

  /* If needed, update gradient w.r.t. the model complex visibilities. */
  if (! is_void(grd)) {
    _mira_t3_post_cost, grd;
  }
  return cost;
}

func _mira_t3_pre_cost(master, db)
{
  extern idx, sgn, re, im;
  extern t3re, t3im, t3idx, t3grd;
  extern re1, im1, re2, im2, re3, im3;
  extern re12, im12;
  eq_nocopy, idx, db.idx;
  eq_nocopy, sgn, db.sgn;
  re = mira_model_re(master, idx);
  im = mira_model_im(master, idx)*sgn;
  re1 = re(1,);
  im1 = im(1,);
  re2 = re(2,);
  im2 = im(2,);
  re3 = re(3,);
  im3 = im(3,);
  re12 = mira_cmult_re(re1,im1, re2,im2);
  im12 = mira_cmult_im(re1,im1, re2,im2);
  t3re = mira_cmult_re(re12,im12, re3,im3);
  t3im = mira_cmult_im(re12,im12, re3,im3);
  if (is_void(grd)) {
    t3idx = [];
    t3grd = [];
  } else {
    n = numberof(t3re);
    t3idx = indgen(n);
    t3grd = array(double, 2, n);
  }
}

func  _mira_t3_post_cost(grd)
/* DOCUMENT _mira_t3_post_cost, grd;

     Private subroutine to update the gradient w.r.t. the model complex
     visibilities, stored in GRD, given G the gradient w.r.t. the complex
     bispectrum.  RE and IM are the real and imaginary parts of the model
     complex visibilities involved in the bispectrum.  SGN is sign of the
     phase (or of the imaginary part) of the involved model complex
     visibilities relative to the global model visibilities (IM is already
     multiplied by SGN).  IDX gives the indices of the involved model complex
     visibilities in the global model visibilities.

     All arguments must be 2D arrays, GRD and G have a leading dimension of 2
     (for the real and imaginary parts respectively), all the others have a
     leading dimension of 3 (for the 3 involved complex visibilities in each
     bispectrum).

   SEE ALSO: _mira_t3_cost, _mira_t3amp_cost, _mira_t3phi_cost.
 */
{
  extern idx, sgn, re, im;
  extern t3re, t3im, t3idx, t3grd;
  extern re1, im1, re2, im2, re3, im3;
  extern re12, im12;
  t3grd_re = t3grd(1,);
  t3grd_im = t3grd(2,);
  _mira_t3_post_cost_helper, idx(1,), sgn(1,),
    mira_cmult_re(re2,im2, re3,im3),
    mira_cmult_im(re2,im2, re3,im3);
  _mira_t3_post_cost_helper, idx(2,), sgn(2,),
    mira_cmult_re(re1,im1, re3,im3),
    mira_cmult_im(re1,im1, re3,im3);
  _mira_t3_post_cost_helper, idx(3,), sgn(3,), re12, im12;
}

func _mira_t3_post_cost_helper(idx, sgn, re, im)
{
  /*
   * For f: ℂ ↦ ℝ, the gradient can be written as:
   *
   *     f'(z) = ∂f(z)/∂z = ∂f(z)/∂Re(z) + i⋅∂f(z)/∂Im(z)
   *
   * where z ∈ ℂ, then for any w ∈ ℂ:
   *
   *     ∂f(w⋅z)/∂z = conj(w)⋅f'(w.z)
   *
   * In the code below:
   *
   *     f'(w⋅z) = t3grd_re + i⋅t3grd_im
   *     w = re + i⋅im
   */
  extern grd, t3grd_re, t3grd_im;
  _mira_update_gradient, grd, idx,
    (re*t3grd_re + im*t3grd_im),
    (re*t3grd_im - im*t3grd_re)*sgn;
}

/*--------------------------------------------------------------------------*/
/* COST FUNCTIONS */

func _mira_visibility_cost(flags, mdl_re, mdl_im, wgt_rr, wgt_ri, wgt_ii,
                           dat_re, dat_im, idx, grd)
/* DOCUMENT _mira_visibility_cost(flags, mdl_re, mdl_im,
                                  wgt_rr, wgt_ri, wgt_ii,
                                  dat_re, dat_im, idx, grd);

     yields the misfit error of complex visibility data.  Arguments MDL_RE, and
     MDL_IM are the real and imaginary parts of the model complex visibilities
     corresponding to the data.  WGT_RR, WGT_RI and WGT_II gives the weights of
     the complex data whose real and imaginary parts are DAT_RE and DAT_IM.
     Optional arguments IDX and GRD gives the indices of the model complex
     visibilities corresponding to the data and the array to integrate the
     gradient of the misfit error with respect to the model complex amplitudes.

   SEE ALSO: mira_polar_to_cartesian.
 */
{
  /*
   * The real and imaginary parts of the (anti-)residuals are respectively
   * given by:
   *
   *     ERR_RE = MDL_RE - DAT_RE
   *     ERR_IM = MDL_IM - DAT_IM
   *
   * then:
   *
   *     COST = sum(WGT_RR⋅ERR_RE^2 + 2⋅WGT_RI⋅ERR_RE⋅ERR_IM +
   *                WGT_II⋅ERR_IM^2)
   *
   * can be rewritten as:
   *
   *     COST = sum((WGT_RR⋅ERR_RE + WGT_RI⋅ERR_IM)⋅ERR_RE +
   *                (WGT_RI⋅ERR_RE + WGT_II⋅ERR_IM)⋅ERR_IM)
   *          = (1/2)⋅sum((∂COST/∂MDL_RE)⋅ERR_RE) +
   *            (1/2)⋅sum((∂COST/∂MDL_IM)⋅ERR_IM)
   *
   * where:
   *
   *     ∂COST/∂MDL_RE = 2⋅(WGT_RR⋅ERR_RE + WGT_RI⋅ERR_IM)
   *     ∂COST/∂MDL_IM = 2⋅(WGT_RI⋅ERR_RE + WGT_II⋅ERR_IM)
   *
   * Introducing common expressions:
   *
   *     TMP_RE = WGT_RR⋅ERR_RE + WGT_RI⋅ERR_IM
   *     TMP_IM = WGT_RI⋅ERR_RE + WGT_II⋅ERR_IM
   *
   * the cost and gradient simplifies to:
   *
   *     COST = sum(TMP_RE⋅ERR_RE) + sum(TMP_IM⋅ERR_IM)  (*)
   *     ∂COST/∂MDL_RE = 2⋅TMP_RE
   *     ∂COST/∂MDL_IM = 2⋅TMP_IM
   *
   * (*) sum(A)+sum(B) where A and B are arrays involves as many operations
   *      as sum(A+B) but does not require temporary storage of A+B.
   */
  err_re = mdl_re - dat_re;
  err_im = mdl_im - dat_im;
  tmp_re = wgt_rr*err_re + wgt_ri*err_im;
  tmp_im = wgt_ri*err_re + wgt_ii*err_im;
  cost = sum(tmp_re*err_re) + sum(tmp_im*err_im);
  if (! is_void(grd)) {
    _mira_update_gradient, grd, idx, tmp_re + tmp_re, tmp_im + tmp_im;
  }
  return cost;
}

func _mira_amplitude_cost(flags, re, im, mdl, wgt, dat, idx, grd)
/* DOCUMENT _mira_amplitude_cost(flags, re, im, mdl, wgt, dat);
         or _mira_amplitude_cost(flags, re, im, mdl, wgt, dat, idx, grd);

     yields the misfit error of amplitude-only data DAT.  Arguments RE, IM and
     MDL are the real part, imaginary part and amplitude of the model complex
     visibilities corresponding to the data.  WGT gives the weights of the
     data.  Optional arguments IDX and GRD gives the indices of the model
     complex visibilities corresponding to the data and the array to integrate
     the gradient of the misfit error with respect to the model complex
     amplitudes.  It is assumed that MDL=abs(RE,IM).

     For amplitude-only data, the misfit error writes:

         COST = sum(WGT⋅ERR²)

     with:

         MDL = sqrt(RE² + IM²)
         ERR = MDL - DAT

     hence the gradient is given by:

         ∂COST/∂RE = FCT⋅RE
         ∂COST/∂IM = FCT⋅IM

     with:

         FCT = 2⋅WGT⋅(sqrt(RE² + IM²) - MDL)/sqrt(RE² + IM²)
             = 2⋅WGT⋅ERR/MDL

     Beware that the case MDL = 0 has to be dealt specifically (only the error
     can be computed, not the gradient).

   SEE ALSO: _mira_phase_cost.
 */
{
  if (min(mdl) > 0) {
    /* All model complex visibilities are nonzero. */
    cost = 0.0;
  } else {
    /* Deal with model visibilities equal to zero (where only the error can be
       computed, not the gradient) and only keep the values where the
       model complex visibilities are nonzero. */
    k = where(! mdl);
    tmp = dat(k);
    cost = sum(wgt(k)*tmp*tmp);
    k = where(mdl);
    if (! is_array(k)) {
      /* There are no nonzero complex visibilities. */
      return cost;
    }
    mdl = mdl(k);
    wgt = wgt(k);
    dat = dat(k);
    if (! is_void(grd)) {
      re = re(k);
      im = im(k);
      idx = idx(k);
    }
  }

  /* Penalty (and gradient) for non-zero model visibilities. */
  err = mdl - dat;
  wgt_err = wgt*res;
  cost += sum(wgt_err*err);
  if (! is_void(grd)) {
    fct = (wgt_err + wgt_err)/mdl;
    _mira_update_gradient, grd, idx, fct*re, fct*im;
  }
  return cost;
}

func _mira_phase_cost(flags, re, im, vis2, mdl, wgt, dat, idx, grd)
/* DOCUMENT _mira_phase_cost(flags, re, im, vis2, mdl, wgt, dat);
         or _mira_phase_cost(flags, re, im, vis2, mdl, wgt, dat, idx, grd);

     yields the misfit error of phase-only data DAT.  Arguments RE, IM, VIS2
     and MDL are the real part, imaginary part, squared amplitude and phase of
     the model complex visibilities corresponding to the data.  WGT gives the
     weights of the data.  Optional arguments IDX and GRD gives the indices of
     the model complex visibilities corresponding to the data and the array to
     integrate the gradient of the misfit error with respect to the model
     complex phases.  It is assumed that VIS2=(RE^2+IM^2) and that
     MDL=atan(IM,RE).

   SEE ALSO: _mira_amplitude_cost.
 */
{
  /* Keep only bits corresponding to the phase cost approximation. */
  bits = (flags & (MIRA_HANIFF_APPROX | MIRA_CONVEX_LIMIT |
                   MIRA_VON_MISES_APPROX));

  /* Deal with model visibilities equal to zero (where model phases are
     unknown). */
  if (min(vis2) > 0.0) {
    /* All model complex visibilities are nonzero. */
    cost = 0.0;
  } else {
    /* Set the cost to be maximal where model visibilities equal to zero and
       only keep the values corresponding to the nonzero complex
       visibilities. */
    cost = sum(wgt(where(! vis2)));
    if (bits == MIRA_VON_MISES_APPROX) {
      cost *= 4.0;
    } else if (bits == MIRA_HANIFF_APPROX) {
      cost *= MIRA_PI^2;
    }
    k = where(vis2);
    if (! is_array(k)) {
      /* There are no nonzero complex visibilities. */
      return cost;
    }
    mdl = mdl(k);
    wgt = wgt(k);
    dat = dat(k);
    if (! is_void(grd)) {
      re = re(k);
      im = im(k);
      vis2 = vis2(k);
      idx = idx(k);
    }
  }

  /* Penalty (and gradient) for non-zero model visibilities. */
  if (bits == MIRA_VON_MISES_APPROX) {
    /* With this option, we use the length of the cord to compute penalty (this
     * also corresponds to the Von Mises distribution a.k.a. circular normal
     * distribution -- http://en.wikipedia.org/wiki/Von_Mises_distribution):
     *
     *   COST = 2⋅sum(WGT⋅(1 - cos(ERR)))
     *        = 4⋅sum(WGT⋅sin(ERR/2)^2)             (*)
     *
     * whose gradient w.r.t. the model is:
     *
     *   ∂COST/∂MDL = 4⋅WGT⋅sin(ERR/2)⋅cos(ERR/2)
     *              = 2⋅WGT⋅sin(ERR)                (*)
     *
     * where ERR = MDL - DAT is the phase (or closure) error and WGT the
     * weight.  Maybe 1 - cos(ERR) is more subject to rounding errors than
     * sin(ERR/2)^2 especially for small residual errors.  For that reason, we
     * use the equations marked with an asterisk.
     *
     * Finally:
     *
     *     ∂COST/∂RE = -FCT⋅IM
     *     ∂COST/∂IM = +FCT⋅RE
     *
     * with:
     *
     *     FCT = (∂COST/∂MDL)/(RE² + IM²)
     */
    err = mdl - dat;
    tmp = sin(0.5*err);
    cost += 4.0*sum(wgt*tmp*tmp);
    if (! is_void(grd)) {
      fct = (wgt + wgt)*sin(err)/vis2;
      _mira_update_gradient, grd, idx, -fct*im, fct*re;
    }
  } else if (bits == MIRA_HANIFF_APPROX) {
    /*
     * With this option and using the same notation as above, the cost is:
     *
     *   COST = sum(WGT⋅arc(ERR)^2)
     *
     * whose gradient w.r.t. the model is:
     *
     *   ∂COST/∂MDL = 2⋅WGT⋅arc(ERR)
     */
    arc_err = arc(mdl - dat);
    wgt_arc_err = wgt*arc_err;
    cost += sum(arc_err*wgt_arc_err);
    if (! is_void(grd)) {
      fct = (wgt_arc_err + wgt_arc_err)/vis2;
      _mira_update_gradient, grd, idx, -fct*im, fct*re;
    }
  } else if (bits == MIRA_CONVEX_LIMIT) {
    /*
     * With this option and using the same notation as above, the cost is:
     *
     *   COST = sum(WGT⋅sin(ERR)^2)
     *
     * whose gradient w.r.t. the model is:
     *
     *   ∂COST/∂MDL = 2⋅WGT⋅sin(ERR)⋅cos(ERR)   (*)
     *              = WGT⋅sin(2⋅ERR)
     *
     * with some compilers (not Yorick!) it may be faster to compute the sine
     * and the cosine of the same angle (in which case the marked equation is
     * to be preferred).
     */
    err = mdl - dat;
    sin_err = sin(err);
    cost += sum(wgt*sin_err*sin_err);
    if (! is_void(grd)) {
      fct = wgt*sin(err + err)/vis2;
      _mira_update_gradient, grd, idx, -fct*im, fct*re;
    }
  } else {
    error, "unknown choice for phase cost";
  }
  return cost;
}

func _mira_update_gradient(grd, idx, grd_re, grd_im, fast)
{
  if (fast) {
    grd(1,idx) += grd_re;
    grd(2,idx) += grd_im;
  } else {
    n = numberof(grd)/2;
    grd(1,*) += histogram(idx, grd_re, top=n);
    grd(2,*) += histogram(idx, grd_im, top=n);
  }
}

/*---------------------------------------------------------------------------*/
