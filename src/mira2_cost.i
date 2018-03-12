/*
 * mira2_cost.i -
 *
 * Objective functions in MiRA.
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

/*--------------------------------------------------------------------------*/
/* TOPLEVEL ROUTINES */

func mira_ndata(master)
/* DOCUMENT mira_ndata(master);

      yields the number of valid data measurements in MiRA instance `master`.

   SEE ALSO: mira_new.
 */
{
  if (master.stage < 1) {
    _mira_select_data, master;
  }
  ndata = 0;
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    ndata += db.ops.ndata(master, db);
  }
  return ndata;
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
  dims = dimsof(mira_model_vis(master));
  if (numberof(dims) != 3 || dims(2) != 2) {
    throw, "unexpected dimensions for model complex visibilities";
  }
  tbl = h_new(nfreqs = dims(3));
  cost = 0.0;
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    cost += db.ops.cost(master, db, tbl);
  }

  /* Convert gradient with respect to the complex visibilities into gradient
     with respect to pixel values. */
  increment = _mira_increment_variable;
  grd_re = h_pop(tbl, "re");
  grd_im = h_pop(tbl, "im");
  grd_amp = h_pop(tbl, "amp");
  if (! is_void(grd_amp)) {
    fct = mira_model_reciprocal_amp(master)*grd_amp;
    increment, grd_re, fct*mira_model_re(master);
    increment, grd_im, fct*mira_model_im(master);
  }
  grd_vis2 = h_pop(tbl, "vis2");
  if (! is_void(grd_vis2)) {
    fct = grd_vis2 + grd_vis2;
    increment, grd_re, fct*mira_model_re(master);
    increment, grd_im, fct*mira_model_im(master);
  }
  grd_phi = h_pop(tbl, "phi");
  if (! is_void(grd_phi)) {
    fct = mira_model_reciprocal_vis2(master)*grd_phi;
    increment, grd_re, -fct*mira_model_im(master);
    increment, grd_im, +fct*mira_model_re(master);
  }
  grd = array(double, dims);
  if (! is_void(grd_re)) {
    grd(1,) = grd_re;
  }
  if (! is_void(grd_im)) {
    grd(2,) = grd_im;
  }
  grd = master.xform(grd, 1);
  return cost;
}

func _mira_update_gradient(tbl, key, idx, val)
{
  if (!is_void(idx)) {
    val = histogram(idx, val, top=tbl.nfreqs);
  }
  if (h_has(tbl, key)) {
    h_set, tbl, key, h_get(tbl, key) + val;
  } else {
    h_set, tbl, key, val;
  }
}

func _mira_increment_variable(&var, inc)
{
  if (is_void(var)) {
    eq_nocopy, var, inc;
  } else {
    var = unref(var) + inc;
  }
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

func _mira_vis2_cost(master, db, tbl)
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
  if (is_hash(tbl)) {
    _mira_update_gradient, tbl, "vis2", idx, wgt_err + wgt_err;
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

func _mira_vis_cost(master, db, tbl)
{
  local idx;
  eq_nocopy, idx, db.idx;
  grd_re = grd_im = compute_gradient = is_hash(tbl);
  cost = _mira_visibility_cost(master.flags,
                               mira_model_re(master, idx),
                               mira_model_im(master, idx),
                               db.wrr, db.wri, db.wii,
                               db.re, db.im, grd_re, grd_im);
  if (compute_gradient) {
    _mira_update_gradient, tbl, "re", idx, grd_re;
    _mira_update_gradient, tbl, "im", idx, grd_im;
  }
  return cost;
}

func _mira_visamp_cost(master, db, tbl)
{
  local idx;
  eq_nocopy, idx, db.idx;
  grd = compute_gradient = is_hash(tbl);
  cost = _mira_amplitude_cost(master.flags,
                              mira_model_amp(master, idx),
                              db.wgt, db.dat, grd);
  if (compute_gradient) {
    _mira_update_gradient, tbl, "amp", idx, grd;
  }
  return cost;
}

func _mira_visphi_cost(master, db, tbl)
{
  local idx;
  eq_nocopy, idx, db.idx;
  grd = compute_gradient = is_hash(tbl);
  cost = _mira_phase_cost(master.flags,
                          mira_model_phi(master, idx),
                          db.wgt, db.dat, grd);
  if (compute_gradient) {
    _mira_update_gradient, tbl, "phi", idx, grd;
  }
  return cost;
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

func _mira_t3_cost(master, db, tbl)
{
  /* Retrieve indices and respective imaginary signs of real and imaginary
     parts of the complex visibilities involved in the measured complex
     bispectrum. */
  local idx, sgn;
  eq_nocopy, idx, db.idx;
  eq_nocopy, sgn, db.sgn;

  /* Retrieve real and imaginary parts of the complex visibilities involved in
     the measured complex bispectrum and compute the model bispectrum. */
  re = mira_model_re(master, idx);
  im = mira_model_im(master, idx)*sgn;
  re1 = re(1,);
  im1 = im(1,);
  re2 = re(2,);
  im2 = im(2,);
  re3 = re(3,); re = [];
  im3 = im(3,); im = [];
  re12 = mira_cmult_re(re1,im1, re2,im2);
  im12 = mira_cmult_im(re1,im1, re2,im2);
  t3mdl_re = mira_cmult_re(re12,im12, re3,im3);
  t3mdl_im = mira_cmult_im(re12,im12, re3,im3);

  /* Compute the cost (and its gradient) with respect to the bispectrum. */
  t3grd_re = t3grd_im = compute_gradient = is_hash(tbl);
  cost = _mira_visibility_cost(master.flags, t3mdl_re, t3mdl_im,
                               db.wrr, db.wri, db.wii,
                               db.re, db.im, t3grd_re, t3grd_im);

  /* If needed, update gradient w.r.t. the model complex visibilities. */
  if (compute_gradient) {
    _mira_t3_post_cost_helper, idx(1,), sgn(1,),
      mira_cmult_re(re2,im2, re3,im3),
      mira_cmult_im(re2,im2, re3,im3);
    _mira_t3_post_cost_helper, idx(2,), sgn(2,),
      mira_cmult_re(re1,im1, re3,im3),
      mira_cmult_im(re1,im1, re3,im3);
    _mira_t3_post_cost_helper, idx(3,), sgn(3,), re12, im12;
  }
  return cost;
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
  extern tbl, t3grd_re, t3grd_im;
  _mira_update_gradient, tbl, "re", idx, (re*t3grd_re + im*t3grd_im);
  _mira_update_gradient, tbl, "im", idx, (re*t3grd_im - im*t3grd_re)*sgn;
}

func _mira_t3amp_cost(master, db, tbl)
{
  /* Retrieve indices of amplitudes involved in the measured bispectrum
     amplitudes. */
  local idx;
  eq_nocopy, idx, db.idx;

  /* Compute the model of the bispectrum amplitude. */
  amp = mira_model_amp(master, idx);
  mdl = amp(1,)*amp(2,)*amp(3,);
  amp = [];

  /* Compute the cost (and its gradient) with respect to the bispectrum. */
  grd = compute_gradient = is_hash(tbl);
  cost = _mira_amplitude_cost(master.flags, mdl, db.wgt, db.dat, grd);

  /* If needed, update gradient. */
  if (compute_gradient) {
    local reciprocal_vis2;
    eq_nocopy, reciprocal_vis2, mira_model_reciprocal_vis2(master);
    n = numberof(reciprocal_vis2);
    fct = reciprocal_vis2*histogram(idx, (mdl*grd)(-:1:3,), top=n);
    _mira_update_gradient, tbl, "re", [], fct*mira_model_re(master);
    _mira_update_gradient, tbl, "im", [], fct*mira_model_im(master);
  }
  return cost;
}

func _mira_t3phi_cost(master, db, tbl)
{
  /* Retrieve indices and respective signs of phases involved in the measured
     phase closures. */
  local idx, sgn;
  eq_nocopy, idx, db.idx;
  eq_nocopy, sgn, db.sgn;

  /* Compute the cost (and its gradient) with respect to the phase closures
     given the model of the phase closures. */
  grd = compute_gradient = is_hash(tbl);
  cost = _mira_phase_cost(master.flags,
                          (sgn*mira_model_phi(master, idx))(sum,),
                          db.wgt, db.dat, grd);

  /* If needed, update gradient.  */
  if (compute_gradient) {
    _mira_update_gradient, tbl, "phi", idx, sgn*grd(-,);
  }
  return cost;
}

/*--------------------------------------------------------------------------*/
/* COST FUNCTIONS */

func _mira_visibility_cost(flags, mdl_re, mdl_im, wgt_rr, wgt_ri, wgt_ii,
                           dat_re, dat_im, &grd_re, &grd_im)
/* DOCUMENT fdata = _mira_visibility_cost(flags, mdl_re, mdl_im,
                                          wgt_rr, wgt_ri, wgt_ii,
                                          dat_re, dat_im);
         or fdata = _mira_visibility_cost(flags, mdl_re, mdl_im,
                                          wgt_rr, wgt_ri, wgt_ii,
                                          dat_re, dat_im, grd_re, grd_im);

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
  if (grd_re) {
    grd_re = tmp_re + tmp_re;
  }
  if (grd_im) {
    grd_im = tmp_im + tmp_im;
  }
  return cost;
}

func _mira_amplitude_cost(flags, mdl, wgt, dat, &grd)
/* DOCUMENT fdata = _mira_amplitude_cost(flags, mdl, wgt, dat);
         or fdata = _mira_amplitude_cost(flags, mdl, wgt, dat, grd);

     yields the data fidelity cost for amplitude-only data.  FLAGS is a bitwise
     combination of options to choose the expression of the data fidelity (for
     now only Gaussian approximation is supported).  MDL is the model value(s)
     of the amplitude(s), DAT specifies the amplitude data and WGT their
     respective weights.  If GRD is given, it must be a caller's variable set
     to true on entry and is replaced by the gradient of the fidelity function
     with respect to the amplitude model:

         GRD = ∂COST/∂MDL

     Then, if MDL = sqrt(RE² + IM²), the gradients with respect to RE and IM
     write:

         ∂COST/∂RE = GRD⋅RE/MDL = cos(PHI)⋅GRD
         ∂COST/∂IM = GRD⋅IM/MDL = sin(PHI)⋅GRD

     with PHI = atan(IM, RE).


   SEE ALSO: mira_cost, _mira_phase_cost, _mira_visibility_cost.
 */
{
  err = mdl - dat;
  wgt_err = wgt*err;
  cost = sum(wgt_err*err);
  if (grd) {
    grd = wgt_err + wgt_err;
  }
  return cost;
}

func _mira_phase_cost(flags, mdl, wgt, dat, &grd)
/* DOCUMENT fdata = _mira_phase_cost(flags, mdl, wgt, dat);
         or fdata = _mira_phase_cost(flags, mdl, wgt, dat, grd);

     yields the data fidelity cost for phase-only data.  FLAGS is a bitwise
     combination of options to choose the expression of the data fidelity
     (MIRA_VON_MISES_APPROX, MIRA_HANIFF_APPROX and MIRA_CONVEX_LIMIT are
     supported).  MDL is the model value(s) of the phase(s), DAT specifies the
     phase data and WGT their respective weights.  If GRD is given, it must be
     a caller's variable set to true on entry and is replaced by the gradient
     of the fidelity function with respect to the phase model:

         GRD = ∂COST/∂MDL

     Then, if MDL = atan(IM, RE), the gradients with respect to RE and IM
     write:

         ∂COST/∂RE = -FCT⋅IM
         ∂COST/∂IM = +FCT⋅RE

     with:

         FCT = (∂COST/∂MDL)/(RE² + IM²)


   SEE ALSO: mira_cost, _mira_amplitude_cost, _mira_visibility_cost.
 */
{
  /* Keep only bits corresponding to the phase cost approximation. */
  bits = (flags & _MIRA_PHASE_ONLY_BITS);

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
     */
    err = mdl - dat;
    tmp = sin(0.5*err);
    cost = 4.0*sum(wgt*tmp*tmp);
    if (grd) {
      grd = (wgt + wgt)*sin(err);
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
    cost = sum(arc_err*wgt_arc_err);
    if (! is_void(grd)) {
      grd = wgt_arc_err + wgt_arc_err;
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
    cost = sum(wgt*sin_err*sin_err);
    if (grd) {
      grd = wgt*sin(err + err);
    }
  } else {
    error, "unknown choice for phase cost";
  }
  return cost;
}

/*---------------------------------------------------------------------------*/
