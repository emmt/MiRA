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
  for (db = _mira_first(master); db; _mira_next(master, db)) {
    n += db.ops.ndata(master, db);
  }
  return n;
}

func mira_cost(master, x)
{
  /* Update model and integrate cost for each datablock. */
  _mira_update_model, master, x;
  cost = 0.0;
  for (db = _mira_first(master); db; _mira_next(master, db)) {
    cost += db.ops.cost(master, db);
  }
  return cost;
}

func mira_cost_and_gradient(master, x, &grd)
{
  /* Update model and integrate cost and gradient with respect to the complex
     visibilities. */
  _mira_update_model, master, x;
  cost = 0.0;
  grd = array(double, dimsof(mira_model_vis(master)));
  for (db = _mira_first(master); db; _mira_next(master, db)) {
    cost += db.ops.cost(master, db, grd);
  }

  /* Convert gradient with respect to the complex visibilities into gradient
     with respect to pixel values. */
  grd = master.xform(grd, 1);
  return cost;
}

func mira_polar_to_cartesian(amp, amperr, phi, phierr, what, goodman=)
/* DOCUMENT ws = mira_polar_to_cartesian(amp, amperr, phi, phierr, what)

     Convert complex data given in polar coordinates `(amp,phi)` with their
     standard deviations `(amperr,phierr)` into Cartesian coordinates `(re,im)`
     and associated noise model.  The result is a hash table:

       ws.re  = real part of complex data
       ws.im  = imaginary part of complex data
       ws.wrr = statistical weight for real part of residuals
       ws.wii = statistical weight for imaginary part of residuals
       ws.wri = statistical weight for real times imaginary parts of residuals

     The quadratic penalty writes:

       err =     ws.wrr*(ws.re - re)^2
             +   ws.wii*(ws.im - im)^2
             + 2*ws.wri*(ws.re - re)*(ws.im - im);

   SEE ALSO: mira_cost.
 */
{
  local re, im, wrr, wri, wii;

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
  test = (amperr > 0.0) & (phierr > 0.0) & (amp > 0.0);
  k = where(test);
  if (goodman) {
    /* Use Goodman approximation with SIGMA such that the area of ellipsoids
     * at one standard deviation are the same:
     *    PI*SIGMA^2 = PI*(RHO*SIGMA_THETA)*SIGMA_RHO
     * hence:
     *    SIGMA = sqrt(RHO*SIGMA_THETA*SIGMA_RHO)
     * If SIGMA_THETA or SIGMA_RHO is invalid, take:
     *    SIGMA = SIGMA_RHO
     *    SIGMA = RHO*SIGMA_THETA
     * accordingly.
     */
    wrr = wri = array(double, dimsof(test));
    if (is_array(k)) {
      /* Valid amplitude and phase data. */
      wrr(k) = 1.0/(amp(k)*phierr(k)*amperr(k));
    }
    eq_nocopy, wii, wrr;
  } else {
    /* Use convex quadratic local approximation. */
    wrr = wri = wii = array(double, dimsof(test));
    if (is_array(k)) {
      /* Valid amplitude and phase data. */
      err1 = amperr(k);
      var1 = err1*err1;
      err2 = amp(k)*phierr(k);
      var2 = err2*err2;
      cs = cos_phi(k);
      sn = sin_phi(k);
      crr = cs*cs*var1 + sn*sn*var2;
      cri = cs*sn*(var1 - var2);
      cii = sn*sn*var1 + cs*cs*var2;
      a = 1.0/(var1*var2);
      wrr(k) =  a*cii;
      wri(k) = -a*cri;
      wii(k) =  a*crr;
    }
  }

  return h_new(re=re, im=im, wrr=wrr, wri=wri, wii=wii);
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
    four_wgt_err = 4.0*unref(wgt_err);
    grd(1, idx) += four_wgt_err*re;
    grd(2, idx) += four_wgt_err*im;
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
  /* RE = real part of (anti-)residuals
   * IM = imaginary part of (anti-)residuals
   * COST = sum(WRR⋅RE^2 + 2⋅WRI⋅RE⋅IM + WII⋅IM⋅IM)
   *      = sum((WRR⋅RE + WRI⋅IM)⋅RE + (WRI⋅RE + WII⋅IM)⋅IM)
   *      = (1/2)⋅sum((∂COST/∂RE)⋅RE) + (1/2)⋅sum((∂COST/∂IM)⋅IM)
   * ∂COST/∂RE = 2⋅(WRR⋅RE + WRI⋅IM)
   * ∂COST/∂IM = 2⋅(WRI⋅RE + WII⋅IM)
   */
  local idx; eq_nocopy, idx, db.idx;
  re = mira_model_re(master, idx) - db.re; /* real part of residuals */
  im = mira_model_im(master, idx) - db.im; /* imaginary part of residuals */
  wr = db.wrr*re + db.wri*im;
  wi = db.wri*re + db.wii*im;
  cost = sum(wr*re) + sum(wi*im);
  if (! is_void(grd)) {
    grd(1, idx) += (wr + wr);
    grd(2, idx) += (wi + wi);
  }
  return cost;
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
  return 2*numberof(db.idx1);
}

func _mira_t3amp_ndata(master, db)
{
  return numberof(db.idx1);
}

func _mira_t3phi_ndata(master, db)
{
  return numberof(db.idx1);
}

func _mira_t3_cost(master, db, grd)
{
  local idx1; eq_nocopy, idx1, db.idx1;
  local idx2; eq_nocopy, idx2, db.idx2;
  local idx3; eq_nocopy, idx3, db.idx3;
  trow("not yet implemented");
}

func _mira_t3amp_cost(master, db, grd)
{
  local idx1; eq_nocopy, idx1, db.idx1;
  local idx2; eq_nocopy, idx2, db.idx2;
  local idx3; eq_nocopy, idx3, db.idx3;
  trow("not yet implemented");
}

func _mira_t3phi_cost(master, db, grd)
{
  local idx1; eq_nocopy, idx1, db.idx1;
  local idx2; eq_nocopy, idx2, db.idx2;
  local idx3; eq_nocopy, idx3, db.idx3;
  trow("not yet implemented");
}

/*--------------------------------------------------------------------------*/
/* COST FUNCTIONS */

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
    grd(1, idx) += fct*re;
    grd(2, idx) += fct*im;
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
      grd(1,idx) -= fct*im;
      grd(2,idx) += fct*re;
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
      grd(1,idx) -= fct*im;
      grd(2,idx) += fct*re;
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
      grd(1,idx) -= fct*im;
      grd(2,idx) += fct*re;
    }
  } else {
    error, "unknown choice for phase cost";
  }
  return cost;
}

/*---------------------------------------------------------------------------*/
