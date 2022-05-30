/*
 * mira2_solver.i -
 *
 * Image reconstruction solver of MiRA.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of MiRA, a "Multi-aperture Image Reconstruction
 * Algorithm", <https://github.com/emmt/MiRA>.
 *
 * Copyright (C) 2001-2021, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

local mira_solve;
/* DOCUMENT img = mira_solve(data, key1=value1, key2=value2, ...);
         or img = mira_solve(data, img_init, key1=value1, key2=value2, ...);
         or img = mira_solve(data, img_init, misc, key1=value1, ...);

     Builds an image from the interferometric data stored into instance DATA
     (see mira_new) which must have been properly initialized for
     reconstruction (see mira_config).  The reconstruction is done by
     minimizing an objective function which is the sum of the likelihood
     penalty (which enforces agreement of the model image with the data) and a
     regularization penalty (which enforces some priors about the image).

     Optional argument IMG_INIT is the initial image.  Most arguments are
     provided by keywords (see below). At least, you want to specify the
     regularization method and level by the keywords REGUL and MU.  Usually,
     you also want to impose a positivity constraint with XMIN=0.0 (use a small
     but strictly positive value if regularization such as entropy with a
     logarithm is used) and a normalization constraint with FLUX=1.0
     (which is suitable for OI-FITS data, otherwise the actual value should be
     equal to the total flux in the image).

     Optional argument MISC is a simple symbol name to store the following
     informations:

         misc.ndata  = number of measurements
         misc.fdata  = data penalty per measurement
         misc.mu     = regularization weight
         misc.fprior = prior penalty
         misc.ftot   = total penalty (ndata*fdata + mu*fprior)
         misc.gpnorm = Euclidean norm of projected gradient
         misc.nevals = number of evaluations of the objective function
         misc.niters = number of iterations of the algorithm

   KEYWORDS:
     XMIN - minimum allowed value in the image; can be a scalar or a pixel-wise
            value; there is no limit if this option is not set.
     XMAX - maximum allowed value in the image; can be a scalar or a pixel-wise
            value; there is no limit if this option is not set.
     FLUX - value of the sum of the image pixels to impose a constraint on the
            total flux.
     FLUXERR - standard deviation of the flux.  Can be 0 to apply a srict
            equality constraint for the flux, or "auto" to set it automatically
            according to the value of FLUX and to the effective number of
            measurements, must be strictly positive otherwise.  If unspecified,
            FLUXERR=0 is assumed, i.e. a strict constraint is applied.
     ZAPDATA - is true to completely ignore the data.
     MU    - regularization level; the higher is MU the more the solution is
             influenced by the priors set by the regularization.
     REGUL - regularization method (see rgl_new).
     CUBIC - is true to use cubic interpolation if initial image needs to be
            resized.
     MAXITER - maximum number of iterations, unlimited if not set.
     MAXEVAL - maximum number of evaluations of the objective function,
            unlimited if not set.
     OUTPUT - output stream/file for iteration informations; default is
            standard text output (see keyword VERB).
     VERB - verbose level: informations get printed out every VERB iterations
            and at convergence.
     MEM  - control the memory usage by the optimizer; the value is the number
            of corrections and gradient differences memorized by the variable
            metric algorithm; by default, MEM=7 (see op_vmlmb).
     FTOL - relative function tolerance for the stopping criterion of the
            optimizer (opl_vmlmb/optm_vmlmb); default value is: FTOL = 1e-15.
     GTOL - relative gradient tolerance for the stopping criterion of the
            optimizer (opl_vmlmb/optm_vmlmb); default value is: GTOL = 0.0.
     XTOL - relative tolerance in the variables for the stopping criterion of
            the optimizer (optm_vmlmb); default value is: XTOL = 1e-5.
     SFTOL, SGTOL, SXTOL - control the stopping criterion of the
            line-search method in the optimizer (opl_vmlmb/optm_vmlmb).

   SEE ALSO:
     mira_new, mira_config.
 */
func mira_solve(master, x, &misc, xmin=, xmax=,
                flux=, fluxerr=, zapdata=, regul=, mu=, cubic=,
                /* options for VMLMB */
                useoptimpacklegacy=, blmvm=,
                mem=, verb=, maxiter=, maxeval=, output=,
                ftol=, gtol=, xtol=, sftol=, sgtol=, sxtol=,
                gpnormconv=)
{
  /* Set default values for optimizer. */
  if (is_void(ftol)) ftol =  1e-15;
  if (is_void(gtol)) gtol = 0.0;
  if (is_void(xtol)) xtol = 1e-5;
  if (is_void(mem)) mem = 7;
  //if (is_void(useoptimpacklegacy)) useoptimpacklegacy = 1n;

  /* Update internal cache (if needed). */
  mira_update, master;

  /* Default solution. */
  dims = mira_image_size(master);
  if (is_void(x)) {
    x = random(dims);
  } else if (! mira_same_dimensions(dimsof(x), dims)) {
    throw, "incompatible image dimensions";
    x = mira_rescale(x, dims, cubic=cubic); // FIXME:
  }

  /* Build functor to evaluate the objective function. */
  fg = mira_objfunc(master, flux=flux, fluxerr=fluxerr,
                    zapdata=zapdata, mu=mu, regul=regul);
  if (verb && fg.flux != 0) {
    inform, "using normalization constraint with FLUX=%g and FLUXERR=%g\n",
      fg.flux, fg.fluxerr;
  }

  // FIXME: fix bug in op_mnb:
  if (structof(x) != double || (is_void(xmin) && is_void(xmax))) {
    x += 0.0; // force copy and conversion
  }

  if (is_func(optm_vmlmb) && !useoptimpacklegacy) {
      inform, "using `%s` (maxeval=%d, maxiter=%d)\n", "optm_vmlmb",
          (is_void(maxeval) ? -1 : maxeval),
          (is_void(maxiter) ? -1 : maxiter);
      local fx, gx, status;
      optm_vmlmb, fg, x, fx, gx, status, lower=xmin, upper=xmax, mem=mem,
          fmin=0, // lnsrch=, delta=, epsilon=, lambda=,
          ftol=[0.0,ftol], gtol=gtol, xtol=xtol, blmvm=blmvm,
          maxiter=maxiter, maxeval=maxeval, verb=verb, cputime=1n,
          observer=_mira_optm_observer,
          printer=_mira_optm_printer, output=output,
          throwerrors=1n;
  } else {
      inform, "using `%s` (maxeval=%d, maxiter=%d)\n", "opl_vmlmb",
          (is_void(maxeval) ? -1 : maxeval),
          (is_void(maxiter) ? -1 : maxiter);
      opl_vmlmb, fg, x,
          xmin=xmin, xmax=xmax, mem=mem,
          verb=verb, viewer=_mira_opl_viewer, printer=_mira_opl_printer,
          maxiter=maxiter, maxeval=maxeval, output=output,
          frtol=ftol, fatol=0.0, //gatol=0.0, grtol=gtol,
          sftol=sftol, sgtol=sftol, sxtol=sftol;
  }

  /* Retrieve the best solution so far (which is already normalized) and the
     various termes of the objective function and its gradient at the
     solution. */
  eq_nocopy, x, fg.best_x;
  gpnorm =  mira_projected_gradient_norm(x, fg.best_g, xmin=xmin, xmax=xmax);
  misc = h_new(ndata  = fg.ndata,
               fdata  = (fg.ndata > 0 ? fg.best_fdata/fg.ndata : 0.0),
               mu     = max(fg.mu, 0.0),
               fprior = (fg.mu > 0.0 ? fg.best_fprior/fg.mu : 0.0),
               ftot   = fg.best_f,
               gpnorm = gpnorm,
               nevals = fg.best_nevals,
               niters = fg.best_niters);
  return x;
}

func mira_objfunc(master, flux=, fluxerr=, zapdata=, regul=, mu=)
/* DOCUMENT fg = mira_objfunc(master, ...);

     Get the objective function used by MiRA to retrieve an image. The
     returned object can be called as a function:

         local gx;
         fx = fg(x, gx);

     where `x` is the image, `gx` is a variable to store the gradient of the
     objective function at `x` and `fx` is the value objective function at `x`.

   SEE ALSO: mira_solve.
 */
{
  /* Check whether to consider the data or not. */
  if (is_array(zapdata) && ! is_scalar(zapdata)) {
    throw, "invalid value for `zapdata`";
  }

  /* Effective number of measurements. */
  ndata = mira_ndata(master);

  /* Flux constraint and its error bar.  By default, a strict constraint is
     applied and a normalized image is assumed. */
  if (! scalar_double(flux, 1) || flux < 0) {
    throw, "invalid value for `flux`";
  }
  if (is_scalar(fluxerr) && is_string(fluxerr) && fluxerr == "auto") {
    fluxerr = flux/sqrt(ndata);
  } else if (! scalar_double(fluxerr, 0) || fluxerr < 0) {
    throw, "invalid value for `fluxerr`";
  }

  /* Setup for regularization. */
  if (is_void(regul)) {
    mu = 0.0;
  } else {
    temp = rgl_get_global_weight(regul);
    if (is_void(mu)) {
      mu = temp;
    } else if (is_scalar(mu) && (is_real(mu) || is_integer(mu)) && mu >= 0.0) {
      mu = double(mu);
      if (mu != temp) {
        rgl_set_global_weight, regul, mu;
      }
    } else {
      error, "global regularization weight MU must be a non negative scalar";
    }
  }

  /* Build a functor. */
  fg = h_new(flux = flux, fluxerr = fluxerr, master = master,
             ndata = (zapdata ? 0 : ndata),
             mu = mu, regul = regul,
             niters = 0, nevals = 0);
  h_evaluator, fg, "_mira_eval_objfunc";
  return fg;
}

func _mira_eval_objfunc(fg, x, &grd)
/* DOCUMENT fx = _mira_eval_objfunc(fg, x, gx);

      Private functions to evaluate the objective function F and its gradient
      at X.  FG is the object implementing the objective function.

   SEE ALSO: mira_objfunc.
*/
{
  /* Normalization constraint. */
  flux = fg.flux;
  fluxerr = fg.fluxerr;
  strict_flux = (flux > 0 && fluxerr == 0);
  loose_flux  = (flux > 0 && fluxerr > 0);

  if (strict_flux && (xsum = sum(x)) > 0) {
    xscl = flux/double(xsum);
    if (xscl != 1) {
      x *= xscl;
    }
  }

  /* Data penalty. */
  if (fg.ndata > 0) {
    fdata = mira_cost_and_gradient(fg.master, x, grd);
  } else {
    fdata = 0.0;
  }

  /* Prior penalty. */
  if (fg.mu > 0) {
    rgl_update, fg.regul, x;
    fprior = rgl_get_penalty(fg.regul, x);
    if (fg.ndata > 0) {
      grd += rgl_get_gradient(fg.regul, x);
    } else {
      eq_nocopy, grd, rgl_get_gradient(fg.regul, x);
    }
  } else {
    fprior = 0.0;
    if (fg.ndata <= 0) {
      grd = array(double, dimsof(x));
    }
  }

  /* Account for the flux constraint. */
  if (strict_flux && (xsum !=  0)) {
    grd = xscl*grd - sum(grd*x)/xsum;
  } else if (loose_flux) {
    fluxres = sum(x) - flux;
    if (fluxres != 0) {
      fluxwgt = 1.0/fluxerr^2;
      fdata += fluxwgt*fluxres*fluxres;
      grd += 2.0*fluxwgt*fluxres;
    }
  }

  /* Number of evaluations of the objective function. */
  h_set, fg, nevals = fg.nevals + 1;

  /* Compute total penalty and update best solution so far. */
  ftot = fdata + fprior;
  if (! h_has(fg, best_f=) || fg.best_f > ftot) {
    best_x = x; // make a copy
    best_g = grd; // make a copy
    h_set, fg, best_f = ftot, best_fdata = fdata,
      best_fprior = fprior, best_g = best_g, best_x = best_x,
      best_nevals = fg.nevals,
      best_niters = fg.niters + 1;
  }
  return ftot;
}

// Use the "viewer" of `opl_vmlmb` to track the number of iterations.
func _mira_opl_viewer(x, fg, ws)
{
  h_set, fg, niters = fg.niters + 1;
}

// Use the "observer" of `optm_vmlmb` to track the number of iterations.
func _mira_optm_observer(iters, evals, rejects, t, x, f, g, gnorm, alpha, fg)
{
    h_set, fg, niters = iters;
}

func _mira_optm_printer(output, iters, evals, rejects, t, x, f, g, gnorm,
                        alpha, fg)
{
  if (evals == 1) {
    write, output, format="# %s\n# %s\n",
      "ITER  EVAL     CPU (ms)        FUNC               <FDATA>     FPRIOR    GNORM   STEPLEN",
      "---------------------------------------------------------------------------------------";
  }
  fdata = (fg.ndata > 0 ? fg.best_fdata/fg.ndata : 0.0);
  fprior = (fg.mu > 0 ? fg.best_fprior/fg.mu : 0.0);
  write, output,
    format=" %5d %5d %12.3f  %+-24.15e%-11.3e%-11.3e%-9.1e%-9.1e\n",
    iters, evals, t*1e3, f, fdata, fprior, gnorm, alpha;
}

func _mira_opl_printer(output, iter, eval, cpu, f, gnorm, steplen, x, fg)
{
  if (eval == 1) {
    write, output, format="# %s\n# %s\n",
      "ITER  EVAL     CPU (ms)        FUNC               <FDATA>     FPRIOR    GNORM   STEPLEN",
      "---------------------------------------------------------------------------------------";
  }
  fdata = (fg.ndata > 0 ? fg.best_fdata/fg.ndata : 0.0);
  fprior = (fg.mu > 0 ? fg.best_fprior/fg.mu : 0.0);
  write, output,
    format=" %5d %5d %12.3f  %+-24.15e%-11.3e%-11.3e%-9.1e%-9.1e\n",
    iter, eval, cpu*1e3, f, fdata, fprior, gnorm, step;
}

func mira_projected_gradient_norm(x, gx, xmin=, xmax=)
{
  // FIXME: strict flux constraint not taken into account
  local gp;
  if (is_void(xmin)) {
    eq_nocopy, gp, gx;
  } else {
    gp = gx*((gx < 0.0)|(x > xmin));
  }
  if (! is_void(xmax)) {
    gp *= ((gx > 0.0)|(x < xmax));
  }
  return sqrt(sum(gp*gp));
}

/*---------------------------------------------------------------------------*/
/* REGULARIZATIONS */

func mira_new_compactness_prior(dat, fwhm)
/* DOCUMENT rgl = mira_new_compactness_prior(dat, fwhm);

     yields a regularization which implements a compactness prior.  The default
     shape is isotropic with a Cauchy profile of full-width at half-maximum
     FWHM (in radians).

   SEE ALSO rgl_quadratic.
 */
{
  s = 2.0/fwhm;
  x = s*mira_image_x(dat);
  y = s*mira_image_y(dat);
  w = 1.0 + (x*x + (y*y)(-,));
  rgl = rgl_new("quadratic");
  rgl_config, rgl, "W", linop_new("diagonal", w);
  return rgl;
}
