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
            optimizer; default value is: FTOL = 1e-15 (see op_vmlmb).
     GTOL - relative gradient tolerance for the stopping criterion of the
            optimizer; default value is: GTOL = 0.0 (see op_vmlmb).
     SFTOL, SGTOL, SXTOL - control the stopping criterion of the
            line-search method in the optimizer (see op_vmlmb).

   SEE ALSO:
     mira_new, mira_config.
 */
func mira_solve(master, x, &misc, xmin=, xmax=,
                flux=, fluxerr=, zapdata=, regul=, mu=, cubic=,
                /* options for OptimPackLegacy */
                mem=, verb=, maxiter=, maxeval=, output=,
                ftol=, gtol=, sftol=, sgtol=, sxtol=,
                gpnormconv=)
{
  /* Set default values for optimizer. */
  if (is_void(ftol)) ftol =  1e-15;
  if (is_void(gtol)) gtol = 0.0;
  if (is_void(mem)) mem = 7;

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
  fg = mira_objective_function(master, flux=flux, fluxerr=fluxerr,
                               zapdata=zapdata, mu=mu, regul=regul);
  if (verb && fg.flux != 0) {
    inform, "using normalization constraint with FLUX=%g and FLUXERR=%g\n",
      fg.flux, fg.fluxerr;
  }

  printer = _mira_solve_printer;
  viewer = _mira_solve_viewer;

  // FIXME: fix bug in op_mnb:
  if (structof(x) != double || (is_void(xmin) && is_void(xmax))) {
    x += 0.0; // force copy and conversion
  }

  opl_vmlmb, _mira_solve_objfunc, x, extra=fg,
    xmin=xmin, xmax=xmax, mem=mem,
    verb=verb, viewer=viewer, printer=printer,
    maxiter=maxiter, maxeval=maxeval, output=output,
    frtol=ftol, fatol=0.0, //gatol=0.0, grtol=gtol,
    sftol=sftol, sgtol=sftol, sxtol=sftol;

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

func mira_objective_function(master, flux=, fluxerr=,
                             zapdata=, regul=, mu=)
/* DOCUMENT fg = mira_objective_function(master, ...);

     yields the objective function used by MiRA to retrieve an image.

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
  h_evaluator, fg, "_mira_evaluate_objfunc";
  return fg;
}

local _mira_solve_objfunc, _mira_evaluate_objfunc;
/* DOCUMENT fx = _mira_evaluate_objfunc(this, x, gx);
         or fx = _mira_solve_objfunc(x, gx, this);

      Private functions to evaluate the objective function and its gradient
      at X.

   SEE ALSO: mira_objective_function.
*/

func _mira_evaluate_objfunc(this, x, &grd) /* DOCUMENTED */
{
  return _mira_solve_objfunc(x, grd, this);
}

func _mira_solve_objfunc(x, &grd, this) /* DOCUMENTED */
{
  /* Normalization constraint. */
  flux = this.flux;
  fluxerr = this.fluxerr;
  strict_flux = (flux > 0 && fluxerr == 0);
  loose_flux  = (flux > 0 && fluxerr > 0);

  if (strict_flux && (xsum = sum(x)) > 0) {
    xscl = flux/double(xsum);
    if (xscl != 1) {
      x *= xscl;
    }
  }

  /* Data penalty. */
  if (this.ndata > 0) {
    fdata = mira_cost_and_gradient(this.master, x, grd);
  } else {
    fdata = 0.0;
  }

  /* Prior penalty. */
  if (this.mu > 0) {
    rgl_update, this.regul, x;
    fprior = rgl_get_penalty(this.regul, x);
    if (this.ndata > 0) {
      grd += rgl_get_gradient(this.regul, x);
    } else {
      eq_nocopy, grd, rgl_get_gradient(this.regul, x);
    }
  } else {
    fprior = 0.0;
    if (this.ndata <= 0) {
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
  h_set, this, nevals = this.nevals + 1;

  /* Compute total penalty and update best solution so far. */
  ftot = fdata + fprior;
  if (! h_has(this, best_f=) || this.best_f > ftot) {
    best_x = x; // make a copy
    best_g = grd; // make a copy
    h_set, this, best_f = ftot, best_fdata = fdata,
      best_fprior = fprior, best_g = best_g, best_x = best_x,
      best_nevals = this.nevals,
      best_niters = this.niters + 1;
  }
  return ftot;
}

/* The following function is just used to store the number of objective
   function evaluations and the number of algorithm iterations. */
func _mira_solve_viewer(x, this, ws)
{
  h_set, this, niters = this.niters + 1;
}

func _mira_solve_printer(output, iter, eval, cpu, fx, gnorm, steplen, x, extra)
{
  if (eval == 1) {
    write, output, format="# %s\n# %s\n",
      "ITER  EVAL   CPU (ms)        FUNC               <FDATA>     FPRIOR    GNORM   STEPLEN",
      "-------------------------------------------------------------------------------------";
  }
  fdata = (extra.ndata > 0 ? extra.best_fdata/extra.ndata : 0.0);
  fprior = (extra.mu > 0 ? extra.best_fprior/extra.mu : 0.0);
  write, output,
    format=" %5d %5d %10.3f  %+-24.15e%-11.3e%-11.3e%-9.1e%-9.1e\n",
    iter, eval, cpu, fx, fdata, fprior, gnorm, step;
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
