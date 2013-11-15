/*
 * fmin.i -
 *
 *	Minimization of an univariate function for Yorick.  The method is
 *	based on original Brent's method modified to allow for different
 *	kind of bounds (both, left, right or none).
 *
 *-----------------------------------------------------------------------------
 *
 *      Copyright (C) 2001, Eric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 *	This file is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License version 2 as
 *	published by the Free Software Foundation.
 *
 *	This file is distributed in the hope that it will be useful, but
 *	WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *-----------------------------------------------------------------------------
 *
 *	$Id: fmin.i,v 1.2 2008/07/12 06:47:24 eric Exp eric $
 *	$Log: fmin.i,v $
 *	Revision 1.2  2008/07/12 06:47:24  eric
 *	 - Added final comment for setting local variables of Emacs.
 *
 *	Revision 1.1  2005/06/24 18:14:35  eric
 *	Initial revision
 *
 *-----------------------------------------------------------------------------
 */

func fmin(f, a, b, lim, tol=, all=, eps=)
/* DOCUMENT fmin(f, a, b)
       -or- fmin(f, a, b, lim)
     Get the location  of the minimum of univariate function  F(X).  F is a
     Yorick  function,  A and  B  are bounds  or  starting  values for  the
     variable X and optional LIM specifies the kind of limits for X:
       If LIM=0 or nil, there is no bounds for X: F is first evaluated at A
         and B, then the interval of  search is enlarged until a minimum is
         bracketed.
       If LIM=1,  the interval is bounded by  A: F is first  evaluated at B
         and  the interval is  enlarged (away  from A)  until a  minimum is
         bracketed  --  i.e.  the location of  the  minimum X is such that:
         A < X, if A < B; or A > X, if A > B.
       If LIM=2, the  interval is bounded by B (same as  with LIM=1 but the
         role of A and B exchanged).
       If LIM=3, the minimum is searched  in the interval (A,B) -- i.e. the
         location of the minimum X is such that: min(A,B) < X < max(A,B).

     Keyword  TOL can be  used to  specify the  relative precision  for the
     solution.  The default value is TOL=sqrt(EPS) (see below).

     If keyword ALL  is true (non-nil and non-zero)  the returned value is:
     [X, FX, XLO, XHI] where X is the approximated location of the minimum,
     FX is F(X) and XLO and XHI are the lower and upper bounds for the true
     minimum; otherwise, only X is returned.

     Keyword EPS  can be  used to specify  the machine  relative precision.
     Default value is EPS~2.22e-16,  which corresponds to IEEE standard for
     double precision.

   NOTES:
     (1) The  minimization routine  should never evaluates  F at  the given
         bounds if any.
     (2) If the function F(X) is  not unimodal, only a local minimum can be
         found.

   REFERENCES:
     The method is  based on original Brent's method  modified to allow for
     different kind of bounds (both, left, right or none).

     [1] Brent, R.P. 1973, "Algorithms for Minimization without Derivatives"
         (Englewood Cliffs, NJ: Prentice-Hall), Chapter 5.


   SEE ALSO: */
{
  /* Make sure A and B are double precision values. */
  a += 0.0;
  b += 0.0;

  /* EPS is approximately the square root of the relative machine
     precision. */
  if (is_void(eps)) eps = 2.2204460492503131e-16; /* assume IEEE double */
  tol1 = eps + 1.0;
  eps = sqrt(eps);
  tol3 = (is_void(tol) ? eps : tol)/3.0;
  /* TOL not used below */

  /* C = (3 - sqrt(5))/2 is the squared inverse of the golden ratio */
  c = 0.3819660112501051517954131656343618822796908201942371378645513772947395;

  /* S = (1 + sqrt(5))/2 = 2 - C is a constant used to increase the width
     of the interval with a golden ratio. */
  s = 1.6180339887498948482045868343656381177203091798057628621354486227052605;

  /* Original Brent's method assumes that the minimum is in (A,B) with
   * A<=B and keeps track of the following variables:
   *   X, FX = least function value found so far
   *   W, FW = previous value of X, FX
   *   V, FV = previous value of W, FW
   *   U, FU = last function evaluation
   * If the interval to consider is not bounded or only left/right bounded,
   * the idea is to find a suitable interval (A,B) where at least one
   * minimum must exists (if the function is continue) and start Brent's
   * algorithm with correct values for X, FX, ... (in order to save some
   * function evaluations).
   */
  if (! lim) {
    /* The interval of search is unlimited, we start with A, B and then
       search for a bracket. */
    x = a;
    fx = f(x);
    w = b;
    fw = f(w);

    /* Make sure X is the best location found so far, we therefore exchange
       W and X if FW <= FX (we exchange the two point in case of equality
       to alternatively search on the other side where the function is,
       numerically, flat) */
    if (fw <= fx) {
      tmp=w; w=x; x=tmp;
      tmp=fw; fw=fx; fx=tmp;
    }

    /* Loop until a bracket is found. Possible improvements: (1) use parabolic
       extrapolation to allow for bigger jumps, (2) keep track of one more
       point in order to be able to slightly reduce the size of the interval
       once a bracket has been found. */
    for (;;) {
      /* Take a golden step in the descent direction. */
      v = x + s*(x - w);
      fv = f(v);

      if (fw > fx) {
        if (fv > fx) {
          /* Bracket found: the minimum is in (V,W).  Set variables for
             Brent's method: set bounds such that A is smaller than B. */
          if (v < w) {
            a = v;
            b = w;
          } else {
            a = w;
            b = v;
          }
          break; /* branch to Brent's method */
        } else {
          /* Continue with golden search with X and V. */
          w=x; fw=fx;
          x=v; fx=fv;
        }
      } else {
        /* We are in a, numerically, flat region (FW=FX). Enlarge interval
           (V,W) by a factor 1+S ~ 2.62 */
        if (fv >= fx) {
          x=w; fx=fw;
          w=v; fw=fv;
        } else {
          w=x; fw=fx;
          x=v; fx=fv;
        }
      }
    }
  } else if (lim == 1 || lim == 2) {
    /* Interval is bounded by A or B.  Possibly exchange A and B, so that
       A is the bound and search until a bound for the other side is found. */
    if (lim == 2) {tmp=a; a=b; b=tmp;}
    w = x = b;
    fw = fx = f(x);
    for (;;) {
      v = x + s*(x - a);
      fv = f(v);
      if (fv > fx) {
        /* We have found a bound for the other side.  Set search interval
           to be (A,B) := (A,V) or (V,A) such that A is smaller than B. */
        if (v > a) {
          b = v;
        } else {
          b = a;
          a = v;
        }
        break; /* branch to Brent's method */
      }
      w=x; fw=fx;
      x=v; fx=fv;
      if (fw > fx) a = w;
    }
  } else if (lim == 3) {
    /* The minimum is to be found in (A,B) -- this is original Brent's
       method.  Make sure that A is smaller than B and set X,FX ... */
    if (a > b) {tmp=a; a=b; b=tmp;}
    v = w = x = a + c*(b - a);
    fv = fw = fx = f(x);
  } else {
    error, "bad value for keyword LIM";
  }

  /*** Brent's method. ***/

  /* Set E and D (note: the golden step instead of the parabolic step is
     taken if abs(E) is too small). */
  e = x - v;
  d = x - w;

  /* main loop starts here */
  for (;;) {
    xm = (a + b)*0.5;
    tol1 = eps*abs(x) + tol3;
    tol2 = tol1 + tol1;

    /* check stopping criterion */
    if (abs(x - xm) <= tol2 - (b - a)*0.5) {
      if (all) return [x, fx, a, b];
      return x;
    }
    if (abs(e) > tol1) {
      /* fit parabola */
      q = (x - v)*(fx - fw);
      r = (x - w)*(fx - fv);
      if (q <= r) {
	p = (x - v)*q - (x - w)*r;
	q = (r - q)*2.0;
      } else {
	p = (x - w)*r - (x - v)*q;
	q = (q - r)*2.0;
      }
      if (abs(p) < abs(0.5*q*e) && p > q*(a - x) && p < q*(b - x)) {
	/* use a parabolic-interpolation step */
        e = d;
        u = x + (d = p/q);

	/* F must not be evaluated too close to A or B */
	if (u - a < tol2 || b - u < tol2) {
	  d = (x < xm ? tol1 : -tol1);
	}
      } else {
	/* use a golden-section step */
	e = (x >= xm ? a : b) - x;
	d = c*e;
      }
    } else {
      /* use a golden-section step */
      e = (x >= xm ? a : b) - x;
      d = c*e;
    }

    /* F must not be evaluated too close to X */
    u = (abs(d) >= tol1 ? x + d : (d > 0.0 ? x + tol1 : x - tol1));
    fu = f(u);

    /* update A, B, V, W, and X */
    if (fx <= fu) {
      if (u >= x) b = u;
      else        a = u;
    }
    if (fu <= fx) {
      if (u >= x) a = x;
      else        b = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else if (fu <= fw || w == x) {
      v = w; fv = fw;
      w = u; fw = fu;
    } else if (fu <= fv || v == x || v == w) {
      v = u; fv = fu;
    }
  } /* end of main loop */
}

/*---------------------------------------------------------------------------*
 * Local Variables:                                                          *
 * mode: Yorick                                                              *
 * tab-width: 8                                                              *
 * fill-column: 75                                                           *
 * c-basic-offset: 2                                                         *
 * coding: latin-1                                                           *
 * End:                                                                      *
 *---------------------------------------------------------------------------*/
