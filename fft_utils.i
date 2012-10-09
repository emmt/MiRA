/*
 * fft_utils.i --
 *
 *	Useful routines for FFT operations in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 *      Copyright (C) 1995, Eric Thiébaut <thiebaut@obs.univ-lyon1.fr>
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
 * History:
 *	$Id: fft_utils.i,v 1.15 2008/12/24 11:04:04 eric Exp $
 *	$Log: fft_utils.i,v $
 *	Revision 1.15  2008/12/24 11:04:04  eric
 *	 - Fixed a bug in fft_best_dim which might not return the smallest
 *	   multiple of 2, 3, and/or 5 greater of equal its argument (thanks
 *	   to Clementine Bechet for the fix).
 *
 *	Revision 1.14  2008/09/24 16:18:30  eric
 *	 - Fixed bug in fft_recenter found by FerrÃ©ol Soulez.
 *
 *	Revision 1.13  2008/07/12 06:38:26  eric
 *	 - Considerable speed-up of fft_gaussian_mtf and fft_gaussian_psf for
 *	   multidimensional array (accounting for the fact that the Gaussian
 *	   is separable).  Also FWHM can now be a scalar or a vector with as
 *	   many values as number of dimensions.
 *	 - New function fft_get_ndims.
 *
 *	Revision 1.12  2007/04/24 07:11:43  eric
 *	 - Function grow_dimlist replaced by make_dimlist.
 *	 - Function fft_paste fixed.
 *
 *	Revision 1.11  2005/10/18 18:20:39  eric
 *	 - Oops, the previous bug fix yields a new bug, OK now fft_dist
 *	   should work again for n-D arrays (with n>1).
 *
 *	Revision 1.10  2005/10/18 18:09:19  eric
 *	 - new function: fft_paste;
 *	 - fixed bug in fft_dist for 1-D dimension lists which affected
 *	   functions such as fft_smooth for 1-D arrays (thanks to Christophe
 *	   Pichon for pointing the bug);
 *
 *	Revision 1.9  2004/10/11 11:18:21  eric
 *	 -  Fix fft_dist() so that it is as flexible as, e.g., array() for
 *	    the dimension list.
 *
 *	Revision 1.8  2004/08/31 16:20:22  eric
 *	 - New routines for Fourier interpolation: fft_fine_shift,
 *	   fft_unphasor, fft_interp, fft_interp_complex, fft_interp_real.
 *
 *	Revision 1.7  2003/08/23 09:55:57  eric
 *	 - new functions: fft_gaussian_mtf and fft_gaussian_psf;
 *	 - *** POSSIBLE INCOMPATIBILITY *** change order of args in
 *	   fft_recenter and use fft_setup to speed up FFT's;
 *
 *	Revision 1.6  2003/01/31 15:58:32  eric
 *	 - Added new routines: fft_recenter, fft_smooth and reverse_dims.
 *
 *	Revision 1.5  2002/11/20 09:20:56  eric
 *	 - new keywords SQUARE and NYQUIST in fft_dist
 *	 - make use of grow_dimlist in "utils.i"
 *
 *	Revision 1.4  2002/11/14 10:59:38  eric
 *	 - new graphics routines: fft_plh, fft_plc, fft_plfc
 *	 - new routines: abs2, fft_roll_1d, fft_roll_2d, __fft & Co
 *	 - introduction documentation by: help, fft_utils;
 *
 *	Revision 1.3  2001/04/23 14:16:15  eric
 *	New routine: fft_symmetric_index.
 *
 *	Revision 1.2  2001/04/06 16:42:59  eric
 *	 - new routine: fft_plg
 *	 - fixed routine: fft_recenter_at_max
 *
 *	Revision 1.1  2001/03/23 16:45:54  eric
 *	Initial revision
 *
 *-----------------------------------------------------------------------------
 */

require, "utils.i";

local fft_utils;
/* DOCUMENT: FFT utility routines in "fft_utils.i"
     This package  is mainly written  to deal with the  particular indexing
     rules in FFT transformed arrays.  The following routines are provided:
       abs2         - squared absolute value.
       fft_best_dim - get best dimension to compute FFT.
       fft_centroid - get centroid in FFT arrays.
       fft_convolve - compute discrete convolution thanks to FFT.
       fft_dist     - compute length of FFT frequencies/coordinates.
       fft_fine_shift, fft_unphasor - shift/roll an array by non-integer
                      offset by means of Fourier interpolation.
       fft_gaussian_mtf - compute Gaussian modulation transfer function.
       fft_gaussian_psf - compute Gaussian point spread function.
       fft_indgen   - generate index of FFT frequencies/coordinates.
       fft_interp, fft_interp_complex, fft_interp_real - interpolate array
                      at non-integer offsets.
       fft_plc      - plot contours of 2D FFT.
       fft_plfc     - plot filled contours of 2D FFT.
       fft_plg      - plot 1D FFT as curve.
       fft_plh      - plot 1D FFT with stairs.
       fft_pli      - plot 2D FFT as image.
       fft_recenter - recenter array with respect to a template.
       fft_recenter_at_max - recenter FFT arrays at their maximum.
       fft_roll_1d  - roll dimension of 1D arrays.
       fft_roll_2d  - roll dimension of 2D arrays.
       fft_shift_phasor - get complex phasor for arbitrary shift.
       fft_smooth   - smooth an array by convolution with a gaussian.
       fft_symmetric_index - get hermitian-symmetry index for FFT arrays.
       reverse_dims - reverse all dimensions of an array.
       __fft        - expert driver for repeated FFT's with same dimensions.
       __fft_init   - initialization for __fft.

   SEE ALSO: fft, fftw. */

func abs2(x)
/* DOCUMENT abs2(x)
     Returns abs(X)^2
   SEE ALSO: abs. */
{
  if (structof(x) != complex) return x*x;
  y = x.im;
  x = double(x);
  return x*x + y*y;
}

func fft_best_dim(len)
/* DOCUMENT fft_best_dim(len);
     Return the smallest integer which is greater or equal LEN and which is
     a multiple of powers of 2, 3 and/or 5
   SEE ALSO fft_indgen, fft. */
{
  len = long(len);
  best = 2*len;
  for (i5 = 1; i5 < best; i5 *= 5) {
    for (i3 = i5; i3 < best; i3 *= 3) {
      /* last loop (power of 2) is exited as soon as N >= LEN */
      for (i2 = i3; i2 < len; i2 *= 2)
        ; /* empty loop body */
      if (i2 == len) return len;
      best = min(best, i2);
    }
  }
  return best;
}

func fft_indgen(dim) { return (u= indgen(0:dim-1)) - dim*(u > dim/2); }
/* DOCUMENT fft_indgen(len)
     Return FFT frequencies along a dimension of length LEN.
   SEE ALSO: indgen, span, fft_dist, fft_freqlist, fft_symmetric_index. */

func fft_dist(.., nyquist=, square=)
/* DOCUMENT fft_dist(dimlist);
       -or- fft_dist(dim1, dim2, ...);
     Returns Euclidian lenght of spatial frequencies in frequel units for a
     FFT of dimensions DIMLIST.

     If keyword  NYQUIST is true,  the frequel coordinates get  rescaled so
     that the Nyquist frequency is  equal to NYQUIST along every dimension.
     This is obtained by using coordinates:
        (2.0*NYQUIST/DIM(i))*fft_indgen(DIM(i))
     along i-th dimension of lenght DIM(i).

     If  keyword  SQUARE is  true,  the square  of  the  Euclidian norm  is
     returned instead.

   SEE ALSO: fft_indgen, fft_symmetric_index. */
{
  /* Build dimension list. */
  local arg, dims;
  while (more_args()) {
    eq_nocopy, arg, next_arg();
    if ((s = structof(arg)) == long || s == int || s == short || s == char) {
      /* got an integer array */
      if (! (n = dimsof(arg)(1))) {
        /* got a scalar */
        grow, dims, arg;
      } else if (n == 1 && (n = numberof(arg) - 1) == arg(1)) {
        /* got a vector which is a valid dimension list */
        if (n) grow, dims, arg(2:);
      } else {
        error, "bad dimension list";
      }
    } else if (! is_void(arg)) {
      error, "unexpected data type in dimension list";
    }
  }
  if (! (n = numberof(dims))) return 0.0; /* scalar array */
  if (min(dims) <= 0) error, "negative value in dimension list";

  /* Build square radius array one dimension at a time, starting with the
     last dimension. */
  if (is_void(nyquist)) {
    r2 = (u = double(fft_indgen(dims(n))))*u;
    while (--n >= 1) {
      r2 = r2(-,..) + (u = double(fft_indgen(dims(n))))*u;
    }
  } else {
    s = 2.0*nyquist;
    dim = dims(n);
    r2 = (u = (s/dim)*fft_indgen(dim))*u;
    while (--n >= 1) {
      dim = dims(n);
      r2 = r2(-,..) + (u = (s/dim)*fft_indgen(dim))*u;
    }
  }
  return (square ? r2 : sqrt(r2));
}

func fft_freqlist(dimlist)
/* DOCUMENT ptr = fft_freqlist(dimlist)
     returns a vector of DIMLIST(1) pointers with normalized FFT
     frequencies along all dimensions of DIMLIST (must be
     dimsof(SOME_ARRAY)) with "adequate" geometry:

       *ptr(1) =   (2*pi/dimlist(2))*fft_indgen(dimlist(2));
       *ptr(2) =  [(2*pi/dimlist(3))*fft_indgen(dimlist(3))];
       *ptr(3) = [[(2*pi/dimlist(4))*fft_indgen(dimlist(4))]];
        ...

   SEE ALSO: fft_indgen, fft_shift_phasor. */
{
  /* Precompute (scaled) Fourier frequencies along every dimensions (the
     trick is to build these as "vectors" with adequate dimension list so
     as to minimize the number of operations during the search). */
  PI = 3.1415926535897932384626433832795029;
  ndims = dimlist(1);
  ptr = array(pointer, ndims);
  for (k=1 ; k<=ndims ; ++k) {
    len = dimlist(k + 1);
    ws = array(1, k + 1); /* to build dimension list of k-th dimension */
    ws(1) = k;
    ws(0) = len;
    (ws = array(double, ws))(*) = (2.0*PI/len)*fft_indgen(len);
    ptr(k) = &ws;
  }
  return ptr;
}

func fft_smooth(a, fwhm, setup=)
/* DOCUMENT fft_smooth(a, fwhm)
       -or- fft_smooth(a, fwhm, setup=workspace)
     Returns array A smoothed along  all its dimensions by convolution by a
     gaussian  with  full  width  at  half maximum  equals  to  FWHM.   See
     fft_setup for the meaning of keyword SETUP.

   SEE ALSO: fft, fft_setup, fft_gaussian_mtf. */
{
  dims = dimsof(a);
  if (is_void(setup)) setup=fft_setup(dims);
  as_double = (structof(a) != complex);
  a = fft((1.0/numberof(a))*fft_gaussian_mtf(dims, fwhm)*
          fft(a, +1, setup=setup), -1, setup=setup);
  return (as_double ? double(a) : a);
}

local fft_gaussian_psf;
local fft_gaussian_mtf;
/* DOCUMENT fft_gaussian_psf(dimlist, fwhm)
       -or- fft_gaussian_mtf(dimlist, fwhm)

     Returns   normalized   Gaussian  point   spread   function  (PSF)   or
     corresponding modulation  transfer function (MTF)  with dimension list
     DIMLIST  and full  width at  half maximum  equals to  FWHM  along each
     dimensions (in the  PSF space).  Up to errors  due to limited support,
     numerical precision and finite sampling, the PSF and the MTF obey:

        sum(PSF) = MTF(1) = 1                       (normalization)
        MTF = fft(PSF, +1)
        PSF = fft(MTF, -1)/numberof(MTF)

     where MTF(1) is the 0-th frequency in the MTF.  The standard deviation
     SIGMA and the FWHM are related by:

        FWHM = sqrt(8*log(2))*SIGMA
             ~ 2.354820045031*SIGMA

     Note that, owing  to the limited size of  the support and/or numerical
     precision, these properties may not be perfectly met; for that reason,
     _always_ compute directly  what you need, e.g. do not  take the FFT of
     the PSF if what  you need is the MTF.  Also note  that the geometry is
     that of  the FFT and that  for unequal dimension lengths,  the PSF has
     the same width (in "pixels") along every dimension but not the MTF.

     FWHM can  be a scalar  or a  vector with as  many values as  number of
     dimensions.

   SEE ALSO: fft_get_ndims, fft_dist, fft_smooth. */

func fft_gaussian_psf(dims, fwhm)
{
  ndims = fft_get_ndims(dims);
  if (ndims <= 0L) {
    if (ndims == 0L) return (1.6651092223153955127063292897904020952612/fwhm);
    error, "bad dimension list";
  }

  //r = (sqrt(log(16.0)/pi)/fwhm) + array(0.0, ndims);
  //s = (sqrt(log(16.0))/fwhm) + array(0.0, ndims);
  r = (0.9394372786996513337723403284101825868414/fwhm) + array(0.0, ndims);
  s = (1.6651092223153955127063292897904020952612/fwhm) + array(0.0, ndims);
  j = ndims;
  u = s(j)*fft_indgen(dims(j));
  p = exp(-u*u);
  q = r(j);
  while (--j >= 1L) {
    u = s(j)*fft_indgen(dims(j));
    p = exp(-u*u)*p(-,..);
    q *= r(j);
  }
  return q*p;
}

func fft_gaussian_mtf(dims, fwhm)
{
  ndims = fft_get_ndims(dims);
  if (ndims <= 0L) {
    if (ndims == 0L) return 1.0;
    error, "bad dimension list";
  }

  // s = (pi/sqrt(log(16)))*fwhm/dim
  s = 1.8867186677527935983734215966417072034386*fwhm/dims;
  j = ndims;
  u = s(j)*fft_indgen(dims(j));
  p = exp(-u*u);
  while (--j >= 1L) {
    u = s(j)*fft_indgen(dims(j));
    p = exp(-u*u)*p(-,..);
  }
  return p;
}

func fft_get_ndims(&dims)
/* DOCUMENT fft_get_ndims(dimlist)
     Returns the number of dimensions in dimension list DIMLIST and modify
     (in-place) DIMLIST to be only the list of dimensions (that is without
     the number of dimensions).  If DIMLIST is invalid, -1 is returned.
     
   SEE ALSO: dimsof.
 */
{
  if (((s = structof(dims)) == long || s == int
       || s == short || s == char) && min(dims) > 0) {
    temp = dimsof(dims)(1);
    if (temp == 0L) {
      dims = long(dims);
      return 1L;
    } else if (temp == 1L && (ndims = dims(1)) == numberof(dims) - 1L) {
      dims = (ndims >= 1L ? long(dims(2:0)) : []);
      return ndims;
    }
    return ndims;
  } else if (is_void(dims)) {
    return 0L;
  }
  return -1L;
}

/*
 * Notes:
 *   The FFT of a Gaussian of given FWHM is:
 *      exp(-pi^2*fwhm^2*(k/dim)^2/log(16))
 *   where K is the FFT index; hence:
 *      Nyquist = sqrt(pi^2*fwhm^2*(dim/2/dim)^2/4/log(2))
 *              = pi*fwhm/4/sqrt(log(2))
 *              ~ 0.9433593338763967992*fwhm
 *   The Gaussian is (N is the number of dimensions):
 *      ((sqrt(log(16)/pi)/fwhm)^n)*exp(-log(16)*(k/fwhm)^2)
 *     ~ ((0.9394372786996513337/fwhm)^n)*exp(...)
 *
 *   Constants (with 40 significant digits):
 *      pi/4/sqrt(log(2)) = 0.9433593338763967991867107983208536017193;
 *      sqrt(log(16)/pi)  = 0.9394372786996513337723403284101825868414;
 *      log(16)           = 2.772588722239781237668928485832706272302;
 */

/*---------------------------------------------------------------------------*/
/* SIMPLE FAST FOURIER TRANSFORM */

func __fft_init(dimlist)
/* DOCUMENT __fft_init, dimlist;
     Initializes  FFT  workspace  for  further  calls to  __fft  (to  see).
     DIMLIST is the dimension list of the arrays to transform.  The routine
     defines 2 external symbols:
       __fft_setup  - used to store the FFT workspace)
       __fft_number - used to keep track of the number of calls to __fft
     In order to avoid namespace pollution/clash, a routine that uses __fft
     should declare these symbols as local before calling __fft_init, e.g.:
       local __fft_setup, __fft_number;
       __fft_init, dimlist;

   SEE ALSO: fft_setup, __fft. */
{
  extern __fft_setup, __fft_number;
  dims = dimlist(2:);
  ndims = numberof(dims);
  __fft_setup = array(pointer, ndims);
  __fft_number = 0;
  for (i=1 ; i<=ndims ; ++i) {
    if (! __fft_setup(i)) {
      dim = dims(i);
      ws = array(double, 6*dim + 15);
      fft_init, dim, ws;
      __fft_setup(where(dims == dim)) = &ws;
    }
  }
}

func __fft(x, dir)
/* DOCUMENT __fft(x);
       -or- __fft(x, dir);
     Replacement for Yorick's fft to speed up fast Fourier transforms (FFT)
     in,  e.g.,  iteratives  algorithms  that  necessitate  computation  of
     several  FFT's  of  arrays  with  same dimension  list.   The  FFT  is
     performed on  all dimensions of  X and DIR  must be a  scalar (default
     +1). FFT workspace must be initialized by __fft_init.

   SEE ALSO: __fft_init, fft. */
{
  extern __fft_setup, __fft_number;

  /* Make a private copy of input array even if it is already complex and
     get its dimension list. */
  if (is_void(dir)) dir = +1;
  x = complex(x);
  dims = dimsof(x);
  ndims = dims(1);
  dims = dims(2:0);
  if (structof(__fft_setup) != pointer ||
      numberof(__fft_setup) != ndims) error, "unitialized FFT workspace";

  /* Do the transform along every dimension of X. */
  len = 6*dims + 15; // expected length of workspace vectors
  std = 1;
  top = numberof(x);
  for (i=1 ; i<=ndims ; i++) {
    dim = dims(i);
    ws = __fft_setup(i);
    if (numberof(*ws) != len(i) || structof(*ws) != double)
      error, swrite(format="bad FFT workspace for dimension length %d", dim);
    top /= dim;
    fft_raw, dir, x, std, dim, top, ws;
    std *= dim;
  }

  /* increment counter and return result */
  if (is_void(__fft_number)) __fft_number= 0;
  ++__fft_number;
  return x;
}

/*---------------------------------------------------------------------------*/

func fft_symmetric_index(..)
/* DOCUMENT fft_symmetric_index(dimlist)
       -or- fft_symmetric_index(dim1, dim2, ...);
     Returns  indices  of  hermitian-symmetry  transform  for  a  FFT  with
     dimension list DIMLIST.  For instance,  if A is a N-dimensional array,
     then:
        AP= A(fft_symmetric_index(dimsof(A)))
     is  equal to array  A with  its coordinates  negated according  to FFT
     convention:
        AP(X1, X2, ..., XN) = A(-X1, -X2, ..., -XN)
     consequently if A is hermitian then:
        AP= conj(A).

   SEE ALSO: fft_indgen. */
{
  /* Build dimension list. */
  dimlist = [0];
  while (more_args()) make_dimlist, dimlist, next_arg();

  /* Compute result starting by last dimension. */
  local u;
  if ((n = numberof(dimlist)) == 1) return 1; /* scalar array */
  for (k=n ; k>1 ; --k) {
    dim = dimlist(k);
    (q = indgen(dim:1:-1))(1) = 0; // neg. of freq along that dim
    u= k<n ? dim*u(-,..) + q : q;
  }
  return u + 1; /* <== indices start at 1 in Yorick */
}

/*---------------------------------------------------------------------------*/

func _fft_centroid(a1, repeat)
/* DOCUMENT _fft_centroid(a1)
       -or- _fft_centroid(a1, repeat)
     Working routine for fft_centroid: return position of centroid of 1D array
     A1 assuming FFT geometry for the coordinates.

   SEE ALSO fft_centroid. */
{
  dim= numberof(a1);
  u= fft_indgen(dim);
  x0= u(a1(mxx));
  x= double(u);
  do {
    w= x0==0 ? a1 : roll(a1, -x0);
    if (dim%2 == 0) w(dim/2+1)= 0.0; // remove Nyquist frequency
    x1= x0+sum(w*x)/sum(w);
    x0p= x0;
    x0= long(floor(x1+0.5));
  } while (--repeat>0 && abs(x0-x0p)>=1);
  return x1;
}

func fft_centroid(a, repeat)
/* DOCUMENT fft_centroid(a)
       -or- fft_centroid(a, repeat)
     Return  the  position  of  centroid  of N-dimensional  array  A  assuming
     coordinates  along  dimensions  of  A  are  wrapped  as  in  a  FFT  (see
     fft_indgen).  The  algorithm proceeds by computing the  center of gravity
     of A around its  central element which is the maximum of  A for the first
     iteration  and  the  closest  to  the previously  computed  centroid  for
     subsequent  iterations.   The  maximum  number  of  iteration  is  REPEAT
     (default: 3;  in any  cases, at least  one iteration is  performed).  The
     Nyquist frequency along each even dimension is omitted to avoid a bias.

   SEE ALSO fft_indgen. */
{
  if (is_void(repeat)) repeat= 3;
  dims= dimsof(a);
  if ((ndims= dims(1)) <= 2) {
    if (ndims==2) return [_fft_centroid(a(,sum), repeat),
			 _fft_centroid(a(sum,), repeat)];
    if (ndims==1) return _fft_centroid(a, repeat);
    return 0.0; // ndims==0
  }
  result= array(double, ndims);
  for (i=1 ; i<=ndims ; ++i) {
    dim= dims(i+1);
    a1= array(double, numberof(a)/dim, dim);
    a1(*)= (i==ndims ? a(*) : transpose(a, [ndims,i])(*));
    a1= a1(sum,);
    result(i)= _fft_centroid(a1, repeat);
  }
  return result;
}

/*---------------------------------------------------------------------------*/
/* CENTERING / ROLLING */

func reverse_dims(a)
/* DOCUMENT reverse_dims(a)
     Returns array A with all its dimensions reversed.

   SEE ALSO: fft_recenter. */
{
  n = dimsof(a)(1);
  r = ::-1;
  if (n ==  1) return a(r);
  if (n ==  2) return a(r, r);
  if (n ==  3) return a(r, r, r);
  if (n ==  4) return a(r, r, r, r);
  if (n ==  5) return a(r, r, r, r, r);
  if (n ==  6) return a(r, r, r, r, r, r);
  if (n ==  7) return a(r, r, r, r, r, r, r);
  if (n ==  8) return a(r, r, r, r, r, r, r, r);
  if (n ==  9) return a(r, r, r, r, r, r, r, r, r);
  if (n == 10) return a(r, r, r, r, r, r, r, r, r, r);
  if (n == 11) return a(r, r, r, r, r, r, r, r, r, r, r);
  if (n == 12) return a(r, r, r, r, r, r, r, r, r, r, r, r);
  if (n == 13) return a(r, r, r, r, r, r, r, r, r, r, r, r, r);
  if (n == 14) return a(r, r, r, r, r, r, r, r, r, r, r, r, r, r);
  if (n == 15) return a(r, r, r, r, r, r, r, r, r, r, r, r, r, r, r);
  if (n == 16) return a(r, r, r, r, r, r, r, r, r, r, r, r, r, r, r, r);
  if (n == 17) return a(r, r, r, r, r, r, r, r, r, r, r, r, r, r, r, r, r);
  if (n == 18) return a(r, r, r, r, r, r, r, r, r, r, r, r, r, r, r, r, r, r);
  error, "too many dimensions";
}

func fft_recenter(x, template, reverse)
/* DOCUMENT fft_recenter(x, template)
       -or- fft_recenter(x, template, reverse)
     Returns array X rolled so that it matches most closely array TEMPLATE.
     X and TEMPLATE may be arrays  of real or complex numbers but must have
     the same dimension lists.  The  returned value is roll(X, S) where the
     offsets S minimize:

         sum(abs(roll(X, S) - TEMPLATE)^2)

     If optional  argument REVERSE is true,  X is also allowed  to have all
     its dimensions reversed in order to math TEMPLATE.

   SEE ALSO: fft, roll, reverse_dims. */
{
  ws = fft_setup(dimsof(x));
  conj_fft_template = conj(fft(template, +1, setup=ws));
  c1 = double(fft(conj_fft_template*fft(x, +1, setup=ws), -1, setup=ws));
  max_c1 = max(c1);
  if (reverse) {
    x2 = reverse_dims(x);
    c2 = double(fft(conj_fft_template*fft(x2, +1, setup=ws), -1, setup=ws));
    if ((max_c2 = max(c2)) > max_c1) {
      eq_nocopy, x, x2;
      i = where(c2 == max_c2);
    } else {
      i = where(c1 == max_c1);
    }
    c2 = [];
  } else {
    i = where(c1 == max_c1);
  }
  c1 = [];
  if (numberof(i) != 1) error, swrite("too many maxima (%d)", numberof(i));
  index = i(1) - 1; /* 0-based position of the maximum of correlation */
  dims = dimsof(x);
  ndims = dims(1);
  dims = dims(2:);
  offset = array(long, ndims);
  for (i=1 ; i<=ndims ; ++i) {
    dim = dims(i);
    offset(i) = index % dim;
    index /= dim;
  }
  //return roll(x, dims - offset);
  return (anyof(offset) ? roll(x, dims - offset) : x);
}

func fft_recenter_at_max(z, middle=)
/* DOCUMENT fft_recenter_at_max(z)

     Return Z  rolled so  that its element  with maximum value  (or maximum
     absolute value if  Z is complex) is at the  origin.  If keyword MIDDLE
     is true  (non-zero and non-nil) the  center is at the  middle of every
     dimension  otherwise the  center is  the first  element of  the output
     array (as assumed by the FFT).

   SEE ALSO: roll. */
{
  index = (structof(z) == complex ? abs(z) : z)(*)(mxx) - 1;
  dims = dimsof(z);
  ndims = dims(1);
  dims = dims(2:);
  offset = array(long, ndims);
  for (i=1 ; i<=ndims ; ++i) {
    dim = dims(i);
    offset(i) = index % dim;
    index /= dim;
  }
  if (middle) dims += (dims+1)/2;
  return roll(z, dims - offset);
}

func fft_roll_1d(a, off) {
  if ((dimlist = dimsof(a))(1) != 1) error, "expecting 1D array";
  n = dimlist(2);  k = (n + (off%n))%n; /* wrap offset in the range [0, n-1] */
  if (! k) return a;
  b = array(structof(a), dimlist);
  b(k+1:n) = a(1:n-k); b(1:k) = a(n-k+1:n); return b; }
func fft_roll_2d(a, off1, off2)
/* DOCUMENT fft_roll_1d(v, off)
       -or- fft_roll_2d(m, off1, off2)
     "rolls" dimensions of the vector V (1D array) or matrix M (2D array)
     and return a result with same data type than original array.

   SEE ALSO: roll. */
{
  if ((dimlist = dimsof(a))(1) != 2) error, "expecting 2D array";
  n1 = dimlist(2);
  k1 = (n1 + (off1%n1))%n1; /* wrap offset in the range [0, n1-1] */
  n2 = dimlist(3);
  k2 = (n2 + (off2%n2))%n2; /* wrap offset in the range [0, n2-1] */
  case = (! k1) + 2*(! k2);
  if (case == 3) return a;
  b = array(structof(a), dimlist);
  if (case == 0) {
    b((r1=k1+1:n1), (r2=k2+1:n2)) = a((s1=1:n1-k1), (s2=1:n2-k2));
    b(r1, (t2=1:k2)) = a(s1, (u2=n2-k2+1:n2));
    b((t1=1:k1), t2) = a((u1=n1-k1+1:n1), u2);
    b(t1, r2) = a(u1, s2);
  } else if (case == 1) {
    b( , k2+1:n2) = a( , 1:n2-k2);
    b( , 1:k2)    = a( , n2-k2+1:n2);
  } else {
    b(k1+1:n1, ) = a(1:n1-k1, );
    b(1:k1, )    = a(n1-k1+1:n1, );
  }
  return b;
}

/*---------------------------------------------------------------------------*/
/* FOURIER INTERPOLATION */

func fft_shift_phasor(off, u)
/* DOCUMENT fft_shift_phasor(off, dimlist)
       -or- fft_shift_phasor(off, fft_freqlist(dimlist))
     returns complex phasor to apply in FFT space for a shift by OFF cells
     in the real space.  DIMLIST is a list of dimensions -- the second
     calling sequence is to allow for computing the normalized FFT
     frequencies only once.  The offset OFF must have as many elements as
     PTR or as many as dimensions in DIMLIST (i.e. a shift for each
     dimension) and may be fractionnal.

     This function is intended for Fourier interpolation.  For instance,
     assuming A is 2-D real array:

         z = fft(a, +1);
         u = fft_freqlist(dimsof(a));
         q = fft_unphasor([0.33, -0.47], u);

     Then:

         fft_interp(z, q, real=1);

     yields the value of A interpolated at coordinate (0.33, -0.47) in FFT
     frame, i.e. center of lower left cell is at (0,0).  The shfited version
     of A by (0.33, -0.47) can be obtained by:

         fft_shift(z, q, real=1);


   SEE ALSO: fft_freqlist, fft_interp, fft_fine_shift, roll. */
{
  if (structof(u) != pointer) u = fft_freqlist(u);
  ndims = numberof(u);
  for (i=1 ; i<=ndims ; ++i) {
    a = (*u(i))*off(i);
    if (i == 1) p  =  cos(a) + 1i*sin(a);
    else        p *= (cos(a) + 1i*sin(a));
  }
  return p;
}

func fft_unphasor(z, phasor, setup=, real=) {
  z = fft(z*conj(phasor), -1, setup=setup);
  return (1.0/numberof(z))*(real ? double(z) : z); }
func fft_fine_shift(a, off, setup=)
/* DOCUMENT fft_fine_shift(a, off, setup=)
       -or- fft_unphasor(z, phasor, setup=, real=)
     The function fft_fine_shift returns array A shifted by offset OFF
     which can be fractionnal.  Alternatively, the function fft_unphasor
     can be used when the forward FFT of A and/or the complex phasor
     corresponding to the shift are already computed:

       fft_unphasor(fft(A, +1), fft_shift_phasor(OFF, dimsof(A)),
                    real=(structof(A) != complex))

     yields the same result as fft_fine_shift(A, OFF).

     These functions can make use of pre-computed FFT workspace specified
     by keyword SETUP (see fft_setup).

     TO-DO: Improve code by not Fourier transforming along direction
            where OFF is zero (or equal to an integer times the length
            of the dimension).

   SEE ALSO: fft, fft_setup, fft_shift_phasor, roll. */
{
  real = (structof(a) != complex);
  dims = dimsof(a);
  if (is_void(setup)) setup = fft_setup(dims);
  a = fft(fft(a, +1, setup=setup)*fft_shift_phasor(-off, dims), -1);
  return (1.0/numberof(a))*(real ? double(a) : a);
}

func fft_interp_real(z, phasor)
{ return (1.0/numberof(z))*double(sum(z*phasor)); }
func fft_interp_complex(z, phasor)
{ return (1.0/numberof(z))*sum(z*phasor); }
func fft_interp(a, off, setup=)
/* DOCUMENT fft_interp(a, off, setup=)
       -or- fft_interp_real(z, phasor)
       -or- fft_interp_complex(z, phasor)
     returns value obtained by Fourier interpolation of A at offset OFF.
     The function fft_interp computes the forward FFT of A and can make use
     of pre-computed FFT workspace specified by keyword SETUP.  The two
     other functions (fft_interp_complex, if A is complex; fft_interp_real
     otherwise) are usefull when the forward FFT of A and/or the complex
     phasor corresponding to the shift are already computed, their
     arguments are:

       Z = fft(A, +1);
       PHASOR = fft_shift_phasor(OFF, dimsof(A));

   SEE ALSO: fft, fft_fine_shift, fft_shift_phasor. */
{
  real = (structof(a) != complex);
  dims = dimsof(a);
  if (is_void(setup)) setup = fft_setup(dims);
  n = numberof(a);
  a = sum(fft(a, +1, setup=setup)*fft_shift_phasor(-off, dims));
  return (1.0/n)*(real ? double(a) : a);
}

/*---------------------------------------------------------------------------*/
/* GRAPHICS */

func fft_plh(y, scale=, legend=, hide=, type=, width=, color=, smooth=,
             marks=, marker=, mspace=, mphase=)
{ fft_plg, y, scale=scale, legend=legend, hide=hide, type=type, width=width,
    color=color, smooth=smooth, marks=marks, marker=marker, mspace=mspace,
    mphase=mphase, stair=1; }
func fft_plg(y, scale=, legend=, hide=, type=, width=, color=, smooth=,
             marks=, marker=, mspace=, mphase=, stair=)
/* DOCUMENT fft_plg, y;
        -or fft_plh, y;
     Plot 1-D FFT array Y as a curve, taking care of "rolling" Y and setting
     correct coordinates.  Keyword SCALE can be used to indicate the
     "frequel" scale along X-axis (SCALE is a scalar); by default,
     SCALE=1.0.

   KEYWORDS legend, hide, type, width, color, closed, smooth
            marks, marker, mspace, mphase.

   SEE ALSO plh, plg, roll. */
{
  if (is_void(scale)) scale= 1.0;
  else if (! is_array(scale) || dimsof(scale)(1)!=0)
    error, "expecting a scalar for SCALE";
  if (! is_array(y) || (dims= dimsof(y))(1)!=1) error, "expecting 1-D array";

  dim1= dims(2);
  min1= (max1= dim1/2) - dim1 + 1;

  if (stair) {
    // just=
    plh, roll(y, -min1), scale*indgen(min1:max1), legend=legend, hide=hide,
      type=type, width=width, color=color,
      marks=marks, marker=marker, mspace=mspace, mphase=mphase;
  } else {
    plg, roll(y, -min1), scale*indgen(min1:max1), legend=legend, hide=hide,
      type=type, width=width, color=color, smooth=smooth,
      marks=marks, marker=marker, mspace=mspace, mphase=mphase;
  }
}

func fft_pli(a, scale=, legend=, hide=, top=, cmin=, cmax=)
/* DOCUMENT fft_pli, a;
     Plot 2-D FFT array A as an image, taking care of "rolling" A and setting
     correct world boundaries.  Keyword SCALE can be used to indicate the
     "frequel" scale along both axis (SCALE is a scalar) or along each axis
     (SCALE is a 2-element vector: SCALE=[XSCALE,YSCALE]); by default,
     SCALE=[1.0, 1.0].

   KEYWORDS legend, hide, top, cmin, cmax.

   SEE ALSO pli, fft_roll_2d. */
{
  local scale1, dim1, min1, max1, scale2, dim2, min2, max2;
  __fft_pl2d_limits, a, scale;
  pli, fft_roll_2d(bytscl(a, top=top, cmin=cmin, cmax=cmax), -min1, -min2),
    scale1*(min1 - 0.5), scale2*(min2 - 0.5),
    scale1*(max1 + 0.5), scale2*(max2 + 0.5),
    legend=legend, hide=hide;
}

func fft_plc(a, scale=, levs=, type=, width=, color=, smooth=,
             legend=, hide=, marks=, marker=, mspace=, mphase=)
/* DOCUMENT fft_plc, a;
     Plot contour levels  of a 2-D FFT array A, taking  care of "rolling" A
     and setting  correct world boundaries.   Keyword SCALE can be  used to
     indicate the  "frequel" scale along both  axis (SCALE is  a scalar) or
     along each axis (SCALE  is a 2-element vector: SCALE=[XSCALE,YSCALE]);
     by default, SCALE=[1.0, 1.0].  Other  keywords have same meaning as in
     plc routine.

   KEYWORDS scale, levs, type, width, color, smooth,
            legend, hide, marks, marker, mspace, mphase.

   SEE ALSO plc, roll, fft_plfc. */
{
  local scale1, dim1, min1, max1, scale2, dim2, min2, max2;
  __fft_pl2d_limits, a, scale;
  u1 = (scale1*indgen(min1:max1))(,-:1:dim2);
  u2 = (scale2*indgen(min2:max2))(-:1:dim1,);
  plc, roll(a, [-min1, -min2]), u2, u1,
    levs=levs, type=type, width=width, color=color, smooth=smooth,
    legend=legend, hide=hide,
    marks=marks, marker=marker, mspace=mspace, mphase=mphase;
}

func fft_plfc(a, scale=, levs=, colors=)
/* DOCUMENT fft_plfc, a;
     Plot  filled contour  levels of  a  2-D FFT  array A,  taking care  of
     "rolling" A  and setting correct world boundaries.   Keyword SCALE can
     be used  to indicate the "frequel"  scale along both axis  (SCALE is a
     scalar)   or  along   each  axis   (SCALE  is   a   2-element  vector:
     SCALE=[XSCALE,YSCALE]); by default,  SCALE=[1.0, 1.0].  Other keywords
     have  same meaning  as in  plfc routine.   As with  plfc  routine, the
     actual level values get saved in external symbol plfc_levs.

   KEYWORDS scale, levs, colors.

   SEE ALSO plc, roll, fft_plc. */
{
  local scale1, dim1, min1, max1, scale2, dim2, min2, max2;
  __fft_pl2d_limits, a, scale;
  u1 = (scale1*indgen(min1:max1))(,-:1:dim2);
  u2 = (scale2*indgen(min2:max2))(-:1:dim1,);
  plfc, roll(a, [-min1, -min2]), u2, u1, levs=levs, colors=colors;
}

func __fft_pl2d_limits(z, scale)
/* DOCUMENT __fft_pl2d_limits, z, scale;
     Private routine used by fft_pli, fft_plc and fft_plfc.

   SEE ALSO fft_pli, fft_plc, fft_plfc. */
{
  extern scale1, dim1, min1, max1, scale2, dim2, min2, max2;
  if (is_void(scale)) {
    scale1 = scale2 = 1.0;
  } else if (dimsof(scale)(1) == 0) {
    scale1 = scale2 = scale;
  } else if (numberof(scale) == 2) {
    scale1 = scale(1);
    scale2 = scale(2);
  } else {
    error, "bad number of elements in SCALE";
  }
  if ((dims = dimsof(z))(1) != 2) error, "expecting 2-D array";
  dim1 = dims(2);
  min1 = (max1 = dim1/2) - dim1 + 1;
  dim2 = dims(3);
  min2 = (max2 = dim2/2) - dim2 + 1;
}

/*---------------------------------------------------------------------------*/

func fft_convolve(orig, psf, do_not_roll)
/* DOCUMENT fft_convolve(orig, psf);
       -or- fft_convolve(orig, psf, do_not_roll);
     Return discrete convolution (computed by FFT) of array ORIG by point
     spread function PSF.  Unless argument DO_NOT_ROLL is true, PSF is
     rolled before.  Note: ORIG and PSF must have same dimension list.

   SEE ALSO: fft, fft_setup, roll. */
{
  real = (structof(orig) != complex && structof(psf) != complex);
  dims = dimsof(orig);
  ws = fft_setup(dims);
  p = fft(orig, -1, setup=ws) * fft((do_not_roll?psf:roll(psf)), -1, setup=ws);
  orig = psf = []; // possibly free some memory
  fft_inplace, p, +1, setup=ws;
  if (real) p = double(p);
  return (1.0/numberof(p))*p;
}

#if 0
func fft_of_two_real_arrays(a, b, &ft_a, &ft_b, ljdir, rjdir, setup=)
/* DOCUMENT fft_of_two_real_arrays, a, b, ft_a, ft_b, direction;
       -or- fft_of_two_real_arrays, a, b, ft_a, ft_b, ljdir, rjdir;
     Computes the FFT  of arrays A and  B and stores them in  TF_A and FT_B
     respectively.  A and B must have same dimension list.  A single FFT is
     needed.  Agrguments  DIRECTION, LJDIR,  RJDIR, and keyword  SETUP have
     the same meaning as for the fft function (which see).

   SEE ALSO: fft_setup, fft_inplace. */
{
  if (structof(a) == complex || structof(b) == complex)
    error, "A and B must be non-complex";
  c = a + 1i*b;
  fft_inplace, c, ljdir, rjdir, setup=setup;
  b = c(fft_symmetric_index(dimsof(c)));
  a = c + b;
  b = c - b;
  ft_a = 0.5*double(a) + 0.5i*(b.im);
  ft_b = 0.5*(a.im)    - 0.5i*double(b);
}
#endif

/*---------------------------------------------------------------------------*/

/* Notes:
 *   there are 1 + n/2 "positive" frequencies
 *   there are (n - 1)/2 "negative" frequencies
 */

func fft_paste(a, b)
/* DOCUMENT fft_paste(a, b)
 *     -or- fft_paste, a, b;
 *   Paste array B into array A in the sense of FFT indexing.  All
 *   dimensions of A must be greater or equal the corresponding dimension
 *   of B.  When called as a subroutine, the operation is done in-place.
 *
 * RESTRICTIONS:
 *   For even dimensions, the Nyquist frequency from B is not pasted
 *   into A.
 *
 * SEE ALSO: fft_indgen.
 */
{
  if (! is_array(a) || ! is_array(b))
    error, "expecting array argument(s)";
  adim = dimsof(a);
  bdim = dimsof(b);
  if ((n = adim(1) - bdim(1)) != 0) {
    if (n > 0) {
      grow, bdim, array(1L, n);
      bdim(1) = adim(1);
    } else {
      grow, adim, array(1L, -n);
      adim(1) = bdim(1);
    }
  }
  if (anyof(adim < bdim)) error, "destination array is too small";
  n = numberof(adim);
  ia = ib = 1L; /* indices start at one in Yorick */
  sa = numberof(a); /* stride in A */
  sb = numberof(b); /* stride in B */
  for (k=n ; k>=2 ; --k) {
    alen = adim(k);
    sa /= alen;
    blen = bdim(k);
    sb /= blen;
    if (blen >= 3) {
      j = (blen - 1)/2; /* maximum absolute frequency */
      ja = jb = array(long, 2*j + 1);
      ja(1:j+1) = indgen(0 : j*sa : sa);
      jb(1:j+1) = indgen(0 : j*sb : sb);
      ja(j+2:) = indgen((alen - j)*sa : (alen - 1)*sa : sa);
      jb(j+2:) = indgen((blen - j)*sb : (blen - 1)*sb : sb);
    } else {
      ja = jb = 0L;
    }
    if (k == n) {
      ia += ja;
      ib += jb;
    } else {
      ia = ja + ia(-,..);
      ib = jb + ib(-,..);
    }
  }
  if (! am_subroutine()) a = a; /* make a copy */
  a(ia) = b(ib);
  return a;
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
