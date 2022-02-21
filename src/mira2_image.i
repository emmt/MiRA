/*
 * mira2_image.i -
 *
 * Management of images and FITS files in MiRA.
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

func mira_soft_threshold(img, lvl, nrm)
/* DOCUMENT mira_soft_threshold(img, lvl);
         or mira_soft_threshold(img, lvl, nrm);

     Perform soft-threshlding of image IMG. Argument LVL is the threshold level
     specified in terms of the fraction of the non-zero pixels.  For instance,
     LVL = 0.05 means that the threshold will be such that 5% of the less
     bright pixels (in absolute value) will be set to zero.

     If optional argument NRM is specified, the result is rescaled so that its
     sum is equal to NRM.  Rescaling is only applied if the sum of
     soft-thresholded pixels is strictly positive.
 */
{
  j = where(img);
  if (is_array(j)) {
    n = numberof(j);
    v = abs(img(j));
    l = quick_select(v, lround(1 + (n - 1)*lvl));
    tmp = abs(img) - l;
    if (max(tmp) <= 0) {
      warn, "Thresholding is not applied to avoid a zero-filled image";
    } else {
      img = sign(img)*max(0.0, tmp);
    }
  }
  if (! is_void(nrm)) {
    s = double(sum(img));
    if (s > 0.0) {
      img = (double(nrm)/s)*unref(img);
    }
  }
  return img;
}

func mira_resample_image(src, pad=, norm=,
                         naxis1=, crpix1=, crval1=, cdelt1=, cunit1=,
                         naxis2=, crpix2=, crval2=, cdelt2=, cunit2=,
                         naxis3=, crpix3=, crval3=, cdelt3=, cunit3=)
/* DOCUMENT mira_resample_image(src, key=val, ...);

     Resample an image (as returned by `mira_read_image`) to given parameters,
     all passed as keywords and with default values taken from SRC.  Number of
     samples, index of reference position, coordinate at reference position and
     increment along n-th axis may be set with keywords NAXISn, CRPIXn, CRVALn
     and CDELTn respectively.  It is assumed that keywords CRVALn and CDELTn
     are in units given by keyword CUNITn if specified, in units given by
     SRC.CUNITn otherwise.  Units (via keywords CUNITn, with n the axis number)
     may be changed but must be consistent (a lenght cannot be converted into
     an angle).

     Keywords PAD and NORM can be used to specify the padding value for
     extrapolating pixels and whether or not to preserve the normalization.

   SEE ALSO: mira_resample_axis.
 */
{
  rdif = mira_rdif;
  eps = 1e-14;
  local arr;
  eq_nocopy, arr, src.arr;
  naxis = src.naxis;
  dst = h_new(naxis=naxis);
  if (naxis >= 1) {
    if (is_void(cunit1)) {
      cunit1 = src.cunit1;
    }
    if (is_scalar(cunit1) && is_string(cunit1)) {
      factor = mira_convert_units(src.cunit1, cunit1);
    } else {
      error, "bad value for CUNIT1";
    }
    if (is_void(naxis1)) {
      naxis1 = src.naxis1;
    }
    if (is_scalar(naxis1) && is_integer(naxis1) && naxis1 >= 1) {
      naxis1 = long(naxis1)
    } else {
      error, "bad value for NAXIS1";
    }
    if (is_void(crpix1)) {
      crpix1 = src.crpix1;
    }
    if (is_scalar(crpix1) && identof(crpix1) <= Y_DOUBLE) {
      crpix1 = double(crpix1);
    } else {
      error, "bad value for CRPIX1";
    }
    if (is_void(crval1)) {
      crval1 = factor*src.crval1;
    }
    if (is_scalar(crval1) && identof(crval1) <= Y_DOUBLE) {
      crval1 = double(crval1);
    } else {
      error, "bad value for CRVAL1";
    }
    if (is_void(cdelt1)) {
      cdelt1 = factor*src.cdelt1;
    }
    if (is_scalar(cdelt1) && identof(cdelt1) <= Y_DOUBLE && cdelt1 != 0) {
      cdelt1 = double(cdelt1);
    } else {
      error, "bad value for CDELT1";
    }
    h_set, dst, naxis1 = naxis1, crpix1 = crpix1, crval1 = crval1,
      cdelt1 = cdelt1, cunit1 = cunit1, ctype1 = src.ctype1;
    if (dst.naxis1 != src.naxis1 ||
        rdif(dst.cdelt1, factor*src.cdelt1) > eps ||
        rdif(dst.crval1 - dst.crpix1*dst.cdelt1,
             factor*(src.crval1 - src.crpix1*src.cdelt1)) > eps) {
      arr = mira_resample_axis(arr, src.crpix1, factor*src.crval1,
                               factor*src.cdelt1, dst.naxis1,
                               dst.crpix1, dst.crval1, dst.cdelt1,
                               1, pad=pad, norm=norm);
    }
  }
  if (naxis >= 2) {
    if (is_void(cunit2)) {
      cunit2 = src.cunit2;
    }
    if (is_scalar(cunit2) && is_string(cunit2)) {
      factor = mira_convert_units(src.cunit2, cunit2);
    } else {
      error, "bad value for CUNIT2";
    }
    if (is_void(naxis2)) {
      naxis2 = src.naxis2;
    }
    if (is_scalar(naxis2) && is_integer(naxis2) && naxis2 >= 2) {
      naxis2 = long(naxis2)
    } else {
      error, "bad value for NAXIS2";
    }
    if (is_void(crpix2)) {
      crpix2 = src.crpix2;
    }
    if (is_scalar(crpix2) && identof(crpix2) <= Y_DOUBLE) {
      crpix2 = double(crpix2);
    } else {
      error, "bad value for CRPIX2";
    }
    if (is_void(crval2)) {
      crval2 = factor*src.crval2;
    }
    if (is_scalar(crval2) && identof(crval2) <= Y_DOUBLE) {
      crval2 = double(crval2);
    } else {
      error, "bad value for CRVAL2";
    }
    if (is_void(cdelt2)) {
      cdelt2 = factor*src.cdelt2;
    }
    if (is_scalar(cdelt2) && identof(cdelt2) <= Y_DOUBLE && cdelt2 != 0) {
      cdelt2 = double(cdelt2);
    } else {
      error, "bad value for CDELT2";
    }
    h_set, dst, naxis2 = naxis2, crpix2 = crpix2, crval2 = crval2,
      cdelt2 = cdelt2, cunit2 = cunit2, ctype2 = src.ctype2;
    if (dst.naxis2 != src.naxis2 ||
        rdif(dst.cdelt2, factor*src.cdelt2) > eps ||
        rdif(dst.crval2 - dst.crpix2*dst.cdelt2,
             factor*(src.crval2 - src.crpix2*src.cdelt2)) > eps) {
      arr = mira_resample_axis(arr, src.crpix2, factor*src.crval2,
                               factor*src.cdelt2, dst.naxis2,
                               dst.crpix2, dst.crval2, dst.cdelt2,
                               2, pad=pad, norm=norm);
    }
  }
  if (naxis >= 3) {
    if (is_void(cunit3)) {
      cunit3 = src.cunit3;
    }
    if (is_scalar(cunit3) && is_string(cunit3)) {
      factor = mira_convert_units(src.cunit3, cunit3);
    } else {
      error, "bad value for CUNIT3";
    }
    if (is_void(naxis3)) {
      naxis3 = src.naxis3;
    }
    if (is_scalar(naxis3) && is_integer(naxis3) && naxis3 >= 3) {
      naxis3 = long(naxis3)
    } else {
      error, "bad value for NAXIS3";
    }
    if (is_void(crpix3)) {
      crpix3 = src.crpix3;
    }
    if (is_scalar(crpix3) && identof(crpix3) <= Y_DOUBLE) {
      crpix3 = double(crpix3);
    } else {
      error, "bad value for CRPIX3";
    }
    if (is_void(crval3)) {
      crval3 = factor*src.crval3;
    }
    if (is_scalar(crval3) && identof(crval3) <= Y_DOUBLE) {
      crval3 = double(crval3);
    } else {
      error, "bad value for CRVAL3";
    }
    if (is_void(cdelt3)) {
      cdelt3 = factor*src.cdelt3;
    }
    if (is_scalar(cdelt3) && identof(cdelt3) <= Y_DOUBLE && cdelt3 != 0) {
      cdelt3 = double(cdelt3);
    } else {
      error, "bad value for CDELT3";
    }
    h_set, dst, naxis3 = naxis3, crpix3 = crpix3, crval3 = crval3,
      cdelt3 = cdelt3, cunit3 = cunit3, ctype3 = src.ctype3;
    if (dst.naxis3 != src.naxis3 ||
        rdif(dst.cdelt3, factor*src.cdelt3) > eps ||
        rdif(dst.crval3 - dst.crpix3*dst.cdelt3,
             factor*(src.crval3 - src.crpix3*src.cdelt3)) > eps) {
      arr = mira_resample_axis(arr, src.crpix3, factor*src.crval3,
                               factor*src.cdelt3, dst.naxis3,
                               dst.crpix3, dst.crval3, dst.cdelt3,
                               3, pad=pad, norm=norm);
    }
  }
  return h_set(dst, arr=arr);
}

func mira_resample_axis(arr, crpix1, crval1, cdelt1,
                        n2, crpix2, crval2, cdelt2,
                        which, pad=, norm=)
/* DOCUMENT mira_resample_axis(arr, crpix1, crval1, cdelt1,
                                n2, crpix2, crval2, cdelt2,
                                which);

     Resample array ARR along one of its axis.  Resampling is performed by
     finite difference of the integration of a piecewise linear approximation.
     This is suitable to preserve positivity and for downsampling or
     upsampling.

     Arguments CRPIX1, CRVAL1 and CDELT1 are respectively the index of the
     reference point, the coordinate of the reference point and the coordinate
     increment along the resampled axis in the source array ARR.

     Arguments N2, CRPIX2, CRVAL2 and CDELT2 are respectively the number of
     samples, the index of the reference point, the coordinate of the reference
     point and the coordinate increment along the resampled axis in the
     destination array.

     Arguments CRVAL1, CDELT1, CRVAL2 and CDELT2 are assumed to be given in the
     same units.

     Optional argument WHICH specify the axis to resample (the first one by
     default).

     Keyword PAD may be used to specify a value for extrapolated values (the
     nearest value is used by default).

     Keyword NORM may be set true to preserve the normalization of the result.

   SEE ALSO: integ.
 */
{
  if (! is_array(arr) || identof(arr) > Y_DOUBLE) {
    error, "expecting an array of reals";
  }
  dims = dimsof(arr);
  rank = numberof(dims) - 1;
  if (is_void(which)) {
    which = 1;
  }
  if (which <= 0) {
    which += rank;
  }
  if (which <= 0 || which > rank) {
    error, "out of bound axis index";
  }
  if (! is_void(pad)) {
    dims(1 + which) += 2;
    tmp = array(double(pad), dims);
    /**/ if (which ==  1) tmp(2:-1,..) = arr;
    else if (which ==  2) tmp(,2:-1,..) = arr;
    else if (which ==  3) tmp(,,2:-1,..) = arr;
    else if (which ==  4) tmp(,,,2:-1,..) = arr;
    else if (which ==  5) tmp(,,,,2:-1,..) = arr;
    else if (which ==  6) tmp(,,,,,2:-1,..) = arr;
    else if (which ==  7) tmp(,,,,,,2:-1,..) = arr;
    else if (which ==  8) tmp(,,,,,,,2:-1,..) = arr;
    else if (which ==  9) tmp(,,,,,,,,2:-1,..) = arr;
    else if (which == 10) tmp(,,,,,,,,,2:-1,..) = arr;
    else error, "too many dimensions";
    eq_nocopy, arr, tmp;
    crpix1 += 1.0;
  }
  n1 = dims(1 + which);
  x1 = (double(indgen(n1)) - crpix1)*cdelt1 + crval1;
  x2 = (double(indgen(n2 + 1)) - (crpix2 + 0.5))*cdelt2 + crval2;
  arr = integ(arr, x1, x2, which);
  a = 1.0/(norm ? cdelt1 : cdelt2);
  if (which ==  1) return a*arr(dif,..);
  if (which ==  2) return a*arr(,dif,..);
  if (which ==  3) return a*arr(,,dif,..);
  if (which ==  4) return a*arr(,,,dif,..);
  if (which ==  5) return a*arr(,,,,dif,..);
  if (which ==  6) return a*arr(,,,,,dif,..);
  if (which ==  7) return a*arr(,,,,,,dif,..);
  if (which ==  8) return a*arr(,,,,,,,dif,..);
  if (which ==  9) return a*arr(,,,,,,,,dif,..);
  if (which == 10) return a*arr(,,,,,,,,,dif,..);
  error, "too many dimensions";
}

func mira_recenter(x, quiet=)
/* DOCUMENT mira_recenter(x);

     Recenter model image X at its photo-centre rounded to nearest pixel.
     Argument X must have at least 2 dimensions, the first 2 dimensions of X
     are considered to be the angular direction.  Extra dimensions are ignored
     for the recentering (they can represent other coordinates for instance
     the wavelength or the time).  Along a dimensions of lenght N, the center
     is at (N - N/2)-th pixel -- with integer division -- which corresponds to
     the model of the Fourier transform assumed by MiRA.

     Unless keyword QUIET is true, the coordinates of the center get printed
     out.

   SEE ALSO: mira_solve.
 */
{
  sx = sum(x);
  if (sx <= 0.0) {
    return x;
  }

  dimlist = dimsof(x);
  n1 = dimlist(2);
  n2 = dimlist(3);
  o1 = n1 - (n1/2);
  o2 = n2 - (n2/2);
  c1 = sum(double(indgen(n1) - o1)     * x)/sx;
  c2 = sum(double(indgen(n2) - o2)(-,) * x)/sx;
  if (! quiet) {
    inform, "Offsets of photo-center: (%+.1f, %+.1f) pixels.\n", c1, c2;
  }
  i1 = lround(c1);
  i2 = lround(c2);
  if (i1 == 0 && i2 == 0) {
    return x;
  }
  if (i1 > 0) {
    dst1 = 1:n1-i1;
    src1 = 1+i1:n1;
  } else {
    dst1 = 1-i1:n1;
    src1 = 1:n1+i1;
  }
  if (i2 > 0) {
    dst2 = 1:n2-i2;
    src2 = 1+i2:n2;
  } else {
    dst2 = 1-i2:n2;
    src2 = 1:n2+i2;
  }
  xp = array(structof(x), dimsof(x));
  xp(dst1, dst2, ..) = x(src1, src2, ..);
  return xp;
}

/*---------------------------------------------------------------------------*/
/* SAVE/LOAD IMAGES */

func mira_read_image(inp, start=, hduname=, extname=)
/* DOCUMENT img = mira_read_image(inp);

     Read an image in input FITS file INP which can be a file name or a FITS
     handle.

     Keyword START specifies from which FITS Header Data Unit (HDU) the search
     of an image should start.  Default is to start at the current HDU.

     If keyword EXTNAME is specified, the image must be in a FITS extension
     whose name is EXTNAME.

     If keyword HDUNAME is specified, the header of the image must have a FITS
     keyword "HDUNAME" matching value HDUNAME.

     The first FITS image matching all the above criteria is returned.

     The returned value is a hash table with the following members:

       img.arr    = array of pixel values (2D or 3D)
       img.naxis  = number of dimensions

       img.naxis1 = length of the 1st axis
       img.crval1 = 1st coordinate of the reference pixel
       img.crpix1 = index along 1st axis of the reference pixel
       img.cdelt1 = coordinate increment along 1st axis
       img.ctype1 = coordinate name along 1st axis
       img.cunit1 = coordinate units along 1st axis

       img.naxis2 = length of the 2nd axis
       img.crval2 = 2nd coordinate of the reference pixel
       img.crpix2 = index along 2nd axis of the reference pixel
       img.cdelt2 = coordinate increment along 2nd axis
       img.ctype2 = coordinate name along 2nd axis
       img.cunit2 = coordinate units along 2nd axis

     If image is 3D:

       img.naxis3 = length of the 3rd axis
       img.crval3 = 3rd coordinate of the reference pixel
       img.crpix3 = index along 3rd axis of the reference pixel
       img.cdelt3 = coordinate increment along 3rd axis
       img.ctype3 = coordinate name along 3rd axis
       img.cunit3 = coordinate units along 3rd axis

     The conventions for MiRA are to use the first 2 axes for spatial
     coordinates while the 3rd axis, if any, is for the spectral channel.  The
     values of members ctype1..3 are upper case strings with leading and
     trailing spaces stripped.  The values of members cunit1..3 are strings
     with leading and trailing spaces stripped.  Other c... members are scalar
     doubles and dimensions are scalar long.

     The pixels of the image may be reversed along some axes so that cdeltN
     (the increment along the N-th axis) is nonnegative for all N.


   SEE ALSO: mira_wrap_image, mira_save_image, fits_read_array. */
{
  /* Open FITS file if necessary. */
  local fh;
  if (is_string(inp)) {
    fh = fits_open(inp, 'r');
    close_on_exit = 1n;
  } else {
    eq_nocopy, fh, inp;
    close_on_exit = 0n;
  }

  /* Move to given HDU/EXTNAME/HDUNAME. */
  mira_find_fits_hdu, fh, "IMAGE", start=start,
    extname=extname, hduname=hduname;

  /* Read image, get coordinates (FIXME: deal with PCij and CDij) and fix
     orientation of axes. */
  naxis = fits_get(fh, "NAXIS");
  img = h_new(naxis=naxis, arr=fits_read_array(fh));
  for (i = 1; i <= naxis; ++i) {
    _mira_get_fits_axis, img, fh, i;
  }

  /* Close FITS file and return data. */
  if (close_on_exit) {
    fits_close, fh;
  }
  return img;
}

func _mira_get_fits_axis(img, fh, n, crval=, crpix=, cdelt=, ctype=, cunit=)
/* DOCUMENT _mira_get_fits_axis, img, fh, n;

     private subroutine to extract coordinate information along N-th axis of
     image IMG read in current HDU of FITS handle FH.


   SEE ALSO: mira_read_image.
 */
{
  sfx = swrite(format="%d", n);

  /* Get NAXISn */
  key = "NAXIS" + sfx;
  naxis = fits_get(fh, key);
  if (! is_scalar(naxis) || ! is_integer(naxis)) {
    error, "value of " + key + " must be an integer";
  }
  naxis += 0;

  /* Get CRVALn */
  key = "CRVAL" + sfx;
  value = fits_get(fh, key);
  if (! is_void(value)) {
    crval = value;
  } else if (is_void(crval)) {
    crval = 0.0;
  }
  if (! is_scalar(crval) || identof(crval) > Y_DOUBLE) {
    error, "value of " + key + " must be a real";
  }
  crval += 0.0;

  /* Get CRPIXn */
  key = "CRPIX" + sfx;
  value = fits_get(fh, key);
  if (! is_void(value)) {
    crpix = value;
  } else if (is_void(crpix)) {
    crpix = double(mira_central_index(naxis));
  }
  if (! is_scalar(crpix) || identof(crpix) > Y_DOUBLE) {
    error, "value of " + key + " must be a real";
  }
  crpix += 0.0;

  /* Get CDELTn */
  key = "CDELT" + sfx;
  value = fits_get(fh, key);
  if (! is_void(value)) {
    cdelt = value;
  } else if (is_void(cdelt)) {
    cdelt = 1.0;
  }
  if (! is_scalar(cdelt) || identof(cdelt) > Y_DOUBLE) {
    error, "value of " + key + " must be a real";
  }
  cdelt += 0.0;

  /* Get CTYPEn */
  key = "CTYPE" + sfx;
  value = fits_get(fh, key);
  if (! is_void(value)) {
    ctype = value;
  } else if (is_void(ctype)) {
    ctype = "";
  }
  if (! is_scalar(ctype) || ! is_string(ctype)) {
    error, "value of " + key + " must be a string";
  }
  ctype = strcase(1, strtrim(ctype, 3));

  /* Get CUNITn */
  key = "CUNIT" + sfx;
  value = fits_get(fh, key);
  if (! is_void(value)) {
    cunit = value;
  } else if (is_void(cunit)) {
    cunit = "";
  }
  if (! is_scalar(cunit) || ! is_string(cunit)) {
    error, "value of " + key + " must be a string";
  }
  cunit = strtrim(cunit, 3);

  /* Fix orientation. */
  if (cdelt < 0) {
    /**/ if (n ==  1) h_set, img, arr=img.arr(::-1,..);
    else if (n ==  2) h_set, img, arr=img.arr(,::-1,..);
    else if (n ==  3) h_set, img, arr=img.arr(,,::-1,..);
    else if (n ==  4) h_set, img, arr=img.arr(,,,::-1,..);
    else if (n ==  5) h_set, img, arr=img.arr(,,,,::-1,..);
    else if (n ==  6) h_set, img, arr=img.arr(,,,,,::-1,..);
    else if (n ==  7) h_set, img, arr=img.arr(,,,,,,::-1,..);
    else if (n ==  8) h_set, img, arr=img.arr(,,,,,,,::-1,..);
    else if (n ==  9) h_set, img, arr=img.arr(,,,,,,,,::-1,..);
    else if (n == 10) h_set, img, arr=img.arr(,,,,,,,,,::-1,..);
    cdelt = -cdelt;
    crpix = naxis + 1 - crpix;
  }

  /* Udate information. */
  return h_set(img,
               "naxis"+sfx, naxis,
               "crval"+sfx, crval,
               "crpix"+sfx, crpix,
               "cdelt"+sfx, cdelt,
               "ctype"+sfx, ctype,
               "cunit"+sfx, cunit);
}

func mira_wrap_image(arr, dat)
/* DOCUMENT img = mira_wrap_image(arr);
         or img = mira_wrap_image(arr, dat);

     wraps array ARR in an image with the same structure as the one returned by
     `mira_read_image`.  Optional argument DAT is MiRA data instance needed to
     retrieve image settings such as the pixel size.

   SEE ALSO: mira_read_image, mira_save_image.
 */
{
  if (! is_array(arr) || identof(arr) > Y_DOUBLE) {
    error, "expecting an array of non-complex numerical values";
  }
  dims = dimsof(arr);
  naxis = dims(1);
  if (naxis != 2 && naxis != 3) {
    error, "expecting a 2D or 3D array";
  }
  naxis1 = dims(2);
  naxis2 = dims(3);
  naxis3 = (naxis >= 3 ? dims(4) : 1);
  if (is_hash(dat)) {
    local wave, x, y;
    eq_nocopy, wave, mira_image_wave(dat);
    eq_nocopy, x, mira_image_x(dat);
    eq_nocopy, y, mira_image_y(dat);
    if (naxis3 == 1) {
      wave = (min(wave) + max(wave))/2.0;
    }
    if (naxis1 != numberof(x) ||
        naxis2 != numberof(y) ||
        naxis3 != numberof(wave)) {
      error, "image dimensions incompatible with MiRA data";
    }
    i1 = mira_central_index(naxis1);
    i2 = mira_central_index(naxis2);
    DEGREE = MIRA_DEGREE;
    NANOMETER = MIRA_NANOMETER;
    img = h_new(arr = arr,
                naxis = naxis,
                naxis1 = naxis1,
                crpix1 = double(i1),
                crval1 = x(i1)/DEGREE,
                cdelt1 = avg(x(dif))/DEGREE,
                ctype1 = "RA---TAN",
                cunit1 = "deg",
                naxis2 = naxis2,
                crpix2 = double(i2),
                crval2 = y(i2)/DEGREE,
                cdelt2 = avg(y(dif))/DEGREE,
                ctype2 = "DEC--TAN",
                cunit2 = "deg");
    if (naxis < 3) {
      return h_set(img, wavelength=wave(1));
    }
    return h_set(img,
                 naxis3 = naxis3,
                 crpix3 = 1.0,
                 crval3 = wave(1)/NANOMETER,
                 cdelt3 = (numberof(wave) > 1 ? avg(wave(dif)) : 1.0)/NANOMETER,
                 ctype3 = "WAVELENGTH",
                 cunit3 = "nm");
  } else if (is_void(dat)) {
    i1 = mira_central_index(naxis1);
    i2 = mira_central_index(naxis2);
    img = h_new(arr = arr,
                naxis = naxis,
                naxis1 = naxis1,
                crpix1 = double(i1),
                crval1 = 0.0,
                cdelt1 = 1.0,
                ctype1 = "X",
                cunit1 = "",
                naxis2 = naxis2,
                crpix2 = double(i2),
                crval2 = 0.0,
                cdelt2 = 1.0,
                ctype2 = "Y",
                cunit2 = "");
    if (naxis < 3) {
      return img;
    }
    return h_set(img,
                 naxis3 = naxis3,
                 crpix3 = 1.0,
                 crval3 = 1.0,
                 cdelt3 = 1.0,
                 ctype3 = "SPECTRAL CHANNEL",
                 cunit3 = "");
  } else {
    error, "optional second argument must be a MiRA data instance";
  }
}

func mira_save_image(img, dest, overwrite=, bitpix=, data=, hook=,
                     comment=, history=, extname=, hduname=,
                     save_visibilities=)
/* DOCUMENT mira_save_image, img, dest;
         or fh = mira_save_image(img, dest);

     Save image IMG in FITS file DEST (can be a file name or a FITS handle).
     Image IMG can be a structured object as returned by `mira_read_image` or
     `mira_wrap_image` or an array of pixel values.  In the latter case,
     keyword DATA can be set with the MiRA data instance from which the image
     has been reconstructed (see `mira_wrap_image`).  When called as a
     function, the FITS handle is returned and can be used, for instance, to
     append more FITS extensions.

     Keywords COMMENT and HISTORY can be set with an array of strings to
     specify comments and history records.

     For more flexibility on the contents of the FITS header, keyword HOOK may
     be set with a callable object called as:

         hook, fh;

     with FH the FITS handle and in order to add custom FITS keywords to the
     header before writing the data.

     If DEST is a string, keyword OVERWRITE can be set true to allow for
     overwriting file DEST if it already exists.

     Keyword EXTNAME can be used to specify the name of the FITS extension.

     Keyword HDUNAME can be used to specify the name of the FITS HDU.

     Keyword BITPIX (by default -32, that is single precision floating point)
     can be used to specify the binary format of pixel values written in the
     file.


   SEE ALSO: fits, mira_read_image, mira_wrap_image.
 */
{
  /* Check input image. */
  local img_arr;
  if (is_array(img)) {

    inform, " Saving image..."
    eq_nocopy, img_arr, img;
    img = mira_wrap_image(unref(img), data);
  } else if (is_hash(img)) {
    if (! is_void(data)) {
      write, format="WARNING - %s\n", "MiRA data ignored";
    }
    eq_nocopy, img_arr, img.arr;
  } else {
    error, "unexpected image type";
  }

  naxis = img.naxis;
  if (naxis != 2 && naxis != 3) {
    error, "expecting a 2D or 3D image";
  }

  /* Get FITS handle. */
  local fh;
  if (is_void(bitpix)) {
    bitpix = -32;
  }
  if (is_string(dest)) {
    fh = fits_open(dest, 'w', overwrite=overwrite);
  } else {
    eq_nocopy, fh, dest;
    fits_new_hdu, fh, "IMAGE", "Image extension";
  }
  hdu = fits_current_hdu(fh);
  if (hdu == 1) {
    fits_set, fh, "SIMPLE", 'T', "true FITS file";
  }
  fits_set, fh, "BITPIX", bitpix, "bits per pixel";
  fits_set, fh, "NAXIS",  naxis,  "number of dimensions";
  for (k = 1; k <= naxis; ++k) {
    sfx = swrite(format="%d", k);
    fits_set, fh, "NAXIS"+sfx, h_get(img, "naxis"+sfx),
      "number of elements along axis";
  }
  fits_set, fh, "EXTEND", 'T', "this file may contain FITS extensions";

  /* Save axis information. Manage to have the image correctly displayed with
     most viewers (East toward left and North toward top).  */
  if ((img.ctype1 == "RA---TAN" && img.cdelt1 > 0.0) ||
      (img.ctype2 == "DEC--TAN" && img.cdelt2 < 0.0)) {
    /* Modify orientation paramaters. */
    local flip1, flip2;
    cdelt1 = img.cdelt1;
    crpix1 = img.crpix1;
    if (img.ctype1 == "RA---TAN" && cdelt1 > 0.0) {
      flip1 = ::-1;
      cdelt1 = -cdelt1;
      crpix1 = img.naxis1 - crpix1 + 1;
    }
    cdelt2 = img.cdelt2;
    crpix2 = img.crpix2;
    if (img.ctype2 == "DEC--TAN" && cdelt2 < 0.0) {
      flip2 = ::-1;
      cdelt2 = -cdelt2;
      crpix2 = img.naxis2 - crpix2 + 1;
    }

    /* Clone image structure before modifying it. */
    keys = h_keys(img);
    cpy = h_new();
    for (k = 1; k <= numberof(keys); ++k) {
      key = keys(k);
      h_set, cpy, key, h_get(img, key);
    }
    img = h_set(cpy, arr=img.arr(flip1,flip2,..),
                cdelt1=cdelt1, crpix1=crpix1,
                cdelt2=cdelt2, crpix2=crpix2);
  }
  for (k = 1; k <= naxis; ++k) {
    sfx = swrite(format="%d", k);
    fits_set, fh, "CRPIX"+sfx, h_get(img, "crpix"+sfx),
      "coordinate system reference pixel";
    fits_set, fh, "CRVAL"+sfx, h_get(img, "crval"+sfx),
      "coordinate system value at reference pixel";
    fits_set, fh, "CDELT"+sfx, h_get(img, "cdelt"+sfx),
      "coordinate increment along axis";
    fits_set, fh, "CTYPE"+sfx, h_get(img, "ctype"+sfx),
      "name of the coordinate axis";
    fits_set, fh, "CUNIT"+sfx, h_get(img, "cunit"+sfx),
      "units of the coordinate axis";
  }
  if (hdu > 1 && ! is_void(extname)) {
    fits_set, fh, "EXTNAME", extname, "Name of this HDU";
  }
  if (! is_void(hduname)) {
    fits_set, fh, "HDUNAME", hduname, "Name of this HDU";
  }
  if (! is_void(comment)) {
    for (k = 1; k <= numberof(comment); ++k) {
      fits_set, fh, "COMMENT", comment(k);
    }
  }
  if (! is_void(history)) {
    for (k = 1; k <= numberof(history); ++k) {
      fits_set, fh, "HISTORY", history(k);
    }
  }

  /* Call the hook, if any, to edit the FITS header. */
  if (! is_void(hook)) {
    hook, fh;
  }

  /* Write the header, the data and pad the HDU. */
  fits_write_header, fh;
  fits_write_array, fh, img.arr;
  fits_pad_hdu, fh;


  return fh;
}

/*---------------------------------------------------------------------------*/
/* UTILITIES FOR FITS FILES */

local _MIRA_FITS_TRUE, _MIRA_FITS_FALSE;
local mira_get_fits_string, mira_get_fits_logical;
local mira_get_fits_integer, mira_get_fits_real;
/* DOCUMENT  str = mira_get_fits_string(fh, kwd);
         or bool = mira_get_fits_logical(fh, kwd);
         or ival = mira_get_fits_integer(fh, kwd);
         or rval = mira_get_fits_real(fh, kwd);

     These functions retrieve the value of FITS keyword KWD from the FITS
     handle FH.  Returned values have a specific type.  If keyword does not
     exist, nil is returned.  If keyword does exist but has an incompatible
     value, an error is thrown.  In the above examples:

     - STR is either [] or a string in upper case letters with trailing spaces
       removed.

     - BOOL is either [], 1n (true) or 0n (false).

     - IVAL is either [] or an integer (long) value.

     - RVAL is either [] or a floating-point (double) value.


   SEE ALSO: fits_get.
 */
func mira_get_fits_string(fh, kwd, lower=) /* DOCUMENTED */
{
  value = fits_get(fh, kwd);
  if (is_scalar(value) && structof(value) == string) {
    value = strcase((lower ? 0n : 1n), strtrim(value, 2));
  } else if (! is_void(value)) {
    throw, "invalid value for FITS keyword " + kwd;
  }
  return value;
}
errs2caller, mira_get_fits_string;

func mira_get_fits_logical(fh, kwd) /* DOCUMENTED */
{
  value = fits_get(fh, kwd);
  if (is_scalar(value) && structof(value) == char &&
      (value == _MIRA_FITS_TRUE || value == _MIRA_FITS_FALSE)) {
    value = (value == _MIRA_FITS_TRUE);
  } else if (! is_void(value)) {
    throw, "invalid value for FITS keyword " + kwd;
  }
  return value;
}
_MIRA_FITS_TRUE = 'T';
_MIRA_FITS_FALSE = 'F';
errs2caller, mira_get_fits_logical;

func mira_get_fits_integer(fh, kwd) /* DOCUMENTED */
{
  value = fits_get(fh, kwd);
  if (is_scalar(value) && identof(value) <= Y_LONG) {
    value = long(value);
  } else if (! is_void(value)) {
    throw, "invalid value for FITS keyword " + kwd;
  }
  return value;
}
errs2caller, mira_get_fits_integer;

func mira_get_fits_real(fh, kwd) /* DOCUMENTED */
{
  value = fits_get(fh, kwd);
  if (is_scalar(value) && identof(value) <= Y_DOUBLE) {
    value = double(value);
  } else if (! is_void(value)) {
    throw, "invalid value for FITS keyword " + kwd;
  }
  return value;
}
errs2caller, mira_get_fits_real;

func mira_find_fits_hdu(fh, xtension, start=, extname=, hduname=)
/* DOCUMENT mira_find_fits_hdu, fh, xtension;
         or bool = mira_find_fits_hdu(fh, xtension);

     moves to the first FITS Header Data Units (HDU) of FITS handle FH which
     has extension XTENSION and which matches all criteria set by the
     following keywords.

     Keyword START specifies from which FITS HDU the search should start.
     Default is to use current HDU.

     If keyword EXTNAME is specified, the HDU must be in a FITS extension
     whose name is EXTNAME.

     If keyword HDUNAME is specified, the header of the HDU must have a FITS
     card "HDUNAME" matching value HDUNAME.

     When called as a subroutine, an error is thrown (in the caller's context)
     if no HDU matching the criteria is found.  When called as a function, a
     boolean value indicating whether the HDU has been found is returned.

     If no matching HDU is found, the FITS current HDU remains unchanged;
     otherwise, the FITS current HDU is the one matchng all the criteria.


   SEE ALSO: fits_eof, fits_next_hdu.
 */
{
  if (is_scalar(xtension) && is_string(xtension)) {
    xtension = strcase(1n, strtrim(xtension, 2));
  } else {
    error, "XTENSION must be a string";
  }
  if (! is_void(hduname)) {
    if (is_scalar(hduname) && is_string(hduname)) {
      hduname = strcase(1n, strtrim(hduname, 2));
    } else {
      error, "HDUNAME must be a string";
    }
  }
  if (! is_void(extname)) {
    if (is_scalar(extname) && is_string(extname)) {
      extname = strcase(1n, strtrim(extname, 2));
    } else {
      error, "EXTNAME must be a string";
    }
  }
  initial_hdu = fits_current_hdu(fh);
  current_hdu = initial_hdu;
  if (! is_void(start)) {
    /* Move to specified HDU number. */
    if (! is_scalar(start) || ! is_integer(start) || start < 1) {
      error, "START must be an integer greater or equal 1";
    }
    if (start > current_hdu) {
      fits_rewind, fh;
      current_hdu = fits_current_hdu(fh);
    }
    while (start < current_hdu) {
      fits_next_hdu, fh;
      if (fits_eof(fh)) {
        error, swrite(format="START=%d is too large", start);
      }
      current_hdu = fits_current_hdu(fh);
    }
  }
  while (! fits_eof(fh)) {
    flag = (fits_get_xtension(fh) == xtension);
    if (flag && ! is_void(extname) && current_hdu > 1) {
      value = fits_get(fh, "EXTNAME");
      flag = (is_string(value) && strcase(1, strtrim(value, 2)) == extname);
    }

    if (flag && ! is_void(hduname)) {
      value = fits_get(fh, "HDUNAME");
      flag = (is_string(value) && strcase(1, strtrim(value, 2)) == hduname);
    }
    if (flag) {
      /* All criteria satisfied. */
      return 1n;
    }
    fits_next_hdu, fh;
    current_hdu = fits_current_hdu(fh);
  }

  /* HDU not found. */
  fits_goto_hdu, fh, initial_hdu;
  if (! am_subroutine()) {
    return 0n;
  }
  msg = "no " + xtension;
  if (! is_void(extname)) {
    msg += ("with EXTNAME='" + extname + "'");
  }
  if (! is_void(hduname)) {
    msg += ((is_void(extname) ? " with " : " and ") +
            "HDUNAME='" + hduname + "'");
  }
  msg += " found in FITS file";
  if (! is_void(hdu)) {
    msg += swrite(format=" after HDU=%d", hdu);
  }
  error, msg;
}
errs2caller, mira_find_fits_hdu;

/*---------------------------------------------------------------------------*/
/* READ/WRITE ALGORITHM PARAMATERS */

func mira_read_input_params(src, plugin=)
/* DOCUMENT tab = mira_read_input_params(src, plugin=...);

      yields a hash table with the input parameters found in FITS extension
      "IMAGE-OI INPUT PARAM".  Argument SRC can be a FITS file name or a FITS
      handle.  The returned table may be empty if no such extension is found or
      if it does not contain any known parameters.

      If keyword PLUGIN is specified, the given plugin is invoked to read its
      specific keywords.

   SEE ALSO: mira_write_input_params.
 */
{
  local fh;
  if (is_string(src)) {
    fh = fits_open(src, "r");
    finally_close = 1n;
  } else {
    eq_nocopy, fh, src;
    finally_close = 0n;
  }
  tab = h_new();
  if (mira_find_fits_hdu(fh, "BINTABLE", start=2,
                         extname="IMAGE-OI INPUT PARAM")) {
    h_set, tab,
      target           = mira_get_fits_string(   fh, "TARGET"),
      wavemin          = mira_get_fits_real(     fh, "WAVE_MIN"),
      wavemax          = mira_get_fits_real(     fh, "WAVE_MAX"),
      use_vis          = mira_get_fits_use_polar(fh, "USE_VIS"),
      use_vis2         = mira_get_fits_logical(  fh, "USE_VIS2"),
      use_t3           = mira_get_fits_use_polar(fh, "USE_T3"),
      initialhdu       = mira_get_fits_string(   fh, "INIT_IMG"),
      threshold        = mira_get_fits_real(     fh, "SOFT_CUT"),
      recenter         = mira_get_fits_logical(  fh, "RECENTER"),
      bootstrap        = mira_get_fits_integer(  fh, "REPEAT"),
      regul            = mira_get_fits_string(   fh, "RGL_NAME", lower=1),
      mu               = mira_get_fits_real(     fh, "RGL_WGT"),
      tau              = mira_get_fits_real(     fh, "RGL_TAU"),
      gamma            = mira_get_fits_real(     fh, "RGL_GAMM"),
      flux             = mira_get_fits_real(     fh, "FLUX"),
      fluxerr          = mira_get_fits_real(     fh, "FLUXERR"),
      "min",             mira_get_fits_real(     fh, "PXL_MIN"),
      "max",             mira_get_fits_real(     fh, "PXL_MAX"),
      xform            = mira_get_fits_string(   fh, "XFORM"),
      smearingfunction = mira_get_fits_string(   fh, "SMEAR_FN"),
      smearingfactor   = mira_get_fits_real(     fh, "SMEAR_FC"),
      maxeval          = mira_get_fits_integer(  fh, "MAXEVAL"),
      maxiter          = mira_get_fits_integer(  fh, "MAXITER"),
      mem              = mira_get_fits_integer(  fh, "OPT_MEM"),
      ftol             = mira_get_fits_real(     fh, "OPT_FTOL"),
      gtol             = mira_get_fits_real(     fh, "OPT_GTOL"),
      sftol            = mira_get_fits_real(     fh, "LNS_FTOL"),
      sgtol            = mira_get_fits_real(     fh, "LNS_GTOL"),
      sxtol            = mira_get_fits_real(     fh, "LNS_XTOL"),
      plugin           = mira_get_fits_string(   fh, "PLUGIN", lower=1);

    if (! is_void(tab.use_vis2)) {
      h_set, tab, use_vis2 = (tab.use_vis2 ? "all" : "none");
    }
    if (! is_void(tab.gamma)) {
      /* Fix units (SI for MiRA, degrees in the FITS file). */
      h_set, tab, gamma = tab.gamma*MIRA_DEGREE;
    }
    if (! is_void(tab.xform)) {
      h_set, tab, xform = strcase(0n, strtrim(tab.xform, 3));
    }
    if (! is_void(tab.smearingfunction)) {
      h_set, tab,
        smearingfunction = strcase(0n, strtrim(tab.smearingfunction, 3));
    }
    if (! is_void(tab.initial)) {
      h_set, tab, initial = mira_read_image(fh, hduname=tab.initial);
    }

    if (! is_void(tab.plugin)) {
      plugin_obj = _mira_load_plugin(tab.plugin);
    } else {
      /* support overriding plugin */
      plugin_obj = plugin
    }

    /* Maybe read plug-in specific keywords. */
    if (is_hash(plugin_obj)) {
      inform, "Reading plugin keywords";
      subroutine = plugin_obj.__vops__.read_keywords;
      subroutine, tab, fh;
    }
  }
  if (finally_close) {
    fits_close, fh;
  }
  return tab;
}

func mira_write_input_params(dest, master, opt)
/* DOCUMENT mira_write_input_params, dest, master, opt;

     writes a FITS binary table with the input parameters used for an image
     reconstruction with MiRA algorithm.  The extension name of the table is
     "IMAGE-OI INPUT PARAM" and it contains all parameters for running the
     algorithm.  Argument DEST is the destination which can be a FITS handle
     or a FITS file name.  Arguments MASTER and OPT are the main MiRA instance
     and the command line options.

   SEE ALSO: mira_read_input_params, mira_write_ouput_params.
 */
{
  /* Get FITS handle and start a new binary table. */
  local fh;
  if (is_string(dest)) {
    fh = fits_open(dest, "w");
    finally_close = 1n;
  } else {
    eq_nocopy, fh, dest;
    finally_close = 0n;
  }
  fits_new_bintable, fh;
  fits_set, fh, "EXTNAME", "IMAGE-OI INPUT PARAM";

  /* INIT_IMG keyword. */
  if (! is_void(opt.init_img)) {
    fits_set, fh, "INIT_IMG", opt.init_img,
      "HDUNAME of initial image";
  }
  if (! is_void(opt.threshold)) {
    fits_set, fh, "SOFT_CUT", opt.threshold,
      "Soft threshold of starting images";
  }
  if (! is_void(opt.recenter)) {
    fits_set, fh, "RECENTER", (opt.recenter ? 'T' : 'F'),
      "Recenter starting images";
  }
  if (! is_void(opt.bootstrap)) {
    fits_set, fh, "REPEAT", opt.bootstrap,
      "Number of algorithm repetitions";
  }

  /* XFORM keyword. */
  fits_set, fh, "XFORM", mira_xform_name(master),
    "Image to complex visibility transform";
  fits_set, fh, "SMEAR_FN", mira_smearing_function(master, 1n),
    "Smearing function";
  fits_set, fh, "SMEAR_FC", mira_smearing_factor(master, 1n),
    "Smearing factor";

  /* WAVE_MIN and WAVE_MAX keywords. */
  if (! is_void(opt.wavemin)) {
    fits_set, fh, "WAVE_MIN", opt.wavemin/MIRA_METER,
      "Minimal wavelength [m]";
  }
  if (! is_void(opt.wavemax)) {
    fits_set, fh, "WAVE_MAX", opt.wavemax/MIRA_METER,
      "Maximal wavelength [m]";
  }

  /* USE_VIS keyword. */
  bits = (flags & (MIRA_FIT_VISAMP|MIRA_FIT_VISPHI));
  if (bits == 0) {
    value = "NONE";
  } else if (bits == MIRA_FIT_VISAMP) {
    value = "AMP";
  } else if (bits == MIRA_FIT_VISPHI) {
    value = "PHI";
  } else {
    value = "ALL";
  }
  fits_set, fh, "USE_VIS", value, "Complex visibility data to consider";

  /* USE_VIS2 keyword. */
  fits_set, fh, "USE_VIS2", ((flags & MIRA_FIT_VIS2) != 0 ? 'T' : 'F'),
    "Powerspectrum data to consider";

  /* USE_T3 keyword. */
  bits = (flags & (MIRA_FIT_T3AMP|MIRA_FIT_T3PHI));
  if (bits == 0) {
    value = "NONE";
  } else if (bits == MIRA_FIT_T3AMP) {
    value = "AMP";
  } else if (bits == MIRA_FIT_T3PHI) {
    value = "PHI";
  } else {
    value = "ALL";
  }
  fits_set, fh, "USE_T3", value, "Bispectrum data to consider";

  /* RGL_NAME keyword. */
  regul = opt.regul;
  mu = 0.0;
  fits_set, fh, "RGL_NAME", regul, "Regularization method";

  /* RGL_TAU keyword. */
  if (regul == "hyperbolic") {
    fits_set, fh, "RGL_TAU", opt.tau, "Spatial gradient threshold";
    mu = opt.mu;
  }

  /* RGL_GAMMA keyword. */
  if (regul == "compactness") {
    fits_set, fh, "RGL_GAMM", opt.gamma/MIRA_DEGREE, "Prior FWHM [deg]";
    mu = opt.mu;
  }

  /* RGL_WGT keyword. */
  if (mu > 0) {
    fits_set, fh, "RGL_WGT", mu, "Regularization weight";
  }

  /* FLUX and FLUXERR keywords. */
  flux = opt.flux;
  if (is_void(flux)) {
    flux = 1.0;
  }
  fits_set, fh, "FLUX", flux, "Assumed flux level";
  if (! is_void(opt.fluxerr)) {
    fits_set, fh, "FLUXERR", opt.fluxerr, "Error bar on the flux";
  }

  /* PXL_MIN and PXL_MAX keywords. */
  if (! is_void(opt("min"))) {
    fits_set, fh, "PXL_MIN", opt("min"), "Minimum allowed pixel value";
  }
  if (! is_void(opt("max"))) {
    fits_set, fh, "PXL_MAX", opt("max"), "Maximum allowed pixel value";
  }

  /* MAXEVAL and MAXITER keywords. */
  if (! is_void(opt.maxeval)) {
    fits_set, fh, "MAXEVAL", opt.maxeval, "Maximum number of evaluations";
  }
  if (! is_void(opt.maxiter)) {
    fits_set, fh, "MAXITER", opt.maxiter, "Maximum number of iteration";
  }

  /* Other optimizer keywords. */
  if (! is_void(opt.mem)) {
    fits_set, fh, "OPT_MEM", opt.mem, "Steps memorized by the optimizer";
  }
  if (! is_void(opt.ftol)) {
    fits_set, fh, "OPT_FTOL", opt.ftol, "Function tolerance for convergence";
  }
  if (! is_void(opt.gtol)) {
    fits_set, fh, "OPT_GTOL", opt.gtol, "Gradient tolerance for convergence";
  }
  if (! is_void(opt.sftol)) {
    fits_set, fh, "LNS_FTOL", opt.sgtol, "Line search function tolerance";
  }
  if (! is_void(opt.sgtol)) {
    fits_set, fh, "LNS_GTOL", opt.sgtol, "Line search gradient tolerance";
  }
  if (! is_void(opt.sxtol)) {
    fits_set, fh, "LNS_XTOL", opt.sxtol, "Line search parameter tolerance";
  }

  /* Maybe add plug-in specific keywords. */
  plugin_obj = mira_plugin(master);
  if (is_hash(plugin_obj)) {
    inform, "Adding plugin keywords";
    fits_set, fh, "PLUGIN",  opt.plugin,  "plugin name";
    subroutine = plugin_obj.__vops__.add_keywords;
    subroutine, master, fh;
  }
  /* Finish the header and, possibly, close the FITS file. */
  fits_write_header, fh;
  if (finally_close) {
    fits_close, fh;
  }
}

func mira_write_output_params(dest, misc)
/* DOCUMENT mira_write_output_params, dest, misc;

     writes a FITS binary table with the output parameters used for an image
     reconstruction with MiRA algorithm.  The extension name of the table is
     "IMAGE-OI OUTPUT PARAM" and it contains all output parameters of the
     algorithm.  Argument DEST is the destination which can be a FITS handle or
     a FITS file name.  Arguments MISC is a hash table with the output
     parameters.

   SEE ALSO: mira_write_input_params.
 */
{
  /* Get FITS handle and start a new binary table. */
  local fh;
  if (is_string(dest)) {
    fh = fits_open(dest, "w");
    finally_close = 1n;
  } else {
    eq_nocopy, fh, dest;
    finally_close = 0n;
  }
  fits_new_bintable, fh;
  fits_set, fh, "EXTNAME", "IMAGE-OI OUTPUT PARAM";

  /* PROCSOFT keyword. */
  procsoft = "MiRA";
  if (!is_void(MIRA_VERSION)) {
      procsoft += " v" + MIRA_VERSION;
  }
  fits_set, fh, "PROCSOFT", procsoft, "Data reduction software and version";

  /* LAST_IMG keyword. */
  if (! is_void(misc.last_img)) {
    fits_set, fh, "LAST_IMG", misc.last_img,
      "HDUNAME of final image";
  }

  /* CONVERGENCE keyword. */
  early_stopping = (mira_less_or_equal(misc.maxiter, misc.niters) ||
                    mira_less_or_equal(misc.maxeval, misc.nevals));
  fits_set, fh, "CONVERGENCE", (early_stopping ? 'F' : 'T'),
      "Algorithm has converged";

  /* NITER keyword. */
  if (! is_void(misc.niters)) {
    fits_set, fh, "NITER", misc.niters,
      "Total number of algorithm iterations";
  }

  /* NEVAL keyword. */
  if (! is_void(misc.nevals)) {
    fits_set, fh, "NEVAL", misc.nevals,
      "Total number of objective function evaluations";
  }

  /* NDATA keyword. */
  if (! is_void(misc.ndata)) {
    fits_set, fh, "NDATA", misc.ndata,
      "Total number of independent measurements";
  }

  /* CHISQ keyword. */
  if (! is_void(misc.fdata)) {
    fits_set, fh, "CHISQ", misc.fdata, "Reduced chi-squared";
  }

  /* RGL_WGT keyword. */
  if (! is_void(misc.mu)) {
    fits_set, fh, "RGL_WGT", misc.mu, "Regularization weight";
  }

  /* FPRIOR keyword. */
  if (! is_void(misc.fprior)) {
    fits_set, fh, "FPRIOR", misc.fprior, "Regularization penalty";
  }

  /* FLUX keyword. */
  if (! is_void(misc.fprior)) {
    fits_set, fh, "FLUX", misc.flux, "Total image flux";
  }

  /* Finish the header and, possibly, close the FITS file. */
  fits_write_header, fh;
  if (finally_close) {
    fits_close, fh;
  }
}
