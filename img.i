/*
 * img.i -
 *
 * Routines for dealing with images (i.e. 2D arrays) in Yorick.
 *
 * ----------------------------------------------------------------------------
 *
 * This file is part of YLib (Yorick Library) which is licensed under the MIT
 * "Expat" License:
 *
 * Copyright (C) 2000, 2014, Éric Thiébaut.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 * ----------------------------------------------------------------------------
 *
 * Routines:
 *   img_dims - returns dimension [WIDTH,HEIGHT] of an image
 *   img_plot - plot an image
 *   img_cbar - add a color to an image plot
 *   img_convolve - convolution/correlation of images
 *   img_crop - crop an image
 *   img_interpolate - bi-linear interpolation of an image
 *   img_extract_parallelogram_as_rectangle - (as its name says)
 *   img_max - get coordinates of maximum in a 2-D array
 *   img_fft_centered_at_max - recenter an image at its maximum according
 *           to fft indexing
 *   img_pad - expand an image
 *   img_paste - copy an image into another one at a given location
 *   img_photometry - measure integrated intensity inside a circular region
 *   img_flt_max - filter an image
 *   img_flt_flac - flipped-local-auto-convolution of an image
 *   img_get_format - get image file format
 *   img_read, img_write - read/write PNM/JPEG/PNG/TIFF/FITS/GIF image
 *   img_tmpnam - get name of temporary file
 *   img_to_gray - convert an image to grayscale
 *   img_resize - resize an image by moving average and/or linear interpolation
 *   img_thumbnail - resize an image by moving average
 *
 * ----------------------------------------------------------------------------
 */

_IMG_VERSION = "$Date: 2014-10-21 09:07:31 +0200$";

require, "utils.i";

func img_dims(img)
/* DOCUMENT img_dims(img)
     Returns dimensions of image IMG: [WIDTH,HEIGHT].  IMG must be a 2
     dimensional array.

   SEE ALSO dimsof. */
{
  if (! is_array(img) || (dims = dimsof(img))(1) != 2)
    error, "expecting 2D array";
  return dims(2:3);
}

/*---------------------------------------------------------------------------*/

func img_plot(img, first=, scale=, cmin=, cmax=, top=, cbar=,
              levels=, labels=, nticks=, vert=, adjust=, vport=,
              format=, color=, width=, height=, ticklen=, font=)
/* DOCUMENT img_plot, img;
     Plot image IMG using pli.  Keyword FIRST can be used to set
     coordinates of center of first (lower left) pixel, default is
     FIRST=1.0 (i.e. same coordinates as Yorick's indexing rules).  Keyword
     SCALE can be used to set the step size between adjacent pixels
     (default is SCALE=1.0).  Keywords FIRST and/or SCALE may have one or
     two values depending whether or not both axis have same value.

     If Keyword CBAR is true, then a color bar is drawn with img_cbar.

   SEE ALSO pli, img_dims. */
{
  if (is_array(img)) {
    dims = dimsof(img);
    ndims = dims(1);
  }
  if (ndims == 2) {
    dims = dims(2:3);
  } else if (ndims == 3 && structof(img) == char) {
    dims = dims(3:4);
  } else {
    error, "expecting 2-D array or RGB image";
  }
  if (is_void(scale)) scale = 1.0;
  if (is_void(first)) first = 1.0;
  ll = first - scale/2.0; // lower-left corner coordinates
  ur = ll + scale*dims;   // upper-right corner coordinates
  pli, img, ll(1), ll(0), ur(1), ur(0), top=top, cmin=cmin, cmax=cmax;
  if (cbar) {
    pl_cbar, img, cmin=cmin, cmax=cmax,
      levels=levels, labels=labels, nlabs=nticks, vert=vert,
      adjust=adjust, vport=vport, format=format, color=color,
      width=width, height=height, ticklen=ticklen, font=font;
  }
}

func img_cbar(z, cmin=, cmax=, vert=, vport=, adjust=,
              nticks=, levels=, labels=,
              color=, font=, height=, opaque=, orient=,
              width=, ticklen=, thickness=, format=)
{
  write, format="WARNING - %s\n",
    "img_cbar is deprecated, use pl_cbar instead";
  pl_cbar, z, cmin=cmin, cmax=cmax, vert=vert, vport=vport, adjust=adjust,
    nlabs=nticks, levels=levels, labels=labels,
    color=color, font=font, height=height, opaque=opaque, orient=orient,
    width=width, ticklen=ticklen, thickness=thickness, format=format;
}

/*---------------------------------------------------------------------------*/
/* RESIZING OF AN IMAGE */

func img_resize(z, dim1, dim2, dbl=)
/* DOCUMENT img_resize(z, s);
         or img_resize(z, dim);
         or img_resize(z, width, height);

      This function returns image Z resized by scale factor S, or resized so
      that its maximum dimension is DIM while preservering the aspect ratio,
      or resized to size WIDTH by DIM pixels.  Z can be a grayscale image
      (i.e., an integer or real 2D array) or a RGB image (i.e., a 3-by-M-by-N
      array of type char).

      For each dimension of the image (not the color axis if any), img_resize
      performs an interpolated moving average to subsample the image or a
      linear interpolation to magnify the image.  Thus img_resize can be used
      to produce accurate thumbnail images.

      When rescaling while preserving the aspect ratio, there may be a slight
      attenuation at the edges of the result.

      Keyword DBL can be set true to keep the result in double precision.  The
      default is to round the result back to the type of the input image.

   SEE ALSO: interp, integ, img_thumbnail.
 */
{
  /* Get image size and type. */
  id = identof(z);
  dims = dimsof(z);
  if (id == Y_CHAR && dims(1) == 3 && dims(2) == 3) {
    rgb = 1n;
    width = dims(3);
    height = dims(4);
  } else if (id <= Y_DOUBLE && dims(1) == 2) {
    rgb = 0n;
    width = dims(2);
    height = dims(3);
  } else {
    error, "bad image";
  }

  /* Compute output dimensions and normalized coordinates in the input and
     output images. */
  local x, xp, y, yp;
  q = 1.0; // flux scaling factor
  if (is_void(dim2)) {
    /* Rescale keeping aspect ratio. */
    if (is_integer(dim1) && is_scalar(dim1) && dim1 > 0) {
      scale = double(dim1)/max(width, height);
    } else if (is_real(dim1) && is_scalar(dim1) && dim1 > 0.0) {
      scale = double(dim1);
    } else {
      error, "bad dimension/scale";
    }
    dim1 = lround(scale*width);
    dim2 = lround(scale*height);
    if (dim1 < width) {
      /* Normalized coordinates after integration by "cum" operator. */
      x = indgen(0:width)*(1.0/width);
      xp = (indgen(0:dim1) - dim1/2.0)*(1.0/(scale*width)) + 0.5;
      q *= scale;
    } else if (dim1 > width) {
      /* Normalized coordinates for interpolation. */
      x = (indgen(width) - (width + 1)/2.0)*(1.0/width) + 0.5;
      xp = (indgen(dim1) - (dim1 + 1)/2.0)*(1.0/(scale*width)) + 0.5;
    }
    if (dim2 < height) {
      /* Normalized coordinates after integration by "cum" operator. */
      y = indgen(0:height)*(1.0/height);
      yp = (indgen(0:dim2) - dim2/2.0)*(1.0/(scale*height)) + 0.5;
      q *= scale;
    } else if (dim2 > height) {
      /* Normalized coordinates for interpolation. */
      y = (indgen(height) - (height + 1)/2.0)*(1.0/height) + 0.5;
      yp = (indgen(dim2) - (dim2 + 1)/2.0)*(1.0/(scale*height)) + 0.5;
    }
  } else {
    /* Rescale to given size. */
    if (! (is_integer(dim1) && is_scalar(dim1) && dim1 > 0) ||
        ! (is_integer(dim2) && is_scalar(dim2) && dim2 > 0)) {
      error, "bad dimensions";
    }
    if (dim1 < width) {
      /* Normalized coordinates after integration by "cum" operator. */
      x = indgen(0:width)*(1.0/width);
      xp = indgen(0:dim1)*(1.0/dim1);
      q *= (double(dim1)/width);
    } else if (dim1 > width) {
      /* Normalized coordinates for interpolation. */
      x = (indgen(width) - (width + 1)/2.0)*(1.0/width) + 0.5;
      xp = (indgen(dim1) - (dim1 + 1)/2.0)*(1.0/dim1) + 0.5;
    }
    if (dim2 < height) {
      /* Normalized coordinates after integration by "cum" operator. */
      y = indgen(0:height)*(1.0/height);
      yp = indgen(0:dim2)*(1.0/dim2);
      q *= (double(dim2)/height);
    } else if (dim2 > height) {
      /* Normalized coordinates for interpolation. */
      y = (indgen(height) - (height + 1)/2.0)*(1.0/height) + 0.5;
      yp = (indgen(dim2) - (dim2 + 1)/2.0)*(1.0/dim2) + 0.5;
    }
  }

  /* Resize the image. */
  if (q != 1.0) {
    z *= q; // also make Z of type double
  } else if (id != Y_DOUBLE) {
    z = double(z);
  }
  if (rgb) {
    /* Resize RGB image along second dimension, then along third dimension. */
    if (dim1 > width) {
      z = interp(z, x, xp, 2);
    } else if (dim1 < width) {
      z = interp(z(,cum,), x, xp, 2)(,dif,);
    }
    if (dim2 > height) {
      z = interp(z, y, yp, 3);
    } else if (dim2 < height) {
      z = interp(z(,,cum), y, yp, 3)(,,dif);
    }
  } else {
    /* Resize grayscale image along first dimension, then along second
       dimension. */
    if (dim1 > width) {
      z = interp(z, x, xp, 1);
    } else if (dim1 < width) {
      z = interp(z(cum,), x, xp, 1)(dif,);
    }
    if (dim2 > height) {
      z = interp(z, y, yp, 2);
    } else if (dim2 < height) {
      z = interp(z(,cum), y, yp, 2)(,dif);
    }
  }
  if (dbl || id == Y_DOUBLE) {
    return z;
  }
  if (id == Y_FLOAT) {
    return float(z);
  }
  z = lround(z);
  if (id == Y_CHAR) {
    return char(z);
  }
  if (id == Y_INT) {
    return int(z);
  }
  if (id == Y_SHORT) {
    return short(z);
  }
  return z;
}

func img_thumbnail(z, dim1, dim2, dbl=)
/* DOCUMENT img_thumbnail(z, s);
         or img_thumbnail(z, dim);
         or img_thumbnail(z, width, height);

      This function has been superceeded by img_resize (which to see) as it
      always uses an nterpolated moving average which is less efficient than a
      simple interpolation for magnifying an image.

      This function returns image Z resized by scale factor S, or resized so
      that its maximum dimension is DIM while preservering the aspect ratio,
      or resized to size WIDTH by DIM pixels.  Z can be a grayscale image
      (i.e., an integer or real 2D array) or a RGB image (i.e., a 3-by-M-by-N
      array of type char).

      This function apply a moving averaging to subsample the image, the
      operation is similar to interpolation when the image is magnified. This
      function can thus be used to magnify an image but is it is more
      efficient to use a bi-linear interpolation for that.  See for instance,
      img_resize.

      When rescaling while preserving the aspect ratio, there may be a slight
      attenuation at the edges.

      Keyword DBL can be set true to keep the result in double precision.  The
      default is to round the result back to the type of the input image.

   SEE ALSO: interp, integ, img_resize.
 */
{
  /* Get image size and type. */
  id = identof(z);
  dims = dimsof(z);
  if (id == Y_CHAR && dims(1) == 3 && dims(2) == 3) {
    rgb = 1n;
    width = dims(3);
    height = dims(4);
  } else if (id <= Y_DOUBLE && dims(1) == 2) {
    rgb = 0n;
    width = dims(2);
    height = dims(3);
  } else {
    error, "bad image";
  }

  /* Normalized coordinates after integration by "cum" operator. */
  x = (1.0/width)*indgen(0:width);
  y = (1.0/height)*indgen(0:height);

  /* Normalized coordinates in the thumbnail. */
  if (is_void(dim2)) {
    /* Rescale keeping aspect ratio. */
    if (is_integer(dim1) && is_scalar(dim1) && dim1 > 0) {
      scale = double(dim1)/max(width, height);
    } else if (is_real(dim1) && is_scalar(dim1) && dim1 > 0.0) {
      scale = double(dim1);
    } else {
      error, "bad dimension/scale";
    }
    dim1 = lround(scale*width);
    dim2 = lround(scale*height);
    xp = (indgen(0:dim1) - dim1/2.0)*(1.0/(scale*width)) + 0.5;
    yp = (indgen(0:dim2) - dim2/2.0)*(1.0/(scale*height)) + 0.5;
    scale = scale^2; // flux scale
  } else {
    /* Rescale to given size. */
    if (! (is_integer(dim1) && is_scalar(dim1) && dim1 > 0) ||
        ! (is_integer(dim2) && is_scalar(dim2) && dim2 > 0)) {
      error, "bad dimensions";
    }
    xp = indgen(0:dim1)*(1.0/dim1);
    yp = indgen(0:dim2)*(1.0/dim2);
    scale = (double(dim1)/width)*(double(dim2)/height); // flux scale
  }
  z *= scale; // also make Z of type double
  if (rgb) {
    z = interp(interp(z(,cum,), x, xp, 2)(,dif,)(,,cum), y, yp, 3)(,,dif);
  } else {
    z = interp(interp(z(cum,), x, xp, 1)(dif,)(,cum), y, yp, 2)(,dif);
  }
  if (dbl || id == Y_DOUBLE) {
    return z;
  }
  if (id == Y_FLOAT) {
    return float(z);
  }
  z = lround(z);
  if (id == Y_CHAR) {
    return char(z);
  }
  if (id == Y_INT) {
    return int(z);
  }
  if (id == Y_SHORT) {
    return short(z);
  }
  return z;
}

/*---------------------------------------------------------------------------*/
/* BI-LINEAR INTERPOLATION */

func img_interpolate(z, x, y)
/* DOCUMENT  img_interpolate(img, x, y)
     Returns image IMG interpolated (by bi-linear interpolation) at
     positions (X,Y).  The coordinates X and Y must be conformable and the
     result has dimension list dimsof(X,Y).  Note that coordinates run like
     Yorick indices, for instance (1,1) is the location of the lower-left
     pixel in the image.

   SEE ALSO: interp. */
{
  if (! is_array(z) || (dims = dimsof(z))(1) != 2) {
    error, "expecting 2-D array";
  }
  if ((width = dims(2)) < 2 || (height = dims(3)) < 2) {
    error, "array to interpolate must have at least 2 elements in every dimension";
  }
  x -= (i = ceil(x));
  i = long(i);
  if (min(i) < 1) {
    k = where(i < 1);
    i(k) = 1;
    x(k) = 0.0;
  }
  if (max(i) >= width) {
    k = where(i >= width);
    i(k) = width - 1;
    x(k) = 1.0;
  }
  y -= (j = ceil(y));
  j = long(j);
  if (min(j) < 1) {
    k = where(j < 1);
    j(k) = 1;
    y(k) = 0.0;
  }
  if (max(j) >= height) {
    k = where(j >= height);
    j(k) = height - 1;
    y(k) = 1.0;
  }
  i += (j - 1)*width; /* univariate index */
  return ((1.0 - x)*((1.0 - y)*z(i) + y*z(i + width)) +
          x*((1.0 - y)*z(i + 1) + y*z(i + (1 + width))));
}

func img_extract_parallelogram_as_rectangle(img, x1, y1, x2, y2, x3, y3, w, h)
/* DOCUMENT img_extract_parallelogram_as_rectangle(img, x1, y1, x2, y2, x3, y3,
                                                   width, height)
      Returns a WIDTH-by-HEIGHT rectangle obtained by bi-linear
      interpolation of a parallelogram region from image IMG.  The
      parallelogram is specified by the coordinates of 3 of its corner:
      (X1,Y1) is the upper-left corner, (X2,Y2) is the lower-left corner
      and (X3,Y3) is the lower-right corner.  Note that coordinates run
      like Yorick indices: (1,1) is the location of the lower-left pixel in
      the image.

   SEE ALSO: img_interpolate, LUsolve. */
{
  a = LUsolve([[1,0,1,0,1,0],
               [1,0,1,0,w,0],
               [h,0,1,0,1,0],
               [0,1,0,1,0,1],
               [0,1,0,1,0,w],
               [0,h,0,1,0,1]],
              [x1, y1, x2, y2, x3, y3]);
  u = double(indgen(w));
  v = double(indgen(h)(-,));
  return img_interpolate(img,
                         a(1) + a(2)*u + a(3)*v,
                         a(4) + a(5)*u + a(6)*v);
}

/*---------------------------------------------------------------------------*/

func img_max(img)
/* DOCUMENT img_max(img)
     Returns coordinates of first image maximum.

   SEE ALSO: img_dims, img_flt_max. */
{
  width = img_dims(img)(1);
  i = img(*)(mxx) - 1;
  return [i%width + 1, i/width + 1];
}

func img_fft_centered_at_max(img)
/* DOCUMENT img_fft_centered_at_max(img)
     Returns image IMG rolled so that its maximum is at coordinates (1,1).

   SEE ALSO: img_dims, img_max, img_flt_max, roll. */
{
  type = structof(img);
  img = roll(img, 1 - img_max(img));
  return type == structof(img) ? img : type(img);
}

/*---------------------------------------------------------------------------*/
/* CONVOLUTION/CORRELATION OF IMAGES */

func img_convolve(a, b, unroll=, pad=, correl=)
/* DOCUMENT img_convolve(a, b)

     Convolve image B by image A using FFT's.  If Yeti FFTW plugin is
     loaded, FFTW is used; otherwise Yorick's FFT is used.

     If keyword UNROLL is true, the result is centered at pixel
     ((WIDTH + 1)/2, (HEIGHT + 1)/2) where WIDTH and HEIGHT are the
     dimensuon of the result (and integere division is applied); the
     default is to have the result centered at pixel (1,1) according
     to FFT conventions.

     If keyword CORREL is true, the correlation of B by A instead of
     the convolution is computed.

     If keyword PAD is true, then A and B are padded with zeroes to
     match good dimensions for the FFT.  As a special case, with PAD=2
     the padding is such that there is no aliasing in the result,
     i.e. dimensions of the result are such that:
       WIDTH  >= DIMSOF(A)(2) + DIMSOF(B)(2) - 1
       HEIGHT >= DIMSOF(A)(2) + DIMSOF(B)(2) - 1

   SEE ALSO: fft_best_dim, fft, fftw.
 */
{
  real = (structof(a) != complex && structof(b) != complex);
  if (pad) {
    if (! is_array(a) || (adims = dimsof(a))(1) != 2 ||
        ! is_array(b) || (bdims = dimsof(b))(1) != 2) {
      error, "expecting 2-D arrays";
    }
    if (! is_func(fft_best_dim)) {
      require, "fft_utils.i";
    }
    if (pad == 2) {
      width = fft_best_dim(adims(2) + bdims(2) - 1);
      height = fft_best_dim(adims(3) + bdims(3) - 1);
    } else {
      width = fft_best_dim(max(adims(2), bdims(2)));
      height = fft_best_dim(max(adims(3), bdims(3)));
    }
    cdims = [2, width, height];
    if (adims(2) != width || adims(3) != height) {
      temp = array((real?double:complex), cdims);
      temp(1:adims(2), 1:adims(3)) = a;
      a = temp;
    }
    if (bdims(2) != width || bdims(3) != height) {
      temp = array((real?double:complex), cdims);
      temp(1:bdims(2), 1:bdims(3)) = b;
      b = temp;
    }
    // FIXME: check unroll offsets
    if (unroll) {
      off0 = ((adims(2:) - 1)/2);
      off1 = ((cdims(2:) - 1)/2) - ((bdims(2:) - 1)/2);
      offset = (correl ? off1 + off0 : off1 - off0);
    }
  } else {
    cdims = dimsof(a, b);
    if (is_void(cdims)) {
      error, "non-conformable arrays";
    }
    if (unroll) {
      offset = ((cdims(2:) - 1)/2);
    }
  }

  if (is_func(fftw)) {
    fwd = fftw_plan(cdims, +1, real=real);
    a = fftw(a, fwd);
    b = fftw(b, fwd);
    p = fftw((correl?conj(a):a)*b, fftw_plan(cdims, -1, real=real));
    a = b = []; // possibly free some memory
  } else {
    ws = fft_setup(cdims);
    a = fft(a, +1, setup=ws);
    b = fft(b, +1, setup=ws);
    p = fft((correl?conj(a):a)*b, -1, setup=ws);
    a = b = []; // possibly free some memory
    if (real) p = double(p);
  }
#if 0
  p = (1.0/numberof(p))*(unroll?roll(p,offset):p);
  i = where(p == max(p))(1) - 1;
  write, format="max at offset: %d %d\n", i%width, i/width;
  return p;
#endif
  return (1.0/numberof(p))*(unroll?roll(p,offset):p);
}

/*---------------------------------------------------------------------------*/

func img_pad(img, .., bg=, just=)
/* DOCUMENT img_pad(img)
       -or- img_pad(img, dimlist)
       -or- img_pad(img, width, height)
     Pad an image to another size, one of:
       - the smallest square 2D array that contains the image
       - an 2D array with dimension list DIMLIST
       - a WIDTH-by-HEIGHT array

     The padding value can be specified with keyword BG (for "background"),
     the default is 0.

     The type of the result depends on the types of IMG and BG.

     The justification is set by keyword JUST:
       JUST =  0 or nil -> lower-left (the default)
               1        -> center
              -1        -> at corners to preserve FFT indexing

   SEE ALSO: img_paste, img_dims. */
{
  /* Get new dimension list. */
  local new1, new2;
  dims = img_dims(img);
  old1 = dims(1);
  old2 = dims(2);
  nargs = 0;
  while (more_args()) {
    arg = next_arg();
    ++nargs;
    if ((s = structof(arg)) != long && s != int && s != short && s != char)
      error, "invalid non-integer dimension list";
    if ((arg_dims = dimsof(arg))(1) == 0) {
      /* got a scalar argument */
      if (nargs == 1)      new1 = long(arg);
      else if (nargs == 2) new2 = long(arg);
      else error, "too many new dimensions";
    } else if (nargs==1 && arg_dims(1)==1 && arg_dims(2)==3 && arg(1)==2) {
      /* got a valid dimlist */
      new1 = arg(2);
      new2 = arg(3);
    } else error, "bad dimension list";
  }
  if (nargs==0) new1 = new2 = max(dims);
  else if (nargs == 1 && is_void(new2)) new2 = new1;
  if (new1<old1 || new2<old2) error, "cannot reduce size";

  /* produce new image */
  if (is_void(bg)) {
    new = array(structof(img), new1, new2);
  } else if (is_array(bg) && dimsof(bg)(1) == 0) {
    new = array(structof(img(1) + bg)(bg), new1, new2);
  } else error, "bad value for keyword BG";
  if (is_void(just) || just==0) {
    /* image will not be centered */
    new(1:old1, 1:old2) = img;
  } else if (just==1) {
    /* image will be centered */
    i1 = (new1-old1)/2;
    i2 = (new2-old2)/2;
    new(i1+1:old1+i1, i2+1:old2+i2) = img;
  } else if (just==-1) {
    /* preserve FFT indexing */
    h1 = old1/2;
    h2 = old2/2;
    if (h1 && h2) new(new1-h1+1:new1, new2-h2+1:new2) = img(1:h1, 1:h2);
    if (h1) new(new1-h1+1:new1, 1:old2-h2) = img(1:h1, h2+1:old2);
    if (h2) new(1:old1-h1, new2-h2+1:new2) = img(h1+1:old1, 1:h2);
    new(1:old1-h1, 1:old2-h2) = img(h1+1:old1, h2+1:old2);
  } else error, "bad value for keyword JUST";

  return new;
}

/*---------------------------------------------------------------------------*/

func img_paste(dst, x, y, src)
/* DOCUMENT img_paste(dst, x, y, src)
     Paste image (2D array) SRC  at location (X,Y) in DST.  The coordinates
     (X,Y) are the  integer indices of the location in  DST where the lower
     left pixel  of SRC will be pasted.   SRC is possibly clipped  so as to
     fit  into DST  (e.g.  X  and/or Y  may be  less or  equal  zero).  The
     operation is  done "in-place"  and the result  is returned.   The data
     type of DST is unchanged.

  SEE ALSO: img_crop, img_pad, img_dims. */
{
  dst_dims = img_dims(dst);
  src_dims = img_dims(src);
  if (! is_integer_scalar(x)) error, "X must be an integer scalar";
  if (! is_integer_scalar(y)) error, "Y must be an integer scalar";
  x1 = (x0 = long(x)) + src_dims(1) - 1;
  y1 = (y0 = long(y)) + src_dims(2) - 1;

  /* insure that SRC will not be completely outside DST */
  if (x0<=dst_dims(1) && x1>=1 && y0<=dst_dims(2) && y1>=1) {
      /* maybe clip SRC along 1st dimension */
      if ((clip= x0<1)) {
        s0 = 2 - x0;
        x0 = 1;
      } else s0 = 1;
      if (x1>dst_dims(1)) {
        x1 = dst_dims(1);
        clip = 1n;
      }
      if (clip) src = src(s0:s0+x1-x0,);

      /* maybe clip SRC along 2nd dimension */
      if ((clip = y0<1)) {
        s0 = 2 - y0;
        y0 = 1;
      } else s0 = 1;
      if (y1>dst_dims(2)) {
        y1 = dst_dims(2);
        clip = 1n;
      }
      if (clip) src = src(,s0:s0+y1-y0);

      /* paste SRC in-place */
      dst(x0:x1, y0:y1) = src;
  }
  return dst;
}

func img_crop(z, bd=, bg=, roi=)
/* DOCUMENT img_crop(z)

     Removes borders from image Z that are the background color, and
     returns the result.  The result is a 2D array of same type as Z or nil
     if all values of Z are equal to the background and no borders (see
     keyword BD below) are requested.

     Keyword BG can be set to specify the background level or "min" to use
     the minimum value in Z, or "max" to use the maximum value in Z;
     otherwise the background level is guessed heuristically.  Note that
     "min" and "max" can be replaced by the min/max Yorick range keywords.

     Keyword BD can be set to specify the border (in pixels) to add around the
     image after cropping.  If not set, BD defaults to 0 (no border).

     If keyword ROI is set, the region of interest as [XMIN,XMAX,YMIN,YMAX]
     is returned where the bounds are 1-based (as for Yorick indices) and
     account for the border (if any).  Beware that if keyword BD is
     specified, the ROI may have out of range bounds (i.e. less or equal
     zero or beyond the size of the image).

   SEE ALSO: img_pad, img_paste.
 */
{
  if (! is_array(z) || (dims = dimsof(z))(1) != 2) {
    error, "expecting a 2-D image";
  }
  width = dims(2);
  height = dims(3);
  type = structof(z);
  if (is_void(bg)) {
    bg = z(1,1);
    if (z(0,1) != bg || z(1,0) != bg || z(0,0) != bg) {
      return z;
    }
  } else if (bg == (dummy = max) || bg == "max") {
    bg = max(z);
  } else if (bg == (dummy = min) || bg == "min") {
    bg = min(z);
  } else {
    bg = type(bg);
  }
  if (is_void(bd)) bd = 0;
  i = where(z != bg);
  if (is_array(i)) {
    --i;
    x = i%width;
    y = i/width;
    xmin = 1 + min(x);
    xmax = 1 + max(x);
    ymin = 1 + min(y);
    ymax = 1 + max(y);
    if (roi) {
      return [xmin - bd, xmax + bd, ymin - bd, ymax + bd];
    }
    if (xmin - bd >= 1 && xmax + bd <= width &&
        ymin - bd >= 1 && ymax + bd <= height) {
      return z(xmin - bd : xmax + bd, ymin - bd : ymax + bd);
    }
    nx = xmax - xmin + 1;
    ny = ymax - ymin + 1;
    c = array(bg, 2*bd + nx, 2*bd + ny);
    c(1 + bd : nx + bd, 1 + bd : ny + bd) = z(xmin:xmax, ymin:ymax);
    return c;
  }
  if (roi) {
    return [-bd, bd, -bd, bd];
  }
  if (bd >= 1) return array(bg, 2*bd, 2*bd);
}

func img_to_gray(img)
/* DOCUMENT img_to_gray(img)
     Convert image IMG to grayscale.  IMG can be a grayscale image or
     an RGB image (3-by-WIDTH-by-HEIGHT array of type char).

   SEE ALSO: bytscl.
 */
{
  if (is_array(img)) {
    if ((rank = (dims = dimsof(img))(1)) == 2) return img;
    if ((structof(img) == char) && (rank == 3) && (dims(2) == 3)) {
      type = structof(img);
      return bytscl(0.33*img(1,..) + 0.5*img(2,..) + 0.17*img(3,..),
                    top=255, cmin=0, cmax=255);
    }
  }
  error, "bad image type or dimensions";
}

/*---------------------------------------------------------------------------*/

func img_flatten(a, n)
{
  if (! is_array(a)) error, "expecting array argument";
  dims = dimsof(a);
  ndims = dims(1);
  if (ndims <= 2) return a;
  width = dims(2);
  height = dims(3);
  depth = numberof(a)/(width*height);
  if (is_void(n)) {
    n = max(1, long(0.5 + sqrt(double(height*depth)/double(width))));
  }
  w = n*width;
  h = ((depth + n - 1)/n)*height;
  out = array(structof(a), w, h);
  for (k=0 ; k<depth ; ++k) {
    i = (k%n)*width;
    j = (k/n)*height;
    out(i+1:i+width, j+1:j+height) = a(,,k+1);
  }
  return out;

}

/*---------------------------------------------------------------------------*/
/* MEASUREMENT */

/*
 * Circle of radius R and center (X0,Y0) rounded to nearest pixel:
 *
 *         floor(sqrt((X - X0)^2 + (Y - Y0)^2) + 0.5) = R
 *   <==>  R <= sqrt((X - X0)^2 + (Y - Y0)^2) + 0.5 < R + 1
 *   <==>  R - 0.5 <= sqrt((X - X0)^2 + (Y - Y0)^2) < R + 0.5
 *
 * Disk of radius R and center (X0,Y0) rounded to nearest pixel:
 *
 *         floor(sqrt((X - X0)^2 + (Y - Y0)^2) + 0.5) <= R
 *   <==>  sqrt((X - X0)^2 + (Y - Y0)^2) + 0.5 < R + 1
 *   <==>  sqrt((X - X0)^2 + (Y - Y0)^2) <  R + 0.5
 *   <==>       (X - X0)^2 + (Y - Y0)^2  < (R + 0.5)^2
 *
 * Corresponding rectangular region of interest:
 *
 *  Y = Y0  ==>  abs(X - X0) < R + 0.5
 *              -(R + 0.5) < X - X0 < R + 0.5
 *              X0 - (R + 0.5) < X < X0 + (R + 0.5)
 *
 */

func img_photometry(img, x, y, r, bg=)
/* DOCUMENT img_photometry(img, x, y, r)
     Returns  sum of pixel  values in  image IMG  whithin circular  area of
     radius R and centered at (X,Y).  Arrays X, Y and R are in pixel units,
     they may have any geometry  but must be conformable.  Coordinates have
     the  same  origin  as  array  indices:  the first  pixel  in  IMG  has
     coordinates (1,1).

   SEE ALSO dimsof, avg. */
{
  /* Get dimensions of image and make X, Y and R arrays conformable. */
  if (! is_array(img) || (dims = dimsof(img))(1) != 2)
    error, "expecting 2D array";
  nx = dims(2);
  ny = dims(3);
  if (is_void((dims = dimsof(x, y, r)))) error, "X, Y and R not conformable";
  zero = array(double, dims);
  if (structof((x += zero)) != double ||
      structof((y += zero)) != double ||
      structof((r += zero)) != double) error, "bad data type";
  zero = [];

  rp = r + 0.5;
  i1 = max(1, long(ceil(x - rp)));
  i2 = min(nx, long(floor(x + rp)));
  j1 = max(1, long(ceil(y - rp)));
  j2 = min(ny, long(floor(y + rp)));
  rp *= rp;
  result = array((structof(img) == complex ? complex : double), dims);
  n = numberof(x);
  for (k=1 ; k<=n ; ++k) {
    if ((i1k = i1(k)) <= (i2k = i2(k)) &&
        (j1k = j1(k)) <= (j2k = j2(k))) {
      irng = i1k:i2k;
      jrng = j1k:j2k;
      u = indgen(irng) - x(k);
      v = indgen(jrng) - y(k);
      if (is_array((i = where(u*u + (v*v)(-,) < rp(k))))) {
        result(k) = avg(img(irng, jrng)(i));
      }
    }
  }
  if (! is_void(bg)) result -= bg;
  PI = 3.14159265358979323848;
  return (2*PI*r*r)*result;
}

/*---------------------------------------------------------------------------*/
/* FILTERS */

func img_flt_max(img, width, uniq=)
/* DOCUMENT img_flt_max(img, width)
     Return indices of pixels in image IMG that have the maximum value in a
     local  WIDTH-by-WIDTH box  (WIDTH must  be odd).   If keyword  UNIQ is
     true, a uniq maximum is  selected in every WIDTH-by-WIDTH box, in this
     case the intensity of the maxima will be in ascending order.

   SEE ALSO: img_flt_flac. */
{
  if (width%2 != 1) error, "width must be odd";
  w = (width-1)/2;
  if (w<=0) {
    if (w<0) error, "width must be non-negative";
    return img;
  }
  dims = img_dims(img);
  dim1 = dims(1);
  dim2 = dims(2);

  /* Filter image to replace pixel value by local maximum value. */
  if (w+1 >= dim1) {
    tmp = img(max,)(-:1:dim1,);
  } else {
    tmp = img;
    for (i=1 ; i<=dim1 ; ++i) tmp(i,) = img(max:max(i-w,1):min(i+w,dim1),);
  }
  if (w+1 >= dim2) {
    lmx = tmp(,max)(,-:1:dim2);
  } else {
    lmx = tmp;
    for (i=1 ; i<=dim2 ; ++i) lmx(,i) = tmp(,max:max(i-w,1):min(i+w,dim2));
  }
  tmp = [];

  /* Get indices of local maxima. */
  i = where(img >= lmx);
  if (! uniq || ! is_array(i)) return i;
  lmx = [];

  /* Sort maxima (because we already know that close maxima must have the same
     value) and remove multiples. */
  z = img(i);
  s = sort(z);
  i = i(s);
  x = 1 + (i-1)%dim1;
  y = 1 + (i-1)/dim1;
  z = z(s);
  if (numberof(z) > 1) {
    j = 1;
    k = where((z(dif)!=0) + (abs(x(dif)) > w) + (abs(y(dif)) > w));
    if (is_array(k)) grow, j, 1+k;
    i = i(j);
  }
  return i;
}

func img_flt_flac(img, width)
/* DOCUMENT img_flt_flac(img, width)
     Compute ``flipped-local-auto-convolution'' (sic!)   of image IMG.  The
     output  pixel  value  is  the   local  convolution  --  in  a  box  of
     WIDTH-by-WIDTH pixels  (WIDTH must  be odd) --  of the input  image by
     itself rotated  by 180 degrees.  This  is usefull to  locate spikes in
     IMG that may have different shapes but are all nearly symmetrical with
     respect to both axis (i.e.   each spike has its shape nearly unchanged
     by a  180 degrees  rotation).  This processing  is a kind  of adaptive
     filtering.   WIDTH should  be roughly  as large  as the  typical spike
     width.  The impact of input noise  in the result is reduced by using a
     larger WIDTH. But  WIDTH should remain smaller than  half the smallest
     separation between  spikes.  The  computation time is  proportional to
     WIDTH^2.

   SEE ALSO: img_flt_max. */
{
  if (width%2 != 1) error, "width must be odd";
  w = (width-1)/2;
  dims = img_dims(img);
  dim1 = dims(1);
  dim2 = dims(2);
  if (max(dim1,dim2)<width) return array(double, dim1, dim2);
  if (structof(img)==string || structof(img(1)+0.0)!=double)
    error, "bad data type for IMG";
  if (structof(img)!=double) img = double(img);
  if (w==0) return img*img;

  // zero offset
  s = img(w+1:-w, w+1:-w);
  s = s*s/2.0;

  // non-zero offsets (use symmetry for speed up by a factor of 2)
  for (w2=0 ; w2<=w ; ++w2) {
    for (w1=(w2>0?-w:1) ; w1<=w ; ++w1) {
      s1 = img(w+1+w1 :  w1-w, w+1+w2 :  w2-w);
      s2 = img(w+1-w1 : -w1-w, w+1-w2 : -w2-w);
      s += s1*s2;
    }
  }

  out = array(double, dim1, dim2);
  out(w+1:-w, w+1:-w) = 2.0*s;
  return out;
}

/*---------------------------------------------------------------------------*/
/* READING/WRITING IMAGES */

local IMG_NONE, IMG_PNM, IMG_JPEG, IMG_PNG, IMG_TIFF, IMG_FITS, IMG_GIF;
func img_get_file_type(filename, type=, reading=, noerr=)
/* DOCUMENT img_get_file_type(filename)
     Returns image type for file FILENAME.  If keyword READING is true,
     then FILENAME must exists and the image type is obtained from the file
     signature (the 4 first bytes in the file).  Otherwise, if keyword TYPE
     is true it must be one of the string: "pnm", "jpeg", "png", "tiff",
     "fits", or "gif". Finally if none of this keywords is set, the image
     type is guessed from the file extension.  Note that the function
     returns IMG_NONE (0) if the file is not readable or if the image type
     cannot be guessed.

     The returned value is one of:
       0 = IMG_NONE   unknown image type;
       1 = IMG_PNM    portable anymap (PBM/PGM/PPM) image;
       2 = IMG_JPEG   JPEG image;
       3 = IMG_PNG,   portable network graphic image;
       4 = IMG_TIFF   TIFF image;
       5 = IMG_FITS   FITS (flexible image transport system) file;
       6 = IMG_GIF    GIF image;

   SEE ALSO: img_read, img_write. */
{
  if (reading) {
    /* Guess file type from 4 byte header, according to magic numbers:
     *   \377 \330      = JPEG
     *   \211 "PNG"     = PNG
     *   "GIF8"         = GIF
     *   "MM" \000 \052 = TIFF image data, big-endian
     *   "II" \052 \000 = TIFF image data, little-endian
     *   "P1" space     = ascii PBM (portable bitmap)
     *   "P2" space     = ascii PGM (portable graymap)
     *   "P3" space     = ascii PPM (portable pixmap)
     *   "P4" space     = raw PBM
     *   "P5" space     = raw PGM
     *   "P6" space     = raw PPM
     */
    magic = array(char, 4);
    if (! (file = open(filename, "rb", 1))) {
      return IMG_NONE;
    }
    n = _read(file, 0, magic);
    close, file;
    z = '\0';
    c1 = (n >= 1 ? magic(1) : z);
    c2 = (n >= 2 ? magic(2) : z);
    c3 = (n >= 3 ? magic(3) : z);
    c4 = (n >= 4 ? magic(4) : z);
    if (c1 == '\377' &&  c2 == '\330')
      return IMG_JPEG;
    if (c1 == '\211' &&  c2 == 'P' && c3 == 'N' && c4 == 'G')
      return IMG_PNG;
    if ((c1 == 'M' &&  c2 == 'M' && c3 == '\000' && c4 == '\052') ||
        (c1 == 'I' &&  c2 == 'I' && c3 == '\052' && c4 == '\000'))
      return IMG_TIFF;
    if (c1 == 'S' &&  c2 == 'I' && c3 == 'M' && c4 == 'P')
      return IMG_FITS;
    if (c1 == 'P' && c2 >= '1' && c2 <= '6' &&
        (c3 == '\n' || c3 == '\r' || c3 == ' ' || c3 == '\t'))
      return IMG_PNM;
    if (c1 == 'G' &&  c2 == 'I' && c3 == 'F' && c4 == '8')
      return IMG_GIF;
    return IMG_NONE;
  }
  if (is_void(type)) {
    /* Guess file type from its extension. */
    if (strpart(filename, -3:-3) == ".") {
      ext = strcase(0, strpart(filename, -2:0));
      if (ext == "png") return IMG_PNG;
      if (ext == "jpg") return IMG_JPEG;
      if (ext == "tif") return IMG_TIFF;
      if (ext == "fit" || ext == "fts") return IMG_FITS;
      if (ext == "gif") return IMG_GIF;
      if (ext == "pnm" || ext == "pbm" ||
          ext == "pgm" || ext == "ppm") return IMG_PNM;
    }
    if (strpart(filename, -4:-4) == ".") {
      ext = strcase(0, strpart(filename, -3:0));
      if (ext == "jpeg") return IMG_JPEG;
      if (ext == "tiff") return IMG_TIFF;
      if (ext == "fits") return IMG_FITS;
    }
    return IMG_NONE;
  }
  if (is_scalar(type)) {
    if (is_string(type)) {
      if (type == "pnm") return IMG_PNM;
      if (type == "jpeg") return IMG_JPEG;
      if (type == "png") return IMG_PNG;
      if (type == "tiff") return IMG_TIFF;
      if (type == "fits") return IMG_FITS;
      if (type == "gif") return IMG_GIF;
      if (! noerr) error, "unknown image type name";
    } else if (is_integer(type)) {
      if ((1 <= type) && (type <= 6)) {
        return long(type);
      }
      if (! noerr) error, "unknown image type identifier";
    }
  }
  if (! noerr) error, "image type must be a name or an integer";
  return IMG_NONE;
}
IMG_NONE = 0;
IMG_PNM = 1;
IMG_JPEG = 2;
IMG_PNG = 3;
IMG_TIFF = 4;
IMG_FITS = 5;
IMG_GIF = 6;

func img_read(filename, tmp=, type=)
/* DOCUMENT img_read(filename)

     This function reads an image in FILENAME and returns it.  Supported
     image formats are: PNM (PBM/PBM/PPM), JPEG, PNG, TIFF, FITS and GIF.
     For some formats, a temporary PNM image needs to be created; the name
     of the temporary file can be specified with keyword TMP.  Whatever is
     the format, the rows of the image are ordered so that the image
     displays consistently in Yorick.

   SEE ALSO: system, pnm_read, fits_read, img_get_file_type,
             img_tmpnam, expand_file_name, protect_file_name,
             jpeg_read, png_read, tiff_read. */
{
  if (! open(filename, "r", 1)) error, "input file does not exists";
  type = img_get_file_type(filename, reading=1n, type=type);
  if (type == IMG_PNM) {
    convert = 0; /* no needs for conversion */
  } else if (type == IMG_JPEG) {
    if (is_func(jpeg_read)) return jpeg_read(filename)(.., ::-1);
    convert = "jpegtopnm -dct float";
  } else if (type == IMG_PNG) {
    if (is_func(png_read)) return png_read(filename)(.., ::-1);
    convert = "pngtopnm";
  } else if (type == IMG_TIFF) {
    if (is_func(tiff_read)) {
      if (strglob("~*", filename)) {
        /* Yeti-TIFF interface does not expand tilde, we must do that
           explicitely. */
        filename = expand_file_name(filename);
      }
      return tiff_read(filename);
    }
    convert = "tifftopnm";
  } else if (type == IMG_FITS) {
    /* FIXME: read compressed FITS file */
    if (! is_func(fits_read)) include, "fits.i", 1;
    return fits_read(filename);
  } else if (type == IMG_GIF) {
    convert = "giftopnm";
  } else {
    error, "unknown image type";
  }

  /* read image as PNM file, possibly after conversion */
  if (! is_func(pnm_read)) include, "pnm.i", 1;
  if (convert) {
    if (strglob("~*", filename)) {
      /* Expand the tilde, if any, for the shell script. */
      filename = expand_file_name(filename);
    }
    if (is_void(tmp)) tmp = img_tmpnam(filename);
    if (catch(-1)) {
      remove, tmp;
      error, catch_message;
    }
    system, swrite(format="%s %s 2>/dev/null >%s", convert,
                   protect_file_name(filename), protect_file_name(tmp));
    img = pnm_read(tmp);
    remove, tmp;
  } else {
    img = pnm_read(filename);
  }
  return img;
}

func img_write(img, filename, type=, cmin=, cmax=, tmp=,
               quality=, optimize=, grey=, gray=, progressive=, comment=,
               smooth=, eps=,
               compression=, interlace=, /* transparent=, gamma=, */
               bitpix=)
/* DOCUMENT img_write, img, filename;
     Writes image IMG into file FILENAME as PNM (PBM/PBM/PPM), JPEG, PNG,
     TIFF or FITS image.  Except for a FITS file, if pixel type of IMG is
     not 'char', the pixel values are scaled to unsigned bytes (0-255) with
     bytscl function (which see).  The image type can be specified by the
     keyword TYPE otherwise it is automatically guessed from FILENAME
     extension (see img_get_file_type).

   KEYWORDS
     Keywords common to all format:
       TYPE - Output image type, one of: "jpeg", "pnm", "png", "tiff"
           or "fits".
       TMP - Name of the temporary file (format PBM, PGM, or PPM, see
           pnm_write) to create.  Default is FILENAME~NUMBER where NUMBER is
           the smallest integer such that no file with the same name already
           exists (note that under race conditions the name of the default
           temporary file is not guaranteed to be unique).
       CMIN, CMAX - Optional lower/upper cutoff (see bytscl) -- not used
           for FITS format.
       GRAY/GREY - Creates grayscale image using NTCS standard.

     Keywords for JPEG images:
       EPS - If true, "jpeg2ps" is used to generate an encapsulated
           PostScript (level 2) image named FILENAME.eps (the JPEG image is
           not removed).
       QUALITY - JPEG quality (default 75).
       OPTIMIZE - Creates optimized JPEG image (on by default if pnmtojpeg
           is used).
       PROGRESSIVE - Creates a progressive JPEG file.
       COMMENT - Text comment.
       SMOOTH=0-100 - Smooth the input image to eliminate dithering
           noise, 0 (the default) means no smoothing.

     Keywords for PNG images:
       COMPRESSION=1-9 - Level of compression (default is 6).
       INTERLACE - Creates an interlaced PNG file.

     Keywords for FITS images:
       BITPIX=n - Bits-per-pixel value.


   SEE ALSO: pnm_write, bytscl, system, img_get_file_type,
             img_tmpnam, expand_file_name, protect_file_name. */
{
  filename = expand_file_name(filename);
  type = img_get_file_type(filename, type=type);

  /* Convert image to bytes (FIXME: PNG can handle 16-bit images but it
     does not work with Yorick-Z extension). */
  if (type == IMG_FITS) {
    rescale = 0n;
  } else {
    rescale = ((structof(img) != char)
               || (! is_void(cmin) && cmin != 0)
               || (! is_void(cmax) && cmax != 255));
  }
  if (gray || grey) {
    dims = dimsof(img);
    if (is_array(img) && (dims = dimsof(img))(1) == 3
        && dims(2) == 3 && structof(img) == char) {
      if (is_void(cmin)) cmin =   0;
      if (is_void(cmax)) cmax = 255;
      img = bytscl(0.33*img(1,..) + 0.50*img(2,..) + 0.17*img(3,..),
                    top=255, cmin=cmin, cmax=cmax);
      rescale = 0n;
    }
  }
  if (rescale) {
    img = bytscl(img, top=255, cmin=cmin, cmax=cmax);
  }

  if (type == IMG_PNM) {
    convert = 0; /* no needs for conversion */
  } else if (type == IMG_JPEG) {
    if (is_void(quality)) quality = 75;
    if (is_func(jpeg_write) && ! progressive && ! optimize
        && is_void(smooth)) {
      jpeg_write, filename, img(.., ::-1), comments, quality;
      if (eps) goto make_eps;
      return;
    }
    convert = swrite(format="pnmtojpeg --dct=float --quality=%d", quality);
    if (is_void(optimize) || optimize) convert += " --optimize";
    if (grey || gray) convert += " --grayscale";
    if (progressive) convert += " --progressive";
    if (! is_void(comment)) convert += swrite(format=" --comment=\"%s\"",
                                              comment);
    if (! is_void(smooth)) convert += swrite(format=" --smooth=%d", smooth);
  } else if (type == IMG_PNG) {
    if (is_func(png_write)) {
      // FIXME: use png_write options...
      png_write, filename, img(.., ::-1);
      return;
    }
    convert = "pnmtopng";
    if (! is_void(compression))
      convert += swrite(format=" -compression=%d", compression);
    if (interlace)
      convert += " -interlace";
  } else if (type == IMG_TIFF) {
    convert = "pnmtotiff";
  } else if (type == IMG_FITS) {
    /* TO DO: read compressed FITS file */
    if (! is_func(fits_write)) include, "fits.i", 1;
    fits_write, filename, img, overwrite=1, bitpix=bitpix;
    return;
  } else if (type == IMG_GIF) {
    error, "GIF format not supported for writing";
  }

  /* Write image as PNM file then convert to output type.  FIXME: the
     format is forced to be 8-bit PGM in part to avoid a bug in "pnm.i"
     with packed bits. */
  if (! is_func(pnm_write)) include, "pnm.i", 1;
  if (convert) {
    if (is_void(tmp)) tmp = img_tmpnam(filename);
    if (catch(-1)) {
      remove, tmp;
      error, catch_message;
    }
    pnm_write, img, tmp, bits=8, noscale=1;
    filename_p = protect_file_name(filename);
    system, swrite(format="%s %s 2>/dev/null >%s",
                   convert, protect_file_name(tmp), filename_p);
    remove, tmp;
  } else {
    pnm_write, img, filename, bits=8, noscale=1;
  }
  if (eps && type == IMG_JPEG) {
    /* make an encapsulated postscript file with JPEG image */
  make_eps:
    if (is_void(filename_p)) filename_p = protect_file_name(filename);
    system, swrite(format="jpeg2ps %s >%s.eps", filename_p, filename_p);
  }
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func img_tmpnam(name)
/* DOCUMENT img_tmpnam(name)
     Return a  string in the form: NAME~#  where # is an  integer chosen so
     that file  NAME~# does not exists.   Beware that there  is no absolute
     warranty that the returned name is not used elsewhere (for instance if
     two programs run at the same time and call the same function) but this
     is highly improbable.  In order to limit the  probabilty of such clash
     to occur, an empty file named NAME~# is created by the function.

   SEE ALSO: open. */
{
  fmt = "%s~%d";
  i = 0;
  do {
    tmp = swrite(format=fmt, name, ++i);
  } while (open(tmp, "r", 1));
  open, tmp, "w"; /* create empty file */
  return tmp;
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 79
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
