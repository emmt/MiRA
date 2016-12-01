/*
 * utils.i --
 *
 * General purpose utility routines for Yorick.
 *
 * ----------------------------------------------------------------------------
 *
 * This file is part of YLib (Yorick Library) which is licensed under the MIT
 * "Expat" License:
 *
 * Copyright (C) 1995-2014, Éric Thiébaut.
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
 *   basename     - get basename of file from path
 *   cast         - change the type and/or the dimensions of an array
 *   depth_of     - get 3rd dimension of an array
 *   dirname      - get directory name from path
 *   dump_text    - dump string array in a text file
 *   eval         - evaluate textual code
 *   expand_file_name - expand leading tilde in file name(s)
 *   filesize     - get size of a file in bytes
 *   get_file_name - get path name of file associated with a stream
 *   glob         - get list of files matching glob-style pattern
 *   guess_compression - guess compression method of existing file
 *   height_of    - get 2nd dimension of an array
 *   is_integer_scalar - check if argument is a scalar of integer type
 *   is_string_scalar - check if argument is a scalar string
 *   lambda       - create anonymous function
 *   load_text    - load all lines of a text file
 *   map          - apply function to every elements of array or list
 *   moments      - get first moments of a sampled distribution
 *   ncols_of     - get number of columns of an array
 *   ndims_of     - get number of dimension of an array
 *   nrows_of     - get number of rows of an array
 *   open_url     - open an URL into a web browser
 *   pdb_list     - lists contents of PDB binary file
 *   pdb_restore_all - restore all non-record variables of a PDB file
 *   protect_file_name - protect special characters in file name
 *   pw_get_gid   - get group numerical identifier from /etc/passwd file
 *   pw_get_home  - get home directory of user from /etc/passwd file
 *   pw_get_name  - get real user name from /etc/passwd file
 *   pw_get_shell - get path to shell of user from /etc/passwd file
 *   pw_get_uid   - get user numerical identifier from /etc/passwd file
 *   pw_get_user  - get user name from /etc/passwd file
 *   pwd          - print/get working directory
 *   raw_read     - read binary data from file
 *   read_ascii   - read ascii numeric data
 *   rescale      - interpolate a multi-dimensional array
 *   smooth       - smooth array along its dimensions
 *   spline_zoom  - resize an array by cubic spline interpolation
 *   stat         - display statistics/info about symbols/expressions
 *   strchr       - forward search for a character in a string
 *   strcut       - cut text to fit into lines of given length
 *   strjoin      - make array of strings into a single string
 *   strlower     - convert string(s) to lower case letters
 *   strrchr      - backward search for a character in a string
 *   strupper     - convert string(s) to upper case letters
 *   strip_file_extension - remove extension from file name
 *   swap_bytes   - swap bytes of an array
 *   tempfile     - get unique file name
 *   timer_elapsed - get/print the elapsed time since timer_start
 *   timer_start  - (re)start the profiling timer
 *   undersample  - shrink array dimension(s)
 *   width_of     - get 1st dimension of an array
 *   xopen        - extended open with (de)compression, primitives, ...
 *
 * ----------------------------------------------------------------------------
 */

_UTILS_VERSION = "$Date: 2014-10-21 09:07:31 +0200$";

local old_eval;
func eval(eval_code, eval_tmp, eval_debug)
/* DOCUMENT eval, code [, tmp [, debug]];
       -or- eval(code [, tmp [, debug]]);
       -or- old_eval(code, tmp=, debug=)

     Evaluates CODE given as a string or as an array of strings (considered
     as different lines in the script).  Since CODE can be dynamically
     build, this routine allows the execution of virtually (see
     hints/restrictions below) any Yorick's code (e.g. dynamic definition
     of structures, of functions, etc.).  For instance, the following
     statement defines a new structure:

       eval, "struct NewStruct {string a; long b; float c, d;}";

     Since the script gets evaluated at the scope level of the "eval"
     routine some local variables of the "eval" routine may be used in
     the script:

       "eval_tmp"    contains the name of the temporary script file and
                     must not be changed by the script;

       "eval_debug"  contains the value of the keyword DEBUG and must
                     not be changed by the script;

       "eval_code"   contains the value of the argument CODE;

       "eval_result" is returned by "eval", its contents may be defined
                     into the script to provide a returned value.

     Note: impredictible results may occur if CODE changes the value of
     symbols "eval_tmp" and "eval_debug".

     Optional second argument TMP can be used to specify the file name of
     the temporary script.  The default file name is:
       "$TMDIR/$USER_XXXXXX" if environment variable "TMDIR" is set;
       "_eval_XXXXXX"        otherwise;
     where the "XXXXXX" are replaced by a random pattern to avoid
     collisions (see (tempfile).

     If optional third argument DEBUG is true (non-zero and non-nil), the
     name of the temporary file is printed out and the file is not removed.


   SEE ALSO: include, unref, tempfile. */
{
  /* Dump script into a temporary file. */
  if (is_void(eval_tmp)) {
    /* Create default name for yorick temporary code. */
    eval_tmp = get_env("TMPDIR");
    if (eval_tmp) eval_tmp += "/" + get_env("USER") + "_XXXXXX.i";
    else eval_tmp = "_eval_XXXXXX.i";
    eval_tmp = tempfile(eval_tmp);
  }
  write, format="%s\n", open(eval_tmp, "w"), eval_code;

  /* Source script and return result. */
  local eval_result;
  include, eval_tmp, 1;
  if (eval_debug) write, format="Yorick code written in \"%s\"\n", eval_tmp;
  else remove, eval_tmp;
  return eval_result;
}

func old_eval(code, tmp=, debug=)
{
  /* Dump script into a temporary file. */
  if (is_void(tmp)) {
    /* Create default name for yorick temporary code. */
    tmp = get_env("YORICK_EVAL_TMP");
    if (! tmp) {
      tmp = get_env("USER");
      tmp = (tmp ? "/tmp/"+tmp+"-" : "~/.") + "eval_tmp.i";
    }
  }
  write, format="%s\n", open(tmp, "w"), code;

  /* Use "eval_" prefix in order to somewhat protect local variables
     from caller's code. */
  local eval_result, eval_tmp, eval_code, eval_debug;
  eq_nocopy, eval_tmp,   tmp;
  eq_nocopy, eval_debug, debug;
  eq_nocopy, eval_code,  code;

  /* Source script and return result. */
  include, eval_tmp, 1;
  if (eval_debug) write, format="Yorick code written in \"%s\"\n", eval_tmp;
  else remove, eval_tmp;
  return eval_result;
}

func undersample(a, nsub, which=, op=)
/* DOCUMENT undersample(a, nsub)

     Returns array A with all (some) dimensions divided by NSUB.  The
     dimensions of interest must be a multiple of NSUB.

     Keyword WHICH can be used to specify the dimension(s) to shrink.  Values
     in WHICH less or equal zero are counted from the last dimension.  By
     default, all dimensions of A get undersampled.

     Keyword OP can be used to specify the range operator to apply to the sets
     of NSUB adjacent values along the considered dimensions:
       OP=sum   to sum the values
       OP=avg   to average values
       OP=min   to keep the smallest value
       OP=max   to keep the largest value
     By default, the median is taken (WARNING: with the median operator, the
     result depends in which order the dimensions of A are considered).

   SEE ALSO: median. */
{
  if (nsub < 1) error, "NSUB must be >= 1";
  if (nsub == 1) return a;
  if (! is_array(a)) error, "expecting an array";
  rank = dimsof(a)(1);
  if (is_void(which)) {
    which = indgen(rank);
  } else {
    which += (which <= 0)*rank;
    if (structof(which) != long)
      error, "bad data type for keyword WHICH";
    if (min(which) < 1 || max(which) > rank)
      error, "out of range dimension in keyword WHICH";
  }
  nw = numberof(which);
  noop = is_void(op); /* take median value */
  if (! noop && typeof(op) != "range") error, "OP must be nil or a range operator";
  dims = array(rank+1, rank+2);
  for (k = 1; k <= nw; ++k) {
    this = which(k);
    if (this != 1) a = transpose(a, [1,this]);
    dims(2:) = dimsof(a);
    dims(2) = nsub;
    dims(3) = dims(3)/nsub;
    (tmp = array(structof(a), dims))(*) = a(*);
    if (noop) {
      a = median(tmp, 1);
    } else {
      a = tmp(op,..);
    }
    tmp = []; /* free some memory */
    if (this != 1) a = transpose(a, [this,1]);
  }
  return a;
}

func resample(a, f)
/* DOCUMENT resample(a, f);
       -or- resample(a, dims);

     Resample input array A so that its dimensions are rescaled by the factor(s)
     F or become the new dimension list DIMS.

     F is a floating point value or vector which specifies the scaling
     factor(s) to apply to the dimensions of the original array A.  If the
     scaling factors are all the same for all dimensions, F may be a scalar;
     otherwise, F must be a vector with as many elements as the rank of A.
     The dimensions of the result are the dimensions of the input array times
     F and rounded to the nearest integer (though the dimension is at least 1).

     DIMS is an integer value or vector which specifies the dimensions of the
     result.  If the result dimensions are all the same, DIMS may be a scalar;
     otherwise, DIMS is the dimension list of the result which must have the
     same rank as the input array A.

     Along dimensions that are reduced, the sampling is decreased by
     integrating values at a lower resolution; along dimensions that are
     augmented, the sampling is increased by linear interpolating at an higher
     resolution; nothing is done for dimensions that remain the same and have
     a scaling factor equals to 1.

     For instance, to resample an RBG image IMG1 of size 3-by-W1-by-H1 to size
     3-by-W2-by-H2:

       img2 = resample(img1, [3,3,w2,h2])

     if you further need to keep the result of type char (i.e. as an RGB image
     for Yorick):

       img2 = char(max(0, min(255, resample(img1, [3,3,w2,h2]) + 0.5)));


   SEE ALSO: bytescl, dimsof, interp, spline_zoom.
 */
{
  if (! is_array(a) || identof(a) >= Y_STRING) {
    error, "expecting a numerical array";
  }
  if (is_scalar(a)) return a;

  local new_dims, scale;
  old_dims = dimsof(a);
  ndims = old_dims(1);
  err = 0n;
  if (is_integer(f)) {
    if (is_scalar(f) && f >= 1) {
      new_dims = array(long, ndims + 1);
      new_dims(1) = ndims;
      new_dims(2:) = f;
    } else if (is_vector(f) && numberof(f) == ndims + 1 && f(1) == ndims
               && (ndims == 0 || min(f) >= 1)) {
      eq_nocopy, new_dims, f;
    } else {
      err = 1n;
    }
    if (! err) {
      f = double(new_dims(2:))/double(old_dims(2:));
    }
  } else if (is_real(f) && min(f) > 0.0) {
    if (is_scalar(f)) {
      f = array(double(f), ndims);
    } else if (numberof(f) == ndims) {
      eq_nocopy, f, double(f);
    } else {
      err = 1n;
    }
    if (! err) {
      new_dims = array(long, ndims + 1);
      new_dims(1) = ndims;
      new_dims(2:) = max(1, long(floor(f*old_dims(2:) + 0.5)));
    }
  }
  if (err) {
    error, "bad new dimension list or scaling factor(s)";
  }

  /* Resample the original array starting by all dimensions that need to be
     reduced to work with as small arrays as possible. */
  flags = (old_dims(2:) != new_dims(2:) | f != 1.0)*((f < 1.0) + 1);
  for (k = 1; k <= ndims; ++k) {
    if (flags(k) == 2) {
      /* Decrease sampling by integrating values at a lower resolution. */
      q = f(k);
      old_len = old_dims(k+1);
      new_len = new_dims(k+1);
      old_pos = (indgen(old_len + 1) - (old_len + 2)/2.0)*q;
      new_pos = (indgen(new_len + 1) - (new_len + 2)/2.0);
      if (k == 1) {
        a = interp(unref(a)(cum,..), old_pos, new_pos, k)(dif,..)*q;
      } else if (k == 2) {
        a = interp(unref(a)(,cum,..), old_pos, new_pos, k)(,dif,..)*q;
      } else if (k == 3) {
        a = interp(unref(a)(,,cum,..), old_pos, new_pos, k)(,,dif,..)*q;
      } else if (k == 4) {
        a = interp(unref(a)(,,,cum,..), old_pos, new_pos, k)(,,,dif,..)*q;
      } else if (k == 5) {
        a = interp(unref(a)(,,,,cum,..), old_pos, new_pos, k)(,,,,dif,..)*q;
      } else if (k == 6) {
        a = interp(unref(a)(,,,,,cum,..), old_pos, new_pos, k)(,,,,,dif,..)*q;
      } else if (k == 7) {
        a = interp(unref(a)(,,,,,,cum,..), old_pos, new_pos, k)(,,,,,,dif,..)*q;
      } else if (k == 8) {
        a = interp(unref(a)(,,,,,,,cum,..), old_pos, new_pos, k)(,,,,,,,dif,..)*q;
      } else if (k == 9) {
        a = interp(unref(a)(,,,,,,,,cum,..), old_pos, new_pos, k)(,,,,,,,,dif,..)*q;
      } else if (k == 10) {
        a = interp(unref(a)(,,,,,,,,,cum,..), old_pos, new_pos, k)(,,,,,,,,,dif,..)*q;
      } else {
        error, "too many dimensions";
      }
    }
    for (k = 1; k <= ndims; ++k) {
      if (flags(k) == 1) {
        /* Increase sampling by linear interpolating at an higher resolution. */
        q = f(k);
        old_len = old_dims(k+1);
        new_len = new_dims(k+1);
        old_pos = (indgen(old_len) - (old_len + 1)/2.0)*q;
        new_pos = (indgen(new_len) - (new_len + 1)/2.0);
        a = interp(unref(a), old_pos, new_pos, k);
        old_pos = (indgen(old_len + 1) - (old_len + 2)/2.0)*q;
        new_pos = (indgen(new_len + 1) - (new_len + 2)/2.0);
      }
    }
  }
  return a;
}

/*---------------------------------------------------------------------------*/

local nrows_of, ncols_of;
func ndims_of(a)  { return (is_array(a) ? dimsof(a)(1) : -1); }
func width_of(a)  { return (is_array(a) ? (numberof((dims = dimsof(a))) >= 2 ? dims(2) : 1) : 0); }
func height_of(a) { return (is_array(a) ? (numberof((dims = dimsof(a))) >= 3 ? dims(3) : 1) : 0); }
func depth_of(a)  { return (is_array(a) ? (numberof((dims = dimsof(a))) >= 4 ? dims(4) : 1) : 0); }
/* DOCUMENT ndims_of(a) - get number of dimensions of array A; returns -1
 *                        for non-array argument;
 *          nrows_of(a) - get number of rows of array A, returns 0
 *                        for non-array or scalar argument;
 *          ncols_of(a) - get number of columns of array A, returns 0
 *                        for non-array or vector argument;
 *          width_of(a) - get length of 1st dimension of A, returns 0
 *                        if A is non-array is a scalar;
 *         height_of(a) - get length of 2nd dimension of A, returns 0
 *                        if A is non-array or has less than 2 dimensions;
 *          depth_of(a) - get length of 3rd dimension of A, returns 0
 *                        if A is non-array or has less than 3 dimensions.
 *
 *   The number of rows of an array is the length of its first dimension;
 *   the number of columns is the length of its second dimension.
 *
 * SEE ALSO: dimsof.
 */
nrows_of = width_of;
ncols_of = height_of;

/*---------------------------------------------------------------------------*/

func spline_zoom(a, factor, rgb=)
/* DOCUMENT spline_zoom(a, fact)
     Return an array obtained by cubic spline interpolation of A with all
     its dimension multiplied by a factor FACT.  If keyword RGB is true the
     first dimsion of A is left unchanged.  If keyword RGB is not
     specified, then it is considered as true if A is a 3 dimensional array
     of 'char' with its first dimension equal to 3.

     The rescale() function is probably more powerful.

   SEE ALSO: spline, transpose, rescale, resample. */
{
  if (! is_func(spline)) require, "spline.i";
  dims = dimsof(a);
  ndims = dims(1);

  if (is_void(rgb)) {
    /* RGB image? */
    rgb = (structof(a) == char && ndims == 3 && dims(2) == 3);
  }
  if (rgb) {
    a = transpose(a, 0);
    k = 1;
  } else {
    k = 0;
  }
  while (++k <= ndims) {
    dims = dimsof(a);
    n0 = dims(2);
    n1 = max(long(factor*n0 + 0.5), 1);
    if (n1 == 1) {
      a = a(avg,..)(..,-);
    } else {
      if (n1 != n0) {
        dims(2) = n1;
        b = array(double, dims);
        x0 = (indgen(n0) - (n0 + 1)/2.0)/n0;
        x1 = (indgen(n1) - (n1 + 1)/2.0)/n1;
        n = numberof(a)/n0;
        for (i = 1; i <= n; ++i) b(,i) = spline(a(,i), x0, x1);
        eq_nocopy, a, b;
      }
      if (ndims > 1) a = transpose(a, 0);
    }
  }
  if (rgb) return char(max(min(floor(a + 0.5), 255.0), 0.0));
  return a;
}

func map(__map__f, __map__x)
/* DOCUMENT map(f, x)
         or map, f, x;
     Map scalar function or subroutine F onto argument X to mimics
     elementwise unary operation.  Argument X must be an array, or
     a list or nil.  When called as a function, the result is: nil
     if X is nil, an array dimsof(F(X(1)))-by-dimsof(X) if X is an
     array, a list if X is a list.

     For example, to kill all windows:
       map, winkill, window_list();

     To print all the elements of a list:
       map, write, list;

   SEE ALSO: lambda, _map. */
{
  /* All locals here must have weird names, since the user's function
     may rely on external variables for arguments not varying in the
     source, or for accumulated outputs. */
  if (is_array(__map__x)) {
    if (am_subroutine()) {
      __map__n = numberof(__map__x);
      for (__map__i = 1; __map__i <= __map__n; ++__map__i) {
        __map__f, __map__x(__map__i);
      }
    } else if (is_scalar(__map__x)) {
      return __map__f(__map__x);
    } else {
      __map__n = numberof(__map__x);
      __map__y1 = __map__f(__map__x(1));
      __map__y = array(structof(__map__y1), dimsof(__map__y1),
                       dimsof(__map__x));
      if ((__map__m = numberof(__map__y1)) == 1) {
        __map__y(1) = __map__y1;
        for (__map__i = 2; __map__i <= __map__n; ++__map__i) {
          __map__y(__map__i) = __map__f(__map__x(__map__i));
        }
      } else {
        __map__y(1:(__map__j1 = __map__m)) = __map__y1(*);
        for (__map__i = 2; __map__i <= __map__n; ++__map__i) {
          __map__j0 = __map__j1 + 1;
          __map__j1 += __map__m;
          __map__y(__map__j0:__map__j1) = __map__f(__map__x(__map__i))(*);
        }
      }
      return __map__y;
    }
  } else if (is_list(__map__x)) {
    if (am_subroutine()) {
      while (__map__x) {
        __map__f, _car(__map__x);
        __map__x = _cdr(__map__x);
      }
    } else {
      __map__y = _lst(__map__f(_car(__map__x)));
      while ((__map__x = _cdr(__map__x))) {
        _cat, __map__y, __map__f(_car(__map__x));
      }
      return __map__y;
    }
  } else if (! is_void(__map__x)) {
    error, "unsupported data type \""+typeof(__map__x)+"\"";
  }
}

func lambda(args, code)
/* DOCUMENT lambda(args, code);
     Return an anonymous function with ARGS its argment list and CODE
     the body of the function.  ARGS and CODE must be scalar strings.
     For instance:

       f1 = lambda("x", "c = x*x; return sqrt(c + abs(x));");
       f2 = lambda("x,y", "return cos(x*y + abs(x));");

     define two functions f1 and f2 which take respectively one and two
     arguments.  When variables f1 and f2 get out of scope the function
     definition is automatically deleted.

     Other example:

       a = _lst(12,34,67);
       b = map(lambda("x", "return sin(x);"), a);

     B is a list with its elements the sines of the elements of A.

   SEE ALSO: map, include.
 */
{
  include,["func __lambda__(" + args + ") { " + code + " }"], 1;
  return __lambda__;
}

func rescale(a, .., scale=, rgb=, cubic=)
/* DOCUMENT rescale(a, dimlist)
       -or- rescale(a, scale=FACT)
     Return an array obtained by interpolation of A with new dimension list
     as given by DIMLIST or with all its dimension multiplied by a scaling
     factor FACT if keyword SCALE is specified.  If keyword RGB is true the
     first dimension of A and of the interpolated array must be 3 and the
     interpolated array is converted into char.  If keyword CUBIC is true
     cubic spline interpolation is used.

   SEE ALSO: interp, spline, transpose, cast, resample. */
{
  /* explicit */ extern spline;

  /* Get dimension lists. */
  if (! is_array(a)) error, "unexpected non-array argument";
  old_dimlist = dimsof(a);
  ndims = old_dimlist(1);
  if (rgb && (ndims != 3 || old_dimlist(2) != 3)) {
    error, "bad dimension list for RGB image";
  }
  if (is_void(scale)) {
    /* Build new dimension list from remaining arguments. */
    new_dimlist = [0];
    while (more_args()) make_dimlist, new_dimlist, next_arg();
    if (new_dimlist(1) != ndims || (rgb && new_dimlist(2) != 3)) {
      error, "incompatible new dimension list";
    }
    scale = double(new_dimlist)/double(old_dimlist);
  } else {
    /* Scale all dimensions. */
    if (more_args()) error, "two many arguments with scale option";
    if (! is_scalar(scale) || ! (is_real(scale) || is_integer(scale)) ||
        scale <= 0) {
      error, "SCALE must be a strictly positive scalar";
    }
    new_dimlist = max(long(scale*old_dimlist + 0.5), 1);
    new_dimlist(1) = ndims;
    scale = array(double(scale), ndims + 1);
    if (rgb) new_dimlist(2) = 3;
  }
  if (rgb) scale(2) = 1.0;

  type = structof(a);
  cmplx = (type == complex);
  if (! cmplx && type != double) a = double(a);
  n1 = (cmplx ? 2 : 1);
  n2 = 1;
  n3 = numberof(a);
  if (rgb) {
    n = old_dimlist(2);
    n1 *= n;
    n3 /= n;
    j = 2;
  } else {
    j = 1;
  }
  for ( ; j <= ndims; ++j) {
    jp1 = j + 1;
    old_dim = old_dimlist(jp1);
    new_dim = new_dimlist(jp1);
    n3 /= old_dim;
    if (new_dim != old_dim) {
      dims = [3, n1, old_dim, n3];
      if (new_dim == 1) {
        a = cast(a, double, dims)(,avg,);
      } else {
        a = cast(a, double, dims);
        old_x = (indgen(old_dim) - (old_dim + 1)/2.0)*scale(jp1);
        new_x = (indgen(new_dim) - (new_dim + 1)/2.0);
        if (cubic) {
          dims(3) = new_dim;
          b = array(double, dims);
          for (i3 = 1; i3 <= n3; ++i3) {
            for (i1 = 1; i1 <= n1; ++i1) {
              b(i1,,i3) = spline(a(i1,,i3), old_x, new_x);
            }
          }
          eq_nocopy, a, b;
        } else {
          a = interp(a, old_x, new_x, 2);
        }
      }
    }
    n1 *= new_dim;
  }

  if (cmplx) {
    a = cast(a, complex, new_dimlist);
  } else if (anyof(new_dimlist != old_dimlist)) {
    a = cast(a, double, new_dimlist);
  }
  return a;
}

func cast(a, type, dimlist)
/* DOCUMENT cast(a, type, dims)
     This function returns array A reshaped to an array with given TYPE and
     dimension list.

   SEE ALSO: reshape.
 */
{
  local r;
  reshape, r, &a, type, dimlist;
  return r;
}

func swap_bytes(a)
/* DOCUMENT swap_bytes(a)
     Swap the bytes of the array A.

   SEE ALSO: swap, cast.
 */
{
  type = structof(a);
  size = sizeof(type);
  if (size == 1) return a;
  dims = dimsof(a);
  a = cast(a, char, [2, size, numberof(a)]);
  for (i = 1, j = size; i < j; ++i, --j) {
    t = a(i,);
    a(i,) = a(j,);
    a(j,) = t;
  }
  return cast(a, type, dims);
}

/*---------------------------------------------------------------------------*/
/* STRING ROUTINES */

func strlower(s) { return strcase(0, s); }
func strupper(s) { return strcase(1, s); }
/* DOCUMENT strlower(s)
       -or- strupper(s)
     Convert (array of) string(s) S to lower/upper case letters.

   SEE ALSO strcase */

func strcut(str, len)
/* DOCUMENT strcut(str, len)
     Cut input scalar string STR in pieces of length less or equal LEN and
     return an array of such pieces.

   SEE ALSO strjoin */
{
  if ((str_len = strlen(str))<=len) return str;
  n = (str_len+len-1)/len;
  result = array(string, n);
  for (i = 1, i1 = 1, i2 = len; i <= n; ++i, i1 += len, i2 += len)
    result(i)= strpart(str, i1:i2);
  return result;
}

func strjoin(str, glue)
/* DOCUMENT strjoin(str)
       -or- strjoin(str, glue)
     Join strings from array STR into a single long string.  The string GLUE
     (default "") is used between each pair of element from STR.

   SEE ALSO strcut */
{
  if (structof(str) != string) error, "expecting a string argument";
  if (glue && ((n = numberof(str)) > 1)) {
    return str(1) + sum(glue + str(2:n));
  }
  return sum(str);
}

local strchr, strrchr; /* required for the documentation */
/* DOCUMENT strchr(s, c)
         or strrchr(s, c)
     The strchr (strrchr) function returns the index of the first (last)
     occurrence of the character C in the string S.  If C is not found an
     index of zero is returned.  S can be an array of strings, in which case
     the result is an array of long's.  C must be a scalar character or a
     scalar string of length 1).

   SEE ALSO strfind.
*/

func strchr(s, c, back)
{
  if (structof(c) == char) c = strchar(c);
  sel = strfind(c, s, back=back);
  if (dimsof(sel)(1) == 1) {
    /* workaround fast scalar "bug" */
    return (sel(2) < 0 ? 0 : sel(1) + 1);
  }
  i = sel(1, ..) + 1;
  test = (sel(2, ..) < 0);
  if (anyof(test)) i(where(test)) = 0;
  return i;
}

func strrchr(s, c)
{
  return strchr(s, c, 1n);
}

/*---------------------------------------------------------------------------*/
/* LOGICAL ROUTINES */

func is_integer_scalar(x) { return (is_integer(x) && is_scalar(x)); }
func is_string_scalar(x) { return (is_string(x) && is_scalar(x)); }
/* DOCUMENT is_integer_scalar(x)
       -or- is_string_scalar(x)
     Check whether or not X is an integer/string scalar.

   SEE ALSO is_scalar, is_integer. */

/*---------------------------------------------------------------------------*/
/* FILE ROUTINES */

func is_absolute_path(name)
/* DOCUMENT is_absolute_path(name)
     Check whether NAME is an absolute path to a file.  Argument NAME must be
     a string or char array.  If NAME is a string array, the result has same
     dimensionality as NAME; otherwise, the result is a scalar.  This function
     should behaves as what is internally assumed by Yorick.

   SEE ALSO: strchar.
 */
{
  if ((type = structof(name)) == string) {
    result = array(int, dimsof(name));
    n = numberof(name);
    for (j = numberof(name); j >= 1; --j) {
      result(j) = __is_absolute_path(strchar(name(j)));
    }
    return result;
  }
  if (type == char) {
    return __is_absolute_path(name);
  }
  error, "expecting a string or char array";
}

func __is_absolute_path(name)
/** DOCUMENT private function, argument must be a char array

   SEE ALSO: is_absolute_path.
 */
{
  n = numberof(name);
  if (n >= 1) {
    if (! (c1 = name(1))) return 0n;
    dirsep = '/';
    if (c1 == dirsep || c1 == '~') return 1n;
    if (n >= 2) {
      if (! (c2 = name(2))) return 0n;
      /* handle MS Windows a:/blahblah */
      if (n >= 3 && c2 == ':' && name(3) == dirsep &&
          ((c1 >= 'A' && c1 <= 'Z') || (c1 >= 'a' && c1 <= 'z'))) return 1n;
      /* handle MS Windows \\server\name */
      if (c1 == '\\' && c2 == '\\') return 1n;
    }
  }
  return 0n;
}

func locate(name, split=, path=)
/* DOCUMENT locate(name)
         or locate, name;

     Returns the full path of an existing file.  NAME is the file name as a
     scalar string.  If no file is found (see below) or if caller has not read
     access to the file, nil is returned.  When called as a subroutine, the
     result, if any, is simply printed on standard output (keyword SPLIT is
     ignored).

     If keyword SPLIT is true and file is readable, the result is the 2
     element array [DIRNAME,BASENAME] where DIRNAME is the expanded directory
     name of the file (with a terminating /) and BASENAME is the name of the
     file without the directory part.

     If NAME is a relative path and keyword PATH is unset, then the file
     loaction is relative to the current working directory; otherwise, PATH
     gives a list of possible relative directories for the file location.
     PATH can be an array of strings with individual paths separated by ':' or
     an integer: PATH=1 to use Yorick search path given by get_path() and
     PATH=2 to use get_env("PATH").

   EXAMPLES
     locate("~/.bashrc")          --> "/home/eric/.bashrc"
     locate("~/.bashrc", split=1) --> ["/home/eric/", ".bashrc"]
     locate("bash", path=2)       --> "/bin/bash"
     locate("dawson.i", path=1)   --> "/usr/local/libexec/yorick/i/dawson.i"

   SEE ALSO: cd, get_env, get_path, open, is_absolute_path.
 */
{
  /* Check name. */
  if (is_void(name) || structof(name) != string || dimsof(name)(1) != 0) {
    error, "expecting a scalar string";
  }

  /* Split name. */
  dirsep = "/";
  j = strfind(dirsep, name, back=1)(2);
  if (j < 0) {
    dirname = "." + dirsep;
    basename = strpart(name, 1:0);
    absolute = 0n;
  } else {
    dirname = strpart(name, 1:j);
    basename = strpart(name, j+1:0);
    absolute = __is_absolute_path(strchar(dirname));
  }

  /* Find file and expand directory name. */
  cwd = get_cwd();
  if (absolute || is_void(path)) {
    temp = cd(dirname);
    if (! is_void(temp)) {
      cd, cwd;
      dirname = temp;
    }
    if (! open(dirname + basename, "r", 1)) {
      return;
    }
  } else {
    /* Use list of search directories. */
    if (structof(path) != string) {
      if (path == 1) path = get_path();
      else if (path == 2) path = get_env("PATH");
      else error, "bad value for keyword PATH";
    }
    /* Split the list. */
    temp = strchar(path);
    j = where(temp == ':');
    if (is_array(j)) {
      temp(j) = '\0';
      path = strchar(temp);
    }
    /* Find the file. */
    dirsep = "/";
    n = numberof(path);
    for (j = 1; j <= n; ++j) {
      dir = path(j);
      if (strpart(dir, 0:0) != dirsep) dir += dirsep;
      if (open(dir + name, "r", 1)) {
        /* File exists and is readable.  Expand the directory
           part of the name. */
        temp = cd(dir + dirname);
        if (is_void(temp))  error, "unexpected result (BUG)";
        cd, cwd;
        dirname = temp;
        break;
      }
    }
    if (j > n) return;
  }
  if (am_subroutine()) {
    write, format="\"%s\"\n", dirname + basename;
  } else {
    if (split) return [dirname, basename];
    return dirname + basename;
  }
}

func load(__c6ce8719)
/* DOCUMENT load, name;
         or load(name);
     Include Yorick script immediately, taking care of pre-inserting the
     directory where is the script to the plug-in directory list.  When called
     as a function, returns the full path to the script or nil is not found.
     When called as a subroutine, an error is raised if the script is not
     found.

   SEE ALSO: include, locate, plug_dir.
 */
{
  /* Since the "included" script will see the local variables here, I must use
   * dummy variable names to avoid namespace corruption:
   * __c6ce8719 = NAME
   * __1f465328 = [DIRNAME, BASENAME], then FULLNAME
   * __5f90ec32 = DIRLIST
   */

  /* Check name. */
  if (is_void(__c6ce8719) || structof(__c6ce8719) != string
      || dimsof(__c6ce8719)(1) != 0 || strpart(__c6ce8719, -1:0) != ".i") {
    error, "expecting a script name like sample.i";
  }

  /* Find the script and source it with plugin directory pre-set. */
  __1f465328 = locate(__c6ce8719, path=2, split=1);
  if (is_void(__1f465328)) {
    if (am_subroutine()) error, "include file \""+__c6ce8719+"\" not found";
    return;
  }
  __5f90ec32 = plug_dir();
  if (is_void(__5f90ec32) || __5f90ec32(1) != __1f465328(1)) {
    plug_dir, grow(__1f465328(1), __5f90ec32);
  } else {
    __5f90ec32 = [];
  }
  __1f465328 = sum(__1f465328);
  include, __1f465328, 1;
  if (is_array(__5f90ec32)) plug_dir, __5f90ec32;
  return __1f465328;
}

local _TEMPFILE_ALPHABET, _TEMPFILE_SEED;
func tempfile(template)
/* DOCUMENT tempfile(template)
     Returns a file name build from TEMPLATE and which did not exists when the
     function checked.  If the string "XXXXXX" is found in the file part of
     TEMPLATE (that is after the last / if any), these characters get replaced
     with a pseudo-random string; otherwise, the pseudo-random string is simply
     appended to TEMPLATE.  The pseudo-random string is chosen so as to make
     the filename unique.  There is however a very small chance that the
     returned filename is not unique.  To limit conflicts, an empty file with
     the same name as the returned value is created; this file can be deleted
     or overwritten by the caller.  The caller is responsible to delete the
     temporary file (with remove) when no longer needed.

     Note that the tempfile function uses its own internal random generator to
     avoid changing the sequence of random values returned by Yorick's builtin
     random generator.  When called as a subroutine, the internal random
     generator is (re)initialized.

   SEE ALSO: open, remove.
 */
{
  extern _TEMPFILE_ALPHABET, _TEMPFILE_SEED;

  if (am_subroutine()) {
    /* Initialization: seed random generator. */
#if 0 /* _read is *not* designed to read from devices */
    file = open("/dev/urandom", "rb", 1);
    if (file) {
      /* plan A: use system random generator */
      t = array(char, 4);
      if (_read(file, 0L, t) == 4) {
        m = 256.0; /* multiplier */
        _TEMPFILE_SEED = t(1) + m*(t(2) + m*(t(3) + m*t(4)));
        return;
      }
    }
#endif
    /* plan B: use elapsed time in microseconds to seed random generator */
    n = 2.0^32; /* number of values */
    r = 1e6; /* microseconds */
    t = array(double, 3);
    timer, t;
    _TEMPFILE_SEED = floor((sum(t) % (n/r))*r + 0.5) % n;
    return;
  }

  /* Locate substitution pattern in template name. */
  pattern = "XXXXXX";
  off = strfind("/", template, back=1)(1);
  sel = strfind(pattern, template, off, back=1);
  if (sel(2) > sel(1)) {
    i1 = sel(1) + 1;
    i2 = sel(2);
    buf = strchar(template);
  } else {
    buf = strchar(template + pattern);
    i1 = strlen(template) + 1;
    i2 = i1 + strlen(pattern) - 1;
  }

  /* Randomly generate a filename from the template. */
  n = numberof(_TEMPFILE_ALPHABET);
  do {
    /* 32-bit pseudo-random number generator */
    _TEMPFILE_SEED = (1664525.0*_TEMPFILE_SEED + 1013904223.0)%4294967296.0;
    x = _TEMPFILE_SEED;
    for (i = i1; i <= i2; ++i) {
      r = x % n;
      buf(i) = _TEMPFILE_ALPHABET(long(r + 1.5));
      x = (x - r)/n;
    }
    filename = strchar(buf);
  } while (open(filename, "r", 1));
  open, filename, "w"; /* create empty file */
  return filename;
}
if (structof(_TEMPFILE_SEED) != double
    || dimsof(_TEMPFILE_SEED)(1) != 0) {
  tempfile; /* setup random generator */
}
_TEMPFILE_ALPHABET = ['0','1','2','3','4','5','6','7','8','9',
                      'A','B','C','D','E','F','G','H','I','J','K','L','M',
                      'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
                      'a','b','c','d','e','f','g','h','i','j','k','l','m',
                      'n','o','p','q','r','s','t','u','v','w','x','y','z'];

func filesize(filename, errmode)
/* DOCUMENT filesize(filename)
 *     -or- filesize(filename, errmode)
 *
 *   Returns the size (in bytes) of file FILENAME which must be readable.
 *   If ERRMODE is non-nil and non-zero, fail by returning -1L, otherwise
 *   failure to open the file in read mode is a runtime error.
 *
 * SEE ALSO: open, _read.
 */
{
  /* Open file in binary mode and get the file size by a bisection method,
     attempting to read a single byte at various offsets.  To avoid integer
     overflows, all offset computations are performed in double
     precision. */
  stream = open(filename, "rb", 1);
  if (! stream) {
    if (errmode) return -1L;
    error, "cannot read file";
  }
  byte = char(0);
  nbits = 8*sizeof(long) - 1;
  max_offset = 2.0^nbits - 1.0;
  min_size =  0.0;
  max_size = -1.0; /* negative until maximum size is detected */
  offset = 16384.0;
  gain = 2.0;
  for (;;) {
    if (_read(stream, long(offset), byte)) {
      min_size = offset + 1.0;
    } else {
      max_size = offset;
    }
    if (max_size >= min_size) {
      /* Maximum size has been set: check for convergence, otherwise
         take the safeguarded bissection step. */
      if (min_size == max_size) {
        return long(max_size);
      }
      offset = min(max_offset, floor((min_size + max_size)/2.0));
    } else {
      /* Maximum size has not yet been set: grow the size. */
      offset = min(max_offset, floor(min_size*gain));
    }
  }
}

func pwd(nil)
/* DOCUMENT pwd
       -or- pwd()
     Prints out (subroutine form) or returns (function form) full path
     of current working directory.

   SEE ALSO: cd, lsdir. */
{
  if (! is_void(nil)) error, "unexpected non-nil argument";
  dir = cd(".");
  if (am_subroutine()) write, format="%s\n", dir;
  else return dir;
}

func glob(pat)
/* DOCUMENT glob(pat)
 *   This function returns a list of files matching glob-style pattern PAT.
 *   Only the 'file' part of PAT can have wild characters (the file part is
 *   after the last '/' in PAT).  If no files match PAT, nil is returned.
 *
 * SEE ALSO: lsdir, strglob.
 */
{
  i = strfind("/", pat, back=1);
  if ((i = i(2)) >= 1) {
    dir = strpart(pat, 1:i);
    pat = strpart(pat, i+1:0);
  } else {
    dir = "./";
  }
  list = lsdir(dir);
  if (structof(list) == string) {
    list = list(where(strglob(pat, list)));
    if (! is_void(list)) return dir + list;
  }
}

func dirname(path)
/* DOCUMENT dirname(path)
     Returns PATH with its  trailing "/component" removed; if PATH contains
     no /'s, returns  "./" (meaning the current directory).   The result is
     always terminated  by a "/", so that  dirname(dirname(PATH)) gives the
     same result as dirname(PATH).

   SEE ALSO: basename, strrchr. */
{
  if (! (i = strrchr(path, '/'))) return "./";
  return strpart(path, 1:i);
}

func basename(path, ext)
/* DOCUMENT basename(path)
       -or- basename(path, suffix)
     Returns  PATH  with  any  leading directory  components  removed.   If
     specified, also remove  a trailing SUFFIX if the  tail of PATH matches
     SUFFIX.   Arguments PATH  and  SUFFIX, if  specified,  must be  scalar
     strings.

   SEE ALSO: dirname, strrchr. */
{
  if ((i = strrchr(path, '/'))) {
    path = strpart(path, i+1:0);
  }
  if (ext) {
    if ((n2 = strlen(path)) >= (n1 = strlen(ext))
        && strpart(path, n2-n1+1:n2) == ext) {
      return strpart(path, 1:n2-n1);
    }
  }
  return path;
}

func strip_file_extension(path, ext, ..)
/* DOCUMENT strip_file_extension(path);
 *     -or- strip_file_extension(path, ext, ..);
 *
 *   Strip file extension from filename PATH (an array of strings).
 *   If no other arguments are specified, the result is PATH with
 *   trailing characters after (and including) the last dot '.'
 *   stripped (the last dot must however occur after the last slash
 *   '/').  Otherwise, there may be any number of extensions EXT which
 *   are tried in turn until one matches the end of PATH, in that case
 *   the result is PATH stripped from that particular extension (which
 *   may or may not contain a dot).  If PATH is an array, the same
 *   processus is applied to every element of PATH (i.e. they can
 *   match a different extension).
 *
 *
 * SEE ALSO: strpart, strgrep, streplace.
 */
{
  if (! is_string(path)) {
    error, "expecting a string";
  }
  while (more_args()) {
    grow, ext, next_arg()(*);
  }
  if (is_void(ext)) {
    return streplace(path, strgrep("\\.[^./]*$", path), "");
  }
  result = path; /* make a private copy */
  path_len = strlen(path);
  n = numberof(ext);
  ext_len = strlen(ext);
  for (j = numberof(path); j >= 1; --j) {
    len1 = path_len(j);
    path1 = path(j);
    for (k = 1; k <= n; ++k) {
      if ((len2 = ext_len(k)) <= len1 && len2 > 0 &&
          strpart(path1, 1 - len2 : 0) == ext(k)) {
        result(j) = strpart(path1, 1 : -len2);
        break;
      }
    }
    if (k > n) {
      result(j) = path1;
    }
  }
  return result;
}

func expand_file_name(path)
/* DOCUMENT expand_file_name(path)
     Expand leading "~" in file name PATH which must be an array of
     strings (or a scalar string).
   SEE ALSO: strpart, cd, get_cwd, get_file_name, protect_file_name. */
{
  select = strglob("~*", path);
  if (noneof(select)) return path;
  result = unref(path); /* avoid a copy if PATH is temporary */
  select = where(select);
  cwd = get_cwd(); /* memorize current working directory */
  head = string();
  for (k = numberof(select); k >= 1; --k) {
    i = select(k);
    name = result(i);
    sread, name, format="%[^/]", head;
    home = cd(head);
    if (home) {
      cd, cwd; /* restore working directory */
      if ((off = strlen(head) + 2) <= strlen(name)) {
        result(i) = home + strpart(name, off:0);
      } else {
        result(i) = home;
      }
    }
  }
  return result;
}

func get_file_name(obj)
/* DOCUMENT get_file_name(obj)
     If OBJ is a stream, returns the path name of the file associated with
     the stream.  Otherwise if OBJ is an array of strings, expand leading
     "~" in the elements of OBJ.

   SEE ALSO: open, print, expand_file_name, protect_file_name. */
{
  /* Expand leading "~" if string object. */
  if (structof(obj) == string) {
    return expand_file_name(obj);
  }

  /* Check input and get description of stream by the print() command. */
  if ((id = typeof(obj)) == "stream") {
    id = 1;
    s = print(obj);
  } else if (id == "text_stream") {
    id = 2;
    s = print(obj)(2:);
  } else {
    error, "unexpected non-string, non-stream argument";
  }

  /* Join backslash terminated lines from print() result (another
     possibility would be to change the line length with `print_format' but
     there is no way to restore the previous line_lenght unles we building
     a wrapper around original `print_format' routine and make a
     substitution). */
  join = (strpart(s, 0:0) == "\\");
  if (anyof(join)) {
    r = array(string, (ns = numberof(s)) - sum(join) + join(0));
    i = j = 0;
    while (i < ns) {
      w = s(++i);
      while (join(i)) {
        w = strpart(w, :-1);
        if (++i > ns) break;
        w += s(i);
      }
      r(++j) = w;
    }
    s = r;
    w = r = [];
  }

  /* Recover the full path of the stream file from the joined lines. */
  if (id == 1) {
    /* Binary stream. */
    if (numberof(s) == 2) {
      w1 = w2 = string(0);
      if (sread(s(1), format="%[^:]", w1) == 1 &&
          sread(s(2), format="%[^/]", w2) == 1) {
        return strpart(s(2), strlen(w2)+1:0) + strpart(s(1), strlen(w1)+3:0);
      }
    }
    error, "unexpected binary stream descriptor";
  } else {
    /* Text stream. */
    if (numberof(s) == 1) {
      w = string(0);
      if (sread(s(1), format="%[^/]", w) == 1) {
        return strpart(s(1), strlen(w)+1:0);
      }
    }
    error, "unexpected text stream descriptor";
  }
}

local _protect_file_name_table, _protect_file_name_list;
func protect_file_name(path)
/* DOCUMENT protect_file_name(path)
     Protect special characters in PATH (apostrophes, quotes, $, *, ?, [,
     etc) by backslashes to avoid them being interpreted by the shell, for
     instance when using the system() builtin function.
   SEE ALSO: system, get_file_name, expand_file_name. */
{
  c = *pointer(path);
  n = numberof(c) - 1; /* same as strlen(path) */
  p = _protect_file_name_table(1 + c);
  r = array(char, 2*n + 1);
  i = 0;
  j = 0;
  while (++i <= n) {
    if (p(i)) r(++j) = '\\';
    r(++j) = c(i);
  }
  return string(&r);
}
_protect_file_name_list = ['$','&','!','#','?','*',
                           '[',']','{','}','<','>',
                           '\\','"','\'',
                           ' ','\t','\r','\n','\v','\f'];
(_protect_file_name_table = array(char, 256))(1 + _protect_file_name_list) = 1;

func read_ascii(file, compress=, maxcols=, skip=)
/* DOCUMENT arr = read_ascii(file_or_name);

     Reads ascii numeric  data in columns from text file.   FILE_OR_NAME is the
     name  of the  file  or an  already  open  file stream.   The  result is  a
     NCOLS-by-NROWS array of doubles.

     Data are read as double values arranged in columns separated by any number
     of spaces  or tabs.  Comments starting  with a "#" or  any other character
     which is not  part of a number  are ignored up to  the end-of-line.  Blank
     lines are ignored.  The first non-blank/commented line gives the number of
     values per column,  for subsequent lines.  Subsequent lines  must have the
     same number  of columns --  blanks in columns  are not permitted,  use 0.0
     instead.  However, minimal error checking is performed, and if the data is
     not really in columns, read_ascii can silently fail to interpret your file
     as you would scanning it by eye.

     Keyword SKIP can be  set with the number of lines  to skip before starting
     the parsing.

     The read operation will be much faster if the number of commented lines is
     relatively small.  Blank lines cost  nothing, while a line containing just
     a "#" is expensive.

     If the file is  specified by its name, it may be  compressed in which case
     it is automatically decompressed while  reading (see xopen).  The value of
     keyword COMPRESS ("auto" by default) is passed to xopen.

     For  very large  data file,  keyword MAXCOLS  can be  used to  specify the
     expected maximum number of columns (10000 by default).

   SEE ALSO: xopen, read, raw_read. */
{
  /* open the file if it's not already open */
  if (structof(file) == string)
    file = xopen(file, compress=(is_void(compress) ? "auto" : compress));

  /* read lines one at a time until the "model" line which
   * determines the number of columns is discovered
   * assume the number of columns is less than MAXCOLS */
  if (! is_void(skip)) {
    rdline, file, skip;
  }
  x = array(double, (is_void(maxcols) ? 10000 : maxcols));
  ncols = 0;
  while ((line = rdline(file))) {
    ncols = sread(line, x);
    if (ncols) break;          /* got a line with numbers */
  }
  if (! ncols) return [];

  nrows = 1;
  list = _lst([x(1:ncols)]);
  x = array(double, ncols, 10000/ncols + 1);
  for(;;) {
    /* try to grab at least 10000 numbers from the file
     * blank lines will be skipped, but any comments will
     * interrupt the read */
    if (! (n = read(file, x))) {
      /* if didn't get any, drop back to reading comments one
       * line at a time until we get some more numbers */
      while ((line = rdline(file))) {
        if ((n = sread(line, x))) break;
      }
      if (! line) break;    /* rdline detected end-of-file, n==0 too */
    }
    if (n%ncols) error, "data is not in columns";
    n /= ncols;

    /* grow the list the fast way, adding new values to its head
     * (adding to the tail would make growth an n^2 proposition,
     *  as would using the grow function) */
    list = _cat(x(,1:n), list);
    nrows += n;
  }

  /* pop chunks off list and reassemble result */
  x = array(0.0, ncols, nrows);
  for (i = nrows; list; list = _cdr(list)) {
    n = numberof(_car(list))/ncols;
    x(,i-n+1:i) = _car(list);
    i -= n;
  }

  return x;
}

func load_text(file, compress=)
/* DOCUMENT load_text(file)
     Returns all lines of text file as a vector of strings.  Returns nil if
     there are no lines to read.  FILE can be a file name or a text stream.
     In the latter case, the lines not yet read get returned.

     By  default,  if  the file  is  specified  by  its  name, it  will  be
     automatically decompressed if its  compressed; it is possible to force
     another  behaviour by  specifying a  value different  from  "auto" for
     keyword COMPRESS (see xopen).

   SEE ALSO: xopen, rdline, dump_text. */
{
  if (structof(file) == string)
    file = xopen(file, compress=(is_void(compress) ? "auto" : compress));
  text = rdline(file, 1000);
  while (text(0)) grow, text, rdline(file, numberof(text));
  if (is_array((j = where(text)))) return text(j);
}

func dump_text(file, text, compress=, level=, preserve=)
/* DOCUMENT dump_text, file, text;
     Dump every  elements of string array  TEXT as individual  lines into a
     text file.  FILE can be a file name or a text stream.

     If the file is specified by  its name, keywords COMPRESS, LEVEL can be
     used to specify  a compression method (see xopen).   The default value
     of COMPRESS is  "auto" (i.e., method guessed on the  basis of the file
     extension).

     If the file  is specified by its name, keyword PRESERVE  can be set to
     true to avoid overwriting an existing file.

   SEE ALSO: xopen, rdline, load_text. */

{
  if ((anything = ! is_void(text)) && structof(text) != string)
    error, "expecting array of strings or nil";
  if (structof(file) == string) {
    file = xopen(file, "w", compress=(is_void(compress) ? "auto" : compress),
                 level=level, preserve=preserve);
  }
  if (anything) write, file, format="%s\n", text;
}

func guess_compression(filename)
/* DOCUMENT guess_compression, filename;
       -or- guess_compression(filename)
     Guess  which compression  program was  used to  produce  file FILENAME
     according to first  bytes of this file.  When  called as a subroutine,
     the file name is printed out  with the name of the compression program
     if any.  If called as a function, the result is an integer:
       1 - if file compressed with "gzip";
       2 - if file compressed with "bzip2";
       3 - if file compressed with "compress";
       4 - if file compressed with "pack";
       5 - if file compressed with "lzma";
       6 - if file compressed with "xz";
       0 - otherwise.

  SEE ALSO: xopen. */
{
  /* according to information in /etc/magic:
   *
   * method    min.            bytes
   *           size      1    2    3    4    5    6
   * --------- ----    ---- ---- ---- ---- ---- ----
   * pack        3     \037 \036
   * compress    3     \037 \235
   * gzip       20     \037 \213   c         (1)
   * bzip2      14      'B'  'Z'  'h'   c    (2)
   * lzma              \135 \000 \000 \200
   * xz                \375 \067 \172 \130 \132 \000
   *                   \375  '7'  'z'  'X'  'Z' \000
   * 7z                 '7'  'z' \274 \257 \047 \034
   * zip                'P'  'K' \003 \004
   *
   *   (1)  if c<8, compression level, else if c==8 deflated
   *   (2)  with '0' <= c <= '9', block size = c*100kB
   *   minimum size has been computed for an empty file:
   *     pack     ->  3 bytes
   *     compress ->  3 bytes
   *     bzip2    -> 14 bytes
   *     gzip     -> 20 bytes, if compression level is specified
   *                 24 bytes, otherwise
   */
  magic = array(char, 6);
  n = _read(open(filename, "rb"), 0, magic);
  if ((c = magic(1)) == '\037') {
    if ((c = magic(2)) == '\213') {
      if (! am_subroutine()) return 1; /* gzip */
      write, format="%s: compressed with \"gzip\"\n", filename;
      return;
    } else if (c == '\235') {
      if (! am_subroutine()) return 3; /* compress */
      write, format="%s: compressed with \"compress\"\n", filename;
      return;
    } else if (c == '\036') {
      if (! am_subroutine()) return 4; /* pack */
      write, format="%s: compressed with \"pack\"\n", filename;
      return;
    }
  } else if (c == 'B' && magic(2) == 'Z' && magic(3) == 'h' &&
             '0' <= (c = magic(4)) && c <= '9') {
    if (! am_subroutine()) return 2; /* bzip2 */
    write, format="%s: compressed with \"bzip2\" (block size = %d kB)\n",
      filename, 100*(c - '0');
    return;
  } else if (c == '\135' && magic(2) == '\000' && magic(3) == '\000' &&
             magic(4) == '\200') {
    if (! am_subroutine()) return 5; /* lzma */
    write, format="%s: compressed with \"lzma\"\n", filename;
    return;
  } else if (n >= 6 && c == '\375' && magic(2) == '7' && magic(3) == 'z' &&
             magic(4) == 'X' && magic(5) == 'Z' && magic(6) ==  '\000') {
    if (! am_subroutine()) return 6; /* xz */
    write, format="%s: compressed with \"xz\"\n", filename;
    return;
  }
  if (! am_subroutine()) return 0;
  write, format="%s: uncompressed?\n", filename;
}

func xopen(filename, filemode, preserve=, nolog=, compress=, level=, prims=)
/* DOCUMENT xopen(filename)
       -or- xopen(filename, filemode)
     Opens the  file FILENAME according  to FILEMODE (both are  strings).  The
     return value  is an IOStream (or  just stream for short).   When the last
     reference to  this return  value is discarded,  the file will  be closed.
     The file  can also  be explicitly closed  with the close  function (which
     see).   The FILEMODE  (default  "r" --  open  an existing  text file  for
     reading) determines whether  the file is to be opened  in read, write, or
     update mode, and whether writes are restricted to the end-of-file (append
     mode).  FILEMODE  also determines  whether the file  is opened as  a text
     file or as a binary file.   FILEMODE can have the following values, which
     are the same as for the ANSI standard fopen function:
         "r"   - read only
         "w"   - write only, random access, existing file overwritten
         "a"   - write only, forced to end-of-file, existing file preserved
         "r+"  - read/write, random access, existing file preserved
         "w+"  - read/write, random access, existing file overwritten
         "a+"  - read/write, reads random access, writes forced to
                 end-of-file, existing file preserved
         "rb"  "wb"  "ab"  "r+b"  "rb+"  "w+b"  "wb+"  "a+b"  "ab+"
                 without b means text file, with b means binary file

     Keyword COMPRESS  can be  used to specify  compression method for  a text
     file open  for reading (FILEMODE="r")  or writing (FILEMODE="w")  only --
     (de)compression is unsupported  in append mode or for  binary files.  The
     value of keyword COMPRESS can be a scalar string or an integer:
          "auto"     - guess compression according to first bytes of file
                       in read  mode, or according to file extension in
                       write mode: ".gz" for gzip, ".bz2" for bzip2 and
                       ".Z" for compress.
       0  "none"     - no (de)compression
       1  "gzip"     - use gzip to (de)compress
       2  "bzip2"    - use bzip2 to (de)compress
       3  "compress" - use compress to (de)compress
       4  "pack"     - use pack to (de)compress
       5  "lzma"     - use lzma to (de)compress
       6  "xz"       - use xz to (de)compress

     The default value for COMPRESS is "auto" in read mode and "none" in write
     mode.   Note that "gzip",  "bzip2", "pack",  "compress", "lzma"  and "xz"
     commands must exists in your $PATH for compressing with the corresponding
     methods.  Decompression of files compressed with "pack" and "compress" is
     done by "gzip".   If keyword COMPRESS is explicitely  0, no decompression
     is ever  applied; if keyword  COMPRESS is explicitely non-zero,  the file
     must have been compressed.  The  compression level for gzip and bzip2 can
     be specified as an integer value thanks to keyword LEVEL.

     Keyword PRIMS can be used  to specify primitives data type different than
     the  native ones  for binary  files (PRIMS  is ignored  for  text files).
     PRIMS  can be  a scalar  string  (i.e., "alpha",  "mac", "sun3",  "cray",
     "macl", "sun", "dec", "pc", "vax", "vaxg", "i86", "sgi64", or "xdr") or a
     32-element vector of long's as taken by set_primitives (which see).

     When a  binary file is created, it  is possible to avoid  the creation of
     the log-file FILENAME+"L" by setting keyword NOLOG to true.

     Keyword PRESERVE can be set to true to avoid overwriting an existing file
     when FILENAME is open for writing (i.e. with a "w" in FILEMODE).


   RESTRICTIONS:
     If (de)compression  is used, FILENAME  must not contain any  double quote
     character (").


   SEE ALSO: close, guess_compression, open, popen, set_primitives. */
{
  if (is_void(filemode) || filemode == "r") {
    /* Open file for reading in text mode. */
    compress = __xopen_get_compress(compress, filename, 1);
    if (! compress) return open(filename, "r");
    if (compress == 2) return popen("bzip2 -dc \"" + filename + "\"", 0);
    if (compress == 5) return popen("lzma -dc \"" + filename + "\"", 0);
    if (compress == 6) return popen("xz -dc \"" + filename + "\"", 0);
    if (compress == -1) error, "bad value for keyword COMPRESS";
    return popen("gzip -dc \"" + filename + "\"", 0);
  }

  if (preserve && strmatch(filemode, "w") && open(filename, "r", 1))
    error, "file \""+filename+"\" already exists";

  if (filemode == "w") {
    /* Open file for writing in text mode. */
    compress = __xopen_get_compress(compress, filename, 0);
    if (! compress) return open(filename, filemode);
    if (compress == 1) {
      if (is_void(level)) command = swrite(format="gzip > \"%s\"", filename);
      else command = swrite(format="gzip -%d > \"%s\"", level, filename);
    } else if (compress == 2) {
      if (is_void(level)) command = swrite(format="bzip2 > \"%s\"", filename);
      else command = swrite(format="bzip2 -%d > \"%s\"", level, filename);
    } else if (compress == 3) {
      command = swrite(format="compress > \"%s\"", filename);
    } else if (compress == 4) {
      command = swrite(format="pack > \"%s\"", filename);
    } else if (compress == 5) {
      if (is_void(level)) command = swrite(format="lzma > \"%s\"", filename);
      else command = swrite(format="lzma -%d > \"%s\"", level, filename);
    } else if (compress == 6) {
      if (is_void(level)) command = swrite(format="xz > \"%s\"", filename);
      else command = swrite(format="xz -%d > \"%s\"", level, filename);
    } else {
      error, "bad value for keyword COMPRESS";
    }
    return popen(command, 1);
  }

  /* Open file for other modes. */
  if (! (is_void(compress) || compress == 0 || compress == "none"))
    error, "(de)compression unsupported in mode \""+filemode+"\"";
  if (binary && nolog) {
    /* will remove log-file unless it already exists */
    logfile = filename + "L";
    if (open(logfile, "r", 1)) logfile = [];
  } else {
    logfile = [];
  }
  file = open(filename, filemode);
  if (logfile) remove, logfile;

  /* Return open file after optionally installing primitive data types. */
  if (is_void(prims) || ! strmatch(filemode, "b")) return file;
  if ((s = structof(prims)) == string) {
    if (! dimsof(prims)(1)) {
      if (prims != "set" && prims != "get" &&
          is_func((sym = symbol_def(prims+"_primitives"))) == 1) {
        sym, file;
        return file;
      } else if (is_array((sym = symbol_def("__"+prims))) &&
                 structof(sym) == long && numberof(sym) == 32 &&
                 dimsof(sym)(1) == 1) {
        set_primitives, file, sym;
        return file;
      }
    } else if (s == long && numberof(prims) == 32 && dimsof(prims)(1) == 1) {
      set_primitives, file, prims;
      return file;
    }
  }
  error, "bad value for keyword PRIMS";
}

func __xopen_get_compress(compress, filename, for_reading)
/* DOCUMENT __xopen_get_compress(compress, filename, for_reading)
     Private function called by xopen.

   SEE ALSO xopen. */
{
  if (is_void(compress)) {
    if (for_reading) return guess_compression(filename);
    return 0;
  } else if (is_array(compress) && ! dimsof(compress)(1)) {
    if ((s = structof(compress)) == string) {
      if (compress == "auto") {
        if (for_reading) return guess_compression(filename);
        if (strpart(filename, -2:0) == ".gz"  ) return 1; /* gzip */
        if (strpart(filename, -3:0) == ".bz2" ) return 2; /* bzip2 */
        if (strpart(filename, -1:0) == ".Z"   ) return 3; /* compress */
        if (strpart(filename, -4:0) == ".lzma") return 5; /* lzma */
        if (strpart(filename, -2:0) == ".xz"  ) return 6; /* xz */
        return 0;
      }
      if (compress == "none"    ) return 0;
      if (compress == "gzip"    ) return 1;
      if (compress == "bzip2"   ) return 2;
      if (compress == "compress") return 3;
      if (compress == "pack"    ) return 4;
      if (compress == "lzma"    ) return 5;
      if (compress == "xz"      ) return 6;
    } else if ((s == long || s == int || s == short || s == char) &&
               compress >= 0 && compress <= 4) {
      return compress;
    }
  }
  return -1;
}

func raw_read(filename, type, .., encoding=, offset=)
/* DOCUMENT raw_read(filename, type, dimlist, ...)

     Read binary array of TYPE  elements with DIMLIST dimension list into file
     FILENAME and return the array.

     Keyword OFFSET can  be used to set  the number of bytes to  skip prior to
     reading the data.

     Keyword ENCODING  can be used  to change the  data encoding of  the file.
     The value of the keyword is a string like:

        "xdr", "sun"  - eXternal Data Representation (IEEE big endian)
         "native"     - native data representation (i.e. no conversion)
         "i86", "pc"  - IEEE little endian machines
         ...

     see  documentation for  "__sun" for  a list  of supported  encodings; the
     default is "native".


   SEE ALSO: open, _read, __sun, make_dimlist, read_ascii.
 */
{
  /* Get dimension list. */
  dims = [0];
  while (more_args()) {
    make_dimlist, dims, next_arg();
  }

  /* Open file. */
  stream = open(filename, "rb");
  if (! is_void(encoding) && encoding != "native") {
    symbol = encoding + "_primitives";
    if (is_func(symbol_exists)) {
      /* Function 'symbol_exists' is provided by Yeti. */
      set_encoding = (symbol_exists(symbol) ? symbol_def(symbol) : -1);
    } else {
      /* symbol_def will raise an error if symbol does not exists */
      set_encoding = symbol_def(symbol);
    }
    if (is_func(set_encoding) != 1) {
      error, "bad encoding \""+encoding+"\"";
    }
    set_encoding, stream;
  }
  save, stream, complex; /* make stream aware of the definition of a complex */

  /* Read data. */
  if (is_void(offset)) {
    offset = 0L;
  }
  data = array(type, dims);
  if (type == char) {
    nbytes = _read(stream, offset, data);
    if (nbytes != numberof(data)) {
      error, "short file";
    }
  } else {
    _read, stream, offset, data;
  }
  return data;
}

/*---------------------------------------------------------------------------*/
/* PARSING OF PASSWORD FILE */

local pw_get_user, pw_get_uid, pw_get_gid;
local pw_get_name, pw_get_home, pw_get_shell;
/* DOCUMENT pw_get_user(id)          // user name
         or pw_get_uid(id)           // user numerical identifier
         or pw_get_gid(id)           // group numerical identifier
         or pw_get_name(id)          // real user name (from GECOS field)
         or pw_get_home(id)          // home directory of user
         or pw_get_shell(id)         // path to shell of user

     These functions return the value(s) of a specific field of password entry
     matching user ID by parsing /etc/passwd file.  If ID is unspecified,
     get_env("USER") is used; otherwise, ID can be a string or a numerical
     identifier.  If ID is specified, it can be a scalar or an array and the
     result has the same geometry as ID.  For a non-existing entry (or empty
     field), the returned value is -1 or the nil-string.

     Keyword PASSWD can be used to specify another location for the password
     file (default: "/etc/passwd").

     Keyword SED can be used to specify the path to the 'sed' command (default:
     "sed").

   SEE ALSO: get_env, popen.
 */
func pw_get_user (id,sed=,passwd=)
{ return _pw_get("s/^\\([^:]*\\).*$/\\1/"); }
func pw_get_uid  (id,sed=,passwd=)
{ return _pw_get("s/^\\([^:]*:\\)\\{2\\}\\([^:]*\\).*$/\\2/", 1); }
func pw_get_gid  (id,sed=,passwd=)
{ return _pw_get("s/^\\([^:]*:\\)\\{3\\}\\([^:]*\\).*$/\\2/", 1); }
func pw_get_name (id,sed=,passwd=)
{ return _pw_get("s/^\\([^:]*:\\)\\{4\\}\\([^,:]*\\).*$/\\2/"); }
func pw_get_home (id,sed=,passwd=)
{ return _pw_get("s/^\\([^:]*:\\)\\{5\\}\\([^:]*\\).*$/\\2/"); }
func pw_get_shell(id,sed=,passwd=)
{ return _pw_get("s/^.*:\\([^:]*\\)$/\\1/"); }
func _pw_get(script, as_integer)
/* DOCUMENT _pw_get(script, as_integer)
     Private function used by pw_get_* functions.
   SEE ALSO: pw_get_real_name. */
{
  extern id, sed, passwd;
  if (is_void(passwd)) passwd = "/etc/passwd";
  if (is_void(sed)) sed = "sed";
  if (is_void(id)) id = get_env("USER");
  if ((s = structof(id)) == string) {
    /* user specified by its name */
    select = swrite(format="^%s:", id);
    result = (as_integer ? array(-1, dimsof(id)) : id);
  } else if (s == long || s == int || s == short || s == char) {
    /* user specified by its numerical ID */
    select = swrite(format="^[^:]*:[^:]*:%d:", id);
    result = array((as_integer ? -1 : string), dimsof(id));
  } else {
    error, "missing or bad user ID";
  }
  if (! open(passwd, "r" , 1)) error, "cannot read password file";
  format = "%s -e '/%s/!d;%s' '%s'";
  n = numberof(id);
  for (i = 1; i <= n; ++i) {
    value = rdline(popen(swrite(format=format, sed,
                                select(i), script, passwd), 0));
    if (strlen(value)) {
      if (as_integer) {
        x = 0;
        if (sread(value, x) == 1) result(i) = x;
      } else {
        result(i) = value;
      }
    }
  }
  return result;
}

/*---------------------------------------------------------------------------*/
/* PDB FILES */

func pdb_list(file)
/* DOCUMENT pdb_list, file;
       -or- pdb_list(file)
     Lists contents of PDB binary file.  FILE can be either a file name or
     a binary stream.

   SEE ALSO: createb, openb, restore, pdb_restore_all. */
{
  if (structof(file) == string) file = openb(file);
  vars = get_vars(file);
  if (! am_subroutine()) return vars;
  title = ["Non-record variables", "    Record variables"];
  for (i = 1; i <= 2; ++i) {
    write, format="%s:", title(i);
    if (numberof(*vars(i))) {
      write, format=" %s", *vars(i);
      write, format="%s\n", "";
    } else {
      write, format="%s\n", " <none>";
    }
  }
}

func pdb_restore_all(_f_i_l_e_)
/* DOCUMENT pdb_restore_all, file;
     Restore all non-record variables of a PDB file.  FILE can be either a
     file name or a binary stream.

   SEE ALSO: createb, openb, restore, pdb_list. */
{
  if (structof(_f_i_l_e_) == string) _f_i_l_e_ = openb(_f_i_l_e_);
  restore, _f_i_l_e_;
}

/*---------------------------------------------------------------------------*/
/* PROFILING ROUTINES */

local _timer_stamp;
func timer_start {
  extern _timer_stamp;
  _timer_stamp = array(double, 3);
  timer, _timer_stamp;
}
func timer_elapsed(count)
/* DOCUMENT timer_start;
       -or- timer_elapsed;
       -or- timer_elapsed, count;
       -or- timer_elapsed()
       -or- timer_elapsed(count)
     The  subroutine  timer_start (re)starts  the  timer  and the  function
     timer_elapsed   computes  the   elapsed  time   since  last   call  to
     timer_start.  If COUNT is given, the elapsed time is divided by COUNT.
     When  called as  a subroutine,  timer_elapsed prints  out  the elapsed
     time; when  called as a  function, it returns  [CPU,SYSTEM,WALL], with
     all three  times measured in seconds.   The two functions  make use of
     external variable _timer_stamp to memorize the initiale times.

     For instance:
       timer_start;
       ...             // some code to be profiled
       timer_elapsed;
       ...             // some more code to be profiled
       timer_elapsed;  // prints out _total_ elapsed time

  SEE ALSO: timer, benchmark. */
{
  extern _timer_stamp;
  elapsed = _timer_stamp;
  timer, elapsed;
  elapsed -= _timer_stamp;
  if (! is_void(count)) elapsed /= count;
  if (am_subroutine()) {
    write, format="cpu=%g, system=%g, wall=%g\n",
      elapsed(1), elapsed(2), elapsed(3);
  } else {
    return elapsed;
  }
}

local __benchmark_proc;
func benchmark(script, repeat)
/* DOCUMENT t = benchmark(script);
         or t = benchmark(script, repeat);
         or benchmark, script;
         or benchmark, script, repeat;

     Measure the  time spent  by executing  Yorick code  in the  SCRIPT string.
     Argument SCRIPT can be a scalar  string or an array of strings (typically,
     one per line of code).  Optional argument REPEAT gives the number of times
     the script is executed; if omitted, its default value is 100.  When called
     as a function,  the returned value is the elapsed  times [CPU,SYS,WALL] in
     seconds  (see `timer`).   When called  as  a subroutine,  the results  are
     printed to standard output.

     The  function `__benchmark_proc`  is created  on  the fly  to perform  the
     benchmark.  A consequence  is that it is not possible  to call `benchmark`
     from the  script.  Also symbols  prefixed by `__benchmark_` should  not be
     used  by  the script  (they  are  reserved  to implement  the  `benchmark`
     function).


   SEE ALSO: timer, timer_start.
 */
{
  extern __benchmark_proc;
  if (is_void(repeat)) repeat = 100;
  if (! is_string(script)) error, "SCRIPT must be a string";
  include, grow("func __benchmark_proc(__benchmark__counter) {",
                "  while (--__benchmark__counter >= 0) { ",
                script,
                "  }",
                "}"), 1;
  t0 = t1 = array(double, 3);
  timer, t0;
  __benchmark_proc, repeat;
  timer, t1;
  __benchmark_proc = [];
  t = (t1 - t0)/repeat;
  if (am_subroutine()) {
    write, format="cpu=%gms, system=%gms, wall=%gms (measured for %d iteration%s)\n",
      t(1)*1e3, t(2)*1e3, t(2)*1e3, repeat, (repeat > 1 ? "s" : "");
  }
  return t;
}

/*---------------------------------------------------------------------------*/

func moments(x, mean, variance)
/* DOCUMENT moments(x)
 *     -or- moments(x, mean)
 *     -or- moments(x, mean, variance)
 *
 *   Returns the first moments of values in array X as:
 *
 *       [MEAN, VARIANCE, SKEWNESS, KURTOSIS]
 *
 *   where MEAN and VARIANCE are the mean and variance of the distribution
 *   sampled by X.  They can be specifed as the 2nd and/or 3rd arguments,
 *   otherwise sample estimates are used:
 *
 *       MEAN = avg(X)
 *
 *       VARIANCE = 1/M * sum((x - MEAN)^2)
 *
 *   where N = numberof(X) and where M = N, if the mean is provided; or
 *   M = N - 1, when the sample mean avg(X) is used.  The other moments
 *   are:
 *
 *       SKEWNESS = 1/N * sum(((x - MEAN)/sqrt(VARIANCE))^3)
 *
 *       KURTOSIS = 1/N * sum(((x - MEAN)/sqrt(VARIANCE))^4) - 3
 *
 *   The skewness is non-dimensional and characterizes the degree of
 *   asymmetry of a distribution around its mean.  For a Gaussian
 *   distribution, the standard deviation of the skewness is sqrt(15/N)
 *   when the true mean is used, and sqrt(6/N) when the mean is estimated
 *   by the sample mean.
 *
 *   The kurtosis is non-dimensional and measures the peakedness or
 *   flatness of a distribution with respect to a Normal distribution.  For
 *   a Gaussian distribution, the standard deviation of the kurtosis is
 *   sqrt(96/N) when the true variance is known and sqrt(24/N) when it is
 *   the sample estimate.  A distribution with positive kurtosis is termed
 *   'leptokurtic' and is more peaked than the Normal distribution.  A
 *   distribution with negative kurtosis is termed 'platykurtic' and is
 *   more flat than the Normal distribution.  An in-between distribution is
 *   termed 'mesokurtic'.
 *
 * SEE ALSO: avg, rms, median.
 */
{
  n = numberof(x);
  if (is_void(mean)) {
    mean = avg(x);
    m = n - 1;
  } else {
    mean += 0.0; /* make sure it is floating point */
    m = n;
  }
  if (mean) {
    x -= mean;
  }
  if (is_void(variance)) {
    variance = (1.0/m)*sum(x*x);
  }
  if (variance != 1) {
    x *= (1.0/sqrt(variance));
  }
  x3 = x*x*x;
  skewness = (1.0/n)*sum(x3);
  kurtosis  = (1.0/n)*sum(x*x3) - 3.0;
  return [mean, variance, skewness, kurtosis];
}

func _stat_worker(x)
/* DOCUMENT _stat_worker(x)
     Private routine used by stat, returns vector of double's:
        [min(X), max(X), avg(X), std(X)]
     where std(X) is the standard deviation of X.

   SEE ALSO stat. */
{
  if (structof(x)!=double) x = double(x);
  avg_x = avg(x);
  dx = x - avg_x;
  return [min(x), max(x), avg_x, sqrt(avg(dx*dx))];
}

func stat(..)
/* DOCUMENT stat, x, ...
     Print out statistics and information for all the arguments. */
{
  ith = 0;
  while (more_args()) {
    ++ith;
    x = next_arg();
    write, format="%2d: ", ith;
    if (is_array(x)) {
      write, format="array(%s", typeof(x);
      dims = dimsof(x);
      n = numberof(dims);
      for (k = 2; k <= n; ++k) write, format=",%d", dims(k);
      type = structof(x);
      is_numerical = (type==double || type==long || type==int || type==char ||
                     type==complex || type==float || type==short);
      write, format=")%s", (is_numerical ? " " : "\n");
      if (is_numerical) {
        fmt = "min=%g max=%g avg=%g std=%g\n";
        if (type == complex) {
          s = _stat_worker(double(x));
          write, format="\n         real part: "+fmt, s(1), s(2), s(3), s(4);
          s = _stat_worker(x.im);
          write, format="    imaginary part: "+fmt, s(1), s(2), s(3), s(4);
          s = _stat_worker(abs(x));
          write, format="           modulus: "+fmt, s(1), s(2), s(3), s(4);
        } else {
          s = _stat_worker(x);
          write, format=fmt, s(1), s(2), s(3), s(4);
        }
      }
    } else {
      write, format="%s, %s\n", typeof(x), strjoin(print(x));
    }
  }
}

/*---------------------------------------------------------------------------*/

func open_url(url, new=, browser=)
/* DOCUMENT open_url, url;
     Open URL into existing browser.  Keyword NEW can be set to "tab" or
     anything else on-false to open URL into a new tab or a new window.
     Keyword BROWSER can be set to the path of the browser to use (default:
     "firefox").

   SEE ALSO: system.
 */
{
  if (is_void(browser)) browser = "firefox";
  if (new) {
    if (new == "tab")  fmt = "%s -remote 'openURL(%s,new-tab)'";
    else fmt = "%s -remote 'openURL(%s,new-window)'";
  } else {
    fmt = "%s -remote 'openURL(%s)'";
  }
  system, swrite(format=fmt, browser, url);
}

/*---------------------------------------------------------------------------*/

func smooth(a, level)
/* DOCUMENT smooth(a)
       -or- smooth(a, level)
     Returns array A smoothed along its dimensions.  I.e. for a 1D array:
       smooth(A) = A(pcen)(zcen)
     for a 2D array:
       smooth(A) = A(pcen,pcen)(zcen,zcen)
     ... (up to 6 dimensions).

     For a greater number of dimensions,  each  direction  is  smoothed and
     transposed in turn: apart from rounding errors, the result is the same
     but the computation time is approximately  multiplied  by  3.   If you
     oftenly smooth arrays with more than 6 dimensions you may  think about
     modifying the source...

     Optional argument  LEVEL  (default  1)  set  the  number  of  time the
     smoothing operation is performed.

   PROPERTIES OF THE SMOOTHING OPERATOR:
     (i)   The smoothing operator is linear and symmetric.  For instance,
           for a vector, A, smooth(A)=S(,+)*A(+) where the matrix S is
           tridiagonal:
                    [3 1         ]
                    [1 2 1       ]
                    [  1 2 1     ]
             0.25 * [   \ \ \    ]    where, to improve readability,
                    [    \ \ \   ]    missing values are all zero.
                    [     1 2 1  ]
                    [       1 2 1]
                    [         1 3]
           You can, in principle, reverse the smoothing operation with
           TDsolve along each dimensions of smooth(A).  Note:For a vector
           A, the operator S-I applied to A (where I is the identity
           matrix) is the finite difference 2nd derivatives of A (but for
           the edges).

     (ii)  The smoothing operator does not change the sum of the element
           values of its argument, i.e.: sum(smooth(A)) = sum(A).

     (iii) Only an array with all elements having the same value is
           invariant by the smoothing operator.  In fact "slopes" along
           dimensions of A are almost invariant, only the values along the
           edges are changed.

     The symmetry of the smoothing operator is important for the
     computation of gradients.  For instance, let Y = smooth(X) and DQ_DY
     be the gradient of a scalar function Q with respect to Y, then the
     gradient of Q with respect to X is simply: DQ_DX = smooth(DQ_DY)

   TO DO:
     By default A is smoothed along all its dimensions, but the list
     of dimensions to smooth can be specified with keyword WHICH.  As
     usual, negative dimensions are taken as offset from the last one.

     If keyword WRAP is true (non-nil and non-zero) a wrapped version
     of the operator (with same properties but no longer tridiagonal)
     is applied instead.  This is suitable for periodic arrays (e.g.
     FFT transformed arrays).

   SEE ALSO: TDsolve. */
{
  n = dimsof(a)(1);
  if (is_void(level) || level == 1) {
    if (n == 1)
      return a(pcen)(zcen);
    if (n == 2)
      return a(pcen,pcen)(zcen,zcen);
    if (n == 3)
      return a(pcen,pcen,pcen)(zcen,zcen,zcen);
    if (n == 4)
      return a(pcen,pcen,pcen,pcen)(zcen,zcen,zcen,zcen);
    if (n == 5)
      return a(pcen,pcen,pcen,pcen,pcen)(zcen,zcen,zcen,zcen,zcen);
    if (n == 6)
      return a(pcen,pcen,pcen,pcen,pcen,pcen)(zcen,zcen,zcen,zcen,zcen,zcen);
    while (n--)
      a = transpose(a(pcen,..)(zcen,..));
    return a;
  }
  if (n == 1) {
    for (i = 1; i <= level; ++i) {
      a = a(pcen)(zcen);
    }
  } else if (n == 2) {
    for (i=1; i<=level; i++) {
      a = a(pcen,pcen)(zcen,zcen);
    }
  } else if (n == 3) {
    for (i = 1; i <= level; ++i) {
      a = a(pcen,pcen,pcen)(zcen,zcen,zcen);
    }
  } else if (n == 4) {
    for (i = 1; i <= level; ++i) {
      a = a(pcen,pcen,pcen,pcen)(zcen,zcen,zcen,zcen);
    }
  } else if (n == 5) {
    for (i = 1; i <= level; ++i) {
      a = a(pcen,pcen,pcen,pcen,pcen)(zcen,zcen,zcen,zcen,zcen);
    }
  } else if (n == 6) {
    for (i = 1; i <= level; ++i) {
      a = a(pcen,pcen,pcen,pcen,pcen,pcen)(zcen,zcen,zcen,zcen,zcen,zcen);
    }
  } else {
    while (n--) {
      for (i = 1; i <= level; ++i) {
        a = a(pcen,..)(zcen,..);
      }
      a = transpose(a);
    }
  }
  return a;
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
