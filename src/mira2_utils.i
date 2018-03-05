/*
 * mira2_utils.i -
 *
 * General purpose functions for MiRA.
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

if (! is_scalar(MIRA_SRCDIR) || ! is_string(MIRA_SRCDIR)) {
  error, "include \"mira2.i\" first";
}
include, MIRA_SRCDIR+"utils.i";

/*---------------------------------------------------------------------------*/
/* UNIVERSAL CONSTANTS AND UNITS */

/* Some constants. */
local MIRA_INF, MIRA_NAN;
local MIRA_PI, MIRA_METER, MIRA_MICRON, MIRA_MICROMETER;
local MIRA_DEGREE, MIRA_DEG, MIRA_ARCSECOND, MIRA_ARCSEC;
local MIRA_MILLIARCSECOND, MIRA_MILLIARCSEC, MIRA_MAS;
/* DOCUMENT Mathematical constants and physical units defined in MiRA.

     +-----------------------------------------------------------+
     | Constants                       Description               |
     +-----------------------------------------------------------+
     | MIRA_PI                         3.1415...                 |
     | MIRA_DEG, MIRA_DEGREE           degree to radian          |
     | MIRA_ARCSEC, MIRA_ARCSECOND     arcsecond in SI units     |
     | MIRA_MAS, MIRA_MILLIARCSEC,     millarcsecond in SI units |
     | MIRA_MILLIARCSECOND                                       |
     | MIRA_MICRON, MIRA_MICROMETER    micrometer in SI units    |
     | MIRA_NANOMETER                  nanometer in SI units     |
     | MIRA_INF                        Inf                       |
     | MIRA_NAN                        quiet NaN                 |
     +-----------------------------------------------------------+

   SEE ALSO: mira, mira_length, mira_angle.
 */
MIRA_PI = 3.141592653589793238462643383279503;
MIRA_TWO_PI = 2.0*MIRA_PI;
MIRA_RAD = MIRA_RADIAN = 1.0;
MIRA_DEG = MIRA_DEGREE = MIRA_PI*MIRA_RADIAN/180.0;
MIRA_ARCSEC = MIRA_ARCSECOND = MIRA_DEGREE/3600.0;
MIRA_MILLIARCSEC = MIRA_MAS = MIRA_MILLIARCSECOND = 1e-3*MIRA_ARCSECOND;
MIRA_METER = 1.0;
MIRA_MICRON = MIRA_MICROMETER = 1e-6*MIRA_METER;
MIRA_NANOMETER = 1e-9*MIRA_METER;
MIRA_LIGHT_YEAR = 9.460730472580800e15; // metres (exactly)
MIRA_ASTRONOMICAL_UNIT = 1.49597870700e11; // metres
MIRA_PARSEC = (648000/MIRA_PI)*MIRA_ASTRONOMICAL_UNIT;

MIRA_INF = ieee_generate(double, 1); // positive infinite
MIRA_NAN = ieee_generate(double, 2); // quiet NAN

local _MIRA_LENGTH_UNITS;
func mira_length(&var, def)
/* DOCUMENT mira_length(var)
         or mira_length(var, def)

     converts value in variable `var` in a length in meters.  If `var` is unset
     on entry, its value is assumed to be `def` (the default value).  Valid
     value can be a real (assumed to be in meters) or a string with human
     readable value and optional units (SI units, that is meters, by default).
     If the value is valid, caller's variable `var` is set to a double
     precision value (in meters) and true is returned.  If the value is
     invalid, an error is thrown if called as a subroutine and false is
     returned if called as a function.

     Examples of valid values:

        "1.67 µm" ---> 1.67e-6
        "3mm" -------> 3.0e-3
        1.3 ---------> 1.3

   SEE ALSO: mira_angle, errs2caller.
 */
{
  if (is_void(var)) eq_nocopy, var, def;
  if (is_scalar(var)) {
    if (is_string(var)) {
      units = string(0);
      value = 0.0;
      n = sread(var, format="%f %[^\a]", value, units);
      if (n == 2) {
        units = strtrim(units, 3, blank=" \t\n\v\f");
        if (h_has(_MIRA_LENGTH_UNITS, units)) {
          var = value*h_get(_MIRA_LENGTH_UNITS, units);
          return 1n;
        }
      } else if (n == 1) {
        /* assume SI units */
        var = double(var);
        return 1n;
      }
    } else if (identof(var) <= Y_DOUBLE) {
      var = double(var);
      return 1n;
    }
  }
  if (am_subroutine()) {
    error, "invalid length";
  }
  return 0n;
}
errs2caller, mira_length;

_MIRA_LENGTH_UNITS = h_new("meter",      1e0,    "m",  1e0,
                           "decimeter",  1e-1,   "dm", 1e-1,
                           "centimeter", 1e-2,   "cm", 1e-2,
                           "millimeter", 1e-3,   "mm", 1e-3,
                           "micrometer", 1e-6,   "µm", 1e-6, "micron", 1e-6,
                           "nanometer",  1e-9,   "nm", 1e-9);
//                         "picometer",  1e-12,  "pm", 1e-12,
//                         "femtometer", 1e-15,  "fm", 1e-15,
//                         "attometer",  1e-18,  "am", 1e-18,
//                         "zeptometer", 1e-21,  "zm", 1e-21,
//                         "yoctometer", 1e-24,  "ym", 1e-24,
//                         "decameter",  1e1,    "dam",1e1,
//                         "hectometer", 1e2,    "hm", 1e2,
//                         "kilometer",  1e3,    "km", 1e3,
//                         "megameter",  1e6,    "Mm", 1e6,
//                         "gigameter",  1e9,    "Gm", 1e9,
//                         "terameter",  1e12,   "Tm", 1e12,
//                         "petameter",  1e15,   "Pm", 1e15,
//                         "exameter",   1e18,   "Em", 1e18,
//                         "zettameter", 1e21,   "Zm", 1e21,
//                         "yottameter", 1e24,   "Ym", 1e24,
//                         "astronomical unit", MIRA_ASTRONOMICAL_UNIT,
//                         "astronomical-unit", MIRA_ASTRONOMICAL_UNIT,
//                         "astronomical_unit", MIRA_ASTRONOMICAL_UNIT,
//                         "au", MIRA_ASTRONOMICAL_UNIT,
//                         "parsec", MIRA_PARSEC,
//                         "pc", MIRA_PARSEC,
//                         "light-year", MIRA_LIGHT_YEAR,
//                         "ly", MIRA_LIGHT_YEAR);

local _MIRA_ANGLE_UNITS;
func mira_angle(&var, def)
/* DOCUMENT mira_angle(var)
         or mira_angle(var, def)

     converts value in variable `var` in a angle in radians.  If `var` is unset
     on entry, its value is assumed to be `def` (the default value).  Valid
     value can be a real (assumed to be in radians) or a string with human
     readable value and optional units (SI units, that is radians, by default).
     If the value is valid, caller's variable `var` is set to a double
     precision value (in radians) and true is returned.  If the value is
     invalid, an error is thrown if called as a subroutine and false is
     returned if called as a function.

     Examples of valid values:

        "3mas" ---------> 1.45444e-08
        "0.7 arcsec" ---> 3.3937e-06
        1.3e-4 ---------> 1.3e-4

   SEE ALSO: mira_length, errs2caller.
 */
{
  if (is_void(var)) eq_nocopy, var, def;
  if (is_scalar(var)) {
    if (is_string(var)) {
      units = string(0);
      value = 0.0;
      n = sread(var, format="%f %[^\a]", value, units);
      if (n == 2) {
        units = strtrim(units, 3, blank=" \t\n\v\f");
        if (h_has(_MIRA_ANGLE_UNITS, units)) {
          var = value*h_get(_MIRA_ANGLE_UNITS, units);
          return 1n;
        }
      } else if (n == 1) {
        /* assume SI units */
        var = double(var);
        return 1n;
      }
    } else if (identof(var) <= Y_DOUBLE) {
      var = double(var);
      return 1n;
    }
  }
  if (am_subroutine()) {
    error, "invalid angle";
  }
  return 0n;
}
errs2caller, mira_angle;

_MIRA_ANGLE_UNITS = h_new("radian", 1e0,
                          "rad", 1e0,
                          "degree", MIRA_DEGREE,
                          "deg", MIRA_DEGREE,
                          "arcsecond", MIRA_ARCSEC,
                          "arcsec", MIRA_ARCSEC,
                          "milliarcsecond", 1e-3*MIRA_ARCSEC,
                          "mas", 1e-3*MIRA_ARCSEC);

/*--------------------------------------------------------------------------*/
/* MISCELLANEOUS */

func mira_sky_coordinates(n, s)
/* DOCUMENT mira_sky_coordinates(n, s);
     Yields a vector of `n` coordinates with step `s` and centered according to
     conventions in `fftshift` and `nfft`.

   SEE ALSO: fftshift, nfft.
 */
{
  return double(s)*indgen(-(n/2):n-1-(n/2))
}

/*--------------------------------------------------------------------------*/
/* HASH TABLES */

local mira_grow_fields;
func mira_grow_field(tbl, key, val)
/* DOCUMENT mira_grow_field, tbl, key, val;
         or mira_grow_fields, tbl, key1, val1, key2, val2, ...;

     The first subroutine, grows (as with the `grow` function) the field KEY of
     the hash table TBL with value VAL.

     The second subroutine, grows (as with the `grow` function) the fields
     KEY1, KEY2, ... of the hash table TBL with corresponding values VAL1,
     VAL2, ...

   SEE ALSO: h_set, h_get, grow.
 */
{
  local oldval;
  eq_nocopy, oldval, h_get(tbl, key);
  h_set, tbl, key, (numberof(oldval) > 0 ? grow(oldval, val) : val);
}

func mira_grow_fields(tbl, ..)
{
  while (more_args()) {
    mira_grow_field, tbl, next_arg(), next_arg();
  }
}

func mira_define_table(args)
/* DOCUMENT mira_define_table, tbl, var1, var2, ..., key1=val1, key2=val2, ...;

     (Re)define a hash table `tbl` from the remaining arguments which are
     interpreted as the contents of a hash table as for the `h_save` function.
     Fist argument `tbl` must be a simple variable.  On entry, if `tbl`
     contains a hash-table, this table is cleared before being filled with the
     new contents; otherwise `tbl` must be undefined and it is set to a fresh
     hash table with the given contents.  This behavior is intended in case
     `tbl` is referenced by other objects and it is desirable to maintain this
     link.

   SEE ALSO: h_save. */
{
  /* Get the keywords and the number of positional arguments. */
  keys = args(-);
  nargs = args(*);
  nkeys = numberof(keys);

  /* Create/fetch the hash table and process positional arguments. */
  local obj, key;
  if (nargs >= 1) {
    eq_nocopy, obj, args(1);
    if (is_hash(obj)) {
      /* Clear current contents. */
      keys = h_keys(obj);
      for (i = numberof(keys); i >= 1; --i) {
        h_pop, obj, keys(i);
      }
    } else if (is_void(obj)) {
      /* Make a fresh table. */
      obj = h_new();
      args, 1, obj;
    }
  }
  if (! is_hash(obj)) {
    error, "expecting a hash table object as first argument or nothing";
  }
  i = 1;
  while (++i <= nargs) {
    eq_nocopy, key, args(-,i);
    if (! key) {
      /* Argument is not a simple variable reference; then, the key is the
         value of the current positional argument and the the value is that of
         the next positional argument. */
      eq_nocopy, key, args(i);
      if (! scalar_string(key)) {
        error, "expecting key name or variable reference";
      }
      if (++i > nargs) {
        error, "missing value after last key name in argument list";
      }
      if (! key) {
        h_evaluator, obj, args(i);
        continue;
      }
    }
    h_set, obj, key, args(i);
  }

  /* Process keywords. */
  for (i = 1; i <= nkeys; ++i) {
    key = keys(i);
    h_set, obj, key, args(key);
  }
  return obj;
}
wrap_args, mira_define_table;

/*--------------------------------------------------------------------------*/
/* GENERAL PURPOSE */

func mira_vector_double(&var, def)
/* DOCUMENT bool = mira_vector_double(var [, def]);
         or mira_vector_double, var [, def];

     Make sure that variable `var` stores a vector of type `double`.  `var`
     should be a variable, not an expression because its contents may be
     redefined by this routine.  On entry, if `var` is undefined, it is set
     with the default value `def`.  On successful return, it is guaranteed that
     `var` stores a vector of type `double`.  When called as a function, a
     boolean result is returned: true in case of success and false in case of
     failure.  When called as a subroutine, an error is thrown in case of
     failure.  See the documentation of `scalar_double` for more details and
     typical usage.

   SEE ALSO: scalar_double.
 */
{
  if (is_void(var)) {
    var = def;
  }
  if ((id = identof(var)) <= Y_DOUBLE) {
    if (is_scalar(var)) {
      var = [double(var)];
      return 1n;
    }
    if (is_vector(var)) {
      if (id != Y_DOUBLE) {
        var = double(var);
      }
      return 1n;
    }
  }
  if (am_subroutine()) error, "expecting an integer";
  return 0n;
}
errs2caller, mira_vector_double;

func mira_same_strings(str1, str2)
/* DOCUMENT mira_same_strings(str1, str2);
     returns whether `str1` and `str2` are identical strings, disregarding the
     leading and trailing spaces and the case of letters.

   SEE ALSO strcase, strtrim.
*/
{
  return (strcase(0n, strtrim(str1)) == strcase(0n, strtrim(str2)));
}

func mira_same_dimensions(adims, bdims)
/* DOCUMENT mira_same_dimensions(adims, bdims);
     returns whether `adims` and `bdims` are identical dimension lists.

   SEE ALSO dimsof.
*/
{
  return (numberof(adims) == numberof(bdims) && allof(adims == bdims));
}

local mira_cast_real_as_complex;
func mira_cast_complex_as_real(z)
/* DOCUMENT z = mira_cast_real_as_complex(x);
         or x = mira_cast_complex_as_real(z);

     The first function converts a 2-by-any real array X into a complex array Z
     such that:

         Z.re = X(1,..)
         Z.im = X(2,..)

     The second function does the inverse operation.

   SEE ALSO: reshape.
 */
{
  if (! is_complex(z)) error, "expecting complex argument";
  reshape, z, &z, double, 2, dimsof(z);
  return z;
}

func mira_cast_real_as_complex(r) /* DOCUMENTED */
{
  if (identof(r) > Y_DOUBLE || (n = numberof((dims = dimsof(r)))) < 2
      || dims(2) != 2) {
    error, "expecting integer or real array with leading dimension equals to 2";
  }
  dims = dims(2:);
  dims(1) = n - 2;
  reshape, r, &double(r), complex, dims;
  return r;
}

local mira_cmult_re, mira_cmult_im;
/* DOCUMENT mira_cmult_re(re1,im1, re2,im2);
         or mira_cmult_im(re1,im1, re2,im2);

     yield the real and imaginary parts of the complex multiplication of
     re1 + 1i*im1 by re2 + 1i*im2.
 */
func mira_cmult_re(re1,im1, re2,im2)
{
  return re1*re2 - im1*im2;
}

func mira_cmult_im(re1,im1, re2,im2)
{
  return re1*im2 + im1*re2;
}

_MIRA_DEBUG_STYLE  = ansi_term(ANSI_TERM_BOLD,ANSI_TERM_FG_CYAN);
func mira_debug(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
{
  str = printf(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
  if (strpart(str, 0:0) == "\n") {
    str = strpart(str, 1:-1);
  }
  write, format="%sDEBUG - %s%s\n", _MIRA_DEBUG_STYLE, str, _MSG_RESET_STYLE;
}

func mira_cat(.., last=)
/* DOCUMENT mira_cat(.., last=0/1);

     yields an array which concatenate all given arguments ignoring void
     arguments.  Non-void input arguments are indexed by the last or first
     index of the returned array depending whether keyword LAST is true or
     false.  Hence, if LAST is true, the result is as with [...], the
     concatenate operator of Yorick.  All non-void input arguments are
     converted to a common type and dimension list.  If all arguments are void,
     a void result is returned.

   SEE ALSO: grow.
 */
{
  local type, dims, rank, arg, args;
  argc = 0;
  while (more_args()) {
    eq_nocopy, arg, next_arg();
    if (is_array(arg) && identof(arg) <= Y_STRING) {
      if (argc == 0) {
        /* Fetch type and dimensions of the first argument. */
        type = structof(arg);
        dims = dimsof(arg);
        rank = numberof(dims);
        args = _lst(arg);
      } else {
        /* Update the type and the dimension list. */
        argtype = structof(arg);
        if (type != argtype) {
          if (type == string) {
            error, "cannot concatenate string and numerical arrays";
          }
          type = structof(type(0) + argtype(0));
        }
        argdims = dimsof(arg);
        argrank = numberof(argdims);
        minrank = min(argrank, rank);
        if (minrank >= 1) {
          /* Common dimensions must be comptatible. */
          a = dims(2:minrank);
          b = argdims(2:minrank);
          if (anyof((a != b) & (min(a,b) != 1))) {
            error, "incompatible dimensions";
          }
          if (argrank > rank) {
            eq_nocopy, dims, argdims;
            rank = argrank;
          }
          dims(2:minrank) = max(a,b);
        } else if (argrank > rank) {
          eq_nocopy, dims, argdims;
          rank = argrank;
        }
        args = _cat(arg, args);
      }
      ++argc;
    } else if (! is_void(arg)) {
      error, "argment(s) must be void or numerical arrays";
    }
  }
  if (argc > 0) {
    if (last) {
      /* Make the argument index the last one. */
      arr = array(type, dims, argc);
      for (k = argc; k > 0; --k) {
        arr(..,k) = _car(args);
        args = _cdr(args);
      }
    } else {
      /* Make the argument index the first one. */
      arr = array(type, argc, dims);
      for (k = argc; k > 0; --k) {
        arr(k,..) = _car(args);
        args = _cdr(args);
      }
    }
    return arr;
  }
}


/*---------------------------------------------------------------------------*/
/* SAFE ATAN(Y,X) */

local mira_atan, _mira_init_atan;
func mira_safe_atan(y, x)
/* DOCUMENT mira_atan(y, x);

     yields the phase of the complex x + i⋅y (both x and y must be real) as
     `atan(y,x)` if (x,y) != (0,0) and yields 0 if (x,y) == (0,0).

     Depending whether the builtin `atan` already has this behavior of if it
     throws an error when (x,y) == (0,0), mira_atan may be defined as atan
     or as mira_safe_atan.  So that the faster version is used if possible.

     According to the man pages of the atan function (in libm, the standard
     mathematical library), it seems guaranteed that atan(y,x) = ±0, whenever x
     = ±0.  So I may be a bit paranoid...

   SEE ALSO: atan.
 */
{
  bad = ((!x)&(!y));
  if (anyof(bad)) {
    dims = dimsof(bad);
    z = array(double, dims);
    nz = where(!bad);
    if (is_array(nz)) {
      if (numberof(x) != numberof(z)) {
        /* broadcast x */
        x += array(double, dims);
      }
      if (numberof(y) != numberof(z)) {
        /* broadcast y */
        y += array(double, dims);
      }
      z(nz) = atan(y(nz), x(nz));
    }
    return z;
  } else {
    return atan(y, x);
  }
}

/* The following is to define mira_atan as atan or as mira_safe_atan depend
   whether atan(0,0)=0 and raises no errors. */
func _mira_init_atan
{
  extern mira_atan;
  if (catch(0x01)) {
    mira_atan = mira_safe_atan;
  } else if (atan(0,0) == 0) {
    mira_atan = atan;
  } else {
    mira_atan = mira_safe_atan;
  }
}
_mira_init_atan;

/*---------------------------------------------------------------------------*/
