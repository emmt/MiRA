/*
 * linop.i --
 *
 * General linear operator class for Yeti.
 *
 * This file is part of IPY package, "Inverse Problems with Yorick".
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2007-2012 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 * Copyright (C) 2013 MiTiV project <http://mitiv.univ-lyon1.fr>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more generally,
 * to use and operate it in the same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 */

/* IPY is required */
if (! is_func(ipy_include)) {
  include, "ipy.i", 1;
}

/* Yeti is required. */
ipy_include, (is_func(h_new) == 2), "yeti.i";

local _linop_identity, _linop_diagonal, _linop_sparse, _linop_full;
local _linop_function, _linop_parametric_function;
local LINOP_DIRECT, LINOP_TRANSPOSE, LINOP_INVERSE;
local LINOP_INVERSE_TRANSPOSE, LINOP_TRANSPOSE_INVERSE;
local LINOP_AUTO, LINOP_IDENTITY, LINOP_SCALAR, LINOP_USERFUNC;
local LINOP_DIAGONAL, LINOP_SPARSE, LINOP_FULL;
local linop_new, linop_apply, _linop_type_table;
/* DOCUMENT obj = linop_new();
 *     -or- obj = linop_new(id);
 *     -or- obj = linop_new(id, a);
 *     -or- obj = linop_new(other);
 *     -or- obj = linop_new(s);
 *     -or- obj = linop_new(f);
 *     -or- obj = linop_new(f, p);
 *     -or- obj(x);
 *     -or- obj(x, job);
 *     -or- linop_apply(obj, x);
 *     -or- linop_apply(obj, x, job);
 *
 *   The function linop_new() creates a new linear operator object OBJ which
 *   can be used as a function (or with the linop_apply() function) to compute
 *   the dot product of the 'vector' X by the 'matrix' (or its transpose if JOB
 *   is LINOP_TRANSPOSE) which corresponds to the linear operator.  The
 *   different possibilities to define a linear operator object are described
 *   in the following.  If keyword VIRTUAL is true, the function implementing
 *   the linear operator is not directly referenced but solved whenever the
 *   linear operator is evalutaed the actual function defintion is used (see
 *   symlink_to_name).
 *
 *   The function linop_apply() applies the linear operator OBJ to a 'vector' X
 *   according to the value of JOB as follows (the symbolic names of job
 *   constants are indicated in parenthesis and are recommended to improve the
 *   readability of the code):
 *
 *      JOB = 0 (LINOP_DIRECT) or unspecified to apply the direct operator,
 *            1 (LINOP_TRANSPOSE) to apply the transpose operator,
 *            2 (LINOP_INVERSE) to apply the inverse operator,
 *            3 (LINOP_INVERSE_TRANSPOSE  or LINOP_TRANSPOSE_INVERSE)
 *              to apply the inverse transpose operator.
 *
 *   Thanks to hash-table evaluator, it is not necessary to use linop_apply()
 *   and the usage of the linear operator object is simplified as follows:
 *
 *       OBJ(x)        is the same as: linop_apply(OBJ, x)
 *       OBJ(x, job)   is the same as: linop_apply(OBJ, x, job)
 *
 *
 * RECURSION
 *
 *   If linop_new() is called with a single argument which is already a linear
 *   operator object, then this object is returned (not a copy).  Hence:
 *
 *       obj = linop_new(other);
 *
 *   where OTHER is already a linear operator object, simply creates a new
 *   reference to the object OTHER. (FIXME: we should have the possibility to
 *   'clone' such an object).
 *
 *
 * IDENTITY MATRIX
 *
 *   A linear operator object implementing the identity operation can be
 *   defined one of:
 *
 *       obj = linop_new();
 *       obj = linop_new(LINOP_IDENTITY);
 *       obj = linop_new("identity");
 *
 *
 * MULTIPLICATION BY A SCALAR
 *
 *   A linear operator object implementing the multiplication by a scalar can
 *   be defined by one of:
 *
 *       obj = linop_new(LINOP_SCALAR, a);
 *       obj = linop_new("scalar", a);
 *
 *   where A gives the scalar parameter of the 'matrix'.
 *
 *
 * DIAGONAL MATRIX
 *
 *   A linear operator object implementing the multiplication by a diagonal
 *   matrix can be defined by one of:
 *
 *       obj = linop_new(LINOP_DIAGONAL, a);
 *       obj = linop_new("diagonal", a);
 *
 *   where A gives the coefficients of the diagonal of the 'matrix'; A can be
 *   multi-dimensional, it must however be conformable with the 'vector' X when
 *   linop_apply() is used.
 *
 *
 * FULL MATRIX
 *
 *   A linear operator object can be implemented given a full matrix as:
 *
 *       obj = linop_new(LINOP_FULL, a);
 *       obj = linop_new("full", a);
 *
 *   where A is the array of matrix coefficients which are used as with
 *   mvmult() function (which to see).
 *
 *
 * SPARSE MATRIX
 *
 *   A linear operator object can be implemented as a sparse matrix:
 *
 *       obj = linop_new(s);
 *
 *   where S is a sparse matrix object (see sparse_matrix).
 *
 *
 * USER DEFINED FUNCTION
 *
 *   The linear operator object functionalities (dot product with the
 *   corresponding 'matrix' or its transpose) can be implemented by a user
 *   defined function F.  There are two different possibilities depending
 *   whether or not the function needs aditional data (for instance to store
 *   the coefficients of the 'matrix').
 *
 *   If argument P is omitted, the pseudo-code for F must be:
 *
 *       func f(x, job) {
 *         if (! job) return A.x;
 *         else if (job == 1) return A'.x;
 *         else if (job == 2) return (1/A).x;
 *         else if (job == 3) return (1/A)'.x;
 *         error, "unsupported value for JOB";
 *       }
 *
 *   where A is the 'matrix' corresponding to the linear operator, where the
 *   dot and the prime indicate dot product and matrix transposition
 *   respectively and where (1/A) indicates matrix inverse.  Note that,
 *   depending on your needs, not all operations must be implemented in the
 *   function F.  For instance, if only direct and matrix transpose products
 *   are implemented, the function can be something like:
 *
 *       func f(x, job) {
 *         if (! job) return A.x;
 *         else if (job == 1) return A'.x;
 *         error, "unsupported value for JOB";
 *       }
 *
 *   If argument P is specified, the pseudo-code for F must be:
 *
 *       func f(p, x, job) {
 *         if (! job) return A(p).x;
 *         else if (job == 1) return A(p)'.x;
 *         else if (job == 2) return (1/A(p)).x;
 *         else if (job == 3) return (1/A(p))'.x;
 *         error, "unsupported value for JOB";
 *       }
 *
 *   where A(p) is the 'matrix' which depends on the 'parameters' P.  Note that
 *   this is purely a notation: P can be anything needed by the user-defined
 *   operator.
 *
 *
 * SEE ALSO:
 *   is_linop, sparse_matrix, mvmult, h_new, linop_new_fft, linop_new_fftw.
 */

/*
 * IMPLEMENTATION NOTES
 *
 *   Layout of a linear operator object as build by the linop_new function:
 *       obj.f   = user-defined function
 *       obj.p   = client data for f (or nil)
 *       obj.a   = coefficients of the matrix
 *       obj.s   = sparse matrix object
 *       obj.w   = wrapper function
 */

func linop_new(a1, a2, virtual=) /* DOCUMENTED */
{
  t1 = identof(a1);
  t2 = identof(a2);

  if (is_scalar(a1)) {
    if (t1 == Y_STRING) {
      id = _linop_type_table(a1);
    } else if (Y_CHAR <= t1 && t1 <= Y_LONG) {
      id = long(a1);
    }
    if (id == LINOP_IDENTITY) {
      if (t2 != Y_VOID) {
        error, "exceeding argument to define identity operator";
      }
      return _linop_init(h_new(class="linop", type=LINOP_IDENTITY),
                         "_linop_identity");
    } else if (id == LINOP_SCALAR) {
      if (! is_scalar(a2)) {
        error, "non-calar coefficient";
      }
      if (t2 < Y_CHAR || t2 > Y_COMPLEX) {
        error, "non-numerical scalar coefficient";
      }
      return _linop_init(h_new(class="linop", type=LINOP_SCALAR, a=a2),
                         "_linop_scalar");
    } else if (id == LINOP_DIAGONAL) {
      // FIXME: optimize if all coefficients are zero or one
      if (t2 < Y_CHAR || t2 > Y_COMPLEX) {
        error, "non-numerical diagonal coefficients";
      }
      return _linop_init(h_new(class="linop", type=LINOP_DIAGONAL, a=a2),
                         "_linop_diagonal");
    } else if (id == LINOP_FULL) {
      if (t2 < Y_CHAR || t2 > Y_COMPLEX) {
        error, "non-numerical matrix coefficients";
      }
      return _linop_init(h_new(class="linop", type=LINOP_FULL, a=a2),
                         "_linop_full");
    } else if (id == LINOP_AUTO) {
      t1 = t2;
      a1 = unref(a2);
      t2 = Y_VOID;
    } else {
      error, "invalid linear operator identifier";
    }
  }

  if (t1 == Y_FUNCTION || t1 == Y_BUILTIN) {
    if (t2 == Y_VOID) {
      return _linop_init(h_new(class="linop", type=LINOP_USERFUNC, f=a1),
                         "_linop_function");
    } else {
      return _linop_init(h_new(class="linop", type=LINOP_USERFUNC, f=a1, p=a2),
                             "_linop_parametric_function");
    }
  } else if (t1 == Y_VOID) {
    if (t2 == Y_VOID) {
      return _linop_init(h_new(class="linop", type=LINOP_IDENTITY),
                         "_linop_identity");
    }
  } else if (t1 == Y_OPAQUE) {
    if (is_sparse_matrix(a1)) {
      if (t2 != Y_VOID) {
        error, "exceeding argument to define sparse linear operator";
      }
      // FIXME: not needed for sparse matrix?
      return _linop_init(h_new(class="linop", type=LINOP_SPARSE, s=a1),
                         "_linop_sparse");
    } else if (is_linop(a1)) {
      if (t2 != Y_VOID) {
        error, "exceeding argument";
      }
      return a1;
    }
  }
  error, "bad argument(s) to define sparse linear operator";
}

/* Job values for linear operators: */
LINOP_DIRECT            = 0;
LINOP_TRANSPOSE         = 1;
LINOP_INVERSE           = 2;
LINOP_INVERSE_TRANSPOSE = 3;
LINOP_TRANSPOSE_INVERSE = 3;

/* Type of linear operators: */
LINOP_AUTO     = 0;
LINOP_IDENTITY = 1;
LINOP_DIAGONAL = 2;
LINOP_SPARSE   = 3;
LINOP_FULL     = 4;
LINOP_SCALAR   = 5;
LINOP_USERFUNC = 6;

_linop_type_table = h_new(auto=LINOP_AUTO,
                          identity=LINOP_IDENTITY,
                          diagonal=LINOP_DIAGONAL,
                          sparse=LINOP_SPARSE,
                          full=LINOP_FULL,
                          userfunc=LINOP_USERFUNC,
                          scalar=LINOP_SCALAR);

func linop_apply(obj, x, job) /* DOCUMENTED */
{
  /* Call the wrapper. */
  return obj.w(obj, x, job);
}

func _linop_init(obj, evalname)
{
  extern virtual;
  h_evaluator, obj, evalname;
  return h_set(obj, w=(virtual
                        ? symlink_to_name(evalname)
                        : symbol_def(evalname)));
}

func _linop_identity(obj, x, job)
{
  return x;
}

func _linop_scalar(obj, x, job)
{
  if (! job || job == 1) {
    return obj.a*x;
  } else if (job == 2 || job == 3) {
    return (1.0/obj.a)*x;
  }
  error, "unexpected JOB";
}

func _linop_diagonal(obj, x, job)
{
  if (! job || job == 1) {
    return obj.a*x;
  } else if (job == 2 || job == 3) {
    /* Speed-up: compute/get fast matrix inverse. */
    local ainv; eq_nocopy, ainv, obj.ainv;
    if (is_void(ainv)) {
      ainv = 1.0/obj.a;
      h_set, obj, ainv = ainv;
    }
    return ainv*x;
  }
  error, "unexpected JOB";
}

func _linop_sparse(obj, x, job)
{
  if (! job || job == 1) {
    return obj.s(x, job);
  }
  error, "unsupported value for JOB in sparse linear operator";
}

func _linop_function(obj, x, job)
{
  return obj.f(x, job);
}

func _linop_parametric_function(obj, x, job)
{
  return obj.f(obj.p, x, job);
}

func _linop_full(obj, x, job)
{
  if (! job || job == 1) {
    return mvmult(obj.a, x, job);
  }
  error, "unsupported value for JOB for full matrix linear operator";
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

local linop_is_identity, linop_is_scalar, linop_is_diagonal, linop_is_sparse;
local linop_is_full, linop_is_userfunc, linop_is_virtual;
func is_linop(obj)
/* DOCUMENT is_linop(obj);
 *       or linop_is_identity(obj);
 *       or linop_is_scalar(obj);
 *       or linop_is_diagonal(obj);
 *       or linop_is_sparse(obj);
 *       or linop_is_full(obj);
 *       or linop_is_userfunc(obj);
 *       or linop_is_virtual(obj);
 *
 *   The function is_linop() checks whether object OBJ is a linear operator.
 *
 *   The function linop_is_identity() checks whether object OBJ is a linear
 *   operator which implements the identity.
 *
 *   The function linop_is_scalar() checks whether object OBJ is a linear
 *   operator which implements a multiplication by a scalar.
 *
 *   The function linop_is_diagonal() checks whether object OBJ is a linear
 *   operator which implements a diagonal operator.
 *
 *   The function linop_is_sparse() checks whether object OBJ is a linear
 *   operator which implements a sparse operator.
 *
 *   The function linop_is_full() checks whether object OBJ is a linear
 *   operator which is a full matrix.
 *
 *   The function linop_is_userfunc() checks whether object OBJ is a linear
 *   operator which is implemented by a user-defined function.
 *
 *   The function linop_is_virtual() checks whether object OBJ is a linear
 *   operator implemented by a virtual evaluator.
 *
 *
 * SEE ALSO: linop_new, linop_get_type, linop_get_coefs.
 */
{
  return (is_sparse_matrix(obj) ||
          (is_hash(obj) && is_string((class = obj.class))
           && is_scalar(class) && class == "linop"));
}

func linop_is_virtual(obj) /* DOCUMENTED */
{
  return (is_hash(obj) && obj.class == "linop" && is_symlink(obj.w));
}

func linop_is_identity(obj) /* DOCUMENTED */
{
  return (is_hash(obj) && obj.class == "linop" && obj.type == LINOP_IDENTITY);
}

func linop_is_scalar(obj) /* DOCUMENTED */
{
  return (is_hash(obj) && obj.class == "linop" && obj.type == LINOP_SCALAR);
}

func linop_is_diagonal(obj) /* DOCUMENTED */
{
  return (is_hash(obj) && obj.class == "linop" && obj.type == LINOP_DIAGONAL);
}

func linop_is_sparse(obj) /* DOCUMENTED */
{
  if (is_hash(obj)) {
    return (obj.class == "linop" && obj.type == LINOP_SPARSE);
  } else {
    return is_sparse_matrix(obj);
  }
}

func linop_is_full(obj) /* DOCUMENTED */
{
  return (is_hash(obj) && obj.class == "linop" && obj.type == LINOP_FULL);
}

func linop_is_userfunc(obj) /* DOCUMENTED */
{
  return (is_hash(obj) && obj.class == "linop" && obj.type == LINOP_USERFUNC);
}

func linop_get_type(obj)
/* DOCUMENT linop_get_type(obj);
 *   The function linop_get_type() returns the type of the linear operator OBJ
 *   (one of: LINOP_IDENTITY, LINOP_SCALAR, LINOP_DIAGONAL, LINOP_SPARSE,
 *   LINOP_FULL or LINOP_USERFUNC) or nil if OBJ is not a linear operator.
 *
 * SEE ALSO: linop_new, is_linop.
 */
{
  if (is_hash(obj) && obj.class == "linop") {
    return obj.type;
  } else if (is_sparse_matrix(obj)) {
    return LINOP_SPARSE;
  }
}

func linop_get_coefs(obj)
/* DOCUMENT linop_get_coefs(obj);
 *   The function linop_get_coefs() returns the coefficient(s) of the linear
 *   operator OBJ is it is of type LINOP_SCALAR, LINOP_DIAGONAL, LINOP_SPARSE,
 *   or LINOP_FULL; nothing is returned otherwise.
 *
 * SEE ALSO: linop_new, is_linop, linop_get_type, linop_make_matrix.
 */
{
  if (is_hash(obj) && obj.class == "linop" &&
      ((type = obj.type) == LINOP_SCALAR ||
       type == LINOP_DIAGONAL ||
       type == LINOP_FULL)) {
    return obj.a;
  } else if (is_sparse_matrix(obj)) {
    return obj.coefs;
  }
}

func linop_make_matrix(op, x, job, multi=)
/* DOCUMENT a = linop_make_matrix(op, x);
 *       or a = linop_make_matrix(op, x, job);
 *
 *  Use linear operator OP (see linop_new) with input "vectors" of same data
 *  type (real or complex) and dimension list as X to build a "matrix" A with
 *  the same coefficients as the linear operator OP(x, JOB) -- see linop_new
 *  for the meaning of optional argument JOB.  The result A is always a regular
 *  array.  By default, A is a real M-by-N array where:
 *
 *      M =   numberof(Y)  if Y is real,
 *          2*numberof(Y)  if Y is complex,
 *
 *  where Y = OP(X, JOB), and
 *
 *      N =   numberof(X)  if X is real,
 *          2*numberof(X)  if X is complex.
 *
 *  If keyword MULTI is true, the dimension list of A is 2, if Y is complex,
 *  followed by dimsof(Y), followed by 2, if X is complex, followed by
 *  dimsof(X).
 *
 *
 * SEE ALSO: mvmult.
 */
{
  xident = identof(x);
  if (xident == Y_COMPLEX) {
    xcast = ipy_cast_real_as_complex;
    xdims = make_dimlist(2, dimsof(x));
    n = 2*numberof(x);
  } else if (xident <= Y_DOUBLE) {
    xcast = double;
    xdims = dimsof(x);
    n = numberof(x);
  } else {
    error, "invalid data type for X";
  }
  x = array(double, xdims);
  for (j = 1; j <= n; ++j) {
    x(j) = 1.0;
    if (j == 1) {
      y = op(xcast(x), job);
      yident = identof(y);
      if (yident == Y_COMPLEX) {
        ycast = ipy_cast_complex_as_real;
      } else if (yident <= Y_DOUBLE) {
        ycast = double;
      } else {
        error, "invalid data type for Y";
      }
      y = ycast(y);
      ydims = dimsof(y);
      m = numberof(y);
      a = array(double, m, n);
      a(,j) = y(*);
    } else {
      a(,j) = ycast(op(xcast(x), job))(*);
    }
    x(j) = 0.0;
  }
  if (multi) {
    return ipy_reshape(unref(a), ydims, xdims);
  }
  return a;
}

/*---------------------------------------------------------------------------*/
/* DEPRECATED STUFF */

func linop_new_fftw(nil, dims=, measure=, real=)
/* DOCUMENT obj = linop_new_fftw(...)
 *
 *   Return a new new linear operator object which can be used to compute
 *   FFT by means of FFTW, the "fastest FFT in the world".
 *
 *   Keyword DIMS can be used to pre-specify the dimension list of the
 *   arrays to be transformed.  If left unspecified, the actual
 *   dimension list will be initialized the first time the linear
 *   operator is applied.  In any cases, a given operator can only be
 *   used onto arrays with same dimension lists.
 *
 *   Keywords REAL and MEASURE have the same meaning as for fftw_plan
 *   (which to see).  Note that fftw_plan is only called as needed and
 *   cached into OBJ to save computation time.  Also note that if you
 *   use REAL=1, you must correctly initialize the dimension list of
 *   array to be transformed either when linop_new_fftw is called
 *   or by computing the first FFTW with JOB=0 or 3 (*not* 1 or 2).
 *
 *   Examples:
 *
 *       OBJ = linop_new_fftw();
 *       linop_apply(OBJ, x)    // compute FFT of X
 *       linop_apply(OBJ, x, 0) // idem
 *       linop_apply(OBJ, x, 1) // apply conjugate transpose FFT to X
 *       linop_apply(OBJ, x, 2) // compute inverse FFT of X
 *       linop_apply(OBJ, x, 3) // apply conjugate transpose inverse FFT to X
 *
 *       OBJ.nevals = number of FFT computed so far by OBJ
 *
 *   If you have a recent version of Yeti (which implements hash-table
 *   evaluators), then it is not necessary to use linop_apply and the
 *   usage of the linear operator object is simplified as follows:
 *
 *       OBJ(x)        is the same as: linop_apply(OBJ, x)
 *       OBJ(x, job)   is the same as: linop_apply(OBJ, x, job)
 *
 *
 * SEE ALSO: linop_new, fftw, fftw_plan, linop_new_fft.
 */
{
  if (! is_func(fftw_plan)) {
    include, "yeti_fftw.i", 1;
  }
  if (! is_void(nil)) error, "no non-keyword argument allowed";
  if (is_void(dims)) {
    state = 0;
  } else {
    state = 1;
    for (k=numberof(dims), number=1; k >= 2; --k) {
      number *= dims(k);
    }
    scl = (1.0/number);
  }
  obj = h_new(class="linop", type=LINOP_USERFUNC,
              w=_linop_fftw_wrapper, real=(real ? 1n : 0n),
              dims=dims, scl=scl, measure=measure, nevals=0, state=state);
  h_evaluator, obj, "_linop_fftw_wrapper";
  return obj;
}

func _linop_fftw_wrapper(obj, x, job)
{
  if (! job || job == 3) {
    /* forward transform or backward conjugate transpose */
    if (! ((state = obj.state) & 2)) {
      /* compute forward FFTW plan */
      if (! (state & 1)) {
        h_set, obj, dims=dimsof(x), scl=(1.0/numberof(x));
      }
      h_set, obj, state=(state |= 3),
        fwd=fftw_plan(obj.dims, +1, real=obj.real, measure=obj.measure);
    }
    z = fftw(x, obj.fwd);
    h_set, obj, nevals = obj.nevals + 1;
    return (job == 3 ? obj.scl*z : z);
  } else if (job == 1 || job == 2) {
    /* forward conjugate transpose or backward transform */
    if (! ((state = obj.state) & 4)) {
      /* compute backward FFTW plan */
      if (! (state & 1)) {
        if (obj.real) {
          error, "you must initialize dimension list first (see doc)";
        }
        h_set, obj, dims=dimsof(x), scl=(1.0/numberof(x));
      }
      h_set, obj, state=(state |= 5),
        bck=fftw_plan(obj.dims, -1, real=obj.real, measure=obj.measure);
    }
    z = fftw(x, obj.bck);
    h_set, obj, nevals = obj.nevals + 1;
    return (job == 2 ? obj.scl*z : z);
  }
  error, "unsupported value for JOB in FFTW linear operator";
}

func linop_new_fft(dims, ldir, rdir, real=)
/* DOCUMENT obj = linop_new_fft(dims, ldir, rdir)
 *
 *   Return a new new linear operator object which can be used to compute
 *   FFT by means of Swarztrauber's FFT.  DIMS is the dimension list of the
 *   arrays to be transformed and optional arguments LDIR and RDIR indicate
 *   the dimensions to transform and in which directions (see fft and
 *   fft_setup for more detailed explanations).  The returned operator can
 *   only be used onto arrays with same dimension lists.
 *
 *   Keyword REAL can be set true to specify a real to complex transform.
 *
 *   The FFT linear operator is more flexible than the FFTW one (can
 *   transform for only a subset of the dimensions and with different
 *   directions) but is slower.  Otherwise the two should behave the same
 *   and you can see the documentation of linop_new_fftw for examples.
 *
 *
 * SEE ALSO: linop_new, fft, linop_new_fftw.
 */
{
  real = (real ? 1n : 0n);
  if (is_void(dims)) {
    dims = [0];
  } else if (! is_integer(dims) || ! is_vector(dims) ||
             numberof(dims) != dims(1) + 1 || min(dims) <= 0) {
    error, "invalid dimension list";
  }
  ndims = numberof(dims) - 1;
  dims = long(dims);
  ltyp = _linop_fft_get_dir(ldir);
  rtyp = _linop_fft_get_dir(rdir);
  if (! ltyp && ! rtyp) {
    ltyp = 1;
    ldir = 1;
  } else if (ltyp < 0 || rtyp < 0) {
    error, "bad FFT directions";
  }
  llen = numberof(ldir);
  rlen = numberof(rdir);
  if (llen + rlen > ndims) {
    error, "more FFT directions than number of dimensions";
  }
  if (noneof(ltyp) && noneof(rdir)) {
    wrapper = _linop_identity;
    ndirs = 0;
    scale = 1.0;
    setup = list = length = dirs = top = [];
  } else {
    /* compute stride and number of elements */
    stride = array(long, ndims);
    length = dims(2:0);
    stride(1) = 1;
    for (j = 1; j < ndims; ++j) {
      stride(j + 1) = stride(j)*length(j);
    }
    number = stride(ndims)*length(ndims);

    /* select directions of transform */
    dirs = array(long, ndims);
    if (ltyp == 1 && ! rtyp) {
      dirs(*) = ldir;
    } else {
      if (llen) dirs(1:llen) = ldir;
      if (rlen) dirs(1-rlen:0) = rdir;
    }
    list = where(dirs);
    ndirs = numberof(list);
    length = length(list);
    stride = stride(list);

    /* compute FFT workspaces */
    top = number/(stride*length);
    number = 1;
    setup = array(pointer, ndirs);
    for (j = 1 ; j <= ndirs; ++j) {
      len = length(j);
      number *= len;
      if (j > 1 && (k = where(len == length)(1)) < j) {
        setup(j) = setup(k);
      } else {
        ws = array(double, 6*len + 15);
        fft_init, len, ws;
        setup(j) = &ws;
      }
    }
    scale = 1.0/number; /* scale for inverse FFT */
    wrapper = _linop_fft_wrapper;
  }
  obj = h_new(class="linop", type=LINOP_USERFUNC,
              w=_linop_fft_wrapper, dims=dims, nevals=0, real=real,
              scale=scale, list=list, ndirs=ndirs, dirs=dirs,
              stride=stride, length=length, top=top, setup=setup);
  h_evaluator, obj, "_linop_fft_wrapper";
  return obj;
}

func _linop_fft_get_dir(dir)
{
  if (is_void(dir)) return 0;
  if (is_integer(dir) && min(dir) >= -1 && max(dir) <= +1) {
    if (is_scalar(dir)) return 1;
    if (is_vector(dir)) return 2;
  }
  return -1;
}

func _linop_fft_wrapper(obj, x, job)
{
  local dims, dirs, setup, length, stride, top;
  if (! job || job == 3) {
    real = 0n;
  } else if (job == 1 || job == 2) {
    real = obj.real;
  } else {
    error, "unsupported value for JOB in FFT linear operator";
  }
  if ((type = identof(x)) > Y_COMPLEX) {
    error, "non-numerical argument";
  }
  eq_nocopy, dims, obj.dims;
  if ((xdims = dimsof(x))(1) != dims(1) || anyof(xdims != dims)) {
    error, "incompatible dimensions of argument";
  }
  h_set, obj, nevals = (obj.nevals + 1);
  if (! (ndirs = obj.ndirs)) {
    if (real) {
      if (type == Y_DOUBLE) {
        x = x; /* make a copy */
        return x;
      }
      return double(x);
    } else {
      if (type == Y_COMPLEX) {
        x = x; /* make a copy */
        return x;
      }
      return complex(x);
    }
  }
  if (type == Y_COMPLEX) {
    x = x; /* make a copy for in-place FFT */
  } else {
    x = complex(x);
  }

  /* do the requested transforms in-place */
  if (job == 1 || job == 2) {
    dirs = -obj.dirs;
  } else {
    eq_nocopy, dirs, obj.dirs;
  }
  eq_nocopy, setup, obj.setup;
  eq_nocopy, length, obj.length;
  eq_nocopy, stride, obj.stride;
  eq_nocopy, top, obj.top;
  for (j = 1; j <= ndirs; ++j) {
    fft_raw, dirs(j), x, stride(j), length(j), top(j), setup(j);
  }
  if (real) {
    x = double(x);
  }
  return ((job == 3 || job == 2) ? obj.scale*x : x);
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * End:
 */
