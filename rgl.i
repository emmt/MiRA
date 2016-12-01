/*
 * rgl.i --
 *
 * Implement regularization operators for iterative optimization and
 * inverse problems.
 *
 * This file is part of IPY package, "Inverse Problems with Yorick".
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2007-20012 Éric Thiébaut <eric.thiebaut@obs.univ-lyon1.fr>
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

/* Yeti and LINOP are required. */
ipy_include, (is_func(h_new) == 2), "yeti.i";
ipy_include, is_func(linop_new), "linop.i";

/*---------------------------------------------------------------------------*/
/* PUBLIC INTERFACE */

local _RGL_CLASS;
func rgl_info(..)
/* DOCUMENT RGL: Regularization Operators for Iterative
 *               Optimization and Inverse Problems
 *
 *   The following calling sequences are implemented:
 *
 *     rgl_info - Lists all implemented regularization methods and their
 *                configurable attributes (e.g. hyper-parameters).
 *
 *     rgl_info, NAME1, NAME2, ...
 *              - Lists attributes for regularization methods NAME1, NAME2, ...
 *
 *     rgl_info(NAME)
 *              - Returns list of attributes for regularization method NAME.
 *
 *     rgl_info()
 *              - Returns sorted list of all implemented regularization
 *                methods.
 *
 *   The RGL package, implements the following "public" functions, (for more
 *   details about function FN, use: help, FN):
 *
 *      rgl_new - create a new instance of regularization method.
 *      rgl_config - configure regularization parameter(s).
 *
 * SEE ALSO: rgl_new, rgl_config.
 */
{
  local arg, argv;
  while (more_args()) {
    eq_nocopy, arg, next_arg();
    if (is_string(arg)) {
      grow, argv, arg(*);
    } else if (! is_void(arg)) {
      error, "expecting nil or string";
    }
  }
  empty = is_void(argv);
  if (empty) {
    argv = h_keys(_RGL_CLASS);
    if (is_array(argv)) {
      argv = argv(sort(argv));
    }
  }
  argc = numberof(argv);

  if (am_subroutine()) {
    for (i = 1; i <= argc; ++i) {
      name = argv(i);
      class = _RGL_CLASS(name);
      if (is_void(class)) {
        error, ("class \"" + name +
                "\" is not an implemented regularization");
      }
      write, format="%s (", name;
      n = numberof(class.attr_list);
      for (k = 1; k <= n; ++k) {
        write, format="%s%s", (k > 1 ? ", " : ""),
          class.attr_list(k);
      }
      write, format="%s\n", ")";
    }
  }

  if (empty) {
    return argv;
  }
  if (argc == 1) {
    name = argv(1);
    class = _RGL_CLASS(name);
    if (is_void(class)) {
      error, ("class \"" + name +
              "\" is not an implemented regularization");
    }
    list = class.attr_list; /* make a copy */
    return list;
  }
  error, "too many arguments in this context";
}

_RGL_CLASS = h_new(); /* used to store implemented regularization classes */
local rgl_new, _rgl_new;
/* DOCUMENT obj = rgl_new(class);
 *       or obj = rgl_new(class, "attr1", value1, "attr2", value2, ...);
 *       or obj = rgl_new(class, attr1=value1, attr2=value2, ...);
 *
 *   Returns new instance for a specific regularization method.  CLASS can be
 *   a regularization class name (such as "quadratic") or a regularization
 *   class definition (OTHER.class where OTHER is an existing regularization
 *   instance).
 *
 *   Subsequent arguments are keywords or pairs of attribute name and value to
 *   set for the regularizer instance (see rgl_config).
 *
 * SEE ALSO: rgl_config, rgl_update, rgl_get_penalty, rgl_get_gradient,
 *           rgl_get_hessian, rgl_get_diagonal_of_hessian, rgl_apply_hessian,
 *           rgl_get_name.
 */

func _rgl_new(class, ..) /* DOCUMENTED */
{
  /* Get class definition. */
  if (is_string(class) && is_scalar(class)) {
    value = _RGL_CLASS(class);
    if (! _rgl_check_class(value)) {
      error, ("no regularization class match name \"" + class + "\"");
    }
    class = value;
  } else if (! _rgl_check_class(class)) {
    error, "expecting a regularization class name or definition";
  }

  /* Create new object instance and initialize it according to its
     class definition. */
  obj = class.setup(h_new(class=class, mu=1.0));

  /* Apply configuration settings if any. */
  if (more_args()) {
    cfg = h_new();
    while (more_args()) {
      local name;
      eq_nocopy, name, next_arg();
      if (rgl_string_scalar(name)) {
        error, "attribute name must be a string";
      }
      if (! more_args()) {
        error, ("missing value for attribute \"" + name + "\"");
      }
      h_set, cfg, name, next_arg();
    }
    _rgl_config_from_hash, obj, cfg;
  }
  return obj;
}

func rgl_new(args) /* DOCUMENTED */
{
  /* Get number of positional arguments and list of keywords. */
  nargs = args(0);
  key_list = args(-);
  nkeys = numberof(key_list);

  /* Get class definition. */
  local class;
  if (nargs >= 1) eq_nocopy, class, args(1);
  if (is_string(class) && is_scalar(class)) {
    value = _RGL_CLASS(class);
    if (! _rgl_check_class(value)) {
      error, ("no regularization class match name \"" + class + "\"");
    }
    class = value;
  } else if (! _rgl_check_class(class)) {
    error, "expecting a regularization class name or definition";
  }

  /* Create new object instance and initialize it according to its class
     definition and, optionally, to configuration settings specified by the
     other arguments. */
  obj = class.setup(h_new(class=class, mu=1.0));
  if (nargs > 1 || nkeys > 0) {
    cfg = h_new();
    /* Process positional arguments. */
    for (i = 2; i <= nargs; ++i) {
      local key;
      eq_nocopy, key, args(i);
      if (rgl_string_scalar(key)) {
        error, "attribute name must be a string";
      }
      if (++i > nargs) {
        error, ("missing value for attribute \"" + key + "\"");
      }
      h_set, cfg, key, args(i);
    }
    /* Process keywords. */
    for (i = 1; i <= nkeys; ++i) {
      key = key_list(i);
      h_set, cfg, key, args(key);
    }
    /* Apply configuration settings. */
    _rgl_config_from_hash, obj, cfg;
  }
  return obj;
}

local _rgl_config_from_hash, _rgl_config_attribute, _rgl_config_query;
local rgl_config, _rgl_config, _rgl_attribute_index;
/* DOCUMENT rgl_config, this, ..., attr, value, ..., attr=value, ...;
 *     -or- rgl_config, this, ..., cfg, ...;
 *     -or- rgl_config(this);
 *     -or- rgl_config(this, attr);
 *
 *   When called as a subroutine, rgl_config is used to set the value(s) of
 *   configurable attribute(s) of the regularization instance THIS.  Settings
 *   consist in ATTR, VALUE pairs or in CFG hash table.  ATTR and VALUE are an
 *   attribute name (or index) and its new value; keys and associated values in
 *   CFG table are used as attribute names and values.  There may be as many
 *   settings as needed.  If you have a recent version of Yorick (which support
 *   wrap_args), attributes may also be specified as keywords, like in
 *   ATTR=VALUE.
 *
 *   When called as a function, rgl_config is used to query the attributes of
 *   the regularization instance THIS.  Without any other argument, the result
 *   is the list of configurable attribute names. Otherwise, ATTR is the name
 *   (or the index) of an attribute and the result is the current attribute
 *   value.  As a special case, if ATTR is "*" or -1, all the attributes are
 *   returned in the form of a hash-table which can be used to configure
 *   another regularization instance.
 *
 *   For instance:
 *     rgl = rgl_new("xsmooth");
 *     rgl_config, rgl, threshold=1e-4, dimlist=[2,512,512];
 *   or:
 *     rgl_config, rgl, "threshold", 1e-4, "dimlist", [2,512,512];
 *   or:
 *     rgl_config, rgl, h_new(threshold=1e-4, dimlist=[2,512,512]);
 *
 *   The threshold value is obtained by:
 *     rgl_config(rgl, "threshold")
 *
 *   Example to print out hyper-parameter names and types:
 *
 *     list = rgl_config(this);
 *     for (k = 1; k <= numberof(list); ++k) {
 *       write, format="%2d: \"%s\" = %s\n", k, list(k),
 *         pr1(rgl_config(this, k));
 *     }
 *
 *
 * SEE ALSO: rgl_new, rgl_get_global_weight, rgl_set_global_weight, h_new.
 */

func rgl_config(args) /* DOCUMENTED */
{
  /* Get number of positional arguments and list of keywords. */
  nargs = args(0);
  key_list = args(-);
  nkeys = numberof(key_list);

  local this, list, attr;
  if (nargs >= 1) eq_nocopy, this, args(1);
  if (! is_hash(this) || ! h_has(this, "class")) {
    error, "expecting an instance of a regularization class";
  }
  class = this.class; /* shortcut */

  if (am_subroutine()) {
    /* Set attributes from positional arguments and then from keywords. */
    for (i = 2; i <= nargs; ++i) {
      eq_nocopy, attr, args(i);
      if (is_hash(attr)) {
        _rgl_config_from_hash, this, attr;
      } else {
        if (++i > nargs) {
          error, "missing attribute value";
        }
       _rgl_config_attribute, this, attr, args(i);
      }
    }
    for (i = 1; i <= nkeys; ++i) {
      key = key_list(i);
      _rgl_config_attribute, this, key, args(key);
    }
  } else {
    /* Called as function: query attribute(s). */
    if (nargs == 1 && nkeys == 0) {
      list = this.class.attr_list; /* make a copy */
      return list;
    } else if (nargs == 2 && nkeys == 0) {
      /* Query attribute value(s). */
      return _rgl_config_query(this, args(2));
    } else {
      error, "too many arguments or keywords in this context";
    }
  }
}

func _rgl_config(this, ..) /* DOCUMENTED */
{
  local attr;
  if (am_subroutine()) {
    /* Set attributes from other arguments. */
    /* Set attribute(s). */
    while (more_args()) {
      eq_nocopy, attr, next_arg();
      if (is_hash(attr)) {
        _rgl_config_from_hash, this, attr;
      } else {
        if (! more_args()) {
          error, "missing attribute value";
        }
        _rgl_config_attribute, this, attr, next_arg();
      }
    }
  } else {
    /* Called as function: query attribute(s). */
    if (more_args()) {
      /* Query attribute value(s). */
      eq_nocopy, attr, next_arg();
      if (more_args()) {
        error, "too many arguments in this context";
      }
      return _rgl_config_query(this, attr);
    } else {
      list = this.class.attr_list; /* make a copy */
      return list;
    }
  }
}

/* Use old versions of rgl_new and rgl_config if Yorick does not support
   wrap_args. */
if (is_func(wrap_args) == 2) {
  //errs2caller, rgl_new;
  //errs2caller, rgl_config;
  wrap_args, rgl_new;
  wrap_args, rgl_config;
  _rgl_new = [];
  _rgl_config = [];
} else {
  rgl_new = _rgl_new;
  rgl_config = _rgl_config;
}

func _rgl_config_from_hash(this, table)
{
  set_attr = this.class.set_attr; /* shortcut */
  for (key = h_first(table); key ; key = h_next(table, key)) {
    index = _rgl_attribute_index(this, key);
    if (index <= 0) {
      error, ("invalid attribute name \"" + key + "\"");
    }
    set_attr, this, index, table(key);
  }
}

func _rgl_config_attribute(this, attr, value)
{
  index = _rgl_attribute_index(this, attr);
  if (index <= 0) {
    if (! is_scalar(attr) || ! (is_integer(attr) || is_string(attr))) {
      error, "expecting an attribute name (or index) or a hash table";
    }
    error, "invalid attribute name or index";
  }
  set_attr = this.class.set_attr; /* needed to overcome Yorck limitations */
  set_attr, this, index, value;
}

func _rgl_config_query(this, attr)
{
  get_attr = this.class.get_attr; /* shortcut */
  index = _rgl_attribute_index(this, attr, 1n);
  if (index > 0) {
    return get_attr(this, index);
  }
  if (index == -1) {
    /* Return all attributes in the form of a hash-table. */
    result = h_new();
    eq_nocopy, list, class.attr_list;
    n = numberof(list);
    for (k = 1; k <= n; ++k) {
      // FIXME: make a copy here?
      h_set, result, list(k), get_attr(this, k);
    }
    return result;
  }
  error, "invalid attribute name or index";
}
if (is_func(errs2caller) == 2) {
  errs2caller, _rgl_config_from_hash;
  errs2caller, _rgl_config_attribute;
  errs2caller, _rgl_config_query;
}

func _rgl_attribute_index(this, attr, query)
/** DOCUMENT  _rgl_attribute_index(this, attr, query)
 *    Private function to get index of attribute ATTR for regularization
 *    instance THIS.
 * SEE ALSO: rgl_config.
 */
{
  if (is_scalar(attr)) {
    if (is_integer(attr)) {
      if ((index = long(attr)) > 0 && index <= this.class.attr_number) {
        return index;
      }
      if (query && index == -1) {
        return -1;
      }
    } else if (is_string(attr)) {
      index = this.class.attr_table(attr);
      if (index) {
        return index;
      }
      if (query && attr == "*") {
        return -1;
      }
    }
  }
  return 0;
}

local rgl_get_global_weight, rgl_set_global_weight;
/* DOCUMENT mu = rgl_get_global_weight(this);
 *     -or- rgl_set_global_weight, this, mu;
 *
 *   Query or set the regularization global weight MU from/in regularization
 *   instance THIS.  Note that, by definition, the global regularization
 *   weight is the first hyper-parameter in THIS, therefore these routines are
 *   simple shortcuts for:
 *
 *       mu = rgl_config(this, 1);
 *       rgl_config, this, 1, mu;
 *
 *
 * SEE ALSO: rgl_new, rgl_config.
 */

func rgl_get_global_weight(this)
{
  return this.class.get_attr(this, 1);
}

func rgl_set_global_weight(this, value)
{
  return this.class.set_attr(this, 1, value);
}

local rgl_update, rgl_get_penalty;
local rgl_get_gradient, rgl_apply_hessian;
local rgl_get_hessian, rgl_get_diagonal_of_hessian;
/* DOCUMENT rgl_update, this, x;
 *     -or- f = rgl_get_penalty(this, x);
 *     -or- g = rgl_get_gradient(this, x);
 *     -or- p = rgl_apply_hessian(this, x, s);
 *     -or- d = rgl_get_diagonal_of_hessian(this, x);
 *     -or- h = rgl_get_hessian(this, x);
 *
 *   The subroutine rgl_update setup internals of regularization instance THIS
 *   to account for a change of parameters, X are the new parameters.  Note
 *   that the setup may be delayed (for efficiency reasons) until the first
 *   call to one of the subsequent routines.  The subroutine rgl_update must
 *   be called *before* rgl_get_penalty, rgl_get_gradient, rgl_apply_hessian,
 *   etc.  This is mandatory to allow the caching of intermediate quantities
 *   and to speed up computations when the penalty and/or the gradient and/or
 *   the Hessian are needed for the same set of parameters X.  It is also
 *   mandatory that the following functions get called with the same
 *   parameters X as rgl_update.
 *
 *   F is the value of the regularization at X.
 *
 *   G is the gradient of the regularization at X.
 *
 *   P is the result of applying the Hessian approximation of the
 *   regularization at X to the parameter step S.
 *
 *   D is the diagonal approximation of the Hessian of the regularization at
 *   X.
 *
 *   H is the approximation of the Hessian of the regularization at X.
 *
 *
 * SEE ALSO:
 *   rgl_new, rgl_config.
 */

func rgl_update(this, x)
{
  return this.class.update(this, x);
}

func rgl_get_penalty(this, x)
{
  return this.class.get_penalty(this, x);
}

func rgl_get_gradient(this, x)
{
  return this.class.get_gradient(this, x);
}

func rgl_get_hessian(this, x)
{
  return this.class.get_hessian(this, x);
}

func rgl_get_diagonal_of_hessian(this, x)
{
  return this.class.get_diagonal_of_hessian(this, x);
}

func rgl_apply_hessian(this, x, s)
{
  return this.class.apply_hessian(this, x, s);
}

func rgl_get_name(this) { return this.class.name; }
/* DOCUMENT rgl_get_name(this)
 *
 *   Returns name of regularization method implemented by
 *   regularization instance THIS.
 *
 * SEE ALSO: rgl_new.
 */

func _rgl_check_class(class)
{
  return (is_hash(class) && is_symlink((setup = class.setup)) &&
          is_func(value_of_symlink(setup)) && is_string(class.name) &&
          is_scalar(class.name));
}

/*---------------------------------------------------------------------------*/
/* CLASS BUILDER */

_RGL_METHODS = ["update", "get_penalty", "get_gradient", "get_hessian",
                "get_diagonal_of_hessian", "apply_hessian",
                "get_attr", "set_attr", "setup"];

func _rgl_class(name, attr)
/* DOCUMENT class = _rgl_class(name, attr)
 *   Build class definition for regularization method.  NAME is a scalar
 *   string which identifies the regularization method.  ATTR is an array of
 *   strings with attribute (or hyper-parameter) names.
 *
 * SEE ALSO: rgl_new.
 */
{
  class = h_new(name=name);
  bogus = "_rgl_bogus_";
  suffix = "_rgl_" + name + "_";

  /* Link to class-methods. */
  get_method = symlink_to_name; // can also be: symbol_def
  for (k = numberof(_RGL_METHODS); k >= 1; --k) {
    method = _RGL_METHODS(k);
    method_name = suffix + method;
    if (! symbol_exists(method_name)) {
      method_name = bogus + method;
    }
    h_set, class, method, get_method(method_name);
  }

  /* List of hyperparameters. */
  table = h_new();
  for (k = numberof(attr); k >= 1; --k) {
    h_set, table, attr(k), k;
  }
  h_set, class, attr_number = numberof(attr),
    attr_table = table, attr_list = attr;

  /* Remember class name. */
  h_set, _RGL_CLASS, name, class;

  return class;
}

func _rgl_bogus_setup(this)
{
  return this;
}

func _rgl_bogus_not_implemented(this, method)
{
  return ("method \"" + method +
          "\" not implemented in regularization class \"" +
          this.class.name + "\"");
}

func _rgl_bogus_update(this, x)
{
  error, _rgl_bogus_not_implemented(this, "update");
}

func _rgl_bogus_get_penalty(this, x)
{
  error, _rgl_bogus_not_implemented(this, "get_penalty");
}

func _rgl_bogus_get_gradient(this, x)
{
  error, _rgl_bogus_not_implemented(this, "get_gradient");
}

func _rgl_bogus_get_hessian(this, x)
{
  error, _rgl_bogus_not_implemented(this, "get_hessian");
}

func _rgl_bogus_get_diagonal_of_hessian(this, x)
{
  error, _rgl_bogus_not_implemented(this, "get_diagonal_of_hessian");
}

func _rgl_bogus_apply_hessian(this, x)
{
  error, _rgl_bogus_not_implemented(this, "apply_hessian");
}

func _rgl_bogus_get_attr(this, x)
{
  error, _rgl_bogus_not_implemented(this, "get_attr");
}

func _rgl_bogus_set_attr(this, x)
{
  error, _rgl_bogus_not_implemented(this, "set_attr");
}

func _rgl_get_symbol(str)
/* DOCUMENT _rgl_get_symbol(str);
 *   Returns symbol with name STR or nil if it doesn't exists.  This function
 *   is mainly a workaround of builtin symbol_def which raises an error if STR
 *   was never defined.
 *
 * SEE ALSO: symbol_exists, symbol_def.
 */
{
  if (symbol_exists(str)) {
    return symbol_def(str);
  }
}

func _rgl_symlink_to_name(str) { return link(str + ""); }
/* DOCUMENT _rgl_symlink_to_name(str);
 *   The _rgl_symlink_to_name function returns symbolic link to value of STR
 *   which must be a scalar string.  This is a workaround the link function
 *   (which see) which returns a link to the symbol's name rather than the
 *   symbol's value when called with a simple symbol as argument (i.e. not an
 *   expression).
 *
 * SEE ALSO: link, _rgl_class.
 */
if (is_func(symlink_to_name) == 2 && is_func(is_symlink) == 2 &&
    is_func(name_of_symlink) == 2 && is_func(value_of_symlink) == 2) {
  /* Yeti 6.2.1 and newer */
  _rgl_symlink_to_name = [];
} else if (is_func(link) == 2 && is_func(is_link) == 2 &&
           is_func(link_name) == 2 && is_func(solve_link) == 2) {
  /* Yeti 6.2.0 */
  symlink_to_name = _rgl_symlink_to_name;
  symlink_to_variable = link;
  name_of_symlink = link_name;
  value_of_symlink = solve_link;
  is_symlink = is_link;
} else {
  error, "symbolic link not implemented (upgrade Yeti)";
}

/*---------------------------------------------------------------------------*/
/* QUADRATIC SMOOTHNESS */

local rgl_smoothness;
/* DOCUMENT this = rgl_new("smoothness");
 *
 *   Creates a regularization instance suitable for quadratic smoothness
 *   regularization.  The regularization penalty for an array X is:
 *
 *       mu*||D.x||^2
 *
 *   where mu is the global regularization weight and D is a linear operator
 *   such that: D.x = x - smooth3(x).
 *
 *   This regularization has one hyper-parameter:
 *
 *     1 - "mu" = global regularization weight.
 *
 * SEE ALSO: rgl_new, smooth3.
 */

func _rgl_smoothness_setup(this)
{
  return h_set(this, state = 0);
}

func _rgl_smoothness_update(this, x)
{
  h_set, this, state = 0, dx = [];
}

func _rgl_smoothness_set_dx(this, x)
{
  extern dx;
  if (this.state < 1) {
    dx = x - smooth3(x);
    h_set, this, state = 1, dx = dx;
  } else {
    eq_nocopy, dx, this.dx;
  }
}

func _rgl_smoothness_get_penalty(this, x)
{
  local dx;
  _rgl_smoothness_set_dx, this, x;
  return this.mu*sum(dx*dx);
}

func _rgl_smoothness_get_gradient(this, x)
{
  local dx;
  _rgl_smoothness_set_dx, this, x;
  return (2.0*this.mu)*(dx - smooth3(dx));
}

func _rgl_smoothness_apply_hessian(this, x, s)
{
  s -= smooth3(s);
  return (2.0*this.mu)*(s - smooth3(s));
}

#if 0
func _rgl_smoothness_get_hessian(this, x) { }
func _rgl_smoothness_get_diagonal_of_hessian(this, x) { }
#endif

func _rgl_smoothness_get_attr(this, index)
{
  if (index == 1) return this.mu;
}

func _rgl_smoothness_set_attr(this, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), w = [], state = 0;
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "smoothness", "mu";

/*---------------------------------------------------------------------------*/
/* WRAPPER FOR RGL_ROUGHNESS */

local rgl_roughness;
/* DOCUMENT obj = rgl_new("roughness");
 *
 *   Creates a regularization instance suitable for minimizing the roughness
 *   of the parameters by various cost functions. The regularization penalty
 *   writes:
 *
 *     fprior(x) = cost([mu, threshold], D.x)
 *
 *   where D is some finite differences operator.
 *
 *   This regularization has the following hyper-parameters:
 *
 *     1 - "mu"        = global regularization weight.
 *     2 - "threshold" = threshold for L2-L1 or L2-L0 norms.
 *     3 - "cost"      = name of cost function: "l1", "l2", "l2l1", "l2l0",
 *                       or "cauchy"; default is "l2" (i.e. quadratic).
 *     4 - "periodic"  = true for periodic roughness.
 *
 * SEE ALSO: rgl_new, rgl_roughness_l2.
 */

func _rgl_roughness_setup(this)
{
  this = h_set(this, state=0, threshold=1.0);
  _rgl_roughness_set_function, this, "l2", 0;
  return this;
}

func _rgl_roughness_update(this, x)
{
  h_set, this, fx=, gx=, state=0;
}

func _rgl_roughness_compute(this, x)
{
  local gx;
  f = value_of_symlink(this.f);
  fx = (f(this.hyper1,   1,     x, gx) +
        f(this.hyper1, [ 0, 1], x, gx) +
        f(this.hyper2, [-1, 1], x, gx) +
        f(this.hyper2, [ 1, 1], x, gx));
  h_set, this, fx=fx, gx=gx, state=1;
}

func _rgl_roughness_get_penalty(this, x)
{
  if (this.state < 1) {
    _rgl_roughness_compute, this, x;
  }
  return this.fx;
}

func _rgl_roughness_get_gradient(this, x)
{
  if (this.state < 1) {
    _rgl_roughness_compute, this, x;
  }
  return this.gx;
}

func _rgl_roughness_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.threshold;
  if (index == 3) return this.cost;
  if (index == 4) return this.periodic;
}

func _rgl_roughness_set_attr(this, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), state = 0;
      _rgl_roughness_set_hyper, this;
    }
  } else if (index == 2) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "threshold must be a strictly positive real";
    }
    if (value != this.threshold) {
      h_set, this, threshold = double(value), state = 0;
      _rgl_roughness_set_hyper, this;
    }
  } else if (index == 3) {
    if (rgl_string_scalar(value)) {
      error, "cost must be a string";
    }
    if (_rgl_roughness_set_function(this, value, this.periodic)) {
      error, "bad cost function name";
    }
  } else if (index == 4) {
    if (rgl_integer_scalar(value)) {
      error, "periodic must be an integer";
    }
    if (_rgl_roughness_set_function(this, this.cost, value)) {
      error, "bad cost/periodic function";
    }
  }
}

func _rgl_roughness_set_function(this, cost, periodic)
{
  periodic = (periodic ? 1n : 0n);
  name = swrite(format="rgl_roughness_%s%s", cost,
                (periodic ? "_periodic" : ""));
  if (! symbol_exists(name) || ! is_func(symbol_def(name))) {
    return -1n;
  }
  if (cost != this.cost || periodic != this.periodic) {
    if (cost == "l2" || cost == "l1") {
      nhyps = 1;
    } else {
      nhyps = 2;
    }
    h_set, this, f = symlink_to_name(name),
      cost = cost, periodic = periodic, state = 0;
    _rgl_roughness_set_hyper, this;
  }
  return 0n;
}

func _rgl_roughness_set_hyper(this)
{
  if (this.cost == "l2" || this.cost == "l1") {
    hyper1 =     this.mu;
    hyper2 = 0.5*this.mu;
  } else {
    hyper1 = [    this.mu, this.threshold];
    hyper2 = [0.5*this.mu, this.threshold];
  }
  h_set, this, hyper1=hyper1, hyper2=hyper2, state=0;
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "roughness", ["mu", "threshold", "cost", "periodic"];


/*---------------------------------------------------------------------------*/
/* L2-L1 SMOOTHNESS */

local rgl_l2l1_smoothness;
/* DOCUMENT obj = rgl_new("l2l1_smoothness");
 *
 *   Creates a regularization instance suitable for l2-l1 smoothness
 *   regularization.
 *
 *   The regularization penalty writes:
 *
 *     fprior(x) = cost_l2l1([mu, threshold], x - smooth3(x))
 *
 *   This regularization has the following hyper-parameters:
 *
 *     1 - "mu"        = global regularization weight.
 *     2 - "threshold" = threshold for L2-L1.
 *
 * SEE ALSO: rgl_new, cost_l2l1, smooth3.
 */

func _rgl_l2l1_smoothness_setup(this)
{
  return h_set(this, state = 0, mu = 1.0, threshold = 1.0);
}

func _rgl_l2l1_smoothness_update(this, x)
{
  local g;
  r = x - smooth3(x);
  f = cost_l2l1([this.mu, this.threshold], r, g);
  g = g - smooth3(g);
  h_set, this, state = 1, f = f, g = g;
}

func _rgl_l2l1_smoothness_get_penalty(this, x)
{
  if (this.state < 1) {
    _rgl_l2l1_smoothness_update, this, x;
  }
  return this.f;
}

func _rgl_l2l1_smoothness_get_gradient(this, x)
{
  if (this.state < 1) {
    _rgl_l2l1_smoothness_update, this, x;
  }
  return this.g;
}

#if 0
func _rgl_l2l1_smoothness_apply_hessian(this, x, s) {}
func _rgl_l2l1_smoothness_get_hessian(this, x) {}
func _rgl_l2l1_smoothness_get_diagonal_of_hessian(this, x) {}
#endif

func _rgl_l2l1_smoothness_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.threshold;
}

func _rgl_l2l1_smoothness_set_attr(this, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), w = [], state = 0;
    }
  } else if (index == 2) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "l2-l1 threshold must be a strictly positive real";
    }
    if (value != this.threshold) {
      h_set, this, threshold = double(value), w = [], state = 0;
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "l2l1_smoothness", ["mu", "threshold"];

/*---------------------------------------------------------------------------*/
/* GENERAL QUADRATIC REGULARIZATION */

local rgl_quadratic;
/* DOCUMENT obj = rgl_new("quadratic");
 *
 *   Create a regularizer instance for general quadratic regularization.
 *
 *   The general expressions of the penalty and its partial derivatives are
 *   (the dot and the prime indicate dot product and matrix transposition
 *   respectively):
 *
 *       f(x) = mu [A.x - b]'.W.[A.x - b]        (penalty)
 *
 *       g(x) = 2 mu A'.W.[A.x - b]              (gradient)
 *
 *       G    = 2 mu A'.W.A                      (Hessian)
 *
 *   where:
 *
 *       mu = 1-st hyper-parameter
 *       A  = 2-nd hyper-parameter (default: identity)
 *       b  = 3-rd hyper-parameter (default: 0)
 *       W  = 4-th hyper-parameter (default: identity)
 *
 *   If specifed, b must be an array conformable with x.
 *
 *   If specified, A can be an array (diagonal weighting matrix) or a sparse
 *   matrix (see sparse_matrix) or a linear operator object (see linop_new).
 *   The default is the identity.
 *
 *   If specified, W can be an array (diagonal weighting matrix) or a sparse
 *   matrix (see sparse_matrix) or a linear operator object (see linop_new);
 *   whatever is W, it must correspond to a symmetric positive semi-definite
 *   matrix operation.  The default is the identity.
 *
 *   Note that W and A are memorized as linear operator objects (see
 *   linop_new) which must not be forgotten if you use `rgl_config` to query
 *   them.
 *
 *   As a convenience, it is possible to query the hyper-parameters by their
 *   names.
 *
 *
 * SEE ALSO: rgl_info, rgl_config, rgl_identity, sparse_matrix, linop_new.
 */

func _rgl_quadratic_setup(this)
{
  return h_set(this, state = 0, A = rgl_identity, b = [], W = rgl_identity);
}

func _rgl_quadratic_update(this, x)
{
  /* compute the 'residuals': r = A.x - b */
  if (is_void(this.b)) {
    r = this.A(x);
  } else {
    r = this.A(x) - this.b;
  }

  /* compute: 2 mu W.(A.x - b) = 2 mu W.r */
  if ((two_mu = 2.0*this.mu) == 1.0) {
    two_mu_W_r = this.W(r);
  } else {
    two_mu_W_r = two_mu*this.W(r);
  }

  /* update internals (with precomputed penalty and gradient) */
  h_set, this, state = 1,
    f = 0.5*sum(two_mu_W_r*r), /* penalty */
    g = this.A(two_mu_W_r, 1); /* gradient */
}

func _rgl_quadratic_get_penalty(this, x)
{
  if (this.state < 1) _rgl_quadratic_update, this, x;
  return this.f;
}

func _rgl_quadratic_get_gradient(this, x)
{
  if (this.state < 1) _rgl_quadratic_update, this, x;
  return this.g;
}

func _rgl_quadratic_apply_hessian(this, x, s)
{
  /* Result is: 2 mu A'.W.A.s */
  return (2.0*this.mu)*this.A(this.W(this.A(s)), 1);
}

#if 0
func _rgl_quadratic_get_hessian(this, x) {}
func _rgl_quadratic_get_diagonal_of_hessian(this, x) {}
#endif

func _rgl_quadratic_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.A;
  if (index == 3) return this.b;
  if (index == 4) return this.W;
}

func _rgl_quadratic_set_attr(this, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), state = 0;
    }
  } else if (index == 2) {
    h_set, this, A = linop_new(value), state = 0;
  } else if (index == 3) {
    if (is_numerical(value) || is_void(value)) {
      if (noneof(value)) {
        value = [];
      }
      h_set, this, b = value, state = 0;
    } else {
      error, "unexpected value for hyper-parameter \"b\"";
    }
  } else if (index == 4) {
    h_set, this, W = linop_new(value), state = 0;
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "quadratic", ["mu", "A", "b", "W"];

/*---------------------------------------------------------------------------*/
/* GENERAL PURPOSE SMOOTHNESS */

func rgl_xsmooth_residuals(this, x)
/* DOCUMENT rgl_xsmooth_residuals(this, x)
 *
 *   Retuns sample of residual values for unknown X.  THIS must be an
 *   instance of a "xsmooth" regularization.
 *
 * SEE ALSO: rgl_new.
 */
{
  if (is_void(this.a)) {
    _rgl_xsmooth_make_matrix, this;
  }
  xdif = this.a(x);
  w = where(xdif);
  return xdif(w);
}

func _rgl_xsmooth_make_matrix(this)
{
  /* build the finite difference matrix */
  mkmx = rgl_make_2d_finite_difference_matrix; /* shortcut */
  if (this.isotropic) {
    h_set, this, a = mkmx(this.dimlist, [1, 2],
                          1.0,       [+1,  0],
                          1.0,       [ 0, +1],
                          sqrt(2.0), [+1, +1],
                          sqrt(2.0), [+1, -1]);
  } else {
    h_set, this, a = mkmx(this.dimlist, [1, 2],
                          1.0, [+1,  0],
                          1.0, [ 0, +1]);
  }
}

func _rgl_xsmooth_setup(this)
{
  return h_set(this,
               state = 0,
               mu = 1.0,
               threshold = 1.0,
               cost = symlink_to_name("cost_l2"),
               isotropic = 1n);
}

func _rgl_xsmooth_update(this, x)
{
  local g;
  if (is_void(this.a)) {
    _rgl_xsmooth_make_matrix, this;
  }
  f = this.cost([this.mu, this.threshold], this.a(x), g);
  h_set, this, f = f, g = g, state = 1;
}

func _rgl_xsmooth_get_penalty(this, x)
{
  if (this.state < 1) {
    _rgl_xsmooth_update, this, x;
  }
  return this.f;
}

func _rgl_xsmooth_get_gradient(this, x)
{
  if (this.state < 1) {
    _rgl_xsmooth_update, this, x;
  }
  if (this.state < 2) {
    h_set, this, g = this.a(this.g, 1), state = 2;
  }
  return this.g;
}


#if 0
func _rgl_xsmooth_apply_hessian(this, x, s) {}
func _rgl_xsmooth_get_hessian(this, x) {}
func _rgl_xsmooth_get_diagonal_of_hessian(this, x) {}
#endif

func _rgl_xsmooth_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.threshold;
  if (index == 3) return name_of_symlink(this.cost);
  if (index == 4) return this.dimlist; /* FIXME: force new copy? */
  if (index == 5) return this.isotropic;
}

func _rgl_xsmooth_set_attr(this, index, value)
{
  if (index == 1) {
    /* 1 is MU */
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), state = 0;
    }
  } else if (index == 2) {
    /* 2 is THRESHOLD */
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "value of \"threshold\" must be a non-negative real";
    }
    if (value != this.threshold) {
      h_set, this, threshold = double(value), state = 0;
    }
  } else if (index == 3) {
    /* 3 is COST */
    if (! is_string(value) || ! is_scalar(value)) {
      error, "expecting scalar string for attribute \"cost\"";
    }
    if (value != name_of_symlink(this.cost)) {
      h_set, this, cost = symlink_to_name(value), state = 0;
    }
  } else if (index == 4) {
    /* 4 is DIMLIST */
    if (! rgl_check_dimlist(value)) {
      error, "invalid value for attribute \"dimlist\"";
    }
    if (numberof(value) != numberof(this.dimlist) ||
        anyof(value != this.dimlist)) {
      h_set, this, dimlist = value, a = [], state = 0;
    }
  } else if (index == 5) {
    /* 5 is ISOTROPIC */
    if (is_void(value)) {
      value = 0n;
    } else if (is_scalar(value)) {
      value = !(!value);
    } else {
      error, "invalid value for attribute \"isotropic\"";
    }
    if (value != this.isotropic) {
      h_set, this, isotropic = value, a = [], state = 0;
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "xsmooth", ["mu", "threshold", "cost", "dimlist", "isotropic"];

func rgl_make_2d_finite_difference_matrix(dimlist, which, ..)
/* DOCUMENT
 *
 *   rgl_make_2d_finite_difference_matrix(dimlist, which,
 *                                        scale, offset, ...)
 *
 *   Return sparse matrix which compute 2-D finite difference of its
 *   argument.
 *
 *   a(..., i1 + 1, ..., i2, ...) - a(..., i1, ..., i2, ...)
 *
 *   For instance to compute isotropic finite differences between the 2-nd
 *   and 4-th dimensions:
 *
 *   rgl_make_2d_finite_difference_matrix(dimlist, [2, 4],
 *                                        1.0,       [+1,  0],
 *                                        1.0,       [ 0, +1],
 *                                        sqrt(2.0), [+1, +1],
 *                                        sqrt(2.0), [+1, -1]);
 */
{
  /* Check/fix dimension list. */
  make_dimlist, dimlist;
  ndims = dimlist(1);
  if (ndims < 2) {
    error, "expecting at least 2-D dimension list";
  }
  if (! is_integer(which) || numberof(which) != 2) {
    error, "WHICH must be a 2-element vector of integers";
  }
  first = which(1);
  if (first <= 0) first += ndims;
  second = which(2);
  if (second <= 0) second += ndims;
  if (first <= 0 || first > ndims || second == first ||
      second <= 0 || second > ndims) {
    error, "invalid indices for dimensions of interest";
  }
  exchange = (first > second);
  if (exchange) {
    swap, first, second;
  }

  /* Build lists of scales and offsets from the argument list. */
  local scale, scale_list, offset, offset_list;
  while (more_args()) {
    eq_nocopy, scale, next_arg();
    if (! more_args()) {
      error, "missing offset";
    }
    eq_nocopy, offset, next_arg();
    if (! is_real(scale) || min(scale) <= 0.0) {
      error, "scales must be strictly non-negative reals";
    }
    if (! is_integer(offset)) {
      error, "offsets must be integers";
    }
    dwgh = dimsof(scale);
    doff = dimsof(offset);
    if (! dwgh(1) && doff(1) == 1 && doff(2) == 2) {
      grow, scale_list, double(scale);
      grow, offset_list, [long(offset)];
    } else if (dwgh(1) == 1 && doff(1) == 2 &&
               doff(2) == 2 && doff(3) == dwgh(2)) {
      grow, scale_list, double(scale);
      grow, offset_list, long(offset);
    } else {
      error, "bad or incompatible dimensions for scales and/or offsets";
    }
  }
  if (exchange) {
    offset_list = offset_list(::-1,);
  }

  /* Compute strides and total number of elements. */
  stride = array(long, ndims);
  number = 1;
  for (k = 1; k <= ndims; ++k) {
    stride(k) = number;
    number *= dimlist(k + 1);
  }

  /* The idea is to map: a(.., i1, .., i2, ..)
   * as:                 a(j0, j1, j2, j3, j4)
   * where:
   *   j0 = 1:STRIDE1
   *   j1 = i1 = 1:NUMBER1
   *   j2 = 1:(STRIDE2/NUMBER1/STRIDE1)
   *   j3 = i2 = 1:NUMBER2
   *   j4 = 1:(NUMBER/NUMBER2/STRIDE2)
   */

  /* Strides and lengths for dimensions of interest. */
  s1 = stride(first);
  n1 = dimlist(first + 1);
  s3 = stride(second);
  n3 = dimlist(second + 1);

  /* Strides and lengths for interleaving dimensions. */
  s0 = 1;
  n0 = s1;
  s2 = s1*n1;
  n2 = s3/s2;
  s4 = s3*n3;
  n4 = number/s4;

  if (n4 > 1) {
    j4 = indgen(0 : s4*(n4 - 1) : s4);
  }
  j3 = indgen(0 : s3*(n3 - 1) : s3);
  if (n2 > 1) {
    j2 = indgen(0 : s2*(n2 - 1) : s2);
  }
  j1 = indgen(0 : s1*(n1 - 1) : s1);
  if (n0 > 1) {
    /* J0 is 1-based offset */
    j0 = indgen(1 : 1 + s0*(n0 - 1) : s0);
  }

  /* Build sparse matrix. */
  local coef_list, index_list;
  count = 0;
  ndifs = numberof(scale_list);
  for (d = 1; d <= ndifs; ++d) {
    factor = 1.0/scale_list(d);
    offset = offset_list(, d);
    if (noneof(offset)) {
      continue;
    }
    k1 = j1 + s1*offset_list(1, d);
    k3 = j3 + s3*offset_list(2, d);
    w1 = where((k1 >= 0)&(k1 < n1*s1));
    w3 = where((k3 >= 0)&(k3 < n3*s3));
    if (! is_array(w1) || ! is_array(w3)) {
      continue;
    }
    j = j1(w1);
    k = k1(w1);
    if (n4 > 1) {
      j += j4(-,);
      k += j4(-,);
    }
    if (n2 > 1) {
      j = j2 + j(-,..);
      k = j2 + k(-,..);
    }
    j = j3(w3) + j(-,..);
    k = k3(w3) + k(-,..);
    if (n0 > 1) {
      j = j0 + j(-,..);
      k = j0 + k(-,..);
    } else {
      ++j;
      ++k;
    }

    n = numberof(k);
    tmp = array(double, 2, n);
    tmp(1,) =  factor;
    tmp(2,) = -factor;
    grow, coef_list, unref(tmp);
    tmp = array(long, 2, n);
    tmp(1,) = k(*);
    tmp(2,) = j(*);
    grow, index_list, unref(tmp);
    count += n; /* number of finite differences */
  }

  return sparse_matrix(coef_list, [1, count], indgen(count)(-:1:2,),
                       dimlist, index_list);
}

/*---------------------------------------------------------------------------*/
/* (ISOTROPIC) TOTAL VARIATION (2D) */

local rgl_totvar;
/* DOCUMENT this = rgl_new("totvar");
 *
 *   Creates a regularizer instance suitable for "Total Variation" (TV)
 *   regularization.  The regularization penalty writes:
 *
 *       mu*sum(dx)
 *
 *   where DX is the length of the local gradient of X along its dimensions:
 *
 *      DX = sqrt(DX1^2 + DX2^2 + ... + EPSILON^2)
 *
 *   where DXn is the partial derivative of X along n-th dimension.
 *
 *   This regularization has the following hyper-parameters:
 *
 *     1 - "mu" = global regularization weight;
 *     2 - "epsilon" = small value to get rid of singularities;
 *     3 - "isotropic" = flag: use isotropic "Total Variation".
 *     4 - "mask" = mask, an array of 0/1 where pixels are to be
 *                  ignored/considered.
 *
 *   Note that by using a small but not negligible EPSILON, edge-preserving
 *   smoothness is achieved.
 *
 *
 * SEE ALSO: rgl_new, smooth3.
 */

func _rgl_totvar_setup(self)
{
  return h_set(self, state=0, epsilon=1e-8, isotropic=1n);
}

func _rgl_totvar_update(self, x)
{
  h_set, self, state=0;
}

func _rgl_totvar_state1(self, x)
{
  local mask; eq_nocopy, mask, self.mask;
  w = self.mu/sqrt(2.0);
  eps = abs(self.epsilon);
  i0 = 1:-1;
  i1 = 2:0;
  if (self.isotropic) {
    x00 = x(i0,i0);
    x10 = x(i1,i0);
    x01 = x(i0,i1);
    x11 = x(i1,i1);
    d1 = (x11 - x00);
    d2 = (x10 - x01);
    d3 = (unref(x00) - unref(x01) - unref(x10) + unref(x11));
    if (! is_void(mask)) {
      m1 = mask(i1,i1)*mask(i0,i0);
      m2 = mask(i1,i0)*mask(i0,i1);
      d1 = m1*unref(d1);
      d2 = m2*unref(d2);
      d3 = unref(m1)*unref(m2)*unref(d3);
    }
    r = sqrt(d1*d1 + d2*d2 + (1.0/3.0)*d3*d3 + 2.0*eps*eps);
    err = w*(sum(r) - numberof(r)*eps);
    h_set, self, state=1, err=err, w=w, r=r, d1=d1, d2=d2, d3=d3;
  } else {
    if (is_void(mask)) {
      d1 = (x(i1,i1) - x(i0,i0));
      d2 = (x(i1,i0) - x(i0,i1));
    } else {
      d1 = (x(i1,i1) - x(i0,i0))*mask(i1,i1)*mask(i0,i0);
      d2 = (x(i1,i0) - x(i0,i1))*mask(i1,i0)*mask(i0,i1);
    }
    r = sqrt(d1*d1 + d2*d2 + 2.0*eps*eps);
    err = w*(sum(r) - numberof(r)*eps);
    h_set, self, state=1, err=err, w=w, r=r, d1=d1, d2=d2;
  }
}

func _rgl_totvar_state2(self, x)
{
  if (self.state < 1) _rgl_totvar_state1, self, x;
  g = array(double, dimsof(x));
  q = self.w/self.r;
  r0 = 1:-1;
  r1 = 2:0;
  d1 = self.d1*q;
  d2 = self.d2*q;
  if (self.isotropic) {
    d3 = (1.0/3.0)*self.d3*q;
    g(r0, r0) -= d1 - d3;
    g(r1, r0) += d2 - d3;
    g(r0, r1) -= d2 + d3;
    g(r1, r1) += d1 + d3;
  } else {
    g(r0, r0) -= d1;
    g(r1, r0) += d2;
    g(r0, r1) -= d2;
    g(r1, r1) += d1;
  }
  h_set, self, state=2, grd=g;
}

func _rgl_totvar_get_penalty(self, x)
{
  if (self.state < 1) _rgl_totvar_state1, self, x;
  return self.err;
}

func _rgl_totvar_get_gradient(self, x)
{
  if (self.state < 2) _rgl_totvar_state2, self, x;
  return self.grd;
}

#if 0
func _rgl_totvar_apply_hessian(self, x, s) {}
func _rgl_totvar_get_hessian(self, x) {}
func _rgl_totvar_get_diagonal_of_hessian(self, x) {}
#endif

func _rgl_totvar_get_attr(self, index)
{
  if (index == 1) return self.mu;
  if (index == 2) return self.epsilon;
  if (index == 3) return self.isotropic;
  if (index == 4) return self.mask;
}

func _rgl_totvar_set_attr(self, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "\"global\" regularization weight must be a non-negative real";
    }
    if (value != self.mu) {
      h_set, self, mu = double(value), state = 0;
    }
  }
  if (index == 2) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "\"epsilon\" must be a strictly non-negative real";
    }
    if (value != self.epsilon) {
      h_set, self, epsilon = double(value), state = 0;
    }
  }
  if (index == 3) {
    if (rgl_boolean(value)) {
      error, "\"isotropic\" must be a boolean value";
    }
    if (value != self.isotropic) {
      h_set, self, isotropic = value, state = 0;
    }
  }
  if (index == 4) {
    if (is_void(value)) {
      h_set, self, mask = [], state = 0;
    } else if (numberof(dimsof(value)) == 3 && identof(value) <= Y_DOUBLE) {
      zero = structof(value)(0);
      h_set, self, mask = double(value != zero), state = 0;
    } else {
      error, "\"mask\" must be a 2-D boolean array";
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "totvar", ["mu", "epsilon", "isotropic", "mask"];

/*---------------------------------------------------------------------------*/
/* SMOOTHNESS W.R.T. AVERAGE ESTIMATE */

func _rgl_qsmooth_setup(this)
{
  return h_set(this,
               state = 0,
               mu = 1.0,
               prior = []);
}

func _rgl_qsmooth_update(this, x)
{
  h_set, this, state = 0;
}

func _rgl_qsmooth_get_penalty(this, x)
{
  if (this.state < 1) {
    if (! is_void(this.prior)) x -= this.prior;
    cut = 1:-1;
    dx1 = x(dif, cut);
    dx2 = x(cut, dif);
    f = this.mu*(sum(dx1*dx1) + sum(dx2*dx2));
    h_set, this, dx1 = dx1, dx2 = dx2, f = f, q = q, state = 1;
  }
  return this.f;
}

func _rgl_qsmooth_get_gradient(this, x)
{
  if (this.state < 2) {
    if (this.state < 1) _rgl_qsmooth_get_penalty, this, x;
    g = array(double, dimsof(x));
    lo = 1:-1;
    hi = 2:0;
    h_set, this, state = 0; // in case of interrupts...
    temp = h_pop(this, dx1=);
    g(hi, lo)  = temp;
    g(lo, lo) -= temp;
    temp = h_pop(this, dx2=);
    g(lo, hi) += temp;
    g(lo, lo) -= temp;
    h_set, this, g = (2.0*this.mu)*g, state = 2;
  }
  return this.g;
}

#if 0
func _rgl_qsmooth_apply_hessian(this, x, s) {}
func _rgl_qsmooth_get_hessian(this, x) {}
func _rgl_qsmooth_get_diagonal_of_hessian(this, x) {}
#endif

func _rgl_qsmooth_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.threshold;
}

func _rgl_qsmooth_set_attr(this, index, value)
{
  if (index == 1) {
    /* 1 is MU */
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), state = 0;
    }
  } else if (index == 2) {
    /* 2 is PRIOR */
    h_set, this, prior = double(value), state = 0;
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "qsmooth", ["mu", "prior"];

/*---------------------------------------------------------------------------*/
/* PREDEFINED QUADRATIC REGULARIZATIONS */

func rgl_new_quadratic_1(.., fwhm=, shape=, normalization=, method=)
{
  /* Get dimension list. */
  dimlist = [0L]; // initial dimension list
  while (more_args()) {
    make_dimlist, dimlist, next_arg();
  }
  if (! dimlist(1)) {
    error, "you must specify a dimension list";
  }

  if (is_void(fwhm)) {
    fwhm = (1.0/3.0)*dimlist(2:);
  }
  //write, fwhm;
  if (shape == "gauss" || shape == "gaussian") {
    shape = 1;
    scale = sqrt(log(16))/fwhm;
  } else if (is_void(shape) || shape == "lorentz" || shape == "lorentzian") {
    shape = 2;
    scale = 2.0/fwhm;
  } else {
    error, "bad value for keyword SHAPE";
  }
  r2 = rgl_squared_distance(dimlist, scale=scale);
  if (shape == 1) {
    prior = exp(-r2);
  } else {
    prior = 1.0/(1.0 + r2);
  }

  rgl = rgl_new("quadratic", dimlist);

  if (normalization) {
    prior *= (double(normalization)/sum(prior));
  }
  if (method == 2) {
    if (min(prior) <= 0.0) {
      error, "try Lorentzian shape instead";
    }
    rgl_config, rgl, "W", linop_new("diagonal", 1.0/prior);
  } else {
    rgl_config, rgl, "b", prior;
  }
  return rgl;
}

/*---------------------------------------------------------------------------*/
/* MAXIMUM ENTROPY REGULARIZATION */

local rgl_entropy;
/* DOCUMENT obj = rgl_new("entropy", ...)
 *
 *  Beware that if you use regularization based on entropy, you must insure
 *  that min(x) > 0 (strictly positive).  For instance:
 *
 *      XMIN = EPSILON/avg(X)   with EPSILON a small positive number
 *
 *  Possible definitions for neg-entropy regularization (and partial
 *  derivatives):
 *
 *    f1(x) = -mu*sum(sqrt(x + eps));
 *    g1(x) = -0.5*mu/sqrt(x + eps);
 *    h1(x) = diag(0.25*mu/((x + eps)*sqrt(x + eps));
 *
 *    f2(x) = -mu*sum(log(x + eps));
 *    g2(x) = -mu/(x + eps);
 *    h2(x) = diag(mu/(x + eps)^2);
 *
 *  Maximum entropy (no prior, X is normalized):
 *    f3(x) = mu*sum((x + eps)*log(x + eps));
 *    g3(x) = mu + mu*log(x + eps);
 *    h3(x) = diag(mu/(x + eps));
 *
 *  Maximum entropy (P is the prior, X is not normalized):
 *    f4(x) = mu*sum(p - x + (x + eps)*log((x + eps)/(p + eps)));
 *    g4(x) = mu*log((x + eps)/(p + eps));
 *    h4(x) = diag(mu/(x + eps));
 *
 *  Maximum entropy (P is the prior, X is normalized):
 *    f5(x) = mu*sum((x + eps)*log((x + eps)/(p + eps)));
 *    g5(x) = mu + mu*log((x + eps)/(p + eps));
 *    h5(x) = diag(mu/(x + eps));
 *
 *  Maximum entropy (prior linearly depends on X: p = A.x,
 *  A is any linear operator, X is not normalized):
 *    f6(x) = mu*sum(A.x - x + x*log(x/A.x));
 *    g6(x) = mu*(log(x/A.x) + A'.(1 - x/A.x));
 *    h6(x) = mu*(diag(1/x) + 2*A'.diag(x/A.x).A
 *                - A[j,k]/(A.x)[j] - A[k,j]/(A.x)[k]);
 *
 *  Maximum entropy (prior linearly depends on X: p = A.x,
 *  A is normalized such that sum(A.x) = sum(x), X is not normalized):
 *    f7(x) = mu*sum(x*log(x/A.x));
 *    g7(x) = mu*(log(x/A.x) - A'.(x/A.x));
 *
 *
 *  Attributes are:
 *    "mu"         (1) - global weight.
 *    "type"       (2) - type: "sqrt", or "log".
 *    "normalized" (3) - X is normalized?
 *    "prior"      (4) - prior p (if an array) or matrix A (if a linear
 *                       operator), or "none".
 *    "epsilon"    (5) - small value to get rid of singularities near zero,
 *                       denoted EPS in equations above (default is 1E-20).
 *
 *  Cases:
 *     Id. Type    Normalized   Prior   Negentropy
 *     -------------------------------------------------------------
 *     1   "sqrt"  0            nil     -mu*sum(sqrt(x))
 *     2   "log"   0            nil     -mu*sum(log(x))
 *     3   "log"   1            nil      mu*sum(x*log(x))
 *     4   "log"   0            p        mu*sum(p - x + x*log(x/p))
 *     5   "log"   1            p        mu*sum(x*log(x/p))
 *     6   "log"   0            A        mu*sum(A.x - x + x*log(x/A.x))
 *     7   "log"   1            A        mu*sum(x*log(x/A.x))
 *
 * SEE ALSO: rgl_new.
 */

func _rgl_entropy_setup(this)
{
  /* state = 0 must call "finalize" to figure out which method to use
   * state = 1 must call "update"
   * state = 2 ok
   */
  return h_set(this, state = 0,
               type = "sqrt", normalized = 0n,
               prior_type = 0, prior = [], epsilon = 1E-20);
}

func _rgl_entropy_finalize(this)
{
  local type; eq_nocopy, type, this.type;
  local normalized; eq_nocopy, normalized, this.normalized;
  local prior; eq_nocopy, prior, this.prior;
  if (this.type == "sqrt") {
    id = 1;
    if (this.prior_type != 0) {
      write, format="WARNING: %s\n", "regularization prior ignored.";
    }
    if (this.normalized) {
      write, format="WARNING: %s\n", "normalization flag ignored.";
    }
  } else /* must be "log" */ {
    if (this.prior_type == 0) {
      id = (this.normalized ? 3 : 2);
    } else if (this.prior_type == 1) {
      id = (this.normalized ? 5 : 4);
    } else /* must be: prior_type = 2 */ {
      id = (this.normalized ? 7 : 6);
    }
  }
  ops = symlink_to_name(swrite(format="_rgl_entropy%d_ops", id));
  h_set, this, id = id, ops = ops, state = 1;
  h_delete, this, "log_x", "log_x_over_p",
    "Ax", "x_over_Ax", "log_x_over_Ax";
}

func _rgl_entropy_update(this, x)
{
  if (this.state < 1) _rgl_entropy_finalize, this;
  op = this.ops.update;
  op, this, x;
  h_set, this, state = 2;
}

func _rgl_entropy_get_penalty(this, x)
{
  if (this.state < 2) _rgl_entropy_update, this, x;
  return this.ops.get_penalty(this, x);
}

func _rgl_entropy_get_gradient(this, x)
{
  if (this.state < 2) _rgl_entropy_update, this, x;
  return this.ops.get_gradient(this, x);
}

func _rgl_entropy_apply_hessian(this, x, s)
{
  if (this.state < 2) _rgl_entropy_update, this, x;
  return this.ops.apply_hessian(this, x, s);
}

func _rgl_entropy_get_diagonal_of_hessian(this, x)
{
  if (this.state < 2) _rgl_entropy_update, this, x;
  return this.ops.get_diagonal_of_hessian(this, x);
}

func _rgl_entropy_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.type;
  if (index == 3) return this.normalized;
  if (index == 4) return this.prior;
  if (index == 5) return this.epsilon;
}

func _rgl_entropy_set_attr(this, index, value)
{
  if (index == 1) {
    /* 1 is MU */
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value);
    }
  } else if (index == 2) {
    /* 2 is TYPE */
    if (rgl_string_scalar(value) || (value != "log" && value != "sqrt")) {
      error, "entropy type must be \"log\" or \"sqrt\"";
    }
    if (value != this.type) {
      h_set, this, type = (value + ""), state = 0;
    }
  } else if (index == 3) {
    /* 3 is NORMALIZED */
    if (! is_scalar(value) || ! is_numerical(value)) {
      error, "attribute \"normalized\" must be a numerical scalar";
    }
    value = !(!value);
    if (value != this.normalized) {
      h_set, this, normalized = value, state = 0;
    }
  } else if (index == 4) {
    /* 4 is PRIOR */
    if (is_hash(value)) {
      /* assume linear operator (FIXME: there should be a is_linop function) */
      prior_type = 2; /* linear operator */
      change = (value != this.prior);
    } else if (is_string(value) && is_scalar(value) && value == "none") {
      prior_type = 0; /* void */
      value = [];
      change = ! is_void(this.prior);
    } else if (is_real(value)) {
      prior_type = 1; /* array */
      if (is_array(this.prior)) {
        a = dimsof(value);
        b = dimsof(this.prior);
        change = (numberof(a) != numberof(b) || anyof(a != b) ||
                  anyof(value != this.prior));
      } else {
        change = 1n;
      }
      if (change) {
        value = double(value);
      }
    } else {
      error, "invalid data type for attribute \"prior\"";
    }
    if (change) {
      h_set, this, prior_type = prior_type, prior = value, state = 0;
    }
  } else if (index == 5) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "\"epsilon\" must be a strictly non-negative real";
    }
    if (value != this.epsilon) {
      h_set, this, epsilon = double(value), state = 0;
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "entropy", ["mu", "type", "normalized", "prior", "epsilon"];

func _rgl_build_ops(update, get_penalty, get_gradient, apply_hessian,
                    get_diagonal_of_hessian)
{
  return h_new(update = symlink_to_name(update),
               get_penalty = symlink_to_name(get_penalty),
               get_gradient = symlink_to_name(get_gradient),
               apply_hessian = symlink_to_name(apply_hessian),
               get_diagonal_of_hessian = symlink_to_name(get_diagonal_of_hessian));
}

func _rgl_entropy1_update(this, x)
{
  nil = [];
  h_set, this, state = 0,
    _sx = nil,
    _fx = nil,
    _gx = nil,
    _hx = nil;
}

func _rgl_entropy1_preamble(this, x)
{
  h_set, this, state = 1, _sx = sqrt(x + this.epsilon);
}

func _rgl_entropy1_get_penalty(this, x)
{
  if ((this.state & 3) != 3) {
    if ((this.state & 1) != 1) _rgl_entropy1_preamble, this, x;
    h_set, this, state = (this.state | 3),
      _fx = -this.mu*sum(this._sx);
  }
  return this._fx;
}

func _rgl_entropy1_get_gradient(this, x)
{
  if ((this.state & 5) != 5) {
    if ((this.state & 1) != 1) _rgl_entropy1_preamble, this, x;
    h_set, this, state = (this.state | 5),
      _gx = (-0.5*this.mu)/this._sx;
  }
  return this._gx;
}

func _rgl_entropy1_apply_hessian(this, x, s)
{
  return _rgl_entropy1_get_diagonal_of_hessian(this, x)*s;
}

func _rgl_entropy1_get_diagonal_of_hessian(this, x)
{
  if ((this.state & 9) != 9) {
    if ((this.state & 1) != 1) _rgl_entropy1_preamble, this, x;
    h_set, this, state = (this.state | 9),
      _hx = (0.25*this.mu)/((x + this.epsilon)*this._sx);
  }
  return this._hx;
}
_rgl_entropy1_ops = _rgl_build_ops("_rgl_entropy1_update",
                                   "_rgl_entropy1_get_penalty",
                                   "_rgl_entropy1_get_gradient",
                                   "_rgl_entropy1_apply_hessian",
                                   "_rgl_entropy1_get_diagonal_of_hessian");

func _rgl_entropy2_update(this, x)
{
}

func _rgl_entropy2_get_penalty(this, x)
{
  return -this.mu*sum(log(x + this.epsilon));
}

func _rgl_entropy2_get_gradient(this, x)
{
  return (-this.mu)/(x + this.epsilon);
}

func _rgl_entropy2_apply_hessian(this, x, s)
{
  x = unref(x) + this.epsilon;
  return (this.mu*s)/(x*x);
}

func _rgl_entropy2_get_diagonal_of_hessian(this, x)
{
  x = unref(x) + this.epsilon;
  return this.mu/(x*x);
}
_rgl_entropy2_ops = _rgl_build_ops("_rgl_entropy2_update",
                                   "_rgl_entropy2_get_penalty",
                                   "_rgl_entropy2_get_gradient",
                                   "_rgl_entropy2_apply_hessian",
                                   "_rgl_entropy2_get_diagonal_of_hessian");

func _rgl_entropy3_update(this, x)
{
  nil = [];
  h_set, this, state = 0,
    _x = nil,
    _lx = nil,
    _fx = nil,
    _gx = nil,
    _hx = nil;
}

func _rgl_entropy3_preamble(this, x)
{
  x = unref(x) + this.epsilon;
  h_set, this, state = 1, _lx = log(x), _x = x;
}

func _rgl_entropy3_get_penalty(this, x)
{
  if ((this.state & 3) != 3) {
    if ((this.state & 1) != 1) _rgl_entropy3_preamble, this, x;
    h_set, this, state = (this.state | 3),
      _fx = this.mu*sum(this._x*this._lx);
  }
  return this._fx;
}

func _rgl_entropy3_get_gradient(this, x)
{
  if ((this.state & 5) != 5) {
    if ((this.state & 1) != 1) _rgl_entropy3_preamble, this, x;
    h_set, this, state = (this.state | 5),
      _gx = this.mu + this.mu*this._lx;
  }
  return this._gx;
}

func _rgl_entropy3_apply_hessian(this, x, s)
{
  return _rgl_entropy3_get_diagonal_of_hessian(this, x)*s;
}

func _rgl_entropy3_get_diagonal_of_hessian(this, x)
{
  if ((this.state & 9) != 9) {
    if ((this.state & 1) != 1) _rgl_entropy3_preamble, this, x;
    h_set, this, state = (this.state | 9),
      _hx = this.mu/this._x;
  }
  return this._hx;
}

_rgl_entropy3_ops = _rgl_build_ops("_rgl_entropy3_update",
                                   "_rgl_entropy3_get_penalty",
                                   "_rgl_entropy3_get_gradient",
                                   "_rgl_entropy3_apply_hessian",
                                   "_rgl_entropy3_get_diagonal_of_hessian");

func _rgl_entropy4_update(this, x)
{
  nil = [];
  h_set, this, state = 0,
    _x = nil,
    _fx = nil,
    _gx = nil,
    _hx = nil;
}

func _rgl_entropy4_preamble(this, x)
{
  x = unref(x) + this.epsilon;
  h_set, this, state = 5,
    _x = x,
    _gx = this.mu*log(x/(this.prior + this.epsilon));
}

func _rgl_entropy4_get_penalty(this, x)
{
  if ((this.state & 3) != 3) {
    if ((this.state & 5) != 5) _rgl_entropy4_preamble, this, x;
    h_set, this, state = (this.state | 3),
      _fx = this.mu*sum(this.prior - x) + sum(this._x*this._gx);
  }
  return this._fx;
}

func _rgl_entropy4_get_gradient(this, x)
{
  if ((this.state & 5) != 5) {
    _rgl_entropy4_preamble, this, x;
  }
  return this._gx;
}

func _rgl_entropy4_apply_hessian(this, x, s)
{
  return _rgl_entropy4_get_diagonal_of_hessian(this, x)*s;
}

func _rgl_entropy4_get_diagonal_of_hessian(this, x)
{
  if ((this.state & 9) != 9) {
    if ((this.state & 1) != 1) _rgl_entropy4_preamble, this, x;
    h_set, this, state = (this.state | 9),
      _hx = this.mu/this._x;
  }
  return this._hx;
}

_rgl_entropy4_ops = _rgl_build_ops("_rgl_entropy4_update",
                                   "_rgl_entropy4_get_penalty",
                                   "_rgl_entropy4_get_gradient",
                                   "_rgl_entropy4_apply_hessian",
                                   "_rgl_entropy4_get_diagonal_of_hessian");

func _rgl_entropy5_update(this, x)
{
  nil = [];
  h_set, this, state = 0,
    _x = nil,
    _lx = nil,
    _fx = nil,
    _gx = nil,
    _hx = nil;
}

func _rgl_entropy5_preamble(this, x)
{
  x = unref(x) + this.epsilon;
  h_set, this, state = 5,
    _x = x,
    _lx = log(x/(this.prior + this.epsilon));
}

func _rgl_entropy5_get_penalty(this, x)
{
  if ((this.state & 3) != 3) {
    if ((this.state & 1) != 1) _rgl_entropy5_preamble, this, x;
    h_set, this, state = (this.state | 3),
      _fx = this.mu*sum(this._x*this._lx);
  }
  return this._fx;
}

func _rgl_entropy5_get_gradient(this, x)
{
  if ((this.state & 5) != 5) {
    if ((this.state & 1) != 1) _rgl_entropy5_preamble, this, x;
    h_set, this, state = (this.state | 5),
      _gx = this.mu*this._lx + this.mu;
  }
  return this._gx;
}

func _rgl_entropy5_apply_hessian(this, x, s)
{
  return _rgl_entropy5_get_diagonal_of_hessian(this, x)*s;
}

func _rgl_entropy5_get_diagonal_of_hessian(this, x)
{
  if ((this.state & 9) != 9) {
    if ((this.state & 1) != 1) _rgl_entropy5_preamble, this, x;
    h_set, this, state = (this.state | 9),
      _hx = this.mu/this._x;
  }
  return this._hx;
}

_rgl_entropy5_ops = _rgl_build_ops("_rgl_entropy5_update",
                                   "_rgl_entropy5_get_penalty",
                                   "_rgl_entropy5_get_gradient",
                                   "_rgl_entropy5_apply_hessian",
                                   "_rgl_entropy5_get_diagonal_of_hessian");

func _rgl_entropy6_update(this, x)
{
  Ax = this.prior(x);
  x_over_Ax = x/Ax;
  h_set, this,
    Ax = Ax,
    x_over_Ax = x_over_Ax,
    log_x_over_Ax = log(x_over_Ax);
}

func _rgl_entropy6_get_penalty(this, x)
{
  return this.mu*sum(this.Ax - x + x*this.log_x_over_Ax);
}

func _rgl_entropy6_get_gradient(this, x)
{
  return this.mu*(this.log_x_over_Ax +
                  this.prior(1.0 - this.x_over_Ax, 1));
}

func _rgl_entropy6_apply_hessian(this, x, s)
{
  error, "method \"apply_hessian\" not yet implemented";
}

func _rgl_entropy6_get_diagonal_of_hessian(this, x)
{
  error, "method \"get_diagonal_of_hessian\" not yet implemented";
}

_rgl_entropy6_ops = _rgl_build_ops("_rgl_entropy6_update",
                                   "_rgl_entropy6_get_penalty",
                                   "_rgl_entropy6_get_gradient",
                                   "_rgl_entropy6_apply_hessian",
                                   "_rgl_entropy6_get_diagonal_of_hessian");

func _rgl_entropy7_update(this, x)
{
  Ax = this.prior(x);
  x_over_Ax = x/Ax;
  h_set, this,
    x_over_Ax = x_over_Ax,
    log_x_over_Ax = log(x_over_Ax);
}

func _rgl_entropy7_get_penalty(this, x)
{
  return this.mu*sum(x*this.log_x_over_Ax);
}

func _rgl_entropy7_get_gradient(this, x)
{
  return this.mu*(this.log_x_over_Ax - this.prior(this.x_over_Ax, 1));
}

func _rgl_entropy7_apply_hessian(this, x, s)
{
  error, "method \"apply_hessian\" not yet implemented";
}

func _rgl_entropy7_get_diagonal_of_hessian(this, x)
{
  error, "method \"get_diagonal_of_hessian\" not yet implemented";
}

_rgl_entropy7_ops = _rgl_build_ops("_rgl_entropy7_update",
                                   "_rgl_entropy7_get_penalty",
                                   "_rgl_entropy7_get_gradient",
                                   "_rgl_entropy7_apply_hessian",
                                   "_rgl_entropy7_get_diagonal_of_hessian");

/*---------------------------------------------------------------------------*/
/* SMOOTHNESS BY REGION */

local rgl_clique;
/* DOCUMENT obj = rgl_new("clique");
 *          rgl_config, obj, "region", region;
 *      or:
 *          obj = rgl_new("clique", "region", region);
 *
 *   Smoothness regularization by "cliques".  The penalty is the quadratic
 *   difference between adjacent pixels which belong to the same region.  The
 *   parameter REGION is an integer valued array with same dimension list as
 *   the image, same values in REGION indicate pixels which belong to the same
 *   clique.
 *
 * SEE ALSO: rgl_info, rgl_config, sparse_matrix, linop_new.
 */

func _rgl_clique_setup(this)
{
  return h_set(this, state = 0);
}

func _rgl_clique_update(this, x)
{
  h_pop, this, "grd";
  if (this.state < 1) {
    _rgl_clique_builder, this;
  }
  r = x(this.i1) - x(this.i0);
  g = (2.0*this.mu)*this.w*r;
  err = 0.5*sum(g*r);
  grd = array(double, dimsof(x));
  grd(*) = (histogram(this.i1, g, top=numberof(x)) -
            histogram(this.i0, g, top=numberof(x)));
  h_set, this, state=2, grd=grd, err=err;
}

func _rgl_clique_get_penalty(this, x)
{
  if (this.state < 2) _rgl_clique_update, this, x;
  return this.err;
}

func _rgl_clique_get_gradient(this, x)
{
  if (this.state < 2) _rgl_clique_update, this, x;
  return this.grd;
}

#if 0
func _rgl_clique_apply_hessian(this, x, s) { }
func _rgl_clique_get_hessian(this, x) { }
func _rgl_clique_get_diagonal_of_hessian(this, x) { }
#endif

func _rgl_clique_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.region;
}

func _rgl_clique_set_attr(this, index, value)
{
  if (index == 1) {
    /* 1 is MU */
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), w = [], state = 0;
    }
  } else if (index == 2) {
    /* 2 is REGION */
    if (! is_array(value) || dimsof(value)(1) != 2 ||
        (id = identof(value)) > Y_LONG) {
      error, "REGION must be a 2-D integer array";
    }
    if (is_void(this.region)
        || is_void(dimsof(value, this.region))
        || anyof(value != this.region)) {
      value = long(value); /* make a private copy */
      h_set, this, region = value, state = 0;
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "clique", ["mu", "region"];


func _rgl_clique_builder(this)
{
  local region, op_i0, op_i1, op_w, op_reg;
  eq_nocopy, region, this.region;
  if (! is_array(region) || (dims = dimsof(region))(1) != 2) {
    error, "REGION must be a 2-D array";
  }
  width = dims(2);
  height = dims(3);
  offsets = [[1,0],[-1,1],[0,1],[1,1]];
  n = numberof(offsets)/2;
  x0 = indgen(width);
  y0 = indgen(height);
  for (k = 1; k <= n; ++k) {
    x1 = x0 + (xoff = offsets(1,k));
    xsel = where((x1 >= 1)&(x1 <= width));
    y1 = y0 + (yoff = offsets(2,k));
    ysel = where((y1 >= 1)&(y1 <= height));
    if (is_array(xsel) && is_array(ysel)) {
      i0 = x0(xsel) + (width*(y0(ysel) - 1))(-,);
      i1 = x1(xsel) + (width*(y1(ysel) - 1))(-,);
      sel = where(region(i0) == region(i1));
      if (is_array(sel)) {
        i0 = i0(sel);
        i1 = i1(sel);
        grow, op_i0, i0;
        grow, op_i1, i1;
        grow, op_reg, region(i0);
        grow, op_w, array(1.0/(xoff*xoff + yoff*yoff), numberof(sel));
      }
    }
  }

  n = numberof(op_i0);
  if (n > 1) {
    j = heapsort(op_i0);
    op_i0 = op_i0(j);
    op_i1 = op_i1(j);
    op_w = op_w(j);
    op_reg = op_reg(j);
  }
  h_set, this, state=1, i0 = op_i0, i1 = op_i1, w = op_w, reg = op_reg;
}

/*---------------------------------------------------------------------------*/
/* SIMPLE L0 REGULARIZATION */

local rgl_simple;
/* DOCUMENT this = rgl_new("simple");
 *
 *   Creates a regularization instance suitable for simple separable L2-L0
 *   regularization.
 *
 *   The regularization penalty for an image X is:
 *
 *      mu * eps^2 * atan(x/eps))^2
*
 *   This regularization has two hyper-parameters:
 *
 *     1 - "mu" = global regularization weight.
 *     2 - "threshold" = value of threshold EPS.
 *
 * SEE ALSO: rgl_new, smooth3.
 */

func _rgl_simple_setup(this)
{
  h_set, this, state = 0, eps = 0.0, mu = 1.0;
  _rgl_simple_setup_l2, this;
  return this;
}

func _rgl_simple_update(this, x)
{
  return this.update(this, x);
}

func _rgl_simple_get_penalty(this, x)
{
  return this.get_penalty(this, x);
}

func _rgl_simple_get_gradient(this, x)
{
  return this.get_gradient(this, x);
}

func _rgl_simple_get_attr(this, index)
{
  if (index == 1) return this.mu;
  if (index == 2) return this.eps;
  if (index == 3) return this.costname;
}

func _rgl_simple_set_attr(this, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight must be a non-negative real";
    }
    if (value != this.mu) {
      h_set, this, mu = double(value), state = 0;
    }
  }
  if (index == 2) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "threshold must be a non-negative real";
    }
    if (value != this.eps) {
      h_set, this, eps = double(value), state = 0;
    }
  }
  if (index == 3) {
    if (rgl_string_scalar(value)) {
      error, "cost must be a scalar string";
    }
    value = strcase(1, value);
    if (value == "L2-L1" || value == "L2_L1" || value == "L2L1" ||
        value == "L1-L2" || value == "L1_L2" || value == "L1L2") {
      value = "L2-L1";
      setup = _rgl_simple_setup_l2l1;
    } else if (value == "L2-L0" || value == "L2_L0" || value == "L2L0" ||
        value == "L0-L2" || value == "L0_L2" || value == "L0L2") {
      value = "L2-L0";
      setup = _rgl_simple_setup_l2l0;
    } else if (value == "L2" || value == "QUADRATIC" || value == "QUAD") {
      value = "L2";
      setup = _rgl_simple_setup_l2;
    } else {
      error, "unknown cost function";
    }
    if (value != this.costname) {
      setup, this;
    }
  }
}

/* Implementation of separable L2 regularization. */

func _rgl_simple_setup_l2(this)
{
  if (! is_void((cleanup = this.cleanup))) cleanup, this;
  h_set, this, costname = "L2", state = 0,
    cleanup = symlink_to_variable(_rgl_simple_cleanup_l2),
    update = symlink_to_variable(_rgl_simple_update_l2),
    get_penalty = symlink_to_variable(_rgl_simple_get_penalty_l2),
    get_gradient = symlink_to_variable(_rgl_simple_get_gradient_l2);
}

func _rgl_simple_cleanup_l2(this)
{
  h_set, this, state = 0;
}

func _rgl_simple_update_l2(this, x)
{
}

func _rgl_simple_get_penalty_l2(this, x)
{
  return this.mu*sum(x*x);
}

func _rgl_simple_get_gradient_l2(this, x)
{
  return (2.0*this.mu)*x;
}


/* Implementation of separable L2-L0 regularization. */

func _rgl_simple_setup_l2l0(this)
{
  if (! is_void((cleanup = this.cleanup))) cleanup, this;
  h_set, this, costname = "L2-L0", state = 0,
    cleanup = symlink_to_variable(_rgl_simple_cleanup_l2l0),
    update = symlink_to_variable(_rgl_simple_update_l2l0),
    get_penalty = symlink_to_variable(_rgl_simple_get_penalty_l2l0),
    get_gradient = symlink_to_variable(_rgl_simple_get_gradient_l2l0);
}

func _rgl_simple_cleanup_l2l0(this)
{
  h_pop, this, "temp1";
  h_pop, this, "temp2";
  h_set, this, state = 0;
}

func _rgl_simple_update_l2l0(this, x)
{
  temp1 = (1.0/this.eps)*x;
  temp2 = this.eps*atan(temp1);
  h_set, this, state = 1, temp1 = temp1, temp2 = temp2;
}

func _rgl_simple_get_penalty_l2l0(this, x)
{
  if (this.state < 1) {
    _rgl_simple_update_l2l0, this, x;
  }
  local temp2;
  eq_nocopy, temp2, this.temp2;
  return this.mu*sum(temp2*temp2);
}

func _rgl_simple_get_gradient_l2l0(this, x)
{
  if (this.state < 1) {
    _rgl_simple_update_l2l0, this, x;
  }
  local temp1;
  eq_nocopy, temp1, this.temp1;

  return (2.0*this.mu)*this.temp2/(1.0 + temp1*temp1);
}


/* Implementation of separable L2-L1 regularization. */

func _rgl_simple_setup_l2l1(this)
{
  if (! is_void((cleanup = this.cleanup))) cleanup, this;
  h_set, this, costname = "L2-L1", state = 0,
    cleanup = symlink_to_variable(_rgl_simple_cleanup_l2l1),
    update = symlink_to_variable(_rgl_simple_update_l2l1),
    get_penalty = symlink_to_variable(_rgl_simple_get_penalty_l2l1),
    get_gradient = symlink_to_variable(_rgl_simple_get_gradient_l2l1);
}

func _rgl_simple_cleanup_l2l1(this)
{
  h_pop, this, "temp";
  h_pop, this, "fx";
  h_pop, this, "gx";
  h_set, this, state = 0;
}

func _rgl_simple_update_l2l1(this, x)
{
  h_set, this, state = 0;
}

func _rgl_simple_get_penalty_l2l1(this, x)
{
  if (this.state < 1) {
    eps = this.eps;
    mu = this.mu;
    if (eps > 0.0) {
      temp = sqrt(x*x + eps*eps);
      h_set, this, state = 1, temp = temp,
        fx = mu*(sum(temp) - numberof(x)*eps);
    } else  {
      h_set, this, state = 1, fx = mu*abs(x);
    }
  }
  return this.fx;
}

func _rgl_simple_get_gradient_l2l1(this, x)
{
  if (this.state < 2) {
    mu = this.mu;
    eps = this.eps;
    if (this.state < 1) {
      if (eps > 0.0) {
        temp = sqrt(x*x + eps*eps);
        h_set, this, state = 2, temp = temp,
          fx = mu*(sum(temp) - numberof(x)*eps),
          gx = mu*x/temp;
      } else  {
        h_set, this, state = 2,
          fx = mu*abs(x),
          gx = mu*(double(x > 0.0) - double(x < 0.0));
      }
    } else {
      if (eps > 0.0) {
        h_set, this, state = 2, gx = mu*x/this.temp;
      } else  {
        h_set, this, state = 2, gx = mu*(double(x > 0.0) - double(x < 0.0));
      }
    }
  }
  return this.gx;
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "simple", ["mu", "threshold", "cost"];

/*---------------------------------------------------------------------------*/
/* L-p NORM */

local rgl_lpnorm;
/* DOCUMENT this = rgl_new("lpnorm");
 *
 *   fprior(x) = sum(sqrt(x*x + eps*eps)^p)^(1/p) - n^(1/p)*eps
 *             = (sum r^(p/2))^(1/p) - n^(1/p)*eps
 *             = (sum q)^(1/p) - n^(1/p)*eps
 *
 *   gprior(p) = (sum r^(p/2))^(1/p - 1)*x*r^(p/2 - 1)
 *             = (sum q)^(1/p - 1)*x*q/r
 *
 *   with: N = numberof(X), R = X*X + EPS*EPS  and Q = R^(P/2)
 *
 *   This regularization has only one hyper-parameter:
 *
 *     1 - "mu" = global regularization weight.
 *     2 - "power" = power of the norm.
 *     3 - "epsilon" = small value.
 *
 * SEE ALSO: rgl_new, smooth3.
 */

func _rgl_lpnorm_setup(self)
{
  return h_set(self, state=0, epsilon=1e-8, power=2.0);
}

func _rgl_lpnorm_update(self, x)
{
  h_set, self, state=0;
}

func _rgl_lpnorm_state1(self, x)
{
  w = self.mu/sqrt(2.0);
  eps = self.epsilon;
  p = self.power;
  r = x*x + eps*eps;
  q = r^(0.5*p);
  sq = sum(q);
  f = self.mu*(sq^(1.0/p) - numberof(x)^(1.0/p)*eps);
  h_set, self, state=1, f=f, r=r, q=q, sq=sq;
}

func _rgl_lpnorm_state2(self, x)
{
  if (self.state < 1) _rgl_lpnorm_state1, self, x;
  p = self.power;
  g = self.mu*(self.sq^(1.0/p - 1.0)*x*self.q/self.r);
  h_set, self, state=2, g=g;
}

func _rgl_lpnorm_get_penalty(self, x)
{
  if (self.state < 1) _rgl_lpnorm_state1, self, x;
  return self.f;
}

func _rgl_lpnorm_get_gradient(self, x)
{
  if (self.state < 2) _rgl_lpnorm_state2, self, x;
  return self.g;
}

func _rgl_lpnorm_get_attr(self, index)
{
  if (index == 1) return self.mu;
  if (index == 2) return self.power;
  if (index == 3) return self.epsilon;
}

func _rgl_lpnorm_set_attr(self, index, value)
{
  if (index == 1) {
    if (rgl_real_scalar(value) || value < 0.0) {
      error, "global regularization weight (MU) must be a non-negative real";
    }
    if (value != self.mu) {
      h_set, self, mu = double(value), w = [], state = 0;
    }
  }
  if (index == 2) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "POWER must be a strictly non-negative real";
    }
    if (value != self.power) {
      h_set, self, power = double(value), state = 0;
    }
  }
  if (index == 3) {
    if (rgl_real_scalar(value) || value <= 0.0) {
      error, "EPSILON must be a strictly non-negative real";
    }
    if (value != self.epsilon) {
      h_set, self, epsilon = double(value), state = 0;
    }
  }
}

/* After having defined all class methods, define the class itself: */
_rgl_class, "lpnorm", ["mu", "power", "epsilon"];

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func rgl_squared_distance(.., scale=)
/* DOCUMENT rgl_squared_distance(dimlist, ...);
 *   Return squared distance with respect to the geometrical center of an
 *   array of dimension list DIMLIST.  Keyword SCALE can be used to specify
 *   a scale along every dimensions (SCALE can be a scalar or a vector).
 *
 * SEE ALSO: make_dimlist.
 */
{
  dimlist = [0L]; // initial dimension list
  while (more_args()) {
    make_dimlist, dimlist, next_arg();
  }
  ndims = dimlist(1);
  if (! ndims) {
    return 0.0;
  }
  if (is_void(scale)) {
    scale = array(1.0, ndims);
  } else {
    scale *= array(1.0, ndims);
    if (numberof(scale) != ndims || structof(scale) != double) {
      error, "bad dimensions or data type for SCALE keyword";
    }
  }
  flag = 1n;
  for (k = numberof(dimlist); k >= 2; --k) {
    len = dimlist(k);
    c = scale(k - 1)*(indgen(len) - 0.5*(len + 1));
    if (flag) {
      flag = 0n;
      r2 = c*c;
    } else {
      r2 = r2(-,..) + c*c;
    }
  }
  return r2;
}

local rgl_integer_scalar, rgl_real_scalar, rgl_string_scalar;
/* DOCUMENT rgl_integer_scalar(var);
 *     -or- rgl_real_scalar(var);
 *     -or- rgl_string_scalar(var);
 *     -or- rgl_boolean(var);
 *
 *   Returns non-zero (-1) if VAR is not a scalar of a given type; otherwise,
 *   fix VAR to be a double or long for the caller and returns zero.
 *
 * SEE ALSO: is_real, is_integer, is_string, is_scalar.
 */

func rgl_integer_scalar(&x) /* DOCUMENTATION IS ABOVE */
{
  if (is_scalar(x) && is_integer(x)) {
    x = long(x);
    return 0n;
  } else {
    return -1n;
  }
}

func rgl_real_scalar(&x) /* DOCUMENTATION IS ABOVE */
{
  if (is_scalar(x) && (is_real(x) || is_integer(x))) {
    x = double(x);
    return 0n;
  } else {
    return -1n;
  }
}

func rgl_boolean(&x) /* DOCUMENTATION IS ABOVE */
{
  if (! is_array(x) || is_scalar(x)) {
    x = (x ? 1n : 0n);
    return 0n;
  } else {
    return -1n;
  }
}

func rgl_string_scalar(x) /* DOCUMENTATION IS ABOVE */
{
  return (is_scalar(x) && is_string(x) ? 0n : -1n);
}

func rgl_check_dimlist(&dimlist)
/* DOCUMENT rgl_check_dimlist(dimlist)
 *   Check dimension list DIMLIST and return number of elements of array with
 *   that dimension list.  Possibly fix DIMLIST in-place so that is it always
 *   like the result of dimsof (which see).  Zero is returned in case of
 *   error.
 *
 * SEE ALSO: dimsof.
 */
{
  if (is_void(dimlist)) {
    /* scalar */
    dimlist = [0];
    return 1;
  }
  if (! is_integer(dimlist)) {
    /* bad data type */
    return 0;
  }
  if ((rank = dimsof(dimlist)(1)) == 0) {
    if ((number = long(dimlist)) <= 0) {
      return 0;
    }
    dimlist = [1, number];
    return number;
  }
  if (rank != 1) {
    return 0; /* not a vector */
  }
  if ((length = numberof(dimlist)) == 1) {
    if (dimlist(1)) {
      return 0;
    }
    if (s != long) {
      dimlist = long(dimlist);
    }
    return 1;
  }
  if (dimlist(1) != length - 1 || min(dimlist) < 1) {
    return 0;
  }
  number = 1;
  for (k = length; k >= 2; --k) {
    number *= dimlist(k);
  }
  if (s != long) {
    dimlist = long(dimlist);
  }
  return number;
}

/*---------------------------------------------------------------------------*/
/* IDENTITY OPERATOR */

local rgl_identity;
/* DOCUMENT rgl_identity(x);
 *       or rgl_identity(x, transp);
 *   Linear operator implementing identity, simply returns X.
 *
 * SEE ALSO: linop_new.
 */
func rgl_identity(x, transp) { return x; }
#if 0
func _rgl_identity_evaluator(this, arg, transp) { return arg; }
rgl_identity = h_new();
h_evaluator, rgl_identity, "_rgl_identity_evaluator";
#endif

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
