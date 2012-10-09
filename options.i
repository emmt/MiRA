/*
 * options.i --
 *
 * Parsing of command line options in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2009 Eric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/ or redistribute the software under the terms of the CeCILL-C license
 * as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 *-----------------------------------------------------------------------------
 */

/* FIXME: check consistency of option list
 * FIXME: add an OPT_USAGE type
 * FIXME: deal with short/long option flags, and use more standard syntax
 */

local opt_init;
local OPT_FLAG, OPT_INTEGER, OPT_REAL, OPT_STRING, OPT_HELP, OPT_VERSION; 
func opt_init(usage, brief, ops)
/* DOCUMENT tab = opt_init(usage, brief, ops)

   Build table for fast parsing of command line options.  USAGE is the syntax
   of the command, BRIEF is a short description of the command.  OPS is a
   ops of accepted options.  The items in OPS are lists of 5 items:
   
     ops = _lst(_lst(name1, defval1, units1, type1, descr1),
                _lst(name2, defval2, units2, type2, descr2), ...);

   with NAMEi the option name, DEFVALi the default value of the option, UNITSi
   the name of units or type of the option, TYPEi the identifier of the option
   type (one of: OPT_FLAG, OPT_INTEGER, OPT_REAL, OPT_STRING, OPT_HELP, or
   OPT_VERSION) and DESCRi the description of the option.

   The options will be used as (with one or two leading dashes):

     --NAMEi              for an OPT_FLAG, OPT_HELP, or OPT_VERSION option
     --NAMEi=VALUE        for an OPT_INTEGER, OPT_REAL, or OPT_STRING option

   There should be at most one OPT_HELP and one OPT_VERSION options.

   OPT_VERSION: If the option is set on the command line, the version will be
   printed.  The default value is the version string, the units are ignored.
   
   OPT_HELP: If the option is set on the command line, a help message will be
   printed with a list of all options.  The default value and the units are
   ignored.

   OPT_FLAG: If the option is set on the command line, the value of the option
   will be true (1n); otherwise, the value is false (0n).  Given default value
   and units are ignored.
   
   
   For instance:

     NULL = [];
     ops = _lst(_lst("input", NULL, "FILE", OPT_STRING, "name of input file"),
                _lst("safe", NULL, NULL, OPT_FLAG, "use safe mode"),
                _lst("level", 2, "INTEGER", OPT_INTEGER,
                     "level of operation"),
                _lst("scale", 0.3, "FACTOR", OPT_REAL,
                     "output scaling factor"),
                _lst("help", NULL, NULL, OPT_HELP, "print this help"),
                _lst("version", "2.1.4", NULL, OPT_VERSION,
                     "print option number"));

     
   SEE ALSO: opt_parse.
 */
{
  n = _len(ops);
  if (n < 1) opt_error, "empty option list";
  options = array(string, n);
  tab = h_new(":options", options, ":usage", usage, ":brief", brief);
  k = 0;
  while (ops) {
    item = _nxt(ops);
    if (_len(item) != 5) {
      opt_error, swrite(format="syntax error in option list (%d)", k+1);
    }
    name  = _car(item, 1);
    value = _car(item, 2);
    units = _car(item, 3);
    type  = _car(item, 4);
    descr = _car(item, 5);
    if (strglob("*[: ]*", name)) {
      opt_error, swrite(format="bad option name \"%s\" in option list (%d)",
                    name, k+1);
    }
    options(++k) = name;
    h_set, tab,
      name + ":defval", value,
      name + ":units", units,
      name + ":type", type,
      name + ":descr", descr;
  }
  return tab;
}

OPT_FLAG = 0;
OPT_INTEGER = 1;
OPT_REAL = 2;
OPT_STRING = 3;
OPT_HELP = 4;
OPT_VERSION = 5;

func opt_parse(tab, &argv)
/* DOCUMENT opt = opt_parse(tab, argv);

     Parse command line options from array of string ARGV and according to
     compiled rules in TAB.  If one of the special "help" or "version" options
     appears in ARGV, then an empty result is returned.  Otherwise, the result
     is a hash-table with members properly set accounting for default values
     in TAB and to given values in ARGV.  Before return, ARGV is set with the
     remaining unprocessed options.

     The returned result is a hash table which can be used as opt("NAME") or
     opt.NAME to get the value of option NAME.
   
     
   SEE ALSO: opt_init, get_argv.
 */
{
  local options;
  eq_nocopy, options, tab(":options");
  opt = h_new();
  for (k = numberof(options); k >= 1; --k) {
    name = options(k);
    h_set, opt, name, tab(name + ":defval");
  }
  nil = string();
  argc = numberof(argv);
  for (k = 1; k <= argc; ++k) {
    arg = argv(k);
    if (strpart(arg, 1:1) != "-") {
      argv = argv(k:);
      break;
    }
    if (strpart(arg, 2:2) == "-") {
      if (arg == "--") {
        argv = (++k < argc ? argv(k:) : []);
      }
      off = 2;
    } else {
      off = 1;
    }
    sel = strfind("=", arg, off);
    if (sel(2) >= 0) {
      name = strpart(arg, off+1:sel(1));
      value = strpart(arg, sel(1)+2:);
    } else {
      name = strpart(arg, off+1:);
      value = [];
    }
    type = tab(name + ":type");
    units = tab(name + ":units");
    if (is_void(type)) {
      opt_error, "unknown option: " + arg;
    }
    if (is_void(value)) {
      if (type == OPT_FLAG) {
        value = 1n;
      } else if (type == OPT_HELP) {
        write, format="%s\n", tab(":usage");
        write, format="%s\n", tab(":brief");
        n = numberof(options);
        for (k = 1; k <= n; ++k) {
          name = options(k);
          value = opt(name);
          type = tab(name + ":type");
          units = tab(name + ":units");
          str = "-" + name;
          if (! is_void(units)) {
            str += "=" + tab(name + ":units");
          }
          if (! is_void(value) && type != OPT_VERSION) {
            def = " (default " + print(value) + ")";
          } else {
            def = "";
          }
          write, format="  %-17s %s%s\n", str, tab(name + ":descr"), def;
        }
        return;
      } else if (type == OPT_VERSION) {
        write, format="version: %s\n", opt(name);
        return;
      } else {
        opt_error, "option \"" + name + "\" takes a value";
      }
    } else {
      if (type == OPT_INTEGER) {
        temp = 0L;
        dummy = nil;
        if (sread(value, temp, dummy) != 1) {
          opt_error, "expecting integer value for option \"" + name + "\"";
        }
        value = temp;
      } else if (type == OPT_REAL) {
        temp = 0.0;
        dummy = nil;
        if (sread(value, temp, dummy) != 1) {
          opt_error, ("expecting floating-point value for option \""
                      + name + "\"");
        }
        value = temp;
      } else if (type == OPT_STRING) {
        if (strlen(value) <= 0) {
          value = nil;
        }
      } else {
        opt_error, "option \"" + name + "\" takes no value";
      }
    }
    h_set, opt, name, value;
  }
  return opt;
}

func opt_usage(tab, msg)
{
  if (msg && strlen(msg)) {
    msg += "\n" + tab(":usage");
  } else {
    msg = tab(":usage");
  }
  if (batch()) {
    write, format="%s\n", msg;
    quit;
  }
}

local opt_error;
/* DOCUMENT opt_error, msg;

     If in batch mode, this subroutine print the error message and exit
     Yorick.  Otherwise, this subroutine is just an alias to the error
     subroutine of Yorick.  Note that batch mode is detected when the source
     is compiled par Yorick parser *not* at runtime when the subroutine is
     called.

   SEE ALSO: error, batch, quit.
 */
func opt_error(msg)
{
  write, format="%s\n", msg;
  quit;
}
if (! batch()) {
  opt_error = error;
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 78
 * coding: utf-8
 * End:
 */
