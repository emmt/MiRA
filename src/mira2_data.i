/*
 * mira2_data.i -
 *
 * Management of data in MiRA.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of MiRA, a "Multi-aperture Image Reconstruction
 * Algorithm", <https://github.com/emmt/MiRA>.
 *
 * Copyright (C) 2001-2018, Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

func mira_vis2_data(  master) { return mira_get_data(master, "vis2"  ); }
func mira_vis_data(   master) { return mira_get_data(master, "vis"   ); }
func mira_visamp_data(master) { return mira_get_data(master, "visamp"); }
func mira_visphi_data(master) { return mira_get_data(master, "visphi"); }
func mira_t3_data(    master) { return mira_get_data(master, "t3"    ); }
func mira_t3amp_data( master) { return mira_get_data(master, "t3amp" ); }
func mira_t3phi_data( master) { return mira_get_data(master, "t3phi" ); }
func mira_get_data(master, class)
/* DOCUMENT db = mira_get_data(master, class);
         or db = mira_vis2_data(master);
         or db = mira_vis_data(master);
         or db = mira_visamp_data(master);
         or db = mira_visphi_data(master);
         or db = mira_t3_data(master);
         or db = mira_t3amp_data(master);
         or db = mira_t3phi_data(master);

      The first function yields the data-block collecting all data of a given
      kind in MASTER.  CLASS is one of "vis2", "vis", "visamp", "visphi", "t3",
      "t3amp" or "t3phi".  An empty result is returned if no such data exists
      in MASTER.

      The other functions are equivalent to calling the first one for a
      specific class.

      As a side effect of calling any of these functions, `mira_update(master)`
      is called to update internal information.

   SEE ALSO: mira_config, mira_update, mira_add_oidata.
 */
{
  /* Make sure everything up to date. */
  mira_update, master;

  /* Search the specific class of data. */
  for (db = _mira_first(master); db != nil; db = _mira_next(master, db)) {
    if (db.ops.class == class) {
      return db;
    }
  }
}

/*---------------------------------------------------------------------------*/
/* ADDING OI-FITS DATA TO AN INSTANCE */

local _mira_add_oidata;
func mira_add_oidata(master, .., quiet=, noise_method=, noise_level=,
                     atol=, rtol=)
/* DOCUMENT mira_add_oidata, master, data, ...;

     Append interferometric OI-FITS data to MiRA handle `master` (as created by
     `mira_new`, which see).  `data` and subsequent arguments are either OI-FITS
     file names or opaque OI-FITS handles as returned by `oifits_load` (which
     see).

     If keyword `quiet` is true, the operation is performed silently.

     Keywords `noise_method` and `noise_level` can be used to add some noise to
     the data.

     Keywords ATOL and RTOL can be used to specify the absolute and relative
     tolerances when comparing numerical values which should be identical
     (i.e., when merging different OI-FITS files).


   SEE ALSO: mira_new, oifits_load.
 */
{
  if (is_void(quiet)) quiet = 0;
  if (is_void(atol)) atol = 0.0;
  if (is_void(rtol)) rtol = 1e-8;
  local arg;
  while (more_args()) {
    eq_nocopy, arg, next_arg();
    if (is_string(arg)) {
      n = numberof(arg);
      for (i = 1; i <= n; ++i) {
        filename = arg(i);
        if (! (quiet & 2)) write, format="Loading file \"%s\"...\n", filename;
        data = oifits_load(filename, quiet=quiet, errmode=1n);
        _mira_add_oidata, master, data, atol, rtol;
      }
    } else if (! is_void(arg)) {
      _mira_add_oidata, master, arg, atol, rtol;
    }
  }
  return master;
}

local _mira_first, _mira_next;
/* DOCUMENT db = _mira_first(master);
         or db = _mira_next(master, db);

      Private routines to iterate over all data-blocks of MiRA instance
      `master`.   Typical usage is:

          for (db = _mira_first(master); db; db = _mira_next(master, db)) {
              ...
          }

*/
func _mira_first(master) { return master.first; }
func _mira_next(master, db) { return db.next; }

func _mira_add_oidata(master, oidata, atol, rtol)
{
  /* Arguments shared with caller. */
  extern quiet, noise_method, noise_level;

  /* Restrict the OI-FITS data to a common target. */
  local oifits_selected_target;
  oidata = oifits_select_target(oidata, master.target);
  h_set, master, target = oifits_selected_target;

  /* Flag the data where values are not finite or where the errors are not
     strictly greater than zero.  All the comparisons must be carried on in a
     specific order to make sure to avoid floating-point errors (i.e., we
     discard non-finite values first).  Since testing IEEE floating-point type
     is rather time consuming, we try to test the most restricted set of values
     as possible.  Performing these tests now also avoid retesting when
     updating the subset of selected data.  This is why the flags are saved
     into each oi-datablock in field "miraflags". */
  for (oidb = oifits_first(oidata); oidb; oidb = oifits_next(oidata, oidb)) {
    flags = [];
    type = oifits_get_type(oidb);
    if (type == OIFITS_TYPE_VIS) {
      flags = _mira_build_flags(oifits_get_flag(oidata, oidb),
                                MIRA_FIT_VISAMP,
                                oifits_get_visamp(oidata, oidb),
                                oifits_get_visamperr(oidata, oidb),
                                MIRA_FIT_VISPHI,
                                oifits_get_visphi(oidata, oidb),
                                oifits_get_visphierr(oidata, oidb));
    } else if (type == OIFITS_TYPE_VIS2) {
      flags = _mira_build_flags(oifits_get_flag(oidata, oidb),
                                MIRA_FIT_VIS2,
                                oifits_get_vis2data(oidata, oidb),
                                oifits_get_vis2err(oidata, oidb));
    } else if (type == OIFITS_TYPE_T3) {
      flags = _mira_build_flags(oifits_get_flag(oidata, oidb),
                                MIRA_FIT_T3AMP,
                                oifits_get_t3amp(oidata, oidb),
                                oifits_get_t3amperr(oidata, oidb),
                                MIRA_FIT_T3PHI,
                                oifits_get_t3phi(oidata, oidb),
                                oifits_get_t3phierr(oidata, oidb));
    }
    if (is_array(flags)) {
      h_set, oidb, miraflags = flags;
    }
  }

  /* Perhaps add noise. */
  if (! is_void(noise_method)) {
    oifits_add_noise, oidata, noise_method, noise_level;
  }

  /* Merge OI-FITS data. */
  oifits_merge, dest=master.oidata, oidata, quiet=quiet, atol=atol, rtol=rtol;

  return h_set(master, stage = 0);
}

func _mira_build_flags(oiflag, ..)
{
  local miraflags, flags, index, bits, dat, err;
  while (more_args()) {
    eq_nocopy, bits, next_arg();
    eq_nocopy, dat,  next_arg();
    eq_nocopy, err,  next_arg();
    flags = (is_void(oiflag) ? array(1n, dimsof(dat)) : !oiflag);
    index = where(flags); // index of pre-selected data
    if (is_array(index)) {
      _mira_build_flags_worker, (ieee_test(dat(index)) != 0n);
    }
    if (is_array(index)) {
      _mira_build_flags_worker, (ieee_test(err(index)) != 0n);
    }
    if (is_array(index)) {
      _mira_build_flags_worker, (err(index) <= 0.0);
    }
    if (is_void(miraflags)) {
      miraflags = array(0n, dimsof(flags));
    }
    if (is_array(index)) {
      miraflags(index) |= bits;
    }
  }
  return miraflags;
}

func _mira_build_flags_worker(test)
{
  extern flags, index; // variables shared with caller
  j = where(test);
  if (is_array(j)) {
    flags(index(j)) = 0n;
    index = where(flags);
  }
}

func _mira_select_data(master)
/* DOCUMENT _mira_select_data, master;

     Private routine to carry on OI-FITS data selection according to the
     settings in MiRA instance `master`.  This routine must be called after
     changing any data selection option.  Attribute `master.stage` is 0 on
     entry and is set to 1 on return.  When called as a function, `master` is
     returned.

   SEE ALSO: mira_update, _mira_define_xform.
 */
{
  /* Check assumptions and reset selected data blocks and coordinates. */
  if (master.stage != 0) {
    throw, "private routine called at wrong stage";
  }
  oidata = master.oidata;
  wavemin = master.wavemin;
  wavemax = master.wavemax;
  h_set, master, first = [], stage = 0,
    img_wavemin = 0.0, img_wavemax = 0.0, img_wave = 0.0,
    coords = h_new(mode = 0, u = [], v = [], wave = [], band = []);

  /* Iterate over OI-FITS data blocks for selection. */
  local wave, band;
  for (oidb = oifits_first(oidata); oidb; oidb = oifits_next(oidata, oidb)) {
    type = oifits_get_type(oidb);
    if (type == OIFITS_TYPE_VIS) {
      /* Differential visibilities or correlated fluxes are not yet supported
         in MiRA. */
      bits = (master.flags & MIRA_FIT_VIS);
      if ((bits & MIRA_FIT_VISAMP) != 0n) {
        amptyp = oifits_get_amptyp(oidata, oidb);
        if (! is_void(amptyp) && ! mira_same_strings(amptyp, "absolute")) {
          bits &= ~MIRA_FIT_VISAMP;
        }
      }
      if ((bits & MIRA_FIT_VISPHI) != 0n) {
        phityp = oifits_get_phityp(oidata, oidb);
        if (! is_void(phityp) && ! mira_same_strings(phityp, "absolute")) {
          bits &= ~MIRA_FIT_VISPHI;
        }
      }
    } else if (type == OIFITS_TYPE_VIS2) {
      bits = (master.flags & MIRA_FIT_VIS2);
    } else if (type == OIFITS_TYPE_T3) {
      bits = (master.flags & MIRA_FIT_T3);
    } else {
      bits = 0n;
    }
    if (! bits) {
      /* This data block is completely ignored. */
      continue;
    }

    /* Collect common parameters before calling the specific constructors. */
    eq_nocopy, wave, oifits_get_eff_wave(oidata, oidb);
    eq_nocopy, band, oifits_get_eff_band(oidata, oidb);
    select = ((oidb.miraflags & bits) != 0n);
    dims = dimsof(select);
    if (numberof(dims) != 3) {
      error, "bug: SELECT must be a 2D array";
    }
    nrows = dims(2);
    ncols = dims(3);
    if (! vector_double(wave) || numberof(wave) != ncols) {
      error, "bug: WAVE and SELECT have incompatible dimensions";
    }
    if (! vector_double(band) || numberof(band) != ncols) {
    error, "bug: BAND and SELECT have incompatible dimensions";
    }
    if (wavemin > 0.0) {
      select &= (wave >= wavemin)(-,);
    }
    if (wavemax < MIRA_HUGE) {
      select &= (wave <= wavemax)(-,);
    }
    select = where(select);
    if (is_array(select)) {
      if (type == OIFITS_TYPE_VIS) {
        _mira_append_vis_data, master, oidata, oidb, wave, band,
          nrows, ncols, select;
      } else if (type == OIFITS_TYPE_VIS2) {
        _mira_append_vis2_data, master, oidata, oidb, wave, band,
          nrows, ncols, select;
      } else if (type == OIFITS_TYPE_T3) {
        _mira_append_t3_data, master, oidata, oidb, wave, band,
          nrows, ncols, select;
      }
    }
  }

  /* Estimate the wavelength, set stage and return. */
  if (numberof(master.coords.wave) > 0) {
    img_wavemin = min(master.coords.wave);
    img_wavemax = max(master.coords.wave);
    h_set, master,
      img_wavemin = img_wavemin,
      img_wavemax = img_wavemax,
      img_wave = (img_wavemin + img_wavemax)/2.0;
  }
  return h_set(master, stage = 1);
}

func _mira_grow_coordinates(master, u, v, wave, band)
/* DOCUMENT idx = _mira_grow_coordinates(master, u, v, wave, band);

      appends coordinates to MiRA instance `master` and return their signed
      indices in the memorized list.  The sign of the result indicates whether
      the baseline coordinates have been reversed while the absolute of the
      result yields the indices of the new coordinates in the list of
      coordinates memorized by `master`.  When all coordinates have been added,
      `_mira_grow_coordinates` should be called to reduce the list of
      coordinates.

   SEE ALSO:
 */
{
  /* Check existing coordinates. */
  if (master.stage != 0) {
    throw, "private routine called at wrong stage";
  }
  coords = master.coords;
  n = numberof(coords.u);
  if (numberof(coords.v) != n || numberof(coords.wave) != n ||
      numberof(coords.band) != n ||
      (n > 0 && !(is_vector(coords.u) && structof(coords.u) == double &&
                  is_vector(coords.v) && structof(coords.v) == double &&
                  is_vector(coords.wave) && structof(coords.wave) == double &&
                  is_vector(coords.band) && structof(coords.band) == double))) {
    error, "corrupted list of coordinates";
  }

  /* Check types and values and broadcast to convert to same type and
     dimensions. */
  if (identof(u) > Y_DOUBLE) {
    error, "`u` must be real";
  }
  if (identof(v) > Y_DOUBLE) {
    error, "`v` must be real";
  }
  if (identof(wave) > Y_DOUBLE || min(wave) <= 0) {
    error, "`wave` must be strictly positive real(s)";
  }
  if (identof(band) > Y_DOUBLE || min(band) <= 0) {
    error, "`band` must be strictly positive real(s)";
  }
  dims = dimsof(u, v, wave, band);
  if (is_void(dims)) {
    error, "`u`, `v`, `wave` and `band` must be conformable";
  }
  zero = array(double, dims);
  u = (u + zero)(*);
  v = (v + zero)(*);
  wave = (wave + zero)(*);
  band = (band + zero)(*);
  zero = [];

  /* Perhaps round `u` and `v`. */
  if ((p = master.baseline_precision) > 0) {
    u = p*round(u/p);
    v = p*round(v/p);
  }

  /* Re-orient the coordinates (and change the sign of the index to indicate
     this). */
  idx = indgen(n+1 : n+numberof(u));
  i = where((u < 0)|(!u & (v < 0)));
  if (is_array(i)) {
    idx(i) = -idx(i);
    u(i) = -u(i);
    v(i) = -v(i);
  }

  /* Append the new coordinates to the existing ones and return the (signed)
     indices. */
  h_set, coords,
    u = grow(coords.u, u),
    v = grow(coords.v, v),
    wave = grow(coords.wave, wave),
    band = grow(coords.band, band);
  return idx;
}

func _mira_check_baselines(&u, &v, nrows)
{
  if (! vector_double(u) || numberof(u) != nrows) {
    error, swrite(format="U should be a vector of %d real value(s)", nrows);
  }
  if (! vector_double(v) || numberof(v) != nrows) {
    error, swrite(format="V should be a vector of %d real value(s)", nrows);
  }
}
errs2caller, _mira_check_baselines;

func _mira_append_vis2_data(master, oidata, oidb, wave, band,
                            nrows, ncols, select)
/* DOCUMENT _mira_append_vis2_data, master, oidata, oidb, wave, band,
                                    nrows, ncols, select;

     Append OI-VIS2 data in datablock OIDB of OI-FITS instance OIDATA to MiRA
     instance MASTER.  WAVE and BAND are the wavelength and spectral bandwidth
     of the data, NROWS and NCOLS are the dimensions of the data, SELECT is the
     list of indices where to keep data.  This private subroutine is only
     called if there are any data to consider.

   SEE ALSO: _mira_select_data, _mira_append_vis_data _mira_append_t3_data.
 */
{
  local u, v;
  eq_nocopy, u, oifits_get_ucoord(oidata, oidb);
  eq_nocopy, v, oifits_get_vcoord(oidata, oidb);
  _mira_check_baselines, u, v, nrows;
  idx = _mira_grow_coordinates(master,
                               u(,-:1:ncols)(select),
                               v(,-:1:ncols)(select),
                               wave(-:1:nrows,)(select),
                               band(-:1:nrows,)(select));
  dat = oifits_get_vis2data(oidata, oidb)(select);
  err = oifits_get_vis2err(oidata, oidb)(select);
  db = _mira_find_datablock(master, MIRA_FIT_VIS2, _MIRA_VIS2_OPS);
  mira_grow_fields, db,
    "idx", abs(idx),
    "dat", dat,
    "wgt", 1.0/(err*err);
  if (MIRA_DEBUG) {
    inform, "will fit %d OI_VIS2 data", numberof(idx);
  }
}

func _mira_append_vis_data(master, oidata, oidb, wave, band,
                           nrows, ncols, select)
/* DOCUMENT _mira_append_vis_data, master, oidata, oidb, wave, band,
                                   nrows, ncols, select;

     Append OI-VIS data in datablock OIDB of OI-FITS instance OIDATA to MiRA
     instance MASTER.  WAVE and BAND are the wavelength and spectral bandwidth
     of the data, NROWS and NCOLS are the dimensions of the data, SELECT is the
     list of indices where to keep data.  This private subroutine is only
     called if there are any data to consider.

   SEE ALSO: _mira_select_data, _mira_append_vis2_data _mira_append_t3_data.
 */
{
  /* Get coordinates of the data and broadcast their dimensions. */
  local u, v;
  eq_nocopy, u, oifits_get_ucoord(oidata, oidb);
  eq_nocopy, v, oifits_get_vcoord(oidata, oidb);
  _mira_check_baselines, u, v, nrows;
  u = u(,-:1:ncols);
  v = v(,-:1:ncols);
  wave = wave(-:1:nrows,);
  band = band(-:1:nrows,);

  /* Sub-select data according to chosen options and to data flags.  If the
     convex approximation is chosen, we first select complex visibility data
     because this requires that the amplitude and the phase be both available
     and imposes to process them together. */
  DEG = MIRA_DEGREE; // to convert degrees into radians
  flags = oidb.miraflags(select);
  if ((master.flags & MIRA_CONVEX_APPROX) == MIRA_CONVEX_APPROX &&
      (master.flags & MIRA_FIT_VIS) == MIRA_FIT_VIS) {
    j = where((flags & MIRA_FIT_VIS) == MIRA_FIT_VIS);
    if (is_array(j)) {
      flags(j) &= ~MIRA_FIT_VIS; // avoid accounting of these data twice
      j = select(j);
      idx = _mira_grow_coordinates(master, u(j), v(j), wave(j), band(j));
      ws = mira_polar_to_cartesian(oifits_get_visamp(oidata, oidb)(j),
                                   oifits_get_visamperr(oidata, oidb)(j),
                                   oifits_get_visphi(oidata, oidb)(j)*DEG,
                                   oifits_get_visphierr(oidata, oidb)(j)*DEG,
                                   "visibility", goodman = 0n);
      db = _mira_find_datablock(master, MIRA_FIT_VIS|MIRA_CONVEX_APPROX,
                                _MIRA_VIS_OPS);
      sgn = double(sign(idx));
      mira_grow_fields, db,
        "idx", abs(idx),
        "re", ws.re, "im", sgn*ws.im,
        "wrr", ws.wrr, "wii", ws.wii, "wri", sgn*ws.wri;
      if (MIRA_DEBUG) {
        inform, "will fit %d OI_VIS data", numberof(idx);
      }
    }
  }

  /* Select amplitude-only data if this option is chosen. */
  if ((master.flags & MIRA_FIT_VISAMP) == MIRA_FIT_VISAMP) {
    j = where((flags & MIRA_FIT_VISAMP) == MIRA_FIT_VISAMP);
    if (is_array(j)) {
      flags(j) &= ~MIRA_FIT_VISAMP; // avoid accounting of these data twice
      j = select(j);
      idx = _mira_grow_coordinates(master, u(j), v(j), wave(j), band(j));
      dat = oifits_get_visamp(oidata, oidb)(j);
      err = oifits_get_visamperr(oidata, oidb)(j);
      db = _mira_find_datablock(master, MIRA_FIT_VISAMP, _MIRA_VISAMP_OPS);
      mira_grow_fields, db,
        "idx", abs(idx),
        "dat", dat,
        "wgt", 1.0/(err*err);
      if (MIRA_DEBUG) {
        inform, "will fit %d OI_VISAMP data", numberof(idx);
      }
    }
  }

  /* Select phase-only data if this option is chosen. */
  if ((master.flags & MIRA_FIT_VISPHI) == MIRA_FIT_VISPHI) {
    j = where((flags & MIRA_FIT_VISPHI) == MIRA_FIT_VISPHI);
    if (is_array(j)) {
      flags(j) &= ~MIRA_FIT_VISPHI; // avoid accounting of these data twice
      j = select(j);
      idx = _mira_grow_coordinates(master, u(j), v(j), wave(j), band(j));
      dat = oifits_get_visphi(oidata, oidb)(j)*DEG;
      err = oifits_get_visphierr(oidata, oidb)(j)*DEG;
      db = _mira_find_datablock(master, MIRA_FIT_VISPHI, _MIRA_VISPHI_OPS);
      mira_grow_fields, db,
        "idx", abs(idx),
        "dat", double(sign(idx))*dat,
        "wgt", 1.0/(err*err);
      if (MIRA_DEBUG) {
        inform, "will fit %d OI_VISPHI data", numberof(idx);
      }
    }
  }
}

func _mira_append_t3_data(master, oidata, oidb, wave, band,
                          nrows, ncols, select)
/* DOCUMENT _mira_append_t3_data, master, oidata, oidb, wave, band,
                                  nrows, ncols, select;

     Append OI-T3 data in datablock OIDB of OI-FITS instance OIDATA to MiRA
     instance MASTER.  WAVE and BAND are the wavelength and spectral bandwidth
     of the data, NROWS and NCOLS are the dimensions of the data, SELECT is the
     list of indices where to keep data.  This private subroutine is only
     called if there are any data to consider.

   SEE ALSO: _mira_select_data, _mira_append_vis2_data _mira_append_vis_data.
 */
{
  /* Get coordinates of the data. */
  local u1, v1, u2, v2;
  eq_nocopy, u1, oifits_get_u1coord(oidata, oidb);
  eq_nocopy, v1, oifits_get_v1coord(oidata, oidb);
  eq_nocopy, u2, oifits_get_u2coord(oidata, oidb);
  eq_nocopy, v2, oifits_get_v2coord(oidata, oidb);
  _mira_check_baselines, u1, v1, nrows;
  _mira_check_baselines, u2, v2, nrows;
  u3 = -(u1 + u2);
  v3 = -(v1 + v2);

  /* Broadcast coordinates dimensions. */
  u1 = u1(,-:1:ncols);
  v1 = v1(,-:1:ncols);
  u2 = u2(,-:1:ncols);
  v2 = v2(,-:1:ncols);
  u3 = u3(,-:1:ncols);
  v3 = v3(,-:1:ncols);
  wave = wave(-:1:nrows,);
  band = band(-:1:nrows,);

  /* Sub-select data according to chosen options and to data flags.  The
     rationale is the same as in _mira_append_vis_data. */
  DEG = MIRA_DEGREE; // to convert degrees into radians
  flags = oidb.miraflags(select);
  if ((master.flags & MIRA_CONVEX_APPROX) == MIRA_CONVEX_APPROX &&
      (master.flags & MIRA_FIT_T3) == MIRA_FIT_T3) {
    j = where((flags & MIRA_FIT_T3) == MIRA_FIT_T3);
    if (is_array(j)) {
      flags(j) &= ~MIRA_FIT_T3; // avoid accounting of these data twice
      j = select(j);
      wtmp = wave(j);
      btmp = band(j);
      idx1 = _mira_grow_coordinates(master, u1(j), v1(j), wtmp, btmp);
      idx2 = _mira_grow_coordinates(master, u2(j), v2(j), wtmp, btmp);
      idx3 = _mira_grow_coordinates(master, u3(j), v3(j), wtmp, btmp);
      ws = mira_polar_to_cartesian(oifits_get_t3amp(oidata, oidb)(j),
                                   oifits_get_t3amperr(oidata, oidb)(j),
                                   oifits_get_t3phi(oidata, oidb)(j)*DEG,
                                   oifits_get_t3phierr(oidata, oidb)(j)*DEG,
                                   "bispectrum", goodman = 0n);
      db = _mira_find_datablock(master, MIRA_FIT_T3|MIRA_CONVEX_APPROX,
                                _MIRA_T3_OPS);
      idx = mira_cat(idx1, idx2, idx3);
      mira_grow_fields, db,
        "idx", abs(idx),
        "sgn", double(sign(idx)),
        "re", ws.re, "im", ws.im,
        "wrr", ws.wrr, "wii", ws.wii, "wri", ws.wri;
      if (MIRA_DEBUG) {
        inform, "will fit %d OI_T3 data", numberof(idx);
      }
    }
  }

  /* Select amplitude-only data if this option is chosen. */
  if ((master.flags & MIRA_FIT_T3AMP) == MIRA_FIT_T3AMP) {
    j = where((flags & MIRA_FIT_T3AMP) == MIRA_FIT_T3AMP);
    if (is_array(j)) {
      flags(j) &= ~MIRA_FIT_T3AMP; // avoid accounting of these data twice
      j = select(j);
      wtmp = wave(j);
      btmp = band(j);
      idx1 = _mira_grow_coordinates(master, u1(j), v1(j), wtmp, btmp);
      idx2 = _mira_grow_coordinates(master, u2(j), v2(j), wtmp, btmp);
      idx3 = _mira_grow_coordinates(master, u3(j), v3(j), wtmp, btmp);
      dat = oifits_get_t3amp(oidata, oidb)(j);
      err = oifits_get_t3amperr(oidata, oidb)(j);
      db = _mira_find_datablock(master, MIRA_FIT_T3AMP, _MIRA_T3AMP_OPS);
      idx = mira_cat(idx1, idx2, idx3);
      mira_grow_fields, db,
        "idx", abs(idx),
        "sgn", double(sign(idx)),
        "dat", dat,
        "wgt", 1.0/(err*err);
      if (MIRA_DEBUG) {
        inform, "will fit %d OI_T3AMP data", numberof(idx);
      }
    }
  }

  /* Select phase-only data if this option is chosen. */
  if ((master.flags & MIRA_FIT_T3PHI) == MIRA_FIT_T3PHI) {
    j = where((flags & MIRA_FIT_T3PHI) == MIRA_FIT_T3PHI);
    if (is_array(j)) {
      flags(j) &= ~MIRA_FIT_T3PHI; // avoid accounting of these data twice
      j = select(j);
      wtmp = wave(j);
      btmp = band(j);
      idx1 = _mira_grow_coordinates(master, u1(j), v1(j), wtmp, btmp);
      idx2 = _mira_grow_coordinates(master, u2(j), v2(j), wtmp, btmp);
      idx3 = _mira_grow_coordinates(master, u3(j), v3(j), wtmp, btmp);
      dat = oifits_get_t3phi(oidata, oidb)(j)*DEG;
      err = oifits_get_t3phierr(oidata, oidb)(j)*DEG;
      db = _mira_find_datablock(master, MIRA_FIT_T3PHI, _MIRA_T3PHI_OPS);
      idx = mira_cat(idx1, idx2, idx3);
      mira_grow_fields, db,
        "idx", abs(idx),
        "sgn", double(sign(idx)),
        "dat", dat,
        "wgt", 1.0/(err*err);
      if (MIRA_DEBUG) {
        inform, "will fit %d OI_T3PHI data", numberof(idx);
      }
    }
  }
}

func _mira_find_datablock(master, bits, ops)
/* DOCUMENT db = _mira_find_datablock(master, bits);
         or db = _mira_find_datablock(master, bits, ops);

     Yields a MiRA datablock identified by BITS (a bitwise combination of
     MIRA_FIT_* constants) in MiRA instance MASTER.  If no such datablock is
     found, a void result is returned unless OPS is provided.  Optional
     argument OPS is a table of virtual operations to work with the contents of
     a new datablock of a given kind (i.e., as specified by BITS).

     ------------------------------------
     Bits             Table of Operations
     ------------------------------------
     MIRA_FIT_VIS2    _MIRA_VIS2_OPS
     MIRA_FIT_VIS     _MIRA_VIS_OPS
     MIRA_FIT_VISAMP  _MIRA_VISAMP_OPS
     MIRA_FIT_VISPHI  _MIRA_VISPHI_OPS
     MIRA_FIT_T3      _MIRA_T3_OPS
     MIRA_FIT_T3AMP   _MIRA_T3AMP_OPS
     MIRA_FIT_T3PHI   _MIRA_T3PHI_OPS
     ------------------------------------

   SEE ALSO: _mira_select_data.
 */
{
  for (db = _mira_first(master); db; db = _mira_next(master, db)) {
    if (db.bits == bits) {
      return db;
    }
  }
  if (! is_void(ops)) {
    db = h_new(bits=bits, next=master.first, ops=ops, idx=where(0));
    h_set, master, first=db;
  }
  return db;
}

/*---------------------------------------------------------------------------*/
