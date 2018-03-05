/*
 * mira2.i -
 *
 * Implement version 2 of MiRA (Multi-aperture Image Reconstruction Algorithm)
 * in Yeti/Yorick.
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

MIRA_VERSION = "2.0.0a";
local MIRA_SRCDIR;
func _mira_init(path)
{
  extern MIRA_SRCDIR;
  fullpath = filepath(path);
  i = strfind("/", path, back=1)(2);
  if (i <= 0) {
    error, "_mira_init must be called with the result of current_include()";
  }
  MIRA_SRCDIR = strpart(fullpath, 1:i);
}
_mira_init, current_include();
require, "TiPi.i";
//include, MIRA_SRCDIR+"TiPi.i";
include, MIRA_SRCDIR+"mira2_utils.i";
include, MIRA_SRCDIR+"mira2_config.i";
include, MIRA_SRCDIR+"mira2_data.i";
include, MIRA_SRCDIR+"mira2_xform.i";
include, MIRA_SRCDIR+"mira2_cost.i";

/*
  MiRA structure:
     master.oidata = OI-FITS data
     master.stage = configuration stage of MiRA (0 initially, ≥ 1 after data
                    selection, ≥ 2 means all settings have been updated)
     master.coords = table of coordinates
     master.coords.u = baseline U-coordinate [m]
     master.coords.v = baseline V-coordinate [m]
     master.coords.wave = effective wavelength [m]
     master.coords.band = effective bandwidth [m]
     master.model.img = current model image
     master.model.vis = model complex visibility
     master.model.phi = phase of model complex visibility
     master.model.amp = amplitude of model complex visibility
     master.model.re = real part of model complex visibility
     master.model.im = imaginary part of model complex visibility
     master.xform = linear transform
     master.smearingfactor = scaling factor for the bandwidth smearing
     master.smearingfunction = function to shape the bandwidth smearing
     master.wavemin = minimal wavelength (in meters)
     master.wavemax = maximal wavelength (in meters)
     master.pixelsize = pixel size (in radians)
     master.dims = image dimensions (in pixels) as returned by `dimsof`
     master.target = name of target (as found in OI-FITS data)
     master.baseline_precision = precision for rounding baselines (in meters)
     master.flags = processiong options

     FIXME:
     master.datalist = list of data blocks
     master.db1 ... master.dbN = selected datablocks of data to fit

     master.dbK.oidb  = link to OI-FITS datablock
     master.dbK.oiref = indices of data in OI-FITS
     master.dbK.ops = table of operations
     master.dbK.sgn = sign for frequencies
     master.dbK.idx = index for frequencies in model
     master.dbK.*


     1. selection of target object

     master = mira_new(data, ..., target=..., wavemin=..., wavemax=...,
                       freqmin=..., freqmax=..., xform=..., dim=...,
                       pixelsize=...);
     mira_save, master, file; // save model to file
     mira_config, master, wavemin=..., wavemax=..., freqmin=..., freqmax=...,
                  xform=...,  dim=..., pixelsize=...;
     mira_update, master; // automatically done
     mira_apply(master, x) // apply direct transform to image x
     mira_apply(master, x, adj) // apply direct transform or its adjoint to x
     mira_cost(master, x) // yield objective function for image x
     mira_cost_and_gradient(master, x, g) // yield objective function and its
                                          // gradient for image x
 */

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func mira_default_image(master)
/* DOCUMENT mira_default_image(master);
     Yields the default image for MiRA instance `master` which is a centered
     point-like source.

   SEE ALSO: mira_new.
 */
{
  img = array(double, mira_image_size(master));
  img(mira_image_size(master, 1)/2 + 1,
      mira_image_size(master, 2)/2 + 1) = 1.0;
  return img;
}

/*--------------------------------------------------------------------------*/
/* TABLES OF VIRTUAL OPERATIONS */

/* Instanciate the tables of virtual operations for the different data blocks.
   If the objects are already hash tables just refresh their contents. */

mira_define_table, _MIRA_VIS_OPS,
  class = "vis",
  init = symlink_to_variable(_mira_vis_init),
  cost = symlink_to_variable(_mira_vis_cost),
  ndata = symlink_to_variable(_mira_vis_ndata);

mira_define_table, _MIRA_VISAMP_OPS,
  class = "visamp",
  init = symlink_to_variable(_mira_visamp_init),
  cost = symlink_to_variable(_mira_visamp_cost),
  ndata = symlink_to_variable(_mira_visamp_ndata);

mira_define_table, _MIRA_VISPHI_OPS,
  class = "visphi",
  init = symlink_to_variable(_mira_visphi_init),
  cost = symlink_to_variable(_mira_visphi_cost),
  ndata = symlink_to_variable(_mira_visphi_ndata);

mira_define_table, _MIRA_VIS2_OPS,
  class = "vis2",
  init = symlink_to_variable(_mira_vis2_init),
  cost = symlink_to_variable(_mira_vis2_cost),
  ndata = symlink_to_variable(_mira_vis2_ndata);

mira_define_table, _MIRA_T3_OPS,
  class = "t3",
  init = symlink_to_variable(_mira_t3_init),
  cost = symlink_to_variable(_mira_t3_cost),
  ndata = symlink_to_variable(_mira_t3_ndata);

mira_define_table, _MIRA_T3AMP_OPS,
  class = "t3amp",
  init = symlink_to_variable(_mira_t3amp_init),
  cost = symlink_to_variable(_mira_t3amp_cost),
  ndata = symlink_to_variable(_mira_t3amp_ndata);

mira_define_table, _MIRA_T3PHI_OPS,
  class = "t3phi",
  init = symlink_to_variable(_mira_t3phi_init),
  cost = symlink_to_variable(_mira_t3phi_cost),
  ndata = symlink_to_variable(_mira_t3phi_ndata);


