/*
 * mira2_plugin_central_star.i -
 *
 * Plugin to add a central star of given relative flux in MiRA image model.
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

_CENTRAL_STAR_OPTS = _lst("\nCentral Star plugin options:",
                          _lst("central_star_flux", 0, "VALUE",
                               OPT_REAL,
                               "Relative brightness of the central star"));

func mira_plugin_central_star_init(nil) {
  inform, "Loading \"Central Star\" plugin...";
  return mira_new_plugin(options = _CENTRAL_STAR_OPTS,
                         parse_options = central_star_parse_options,
                         tweak_visibilities = central_star_tweak_visibilities,
                         add_keywords = central_star_add_keywords);
}

func central_star_parse_options(plugin, opt)
{
  /* Force some options. */
  if (is_void(opt.normalization)) {
    h_set, opt, normalization = 1.0;
  }
  norm = double(opt.normalization);

  /* Check plugin specific options. */
  flux = double(opt.central_star_flux);
  if (flux < 0 || flux > 1) {
    throw, "Value of `-central_star_flux` is out of range";
  }
  h_set, plugin, flux = flux,
    a0 = norm*flux,
    a1 = norm*(1 - flux);
}

func central_star_tweak_visibilities(master, vis)
{
  plugin = mira_plugin(master);
  if (plugin.flux > 0) {
    vis = plugin.a0 + plugin.a1*unref(vis);
  }
  return vis;
}

func central_star_tweak_gradient(master, grd)
{
  plugin = mira_plugin(master);
  if (plugin.flux > 0) {
    grd = plugin.a1*unref(grd);
  }
  return grd;
}

func central_star_add_keywords(master, fh)
{
  plugin = mira_plugin(master);
  fits_set, fh, "CENSTAR", plugin.flux,
    "Relative brightness of the central star";
}
