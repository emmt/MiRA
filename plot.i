/*
 * plot.i --
 *
 * Additional routines for plotting in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 1996-2010, Eric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 * This file is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * This file is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 *-----------------------------------------------------------------------------
 *
 * Routines:
 *	color_bar: add a color bar to a plot;
 *	pl3dj: plot disjoint lines in 3D space;
 *	pl3s: plot 3D surface;
 *	pl3t: plot text in 3D space;
 *	pl_box: draw a rectangle;
 *	pl_cbar: add a color bar to a plot;
 *	pl_cbox: draw a centered rectangle;
 *	pl_circle: draw a circle;
 *	pl_ellipse: draw an ellipse;
 *	pl_get_axis_flags: parse axis settings;
 *	pl_get_color: parse color value/name;
 *	pl_get_font: parse font value/name;
 *	pl_get_symbol: parse symbol value/name;
 *	pl_hinton: plot Hinton diagram;
 *	pl_img: plot an image with many options;
 *	pl_map: apply a function to all the elements of an array;
 *	pla: plot several curves at the same time;
 *	plh: plot in an "histogram" style;
 *	plhline: plot horizontal line across viewport;
 *	plp: plot points/error bars;
 *	pls: plot surface as a filled mesh with contours;
 *	plvline: plot vertical line across viewport;
 *	ps2png, ps2jpeg: convert PostScript file into bitmap image;
 *	win2png, win2jpeg: dump graphical window into bitmap image;
 *	win_copy_lim: copy plot limits between graphic windows;
 *	xbtn_plot: plot pseudo-GUI buttons;
 *	xbtn_which: manage pseudo-GUI buttons;
 *	xmouse: extended mouse interface;
 *	xmouse_box: interactively select and draw a rectangular region;
 *	xmouse_demo: simple demonstration of xmouse capabilities;
 *	xmouse_length: interactively measure a distance;
 *	xmouse_line: interactively select and draw a line;
 *	xmouse_point: interactively select and draw point(s);
 *	xwindow: setup viewport and graphic style of graphic windows;
 *
 * History:
 *	$Id: plot.i,v 1.27 2012/06/11 16:15:38 eric Exp eric $
 *	$Log: plot.i,v $
 *	Revision 1.27  2012/06/11 16:15:38  eric
 *	 - New routine pl_hinton to plot Hinton diagram.
 *
 *	Revision 1.26  2010/08/19 14:16:00  eric
 *	 - New function plfg to plot filled graph.
 *
 *	Revision 1.25  2010/06/10 14:28:50  eric
 *	 - Changed default values for PIXELBIAS and PIXELREF in pl_img.
 *
 *	Revision 1.24  2010/03/08 09:29:33  eric
 *	 - Bug in win_copy_lim fixed (thanks to Michel Tallon).
 *
 *	Revision 1.23  2010/03/05 16:23:58  eric
 *	 - Bug in win_copy_lim fixed (thanks to Isabelle Tallon-Bosc).
 *	 - win_copy_lim can optionally copy the limits of just one of the axis.
 *
 *	Revision 1.22  2008/09/04 07:19:30  eric
 *	 - Some documentation fixes.
 *
 *	Revision 1.21  2008/07/12 06:41:25  eric
 *	 - Added c-basic-offset for Emacs.
 *
 *	Revision 1.20  2008/07/11 20:43:17  eric
 *	 - New routines: pl_img and pl_cbar.
 *
 *	Revision 1.19  2007/07/25 17:02:19  eric
 *	 - Some changes in pl_get_palette.
 *
 *	Revision 1.18  2007/07/02 22:16:42  eric
 *	 - Added "utils.i" as a required file.
 *
 *	Revision 1.17  2007/07/02 22:15:18  eric
 *	 - New routines: ps2png, ps2jpeg, win2png and win2jpeg.
 *
 *	Revision 1.16  2007/05/31 07:14:09  eric
 *	 - Fixed pl_get_color when argument is already a RGB triplet.
 *
 *	Revision 1.15  2007/03/19 16:52:48  eric
 *	 - Automatically include "style.i" is function `xwindow`.
 *
 *	Revision 1.14  2007/02/06 12:05:18  eric
 *	Heavy changes in pl_get_color:
 *	 a. allow for named colors from the X11 RGB database (more than 600
 *	    color names);
 *	 b. always return a color as a char array (either indexed color
 *	    or [R,G,B] triplet); this should fix a bug which prevent to
 *	    choose the inside color in plp;
 *	 c. allow for multiple colors (array of colors);
 *	 d. allow for packed RGB colors.
 *
 *	Revision 1.13  2007/02/01 10:20:16  eric
 *	Quick reference for all routines added in top comment of this file.
 *
 *	Revision 1.12  2007/02/01 10:09:22  eric
 *	There are a lot of changes this time:
 *
 *	 - 'plp' improved so that it is possible to draw 'transparent'
 *	   symbols and to specify a different color for the inside and
 *	   outline of the symbols;
 *
 *	 - new routines: 'pl_get_symbol', 'pl_get_color', 'pl_get_font',
 *	   'pl_get_axis_flags' to parse and check values of standard graphic
 *	   keywords;
 *
 *	 - new routine: 'xwindow' to manage graphic windows and setup
 *	   viewport and graphic style;
 *
 *	 - new routines: 'xbtn_plot', 'xbtn_which' to create and manage
 *	   GUI buttons in graphic windows;
 *
 *	 - new interactive routines: 'xmouse', 'xmouse_point', 'xmouse_box',
 *	   'xmouse_line', 'xmouse_length' and 'xmouse_demo';
 *
 *	 - new plotting routines to draw simple shapes: 'pl_box', 'pl_cbox',
 *	   'pl_circle', 'pl_ellipse';
 *
 *	 - new utility routine: 'pl_map' to apply a function to all the
 *	   elements of an array.
 *
 *	Revision 1.11  2007/01/16 08:15:01  eric
 *	 - Fixed coding style and (some) documentation.
 *	 - New functions plhline and plvline.
 *
 *	Revision 1.10  2002/05/14 10:03:49  eric
 *	 - plp: use directly plfp instead of plmk to plot symbols and ticks.
 *	 - plp: draw each symbol inside a unit circle so that they look the
 *	   same size.
 *	 - plp: new "star" symbol (SYMBOL=8) and draw polygon for SYMBOL>=9.
 *	 - plp: new keywords XLO, XHI, YLO, YHI to draw non-symmetrical error
 *	   bars and FILL to draw filled symbols.
 *	 - plp: removed unused keyword HIDE.
 *
 *	Revision 1.9  2001/11/15 09:46:19  eric
 *	 - Mock check-in to synchronize version numbers because RCS
 *	   master file was lost.
 *
 *	Revision 1.8  2001/06/21 12:45:17  eric
 *	- fix indentation and floating point representation.
 *
 *	Revision 1.7  1997/04/03 15:52:39  eric
 *	Added routines color_bar and pl_fc.
 *
 *	Revision 1.6  1996/11/22 11:42:02  eric
 *	 - Add keyword ASPECT for plp.
 *
 *	02/20/96 Eric THIEBAUT: plp;
 *	02/20/96 release 1.1.
 *	02/21/96 Eric THIEBAUT: plh;
 *	02/21/96 release 1.2.
 *	02/21/96 Christophe PICHON and Eric THIEBAUT: pla;
 *	03/04/96 Eric THIEBAUT (from routines written by Christophe PICHON
 *	  and David MUNRO): pl3s, pl3dj pl3t;
 *	03/04/96 release 1.3.
 *	03/04/96 Eric THIEBAUT: added more symbols in plp;
 *	24/03/96 Eric THIEBAUT: added "box" keyword in pl3s and fix some bugs;
 *	24/03/96 release 1.4.
 *	May 3, 1996 by Eric THIEBAUT: added point symbol in routine plp;
 *
 *-----------------------------------------------------------------------------
 */

require, "utils.i";

func pl_arrow(x0, y0, x1, y1, head=, size=, width=, color=, angle=)
/* DOCUMENT pl_arrow, x0, y0, x1, y1;

     Plot an arrow starting at (X0,Y0), ending at (X1,Y1).


   RESTRICTIONS
     Works better with square limits.
     
   SEE ALSO:
 */
{
  color = pl_get_color(color);
  pldj, x0,y0,x1,y1, width=width, color=color, type=type;
  if (head) {
    if (is_void(angle)) angle = 50.0;
    if (is_void(size)) size = 1.0;
    u = x1 - x0;
    v = y1 - y0;
    s = size*0.01/abs(u, v);
    alpha = angle*(pi/360.0); // half angle with respect to direction
    sn = s*sin(alpha);
    cs = s*cos(alpha);
    xm = [cs*u - sn*v, 0.0, cs*u + sn*v, 0.0];
    ym = [cs*v + sn*u, 0.0, cs*v - sn*u, 0.0];
    m = numberof(xm);
    n = array(1, 1 + numberof(y1));
    n(1) = m;
    if ((head&1) != 0) {
      plfp, z, grow(ym,y0(*)), grow(xm,x0(*)), n,
        edges=1, ewidth=width, ecolor=color;
    }
    if ((head&2) != 0) {
      plfp, z, grow(-ym,y1(*)), grow(-xm,x1(*)), n,
        edges=1, ewidth=width, ecolor=color;
    }
  }
}

func pl_fc(z, y, x, ireg, levs=, legend=, hide=, type=, width=, color=,
	   colors=, smooth=, marks=, marker=, mspace=, mphase=,
	   triangle=, region=)
{
  d = dimsof(z);
  if (d(1) != 2)
    error, "expecting a 2D array for Z";
  nx = d(2);
  ny = d(3);
  if (is_void(x))
    x = span(1, nx, nx)(,-:1:ny);
  else if (dimsof(x)(1) == 1)
    x = x(,-:1:ny);
  if (is_void(y))
    y = span(1, ny, ny)(-:1:nx,);
  else if (dimsof(y)(1) == 1)
    y = y(-:1:nx,);
  if (is_void(levs)) {
    zmin = min(z);
    zmax = max(z);
    if (zmin >= zmax) levs = zmin;
    else levs = zmin+indgen(9)*0.1*(zmax-zmin);
  }
  plfc, z, y, x, ireg, levs=levs, colors=colors,
    triangle=triangle, region=region;
  plc, z, y, x, ireg, levs=levs, legend=legend, hide=hide, type=type,
    width=width, color=color, smooth=smooth, marks=marks, marker=marker,
    mspace=mspace, mphase=mphase, triangle=triangle, region=region;
}

/*---------------------------------------------------------------------------*/

func pl_img(img, clear=, cmin=, cmax=, zformat=, keeplimits=,
            normalize=, pixelbias=, pixelsize=, pixelref=, pixelunits=)
/* DOCUMENT pl_img, img;
 *
 *   Plot image IMG in current graphics window.
 *
 *   Keyword CLEAR can be used to call fma command (which see): if CLEAR >
 *   0, fma is called prior to drawing the image; if CLEAR < 0, fma is
 *   called after drawing the image (useful in animate mode); if CLEAR = 0
 *   or undefined, then fma is not called.
 *
 *   Keyword ZFORMAT can be used to specify the format for the color bar
 *   labels (see pl_cbar).
 *
 *   Keywords PIXELSIZE, PIXELBIAS and PIXELREF can be set to specify the
 *   pixel coordinates.  The coordinate of the center of j-th pixel is
 *   given by:
 *
 *       PIXELBIAS + PIXELSIZE*(j - PIXELREF)
 *
 *   The default settings are PIXELBIAS=0.0, PIXELSCALE=1.0 and PIXELREF=0.0
 *   (i.e. to yield the same coordinates as Yorick's indexing rules).
 *
 *   Keyword PIXELUNITS can be used to specify the label(s) of the axis.
 *
 *   PIXELSIZE, PIXELBIAS, PIXELREF and PIXELUNITS can have 1 or 2 elements
 *   to specify same or different settings for the X and Y axis.
 *
 *   Keywords CMIN and CMAX can be used to specify cut levels for the
 *   display (see pli).  If keyword NORMALIZE (see below) is true, CMIN and
 *   CMAX are in units of normalized intensity.
 *
 *   If keyword NORMALIZE is true, the flux is normalized --- divided by
 *   the pixel area that is: PIXELSIZE(1)*PIXELSIZE(0).
 *
 *
 * SEE ALSO:
 *   pli, fma, animate, pl_cbar, xytitles.
 */
{
  dims = dimsof(img);
  if (is_void(dims) || dims(1) != 2) {
    error, "expecting a 2-D image";
  }
  width = dims(2);
  height = dims(3);

  if (is_void(pixelref)) {
    pixelref = 0.0;
  }
  if (is_void(pixelbias)) {
    pixelbias = 0.0;
  }
  if (is_void(pixelsize)) {
    pixelsize = 1.0;
  }
  if (is_void(normalize)) {
    scl = 1.0;
  } else {
    scl = 1.0/(pixelsize(1)*pixelsize(0));
    if (scl != 1.0) {
      img *= scl;
    }
  }

  if (is_void(cmin)) cmin = min(img);
  if (is_void(cmax)) cmax = max(img);

  xsize = pixelsize(1);
  ysize = pixelsize(0);
  xbias = pixelbias(1);
  ybias = pixelbias(0);
  xref = pixelref(1);
  yref = pixelref(0);

  x0 = xbias + xsize*(0.5 - xref);
  x1 = xbias + xsize*(width + 0.5 - xref);
  y0 = ybias + ysize*(0.5 - yref);
  y1 = ybias + ysize*(height + 0.5 - yref);

  if (clear && clear > 0) {
    fma;
  }
  pli, img, x0, y0, x1, y1, cmin=cmin, cmax=cmax;
  local red, green, blue;
  palette, red, green, blue, query=1;
  ncolors = numberof(red);
  levs = span(cmin, cmax, ncolors + 1);
  colors = indgen(ncolors);
  pl_cbar, cmin=cmin, cmax=cmax, vert=1, nlabs=11, format=zformat;
  if (! is_void(pixelunits)) {
    xytitles, pixelunits(1), pixelunits(0);
  }
  if (clear && clear < 0) {
    fma;
  }
}

func pl_cbar(z, cmin=, cmax=, vert=, nlabs=, adjust=,
             color=, font=, height=, opaque=, orient=,
             labels=, levels=,
             width=, ticklen=, thickness=, vport=, format=)
/* DOCUMENT pl_cbar, z;
       -or- pl_cbar, cmin=CMIN, cmax=CMAX;

     Draw a color bar below the current coordinate system.  The colors and the
     default associated label values are from min(Z) to max(Z) --
     alternatively, keywords CMIN and CMAX can be specified.

     If keywords VERT is set true, the color bar appears to the left of the
     current coordinate system; by default, an horizontal color bar is drawn.

     By default, the colorbar is drawn next to the current viewport; other
     viewport coordinates can be given by VPORT=[xmin,xmax,ymin,ymax].
     Keyword ADJUST can be used to move the bar closer to (adjust<0) or
     further from (adjust>0) the viewport.

     Keyword LEVELS can be used to specify the positions of the labels and
     ticks to draw in units of Z.

     Keyword LABELS can be used to specify the labels to print along the color
     bar.  If keyword LABELS is not specified, the labels are the textual
     values of the levels and the format of the labels can be specified with
     keyword FORMAT; by default FORMAT="%.3g".  The font type, font height
     and text orientation for the labels can be set with keywords FONT
     (default FONT="helvetica"), HEIGHT (default HEIGHT=14 points) and ORIENT
     respectively.

     If neither LEVELS nor NLABS are specified, keyword NLABS can be used to
     choose the number of displayed labels; by default, NLABS=11 which
     correspond to a label every 10% of the dynamic; use NLABS=0 to suppress
     all labels.

     Keyword COLOR can be used to specify the color of the labels, the
     ticks and the frame of the colorbar.  Default is foreground color.

     Keyword WIDTH can be used to set the width of the lines used to draw
     the frame and the ticks of the colorbar.

     Keyword TICKLEN can be used to set the length (in NDC units) of the
     ticks.  Default is 0.007 NDC.

     Keyword THICKNESS can be used to set the thickness of the colorbar (in
     NDC units).  Default is 0.020 NDC.


    SEE ALSO: pl_img, pl_span, pli, plt, pldj, plg, viewport.
 */
{
  if (is_void(cmin)) cmin = min(z);
  if (is_void(cmax)) cmax = max(z);
  if (is_void(vport)) vport = viewport();
  if (is_void(adjust)) adjust = 0.0;
  if (is_void(ticklen)) ticklen = 0.007;
  if (is_void(thickness)) thickness = 0.020;
  nlevels = numberof(levels);
  nlabels = numberof(labels);
  if (nlabels != nlevels && nlevels != 0 && nlabels != 0) {
    error, "LABELS and LEVELS must have same number of elements";
  }
  n = max(nlabels, nlevels);
  if (is_void(nlabs)) {
    nlabs = (n != 0 ? n : 11);
  } else if (n != 0) {
    error, "NLABS must be unspecified when LEVELS or LABELS are given";
  }

  /* Make the color gradient image to display. */
  local red, green, blue;
  palette, red, green, blue, query=1;
  ncolors = numberof(red);
  if (ncolors < 2) {
    ncolors = 240;
  }
  cells = char(indgen(0 : ncolors - 1));

  linetype = 1; /* "solid" */

  if (vert) {
    x0 = vport(2) + adjust + 0.022;
    x1 = x0 + thickness;
    z0 = y0 = vport(3);
    z1 = y1 = vport(4);
    cells = cells(-,);
  } else {
    z0 = x0 = vport(1);
    z1 = x1 = vport(2);
    y0 = vport(3) - adjust - 0.045;
    y1 = y0 - thickness;
    cells = cells(,-);
  }
  sys = plsys(0);
  pli, cells, x0, y0, x1, y1;
  if (is_void(width) || width != 0) {
    plg, [y0,y0,y1,y1], [x0,x1,x1,x0], closed=1,
      color=color, width=width, type=linetype, marks=0;
  }

  if (nlabs > 0) {
    cmin += 0.0; // make sure CMIN is floating point
    if (is_void(levels)) {
      if (cmin == cmax) {
        levels = cmin;
        lz = (z0 + z1)/2.0;
        if (nlabels > 1) labels = labels((nlabels + 1)/2);
        nlabs = 1;
      } else if (nlabs == 1) {
        levels = (cmin + cmax)/2.0;
        lz = (z0 + z1)/2.0;
      } else {
        levels = pl_span(cmin, cmax, nlabs);
        lz = pl_span(z0, z1, nlabs);
      }
    } else {
      /* Select levels inside the displayed range. */
      j = where((levels >= min(cmin, cmax))&(levels <= max(cmin, cmax)));
      if ((n = numberof(j)) != nlabs) {
        if (is_array(j)) {
          levels = levels(j);
          if (! is_void(labels)) labels = labels(j);
        }
        nlabs = n; /* may be zero here */
      }
      if (nlabs > 0) {
        if (cmin == cmax) {
          levels = cmin;
          lz = (z0 + z1)/2.0;
          if ((n = numberof(labels)) > 1) labels = labels((n + 1)/2);
          nlabs = 1;
        } else {
          lz = interp([z0, z1], [cmin, cmax], levels);
        }
      }
    }
  }

  if (nlabs > 0) {
    if (is_void(labels)) {
      if (is_void(format)) format = "%.3g";
      labels = swrite(format=format, levels);
    }

    local lx0, lx1, lx2, ly0, ly1, ly2;
    if (vert) {
      lx0 = array(x1, nlabs);
      lx1 = array(x1 + ticklen, nlabs);
      lx2 = array(x1 + 1.67*ticklen, nlabs);
      eq_nocopy, ly0, lz;
      eq_nocopy, ly1, lz;
      eq_nocopy, ly2, lz;
      justify = "LH";
    } else {
      ly0 = array(y1, nlabs);
      ly1 = array(y1 - ticklen, nlabs);
      ly2 = array(y1 - 1.67*ticklen, nlabs);
      eq_nocopy, lx0, lz;
      eq_nocopy, lx1, lz;
      eq_nocopy, lx2, lz;
      justify = "CT";
    }
    if (ticklen && (is_void(width) || width != 0)) {
      pldj, lx0, ly0, lx1, ly1,
        color=color, width=width, type=linetype;
    }
    for (i = 1; i <= nlabs; ++i) {
      plt, labels(i), lx2(i), ly2(i), tosys=0, color=color, font=font,
        height=height, opaque=opaque, orient=orient, justify=justify;
    }
  }
  plsys, sys;
}

func color_bar(levs, colors, vert=, labs=, adjust=, color=, width=,
               height=, ticklen=, vport=, format=, font=)
/* DOCUMENT color_bar;
         or color_bar, levs, colors;

     Note: pl_cbar (which see) is probably a better routine.

     Draw a color bar below the current coordinate system.  If LEVS is not
     specified uses plfc_levs (set by previous call to plfc).  If COLORS is
     specified, it should have one more value than LEVS, otherwise equally
     spaced colors are chosen, or plfc_colors if plfc_levs was used.  With
     the VERT=1 keyword the color bar appears to the left of the current
     coordinate system (vert=0 is default).  By default, color_bar will
     attempt to label some of the color interfaces.  With the LABS keyword,
     you can force the labelling algorithm as follows: LABS=0 supresses all
     labels, LABS=n forces a label at every n-th interface, LABS=[i,n]
     forces a label every n-th interface starting at i-th interface
     (0<=i<=numberof(LEVS)).

     You can specify the viewport coordinates by keyword
     VPORT=[xmin,xmax,ymin,ymax]; by default the colorbar is drawn next to
     the current viewport.  You can use the ADJUST keyword to move the bar
     closer to (adjust<0) or further from (adjust>0) the viewport.

     You can specify the string format for labels with keyword FORMAT
     (default "%g"), the font type with keyword FONT (default "helvetica")
     and the font height with keyword HEIGHT (default 14 points).

     Keyword COLOR can be used to specify the color of the labels, the
     ticks and the frame of the colorbar.  Default is foreground color.

     Keyword WIDTH can be used to set the width of the lines used to draw
     the frame and the ticks of the colorbar.

     Keyword TICKLEN can be used to set the length (in NDC units) of the
     ticks.  Default is 0.005 NDC.

   SEE ALSO: pl_cbar, plfc. */
{
  nil = string(0);
  if (is_void(levs)) {
    if (is_void(plfc_levs)) error, "no levels specified";
    levs = plfc_levs;
    n = numberof(levs)+1;
    if (is_void(colors)) colors = plfc_colors;
  } else {
    n = numberof(levs) + 1;
    if (is_void(colors)) colors = bytscl(span(1,n,n),cmin=0.5,cmax=n+0.5);
  }
  if (n != numberof(colors))
    error, "numberof(colors) must be one more than numberof(levs)";

  if (is_void(vport)) vport = viewport();
  if (is_void(adjust)) adjust = 0.0;
  if (is_void(ticklen)) ticklen = 0.005;
  dx = dy = 0.0;
  if (vert) {
    x = (vport(2)+adjust+[0.022,0.042])(-:1:n+1,);
    dx = ticklen;
    y = span(vport(3),vport(4),n+1)(,-:1:2);
  } else {
    y = (vport(3)-adjust-[0.045,0.065])(-:1:n+1,);
    dy = -ticklen;
    x = span(vport(1),vport(2),n+1)(,-:1:2);
  }
  sys = plsys(0);
  plf,[colors],y,x,edges=1,ecolor=color,legend=nil;
  plsys, sys;

  if (is_void(labs) || labs(0) > 0) {
    if (numberof(levs) > 1) {
      dz = levs(dif);
      if (numberof(dz) != numberof(levs) - 1 ||
          anyof((dz > 0.0) != (dz(1) > 0.0)) || !dz(1))
        error, "levs must be monotone 1D";
      levs = levs(1:0);
      levs = grow([2*levs(1)-levs(2)],levs,[2*levs(0)-levs(-1)]);
    } else {
      levs = double(levs(1));
      if (!levs) levs = [-1.0,levs,1.0];
      else levs = [0.0,levs,2*levs];
    }
    if (numberof(labs)<2) {
      if (is_void(labs)) labs = (n-1)/4 + 1;
      orig = where(levs<1.0e-9*max(levs(dif)));
      if (numberof(orig)==1) labs = [orig(1)%labs,labs];
      else labs = [(n%labs)/2,labs];
    }
    list = where(indgen(0:n)%labs(2)==labs(1));
    x = x(list,);
    y = y(list,);
    if (is_void(format)) format = "%g";
    labs = swrite(format=format,levs(list));
    plsys, 0;
    pldj, x(,2),y(,2),x(,2)+dx,y(,2)+dy, legend=nil, color=color, width=width;
    plsys, sys;
    if (is_void(font)) font = "helvetica";
    plt1, labs,x(,2)+2*dx,y(,2)+2*dy, justify=(vert?"LH":"CT"),
      height=height, font=font, color=color;
  }
}

/*---------------------------------------------------------------------------*/

func pla(y, x, every=, legend=, hide=, type=, width=, color=, closed=, smooth=,
	 marks=, marker=, mspace=, mphase=, rays=, arrowl=, arroww=, rspace=,
	 rphase=)
/* DOCUMENT pla, y, x
         or pla, y

     Plot the buddle of curves Y versus X labelled by their last dimension.
     Y must be 2-dimensional, and X may be 2-dimensional, 1-dimensional or
     omitted.  If X is 2-dimensional, it must have the same dimensions as Y
     and Y(,i) versus X(,i) is plotted for each last indice i.  If X is
     1-dimensional, it must have the same length as the 1st dimension of Y
     and Y(,i) versus X is plotted for each last indice i.  If X is
     omitted, it defaults to [1, 2, ..., numberof(Y(,1))].

     The plotting keywords of plg are accepted plus the optional keyword
     EVERY=N which can be used to plot every N curves in the bundle
     (default N=1).

   EXAMPLE
     x = span(0,1,25)(,-:1:25);
     pla, x*transpose(x), marks=0, every=3;

   SEE ALSO
     plg, plp.
*/
{
  if (is_void(every)) {
    n = 1;
  } else if (! is_array(every) || dimsof(every)(1) || (n = long(every)) <= 0) {
    error, "EVERY must be a scalar >= 1";
  }
  if (! is_array(y) || dimsof(y)(1) != 2) {
    error, "Y must be 2-dimensional";
  }
  imax = dimsof(y)(3);
  if (is_void(x)) {
    x2d = 0N;
  } else {
    x2d = dimsof(x)(1) >= 2;
  }
  for(i = (n+1)/2; i <= imax; i += n) {
    px = (x2d ? &x(,i) : &x);
    plg, y(, i), *px, legend=legend, hide=hide, type=type, width=width,
      color=color, closed=closed, smooth=smooth, marks=marks, marker=marker,
      mspace=mspace, mphase=mphase, rays=rays, arrowl=arrowl, arroww=arroww,
      rspace=rspace, rphase=rphase;
  }
}

/*---------------------------------------------------------------------------*/

func pls_mesh(&x, &xx, d, which=, inhibit=)
/* DOCUMENT err_msg = pls_mesh(x, xx, dimsof(z), which=1/2, inhibit=1/2)

     build X and/or XX arrays of coordinates (abscissa if last argument is
     0/nil; otherwise ordinate) for 2-D array Z.  Normally, the returned
     value is string(0) otherwise it is an error message.

     X is input and output, it will have the same shape as Z and will be
     suitable for contour plots.  XX is purely output, it will have 1 more
     element than Z in each dimension and will be suitable for mesh plots.
     In other words, X(i,j) will furnish the coordinate of the centre of
     cell Z(i,j) whereas XX(i,j), XX(i,j+1), XX(i+1,j) and XX(i+1,j+1)
     will give the coordinates of the corners of cell Z(i,j).

     Assuming the length of Z along the considered dimension is N
     (N must be >= 2) there are 3 possibilities:
       (1) if X is a vector with N elements or has the same shape as Z,
           then X is considered to give the coordinates at the centre of Z
           cells: X is unchanged and output XX is build by interpolating
	   (and extrapolating at the edges) X ;
       (2) if X is a vector with N+1 elements or has 1 more element than Z
           in each dimension, then X is considered to give the coordinates
           at the corners of Z cells: output XX is set to input X and
	   output X is build by interpolating output XX;
       (3) if X is nil, it defaults to [0.5, 1.5, ..., N-0.5] and XX
           defaults to [0, 1, ..., N] along the considered dimension.
     Finally, if X is 1-D, it is expanded in the other direction.

     If keyword WHICH is 1 (the default), abscissa is the dimension of
     interest; otherwise WHICH must be 2 and ordinate is the dimension
     of interest.

     If keyword INHIBIT is 1, then only X output is computed; if INHIBIT
     is 2 then only XX output is computed.

   SEE ALSO: pls, pl3s, plmesh.
 */
{
  xx = [];
  if (is_void(which))
    which = 1;
  do_x = inhibit != 1;
  do_xx = inhibit != 2;
  expand=1;
  if (d(1) != 2 || anyof(d < 2))
    return "Z must be 2-dimensional and have at least 2-by-2 elements";
  n1 = d(2);
  n2 = d(3);
  n = d(which+1);
  if (is_void((dx = dimsof(x)))) {
    if (do_x)
      x = span(0.5, n-0.5, n);
    if (do_xx)
      xx = span(0, n, n+1);
  } else if (dx(1) == 1) {
    if (dx(2) == n) {
      if (do_xx) {
	xx = x(pcen);
	xx(1) = 2.0 * x(1) - x(2);
	xx(0) = 2.0 * x(0) - x(-1);
      }
    } else if (dx(2) == n+1) {
      xx = x;
      x = do_x ? xx(zcen) : [];
    }
  } else if (dx(1) == 2) {
    expand = 0;
    if (allof(dx == d)) {
      if (do_xx) {
	t = x(pcen,);
	t(1,) = 2.0 * x(1,) - x(2,);
	t(0,) = 2.0 * x(0,) - x(-1,);
	xx = t(,pcen);
	xx(,1) = 2.0 * t(,1) - t(,2);
	xx(,0) = 2.0 * t(,0) - t(,-1);
	t = [];
      }
    } else if (allof(dx == d + [0,1,1])) {
      xx = x;
      x = do_x ? xx(zcen,zcen) : [];
    }
  }
  if (is_void(xx) && is_void(x)) {
    return "X, Y and Z are not compatible";
  }
  if (expand) {
    if (which == 1) {
      if (do_x)
	x = x(,-:1:n2);
      if (do_xx)
	xx = xx(,-:1:n2+1);
    } else {
      if (do_x)
	x = x(-:1:n1,);
      if (do_xx)
	xx = xx(-:1:n1+1,);
    }
  }
  return string(0);
}

func pls(z, y, x, cbar=, viewport=, title=, xtitle=, ytitle=,
	 legend=, hide=, top=, cmin=, cmax=, edges=, ecolor=, ewidth=,
         height=, font=, levs=, nlevs=, type=, width=, color=,
         marks=, marker=, mspace=, mphase=, smooth=)
/* DOCUMENT pls, z, y, x
         or pls, z
     draws surface plot of Z versus (X,Y) as a filled mesh with
     optional contours.  The Z array must be a 2-dimensional array,
     see documentation of pls_mesh for the meaning of X and Y.

     If keyword CBAR is set to non-zero, a color bar is drawn on the
     right of the plot.  The current viewport (in NDC) may be
     specified with keyword VIEWPORT, default is:
       [0.19, 0.60, 0.44, 0.85].

     The appearance of the filled mesh can be modified by means of
     keywords: LEGEND, HIDE, TOP, CMIN, CMAX, EDGES, ECOLOR and EWIDTH
     (see plf documentation).

     Optional contour plot of Z may be superimposed by either keyword
     NLEVS to set the number of contours or by with keyword LEVS to
     specify the level values.  The appearance of the contour plot can
     be modified by means of keywords: LEGEND, HIDE, TYPE, WIDTH,
     COLOR, MARKS, MARKER, MSPACE, MPHASE and SMOOTH (see plc
     documentation).

   SEE ALSO: pls_mesh, pl3s, plc, plf, plmesh.
*/
{
  local r, g, b;	// these variables are used to query colors
  local xx, yy;

  /*
   * Set some defaults.
   */
  if (is_void(edges))
    edges = 0;
  if (is_void(height)) {
    height = 12;
    small=10;
  } else {
    s = [8,10,12,14,18,24];
    i = where(height == s);
    if (numberof(i) != 1)
      error, "bad font HEIGHT";
    i = i(1);
    small = i > 1 ? s(i-1) : height;
  }
  if (numberof(levs)) {
    nlevs = numberof(levs);
  } else if (is_void(nlevs)) {
      nlevs = 8;
  }

  /*
   * Compute mesh coordinates.
   */
  i = nlevs >= 1 ? 0 : 1;
  if ((msg = pls_mesh(x, xx, dimsof(z), which=1, inhibit=i)) != string(0) ||
      (msg = pls_mesh(y, yy, dimsof(z), which=2, inhibit=i)) != string(0))
    error, msg;

  /*
   * Plot color bar and titles.
   */
  vpmax = [0.127, 0.672, 0.363, 0.908];
  if (numberof(viewport) != 4)
    viewport = [0.19, 0.60, 0.44, 0.85];		// standard viewport
  if (cbar) {
    local r, g, b;
    plsys, 0;
    margin = vpmax(2)-viewport(2);
    x0 = viewport(2) + 0.7 * margin;
    x1 = viewport(2) + 0.9 * margin;
    y0 = viewport(3);
    y1 = viewport(4);
    palette, r, g, b, query=1;
    n = numberof(r);
    r = g = b = [];
    pli, char(indgen(n)-1)(-,), legend=string(0), x0, y0, x1, y1;
    plg, [y0,y0,y1,y1,y0], [x0,x1,x1,x0,x0], legend=string(0), marks=0,
      width=1;
    plsys, 1;
  }
  xc = 0.5*(viewport(1)+viewport(2));
  yc = 0.5*(viewport(3)+viewport(4));
  if (!is_void(title)) {
    plt, title, xc, viewport(4) + 0.9 * (vpmax(4) - viewport(4)), tosys=0,
      legend=string(0), justify="CT", path=0,
      font=font, height=height, opaque=opaque;
  }
  if (!is_void(xtitle)) {
    plt, xtitle, xc, vpmax(3) + 0.05 * (viewport(3) - vpmax(3)), tosys=0,
      legend=string(0), justify="CB", path=0,
      font=font, height=small, opaque=opaque;
  }
  if (!is_void(ytitle)) {
    plt, ytitle, vpmax(1) + 0.05 * (viewport(1) - vpmax(1)), yc, tosys=0,
      legend=string(0), justify="LH", path=1,
      font=font, height=small, opaque=opaque;
  }

  /*
   * Plot filled mesh.
   */
  plf, z, yy, xx, legend=legend, hide=hide,
    top=top, cmin=cmin, cmax=cmax, edges=edges, ecolor=ecolor, ewidth=ewidth;
  xx = yy = [];

  /*
   * Plot contours.
   */
  if (nlevs) {
    if (is_void(levs)) {
      zmax = double(max(z));
      zmin = double(min(z));
      levs = zmin + (zmax-zmin) / double(nlevs+1) * indgen(nlevs);
    }
    plc, z, y, x, levs=levs, legend=legend, hide=hide, type=type, width=width,
      color=color, marks=marks, marker=marker, mspace=mspace, mphase=mphase,
      smooth=smooth;
  }
}

/*---------------------------------------------------------------------------*/

func plh(y, x, just=, legend=, hide=, type=, width=, color=, marks=, marker=,
	 mspace=, mphase=)
/* DOCUMENT plh, y, x
         or plh, y
	plots a graph of Y versus X in an "histogram" style (i.e., with
	steps).  Y and X must be 1-D arrays of equal length; if X is
	omitted, it defaults to [1, 2, ..., numberof(Y)].

	The optional keyword JUST set justification of the histogram:
	JUST=1, 2 or 3 makes the graph be left justified, centered or
	right justified respectively along X axis.  Default is centered.

	Other plotting keywords (legend, hide, type, width, color, marks,
	marker, mspace, and mphase) are passed to the plg routine.

   SEE ALSO: plg, plm, plc, plv, plf, pli, plt, pldj, plfp
             limits, logxy, range, fma, hcp
*/
{
  /* parse/check arguments */
  if (!is_array(y) || dimsof(y)(1)!=1 || (n = numberof(y)) < 2)
    error, "Y must be a vector of at least 2 elements";
  if (is_void(x))
    x = double(indgen(numberof(y)));
  else if (!is_array(x) || dimsof(x)(1)!=1 || numberof(x) != n)
    error, "X must be a vector of same length as Y";
  if (is_void(just))
    just = 2;
  if (! is_void(color)) color = pl_get_color(color);

  /* build new X vector */
  n2 = 2 * n;
  x2 = array(double, n2);
  if (just == 1) {
    /* left justify */
    x2(1::2) = x;
    x2(2:-1:2) = x(2:);
    x2(0) = 2 * x(0) - x(-1);
  } else if (just == 2) {
    /* center */
    d = 0.5 * x(dif);
    dx = d(1);
    grow, dx, d, d(0);
    d = [];
    x2(1::2) = x - dx(:-1);
    x2(2::2) = x + dx(2:);
    dx = [];
  } else if (just == 3) {
    /* right justify */
    x2(1) = 2 * x(1) - x(2);
    x2(2::2) = x;
    x2(3::2) = x(:-1);
  } else {
    error, "bad value for JUST";
  }

  /* build new Y vector */
  y2 = array(double, n2);
  y2(1::2) = y2(2::2) = y;

  /* plot the graph */
  plg, y2, x2,
    legend=legend, hide=hide, type=type, width=width, color=color,
    marks=marks, marker=marker, mspace=mspace, mphase=mphase;
}

func plfg(y, x, base=, vert=, color=, edges=, ecolor=, ewidth=, etype=,
          legend=, hide=)
/* DOCUMENT plfg, y;
         or plfg, y, x;

     This subroutine plots a filled curve of Y versus X.  Y and X must be 1-D
     arrays of equal length; if X is omitted, it defaults to [1, 2, ...,
     numberof(Y)].  The filled area is delimited by the curve and a baseline
     which is horizontal unless keyword VERT is set true, to indicate that the
     baseline is vertical.  Keyword BASE can be used to specify the coordinate
     of the baseline (BASE = 0.0 by default).  The fill color (default
     foreground) can be specified with keyword COLOR.

     If keyword EDGES is true, the curve is also drawn.  The color, type and
     width of the line used to draw the curve can be specified by keywords
     ECOLOR, ETYPE and EWIDTH.  If keyword EDGES is set with the special value
     -1, the full outline of the filled area is drawn; otherwise only the
     curve itself is drawn.  By default, the curve is not drawn unless at
     least one of the keywords ECOLOR, ETYPE or EWIDTHis non nil.

   SEE ALSO: plg, plfp, pl_get_color.
 */
{
  n = numberof(y);
  if (is_void(base)) base = 0.0;
  if (is_void(x)) x = double(indgen(n));
  if (is_void(color)) {
    z = []; // will draw with background color
  } else {
    z = pl_get_color(color)(,-);
  }
  if (! is_void(ecolor)) {
    ecolor = pl_get_color(ecolor);
  }
  if (is_void(edges)) {
    edges = ! (is_void(ecolor) && is_void(ewidth) && is_void(etype));
  }
  local ex, ey;
  if (! vert) {
    x1 = x(0);
    x2 = x3 = x(1);
    y1 = y2 = base;
    y3 = y(1);
  } else {
    x1 = x2 = base;
    x3 = x(1);
    y1 = y(0);
    y2 = y3 = y(1);
  }
  ey = grow(y, y1, y2, y3);
  ex = grow(x, x1, x2, x3);
  plfp, z, ey, ex, n + 3, edges=0, legend=legend, hide=hide;
  if (edges) {
    if (edges != -1) {
      eq_nocopy, ex, x;
      eq_nocopy, ey, y;
    }
    plg, ey, ex, color=ecolor, width=ewidth, type=etype,
      legend=legend, hide=hide, marks=0;
  }
}

/*---------------------------------------------------------------------------*/

func plhline(y, x0, x1, color=, width=, type=) {
  lim = limits(); one = array(1.0, dimsof(y));
  pldj, one*(is_void(x0) ? lim(1) : x0), y, one*(is_void(x1) ? lim(2) : x1), y,
    color=color, width=width, type=type; }
func plvline(x, y0, y1, color=, width=, type=) {
  lim = limits(); one = array(1.0, dimsof(x));
  pldj, x, one*(is_void(y0) ? lim(3) : y0), x, one*(is_void(y1) ? lim(4) : y1),
    color=color, width=width, type=type; }
/* DOCUMENT plhline, y;
       -or- plhline, y, x0, x1;
       -or- plvline, x;
       -or- plhline, x, y0, y1;
     Plots an horizontal/vertical line.

   KEYWORDS color, type, width.
   SEE ALSO pldj. */

/*---------------------------------------------------------------------------*/
/* PLP - PLOTTING POINTS */

func plp(y, x, dx=, xlo=, xhi=, dy=, ylo=, yhi=, size=, symbol=, ticks=,
	 legend=, type=, width=, color=, fill=)
/* DOCUMENT plp, y, x;
 *     -or- plp, y, x, dx=sigma_x, dy=sigma_y;
 *
 *   Plots points (X,Y) with symbols and/or  error bars.  X, and Y may have
 *   any dimensionality, but  must have the same number  of elements.  If X
 *   is nil, it defaults to indgen(numberof(Y)).
 *
 *   Keyword SYMBOL  may be used  to choose the  shape of each  symbol (see
 *   pl_get_symbol).
 *
 *   Keyword SIZE  may be used to change  the size of the  symbols and tick
 *   marks (SIZE acts as a multiplier, default value is 1.0).
 *
 *   If value of  keyword FILL is true (non-nil  and non-zero), symbols are
 *   filled  with COLOR.  The  default is  to draw  open symbols.   You can
 *   specify different colors for the  outline and the inside of the symbol
 *   by using keyword COLOR = [EDGE_COLOR, FILL_COLOR].
 *
 *   Keywords XLO, XHI, YLO, and/or YHI  can be used to indicate the bounds
 *   of the optional  error bars (default is to draw  no error bars).  Only
 *   specified bounds get plotted as  error bars. If value of keyword TICKS
 *   is true  (non-nil and non-zero), ticks  get drawn at  the endpoints of
 *   the error bars.   Alternatively, keywords DX and/or DY  can be used to
 *   plot  error bars  as segments  from XLO=X-DX  to XHI=X+DX  and/or from
 *   YLO=Y-DY to  YHI=Y+DY.  If keyword  DX (respectively DY) is  used, any
 *   value of XLO and XHI (respectively YLO and YHI) is ignored.
 *
 *   The  other keywords are  the same  as for  pldj: LEGEND,  TYPE, WIDTH,
 *   COLOR (TYPE is only used to draw error bars).
 *
 *
 * EXAMPLES:
 *   plp, y, x, symbol='*', color=["orange","DarkSlateBlue"], fill=1,
 *        size=3, width=5;
 *
 *
 * KEYWORDS:
 *   legend, type, width, color.
 *
 *
 * SEE ALSO:
 *   pl_get_symbol, pl_get_color,
 *   pldj, plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmk,
 *   limits, logxy, range, fma, hcp.
 */
{
  /* NDC units for symbols/ticks (one pixel = 0.00125268 NDC at 75 DPI) */
  u0 = 0.0;       // zero
  u1 = 0.0100214; // radius of about 8 pixels at 75 DPI
  if (! is_void(size)) u1 *= size;

  /* parse color(s) */
  if (is_array(color) && dimsof(color)(0) == 2) {
    edge_color = pl_get_color(color(..,1));
    if (fill) {
      fill_color = pl_get_color(color(..,2));
    }
  } else {
    edge_color = pl_get_color(color);
    if (fill) {
      fill_color = edge_color;
    }
  }

  /* default X */
  if (is_void(x)) (x = array(double, dimsof(y)))(*) = indgen(numberof(y));

  /* error bars */
  if (is_void(dx)) {
    err = (! is_void(xlo)) + 2*(! is_void(xhi));
  } else {
    xlo = x - dx;
    xhi = x + dx;
    err = 3;
  }
  if (err) {
    pldj, (is_void(xlo) ? x : xlo), y, (is_void(xhi) ? x : xhi), y,
      type=type, width=width, color=color;
    if (ticks) {
      xm = [ u0, u0];
      ym = [-u1, u1];
      if      (err == 1) __plp,   y,       xlo;
      else if (err == 2) __plp,   y,       xhi;
      else               __plp, [y, y], [xlo, xhi];
    }
    xhi = xlo = [];
  }
  if (is_void(dy)) {
    err = (! is_void(ylo)) + 2*(! is_void(yhi));
  } else {
    ylo = y - dy;
    yhi = y + dy;
    err = 3;
  }
  if (err) {
    pldj, x, (is_void(ylo) ? y : ylo), x, (is_void(yhi) ? y : yhi),
      type=type, width=width, color=color;
    if (ticks) {
      xm = [-u1, u1];
      ym = [ u0, u0];
      if      (err == 1) __plp,    ylo,       x;
      else if (err == 2) __plp,    yhi,       x;
      else               __plp, [ylo, yhi], [x, x];
    }
    yhi = ylo = [];
  }

  /* symbols */
  symbol = pl_get_symbol(symbol);
  if (! symbol) {
    return;
  }
  if (symbol == 1) {
    /* square */
    fillable = 1n;
    u2 = u1*sqrt(0.5);
    xm = [-u2, u2, u2,-u2];
    ym = [ u2, u2,-u2,-u2];
  } else if (symbol == 2) {
    /* + cross */
    fillable = 0n;
    xm = [-u1, u1, u0, u0, u0, u0];
    ym = [ u0, u0, u0, u1,-u1, u0];
    fill = 0;
  } else if (symbol == 3) {
    /* triangle */
    fillable = 1n;
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [u0, u3,-u3];
    ym = [u1,-u2,-u2];
  } else if (symbol == 4) {
    /* hexagon */
    fillable = 1n;
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [ u1, u2,-u2,-u1,-u2, u2];
    ym = [ u0, u3, u3, u0,-u3,-u3];
  } else if (symbol == 5) {
    /* diamond */
    fillable = 1n;
    xm = [u1, u0,-u1, u0];
    ym = [u0, u1, u0,-u1];
  } else if (symbol == 6) {
    /* x cross (rotated 45 degrees) */
    fillable = 0n;
    u2 = u1*sqrt(0.5);
    xm = [u2,-u2, u0, u2,-u2, u0];
    ym = [u2,-u2, u0,-u2, u2, u0];
    fill = 0;
  } else if (symbol == 7) {
    /* triangle (upside down) */
    fillable = 1n;
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [ u0, u3,-u3];
    ym = [-u1, u2, u2];
  } else if (symbol == 8) {
    /* 5 branch star
     *   C18 = cos(18*ONE_DEGREE)
     *   S18 = sin(18*ONE_DEGREE)
     *   C54 = cos(54*ONE_DEGREE)
     *   S54 = sin(54*ONE_DEGREE)
     */
    fillable = 1n;
    u2 = 0.224514*u1; // C54*S18/S54
    u3 = 0.309017*u1; // S18
    u4 = 0.951057*u1; // C18
    u5 = 0.363271*u1; // C18*S18/S54
    u6 = 0.118034*u1; // S18*S18/S54
    u7 = 0.587785*u1; // C54
    u8 = 0.809017*u1; // S54
    u9 = 0.381966*u1; // S18/S54
    xm = [ u0, u2, u4, u5, u7, u0,-u7,-u5,-u4,-u2];
    ym = [ u1, u3, u3,-u6,-u8,-u9,-u8,-u6, u3, u3];
  } else {
    /* N-side polygon in unit circle */
    fillable = 1n;
    PI = 3.141592653589793238462643383279503;
    a = (2.0*PI/symbol)*indgen(0:symbol-1);
    xm = u1*cos(a);
    ym = u1*sin(a);
  }
  __plp, y, x;
}

func __plp(y, x)
/* DOCUMENT __plp, x, y;
     Private routine used by plp. */
{
  extern xm, ym, edge_color, fill_color, fill, legend, width;
  local z;
  n = array(1, 1 + numberof(y));
  m = numberof(ym);
  if (m > 2) {
    if (fill) {
      z = array(fill_color, numberof(n));
    } else if (fillable) {
      /* Draw inside and ouside edges to emulate 'open' (transparent)
         symbols. */
      m += m;
      grow, xm, xm(1), xm(0:2:-1);
      grow, ym, ym(1), ym(0:2:-1);
    }
  }
  n(1) = m;
  plfp, z, grow(ym,y(*)), grow(xm,x(*)), n,
    legend=legend, edges=1, ewidth=width, ecolor=edge_color;
}

/*---------------------------------------------------------------------------*/
/* PARSING OF PLOTTING KEYWORDS */

local PL_BG, PL_FG, PL_BLACK, PL_WHITE, PL_RED, PL_GREEN, PL_BLUE, PL_CYAN;
local PL_MAGENTA, PL_YELLOW, PL_GRAYD, PL_GRAYC, PL_GRAYB, PL_GRAYA;
local PL_XOR, PL_EXTRA;
local _PL_COLOR_NAMES, _PL_COLOR_RGB, _PL_COLOR_PACKED;
func pl_get_color(color, flags)
/* DOCUMENT pl_get_color(color);
 *     -or- pl_get_color(color, flags);
 *
 *   The function pl_get_color parses its argument and returns the
 *   corresponding numerical color(s). Input COLOR can be specified by
 *   color name, color index, packed RGB triplet or non-packed [R,G,B]
 *   triplet.  In any case, the output is an array of char's and either
 *   indexed color(s) or [R,G,B] triplet(s).
 *
 *   The default is to accept only a single input color and to return
 *   the corresponding color index or [R,G,B] triplet.  Optional argument
 *   FLAGS can be specified to change this default behaviour:
 *     FLAGS  = 1 - to force to RGB triplet
 *     FLAGS  = 2 - force to indexed color
 *     FLAGS |= 4 - to accept multiple colors (can be OR'ed with one of the
 *                  other values)
 *   hence:
 *     0 - (default) means: only scalar color accepted, result is either
 *         an index or a triplet;
 *     1 - means: only scalar color accepted, result is always a triplet;
 *     2 - means: only scalar color accepted, result is always an index;
 *     4 - means: multiple colors accepted, result is either an index or
 *         a triplet;
 *     5 - means: multiple colors accepted, result is always a triplet;
 *     6 - means: multiple colors accepted, result is always an index;
 *
 *   If COLOR is nil, the returned value is PL_FG (index for "fg" color).
 *
 *   The following predefined color names, vales and constants (as Yorick
 *   global variables) are supported:
 *
 *      NAME    INTEGER   CONSTANT       NAME      INTEGER    CONSTANT
 *     ---------------------------       -------------------------------
 *     "bg"     255   -1  PL_BG          "magenta"  247   -9  PL_MAGENTA
 *     "fg"     254   -2  PL_FG          "yellow"   246  -10  PL_YELLOW
 *     "black"  253   -3  PL_BLACK       "grayd"    245  -11  PL_GRAYD
 *     "white"  252   -4  PL_WHITE       "grayc"    244  -12  PL_GRAYC
 *     "red"    251   -5  PL_RED         "grayb"    243  -13  PL_GRAYB
 *     "green"  250   -6  PL_GREEN       "graya"    242  -14  PL_GRAYA
 *     "blue"   249   -7  PL_BLUE        "xor"      241  -15  PL_XOR
 *     "cyan"   248   -8  PL_CYAN        "extra"    240  -16  PL_EXTRA
 *     ---------------------------       -------------------------------
 *
 *   Color names are case insensitive and a most X11 color names are
 *   available (thanks to a builtin database).  For instance
 *   "darkslateblue" and "DarkSlateBlue" are valid color names.
 *
 *
 * SEE ALSO:
 *  color, pl_get_symbol, pl_get_font, pl_get_palette.
 */
{
  /* Parse flags (default is to get a single color in a numerical form
     which is usable as the value of the COLOR keyword in Yorick plotting
     functions). */
  if (is_void(flags)) {
    mode = 0; /* any */
    single = 1n;
  } else {
    mode = (flags & 3);
    single = ! (flags & 4);
  }

  /* Parse color value(s). */
  /**/extern _PL_COLOR_NAMES, _PL_COLOR_RGB;
  if (is_void(color)) {
    return (mode == 1 ? _PL_COLOR_RGB(,2) : PL_FG);
  }
  if ((s = structof(color)) == string) {
    dims = dimsof(color);
    if (single && dims(1)) error, "color name must be a scalar";
    ndx = array(long, dims);
    len = strlen(color);
    n = numberof(color);
    for (k = 1; k <= n; ++k) {
      /* Timing Issues: It is slightly faster to use strfind (no string
         duplication).  Searching the color database, takes ~ 10
         microseconds per color on my 2GHz Centrino. */
      sel = strfind(color(k), _PL_COLOR_NAMES, case=0);
      if (is_array((i = where(! sel(1,..)))) &&
          is_array((j = where(sel(2, i) == len(k))))) {
        ndx(k) = i(j(1));
      } else {
        error, ("unrecognized color name: \"" + color(k) + "\"");
      }
    }
    if (! mode) {
      return (max(ndx) <= 16 ? char(256 - ndx) : _PL_COLOR_RGB(, ndx));
    } else if (mode == 1) {
      return _PL_COLOR_RGB(, ndx);
    } else {
      if (max(ndx) > 16) {
        error, "named color is not an indexed one";
      }
      return char(256 - ndx);
    }
  }
  if (s == char) {
    /* Can be an indexed color or a triplet. */
    dims = dimsof(color);
    n = dims(1);
    is_rgb = (n >= 1 && dims(2) == 3);
    if (single && n != is_rgb) {
      error, "expecting a single color";
    }
    if (is_rgb) {
      if (mode == 2) {
        /* FIXME: cannot force triplet to index */
        error, "unsupported conversion of [R,G,B] to indexed color";
      }
    } else if (mode == 1) {
      return pl_get_palette()(, 1 + color);
    }
    return color;
  }
  if (s == long || s == short || s == int) {
    dims = dimsof(color);
    if (single && dims(1)) {
      error, "expecting a single color";
    }
    if ((cmin = min(color)) < 0) {
      if (cmin < -16) {
        error, "out of range indexed color";
      }
      color = long(color); /* force copy/conversion */
      color(where(color < 0)) += 256;
      cmin = min(color);
    }
    cmax = max(color);
    if (cmax <= 255) {
      /* Assume all colors are indexed ones. */
      if (mode == 1) {
        return pl_get_palette()(, 1 + color);
      }
      return char(color);
    }
    if (mode == 2) {
      /* FIXME: cannot force triplet to index */
      error, "unsupported conversion of [R,G,B] to indexed color";
    }
    rgb = array(char, 3, dims);
    if (cmin >= 256) {
      /* Assume all colors are packed RGB triplets. */
      rgb(1, ..) = (color & 0xff);
      rgb(2, ..) = ((color >> 8) & 0xff);
      rgb(3, ..) = ((color >> 16) & 0xff);
      return rgb;
    } else {
      /* There is a mixture of indexed and packed RGB colors. */
      lut = pl_get_palette();
      ndx = (color <= 255);
      if (is_array((i = where(nxd)))) {
        rgb(, i) = lut(, 1 + color(i));
      }
      if (is_array(( i = where(! nxd)))) {
        color = color(i);
        rgb(1, i) = (color & 0xff);
        rgb(2, i) = ((color >> 8) & 0xff);
        rgb(3, i) = ((color >> 16) & 0xff);
      }
    }
    return rgb;
  }
  error, "bad data type for color";
}

/* Color codes as defined in 'play.h': */
PL_BG      = char(255);
PL_FG      = char(254);
PL_BLACK   = char(253);
PL_WHITE   = char(252);
PL_RED     = char(251);
PL_GREEN   = char(250);
PL_BLUE    = char(249);
PL_CYAN    = char(248);
PL_MAGENTA = char(247);
PL_YELLOW  = char(246);
PL_GRAYD   = char(245);
PL_GRAYC   = char(244);
PL_GRAYB   = char(243);
PL_GRAYA   = char(242);
PL_XOR     = char(241);
PL_EXTRA   = char(240);

local PL_RGB_BG, PL_RGB_FG, PL_RGB_XOR, PL_RGB_EXTRA;
local PL_RGB_BLACK, PL_RGB_WHITE, PL_RGB_RED, PL_RGB_GREEN, PL_RGB_BLUE;
local PL_RGB_CYAN, PL_RGB_MAGENTA, PL_RGB_YELLOW;
local PL_RGB_GRAYD, PL_RGB_GRAYC, PL_RGB_GRAYB, PL_RGB_GRAYA;
func pl_get_palette(win)
/* DOCUMENT rgb = pl_get_palette();
 *     -or- rgb = pl_get_palette(win);
 *
 *   Returns color table for current graphics window or for winddow WIN if
 *   it is specified.  The color table is a 3-by-256 array of char's such
 *   that: RGB(1,i), RGB(2,i) and RGB(3,i) are the red, green and blue
 *   levels for color index i respectively.
 *
 *   Note that some color index cannot be obtained ("bf", "fg", "xor" and
 *   "extra") and are left as black/white in the returned color table.  It
 *   is possible to change the predefined colors by setting the global
 *   variables such as PL_RGB_BG with a given [r,g,b] triplet.
 *
 * SEE ALSO: palette, pl_get_color.
 */
{
  local r, g, b, old_win;
  if (! is_void(win)) {
    old_win = current_window();
    window, win;
  }
  palette, r, g, b, query=1;
  lut = array(char, 3, 256);
  if (is_array(r)) {
    lut(1, 1:numberof(r)) = r;
    lut(2, 1:numberof(g)) = g;
    lut(3, 1:numberof(b)) = b;
  } else {
    /* Default is a gray colormap with K=240 colors and there are M=256
     * levels.  The formula is:
     *
     *    y = (2*M*x + M)/(2*K)
     *
     * where integer arithmetic must be used and x = 0, ..., K - 1 is the
     * 0-based color index and y = 0, ..., M - 1 is the color level.
     * After simplification (with K=240 and M=256):
     */
    lut( , 1:240) = (indgen(8:3832:16)/15)(-, );
  }
  lut(, 1 + PL_BG)      = PL_RGB_BG;
  lut(, 1 + PL_FG)      = PL_RGB_FG;
  lut(, 1 + PL_BLACK)   = PL_RGB_BLACK;
  lut(, 1 + PL_WHITE)   = PL_RGB_WHITE;
  lut(, 1 + PL_RED)     = PL_RGB_RED;
  lut(, 1 + PL_GREEN)   = PL_RGB_GREEN;
  lut(, 1 + PL_BLUE)    = PL_RGB_BLUE;
  lut(, 1 + PL_CYAN)    = PL_RGB_CYAN;
  lut(, 1 + PL_MAGENTA) = PL_RGB_MAGENTA;
  lut(, 1 + PL_YELLOW)  = PL_RGB_YELLOW;
  lut(, 1 + PL_GRAYD)   = PL_RGB_GRAYD;
  lut(, 1 + PL_GRAYC)   = PL_RGB_GRAYC;
  lut(, 1 + PL_GRAYB)   = PL_RGB_GRAYB;
  lut(, 1 + PL_GRAYA)   = PL_RGB_GRAYA;
  lut(, 1 + PL_XOR)     = PL_RGB_XOR;
  lut(, 1 + PL_EXTRA)   = PL_RGB_EXTRA;
  if (win != old_win) {
    window, old_win;
  }
  return lut;
}

PL_RGB_BG      = char([214, 214, 214]);
PL_RGB_FG      = char([  0,   0,   0]);
PL_RGB_BLACK   = char([  0,   0,   0]);
PL_RGB_WHITE   = char([255, 255, 255]);
PL_RGB_RED     = char([255,   0,   0]);
PL_RGB_GREEN   = char([  0, 255,   0]);
PL_RGB_BLUE    = char([  0,   0, 255]);
PL_RGB_CYAN    = char([  0, 255, 255]);
PL_RGB_MAGENTA = char([255,   0, 255]);
PL_RGB_YELLOW  = char([255, 255,   0]);
PL_RGB_GRAYD   = char([100, 100, 100]);
PL_RGB_GRAYC   = char([150, 150, 150]);
PL_RGB_GRAYB   = char([190, 190, 190]);
PL_RGB_GRAYA   = char([214, 214, 214]);
PL_RGB_XOR     = char([  0,   0,   0]);
PL_RGB_EXTRA   = char([  0,   0,   0]);


local PL_COURIER, PL_TIMES, PL_HELVETICA, PL_SYMBOL, PL_NEWCENTURY;
local PL_GUI_FONT, PL_BOLD, PL_ITALIC, PL_OPAQUE;
local _PL_FONT_TABLE;
func pl_get_font(value, default)
/* DOCUMENT pl_get_font(value, default);
 *   Parse font VALUE which can have any value recognized by Yorick
 *   for the "font" keyword in builtin plotting functions and return the
 *   corresponding integer value.  In addition, if VALUE is void, DEFAULT
 *   is returned.
 *
 * SEE ALSO:
 *   font, pl_get_color, , pl_get_symbol, xwindow.
 */
{
  extern _PL_FONT_TABLE;
  if (is_void(value)) return default;
  if (is_array(value) && ! dimsof(value)(1)) {
    if ((s = structof(value)) == long) return value;
    if (s == string) {
      n = strmatch(value, "B") | 2*strmatch(value, "I");
      fn = (n==3 ? strpart(value, 1:-2) : (n ? strpart(value, 1:-1) : value));
      if (is_func(h_new)) {
        if (! is_hash(_PL_FONT_TABLE)) {
          _PL_FONT_TABLE = h_new("courier",1, "times",2, "helvetica",3,
                                 "symbol",4,"schoolbook",5);
        }
        index = _PL_FONT_TABLE(fn);
        if (index) return n + 4*(index - 1);
      } else {
        if (fn == "courier")    return n;
        if (fn == "times")      return n+4;
        if (fn == "helvetica")  return n+8;
        if (fn == "symbol")     return n+12;
        if (fn == "schoolbook") return n+16;
      }
      error, "bad font name \""+value+"\"";
    }
    if (s == char || s == short || s == int) return long(value);
  }
  error, "bad font value";
}

/* Font codes/flags as defined in 'play.h': */
PL_COURIER    =  0;
PL_TIMES      =  4;
PL_HELVETICA  =  8;
PL_SYMBOL     = 12;
PL_NEWCENTURY = 16;
PL_GUI_FONT   = 20;
PL_BOLD       =  1;
PL_ITALIC     =  2;
PL_OPAQUE     = 32;

local PL_SQUARE, PL_PLUS, PL_TRIANGLE, PL_UP_TRIANGLE, PL_CIRCLE;
local PL_DIAMOND, PL_CROSS, PL_DOWN_TRIANGLE, PL_STAR;
local _PL_SYMBOL_TABLE;
func pl_get_symbol(symbol)
/* DOCUMENT pl_get_symbol(symbol);
 *
 *   Get symbol value as an integer, SYMBOL must be a scalar and may be either
 *   an integer, a character or a string:
 *
 *     INT CHAR  STRING                 DESCRIPTION
 *     ----------------------------------------------------------------------
 *      0                               nothing (just draw error bars if any)
 *      1   #  "square"                 a square
 *      2   +  "plus"                   a plus sign
 *      3   ^  "triangle" "uptriangle"  a triangle
 *      4   o  "circle"                 a circle (actually an hexagon)
 *      5   @  "diamond"                a square rotated by 45 degrees
 *      6   x  "cross"                  an X-cross    <- this is the default
 *      7   v  "downtriangle"           an upside down triangle
 *      8   *  "star"                   a star
 *     >=9                              a polygon with SYMBOL sides
 *     ----------------------------------------------------------------------
 *
 *   The one-character symbol may given as lower/upper case and as a string
 *   or a  char; e.g. 'v', 'V',  "v" and "V"  all stand for an  upside down
 *   triangle.
 *
 *   For  convenience, global  variables  PL_SQUARE, PL_PLUS,  PL_TRIANGLE,
 *   PL_UP_TRIANGLE,     PL_CIRCLE;     local     PL_DIAMOND,     PL_CROSS,
 *   PL_DOWN_TRIANGLE and PL_STAR are defined with the corresponding symbol
 *   code.
 *
 * SEE ALSO: plp, pl_get_color.
 */
{
  if (is_void(symbol)) {
    return 6;
  }
  if (! is_array(symbol) || dimsof(symbol)(1)) {
    error, "symbol must be a scalar";
  }
  s = structof(symbol);
  if (s == string || s == char) {
    if (is_func(h_new)) {
      /* Use Yeti hash-table to speed-up symbol identification. */
      extern _PL_SYMBOL_TABLE;
      if (! is_hash(_PL_SYMBOL_TABLE)) {
        _PL_SYMBOL_TABLE = h_new("square",1, "#",1,
                                 "plus",2, "+",2,
                                 "triangle",3, "uptriangle",3, "^",3,
                                 "circle",4, "o",4, "O",4,
                                 "diamond",5, "@",5,
                                 "cross",6, "x",6, "X",6,
                                 "downtriangle",7, "v",7, "V",7,
                                 "star",8, "*",8);
      }
      if (s == char) {
        symbol = strchar(symbol);
      }
      symbol = _PL_SYMBOL_TABLE(symbol);
      if (symbol) {
        return symbol;
      }
    } else {
      /* Use vanilla Yorick. */
      local c;
      if (s == char) {
        len = 1;
        eq_nocopy, c, symbol;
      } else {
        len = strlen(symbol);
        c = strchar(symbol)(1);
      }
      if (len == 1) {
        if (c=='#') return 1;
        if (c=='+') return 2;
        if (c=='^') return 3;
        if (c=='o' || c =='O') return 4;
        if (c=='@') return 5;
        if (c=='x' || c=='X') return 6;
        if (c=='v' || c=='V') return 7;
        if (c=='*') return 8;
      } else {
        /* must be a string */
        if (c=='s') {
          if (symbol=="square") return 1;
          if (symbol == "star") return 8;
        } else if (c=='p') {
          if (symbol=="plus") return 2;
        } else if (c=='t') {
          if (symbol == "triangle") return 3;
        } else if (c=='c') {
          if (symbol == "circle") return 4;
          if (symbol == "cross") return 6;
        } else if (c=='d') {
          if (symbol == "diamond") return 5;
          if (symbol == "downtriangle") return 7;
        } else if (c=='u') {
          if (symbol=="uptriangle") return 3;
        }
      }
    }
  } else if ((s == long || s == int || s == short) && symbol >= 0) {
    return long(symbol);
  }
  error, "bad symbol value";
}

/* Symbol codes as used by 'plp': */
PL_SQUARE        = 1;
PL_PLUS          = 2;
PL_TRIANGLE      = 3;
PL_UP_TRIANGLE   = 3;
PL_CIRCLE        = 4;
PL_DIAMOND       = 5;
PL_CROSS         = 6;
PL_DOWN_TRIANGLE = 7;
PL_STAR          = 8;

func pl_get_axis_flags(value, default)
/* DOCUMENT pl_get_axis_flags(value, default);
 *   Parse axis  specification VALUE  which can be  an integer or  a string
 *   where each bits/character toggle an option (see table below). If VALUE
 *   is void, DEFAULT is returned.
 *
 *   char  bit    option
 *   ----  -----  ----------------------------------------------------
 *   t     0x001  Draw ticks on bottom or left edge of viewport
 *   T     0x002  Draw ticks on top or right edge of viewport
 *   c     0x004  Draw ticks centered on origin in middle of viewport
 *   i     0x008  Ticks project inward into viewport
 *   o     0x010  Ticks project outward away from viewport (0x18 for both)
 *   l     0x020  Draw tick label numbers on bottom or left edge of viewport
 *   L     0x040  Draw tick label numbers on top or right edge of viewport
 *   g     0x080  Draw all grid lines down to gridLevel
 *   z     0x100  Draw single grid line at origin
 *
 * SEE ALSO:
 *   xwindow.
 */
{
  if (is_void(value)) return default;
  if (is_scalar(value)) {
    if (is_integer(value)) return long(value);
    if (is_string(value)) {
      flags = 0;
      if (strmatch(value, "t")) flags |= 0x001;
      if (strmatch(value, "T")) flags |= 0x002;
      if (strmatch(value, "c")) flags |= 0x004;
      if (strmatch(value, "i")) flags |= 0x008;
      if (strmatch(value, "o")) flags |= 0x010;
      if (strmatch(value, "l")) flags |= 0x020;
      if (strmatch(value, "L")) flags |= 0x040;
      if (strmatch(value, "g")) flags |= 0x080;
      if (strmatch(value, "z")) flags |= 0x100;
      return flags;
    }
  }
  error, "bad axis flag value";
}

/* Line types as defined in 'play.h': */
PL_SOLID      = 0;
PL_DASH       = 1;
PL_DOT        = 2;
PL_DASHDOT    = 3;
PL_DASHDOTDOT = 4;
PL_SQUARE     = 8;


func pl_set_ndigits(nx, ny)
/* DOCUMENT pl_set_ndigits, nx, ny;
     Set the number of digits for the ticks of a plot. If specified, NX and NY
     are the number of digits for the horizontal and vertical axis respectively.
   SEE ALSO:
 */
{
  local landscape, systems, legends, clegends;
  flags = ((is_void(nx) ? 0 : 1) | (is_void(ny) ? 0 : 2));
  if (flags != 0) {
    get_style, landscape, systems, legends, clegends;
    if ((flags & 1) != 0) systems.ticks.horiz.nDigits = nx;
    if ((flags & 2) != 0) systems.ticks.vert.nDigits = ny;
    set_style, landscape, systems, legends, clegends;
  }
}


/*---------------------------------------------------------------------------*/

/*
 * 3D TRANSFORM:
 * -------------
 * Let (X,Y,Z)  be  the  data  coordinates  (in  a  direct  frame) and
 * (XP,YP,ZP) be the coordinates in the view frame (also direct, XP is
 * from left to right, YP is from bottom to top and  ZP  points toward
 * the observer), then  for  an  altitude  ALT  (angle  of  view above
 * XY-plane) and an azimuth  AZ  (angle  of  view  around  Z-axis) the
 * coordinates transform is obtained  by  a  rotation  around  Oz with
 * angle AZ (azimuth), followed by a rotation around Ox with  angle AX
 * (AX = ALT - 90deg):
 *   XP = X cos(AZ) - Y sin(AZ)
 *   YP = U cos(AX) - Z sin(AX) =   U sin(ALT) + Z cos(ALT)
 *   ZP = U sin(AX) + Z cos(AX) = - U cos(ALT) + Z sin(ALT)
 * where:
 *   U  = X sin(AZ) + Y cos(AZ)
 */

func _pl3xyz(&xp, &yp, &zp, x, y, z)
{
/* DOCUMENT _pl3xyz, xp, yp, zp, x, y, z
     transform data coordinates (X,Y,Z) into viewer coordinates
     (XP,YP,ZP) for an externally defined altitude ALT and
     azimuth AZ (in degrees).
 */
  extern alt, az;
  if (is_void(az))  az = 30.0;	// angle of view around z-axis
  if (is_void(alt)) alt = 45.0;	// angle of view above xy-plane
  d2r = pi / 180.0;
  xp = x * (c = cos(az * d2r)) - y * (s = sin(az * d2r));
  zp = x * s + y * c;
  yp = z * (c = cos(alt * d2r)) + zp * (s = sin(alt * d2r));
  zp = z * s - zp * c;
}

func _pl3xy(&xp, &yp, x, y, z)
{
/* DOCUMENT _pl3xy, xp, yp, x, y, z
     transform data coordinates (X,Y,Z) into viewer coordinates
     (XP,YP) for an externally defined altitude ALT and
     azimuth AZ (in degrees).
 */
  extern alt, az;
  if (is_void(az))  az = 30.0;	// angle of view around z-axis
  if (is_void(alt)) alt = 45.0;	// angle of view above xy-plane
  d2r = pi / 180.0;
  xp = x * (c = cos(az * d2r)) - y * (s = sin(az * d2r));
  yp = z * cos(alt * d2r) + (x * s + y * c) * sin(alt * d2r);
}

/*---------------------------------------------------------------------------*/

func pl3t(text, x, y, z, alt=, az=, legend=, hide=, color=, font=, height=,
	  opaque=, path=, justify=, tosys=)
{
/* DOCUMENT pl3t, text, x, y, z, alt=alt, az=az, tosys=0/1
     plots TEXT (a string) at the point (X,Y,Z) in a 3-dimensional
     graph view from altitude ALT (default 45) and azimuth AZ
     (default 30) both in degrees.   TEXT, X, Y and Z may be arrays
     with the same number of elements.

     Other optional keywords are:
       legend, hide, color, font, height, opaque, path, justify and tosys
     and have the same meaning as in plt.

   SEE ALSO: plt.
 */
  local xp, yp;
  _pl3xy, xp, yp, x, y, z;
  n = numberof(text);
  if (n == 1) {
    plt, text, xp, yp,
      legend=legend, hide=hide, color=color,font=font, height=height,
      opaque=opaque, path=path, justify=justify, tosys=tosys;
  } else {
    for (i=1; i<=n; i++) {
      plt, text(i), xp(i), yp(i),
	legend=legend, hide=hide, color=color,font=font, height=height,
	opaque=opaque, path=path, justify=justify, tosys=tosys;
    }
  }
}

/*---------------------------------------------------------------------------*/

func pl3dj(x0, y0, z0, x1, y1, z1, alt=, az=,
	   legend=, hide=, type=, width=, color=)
{
/* DOCUMENT pl3dj, x0, y0, z0, x1, y1, z1, alt=alt, az=az
     plots disjoint lines from (X0,Y0,Z0) to (X1,Y1,Z1) in a 3-dimensional
     graph view from altitude ALT (default 45) and azimuth AZ (default 30)
     both in degrees.  X0, Y0, Z0, X1, Y1 and Z1 must have the same shapes.

     Additional keywords are those accepted by pldj: legend, hide, type,
     width, and color.

   SEE ALSO: pldj, pl3t, pl3s.
 */
  local x0p, y0p, x1p, y1p;
  _pl3xy, x0p, y0p, x0, y0, z0;
  _pl3xy, x1p, y1p, x1, y1, z1;
  pldj, x0p, y0p, x1p, y1p,
    legend=legend, hide=hide, type=type, width=width, color=color;
}

/*---------------------------------------------------------------------------*/

func pl3s(z, y, x, alt=, az=, axis=, box=, acolor=, fill=, legend=, hide=,
          edges=, ecolor=, ewidth=, height=, font=)
/* DOCUMENT pl3s, z, y, x, fill=0/1/2
         or pl3s, z, fill=0/1/2
      draws 3-D surface plot of Z versus (X,Y).   The  Z  array  must  be a
      2-dimensional array, say NX-by-NY and X  and  Y  must  have  the same
      shape as Z or be vectors  of  length  NX  and  NY  respectively.   If
      omitted, X and Y are set to the first and second  indice  value  of Z
      respectively.

      The FILL keyword indicates the kind of plot: 0 (default) for  3D wire
      frames, 1 for 3D mesh filled with intensity,  2  for  3D  mesh shaded
      with light source aligned with observer.

      The altitude and azimuth angles (in degrees) can be set with keywords
      ALT and AZ, their default values are 30 and 45 deg.

      A solid edge can optionally be drawn around each zone by  setting the
      EDGES keyword non-zero.  ECOLOR and EWIDTH determine  the  edge color
      and width.

      Frame axis can optionally be drawn around  the  plot  by  setting the
      AXIS keyword non-zero.  The  color  of  the  axis  and  label  can be
      modified with keyword ACOLOR.

      If BOX keyword non-zero, the 3-D box borders are drawn (with the same
      color as the axis, i.e., ACOLOR).

   EXAMPLE
     It is usually better to select an eventually new window and choose
     the "nobox" style:
       window, max(0, current_window()), wait=1, style="nobox.gs";
       x = span(-3,3,50);
       y = span(-2,2,40);
       z = cos(x(,-) * y(-,));
       pl3s, z, y, x, axis=1, fill=2, edges=1, font="timesBI", height=10;

      The following keywords are legal (each has a separate help entry):
   KEYWORDS: legend, hide, region, edges, ecolor, ewidth, font, height.
   SEE ALSO: pl3dj, pl3t, plg, plm, plc, plv, plf, pli, plt, pldj, plfp,
     plmesh, limits, range, fma, hcp, palette, bytscl, ...
*/
{
  /*
   * Check dimensions of input arrays.
   */
  local tmp;
  if ((msg = pls_mesh(x, tmp, dimsof(z), which=1, inhibit=2)) != string(0) ||
      (msg = pls_mesh(y, tmp, dimsof(z), which=2, inhibit=2)) != string(0))
    error, msg;
  tmp = [];

  /*
   * Rescale arrays.
   */
  xspan = (xmax = max(x)) - (xmin = min(x));
  yspan = (ymax = max(y)) - (ymin = min(y));
  zspan = (zmax = max(z)) - (zmin = min(z));
  if (xspan <= 0.0 || yspan <= 0.0)
    error, "X and/or Y are constant";
  if (zspan <= 0.0) {
    zmin = -1.0;
    zmax = +1.0;
    zspan = 2.0;
    z(*) = 0.0;
  } else {
    z = (z - zmin) / zspan;
  }
  x = (x - xmin) / xspan;
  y = (y - ymin) / yspan;

  /*
   * Insure that angles are in the range [-180,180].
   */
  if (is_void(alt))                alt  =  45.0;
  else if ((alt %= 360.0) > 180.0) alt -= 360.0;
  else if (alt < -180.0)           alt += 360.0;
  if (is_void(az))                 az   =  30.0;
  else if ((az %= 360.0) > 180.0)  az  -= 360.0;
  else if (az < -180.0)            az  += 360.0;

  /*
   * Plot an invisible box around the plot to left some space
   * around for the axis labels.
   */
  if (axis) {
    local bxp, byp;
    if (is_void(height))
      height = 10.0;
    q = height / 75.0;
    _pl3xy, bxp, byp,
      [0, 1, 1, 0, 0, 1, 1, 0] + q * [-1, 1, 1,-1,-1, 1, 1,-1],
      [0, 0, 1, 1, 0, 0, 1, 1] + q * [-1,-1, 1, 1,-1,-1, 1, 1],
      [0, 0, 0, 0, 1, 1, 1, 1] + q * [-1,-1,-1,-1, 1, 1, 1, 1];
    bxmin = min(bxp); bxmax = max(bxp);
    bymin = min(byp); bymax = max(byp);
    pldj, bxmin, bymin, bxmax, bymax, type="none";
    pldj, bxmin, bymax, bxmax, bymin, type="none";
  }

  /*
   * Plot the rear of the 3-D box: must figure out which box faces
   * are seen, and then which box edges are seen.  The plotting order
   * is (1) rear part of the box, (2) surface mesh and (3) front
   * part of the box and axis.
   */
  if (box) {
    local bxp0, byp0, bxp1, byp1;	// end-points of box edges
    local nxp, nyp, nzp;		// vectors normal to box faces
    _pl3xy, bxp0, byp0,
      [0,1,1,0,0,1,1,0,0,1,1,0],
      [0,0,1,1,0,0,1,1,0,0,1,1],
      [0,0,0,0,1,1,1,1,1,1,1,1];
    _pl3xy, bxp1, byp1,
      [1,1,0,0,1,1,0,0,0,1,1,0],
      [0,1,1,0,0,1,1,0,0,0,1,1],
      [0,0,0,0,1,1,1,1,0,0,0,0];
    _pl3xyz, nxp, nyp, nzp,
      [ 0, 0, 0, 1, 0,-1],
      [ 0, 0,-1, 0, 1, 0],
      [-1, 1, 0, 0, 0, 0];
    face_edges=[[1,2,3,4], [5,6,7,8], [1,5,9,10], [2,6,10,11],
		[3,7,11,12], [4,8,9,12]];
    visible = array(0, 12);
    visible(face_edges(, where(nzp >= 0.0))(*)) = 1;
    fore = where(visible);
    back = where(!visible);
    pldj, bxp0(back), byp0(back), bxp1(back), byp1(back), type=1, color=acolor;
  }

  /*
   * Rotate the surface so as to have the drawing starting at back
   * end (i.e., hidden surfaces are drawn first).  To this end, the
   * first thing to do is to localize the corner of surface z=0
   * which is the farest from the observer.  Depending on the
   * position of the farest corner, X, Y and Z arrays
   * may have to be scrambled.  After what, the 3-D projection
   * of the surface can be computed.
   */
  local xp, yp, zp;
  _pl3xyz, xp, yp, zp,
    [x(0,1), x(1,0), x(0,0), x(1,1)],
    [y(0,1), y(1,0), y(0,0), y(1,1)],
    0;
  far = zp(mnx)(1);		// index of farest point
  if (far == 1) {
    x = x(::-1,);
    y = y(::-1,);
    z = z(::-1,);
  } else if (far == 2) {
    x = x(,::-1);
    y = y(,::-1);
    z = z(,::-1);
  } else if (far == 3) {
    x = x(::-1,::-1);
    y = y(::-1,::-1);
    z = z(::-1,::-1);
  }
  _pl3xyz, xp, yp, zp, x, y, z;

  /*
   * Plot surface.
   */
  if (!fill) {
    colors = [];
    edges = 1;
  } else if (fill == 1) {
    colors = bytscl(z, cmin=0.0, cmax=1.0);
  } else {
    /* compute the two median vectors for each cell */
    m0x = xp(dif,zcen);
    m0y = yp(dif,zcen);
    m0z = zp(dif,zcen);
    m1x = xp(zcen,dif);
    m1y = yp(zcen,dif);
    m1z = zp(zcen,dif);
    /* define the normal vector to be their cross product */
    nx = m0y*m1z - m0z*m1y;
    ny = m0z*m1x - m0x*m1z;
    nz = m0y*m1x - m0x*m1y;
    m0x = m0y = m0z = m1x = m1y = m1z = [];
    colors = bytscl(nz);
    //colors = bytscl(nz / abs(nx, ny, nz), cmin=0.0, cmax=1.0);
    nx = ny = nz = [];
  }
  plf, colors, yp, xp, legend=legend, hide=hide,
    edges=edges, ecolor=ecolor, ewidth=ewidth;
  xp = yp = zp = colors = [];

  /*
   * Plot the axis and the front of the box.
   */
  if (axis) {
    if (far == 1) {
      px = [ 0, 1, 0]; vx = [ 0, 1, 0];
      py = [ 0, 0, 0]; vy = [-1, 0, 0];
      pz = [ 1, 1, 0]; vz = [ 1, 0, 0];
    } else if (far == 2) {
      px = [ 0, 0, 0]; vx = [ 0,-1, 0];
      py = [ 1, 0, 0]; vy = [ 1, 0, 0];
      pz = [ 0, 0, 0]; vz = [-1, 0, 0];
    } else if (far == 3) {
      px = [ 0, 0, 0]; vx = [ 0,-1, 0];
      py = [ 0, 0, 0]; vy = [-1, 0, 0];
      pz = [ 0, 1, 0]; vz = [ 0, 1, 0];
    } else {
      px = [ 0, 1, 0]; vx = [ 0, 1, 0];
      py = [ 1, 0, 0]; vy = [ 1, 0, 0];
      pz = [ 1, 0, 0]; vz = [ 0,-1, 0];
    }
    _pl3tick, xmin, xmax, px, [1,0,0], 0.02 * vx,
      height=height, font=font, color=acolor;
    _pl3tick, ymin, ymax, py, [0,1,0], 0.02 * vy,
      height=height, font=font, color=acolor;
    _pl3tick, zmin, zmax, pz, [0,0,1], 0.02 * vz,
      height=height, font=font, color=acolor;
  }
  if (box) {
    pldj, bxp0(fore), byp0(fore), bxp1(fore), byp1(fore), type=1, color=acolor;
  }
}

/*---------------------------------------------------------------------------*/

func _pl3tick(tmin, tmax, p, d, v, color=, font=, height=, opaque=, path=)
/* DOCUMENT _pl3tick, p, d, v
     draw axis between (P(1),P(2),P(3)) and ((P+D)(1),(P+D)(2),(P+D)(3))
     and ticks with vectorial length (V(1),V(2),V(3)).
 */
{
  extern alt, az;
  local px, py, dx, dy, vx, vy;

  /*
   * Compute projections, draw axis and figure out which
   * justification is the best one..
   */
  _pl3xy, px, py, p(1), p(2), p(3);
  _pl3xy, dx, dy, d(1), d(2), d(3);
  _pl3xy, vx, vy, v(1), v(2), v(3);
  pldj, px, py, px+dx, py+dy, color=color;
  if (where(d != 0)(1) == 3) {
    // horizontal ticks for z-axis
    vx = vx <= 0.0 ? -max(abs(v)) : max(abs(v));
    vy = 0.0;
  }
  justify = vx > 0.0 ? "L" : (vx ? "R" : "C");
  justify += vy > 0.0 ? "B" : (vy ? "T" : "H");

  /*
   * Compute step size in axis direction to get approximatively 6 ticks
   * and plot ticks.
   */
  tspan = tmax - tmin;
  tstep = 10.0^floor(log10(0.2 * tspan) + 0.5);
  if (tspan / tstep < 4)
    tstep *= 0.5;
  else if (tspan / tstep > 8)
    tstep *= 2.0;
  imin = ceil(tmin / tstep);
  imax = floor(tmax / tstep);
  ni = long(imax - imin + 1.0);
  t = tstep * span(imin, imax, ni);
  tn = (t - tmin) / tspan;
  x = px + dx * tn;
  y = py + dy * tn;
  pldj, x, y, x+vx, y+vy, legend=string(0), color=color;

  /*
   * Write tick labels.
   */
  x += 2 * vx;
  y += 2 * vy;
  t = swrite(format="%.3g", t);
  for (i = 0; i <= ni; i++) {
    plt, t(i), x(i), y(i), legend=string(0), justify=justify, tosys=1,
      color=color, font=font, height=height, opaque=opaque;
  }
}

/*---------------------------------------------------------------------------*/
/* EXTENDED WINDOW/VIEWPORT/STYLE INTERFACE */

func xwindow(win, width=, height=, dpi=, display=, private=, dump=, hcp=,
             parent=, xpos=, ypos=,
	     landscape=, viewport=, units=, font=, size=, color=, frame=,
	     keep=, xopt=, yopt=, xmargin=, ymargin=, debug=)
/* DOCUMENT xwindow, win;
       -or- xwindow;
       -or- xwindow(...);
     Create/switch to graphic window WIN.  This routine allows one to
     create a window of given size with viewport(s) different from the
     default one.  Otherwise its behaviour should mimics that of the
     "window" builtin routine.  If called as a function, the return
     value is an array of 6 double values:
       [WIN, ONE_PIXEL, XMIN, XMAX, YMIN, YMAX]
     where WIN is the window number, ONE_PIXEL is the pixel size in NDC
     units and [XMIN,XMAX,YMIN,YMAX] is the bounding box of the window
     in NDC units.

   KEYWORDS
     DPI, DISPLAY, PRIVATE, DUMP, HCP, PARENT, XPOS, YPOS - Same options
                 as for the window builtin routine.
     WIDTH, HEIGHT - Window width/height in pixels; default: 450 at 75dpi
                 and 600 at 100dpi.
     LANDSCAPE - If true, use lanscape orientation; else portrait.
     VIEWPORT  - Viewport coordinates: [XMIN, XMAX, YMIN, YMAX].  Several
                 viewports can be specified at the same time, in this
                 case, first dimension of VIEWPORT must be a 4 and
                 XMIN is VIEWPORT(1,..), XMAX is VIEWPORT(2,..) and so on.
     UNITS     - Units for the viewport keyword: 0 for NDC (default),
                 1 for relative, and 2 for pixels.
     FONT      - Font to use to label axes (see pl_get_font); default is
                 Helvetica.
     SIZE      - Text size for axis labels in points; default is 12.
     COLOR     - Axis color (see pl_get_color); default is foreground.
     KEEP      - Keep WIN if it already exists? Otherwise the window is
                 forced to be recreated (this is the default behaviour).
     XOPT      - Options for X-axis of viewport (see pl_get_axis_flags).
     YOPT      - Options for Y-axis of viewport (see pl_get_axis_flags).
     FRAME     - Line width for frame. Can be zero to disable drawing a
                 frame. Default is 1.0, the normal width of a line.
     XMARGIN, YMARGIN -

     Note: If you specify multiple viewports, FONT, SIZE, COLOR, XOPT and
           YOPT can be arrays to specify different options for each
           viewport (except that multiple colors cannot be specifed as RGB
           triplets).

   SEE ALSO pl_get_axis_flags, pl_get_color, pl_get_font, window. */
{
  if (! is_struct(GfakeSystem)) {
    include, "style.i", 1;
  }

  /* DPI can be set by pldefault but we have no way to know this value
     -- see graph.c */
  if (is_void(dpi)) dpi = 100;
  else              dpi = (dpi < 25 ? 25 : (dpi > 300 ? 300 : long(dpi)));

  /* Define some constants -- see gist.h */
  ONE_POINT = 0.0013000;       // one point in NDC units
  ONE_INCH  = 72.27*ONE_POINT; // one inch in NDC units
  ONE_PIXEL = ONE_INCH/dpi;    // one pixel in NDC units

  /* Figure out window size in pixels (without top text line) -- called
     "topWidth" and "topHeight" in xfancy.c */
  if (is_void(width)) width = 6*dpi;
  else width = (width <= 31 ? 31 : long(width));
  if (is_void(height)) height = 6*dpi;
  else height = (height <= 31 ? 31 : long(height));

  /* Page size (8.5 x 11 inches) in NDC and pixels -- see xfancy.c */
  if (is_void(landscape)) landscape = (width > height);
  if (landscape) {
    pageWidthNDC = 1.033461;
    pageHeightNDC = 0.798584;
  } else {
    pageWidthNDC = 0.798584;
    pageHeightNDC = 1.033461;
  }
  pageWidth = long(pageWidthNDC/ONE_PIXEL);
  pageHeight = long(pageHeightNDC/ONE_PIXEL);
  if (width > pageWidth) width = pageWidth;
  if (height > pageHeight) height = pageHeight;

  /* Compute offsets (really tricky to figure out!) and bounding viewport
     -- see xfancy.c */
  xoff = (pageWidth - width)/2;
  if (landscape) {
    //yoff = (pageHeight - height)/2;
    yoff = (pageHeight - height + 3)/2;
  } else {
    //yoff = (pageHeight - height) - (pageWidth - height)/2;
    yoff = pageHeight - (pageWidth + height)/2;
  }
  if (xoff < 0) xoff = 0;
  if (yoff < 0) yoff = 0;
  vp0_offsets = [ONE_PIXEL*xoff, ONE_PIXEL*xoff,
                 ONE_PIXEL*yoff, ONE_PIXEL*yoff];
  vp0_pixel = [width, width, height, height] - 1.0;
  vp0_ndc = vp0_offsets + ONE_PIXEL*vp0_pixel;

  /* Fix viewport coordinates and figure out how many viewports to make. */
  if (is_void(viewport)) {
    viewport = [0.15, 0.85, 0.1, 0.8];
    units = 1;
  } else if (is_void(units)) {
    units = 0;
  } else if (is_array(units) && ! dimsof(units)(1)) {
    if ((s = structof(units)) == string) {
      if (units == "ndc") units = 0;
      else if (units == "pixel") units = 2;
      else if (units == "relative") units = 1;
      else units = -1;
    } else if (s != long && s != char && s != short && s != int) {
      units = -1;
    }
  } else {
    units = -1;
  }
  if (units) {
    /* Compute viewport size in pixels, then in NDC. */
    if (units == 1) {
      /* relative units: 0 is min and 1 is max */
      viewport *= vp0_pixel;
    } else if (units != 2) {
      error, "bad value for keyword UNITS (0 or \"ndc\", 1 or \"relative\", 2 or \"pixel\")";
    }
    viewport = vp0_offsets + ONE_PIXEL*viewport;
  }
  if (! is_array(viewport) || (dims = dimsof(viewport))(1) < 1 || dims(2) != 4)
    error, "bad VIEWPORT array";
  nvps = numberof(viewport)/4;

  /* Parse other parameters. */
  if (is_void(xmargin)) xmargin = 2.0;
  if (is_void(ymargin)) ymargin = 2.0;
  if (is_void(frame)) frame = 1.0;
  if (is_void(size)) size = 12;
  font = pl_map(pl_get_font, font, 8);
  color = pl_map(pl_get_color, color, 254 /* fg */);
  xopt = pl_map(pl_get_axis_flags, xopt, 0x33);
  yopt = pl_map(pl_get_axis_flags, yopt, 0x33);
  if (nvps>1) {
    /* Multiple viewports: fix per-viewport settings. */
    dims = dims(2:);
    dims(1) = numberof(dims)-1;
    zero = array(long, dims);
    if (is_void(dimsof(xopt, zero))) error, "bad XOPT dimensions";
    xopt += zero;
    if (is_void(dimsof(yopt, zero))) error, "bad YOPT dimensions";
    yopt += zero;
    if (is_void(dimsof(size, zero))) error, "bad SIZE dimensions";
    size += zero;
    if (is_void(dimsof(font, zero))) error, "bad FONT dimensions";
    font += zero;
    if (is_void(dimsof(color, zero))) error, "bad COLOR dimensions";
    color += zero;
    if (is_void(dimsof(frame, zero))) error, "bad FRAME dimensions";
    frame += zero;
  }
  if (numberof(xopt)  != nvps) error, "bad XOPT dimensions";
  if (numberof(yopt)  != nvps) error, "bad YOPT dimensions";
  if (numberof(size)  != nvps) error, "bad SIZE dimensions";
  if (numberof(font)  != nvps) error, "bad FONT dimensions";
  if (numberof(color) != nvps) error, "bad COLOR dimensions";
  if (numberof(frame) != nvps) error, "bad FRAME dimensions";

  /* Build up systems. */
  nHDigits = 5;
  nVDigits = 6;
  tickLen = ONE_PIXEL*[5.0, 3.0, 1.0, 0.0, 0.0];
  system = array(GfakeSystem, nvps);
  viewport = viewport(,*);  // concatenate extra dims
  for (i = 1; i <= nvps; ++i) {
    fontHeight = size(i)*ONE_POINT;
    xFlags = xopt(i);
    yFlags = yopt(i);
    currentColor = color(i);
    tickStyle =  GpLineAttribs(color=currentColor, type=1, width=1);
    frameLineWidth = frame(i);
    if (frameLineWidth) {
      frameFlag = 1n;
      frameStyle = GpLineAttribs(color=currentColor, type=1, width=frameLineWidth);
    } else {
      frameFlag = 0n;
      frameStyle = GpLineAttribs(color=currentColor, type=1, width=0.0);
    }
    gridStyle =  GpLineAttribs(color=currentColor, type=3, width=1);
    textStyle =  GpTextAttribs(color=currentColor, font=font(i), height=fontHeight,
                               alignH=0, alignV=0, opaque=0);
    xLabelOff = (xFlags&0x08 ? max(tickLen) : 0) + fontHeight/1.5;
    yLabelOff = (yFlags&0x08 ? max(tickLen) : 0) + 0.75*fontHeight;
    xTickOff = (xmargin
               ? xmargin*ONE_PIXEL + ((xFlags&0x08) ? max(tickLen) : 0.0)
               : 0.0);
    yTickOff = (ymargin
               ? ymargin*ONE_PIXEL + ((yFlags&0x08) ? max(tickLen) : 0.0)
               : 0.0);
    xTickLen = ((xFlags&0x018)==0x018 ? 2.0*tickLen : tickLen);
    yTickLen = ((yFlags&0x018)==0x018 ? 2.0*tickLen : tickLen);
    horiz = GaAxisStyle(nMajor=7.5, nMinor=50, logAdjMajor=1.2,
                        logAdjMinor=1.2, nDigits=nHDigits, gridLevel=1,
                        flags=xFlags, tickOff=xTickOff, labelOff=xLabelOff,
                        tickLen=xTickLen, tickStyle=tickStyle,
                        gridStyle=gridStyle, textStyle=textStyle,
                        xOver=viewport(2,i), yOver=viewport(3,i)-fontHeight);
    vert =  GaAxisStyle(nMajor=7.5, nMinor=50, logAdjMajor=1.2,
                        logAdjMinor=1.2, nDigits=nVDigits, gridLevel=1,
                        flags=yFlags,  tickOff=yTickOff, labelOff=yLabelOff,
                        tickLen=yTickLen, tickStyle=tickStyle,
                        gridStyle=gridStyle, textStyle=textStyle,
                        xOver=viewport(1,i), yOver=viewport(4,i)+fontHeight);
    system(i) = GfakeSystem(viewport=viewport(,i),
                           ticks=GaTickStyle(horiz=horiz, vert=vert, frame=frameFlag,
                                             frameStyle=frameStyle),
                           legend=swrite(format="System %d", i));
  }

  /* layout of the plot legends */
  legends = GeLegendBox();

  /* layout of the contour legends */
  clegends = GeLegendBox();

  /* create window and apply plot-style settings */
  if (keep) {
    dpi = width = height = display = private = parent = xpos = ypos = dump = [];
  } else {
    wn = (is_void(win) ? current_window() : win);
    if (wn >= 0) {
      window, wn, display="", hcp="";
    } else {
      wn = [];
    }
  }
  window, wn, wait=1, width=width, height=height,
    display=display, private=private, dump=dump, hcp=hcp, dpi=dpi,
    parent=parent, xpos=xpos, ypos=ypos;
  set_style, landscape, system, legends, clegends;
  fma;
  if (debug) {
    for (i = 1; i <= nvps; ++i) {
      pl_box, viewport(1,i), viewport(3,i), viewport(2,i),
        viewport(4,i), color="red", system=0;
    }
  }
  if (! am_subroutine()) {
    result = array(double, 6);
    result(1) = current_window();
    result(2) = ONE_PIXEL;
    result(3:6) = vp0_ndc;
    return result;
  }
}

/*---------------------------------------------------------------------------*/
/* SIMPLE BUTTON WIDGETS (based on button.i) */

struct XButton {
  double x, y;      /* NDC coordinates of button center */
  double dx, dy;    /* button half widths in NDC */
  string text;      /* button text */
  string font;      /* text font (0 for helvetica) */
  string color;     /* text and border color (0 for fg) */
  double height;    /* text height */
  double width;     /* width of line around button (0 is 1.0, <0 no box) */
}

func xbtn_plot(..)
/* DOCUMENT xbtn_plot, button1, button2, ...
     plot the specified BUTTONs.  Each button in the list may be an array
     of XButton structure.  Void arguments are no-ops.

   SEE ALSO: XButton, xbtn_which. */
{
  while (more_args()) {
    btns = next_arg();
    nbtns = numberof(btns);
    for (i=1 ; i<=nbtns ; ++i) {
      btn = btns(i);
      if (!(font = btn.font)) font = "helvetica";
      if (!(color = btn.color)) color = "fg";
      if ((h = btn.height)<=0.0) h = 14.0;
      x = btn.x;
      y = btn.y;
      plt, btn.text, x, y, justify="CH", font=font, height=h,
	color=color, opaque=1;
      dx = btn.dx;
      dy = btn.dy;
      if (!(w = btn.width)) w = 1.0;
      if (w>0.0) {
	oldSys = plsys(0);
	plg, [y-dy,y-dy,y+dy,y+dy],[x-dx,x+dx,x+dx,x-dx],
          closed=1, width=w, type=1, marks=0, color=color;
	plsys, oldSys;
      }
    }
  }
}

func xbtn_which(button, x, y)
/* DOCUMENT xbtn_which(button, x, y)
     Return index of element in button array BUTTON that contains NDC
     coordinates (X,Y).  Return -1 if coordinates do not match any
     buttons and -N if coordinates match N>1 buttons.

   SEE ALSO: XButton, xbtn_plot. */
{
  i = where((abs(y-button.y) < button.dy) & (abs(x-button.x) < button.dx));
  if ((n = numberof(i)) == 1) return i(1);
  if (n == 0) return -1;
  return -n;
}

/*---------------------------------------------------------------------------*/
/* EXTENDED MOUSE ROUTINES */

func xmouse_point(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_point()
       -or- xmouse_point
     Interactively  define  a  point  as with  xmouse("point").   The  same
     keywords as xmouse (which see) can be specified.

   SEE  ALSO: xmouse. */
{
  m = __xmouse(0, "click mouse to choose a position");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="position is (%g,%g)\n", m(1), m(2);
  } else if (is_array(m)) return m(1:2);
}

func xmouse_box(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_box()
       -or- xmouse_box
     Interactively  define a  rectangular box  as with  xmouse("box").  The
     same keywords as xmouse (which see) can be specified.

   SEE  ALSO: xmouse. */
{
  m = __xmouse(1, "click and drag mouse to select a region");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="%g x %g region is [%g : %g]x[%g : %g]\n",
           abs(m(3) - m(1)), abs(m(4) - m(2)),
           min(m(1), m(3)), max(m(1), m(3)),
           min(m(2), m(4)), max(m(2), m(4));
  } else if (is_array(m)) return [min(m(1), m(3)), max(m(1), m(3)),
                                  min(m(2), m(4)), max(m(2), m(4))];
}

func xmouse_line(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_line()
       -or- xmouse_line
     Interactively define a line as with xmouse("line").  The same keywords
     as xmouse (which see) can be specified.

   SEE  ALSO: xmouse. */
{
  m = __xmouse(2, "click and drag mouse to choose a line");
  if (am_subroutine()) {
    if (is_void(m)) {
      write, "aborted";
    } else {
      x0 = m(1);
      y0 = m(2);
      x1 = m(3);
      y1 = m(4);
      dx = x1 - x0;
      dy = y1 - y0;
      write,
        format="line end-points are (%g,%g) and (%g,%g), length = %g, angle=%g deg\n",
        x0, y0, x1, y1, abs(dx, dy), atan(dy, dx)*(180/pi);
    }
  } else if (is_array(m)) {
    return m(1:4);
  }
}

func xmouse_length(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_length()
       -or- xmouse_length
     Interactively  measure a  length as  with xmouse("length").   The same
     keywords as xmouse (which see) can be specified.

   SEE  ALSO: xmouse. */
{
  m = __xmouse(2, "click and drag mouse to measure a distance");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="distance from (%g,%g) to (%g,%g) = %g\n",
           m(1), m(2), m(3), m(4), abs(m(1) - m(3), m(2) - m(4));
  } else if (is_array(m)) return abs(m(1) - m(3), m(2) - m(4));
}

func xmouse(type,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse();
 *     -or- xmouse(type);
 *     -or- xmouse;
 *     -or- xmouse, type;
 *   Asks the  user to click into  a graphic window to  indicate a position
 *   (the default  if TYPE is not specified),  or to define a  segment or a
 *   rectangular box,  or to measure  a distance.  The possible  values for
 *   TYPE are:
 *
 *       TYPE          RESULT
 *     ---------   ----------------
 *     0 "point"   [X,Y]
 *     1 "box"     [XMIN, XMAX, YMIN, YMAX]
 *     2 "line"    [X0, Y0, X1, Y1]
 *     3 "length"  LENGTH
 *
 *   When called as a function, the coordinates of the result get returned;
 *   when called as a subroutine, the position is printed out.  If the user
 *   cancel the operation (see mouse function) or click into another window
 *   than the target one, no position is selected (nil get returned).  This
 *   behaviour can be  changed by setting keyword FOREVER  to true in which
 *   case the function loop until a valid selection is done.
 *
 *   Keyword WIN can be used to  specify the target graphic window which is
 *   by default the curent graphic window.
 *
 *   Keyword PROMPT can  be used to specify a  prompt string different from
 *   the default one.
 *
 *   Keyword SYSTEM  can be used to  specify the coordinate  system to use,
 *   the  default being  to use  the coordinate  system that  is  under the
 *   mouse.   The   returned  coordinates/length  are  in   units  of  that
 *   coordinate system.
 *
 *
 * SEE ALSO:
 *   mouse, xmouse_demo, xmouse_point, xmouse_box, xmouse_line,
 *   xmouse_length.
 */
{
  if (is_void(type)) {
    op = xmouse_point;
  } else if (is_array(type) && ! dimsof(type)(1)) {
    if ((s = structof(type)) == string) {
      op = symbol_def("xmouse_"+type);
      if (is_func(op) != 1) error, "unrecognized type name: \""+type+"\"";
    } else if (s == long || s == int || s == char || s == short) {
      /**/ if (type == 0) op = xmouse_point;
      else if (type == 1) op = xmouse_box;
      else if (type == 2) op = xmouse_line;
      else if (type == 3) op = xmouse_length;
      else error, "bad mouse selection type";
    } else {
      error, "bad mouse selection type must be a scalar integer of string";
    }
  } else {
    error, "bad mouse selection type";
  }
  if (am_subroutine()) op, win=win, system=system, prompt=prompt,
                         forever=forever;
  else return op(win=win, system=system, prompt=prompt, forever=forever);
}

func __xmouse(style, default_prompt)
/* DOCUMENT __xmouse()
     Private function used by xmouse and xmouse_[...] routines.
   SEE  ALSO: xmouse. */
{
  extern win, prompt, system, forever;
  oldwin = current_window();
  if (! is_void(win)) window, win;
  if (is_void(system)) system = -1;
  if (is_void(prompt)) prompt = default_prompt;
  for (;;) {
    m = mouse(system, style, prompt);
    if (is_array(m) && m(10)) break; /* m(10) is the  depressed button */
    if (! forever) {
      m = [];
      break;
    }
  }
  if (oldwin >= 0) window, oldwin;
  return m;
}

func xmouse_demo
/* DOCUMENT xmouse_demo
 *   Run a simple demonstration of the 'xmouse' functions.
 *
 * SEE ALSO:xmouse_point, xmouse_line, xmouse_box, xmouse_length.
 */
{
  t=xmouse_point();
  plp,t(2),t(1),symbol=PL_STAR,color=PL_BLUE,width=5,fill=0;

  t=xmouse_line();
  pldj,t(1),t(2),t(3),t(4),color=PL_GREEN;

  t=xmouse_box();
  pl_box,t,color=PL_GREEN;

  xmouse_length;
}

/*---------------------------------------------------------------------------*/
/* PLOTTING OF SIMPLE SHAPES */

func pl_box(xmin,xmax,ymin,ymax,color=,width=,type=)
/* DOCUMENT pl_box, xmin, xmax, ymin, ymax;
       -or- pl_box, [xmin, xmax, ymin, ymax];
     Plots a rectangular  box onto the current graphic  device.  As for plg
     (which see), keywords COLOR, WIDTH and TYPE can be specified.

   SEE ALSO: plg, pl_get_color, pl_cbox, pl_circle, pl_ellipse. */
{
  if (is_void(xmax) && numberof(xmin)==4) {
    ymax = xmin(4);
    ymin = xmin(3);
    xmax = xmin(2);
    xmin = xmin(1);
  }
  plg, [ymin,ymin,ymax,ymax], [xmin,xmax,xmax,xmin], closed=1,
    color=pl_get_color(color), width=width, type=type;
}

func pl_cbox(x0, y0, xsize, ysize, color=, width=, type=, legend=)
/* DOCUMENT pl_cbox, x0, y0, size;
       -or- pl_cbox, x0, y0, xsize, ysize;
     Draw a  SIZE by SIZE  square box or  a XSIZE by YSIZE  rectangular box
     centered around (X0,Y0).   Keywords COLOR, WIDTH and TYPE  can be used
     and have  the same  meaning as for  builtin routine "plg".  If keyword
     LEGEND is not set, an empty legend will be used.

   SEE ALSO plg, pl_box, pl_circle, pl_ellipse. */
{
  if (is_void(ysize)) ysize = xsize;
  if (is_void(legend)) legend=string(0);
  plg,
    y0 + ysize * [-0.5, +0.5, +0.5, -0.5],
    x0 + xsize * [-0.5, -0.5, +0.5, +0.5],
    color=color, width=width, type=type, marks=0, closed=1, legend=legend;
}

func pl_circle(x0, y0, r, color=, width=, number=, type=, legend=)
/* DOCUMENT pl_circle, x0, y0, r;
     Draw circle(s)  of radius R  around (X0,Y0).  Value of  keyword NUMBER
     tells how many segments to use to approximate the circle (default 20).
     Keywords COLOR, WIDTH  and TYPE can be used and  have the same meaning
     as for  "plg" routine (which  see). If keyword  LEGEND is not  set, an
     empty legend  will be  used.  Arguments may  be conformable  arrays to
     draw several circles in one call (but keywords must have scalar values
     if any).

   SEE ALSO plg, pl_box, pl_cbox, pl_ellipse. */
{
  if (is_void(legend)) legend = string();
  if (is_void(number)) number = 20;
  PI = 3.14159265358979323848;
  t = (2.0*PI/number)*indgen(number);
  cos_t = cos(t);
  sin_t = sin(t);
  if (is_void((dims = dimsof(x0, y0, r))))
    error, "non conformable arguments";
  if (dims(1)) {
    /* Draw several circles. */
    dummy = array(double, dims);
    x0 += dummy;
    y0 += dummy;
    r  += dummy;
    n = numberof(dummy);
    for (i=1 ; i<=n ; ++i) {
      plg, y0(i) + r(i)*cos_t, x0(i) + r(i)*sin_t,
        width=width, color=color, type=type, marks=0, closed=1, legend=legend;
    }
  } else {
    /* Draw a single circle. */
    plg, y0 + r*cos_t, x0 + r*sin_t,
      width=width, color=color, type=type, marks=0, closed=1, legend=legend;
  }
}

func pl_ellipse(x0, y0, a, b, theta, color=, width=, number=, type=, legend=)
/* DOCUMENT pl_ellipse, x0, y0, a, b, theta;
     Draw ellipse(s) centered at (X0,Y0)  with semi-axis A and B.  THETA is
     the angle (in  degrees counterclockwise) of the axis  of semi-length A
     with horizontal axis.  Value of keyword NUMBER tells how many segments
     to use to approximate the  circle (default 20).  Keywords COLOR, WIDTH
     and TYPE  can be used and have  the same meaning as  for "plg" routine
     (which see).   If keyword LEGEND is  not set, an empty  legend will be
     used.  Arguments may be conformable arrays to draw several ellipses in
     one call (but keywords must have scalar values if any).

   SEE ALSO plg, pl_box, pl_cbox, pl_circle. */
{
  if (is_void(legend)) legend = string();
  if (is_void(number)) number = 20;
  PI = 3.14159265358979323848;
  t = (2.0*PI/number)*indgen(number);
  cos_t = cos(t);
  sin_t = sin(t);
  if (is_void((dims = dimsof(x0, y0, a, b, theta))))
    error, "non conformable arguments";
  if (dims(1)) {
    /* Draw several ellipses. */
    dummy = array(double, dims);
    x0 += dummy;
    y0 += dummy;
    a  += dummy;
    b  += dummy;
    theta = (PI/180.0)*theta + dummy;
    n = numberof(dummy);
    for (i=1 ; i<=n ; ++i) {
      u = a(i)*cos_t;
      v = b(i)*sin_t;
      t = theta(i);
      cs = cos(t);
      sn = sin(t);
      plg, y0(i) + u*sn + v*cs, x0(i) + u*cs - v*sn,
        width=width, color=color, type=type, marks=0, closed=1, legend=legend;
    }
  } else {
    /* Draw a single ellipse. */
    theta *= (PI/180.0);
    u = a*cos_t;
    v = b*sin_t;
    cs = cos(theta);
    sn = sin(theta);
    plg, y0 + u*sn + v*cs, x0 + u*cs - v*sn,
      width=width, color=color, type=type, marks=0, closed=1, legend=legend;
  }
}

/*---------------------------------------------------------------------------*/
/* HINTON DIAGRAM */

func pl_hinton(a, x0, y0, x1, y1, cmin=, cmax=,
               bg=, neg=, pos=,
               pad=, debug=)
/* DOCUMENT pl_hinton, a;
         or pl_hinton, a, x0, y0, x1, y1;

     This subroutine plots a Hinton diagram which is a way of visualizing
     numerical values in the matrix/vector A, popular in the neural networks
     and machine learning literature.  Each element of A is represented by a
     square, the area occupied by the square is proportional to the magnitude
     of the value, and the colour (black/white by default) indicates its sign
     (negative/positive).

     Optional arguments X0, Y0, X1 and Y1 can be used to specify the
     coordinates of the centers of the first and last cells.

     Keywords BG, NEG and POS can be used to specify the colors of the
     background, the negative cells and the positive cells (respectively).
     See pl_get_color for means to specify a color.

     Keywodr PAD can be used to specify the minimum amount of padding between
     cells.  PAD is a relative value in the range 0 (no padding) to 1 (full
     padding).  The default value is PAD = 0.05.

     Keywords CMIN and CMAX can be used to specify the minimum and maximum
     values to consider.
     
   SEE ALSO: plf, pl_get_color.
 */
{
  if (is_void(pad)) {
    pad = 0.05;
  } else if (is_scalar(pad) && identof(pad) <= Y_DOUBLE && pad >= 0.0 && pad <= 1.0) {
    pad = double(pad);
  } else {
    error, "value of keyword PAD must be in (0,1)";
  }
  if (identof(a) > Y_DOUBLE || (d = dimsof(a))(1) != 2) {
    error, "expecting a 2-D real array";
  }
  nx = d(2);
  ny = d(3);  
  if (is_void(cmin)) cmin = min(a);
  if (is_void(cmax)) cmax = max(a);
  if (is_void(bg)) {
    /* default background color is dark gray */
    bg = array(char(100), 3);
  } else {
    bg = pl_get_color(bg, 1);
  }
  if (is_void(neg)) {
    /* default color of negative cells is black */
    neg = array(char(0), 3);
  } else {
    neg = pl_get_color(neg, 1);
  }
  if (is_void(pos)) {
    /* default color of positive cells is white */
    pos = array(char(255), 3);
  } else {
    pos = pl_get_color(pos, 1);
  }

  q = max(abs(cmin), abs(cmax));
  if (q == 0.0) {
    return;
  }
  q = (1.0 - pad)/(2.0*sqrt(q));
  t = sqrt(abs(min(max(cmin, a), cmax)));

  /* Fill the z-cells with the background color (gray) and then with
     foreground color (black/white) depending on the sign. */
  z = array(bg, 2*nx + 1, 2*ny + 1);
  i = where(a > 0.0);
  if (is_array(i)) {
    z(, ((--i)%nx)*2 + (i/nx)*(4*nx + 2) + (2*nx + 3)) = pos;
  }
  i = where(a < 0.0);
  if (is_array(i)) {
    z(, ((--i)%nx)*2 + (i/nx)*(4*nx + 2) + (2*nx + 3)) = neg;
  }
  //pli, z;

  if (is_void(x0)) x0 = 1.0;
  if (is_void(y0)) y0 = 1.0;
  if (is_void(x1)) x1 = x0 + nx - 1.0;
  if (is_void(y1)) y1 = y0 + ny - 1.0;
  xstp = (nx > 1 ? (x1 - x0)/(nx - 1.0) : 1.0);
  ystp = (ny > 1 ? (y1 - y0)/(ny - 1.0) : 1.0);
  
  xmin = x0 - 0.5*xstp;
  xmax = x1 + 0.5*xstp;
  ymin = y0 - 0.5*ystp;
  ymax = y1 + 0.5*ystp;
  
  lx = 2*nx + 2;
  ly = 2*ny + 2;

  i0 = 2:lx-2:2;
  i1 = 3:lx-1:2;
  j0 = 2:ly-2:2;
  j1 = 3:ly-1:2;
  
  s = (xstp*q)*t;
  c = pl_span(x0, x1, nx)(,-);
  a = c - s;
  b = c + s;
  x = array(double, lx, ly);
  x(1,) = xmin;
  x(0,) = xmax;
  x(i0, j0) = a;
  x(i0, j1) = a;
  x(i1, j0) = b;
  x(i1, j1) = b;
  x(,1) = x(,2);
  x(,0) = x(,-1);
 
  s = (ystp*q)*t;
  c = pl_span(y0, y1, ny)(-,);
  a = c - s;
  b = c + s;
  y = array(double, lx, ly);
  y(,1) = ymin;
  y(,0) = ymax;
  y(i0, j0) = a;
  y(i0, j1) = b;
  y(i1, j0) = a;
  y(i1, j1) = b;
  y(1,) = y(2,);
  y(0,) = y(-1,);

  if (debug) {
    plf, z, y, x, edges=1n;
  } else {
    plf, z, y, x, edges=0n;
  }
}

/*---------------------------------------------------------------------------*/
/* CONVERT POSTSCRIPT FILE INTO BITMAP IMAGE */

local ps2png, ps2jpeg, _ps2any_worker;
/* DOCUMENT ps2png, inp;
 *     -or- ps2png, inp, out;
 *     -or- ps2png(inp);
 *     -or- ps2png(inp, out);
 *     -or- ps2jpeg, inp;
 *     -or- ps2jpeg, inp, out;
 *     -or- ps2jpeg(inp);
 *     -or- ps2jpeg(inp, out); *
 *
 *   Convert PostScript or PDF file INP into a PNG or JPEG image.  If
 *   the name of the output file OUT is not specified, it get
 *   automatically derived from INP by replacing extension ".eps",
 *   ".ps" or ".pdf" by ".png" or ".jpg" (if INP do not have any of
 *   the extensions ".eps", ".ps" nor ".pdf", the ".png" or ".jpg"
 *   extension is simply appended).
 *
 *   When called as a subroutine, the conversion script is executed
 *   (by the system built-in function); when called as a function, the
 *   command is returned but not executed.
 *
 *   The conversion is handled by GhostScript and by various commands
 *   from the netpbm package.  The gs command and the netpbm commands
 *   must be found in the shell PATH.
 *
 *   The following keywords are supported:
 *
 *       GRAY = flag: convert to grayscale (default is color)?
 *
 *       DPI = resolution in pixel per inch (default 300).
 *
 *       BORDER = border in pixels (default 3).
 *
 *       ALPHABITS = number of antialiasing bits for both text and
 *           graphics: 0 (for none), 1, 2 or 4 (default is 0).
 *
 *       TEXTALPHABITS = same as ALPHABITS but for text only (default
 *           is ALPHABITS).
 *
 *       GRAPHICSALPHABITS = same as ALPHABITS but for graphics only
 *           (default is ALPHABITS).
 *
 *       SMOOTH = flag: smooth image (default is false)?
 *
 *       DEVICE = GhostScript output device to use (default is "ppmraw");
 *           it is probably better to not change this option.
 *
 *   The following keywords are only for conversion to PNG images:
 *
 *       TRANSPARENT = name of transparent color COLOR (default none);
 *
 *
 *   The following keywords are only for conversion to JPEG images:
 *
 *       QUALITY = quality (in percent) for JPEG image (default 85%).
 *
 *
 * SEE ALSO: system.
 */
func ps2png(inp, out, gray=, smooth=, dpi=, border=, device=,
            alphabits=, graphicsalphabits=, textalphabits=,
            transparent=)
{
  local status;
  command = _ps2any_worker(inp, out, 0);
  if (status) {
    error, command;
  }
  if (am_subroutine()) {
    system, command;
  } else {
    return command;
  }
}

func ps2jpeg(inp, out, gray=, smooth=, dpi=, border=, device=,
             alphabits=, graphicsalphabits=, textalphabits=,
             quality=)
{
  local status;
  command = _ps2any_worker(inp, out, 1);
  if (status) {
    error, command;
  }
  if (am_subroutine()) {
    system, command;
  } else {
    return command;
  }
}
func _ps2any_worker(inp, out, jpg)
{
  /* Status and all options as external variables. */
  extern status;
  extern alphabits, graphicsalphabits, textalphabits, quality;
  extern gray, smooth, dpi, border, device, transparent;

  status = -1; /* assume we have an error */

  /* Parse common options. */
  if (is_void(gray)) {
    gray = 0n;
  } else if (is_scalar(gray)) {
    gray = !(!gray);
  } else {
    return "bad value for GRAY keyword";
  }

  if (is_void(smooth)) {
    smooth = 0n;
  } else if (is_scalar(smooth)) {
    smooth = !(!smooth);
  } else {
    return "bad value for SMOOTH keyword";
  }

  if (is_void(dpi)) {
    dpi = 300.0;
  } else if (! is_scalar(dpi) ||
             ! (is_integer(dpi) || is_real(dpi)) ||
             dpi <= 0.0) {
    return "bad value for DPI keyword";
  } else {
    dpi = double(dpi);
  }

  if (is_void(border)) {
    border = 3.0;
  } else if (! is_scalar(border) ||
             ! (is_integer(border) || is_real(border)) ||
             border < 0.0) {
    return "bad value for BORDER keyword";
  } else {
    border = double(border);
  }

  if (is_void(device)) {
    device = "ppmraw";
  } else if (! is_scalar(device) || ! is_string(device)) {
    return "bad value for DEVICE keyword";
  }

  if (is_scalar(alphabits) && is_integer(alphabits) &&
      0 <= alphabits && alphabits <= 4) {
    alphabits = long(alphabits);
  } else if (is_void(alphabits)) {
    alphabits = 0;
  } else {
    return "bad value for keyword ALPHABITS";
  }
  if (is_scalar(textalphabits) && is_integer(textalphabits) &&
      0 <= textalphabits && textalphabits <= 4) {
    textalphabits = long(textalphabits);
  } else if (is_void(textalphabits)) {
    textalphabits  = alphabits;
  } else {
    return "bad value for keyword TEXTALPHABITS";
  }
  if (is_scalar(graphicsalphabits) && is_integer(graphicsalphabits) &&
      0 <= graphicsalphabits && graphicsalphabits <= 4) {
    graphicsalphabits = long(graphicsalphabits);
  } else if (is_void(graphicsalphabits)) {
    graphicsalphabits = alphabits;
  } else {
    return "bad value for keyword GRAPHICSALPHABITS";
  }

  /* Buil up the code. */
  command = swrite(format="gs -dSAFER -dNOPAUSE -dBATCH -q -sDEVICE=%s -r%.0f",
                   device, dpi);
  if (graphicsalphabits > 0) {
    command += swrite(format=" -dGraphicsAlphaBits=%d", graphicsalphabits);
  }
  if (textalphabits > 0) {
    command += swrite(format=" -dTextAlphaBits=%d", textalphabits);
  }
  command += swrite(format=" -sOutputFile=- \'%s\'", inp);
  if (smooth) {
    command += " | pnmalias";
  }
  command += " | pnmcrop";
  if (border > 0) {
    command += swrite(format=" | pnmpad -white -left %.0f -right %.0f -top %.0f -bottom %.0f", border, border, border, border);
  }
  if (is_void(out)) {
    out = strip_file_extension(inp, ".ps", ".eps") + (jpg ? ".jpg" : ".png");
  }
  if (jpg) {
    /* Convert to JPEG file. */
    if (is_void(quality)) {
      quality = 85.0;
    } else if (! is_scalar(quality) ||
               ! (is_integer(quality) || is_real(quality)) ||
               quality < 0.0 || quality > 100.0) {
      return "bad value for QUALITY keyword";
    } else {
      quality = double(quality);
    }
    command += swrite(format=" | pnmtojpeg > '%s' --optimize --quality=%.0f",
                      out, quality);
  } else {
    /* Convert to PNG file. */
    command += swrite(format=" | pnmtopng > '%s' -compression 9", out);
    if (is_scalar(transparent) && is_string(transparent)) {
      command += swrite(format=" -transparent '%s'", transparent);
    } else if (! is_void(transparent)) {
      return "bad value for TRANSPARENT keyword";
    }
  }

  /* Clear error status and return code. */
  status = 0;
  return command;
}

/*---------------------------------------------------------------------------*/
/* DUMP GRAPHICAL WINDOW AS BITMAP IMAGE */

local win2png, win2jpeg, _win2any_worker;
/* DOCUMENT win2png, filename;
 *     -or- win2png, format, value;
 *     -or- win2jpeg, filename;
 *     -or- win2jpeg, format, value;
 *
 *   Dump contents of current graphical window into a PNG or JPEG
 *   image.  FILENAME is the name of the output file.  If can also be
 *   specified by the two arguments FORMAT and VALUE, in which case
 *   FILENAME is computed by:
 *
 *       FILENAME = swrite(format=FORMAT, VALUE);
 *
 *   for instance:
 *
 *       win2png, "img-%04d.png", 3;
 *
 *   to save image into file "img-0003.png".
 *
 *   Keyword WIN can be set to specify another graphical window
 *   than the current one.
 *
 *   Keyword TEMPS can be set to specify the name of the temporary
 *   PostScript file.
 *
 *   If keyword DEBUG is true, the conversion script get printed out
 *   and the temporary PostScript file is not deleted (you can guess
 *   its name from the conversion script).
 *
 *   In addition, the same keywords as ps2png and ps2jpeg (which see)
 *   are available.  Note that, if not specified, a default value is
 *   provided for the DPI keyword of ps2png and ps2jpeg so that the
 *   resulting image has the same resolution as the graphical window
 *   (this however require Yeti plugin).
 *
 *
 * SEE ALSO: ps2png, ps2jpeg, hcps, tempfile, current_window, window.
 */
func win2png(filename, counter, win=, debug=, temp=,
             gray=, smooth=, dpi=, border=, device=,
             alphabits=, graphicsalphabits=, textalphabits=,
             transparent=)
{
  result = _win2any_worker(filename, counter, 0);
  if (result) {
    error, result;
  }
}

func win2jpeg(filename, counter, win=, debug=, temp=,
              gray=, smooth=, dpi=, border=, device=,
              alphabits=, graphicsalphabits=, textalphabits=,
              quality=)
{
  result = _win2any_worker(filename, counter, 1);
  if (result) {
    error, result;
  }
}

func _win2any_worker(filename, counter, jpg)
{
  /* All options as external variables. */
  extern temp, dpi, win;
  local inp, out, status;

  status = -1; /* assume we have an error */

  /* Prepare name of output file. */
  if (is_string(filename) && is_scalar(filename)) {
    if (is_void(counter)) {
      eq_nocopy, output, filename;
    } else {
      output = swrite(format=filename, counter);
    }
  } else {
    return "expecting a scalar string for FILENAME";
  }

  /* Figure out which window to dump. */
  old_win = current_window();
  if (is_void(win)) {
    tgt_win = old_win;
  } else if (window_exists(win)) {
    tgt_win = long(win);
  } else {
    tgt_win = -1;
  }
  if (tgt_win < 0) {
    return "bad graphical window";
  }
  if (is_void(dpi)) {
    if (! is_func(window_geometry)) {
      require, "yeti.i";
    }
    geom = window_geometry(tgt_win);
    dpi = long(floor(geom(1) + 0.5));
  }

  /* Temporary input PS file. */
  if (is_void(temp)) {
    temp = tempfile("/tmp/png_dump-XXXXXX.ps");
  } else if (! is_string(temp) || ! is_scalar(temp)) {
    return "expecting a scalar string for TEMP";
  }

  /* Build command. */
  script = _ps2any_worker(temp, output, jpg);
  if (status) {
    return script;
  }

  /* Dump graphics window. */
  window, tgt_win, hcp=temp, dump=1, legends=0, wait=1;
  hcp;
  window, tgt_win, hcp="";
  if (tgt_win != old_win && old_win != -1) {
    window, old_win;
  }
  pause, 1;
  if (debug) {
    write, format="%s\n", script;
  }
  system, script;
  if (! debug) {
    remove, temp;
  }
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func pl_span(a, b, n)
/* DOCUMENT pl_span(a, b, n)
     The function pl_span() returns an array of N doubles equally spaced from
     A to B.  This function is more simple than built-in span() -- in
     particular, A and B must be scalars -- but avoids rounding errors.

   SEE ALSO: span.
 */
{
  a += 0.0; // make sure A is floating point
  if (a == b || n == 1) return array(a, n);
  return ((b - a)/(n - 1.0))*(indgen(n) - (n + 1)/2.0) + (a + b)/2.0;
}

func pl_map(op, arg, default)
/* DOCUMENT pl_map(op, arg)
       -or- pl_map(op, arg, default)
     Maps scalar function OP onto array argument ARG to mimics element-wise
     unary operation.  Returns DEFAULT if ARG is void.
  
   SEE ALSO xwindow.
 */
{
  if (is_array(arg)) {
    /* use structof to avoid unecessary string duplication for string result */
    out = array(structof((out1 = op(arg(1)))), dimsof(arg));
    out(1) = out1;
    n = numberof(arg);
    for (i=2 ; i<=n ; ++i) out(i) = op(arg(i));
    return out;
  }
  if (is_void(arg)) return default;
  error, "unexpected argument";
}

func win_copy_lim(src, dst, axis)
/* DOCUMENT win_copy_lim, src, dst;
         or win_copy_lim, src, dst, axis;
  
     Make limits of window DST the same as those in window SRC.
  
     AXIS can be used to specify which limits to copy: 1st and 2nd bits
     respectively indicate whether X and/or Y limits should be copied. If
     unspecified, the limits of the two axis are copied (same as with AXIS=3).
  
   SEE ALSO: current_window, limits, window.
 */
{
  win = current_window();
  window, src;
  l = limits();
  window, dst;
  if (is_void(axis) || (axis & 3) == 3) {
    /* copy X and Y limits */
    limits, l(1), l(2), l(3), l(4);
  } else if ((axis & 3) == 1) {
    /* only copy X limits */
    limits, l(1), l(2);
  } else if ((axis & 3) == 2) {
    /* copy Y limits */
    limits, , , l(3), l(4);
  }
  if (win >= 0) window, win; /* restore old current window */
}

/*---------------------------------------------------------------------------*/
/* COLOR DATABASE */

func pl_database_index_to_rgb(color)
{
  return _PL_COLOR_RGB(, color);
}
func pl_database_index_to_packed(color)
{
  return _PL_COLOR_PACKED(color);
}
func pl_rgb_to_packed(color)
{
  return (0x01000000 | long(color(1,..)) |
          (long(color(2,..)) << 8) | (long(color(3,..)) << 16));
}
func pl_packed_to_rgb(color)
{
  rgb = array(char, 3, dimsof(color));
  rgb(1, ..) = (color & 0xff);
  rgb(2, ..) = ((color >> 8) & 0xff);
  rgb(3, ..) = ((color >> 16) & 0xff);
  return rgb;
}

/* Names of known colors (the 16 indexed ones must come first). */
_PL_COLOR_NAMES = \
["bg", "fg", "black", "white", "red", "green", "blue", "cyan", "magenta",
 "yellow", "grayd", "grayc", "grayb", "graya", "extra", "xor",
 "aliceblue", "antiquewhite", "antiquewhite1", "antiquewhite2",
 "antiquewhite3", "antiquewhite4", "aquamarine", "aquamarine1", "aquamarine2",
 "aquamarine3", "aquamarine4", "azure", "azure1", "azure2", "azure3", "azure4",
 "beige", "bisque", "bisque1", "bisque2", "bisque3", "bisque4",
 "blanchedalmond", "blue1", "blue2", "blue3", "blue4", "blueviolet", "brown",
 "brown1", "brown2", "brown3", "brown4", "burlywood", "burlywood1",
 "burlywood2", "burlywood3", "burlywood4", "cadetblue", "cadetblue1",
 "cadetblue2", "cadetblue3", "cadetblue4", "chartreuse", "chartreuse1",
 "chartreuse2", "chartreuse3", "chartreuse4", "chocolate", "chocolate1",
 "chocolate2", "chocolate3", "chocolate4", "coral", "coral1", "coral2",
 "coral3", "coral4", "cornflowerblue", "cornsilk", "cornsilk1", "cornsilk2",
 "cornsilk3", "cornsilk4", "cyan1", "cyan2", "cyan3", "cyan4", "darkblue",
 "darkcyan", "darkgoldenrod", "darkgoldenrod1", "darkgoldenrod2",
 "darkgoldenrod3", "darkgoldenrod4", "darkgray", "darkgreen", "darkgrey",
 "darkkhaki", "darkmagenta", "darkolivegreen", "darkolivegreen1",
 "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkorange",
 "darkorange1", "darkorange2", "darkorange3", "darkorange4", "darkorchid",
 "darkorchid1", "darkorchid2", "darkorchid3", "darkorchid4", "darkred",
 "darksalmon", "darkseagreen", "darkseagreen1", "darkseagreen2",
 "darkseagreen3", "darkseagreen4", "darkslateblue", "darkslategray",
 "darkslategray1", "darkslategray2", "darkslategray3", "darkslategray4",
 "darkslategrey", "darkturquoise", "darkviolet", "deeppink", "deeppink1",
 "deeppink2", "deeppink3", "deeppink4", "deepskyblue", "deepskyblue1",
 "deepskyblue2", "deepskyblue3", "deepskyblue4", "dimgray", "dimgrey",
 "dodgerblue", "dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4",
 "firebrick", "firebrick1", "firebrick2", "firebrick3", "firebrick4",
 "floralwhite", "forestgreen", "gainsboro", "ghostwhite", "gold", "gold1",
 "gold2", "gold3", "gold4", "goldenrod", "goldenrod1", "goldenrod2",
 "goldenrod3", "goldenrod4", "gray", "gray0", "gray1", "gray10", "gray100",
 "gray11", "gray12", "gray13", "gray14", "gray15", "gray16", "gray17",
 "gray18", "gray19", "gray2", "gray20", "gray21", "gray22", "gray23", "gray24",
 "gray25", "gray26", "gray27", "gray28", "gray29", "gray3", "gray30", "gray31",
 "gray32", "gray33", "gray34", "gray35", "gray36", "gray37", "gray38",
 "gray39", "gray4", "gray40", "gray41", "gray42", "gray43", "gray44", "gray45",
 "gray46", "gray47", "gray48", "gray49", "gray5", "gray50", "gray51", "gray52",
 "gray53", "gray54", "gray55", "gray56", "gray57", "gray58", "gray59", "gray6",
 "gray60", "gray61", "gray62", "gray63", "gray64", "gray65", "gray66",
 "gray67", "gray68", "gray69", "gray7", "gray70", "gray71", "gray72", "gray73",
 "gray74", "gray75", "gray76", "gray77", "gray78", "gray79", "gray8", "gray80",
 "gray81", "gray82", "gray83", "gray84", "gray85", "gray86", "gray87",
 "gray88", "gray89", "gray9", "gray90", "gray91", "gray92", "gray93", "gray94",
 "gray95", "gray96", "gray97", "gray98", "gray99", "green1", "green2",
 "green3", "green4", "greenyellow", "grey", "grey0", "grey1", "grey10",
 "grey100", "grey11", "grey12", "grey13", "grey14", "grey15", "grey16",
 "grey17", "grey18", "grey19", "grey2", "grey20", "grey21", "grey22", "grey23",
 "grey24", "grey25", "grey26", "grey27", "grey28", "grey29", "grey3", "grey30",
 "grey31", "grey32", "grey33", "grey34", "grey35", "grey36", "grey37",
 "grey38", "grey39", "grey4", "grey40", "grey41", "grey42", "grey43", "grey44",
 "grey45", "grey46", "grey47", "grey48", "grey49", "grey5", "grey50", "grey51",
 "grey52", "grey53", "grey54", "grey55", "grey56", "grey57", "grey58",
 "grey59", "grey6", "grey60", "grey61", "grey62", "grey63", "grey64", "grey65",
 "grey66", "grey67", "grey68", "grey69", "grey7", "grey70", "grey71", "grey72",
 "grey73", "grey74", "grey75", "grey76", "grey77", "grey78", "grey79", "grey8",
 "grey80", "grey81", "grey82", "grey83", "grey84", "grey85", "grey86",
 "grey87", "grey88", "grey89", "grey9", "grey90", "grey91", "grey92", "grey93",
 "grey94", "grey95", "grey96", "grey97", "grey98", "grey99", "honeydew",
 "honeydew1", "honeydew2", "honeydew3", "honeydew4", "hotpink", "hotpink1",
 "hotpink2", "hotpink3", "hotpink4", "indianred", "indianred1", "indianred2",
 "indianred3", "indianred4", "ivory", "ivory1", "ivory2", "ivory3", "ivory4",
 "khaki", "khaki1", "khaki2", "khaki3", "khaki4", "lavender", "lavenderblush",
 "lavenderblush1", "lavenderblush2", "lavenderblush3", "lavenderblush4",
 "lawngreen", "lemonchiffon", "lemonchiffon1", "lemonchiffon2",
 "lemonchiffon3", "lemonchiffon4", "lightblue", "lightblue1", "lightblue2",
 "lightblue3", "lightblue4", "lightcoral", "lightcyan", "lightcyan1",
 "lightcyan2", "lightcyan3", "lightcyan4", "lightgoldenrod", "lightgoldenrod1",
 "lightgoldenrod2", "lightgoldenrod3", "lightgoldenrod4",
 "lightgoldenrodyellow", "lightgray", "lightgreen", "lightgrey", "lightpink",
 "lightpink1", "lightpink2", "lightpink3", "lightpink4", "lightsalmon",
 "lightsalmon1", "lightsalmon2", "lightsalmon3", "lightsalmon4",
 "lightseagreen", "lightskyblue", "lightskyblue1", "lightskyblue2",
 "lightskyblue3", "lightskyblue4", "lightslateblue", "lightslategray",
 "lightslategrey", "lightsteelblue", "lightsteelblue1", "lightsteelblue2",
 "lightsteelblue3", "lightsteelblue4", "lightyellow", "lightyellow1",
 "lightyellow2", "lightyellow3", "lightyellow4", "limegreen", "linen",
 "magenta1", "magenta2", "magenta3", "magenta4", "maroon", "maroon1",
 "maroon2", "maroon3", "maroon4", "mediumaquamarine", "mediumblue",
 "mediumorchid", "mediumorchid1", "mediumorchid2", "mediumorchid3",
 "mediumorchid4", "mediumpurple", "mediumpurple1", "mediumpurple2",
 "mediumpurple3", "mediumpurple4", "mediumseagreen", "mediumslateblue",
 "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue",
 "mintcream", "mistyrose", "mistyrose1", "mistyrose2", "mistyrose3",
 "mistyrose4", "moccasin", "navajowhite", "navajowhite1", "navajowhite2",
 "navajowhite3", "navajowhite4", "navy", "navyblue", "oldlace", "olivedrab",
 "olivedrab1", "olivedrab2", "olivedrab3", "olivedrab4", "orange", "orange1",
 "orange2", "orange3", "orange4", "orangered", "orangered1", "orangered2",
 "orangered3", "orangered4", "orchid", "orchid1", "orchid2", "orchid3",
 "orchid4", "palegoldenrod", "palegreen", "palegreen1", "palegreen2",
 "palegreen3", "palegreen4", "paleturquoise", "paleturquoise1",
 "paleturquoise2", "paleturquoise3", "paleturquoise4", "palevioletred",
 "palevioletred1", "palevioletred2", "palevioletred3", "palevioletred4",
 "papayawhip", "peachpuff", "peachpuff1", "peachpuff2", "peachpuff3",
 "peachpuff4", "peru", "pink", "pink1", "pink2", "pink3", "pink4", "plum",
 "plum1", "plum2", "plum3", "plum4", "powderblue", "purple", "purple1",
 "purple2", "purple3", "purple4", "red1", "red2", "red3", "red4", "rosybrown",
 "rosybrown1", "rosybrown2", "rosybrown3", "rosybrown4", "royalblue",
 "royalblue1", "royalblue2", "royalblue3", "royalblue4", "saddlebrown",
 "salmon", "salmon1", "salmon2", "salmon3", "salmon4", "sandybrown",
 "seagreen", "seagreen1", "seagreen2", "seagreen3", "seagreen4", "seashell",
 "seashell1", "seashell2", "seashell3", "seashell4", "sienna", "sienna1",
 "sienna2", "sienna3", "sienna4", "skyblue", "skyblue1", "skyblue2",
 "skyblue3", "skyblue4", "slateblue", "slateblue1", "slateblue2", "slateblue3",
 "slateblue4", "slategray", "slategray1", "slategray2", "slategray3",
 "slategray4", "slategrey", "snow", "snow1", "snow2", "snow3", "snow4",
 "springgreen", "springgreen1", "springgreen2", "springgreen3", "springgreen4",
 "steelblue", "steelblue1", "steelblue2", "steelblue3", "steelblue4", "tan",
 "tan1", "tan2", "tan3", "tan4", "thistle", "thistle1", "thistle2", "thistle3",
 "thistle4", "tomato", "tomato1", "tomato2", "tomato3", "tomato4", "turquoise",
 "turquoise1", "turquoise2", "turquoise3", "turquoise4", "violet", "violetred",
 "violetred1", "violetred2", "violetred3", "violetred4", "wheat", "wheat1",
 "wheat2", "wheat3", "wheat4", "whitesmoke", "yellow1", "yellow2", "yellow3",
 "yellow4", "yellowgreen"];
_PL_COLOR_RGB = \
[[255,255,255], [  0,  0,  0], [  0,  0,  0], [255,255,255], [255,  0,  0],
 [  0,255,  0], [  0,  0,255], [  0,255,255], [255,  0,255], [255,255,  0],
 [100,100,100], [150,150,150], [190,190,190], [214,214,214], [  0,  0,  0],
 [  0,  0,  0], [240,248,255], [250,235,215], [255,239,219], [238,223,204],
 [205,192,176], [139,131,120], [127,255,212], [127,255,212], [118,238,198],
 [102,205,170], [ 69,139,116], [240,255,255], [240,255,255], [224,238,238],
 [193,205,205], [131,139,139], [245,245,220], [255,228,196], [255,228,196],
 [238,213,183], [205,183,158], [139,125,107], [255,235,205], [  0,  0,255],
 [  0,  0,238], [  0,  0,205], [  0,  0,139], [138, 43,226], [165, 42, 42],
 [255, 64, 64], [238, 59, 59], [205, 51, 51], [139, 35, 35], [222,184,135],
 [255,211,155], [238,197,145], [205,170,125], [139,115, 85], [ 95,158,160],
 [152,245,255], [142,229,238], [122,197,205], [ 83,134,139], [127,255,  0],
 [127,255,  0], [118,238,  0], [102,205,  0], [ 69,139,  0], [210,105, 30],
 [255,127, 36], [238,118, 33], [205,102, 29], [139, 69, 19], [255,127, 80],
 [255,114, 86], [238,106, 80], [205, 91, 69], [139, 62, 47], [100,149,237],
 [255,248,220], [255,248,220], [238,232,205], [205,200,177], [139,136,120],
 [  0,255,255], [  0,238,238], [  0,205,205], [  0,139,139], [  0,  0,139],
 [  0,139,139], [184,134, 11], [255,185, 15], [238,173, 14], [205,149, 12],
 [139,101,  8], [169,169,169], [  0,100,  0], [169,169,169], [189,183,107],
 [139,  0,139], [ 85,107, 47], [202,255,112], [188,238,104], [162,205, 90],
 [110,139, 61], [255,140,  0], [255,127,  0], [238,118,  0], [205,102,  0],
 [139, 69,  0], [153, 50,204], [191, 62,255], [178, 58,238], [154, 50,205],
 [104, 34,139], [139,  0,  0], [233,150,122], [143,188,143], [193,255,193],
 [180,238,180], [155,205,155], [105,139,105], [ 72, 61,139], [ 47, 79, 79],
 [151,255,255], [141,238,238], [121,205,205], [ 82,139,139], [ 47, 79, 79],
 [  0,206,209], [148,  0,211], [255, 20,147], [255, 20,147], [238, 18,137],
 [205, 16,118], [139, 10, 80], [  0,191,255], [  0,191,255], [  0,178,238],
 [  0,154,205], [  0,104,139], [105,105,105], [105,105,105], [ 30,144,255],
 [ 30,144,255], [ 28,134,238], [ 24,116,205], [ 16, 78,139], [178, 34, 34],
 [255, 48, 48], [238, 44, 44], [205, 38, 38], [139, 26, 26], [255,250,240],
 [ 34,139, 34], [220,220,220], [248,248,255], [255,215,  0], [255,215,  0],
 [238,201,  0], [205,173,  0], [139,117,  0], [218,165, 32], [255,193, 37],
 [238,180, 34], [205,155, 29], [139,105, 20], [190,190,190], [  0,  0,  0],
 [  3,  3,  3], [ 26, 26, 26], [255,255,255], [ 28, 28, 28], [ 31, 31, 31],
 [ 33, 33, 33], [ 36, 36, 36], [ 38, 38, 38], [ 41, 41, 41], [ 43, 43, 43],
 [ 46, 46, 46], [ 48, 48, 48], [  5,  5,  5], [ 51, 51, 51], [ 54, 54, 54],
 [ 56, 56, 56], [ 59, 59, 59], [ 61, 61, 61], [ 64, 64, 64], [ 66, 66, 66],
 [ 69, 69, 69], [ 71, 71, 71], [ 74, 74, 74], [  8,  8,  8], [ 77, 77, 77],
 [ 79, 79, 79], [ 82, 82, 82], [ 84, 84, 84], [ 87, 87, 87], [ 89, 89, 89],
 [ 92, 92, 92], [ 94, 94, 94], [ 97, 97, 97], [ 99, 99, 99], [ 10, 10, 10],
 [102,102,102], [105,105,105], [107,107,107], [110,110,110], [112,112,112],
 [115,115,115], [117,117,117], [120,120,120], [122,122,122], [125,125,125],
 [ 13, 13, 13], [127,127,127], [130,130,130], [133,133,133], [135,135,135],
 [138,138,138], [140,140,140], [143,143,143], [145,145,145], [148,148,148],
 [150,150,150], [ 15, 15, 15], [153,153,153], [156,156,156], [158,158,158],
 [161,161,161], [163,163,163], [166,166,166], [168,168,168], [171,171,171],
 [173,173,173], [176,176,176], [ 18, 18, 18], [179,179,179], [181,181,181],
 [184,184,184], [186,186,186], [189,189,189], [191,191,191], [194,194,194],
 [196,196,196], [199,199,199], [201,201,201], [ 20, 20, 20], [204,204,204],
 [207,207,207], [209,209,209], [212,212,212], [214,214,214], [217,217,217],
 [219,219,219], [222,222,222], [224,224,224], [227,227,227], [ 23, 23, 23],
 [229,229,229], [232,232,232], [235,235,235], [237,237,237], [240,240,240],
 [242,242,242], [245,245,245], [247,247,247], [250,250,250], [252,252,252],
 [  0,255,  0], [  0,238,  0], [  0,205,  0], [  0,139,  0], [173,255, 47],
 [190,190,190], [  0,  0,  0], [  3,  3,  3], [ 26, 26, 26], [255,255,255],
 [ 28, 28, 28], [ 31, 31, 31], [ 33, 33, 33], [ 36, 36, 36], [ 38, 38, 38],
 [ 41, 41, 41], [ 43, 43, 43], [ 46, 46, 46], [ 48, 48, 48], [  5,  5,  5],
 [ 51, 51, 51], [ 54, 54, 54], [ 56, 56, 56], [ 59, 59, 59], [ 61, 61, 61],
 [ 64, 64, 64], [ 66, 66, 66], [ 69, 69, 69], [ 71, 71, 71], [ 74, 74, 74],
 [  8,  8,  8], [ 77, 77, 77], [ 79, 79, 79], [ 82, 82, 82], [ 84, 84, 84],
 [ 87, 87, 87], [ 89, 89, 89], [ 92, 92, 92], [ 94, 94, 94], [ 97, 97, 97],
 [ 99, 99, 99], [ 10, 10, 10], [102,102,102], [105,105,105], [107,107,107],
 [110,110,110], [112,112,112], [115,115,115], [117,117,117], [120,120,120],
 [122,122,122], [125,125,125], [ 13, 13, 13], [127,127,127], [130,130,130],
 [133,133,133], [135,135,135], [138,138,138], [140,140,140], [143,143,143],
 [145,145,145], [148,148,148], [150,150,150], [ 15, 15, 15], [153,153,153],
 [156,156,156], [158,158,158], [161,161,161], [163,163,163], [166,166,166],
 [168,168,168], [171,171,171], [173,173,173], [176,176,176], [ 18, 18, 18],
 [179,179,179], [181,181,181], [184,184,184], [186,186,186], [189,189,189],
 [191,191,191], [194,194,194], [196,196,196], [199,199,199], [201,201,201],
 [ 20, 20, 20], [204,204,204], [207,207,207], [209,209,209], [212,212,212],
 [214,214,214], [217,217,217], [219,219,219], [222,222,222], [224,224,224],
 [227,227,227], [ 23, 23, 23], [229,229,229], [232,232,232], [235,235,235],
 [237,237,237], [240,240,240], [242,242,242], [245,245,245], [247,247,247],
 [250,250,250], [252,252,252], [240,255,240], [240,255,240], [224,238,224],
 [193,205,193], [131,139,131], [255,105,180], [255,110,180], [238,106,167],
 [205, 96,144], [139, 58, 98], [205, 92, 92], [255,106,106], [238, 99, 99],
 [205, 85, 85], [139, 58, 58], [255,255,240], [255,255,240], [238,238,224],
 [205,205,193], [139,139,131], [240,230,140], [255,246,143], [238,230,133],
 [205,198,115], [139,134, 78], [230,230,250], [255,240,245], [255,240,245],
 [238,224,229], [205,193,197], [139,131,134], [124,252,  0], [255,250,205],
 [255,250,205], [238,233,191], [205,201,165], [139,137,112], [173,216,230],
 [191,239,255], [178,223,238], [154,192,205], [104,131,139], [240,128,128],
 [224,255,255], [224,255,255], [209,238,238], [180,205,205], [122,139,139],
 [238,221,130], [255,236,139], [238,220,130], [205,190,112], [139,129, 76],
 [250,250,210], [211,211,211], [144,238,144], [211,211,211], [255,182,193],
 [255,174,185], [238,162,173], [205,140,149], [139, 95,101], [255,160,122],
 [255,160,122], [238,149,114], [205,129, 98], [139, 87, 66], [ 32,178,170],
 [135,206,250], [176,226,255], [164,211,238], [141,182,205], [ 96,123,139],
 [132,112,255], [119,136,153], [119,136,153], [176,196,222], [202,225,255],
 [188,210,238], [162,181,205], [110,123,139], [255,255,224], [255,255,224],
 [238,238,209], [205,205,180], [139,139,122], [ 50,205, 50], [250,240,230],
 [255,  0,255], [238,  0,238], [205,  0,205], [139,  0,139], [176, 48, 96],
 [255, 52,179], [238, 48,167], [205, 41,144], [139, 28, 98], [102,205,170],
 [  0,  0,205], [186, 85,211], [224,102,255], [209, 95,238], [180, 82,205],
 [122, 55,139], [147,112,219], [171,130,255], [159,121,238], [137,104,205],
 [ 93, 71,139], [ 60,179,113], [123,104,238], [  0,250,154], [ 72,209,204],
 [199, 21,133], [ 25, 25,112], [245,255,250], [255,228,225], [255,228,225],
 [238,213,210], [205,183,181], [139,125,123], [255,228,181], [255,222,173],
 [255,222,173], [238,207,161], [205,179,139], [139,121, 94], [  0,  0,128],
 [  0,  0,128], [253,245,230], [107,142, 35], [192,255, 62], [179,238, 58],
 [154,205, 50], [105,139, 34], [255,165,  0], [255,165,  0], [238,154,  0],
 [205,133,  0], [139, 90,  0], [255, 69,  0], [255, 69,  0], [238, 64,  0],
 [205, 55,  0], [139, 37,  0], [218,112,214], [255,131,250], [238,122,233],
 [205,105,201], [139, 71,137], [238,232,170], [152,251,152], [154,255,154],
 [144,238,144], [124,205,124], [ 84,139, 84], [175,238,238], [187,255,255],
 [174,238,238], [150,205,205], [102,139,139], [219,112,147], [255,130,171],
 [238,121,159], [205,104,137], [139, 71, 93], [255,239,213], [255,218,185],
 [255,218,185], [238,203,173], [205,175,149], [139,119,101], [205,133, 63],
 [255,192,203], [255,181,197], [238,169,184], [205,145,158], [139, 99,108],
 [221,160,221], [255,187,255], [238,174,238], [205,150,205], [139,102,139],
 [176,224,230], [160, 32,240], [155, 48,255], [145, 44,238], [125, 38,205],
 [ 85, 26,139], [255,  0,  0], [238,  0,  0], [205,  0,  0], [139,  0,  0],
 [188,143,143], [255,193,193], [238,180,180], [205,155,155], [139,105,105],
 [ 65,105,225], [ 72,118,255], [ 67,110,238], [ 58, 95,205], [ 39, 64,139],
 [139, 69, 19], [250,128,114], [255,140,105], [238,130, 98], [205,112, 84],
 [139, 76, 57], [244,164, 96], [ 46,139, 87], [ 84,255,159], [ 78,238,148],
 [ 67,205,128], [ 46,139, 87], [255,245,238], [255,245,238], [238,229,222],
 [205,197,191], [139,134,130], [160, 82, 45], [255,130, 71], [238,121, 66],
 [205,104, 57], [139, 71, 38], [135,206,235], [135,206,255], [126,192,238],
 [108,166,205], [ 74,112,139], [106, 90,205], [131,111,255], [122,103,238],
 [105, 89,205], [ 71, 60,139], [112,128,144], [198,226,255], [185,211,238],
 [159,182,205], [108,123,139], [112,128,144], [255,250,250], [255,250,250],
 [238,233,233], [205,201,201], [139,137,137], [  0,255,127], [  0,255,127],
 [  0,238,118], [  0,205,102], [  0,139, 69], [ 70,130,180], [ 99,184,255],
 [ 92,172,238], [ 79,148,205], [ 54,100,139], [210,180,140], [255,165, 79],
 [238,154, 73], [205,133, 63], [139, 90, 43], [216,191,216], [255,225,255],
 [238,210,238], [205,181,205], [139,123,139], [255, 99, 71], [255, 99, 71],
 [238, 92, 66], [205, 79, 57], [139, 54, 38], [ 64,224,208], [  0,245,255],
 [  0,229,238], [  0,197,205], [  0,134,139], [238,130,238], [208, 32,144],
 [255, 62,150], [238, 58,140], [205, 50,120], [139, 34, 82], [245,222,179],
 [255,231,186], [238,216,174], [205,186,150], [139,126,102], [245,245,245],
 [255,255,  0], [238,238,  0], [205,205,  0], [139,139,  0], [154,205, 50]];
#if 0
_PL_COLOR_PACKED = (0x01000000 | _PL_COLOR_RGB(1,) |
                    (_PL_COLOR_RGB(2,) << 8) | (_PL_COLOR_RGB(3,) << 16));
#endif
_PL_COLOR_RGB = char(_PL_COLOR_RGB);

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
