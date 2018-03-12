/*
 * mira2_tests.i -
 *
 * Test suite for MiRA.
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

TEST_DIR = fulldirname(current_include());
include, TEST_DIR+"mira2.i";

main = mira_new(TEST_DIR+"../Contest1_J.oifits");
pixelsize = "0.3mas";
//h = collect_coords(main);
//u = uniquevalues(h.u);
//v = uniquevalues(h.v);

func model(main, force=, pixelsize=, xform=, smearingfunction=,
           smearingfactor=)
{
  MIRA_DEBUG = 0n;
  MIRA_QUIET = 1n;
  if (is_void(smearingfactor)) smearingfactor = 1;
  mira_config, main, pixelsize=pixelsize, xform=xform,
    smearingfunction=smearingfunction, smearingfactor=smearingfactor;
  mira_update, main, force=force;
  return main.xform;
}

func reldif(a, b)
{
  return (a == b ? 0.0 : 2.0*abs(a - b)/(abs(a) + abs(b)));
}

_INFO_STYLE    = ansi_term(ANSI_TERM_BOLD,ANSI_TERM_FG_BLUE);
_SUCCESS_STYLE = ansi_term(ANSI_TERM_BOLD,ANSI_TERM_FG_GREEN);
_FAILURE_STYLE = ansi_term(ANSI_TERM_BOLD,ANSI_TERM_FG_RED);
_RESET_STYLE   = ansi_term(ANSI_TERM_RESET);
nerrors = 0;
func report(success,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
{
  str = printf(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
  if (strpart(str, 0:0) == "\n") {
    str = strpart(str, 1:-1);
  }
  if (success) {
    write, format="%sSUCCESS - %s%s\n", _SUCCESS_STYLE, str, _RESET_STYLE;
  } else {
    write, format="%sFAILURE - %s%s\n", _FAILURE_STYLE, str, _RESET_STYLE;
    ++nerrors;
  }
}

func start(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
{
  str = printf(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
  write, format="%s%s", _INFO_STYLE, str;
  tic;
}

func stop(n)
{
  t = toc(n);
  write, format=" %.5g µs/evaluation%s\n", t(3)*1e6, _RESET_STYLE;
}

func testsymbols
{
  names = symbol_names(-1);
  sel = strgrep("^_?(mira|MIRA)_", names);
  names = names(where(sel(1,..) <= sel(2,..)));
  n = numberof(names);
  if (n > 0) {
    names = names(sort(names));
    for (i = 1; i <= n; ++i) {
      name = names(i);
      if (! symbol_exists(name)) {
        warn, "symbol \"%s\" is undefined", name;
      } else if (is_void(symbol_def(name))) {
        warn, "symbol \"%s\" is void", name;
      }
    }
  }
}
testsymbols;

func benchmark_build(main, pixelsize=)
{
  repeat = 100;

  k = repeat+1;
  start, "build NFFT operator:                          ";
  while(--k) {
    H = model(main, force=1, pixelsize=pixelsize, xform="nfft",
              smearingfunction="none", smearingfactor=1);
  }
  stop, repeat;

  k = repeat+1;
  start, "build separable operator:                     ";
  while(--k) {
    H = model(main, force=1, pixelsize=pixelsize, xform="separable",
              smearingfunction="none", smearingfactor=1);
  }
  stop, repeat;

  repeat = 20;
  k = repeat+1;
  start, "build nonseparable operator (smearing=none):  ";
  while(--k) {
    H = model(main, force=1, pixelsize=pixelsize, xform="nonseparable",
              smearingfunction="none", smearingfactor=1);
  }
  stop, repeat;

  k = repeat+1;
  start, "build nonseparable operator (smearing=sinc):  ";
  while(--k) {
    H = model(main, force=1, pixelsize=pixelsize, xform="nonseparable",
              smearingfunction="sinc", smearingfactor=1);
  }
  stop, repeat;

  k = repeat+1;
  start, "build nonseparable operator (smearing=gauss): ";
  while(--k) {
    H = model(main, force=1, pixelsize=pixelsize, xform="nonseparable",
              smearingfunction="gauss", smearingfactor=1);
  }
  stop, repeat;
}

func test_xform(main, pixelsize=)
{
  /* Without smearing. */
  write, format="\n%s%s%s\n", _INFO_STYLE,
    "Tests without smearing", _RESET_STYLE;

  H0 = model(main, pixelsize=pixelsize, xform="nfft",
             smearingfunction="none", smearingfactor=1);
  H1 = model(main, pixelsize=pixelsize, xform="nonseparable",
             smearingfunction="none", smearingfactor=1);
  H2 = model(main, pixelsize=pixelsize, xform="separable",
             smearingfunction="none", smearingfactor=1);

  img = random(mira_image_size(main));
  img *= 1.0/sum(img);

  z0 = H0(img);
  z1 = H1(img);
  z2 = H2(img);

  z = random_n(dimsof(z2));
  g0 = H0(z, 1);
  g1 = H1(z, 1);
  g2 = H2(z, 1);

  x = random_n(dimsof(img));
  y = random_n(dimsof(H0(x)));

  res01 = sum(y*H0(x));
  res02 = sum(H0(y,1)*x);
  err0 = reldif(res01, res02);
  report, err0 < 1e-10, "direct/adjoint inner product test for NFFT operator (%g)", err0;

  res11 = sum(y*H1(x));
  res12 = sum(H1(y,1)*x);
  err1 = reldif(res11, res12);
  report, err1 < 1e-10, "direct/adjoint inner product test for nonseparable operator (%g)", err1;

  res21 = sum(y*H2(x));
  res22 = sum(H2(y,1)*x);
  err2 = reldif(res21, res22);
  report, err2 < 1e-10, "direct/adjoint inner product test for separable operator (%g)", err2;

  err = max(abs(z0 - z1));
  report, err < 1e-10, "direct nfft vs. nonseparable (%g)", err;

  err = max(abs(g0 - g1));
  report, err < 1e-10, "adjoint nfft vs. nonseparable (%g)", err;

  err = max(abs(z2 - z1));
  report, err < 1e-10, "separable vs. nonseparable (%g)", err;

  err = max(abs(g2 - g1));
  report, err < 1e-10, "adjoint separable vs. nonseparable (%g)", err;

  repeat = 100;

  k = repeat+1; start,"apply direct NFFT operator:         "; while(--k)z0=H0(img);stop,repeat;
  k = repeat+1; start,"apply direct nonseparable operator: "; while(--k)z1=H1(img);stop,repeat;
  k = repeat+1; start,"apply direct separable operator:    "; while(--k)z2=H2(img);stop,repeat;

  k = repeat+1; start,"apply direct NFFT operator:         "; while(--k)z0=H0(img);stop,repeat;
  k = repeat+1; start,"apply direct nonseparable operator: "; while(--k)z1=H1(img);stop,repeat;
  k = repeat+1; start,"apply direct separable operator:    "; while(--k)z2=H2(img);stop,repeat;

  k = repeat+1; start,"apply adjoint NFFT operator:        "; while(--k)g0=H0(z,1);stop,repeat;
  k = repeat+1; start,"apply adjoint nonseparable operator:"; while(--k)g1=H1(z,1);stop,repeat;
  k = repeat+1; start,"apply adjoint separable operator:   "; while(--k)g2=H2(z,1);stop,repeat;

  /* With smearing. */
  write, format="\n%s%s%s\n", _INFO_STYLE,
                "Tests with smearing", _RESET_STYLE;

  S1 = model(main, pixelsize=pixelsize, xform="nonseparable",
             smearingfunction="sinc", smearingfactor=1);
  S2 = model(main, pixelsize=pixelsize, xform="nonseparable",
             smearingfunction="gauss", smearingfactor=1);

  z1 = S1(img);
  z2 = S2(img);

  err = max(abs(z2 - z1));
  report, err < 0.05, "direct nonseparable, sinc vs. gauss (%g)", err;

  z = (z1 + z2)/2;
  g1 = S1(z, 1);
  g2 = S2(z, 1);

  err = max(abs(g2 - g1));
  report, err < 0.05, "adjoint nonseparable, sinc vs. gauss (%g)", err;

  res31 = sum(y*S1(x));
  res32 = sum(S1(y,1)*x);
  err3 = reldif(res31, res32);
  report, err3 < 1e-10, "direct/adjoint inner product test for nonseparable + sinc (%g)", err3;

  res41 = sum(y*S2(x));
  res42 = sum(S2(y,1)*x);
  err4 = reldif(res41, res42);
  report, err4 < 1e-10, "direct/adjoint inner product test for nonseparable + gauss (%g)", err4;
}

if (0n) benchmark_build, main, pixelsize=pixelsize;
if (0n) test_xform, main, pixelsize=pixelsize;

/* Checking gradients. */
x = random(mira_image_size(main));
x /= sum(x);


func test_gradient(main, x, pixelsize=, flags=)
{
  local grd;
  write, format="%s\n", "";
  inform, "Checking gradients with flags="+mira_format_flags(flags);
  mira_config, main, pixelsize=pixelsize, flags=flags;
  mira_cost_and_gradient, main, x, grd;
  checkgradient, closure("mira_cost", main), grd, x, number=1000;
}

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_VIS2;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_T3PHI|MIRA_HANIFF_APPROX;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_T3PHI|MIRA_VON_MISES_APPROX;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_T3PHI|MIRA_CONVEX_LIMIT;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_T3AMP;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_T3AMP|MIRA_FIT_T3PHI;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_T3AMP|MIRA_FIT_T3PHI|MIRA_CONVEX_APPROX;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_VIS2|MIRA_FIT_T3;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_VIS2|MIRA_FIT_T3|MIRA_CONVEX_APPROX;

if (1n) test_gradient, main, x, pixelsize=pixelsize,
                 flags=MIRA_FIT_VIS2|MIRA_FIT_T3PHI;

func phase_cost_checker(dims, flags)
{
  if (is_void(dims)) dims = 100;
  if (is_void(flags)) flags = MIRA_VON_MISES_APPROX;
  obj = h_new(dat = random_n(dims)*pi,
              wgt = 1 + random(dims),
              mdl = random_n(dims)*pi,
              flags = flags);
  h_evaluator, obj, "_eval_phase_cost";
  return obj;
}

func _eval_phase_cost(this, x, &g)
{
  return _mira_phase_cost(this.flags, x, this.wgt, this.dat, g);
}

func amplitude_cost_checker(dims, flags)
{
  if (is_void(dims)) dims = 100;
  if (is_void(flags)) flags = 0;
  obj = h_new(dat = max(0, random_n(dims)+3),
              wgt = 1 + random(dims),
              mdl = random_n(dims)*pi,
              flags = flags);
  h_evaluator, obj, "_eval_amplitude_cost";
  return obj;
}

func _eval_amplitude_cost(this, x, &g)
{
  return _mira_amplitude_cost(this.flags, x, this.wgt, this.dat, g);
}

if (0n) {
  f1 = phase_cost_checker(100, MIRA_VON_MISES_APPROX);
  x1 = f1.mdl; g1 = 1n; f1,x1,g1; checkgradient,f1,g1,x1;

  f2 = phase_cost_checker(100, MIRA_HANIFF_APPROX);
  x2 = f2.mdl; g2 = 1n; f2,x2,g2; checkgradient,f2,g2,x2;

  f3 = phase_cost_checker(100, MIRA_CONVEX_LIMIT);
  x3 = f3.mdl; g3 = 1n; f3,x3,g3; checkgradient,f3,g3,x3;

  f4 = amplitude_cost_checker(100);
  x4 = f4.mdl; g4 = 1n; f4,x4,g4; checkgradient,f4,g4,x4;
}
