/*
 * mira-demo.i -
 *
 * Demonstration of image reconstruction session with MiRA (Multi-aperture
 * Image Reconstruction Algorithm) using the first dataset from the "The 2004
 * Optical/IR Interferometry Imaging Beauty Contest" (Lawson et al., 2004).
 */

if (is_void(MIRA_HOME)) {
  /* Load MiRA software. */
  include, "../src/mira.i";
}

/* Load OI-FITS data file ('mh1' will be our MIRA instance for this
   data file; if there are several spectral channels in the data file,
   you must choose one with keyword EFF_WAVE or choose a spectral
   range with keywords EFF_WAVE and EFF_BAND): */
mh1 = mira_new(MIRA_HOME+"../data/data1.oifits");

/* Configure data instance for image reconstruction parameters (DIM is
   the number of pixels along the width and height of the restored
   image; FOV is the size of the corresponding field of view in
   radians; XFORM is the name of the method to approximate the Fourier
   transform, can be "exact" or "fft", default is "exact"): */
mira_config, mh1, dim=256, pixelsize=0.1*MIRA_MILLIARCSECOND, xform="fft";

/* Smooth support. */
dim = mira_get_dim(mh1);
r = abs(mira_get_ra(mh1), mira_get_dec(mh1)(-,));
prior = 1.0/(1.0 + (2.0*r/(5.0*MIRA_MILLIARCSECOND))^2);
prior *= 1.0/sum(prior);
rgl_config, (rgl = rgl_new("quadratic")), "W", linop_new("diagonal", 1.0/prior);

if (! window_exists(0)) {
  window, 0, style="work.gs", dpi=75, width=550, height=450;
} else {
  window, 0, style="work.gs";
}
limits,10,-10,-10,10;
palette, "heat.gp";

if (! window_exists(1)) {
  window, 1, style="work.gs", dpi=75, width=450, height=450;
} else {
  window, 1, style="work.gs";
}
palette, "heat.gp";

/*---------------------------------------------------------------------------*/
/* RECONSTRUCTION WITH PHASE CLOSURE */

/*  seed  rotated  final penalty
 *   0.1  yes      4.94723e3 (viusaly the best)
 *   0.2  bad      6.95202e3
 *   0.3  bad      6.96364e3
 *   0.4  yes      4.96052e3
 *   0.5  bad      6.95272e3
 *   0.6  no       7.27172e3
 *   0.7  bad      6.67580e3
 *   0.8  no       4.94379e3
 *   0.9  bad      5.62120e3
 *   1.0  yes      4.91903e3
 */
random_seed, (is_func(fftw) ? 0.2 : 0.7);
img0 = random(dim, dim);

rdline, prompt="hit [Return] to start reconstruction";
img1 = mira_solve(mh1, img0, maxeval=1, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e5, view=3, title="Reconstruction with phase closure");
rdline, prompt="hit [Return] to continue reconstruction";
img1 = mira_solve(mh1, img1, maxeval=100, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e5, view=3, title="Reconstruction with phase closure");
img1 = mira_solve(mh1, img1, maxeval=100, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e4, view=3, title="Reconstruction with phase closure");
img1 = mira_solve(mh1, img1, maxeval=200, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e3, view=3, title="Reconstruction with phase closure");
rdline, prompt="hit [Return] to start reconstruction";

/*---------------------------------------------------------------------------*/
/* RECONSTRUCTION WITHOUT PHASE INFORMATION */

/* The following results were with the old (bogus) optimizer.
 *
 * seed  rotated  final penalty
 *   0.1  yes      4.94723e3
 *   0.2  bad      6.95202e3
 *   0.3  bad      6.96364e3
 *   0.4  yes      4.96052e3
 *   0.5  bad      6.95272e3     1.173965308004727e+04
 *   0.6  no       7.27172e3     1.363081123728543e+04
 *   0.7  bad      6.67580e3     1.176493065068268e+04
 *   0.8  no       4.94379e3     1.175854673628141e+04
 *   0.9  bad      5.62120e3     1.241042996781854e+04
 *   1.0  yes      4.91903e3     1.368804722873744e+04
 *   1.1  yes                    1.368804722873744e+04
 */
random_seed, (is_func(fftw) ? 0.5 : 0.7);
ROTATE = 0;
img0 = random(dim, dim)//(::-1,::-1);

img4 = mira_solve(mh1, img0, maxeval=0, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e5, zap_phase=1, view=3, title="Reconstruction with no phase data");
rdline, prompt="hit [Return] to start reconstruction";
img4 = mira_solve(mh1, img0, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e5, zap_phase=1, view=3, title="Reconstruction with no phase data");
img4 = mira_solve(mh1, img4, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=5e4, zap_phase=1, view=3, title="Reconstruction with no phase data");
img4 = mira_solve(mh1, img4, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=2e4, zap_phase=1, view=3, title="Reconstruction with no phase data");
img4 = mira_solve(mh1, img4, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e4, zap_phase=1, view=3, title="Reconstruction with no phase data");
img4 = mira_solve(mh1, img4, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=5e3, zap_phase=1, view=3, title="Reconstruction with no phase data");
img4 = mira_solve(mh1, img4, maxeval=70, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=3e3, zap_phase=1, view=3, title="Reconstruction with no phase data");
//img4 = mira_solve(mh1, img4, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=2e3, zap_phase=1, view=3, title="Reconstruction with no phase data");
//img4 = mira_solve(mh1, img4, maxeval=50, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e3, zap_phase=1, view=3, title="Reconstruction with no phase data");
if (ROTATE) {
  rdline, prompt="hit [Return] to turn the image";
  img4 = mira_solve(mh1, img4(::-1,::-1), maxeval=70, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=3e3, zap_phase=1, view=3, title="Reconstruction with no phase data");
 }

