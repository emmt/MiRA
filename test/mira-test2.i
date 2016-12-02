/*
 * mira-test2.i -
 *
 * Example of image reconstruction session with MiRA (Multi-aperture Image
 * Reconstruction Algorithm) using the second dataset from the "The 2004
 * Optical/IR Interferometry Imaging Beauty Contest" (Lawson et al., 2004).
 */

if (is_void(MIRA_HOME)) {
  /* Load MiRA software. */
  include, "../src/mira.i", 1;
}

/* Load OI-FITS data file ('db2' will be our MiRA instance for this data
   file; if there are several spectral channels in the data file, you
   must choose one with keyword EFF_WAVE or choose a spectral range
   with keywords EFF_WAVE and EFF_BAND): */
db2 = mira_new(MIRA_HOME+"../data/data2.oifits");

/* Configure data instance for image reconstruction parameters (DIM is
   the number of pixels along the width and height of the restored
   image; FOV is the size of the corresponding field of view in
   radians; XFORM is the name of the method to approximate the Fourier
   transform, can be "exact" or "fft", default is "exact"): */
mira_config, db2, dim=100, pixelsize=0.5*MIRA_MILLIARCSECOND, xform="fft";

/* Choose a suitable regularization method: */
rgl = rgl_new("smoothness");

/* Attempt an image reconstruction (from scratch): */
dim = mira_get_dim(db2);
img0 = array(double, dim, dim);
img0(dim/2, dim/2) = 1.0;
img2 = mira_solve(db2, img0, maxeval=500, verb=1, xmin=0.0, normalization=1,
                  regul=rgl, mu=1e6);

/* Iterate until convergence. */
img2 = mira_solve(db2, img2, maxeval=500, verb=1, xmin=0.0, normalization=1,
                  regul=rgl, mu=1e6);
img2 = mira_solve(db2, img2, maxeval=500, verb=1, xmin=0.0, normalization=1,
                  regul=rgl, mu=1e6);
