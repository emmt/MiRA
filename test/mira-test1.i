/*
 * mira-test1.i -
 *
 * Example of image reconstruction session with MiRA (Multi-aperture Image
 * Reconstruction Algorithm) using the first dataset from the "The 2004
 * Optical/IR Interferometry Imaging Beauty Contest" (Lawson et al., 2004).
 */

if (is_void(MIRA_HOME)) {
  /* Load MiRA software. */
  include, "../src/mira.i";
}

/* Load OI-FITS data file ('db1' will be our MiRA instance for this
   data file; if there are several spectral channels in the data file,
   you must choose one with keyword EFF_WAVE or choose a spectral
   range with keywords EFF_WAVE and EFF_BAND): */
db1 = mira_new(MIRA_HOME+"../data/data1.oifits");

/* Configure data instance for image reconstruction parameters (DIM is
   the number of pixels along the width and height of the restored
   image; FOV is the size of the corresponding field of view in
   radians; XFORM is the name of the method to approximate the Fourier
   transform, can be "exact" or "fft", default is "exact"): */
mira_config, db1, dim=100, pixelsize=0.5*MIRA_MILLIARCSECOND,
  xform="exact" /*"fft"*/;

/* Choose a suitable regularization method: */
rgl = rgl_new("smoothness");

/* Attempt an image reconstruction (from scratch): */
dim = mira_get_dim(db1);
img0 = array(double, dim, dim);
img0(dim/2, dim/2) = 1.0;
img1 = mira_solve(db1, img0, maxeval=500, verb=1, xmin=0.0, normalization=1,
                  regul=rgl, mu=1e6);

/* Continue reconstruction with recentered image: */
img1 = mira_solve(db1, mira_recenter(img1), maxeval=500, verb=1, xmin=0.0,
                  normalization=1, regul=rgl, mu=1e6);
img1 = mira_solve(db1, mira_recenter(img1), maxeval=500, verb=1, xmin=0.0,
                  normalization=1, regul=rgl, mu=1e6);
img1 = mira_solve(db1, mira_recenter(img1), maxeval=500, verb=1, xmin=0.0,
                  normalization=1, regul=rgl, mu=1e6);


/* Extract central part of image and restart with higher resolution
   image and smaller field of view: */
dim = mira_get_dim(db1);
cut = (dim + 2)/4;
scale = 4.0;
new_img1 = mira_rescale(img1(1+cut:-cut, 1+cut:-cut), scale=scale);
new_dim = dimsof(new_img1)(2);
new_pixelsize = mira_get_pixelsize(db1)/scale;
mira_config, db1, dim=new_dim, pixelsize=new_pixelsize;
rgl = rgl_new("smoothness");
new_img1 = mira_solve(db1, new_img1, maxeval=500, verb=1, xmin=0.0,
                      normalization=1, regul=rgl, mu=1e6);

/* Choose a l2-l1 smoothness. */
rgl = rgl_new("xsmooth", "cost","cost_l2l1", "threshold",2e-5,
              "dimlist",dimsof(new_img1));

#if 0
/* Smooth support. */
mira_config, db1, dim=256, pixelsize=0.1*MIRA_MILLIARCSECOND, xform="fft";
dim = mira_get_dim(db1);
r = abs(mira_get_x(db1), mira_get_x(db1)(-,));
prior = 1.0/(1.0 + (2.0*r/(5.0*MIRA_MILLIARCSECOND))^2);
prior *= 1.0/sum(prior);
rgl_config, (rgl = rgl_new("quadratic")), "W", linop_new("diagonal", 1.0/prior);
img0 = (prior == max(prior)); img0 *= 1.0/sum(img0);
img1 = mira_solve(db1, img0, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e4);
img1 = mira_solve(db1, img1, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e4);
img2 = mira_solve(db1, img1, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=5e3);
img2 = mira_solve(db1, img2, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=5e3);
img3 = mira_solve(db1, img2, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=2e3);
img3 = mira_solve(db1, img3, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=2e3);
img4 = mira_solve(db1, img3, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e3);
img4 = mira_solve(db1, img4, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=1e3);
#endif

#if 0
/* maximum entropy. */
SOFT = "/home/eric/work/mira-pub/mira-0.9/";
include, "histo.i";
include, SOFT+"mira.i";
include, "~/devel/nelder-mead.i";
include, SOFT+"fit_limb.i";
include, SOFT+"azimuthal.i";
DIM = 256;
PIXELSIZE = 0.1*MIRA_MILLIARCSECOND;

mira_config, db1, dim=DIM, pixelsize=PIXELSIZE, xform="fft";

// fir a 2-D Gaussian as the prior:
gauss = mira_fit_gaussian_2d(db1, 3*MIRA_MILLIARCSECOND, 5*MIRA_MILLIARCSECOND, 150*(MIRA_PI/180));
gauss = mira_fit_gaussian_2d(db1, gauss.xbest);
prior = gauss.model; prior *= 1.0/sum(prior);
rgl_config, (ent = rgl_new("entropy")), "type","log", "normalized",0n, "prior",prior;

img_ent = mira_solve(db1, prior, maxeval=500, verb=1, xmin=1e-8, normalization=1, regul=ent, mu=4e3, view=2, mem=3);

/* Tikhonov. */
img_l2 = mira_solve(db1, prior, maxeval=500, verb=1, xmin=0.0, normalization=1, regul=rgl, mu=5e7, view=2);

#endif
