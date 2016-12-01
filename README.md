# MIRA: Multi-aperture Image Reconstruction Algorithm

MiRA (Multi-aperture Image Reconstruction Algorithm) is an algorithm for
image reconstruction from interferometric data.


## Quick instructions for using MiRA

You must have [Yorick](http://dhmunro.github.io/yorick-doc/) (at least
version 2.1), [Yeti](https://github.com/emmt/Yeti) (at least version 6.2.0) and
[OptimPack](https://cral.univ-lyon1.fr/labo/perso/eric.thiebaut/?Software/OptimPack)
(version 1.2) installed.  For faster operations, you may want to use
[YNFFT](https://github.com/emmt/ynfft), a Yorick plugin for the *Nonequispaced
Fast Fourier Transform*.


Launch Yorick and load `"mira.i"` (this should also load Yeti plugin):

    include, "mira.i";

Load OI-FITS data file (`db` will be your MiRA instance for this data file; in
the data file, you may select a spectral range with keywords `eff_wave` and
`eff_band` or with keywords `wavemin` and `wavemax`):

    db = mira_new("data/data1.oifits");

Configure data instance for image reconstruction parameters (keyword `dim` is
the number of pixels along the width and height of the restored image; keyword
`pixelsize` is the size of the pixel in radians; keyword `xform` is the name of
the method to approximate the Fourier transform, can be `"exact"`, `"fft"` or
`"nfft"`, default is `"exact"`):

    mira_config, db, dim=50, pixelsize=0.5*MIRA_MILLIARCSECOND,
            xform="nfft";

Choose a suitable regularization method:

    rgl = rgl_new("smoothness");

Attempt an image reconstruction (from scratch):

    img1 = mira_solve(db, maxeval=500, verb=1, xmin=0.0, normalization=1,
                      regul=rgl, mu=1e6);

where `img1` is the output image, `db` is the data instance, `maxeval` is the
maximum number of evaluations of the cost function, `verb` is set to one to
display convergence information at every successful step (`verb=10` to display
this information every 10 steps and `verb=0` or nil to display no information),
`xmin=0` to enforce a positivity constraint (`xmin` is the minimum allowed
value in the restored image, you certainly do not to want to omit this option),
`regul` set the regularization method, `mu` is the global weight of
regularization (the higher its value the smoother the result), `ftol` controls
the stopping criterion of OptimPack1 (which to see).

You can also try a reconstruction given an initial image `img0`:

    img2 = mira_solve(db, img0, maxeval=500, verb=1, xmin=0.0,
                        normalization=1, regul=rgl, mu=1e6);

Note that if `img0` is not of size `dim×dim` it will be resampled to that size
(*i.e.*, assuming the field of view is the same).

With keyword `zap_data`, you can just build the default image as imposed by
the regularization:

    img0 = mira_solve(db, maxeval=500, verb=1, xmin=0.0, normalization=1,
                      regul=rgl, mu=1e6, zap_data=1);


## Caveats

* *All* units are in SI (*i.e.*, angles are in radians, wavelengths in meters,
  *etc.*).  A very common error in parameter settings is to use completely out
  of range values because you assume the wrong units.  To make things a little
  bit easier, some constants are pre-defined by MiRA package, and you can use
  for instance: `3*MIRA_MILLIARCSECOND` (instead of 1.45444e-08 radians) or
  `5*MIRA_MICRON` (instead of 5e-6 meters).

* Positivity is a *must*, you cannot expect a good image reconstruction (at
  least with any practical u-v coverage) without option `xmin=0` in
  `mira_solve`.  If you use certain regularizations such as the entropy, you
  must specify a strictly positive value for `xmin` (choose a very small value,
  for instance: `1e-50*nrm/npix` where `nrm` is the normalization level and
  `npix` the total number of pixels).

* Likewise, `normalization=1` must not be omitted if your data set obeys
  OI-FITS standard (*i.e.*, visibilities are normalized) and has no explicit
  measurement at frequency `(0,0)`.

* The cost function is highly non-quadratic and may cause difficulties for
  OptimPack to converge.  To overcome this, it is sometimes very effective to
  simply restart `mira_solve` with the an initial image provided by a previous
  reconstruction (possibly after recentering by `mira_recenter`).  For
  instance:

  ```
    img2 = mira_solve(db, img0, maxeval=500, verb=1, xmin=0.0,
                      normalization=1, regul=rgl, mu=1e6);
    img2 = mira_recenter(img2);
    img2 = mira_solve(db, img2, maxeval=500, verb=1, xmin=0.0,
                      normalization=1, regul=rgl, mu=1e6);
    img2 = mira_recenter(img2);
    img2 = mira_solve(db, img2, maxeval=500, verb=1, xmin=0.0,
                      normalization=1, regul=rgl, mu=1e6);
    ...
  ```

* It is usually better to work with a purposely too high regularization level
  and then lower the value of `mu` as the image reconstruction converges.


# References

* Thiébaut, É.: *"MiRA: an effective imaging algorithm for optical
  interferometry"* in SPIE Proc. Astronomical Telescopes and Instrumentation
  **7013**, 70131I-1-70131I-12 (2008)
  [DOI](http://dx.doi.org/10.1117/12.788822).

* Thiébaut, É. & Giovannelli, J.-F.: *"Image Reconstruction in Optical
  Interferometry"* in IEEE Signal Processing Magazine **27**, pp. 97-109 (2010)
  [DOI](http://dx.doi.org/10.1109/MSP.2009.934870).

* Thiébaut, É.: *"Image reconstruction with optical interferometers"* in New
  Astronomy Reviews **53**, pp. 312-328 (2009)
  [DOI](http://dx.doi.org/10.1016/j.newar.2010.07.011).
