# Quick instructions for using MiRA

It is assumed that MiRA is properly installed (see [`INSTALL.md`](INSTALL.md)).


## Using MiRA in Yorick interpreter

Launch Yorick and load `"mira.i"` (this should also load Yeti plugin):

    include, "mira2.i";

Load OI-FITS data file (`db` will be your MiRA instance for this data file; in
the data file, you may select a spectral range with keywords `eff_wave` and
`eff_band` or with keywords `wavemin` and `wavemax`):

    db = mira_new("data/data1.oifits");

Configure data instance for image reconstruction parameters (keyword `dim` is
the number of pixels along the width and height of the restored image; keyword
`pixelsize` is the size of the pixel in radians; keyword `xform` is the name of
the method to approximate the Fourier transform, can be `"exact"`, `"fft"` or
`"nfft"`, default is `"exact"`):

    mira_config, db, dims=50, pixelsize=0.5*MIRA_MILLIARCSECOND,
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


## Using MiRA from the command line

MiRA can be run from the command line:

    mira [OPTIONS] INPUT [...] OUTPUT

where `[OPTIONS]` are optional settings, `INPUT [...]` are any number of OIFITS
input data files and `OUTPUT` is the name of the output FITS file to save the
resulting image.  Option `-help` can be used for a short description of all
options.


### Data selection

MiRA reconstructs a gray image from the interferometric data.  The wavelengths
of the selected data can be specified as the end points of the spectral range:

    ... -wavemin=MINVAL -wavemax=MAXVAL

or as the central wavelength and bandwidth:

    ... -effwave=CENTER -effband=WIDTH

The arguments of these options have length units.  For instance:

    ... -wavemin=1.6µm -wavemax=1.8microns

is the same as:

    ... -effwave=1700nm -effband=200nm

If the input data files contain observations for more than one object, the
target to consider can be specified as:

    ... -target=NAME


### Image parameters

An initial image for the reconstruction can be provided as:

    ... -initial=FILENAME

where `FILENAME` is the name of the FITS file with the image to start with.

If no initial image is provided, the reconstruction starts with a random image
and option `-seed=NUMBER` can be used to seed the random generator.

By default, the reconstructed image will have the same pixel size and
dimensions as the initial image if provided.  Otherwise, the pixel size can be
specified with option `-pixelsize=PIXSIZ` and the image dimensions can be chosen
with `-dim=NUMBER` or `-fov=ANGLE`, where `PIXSIZ` and `ANGLE` are in angular units and `NUMBER` is a number of pixels.  For instance:

    ... -pixelsize=0.25mas -fov=100mas

or

    ... -pixelsize=250e-6arcsec -dim=400


both yield a 400×400 image with a pixel size of 0.25 milliarcsecond.


### Non-uniform Fourier transform

The nonequispaced Fourier transform of the pixels can be computed by different
methods: `-xform=exact` uses an exact transform, `-xform=nfft` use a precise
approximation by the NFFT algorithm while `-xform=fft` uses a built-in
algorithm which is less precise.  The exact transform can be very slow if the
image is large and/or if there are many data.  The two others exploit an FFT
algorithm and are faster.  If you have installed
[Yorick NFFT plug-in](https://github.com/emmt/ynfft), `-xform=nfft` is
certainly the method of choice.


### Image constraints

The total flux of the sought image, and the lower and upper bounds for pixel
values can be specified with options `-normalization=SUM`, `-min=LOWER` and
`-max=UPPER` respectively.  For instance:

    ... -normalization=1 -min=0

is a must for processing OIFITS data.


### Regularization

The regularization is the kind of prior imposed to the reconstructed image.
For all implemented regularizations, the relative strength of the prior
(compared to the data) is specified with option `-mu=µ` where `µ ≥ 0`.

The following regularizations are available:

* **Edge-preserving smoothness** is selected with the following options:

  ````
  -regul=hyperbolic -mu=µ -tau=τ -eta=η
  ````

  where `τ` is the edge threshold and `η` the scale of the finite differences
  to estimate the local gradient of the image.

  Different scales can be set for different dimensions by providing a list of
  values to `-eta`, *e.g.* `-eta=1,1,0.3`.  By default, `-eta=1` which implies
  that all finite differences are scaled the same.  If any of the scales set by
  `-eta=...` is zero, then no regularity is imposed along the corresponding
  dimensions.

  Using a very small edge threshold `τ`, compared to the norm of the local
  gradients, mimics the effects of *total variation* (TV) regularizations.
  Conversely, using a very small edge threshold yields a regularization
  comparable to *quadratic smoothness*.


* **Quadratic compactness** is selected with the following options:

  ````
  -regul=compactness -mu=µ -gamma=γ
  ````

  where `γ` is the full width at half maximum (FWHM) of the prior distribution
  of light.  This parameter has angular units.  For instance `-gamma=15mas`.


### Tuning the reconstruction

The image reconstruction amounts to minimizing a criterion which is the sum of
two terms: a data fidelity term and a regularization term.  Given these two
terms and an initial image, the algorithm proceeds by iteratively improving the
solution so as to reduce the criterion.


There are parameters to control the convergence of the algorithm and the amount
of computation.

* Options `--maxiter=NUMBER` and `--maxeval=NUMBER` control the maximum number
  of iterations and evaluations of the criterion respectively.  By default,
  there are no limits on these numbers.

* Options `--sftol=SFTOL`, `--sgtol=SGTOL` and `--sxtol=SXTOL` control the line
  search convergence.  `SFTOL` and `SGTOL` correspond to the first Wolfe
  condition (Armijo's criterion) and to the second Wolfe condition.

* Options `--ftol=FTOL` and `--gtol=GTOL` control the global convergence of the
  reconstruction.  The algorithm stops when the relative function change is
  less than `FTOL` between two successive iterations or when the norm of the
  gradient becomes smaller than `GTOL`.


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

