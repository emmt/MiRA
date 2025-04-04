# MIRA: Multi-aperture Image Reconstruction Algorithm

MiRA (Multi-aperture Image Reconstruction Algorithm) is an algorithm for image
reconstruction from interferometric data. The software is written in
[Yorick](https://github.com/LLNL/yorick/) and expects data file in the OIFITS format.

A quick way to install MiRA is to use [EasyYorick](https://github.com/emmt/EasyYorick)
script which takes care of downloading, compiling and installing all needed code.
Otherwise, see [*Installation*](doc/INSTALL.md) for prerequisites and instructions to
install MiRA.

The man pages of MiRA give all options. See [*Usage*](doc/USAGE.md) for a short
introduction about using MiRA.

## Links

MiRA is used in a number of projects:

* [OImaging](https://www.jmmc.fr/english/tools/data-analysis/oimaging/) provides a GUI to
  assist in image reconstruction form optical interferometric data. MiRA is one of the
  image reconstruction software available in OImaging.

* [PYRA](https://github.com/jdrevon/PYRA) is a Python script by Julien Drevon to
  automatically sample and tune the parameters of MiRA.

## References

* Thiébaut, É.: *"MiRA: an effective imaging algorithm for optical interferometry"* in
  SPIE Proc. Astronomical Telescopes and Instrumentation **7013**, 70131I-1-70131I-12
  (2008) [DOI](http://dx.doi.org/10.1117/12.788822).

* Thiébaut, É. & Giovannelli, J.-F.: *"Image Reconstruction in Optical Interferometry"* in
  IEEE Signal Processing Magazine **27**, pp. 97-109 (2010)
  [DOI](http://dx.doi.org/10.1109/MSP.2009.934870).

* Thiébaut, É.: *"Image reconstruction with optical interferometers"* in New Astronomy
  Reviews **53**, pp. 312-328 (2009) [DOI](http://dx.doi.org/10.1016/j.newar.2010.07.011).
