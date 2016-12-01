.SUFFIXES: .i
.PHONY: clean default distclean distrib update install

#------------------------------------------------------------------------------
#

SOURCES = fft_utils.i \
          ipy.i \
          fmin.i \
          img.i \
          linop.i \
          mira.i \
          oifits.i \
          options.i \
          plot.i \
          rgl.i \
          utils.i \
          mira mira-batch.i \
          mira-demo.i mira-test1.i mira-test2.i

IPY_DIR = /home/eric/work/code/yorick/ipy
YOR_DIR = /home/eric/yorick
YLIB_DIR = /home/eric/git/ylib
YOIFITS_DIR = /home/eric/git/YOIFITS

default:
	@echo "There is no default target.  Try one of:"
	@echo "    make update"
	@echo "    make clean"
	@echo "    make distclean"
	@echo "    make distrib [VERSION=###]"
	@echo "    make test"
#	@echo "    make install"

update: $(SOURCES)

clean:
	rm -f core *~


distclean: clean
	rm -rf old

TEST_FLAGS=-pixelsize=0.2 -dim=100 --regul=totvar -regul_isotropic -ftol=0 -gtol=0 -maxeval=2000 --overwrite -normalization=1.0 -xmin=0.0 --regul_mu=1E4
test:
	./MiRA $(TEST_FLAGS)                     data/data1.oifits test1.fits
	./MiRA $(TEST_FLAGS) -initial=test1.fits data/data1.oifits test2.fits

fft_utils.i: $(YLIB_DIR)/fft_utils.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

#fits.i: $(YOR_DIR)/fits.i
#	rm -f $@.bak
#	test -f $@ -o -h $@ && mv $@ $@.bak || true
# 	cp -a $< $@

ipy.i: $(IPY_DIR)/ipy.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

rgl.i: $(IPY_DIR)/rgl.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

linop.i: $(IPY_DIR)/linop.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

fmin.i: $(YLIB_DIR)/fmin.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

img.i: $(YLIB_DIR)/img.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@


options.i: $(YOR_DIR)/options.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

optimpack.i: $(YOR_DIR)/optimpack.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

plot.i: $(YLIB_DIR)/plot.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

#oifits.i: /home/eric/work/mira/oifits.i
#	rm -f $@.bak
#	test -f $@ -o -h $@ && mv $@ $@.bak || true
#	cp -a $< $@

utils.i: $(YLIB_DIR)/utils.i
	rm -f $@.bak
	test -f $@ -o -h $@ && mv $@ $@.bak || true
	cp -a $< $@

DISTRIB_FILES = $(SOURCES) Makefile AUTHOR COPYING README NEWS
DISTRIB_SUBDIRS = data

distrib:
	@if test "x$(VERSION)" = "x"; then \
	  version=`grep '^MIRA_VERSION *= *"' mira.i | sed 's/^.*= *"\([^"]*\).*/\1/'`; \
	  if test "x$$version" = "x"; then \
	    echo >&2 "error: bad MIRA_VERSION in file \"mira.i\""; \
	    return 1; \
	  fi; \
	else \
	  version=$(VERSION); \
	fi; \
	pkgdir=mira-$$version; \
	archive=$$pkgdir.tar.bz2; \
	if test -e "$$pkgdir"; then \
	  echo >&2 "error: $$pkgdir already exists"; \
	  return 1; \
	fi; \
	if test -e "$$archive"; then \
	  echo >&2 "error: $$archive already exists"; \
	  return 1; \
	fi; \
	mkdir "$$pkgdir"; \
	for file in $(DISTRIB_FILES); do \
	  cp -p "$$file" "$$pkgdir/."; \
	done; \
	for dir in $(DISTRIB_SUBDIRS); do \
	  cp -a "$$dir" "$$pkgdir/$$dir"; \
	  (cd "$$pkgdir/$$dir/."; rm -f *~); \
	done; \
	if test "x$(VERSION)" != "x"; then \
	  rm -f "$$pkgdir/mira.i"; \
	  sed <mira.i >"$$pkgdir/mira.i" -e 's/^MIRA_VERSION *= *".*/MIRA_VERSIOn = "$(VERSION)";/'; \
	  touch -r mira.i "$$pkgdir/mira.i"; \
	fi; \
	echo "$$version" > "$$pkgdir/VERSION"; \
	tar cf - "$$pkgdir" | bzip2 -9 > "$$archive"; \
	rm -rf "$$pkgdir"; \
	echo "archive $$archive created"; \
	return 0

