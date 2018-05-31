.SUFFIXES: .i
.PHONY: clean default distclean test distrib install
srcdir = .

# DESTDIR can be overriden by package managers.
DESTDIR=

#------------------------------------------------------------------------------
#

MAN_PAGES = ymira.1
BIN_FILES = ymira
MIRA_FILES = mira2.i mira2_batch.i mira2_config.i mira2_cost.i mira2_data.i \
    mira2_image.i mira2_solver.i mira2_tests.i mira2_utils.i mira2_xform.i
#IPY_FILES = linop.i rgl.i
#YLIB_FILES = options.i utils.i
#YOIFITS_FILES = oifits.i
TEST_FILES = mira-demo.i mira-test1.i mira-test2.i
DATA_FILES = data1.oifits data2.oifits README
OTHER_FILES = AUTHOR LICENSE Makefile configure \
    README.md INSTALL.md USAGE.md NEWS.md
DOC_FILES = $(MAN_PAGES)

MIRA_SRC = $(srcdir)/src
IPY_SRC = $(srcdir)/lib/ipy
YLIB_SRC = $(srcdir)/lib/ylib
YOIFITS_SRC = $(srcdir)/lib/yoifits

default:
	@echo "There is no default target.  Try one of:"
	@echo "    make clean"
	@echo "    make distclean"
	@echo "    make distrib [VERSION=...]"
	@echo "    make install [INCDIR=...] [YORICK=...] [BINDIR=...]"
	@echo "    make test"

clean:
	rm -f core *~ $(srcdir)/src/*~

distclean: clean
	rm -f $(CONFIG)

TEST_FLAGS=-pixelsize=0.1mas -fov=20mas -regul=hyperbolic -bootstrap=1 -recenter -mu=3e3 -tau=5e-5 -ftol=0 -gtol=0 -maxeval=2000 -overwrite -save_visibilities -save_initial -normalization=1 -min=0 -verb=10
test:
	$(srcdir)/bin/ymira $(TEST_FLAGS) -initial=Dirac \
	    $(srcdir)/data/data1.oifits test1.fits
	$(srcdir)/bin/ymira $(TEST_FLAGS) -initial=test1.fits \
	    $(srcdir)/data/data1.oifits test2.fits

# Installation parameters are variables so that they can be overwritten when
# calling make (the installation script copies the make variables to avoid
# interpreting the make variables more than necessary).
CONFIG=install.cfg
YORICK=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *YORICK *=/!d;s/^[^=]*= *//;s/ *$$//'`
BINDIR=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *BINDIR *=/!d;s/^[^=]*= *//;s/ *$$//'`
INCDIR=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *INCDIR *=/!d;s/^[^=]*= *//;s/ *$$//'`
MANDIR=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *MANDIR *=/!d;s/^[^=]*= *//;s/ *$$//'`

$(CONFIG):
	@echo "You must run the configuration script first"
	@echo "Try \"configure --help\" for options."

install:
	@ INCDIR=$(INCDIR); \
	  BINDIR=$(BINDIR); \
	  MANDIR=$(MANDIR); \
	  YORICK=$(YORICK); \
	  if test "x$$INCDIR" = "x" -o \( "x$$BINDIR" != "x" -a "x$$YORICK" = "x" \); then \
	      echo >&2 "Run the configuration script \"configure\" first, or specify"; \
	      echo >&2 "installation parameters on the command line:"; \
	      echo >&2 ""; \
	      echo >&2 "    make install YORICK=... BINDIR=... INCDIR=... MANDIR=... DESTDIR=..."; \
	      echo >&2 ""; \
	      return 1; \
	  fi; \
	  echo "Installation settings:"; \
	  echo "  BINDIR = $$BINDIR"; \
	  echo "  MANDIR = $$MANDIR"; \
	  echo "  INCDIR = $$INCDIR"; \
	  echo "  YORICK = $$YORICK"; \
	  if test "x$$INCDIR" != "x"; then \
	      if test "x$(DESTDIR)" != "x"; then \
	          INCDIR=$(DESTDIR)/$$INCDIR; \
	      fi; \
	      mkdir -p "$$INCDIR"; \
	      for file in $(MIRA_FILES); do \
	          echo "Installing $$file in $$INCDIR"; \
	          dst=$$INCDIR/$$file; \
	          src=$(srcdir)/src/$$file; \
	          if ! test -f "$$src"; then \
	              echo >&2 "Missing file \"$$file\""; \
	              return 1; \
	          fi; \
	          cp -p "$$src" "$$dst"; \
	          chmod 644 "$$dst"; \
	      done; \
	      for file in $(IPY_FILES); do \
	          echo "Installing $$file in $$INCDIR"; \
	          dst=$$INCDIR/$$file; \
	          src=$(IPY_SRC)/$$file; \
	          if ! test -f "$$src"; then \
	              src=$(srcdir)/src/$$file; \
	              if ! test -f "$$src"; then \
	                  echo >&2 "Missing file \"$$file\""; \
	                  return 1; \
	              fi; \
	          fi; \
	          cp -p "$$src" "$$dst"; \
	          chmod 644 "$$dst"; \
	      done; \
	      for file in $(YLIB_FILES); do \
	          echo "Installing $$file in $$INCDIR"; \
	          dst=$$INCDIR/$$file; \
	          src=$(YLIB_SRC)/$$file; \
	          if ! test -f "$$src"; then \
	              src=$(srcdir)/src/$$file; \
	              if ! test -f "$$src"; then \
	                  echo >&2 "Missing file \"$$file\""; \
	                  return 1; \
	              fi; \
	          fi; \
	          cp -p "$$src" "$$dst"; \
	          chmod 644 "$$dst"; \
	      done; \
	      for file in $(YOIFITS_FILES); do \
	          echo "Installing $$file in $$INCDIR"; \
	          dst=$$INCDIR/$$file; \
	          src=$(YOIFITS_SRC)/$$file; \
	          if ! test -f "$$src"; then \
	              src=$(srcdir)/src/$$file; \
	              if ! test -f "$$src"; then \
	                  echo >&2 "Missing file \"$$file\""; \
	                  return 1; \
	              fi; \
	          fi; \
	          cp -p "$$src" "$$dst"; \
		  chmod 644 "$$dst"; \
	      done; \
	  fi; \
	  if test "x$$BINDIR" != "x"; then \
	      if test "x$(DESTDIR)" != "x"; then \
	          BINDIR=$(DESTDIR)/$$BINDIR; \
	      fi; \
	      mkdir -p "$$BINDIR"; \
	      for file in $(BIN_FILES); do \
	          echo "Installing $$file in $$BINDIR"; \
	          dst=$$BINDIR/$$file; \
	          src=$(srcdir)/bin/$$file; \
	          if ! test -f "$$src"; then \
	              echo >&2 "Missing file \"$$file\""; \
	              return 1; \
	          fi; \
	          sed <"$$src" >"$$dst" \
	            -e "s,^YORICK=.*,YORICK=$$YORICK," \
	            -e "s,^INCDIR=.*,INCDIR=$$INCDIR,"; \
	          chmod 755 "$$dst"; \
	      done; \
	      if test "x$$MANDIR" != "x"; then \
	          if test "x$(DESTDIR)" != "x"; then \
	              MANDIR=$(DESTDIR)/$$MANDIR; \
	          fi; \
	          mkdir -p "$$MANDIR/man1"; \
	          for file in $(MAN_PAGES); do \
	              echo "Installing $$file.gz in $$MANDIR/man1"; \
	              dst=$$MANDIR/man1/$$file.gz; \
	              src=$(srcdir)/doc/$$file; \
	              if ! test -f "$$src"; then \
	                  echo >&2 "Missing file \"$$file\""; \
	                  return 1; \
	              fi; \
	              gzip <"$$src" >"$$dst"; \
	              chmod 644 "$$dst"; \
	          done; \
	      fi; \
	  fi

distrib:
	@if test "x$(VERSION)" = "x"; then \
	  version=`grep '^MIRA_VERSION *= *"' "$(MIRA_SRC)/mira.i" | sed 's/^.*= *"\([^"]*\).*/\1/'`; \
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
	dstdir=$$pkgdir/src; \
	mkdir -p "$$dstdir"; \
	for file in $(MIRA_FILES); do \
	  if test "$$file" != "mira.i" -o "x$(VERSION)" != "x"; then \
	    cp -p "$(srcdir)/src/$$file" "$$dstdir/."; \
	  else \
	    sed <"$(srcdir)/src/$$file" >"$$dstdir/$$file" \
	      -e 's/^MIRA_VERSION *= *".*/MIRA_VERSIOn = "$(VERSION)";/'; \
	    touch -r "$(srcdir)/src/$$file" "$$dstdir/$$file"; \
	  fi; \
	done; \
	for file in $(IPY_FILES); do \
	  cp -p "$(IPY_SRC)/$$file" "$$dstdir/."; \
	done; \
	for file in $(YLIB_FILES); do \
	  cp -p "$(YLIB_SRC)/$$file" "$$dstdir/."; \
	done; \
	for file in $(YOIFITS_FILES); do \
	  cp -p "$(YOIFITS_SRC)/$$file" "$$dstdir/."; \
	done; \
	dstdir=$$pkgdir; \
	mkdir -p "$$dstdir"; \
	for file in $(OTHER_FILES); do \
	  cp -p "$(srcdir)/$$file" "$$dstdir/."; \
	done; \
	dstdir=$$pkgdir/data; \
	mkdir -p "$$dstdir"; \
	for file in $(DATA_FILES); do \
	  cp -p "$(srcdir)/data/$$file" "$$dstdir/."; \
	done; \
	dstdir=$$pkgdir/doc; \
	mkdir -p "$$dstdir"; \
	for file in $(DOC_FILES); do \
	  cp -p "$(srcdir)/doc/$$file" "$$dstdir/."; \
	done; \
	dstdir=$$pkgdir/test; \
	mkdir -p "$$dstdir"; \
	for file in $(TEST_FILES); do \
	  cp -p "$(srcdir)/test/$$file" "$$dstdir/."; \
	done; \
	dstdir=$$pkgdir/bin; \
	mkdir -p "$$dstdir"; \
	for file in $(BIN_FILES); do \
	  cp -p "$(srcdir)/bin/$$file" "$$dstdir/."; \
	done; \
	echo "$$version" > "$$pkgdir/VERSION"; \
	tar cf - "$$pkgdir" | bzip2 -9 > "$$archive"; \
	rm -rf "$$pkgdir"; \
	echo "archive $$archive created"; \
	return 0
