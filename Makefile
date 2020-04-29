.SUFFIXES: .i
.PHONY: clean all default distclean test distrib install help
srcdir = .

# DESTDIR can be overriden by package managers.
DESTDIR=

#------------------------------------------------------------------------------

MAN_PAGES = ymira.1
BIN_FILES = ymira
MIRA_FILES = \
    mira2.i \
    mira2_batch.i \
    mira2_config.i \
    mira2_cost.i \
    mira2_data.i \
    mira2_dirty.i \
    mira2_image.i \
    mira2_plugin_central_star.i \
    mira2_solver.i \
    mira2_tests.i \
    mira2_utils.i \
    mira2_xform.i
TEST_FILES = \
    mira-demo.i \
    mira-test1.i \
    mira-test2.i
DATA_FILES = data1.oifits data2.oifits README
OTHER_FILES = AUTHOR LICENSE Makefile configure \
    README.md INSTALL.md USAGE.md NEWS.md
DOC_FILES = $(MAN_PAGES)

MIRA_SRC = $(srcdir)/src

default: all

all:
	@echo "Nothing to do, execute 'make install' to install"

help:
	@echo "Try one of:"
	@echo " - 'make' or 'make all' to build"
	@echo " - 'make clean' to cleanup source tree but keep configuration"
	@echo " - 'make distclean' to cleanup source tree"
	@echo " - 'make distrib [VERSION=...]' to build archive"
	@echo " - 'make install [INCDIR=...] [YORICK=...] [BINDIR=...]' to install"
	@echo " - 'make test' to run some tests"

clean:
	rm -f core *~ $(srcdir)/src/*~

distclean: clean
	rm -f $(CONFIG)

TEST_FLAGS=-pixelsize=0.1mas -fov=20mas -regul=hyperbolic -mu=3e3 -tau=5e-5 -ftol=0 -gtol=0 -maxeval=1000 -overwrite -save_visibilities -save_initial -flux=1 -min=0 -verb=10
test:
	$(srcdir)/bin/ymira $(TEST_FLAGS) -initial=Dirac -bootstrap=1 -recenter \
	    $(srcdir)/data/data1.oifits test1.fits
	$(srcdir)/bin/ymira $(TEST_FLAGS) -initial=test1.fits -recenter \
	    $(srcdir)/data/data1.oifits test2.fits

test-plugin:
	$(srcdir)/bin/ymira -plugin=example $(TEST_FLAGS) -example_option=42 -initial=Dirac \
	    $(srcdir)/data/data1.oifits test3.fits

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
	            -e "s,^\(YORICK=.{MIRA_YORICK:\).*,\1-$$YORICK}," \
	            -e "s,^\(SRCDIR=.{MIRA_SRCDIR:\).*,\1-$$INCDIR},"; \
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
