.SUFFIXES: .i
.PHONY: clean default distclean test distrib install
srcdir = .

# DESTDIR can be overriden by package managers.
DESTDIR=

#------------------------------------------------------------------------------
#

BIN_FILES = mira
MIRA_FILES = mira.i mira-batch.i
IPY_FILES = linop.i rgl.i
YLIB_FILES = options.i xplot.i xplot0.i
YOIFITS_FILES = oifits.i
TEST_FILES = mira-demo.i mira-test1.i mira-test2.i
DATA_FILES = data1.oifits data2.oifits README
OTHER_FILES = AUTHOR LICENSE Makefile configure \
    README.md INSTALL.md USAGE.md NEWS.md

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

TEST_FLAGS=-pixelsize=0.1mas -fov=20mas -regul=hyperbolic -bootstrap=1 -recenter -mu=3e3 -tau=5e-5 -ftol=0 -gtol=0 -maxeval=2000 -overwrite -normalization=1 -min=0 -verb=10
test:
	$(srcdir)/bin/mira $(TEST_FLAGS) \
	    $(srcdir)/data/data1.oifits test1.fits
	$(srcdir)/bin/mira $(TEST_FLAGS) -view=-1 -initial=test1.fits \
	    $(srcdir)/data/data1.oifits test2.fits

# Installation parameters are variables so that they can be overwritten when
# calling make (the installation script copies the make variables to avoid
# interpreting the make variables more than necessary).
CONFIG=install.cfg
YORICK=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *YORICK *=/!d;s/^[^=]*= *//;s/ *$$//'`
BINDIR=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *BINDIR *=/!d;s/^[^=]*= *//;s/ *$$//'`
INCDIR=`test -f $(CONFIG) && sed <$(CONFIG) '/^ *INCDIR *=/!d;s/^[^=]*= *//;s/ *$$//'`

$(CONFIG):
	@echo "You must run the configuration script first"
	@echo "Try \"configure --help\" for options."

install:
	@ INCDIR=$(INCDIR); \
	  BINDIR=$(BINDIR); \
	  YORICK=$(YORICK); \
	  if test "x$$INCDIR" = "x" -o \( "x$$BINDIR" != "x" -a "x$$YORICK" = "x" \); then \
	      echo >&2 "Run the configuration script \"configure\" first, or specify"; \
	      echo >&2 "installation parameters on the command line:"; \
	      echo >&2 ""; \
	      echo >&2 "    make install YORICK=... BINDIR=... INCDIR=..."; \
	      echo >&2 ""; \
	      return 1; \
	  fi; \
	  echo "Installation settings:"; \
	  echo "  BINDIR = $$BINDIR"; \
	  echo "  INCDIR = $$INCDIR"; \
	  echo "  YORICK = $$YORICK"; \
	  if test "x$$INCDIR" != "x"; then \
	      mkdir -p "$$INCDIR"; \
	      for file in $(MIRA_FILES); do \
	          echo "Installing $$file in $$INCDIR"; \
	          dst=$(DESTDIR)$$INCDIR/$$file; \
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
	          dst=$(DESTDIR)$$INCDIR/$$file; \
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
	          dst=$(DESTDIR)$$INCDIR/$$file; \
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
	          dst=$(DESTDIR)$$INCDIR/$$file; \
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
	      mkdir -p "$$BINDIR"; \
	      for file in $(BIN_FILES); do \
	          echo "Installing $$file in $$BINDIR"; \
	          dst=$(DESTDIR)$$BINDIR/$$file; \
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

