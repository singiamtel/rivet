ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext include bin analyses test doc

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-yoda=$(YODAPATH) --with-hepmc=$(HEPMCPATH) --with-fastjet=$(FASTJETPATH)

pkgconfigdir = $(libdir)/pkgconfig
dist_pkgconfig_DATA = rivet.pc

doc:
	cd doc && $(MAKE) doc

.PHONY : doc dox pyclean analyses-indiv

clean-local:
	@rm -rf a.out

if ENABLE_ANALYSES
if ENABLE_PYEXT

install-data-local:
	$(MAKE) -C doc dat json
	$(mkdir_p) $(DESTDIR)$(pkgdatadir)
	$(install_sh_DATA) doc/analists/analyses.dat $(DESTDIR)$(pkgdatadir)/ || true
	$(install_sh_DATA) doc/analists/analyses.json $(DESTDIR)$(pkgdatadir)/ || true

endif	
endif

## Upload tarballs to HepForge
DEST = login.hepforge.org:rivet/downloads/
upload: dist
	$(RSH) $(DIST_ARCHIVES) $(DEST)

ANASOURCES ?= $(shell cat doc/anaindex/analyses.dat)


# Check that analyses can be built individually
%.anacheck:
	@ANA=$(basename $@); \
  SRC=$(shell find analyses -type f -name "$(basename $@).cc"); \
  rivet-build Rivet$${ANA}.so $${SRC} &>/dev/null || echo "Build problem in $${ANA}"; \
  rm -rf "Rivet$${ANA}.so";

analyses-indiv: $(addsuffix .anacheck, ${ANASOURCES})
	@echo "Individual analysis checks done."
