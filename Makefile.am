ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext include bin analyses test doc

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-yoda=$(YODAPATH) --with-hepmc=$(HEPMCPATH) --with-fastjet=$(FASTJETPATH)

pkgconfigdir = $(libdir)/pkgconfig
dist_pkgconfig_DATA = rivet.pc

if CODE_COVERAGE_ENABLED

REPORTDIR=$(top_srcdir)/test/coverage

REPORT=$(REPORTDIR)/coverage.info

COVERAGE_EXCL = '*/eigen3/*' '*/HepMC3/*' '*/YODA/*' '/usr/*' '*/test/*'

LCOV_FLAGS = --quiet

if WITH_OSX

LCOV_FLAGS += --ignore-errors inconsistent,inconsistent,usage

endif

cov-report:
	@lcov --compat-libtool --capture --directory $(top_srcdir) --output-file $(REPORT) $(LCOV_FLAGS)
	@lcov --remove $(REPORT) --output-file $(REPORT) $(LCOV_FLAGS) $(COVERAGE_EXCL)
	@genhtml -o $(REPORTDIR)/html $(REPORT) $(LCOV_FLAGS)

cov-reset:
	@rm -fr $(REPORTDIR)/html $(REPORT)
	@find . \( -name "*.gcda" -o -name "*.gcno" \) -exec rm {} \;

check:
	@make cov-report

clean-local:
	@make cov-reset
	@rm -rf a.out

else

clean-local:
	@rm -rf a.out

endif

doc:
	cd doc && $(MAKE) doc

.PHONY : doc dox pyclean analyses-indiv

if ENABLE_ANALYSES
if ENABLE_PYEXT

install-data-local:
	$(MAKE) -C doc dat json
	$(mkdir_p) $(DESTDIR)$(pkgdatadir)
	$(install_sh_DATA) doc/anaindex/analyses.dat $(DESTDIR)$(pkgdatadir)/ || true
	$(install_sh_DATA) doc/anaindex/analyses.json $(DESTDIR)$(pkgdatadir)/ || true

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
