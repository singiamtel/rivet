check_PROGRAMS = testMath testMatVec testCmp testAPI testNaN testBeams testMerging

AM_LDFLAGS = -L$(top_srcdir)/src $(YAMLCPP_LDFLAGS) -L$(YODALIBPATH)
if WITH_OSX
AM_LDFLAGS += -Wl,-rpath,${HEPMC3LIBPATH}
endif
LIBS = -lm -lYODA

TEST_LDADD = $(top_builddir)/src/libRivet.la $(HEPMCLDFLAGS) $(HEPMC3LDFLAGS) $(YODALDFLAGS) $(HEPMCLDLIBS) $(HEPMC3LDLIBS) $(YODALDLIBS) $(FASTJETLIBADD)

testMath_SOURCES = testMath.cc
testMath_LDADD = $(TEST_LDADD)
testMatVec_SOURCES = testMatVec.cc
testMatVec_LDADD = $(TEST_LDADD)
testCmp_SOURCES = testCmp.cc
testCmp_LDADD = $(TEST_LDADD)
testAPI_SOURCES = testAPI.cc
testAPI_LDADD = $(TEST_LDADD)
testNaN_SOURCES = testNaN.cc
testNaN_LDADD = $(TEST_LDADD)
testBeams_SOURCES = testBeams.cc
testBeams_LDADD = $(TEST_LDADD)
testMerging_SOURCES = testMerging.cc
testMerging_LDADD = $(TEST_LDADD)

TESTS_ENVIRONMENT = \
  RIVET_ANALYSIS_PATH=$(top_builddir)/analyses \
  RIVET_DATA_PATH=$(top_builddir)/analyses/data \
  RIVET_INFO_PATH=$(top_builddir)/analyses/data \
  RIVET_REF_PATH=$(top_builddir)/analyses/data \
  LIBLOCATION=$(top_builddir)/src/.libs \
  LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(HEPMCLIBPATH):$(HEPMC3LIBPATH):$(FASTJETLIBPATH):$(FJCONTRIBLIBPATH):$(YODALIBPATH) \
  DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(HEPMCLIBPATH):$(HEPMC3LIBPATH):$(FASTJETLIBPATH):$(FJCONTRIBLIBPATH):$(YODALIBPATH) \
  PYTHONPATH=$(top_builddir)/pyext/build:$(YODA_PYTHONPATH):$(PYTHONPATH) \
  PATH=$(top_builddir)/bin:$(top_srcdir)/bin:$(YODABINPATH):$(PATH) \
  PYTHON=$(PYTHON) \
  RIVET_TESTS_SRC=$(srcdir)

TESTS = \
testMath testMatVec testCmp testAPI.sh testNaN.sh testBeams testMerging.sh

if ENABLE_PYEXT
TESTS += testImport.sh

if ENABLE_ANALYSES

TESTS += testCmdLine.sh

endif
endif

EXTRA_DIST = testAPI.hepmc testCmdLine.sh testImport.sh testAPI.sh testMerging.sh testNaN.sh

if CODE_COVERAGE_ENABLED
  AM_LIBS = $(CODE_COVERAGE_LIBS)
  AM_CFLAGS = $(CODE_COVERAGE_CFLAGS)
  AM_CPPFLAGS += $(CODE_COVERAGE_CPPFLAGS)
  AM_CXXFLAGS += $(CODE_COVERAGE_CXXFLAGS)
endif

CLEANFILES = log a.out fifo.hepmc file2.hepmc out*.yoda Rivet.yoda
