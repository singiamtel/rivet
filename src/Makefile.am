SUBDIRS = Core Tools Projections AnalysisTools

lib_LTLIBRARIES = libRivet.la

libRivet_la_SOURCES =

libRivet_la_LDFLAGS = -export-dynamic -avoid-version -L$(YODALIBPATH)

libRivet_la_LIBADD = \
Core/libRivetCore.la \
Projections/libRivetProjections.la \
Tools/libRivetTools.la \
AnalysisTools/libRivetAnalysisTools.la \
-lYODA -ldl -lm  \
$(FASTJETLIBADD)

libRivet_la_LIBADD += -lHepMC3 -lHepMC3search
libRivet_la_LDFLAGS += -L$(HEPMC3LIBPATH)

if CODE_COVERAGE_ENABLED
  AM_LIBS = $(CODE_COVERAGE_LIBS)
  AM_CFLAGS = $(CODE_COVERAGE_CFLAGS)
  AM_CPPFLAGS += $(CODE_COVERAGE_CPPFLAGS)
  AM_CXXFLAGS += $(CODE_COVERAGE_CXXFLAGS)
endif
