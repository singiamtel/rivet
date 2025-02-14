if WITH_YAML_CPP
SUBDIRS = .
else
SUBDIRS = yamlcpp .
endif


EXTRA_DIST = zstr

noinst_LTLIBRARIES  = libRivetCore.la

libRivetCore_la_SOURCES = \
  Run.cc Event.cc \
  Jet.cc Particle.cc DressedLepton.cc \
  ProjectionApplier.cc Projection.cc \
  Analysis.cc AnalysisLoader.cc AnalysisInfo.cc \
  AnalysisHandler.cc ProjectionHandler.cc

if WITH_YAML_CPP
libRivetCore_la_CPPFLAGS = $(AM_CPPFLAGS)
else
libRivetCore_la_CPPFLAGS = -I$(srcdir)/yamlcpp -DYAML_NAMESPACE=RIVET_YAML $(AM_CPPFLAGS)
libRivetCore_la_LIBADD = $(builddir)/yamlcpp/librivet-yaml-cpp.la
endif
