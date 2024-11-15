noinst_LTLIBRARIES = libRivetTools.la

libRivetTools_la_SOURCES = \
  Beams.cc \
  Correlators.cc \
  Cuts.cc \
  JetUtils.cc \
  Random.cc \
  Logging.cc \
  ParticleUtils.cc \
  ParticleName.cc \
  Percentile.cc \
  ProjectionTreePlotter.cc \
  RivetHepMC_3.cc \
  RivetMT2.cc \
  RivetPaths.cc \
  RivetYODA.cc \
  Utils.cc \
  Vector4.cc \
  binreloc.c

dist_noinst_HEADERS = binreloc.h lester_mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -DENABLE_BINRELOC \
  -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"

if ENABLE_H5
  libRivetTools_la_SOURCES += RivetHDF5.cc
endif
