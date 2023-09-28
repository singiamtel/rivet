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
  RivetHDF5.cc \
  RivetMT2.cc \
  RivetPaths.cc \
  RivetYODA.cc \
  TypeRegistry.cc \
  Utils.cc \
  Vector4.cc \
  binreloc.c

dist_noinst_HEADERS = binreloc.h lester_mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -DENABLE_BINRELOC \
  -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"

libRivetTools_la_CPPFLAGS += $(HDF5_CPPFLAGS)
libRivetTools_la_LDFLAGS = $(HDF5_LDFLAGS)
libRivetTools_la_LIBADD = $(HDF5_LIBS)
