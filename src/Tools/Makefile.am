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
  RivetHDF5.cc \
  RivetMT2.cc \
  RivetPaths.cc \
  RivetYODA.cc \
  TypeRegistry.cc \
  Utils.cc \
  Vector4.cc \
  binreloc.c

if ENABLE_HEPMC_3
libRivetTools_la_SOURCES += RivetHepMC_3.cc
else
libRivetTools_la_SOURCES += RivetHepMC_2.cc
endif

dist_noinst_HEADERS = binreloc.h lester_mt2_bisect.hh

libRivetTools_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -DENABLE_BINRELOC \
  -DDEFAULTDATADIR=\"$(datadir)\" -DDEFAULTLIBDIR=\"$(libdir)\"

libRivetTools_la_CPPFLAGS += $(HDF5_CPPFLAGS)
libRivetTools_la_LDFLAGS = $(HDF5_LDFLAGS)
libRivetTools_la_LIBADD = $(HDF5_LIBS)
