bin_SCRIPTS = rivet-config rivet-build
dist_bin_SCRIPTS = make-plots

EXTRA_DIST =

RIVETPROGS = \
    rivet \
    rivet-mkanalysis \
    rivet-findid rivet-which \
    rivet-diffhepdata rivet-diffhepdata-all \
    rivet-cmphistos rivet-mkhtml rivet-mkhtml-tex rivet-merge rivet-mkvaldir
if ENABLE_PYEXT
dist_bin_SCRIPTS += $(RIVETPROGS)
else
EXTRA_DIST += $(RIVETPROGS)
endif

## bash completion
if ENABLE_PYEXT
dist_pkgdata_DATA = rivet-completion
bashcomp_dir = $(DESTDIR)$(sysconfdir)/bash_completion.d
install-data-local:
	$(mkdir_p) $(bashcomp_dir)
	$(install_sh_DATA) $(srcdir)/rivet-completion $(bashcomp_dir)/
uninstall-local:
	rm -f $(bashcomp_dir)/rivet-completion
else
EXTRA_DIST += rivet-completion
endif

## No-Python Rivet program
noinst_PROGRAMS = rivet-nopy
rivet_nopy_SOURCES = rivet-nopy.cc
rivet_nopy_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS) $(CPPFLAGS)
rivet_nopy_LDADD = $(top_builddir)/src/libRivet.la
rivet_nopy_LDFLAGS = -L$(HEPMC3LIBPATH)
rivet_nopy_LDADD += $(top_builddir)/src/libRivet.la -lHepMC3 -lHepMC3search
