SUBDIRS = rivet

if ENABLE_PYEXT

PYEXT_ENV = LDFLAGS="$(LDFLAGS)" CC="$(CXX)" CXX="$(CXX)" CXXFLAGS="$(PYEXT_CXXFLAGS)" ARCHFLAGS="" PYTHONPATH="$(DESTDIR):$(PYTHONPATH)"

PYCORE_SRC = \
    rivet/*.py \
    rivet/pdg/*.py \
    rivet/plotting/*.py

all-local: build/rivet/core.so

if WITH_OSX

build/rivet/core.so: build.py rivet/core.cpp $(PYCORE_SRC)
	$(MKDIR_P) $(abs_builddir)/build/rivet
	cd $(abs_builddir) && $(PYEXT_ENV) RIVET_LOCAL=1 $(PYTHON) $(abs_builddir)/build.py
	touch $@
	install_name_tool -change \
           $(libdir)/libRivet.dylib \
           $(abs_top_srcdir)/src/.libs/libRivet.dylib \
           $(abs_builddir)/build/rivet/core.so

# TODO: always update the .py files, without triggering an extension rebuild
install-exec-local: build/rivet/core.so
	$(MKDIR_P) $(DESTDIR)$(RIVET_PYTHONPATH)
	cp -r $(abs_builddir)/build/rivet $(DESTDIR)$(RIVET_PYTHONPATH)/
	install_name_tool -change \
         $(abs_top_srcdir)/src/.libs/libRivet.dylib \
         $(libdir)/libRivet.dylib \
         $(abs_builddir)/build/rivet/core.so

else

build/rivet/core.so: build.py rivet/core.cpp $(PYCORE_SRC)
	$(MKDIR_P) $(abs_builddir)/build/rivet
	cd $(abs_builddir) && $(PYEXT_ENV) RIVET_LOCAL=1 $(PYTHON) $(abs_builddir)/build.py
	touch $@

# TODO: always update the .py files, without triggering an extension rebuild
install-exec-local: build/rivet/core.so
	$(MKDIR_P) $(DESTDIR)$(RIVET_PYTHONPATH)
	cp -r $(abs_builddir)/build/rivet $(DESTDIR)$(RIVET_PYTHONPATH)/

endif

uninstall-local:
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/rivet

clean-local:
	rm -rf $(builddir)/build

endif
