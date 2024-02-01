SUBDIRS = rivet

if ENABLE_PYEXT

PYEXT_ENV = LDFLAGS="$(LDFLAGS)" CC="$(CXX)" CXX="$(CXX)" CXXFLAGS="$(PYEXT_CXXFLAGS)" ARCHFLAGS="" PYTHONPATH="$(DESTDIR):$(PYTHONPATH)"

.build: $(wildcard rivet/*.py)
	$(MKDIR_P) $(abs_builddir)/build/rivet
	cd $(abs_builddir) && $(PYEXT_ENV) RIVET_LOCAL=1 $(PYTHON) $(abs_builddir)/build.py
	touch $@ 

all-local: .build

install-exec-local:
	$(MKDIR_P) $(DESTDIR)$(RIVET_PYTHONPATH)
	cp -r $(abs_builddir)/build/rivet $(DESTDIR)$(RIVET_PYTHONPATH)/

uninstall-local:
	rm -rf $(DESTDIR)$(RIVET_PYTHONPATH)/rivet

clean-local:
	rm -rf $(builddir)/build .build

endif
