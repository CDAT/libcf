#
# $Id: Makefile.am 930 2012-07-10 21:13:22Z pletzer $
#

libcfpydir = $(exec_prefix)/pycf

if HAVE_PYTHON
libCFConfig.py: $(top_builddir)/config.log
	@PYTHON@ $(top_srcdir)/pycf/generateLibCFConfig.py -b $(top_builddir) -s $(top_srcdir) -c $@

dist_libcfpy_SCRIPTS = libCFConfig.py __init__.py
endif

EXTRA_DIST = generateLibCFConfig.py

clean-local:
	-rm -rf libCFConfig.py


