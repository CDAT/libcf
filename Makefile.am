## This is a automake file, part of Unidata's NetCDF CF Library package.
# Copyright 2006-2011, see the COPYRIGHT file for more information.

# This is the main automake file for the libCF library. 

# These files get added to the distribution.
EXTRA_DIST = COPYRIGHT
AM_DISTCHECK_CONFIGURE_FLAGS = "LDFLAGS=-L/home/pletzer/software/netcdf-4.1.1-ser/lib CPPFLAGS=-I/home/pletzer/software/netcdf-4.1.1-ser/include"


ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I $(top_srcdir)/src -I $(top_srcdir)/gridspec_api

SUBDIRS = include src gridspec_api examples pycf

include_HEADERS = cf_config.h

lib_LTLIBRARIES = libcf.la
libcf_la_SOURCES = 
libcf_la_LIBADD = gridspec_api/host/libnccf_host.la			\
gridspec_api/mosaic/libnccf_mosaic.la					\
gridspec_api/regrid/libnccf_regrid.la					\
gridspec_api/data/libnccf_data.la                                       \
gridspec_api/grid/libnccf_grid.la	                                \
gridspec_api/coord/libnccf_coord.la					\
gridspec_api/axis/libnccf_axis.la					\
gridspec_api/global/libnccf_global.la                                   \
src/libcf_src.la -lnetcdf

# What needs to go in the binary dist?
BINFILES = README_BINARIES.txt include/libcf.h lib/libcf.a

# At Unidata, package up binaries.
ftpbin: install
	echo "Getting binaries from ${prefix}"
	ls ${prefix}
	echo "These are the binaries for libCF @PACKAGE_VERSION@." > ${prefix}/README_BINARIES.txt
	echo "For this build: CC=$(CC) CXX=$(CXX) FC=$(FC) F77=$(F77) $F90=$(F90)" >> ${prefix}/README_BINARIES.txt
	echo "CFLAGS=$(CFLAGS) CXXFLAGS=$(CXXFLAGS)" >> ${prefix}/README_BINARIES.txt
	echo "FCFLAGS=$(FCFLAGS) F77FLAGS=$(F77FLAGS) $F90FLAGS=$(F90FLAGS)" >> ${prefix}/README_BINARIES.txt
	tar cf binary-libcf-@PACKAGE_VERSION@.tar -C ${prefix} ${BINFILES}
	gzip -f binary-libcf-@PACKAGE_VERSION@.tar
