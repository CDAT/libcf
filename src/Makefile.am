## This is an automake file, part of Unidata's NetCDF CF Library package.
# Copyright 2006-2011, see the COPYRIGHT file for more information.

# This automake file is in charge of building the src directory.

# Test programs will need to link to this library (libcf) and also
# netcdf.
LDADD = ${lib_LTLIBRARIES}
LDADD += -lnetcdf

# Initalize these.
AM_CPPFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS =

include_HEADERS =          \
nccf_find_next_indices.h   \
nccf_keyvalue.h            

# This is our output, the NetCDF CF Library.
lib_LTLIBRARIES = libcf_src.la

# I think we always need this one
LDADD += -lm

# These linker flags specify libtool version info.
libcf_src_la_LDFLAGS = -version-info 0:0:0

check_PROGRAMS = tst_coords tst_files tst_vars tst_cvars tst_vert	\
tst_subsets tst_princeton tst_listitem tst_nccf_get_flat_index		\
tst_nccf_index2vector tst_nccf_varObj tst_solve tst_linear_interp	\
tst_nccf_get_bipolar_cap tst_nccf_get_tripolar_grid			\
tst_nccf_generate_id tst_nccf_get_tripolar_halfgrid			\
tst_nccf_get_cubedsphere_grid tst_nccf_find_adjacent			\
tst_nccf_linear_interp_mask tst_nccf_keyvalue \
tst_nccf_get_flat_index tst_nccf_index2vector		\
tst_nccf_varObj tst_solve tst_linear_interp tst_nccf_get_bipolar_cap	\
tst_nccf_get_tripolar_grid tst_nccf_generate_id				\
tst_nccf_get_tripolar_halfgrid tst_nccf_get_cubedsphere_grid		\
tst_nccf_find_adjacent tst_nccf_linear_interp_mask tst_nccf_keyvalue

#noinst_PROGRAMS = tst_nccf_get_flat_index tst_nccf_index2vector		

if BUILD_F90
lib_LTLIBRARIES += libcff_src.la
libcff_src_la_DEPENDENCIES = libcf_api.f90
EXTRA_DIST = libcf_api.f90
libcff_src_la_SOURCES = cffcvars.f90
# if USE_NETCDF4
# check_PROGRAMS += tst_f03_file
# tst_f03_file_SOURCES = tst_f03_file.f90
# endif # USE_NETCDF4
LDADD += libcf_src.la
endif # BUILD_F90

# This is the coordinate systems stuff.
libcf_src_la_SOURCES = nccoord.c libcf_int.h cffiles.c cfvars.c		\
cfcvars.c cfsubsets.c cferror.c cdTimeConv.c cdmsint.h timeArith.c	\
timeConv.c libcfAddOn.c cflistitem.c cdms.h cddrs.h drscdf.h cdunif.h

# in memory variable, interpolation, indexing, ...
libcf_src_la_SOURCES += nccf_varObj.c nccf_varSetData_generic.h		\
nccf_strerror.c nccf_solve.h nccf_solve.c nccf_put_var_slice_text.c	\
nccf_linesearch_indices.h nccf_linesearch_indices.c			\
nccf_linear_interp.h nccf_linear_interp.c				\
nccf_index_to_corner_vector.c nccf_index2vector.c nccf_handle_error.c	\
nccf_get_tripolar_halfgrid.c nccf_get_tripolar_grid.c			\
nccf_get_start_end_bound_indices.c nccf_get_position.h			\
nccf_get_position.c nccf_get_multi_index.c nccf_get_linear_weights.h	\
nccf_get_linear_weights.c nccf_get_flat_index.c				\
nccf_get_cubedsphere_grid.c nccf_get_bipolar_cap.c nccf_generate_id.c	\
nccf_find_next_indices.h nccf_find_next_indices.c			\
nccf_find_next_adjcnt_index.c nccf_find_next_adjacent.c			\
nccf_find_indices.h nccf_find_indices.c nccf_add_id_to_files.c		\
nccf_keyvalue.h nccf_keyvalue.c                                         \
nccf_get_distance_in_coord.c nccf_get_distance_in_coord.h               \
nccf_linear_interp_mod.h nccf_linear_interp_mod.c			

# These files are created by the tests.
CLEANFILES = tst_*.nc libcf.mod

# created by one of the unit tests
clean-local:
	-rm -rf tst_linear_interp_tripolar.py

if EXTRA_DATA_TESTS
check_PROGRAMS += tst_example_data
endif # EXTRA_DATA_TESTS

# Rule to build Fortran 2003 code.
#cfcvars.o:cfcvsars.f90
#	$(FC) -c $(AM_FCFLAGS) $(FCFLAGS) -c $? -o $@

tst_f03_file.o: tst_f03_file.f90
	$(FC) -c $(AM_FCFLAGS) $(FCFLAGS) -c $? -o $@

TESTS = $(check_PROGRAMS)



