!
! Example showing how to call a LibCF function from Fortran 2003
! $Id: $

program test
   use iso_c_binding
   implicit none

   integer, parameter            :: r8 = selected_real_kind(12,100)
   real(r8), allocatable, target :: lon(:)
   real(r8), parameter           :: lonmin = 0.0, lonmax = 360.0
   real(r8)                      :: dlon
   integer                       :: nlon, nlat, i, j, k, stat, sav, lonid
   integer, parameter            :: ndims = 2
   integer                       :: dims(ndims)
   character(len=128), target    :: dimnames(ndims)
   type(c_ptr)                   :: sa_dimnames(ndims)

   ! interface for calling C routine
   interface 
      function nccf_def_lon_coord(ndims, dims, sa_dimnames, coorddata, &
         &                           sav, coordid) result(stat) &
         &                           bind(C, name="nccf_def_lon_coord")
         use iso_c_binding
         implicit none
         integer(c_int), value     :: ndims
	     integer(c_int)            :: dims(*)
         type(c_ptr)               :: sa_dimnames(*)
         real(c_double)            :: coorddata(*)
         integer(c_int), value     :: sav
         integer(c_int)            :: coordid
         integer(c_int)            :: stat
       end function 

       function nccf_free_coord(coordid) result(stat) &
         &                           bind(C, name="nccf_free_coord")
          use iso_c_binding
          implicit none
          integer(c_int), value     :: coordid
          integer(c_int)            :: stat
       end function
    end interface

    ! create a longitude coordinate
    nlon = 11
    nlat = 12
    dlon = (lonmax - lonmin)/real(nlon-1, r8)
    allocate(lon(nlon*nlat))
    do j = 1, nlat
       do i = 1, nlon
          k = i + nlon*(j-1)
          lon(k) = lonmin + dlon*(i-1)
       enddo
    enddo
	
    ! note: index order is reversed between fortran and C
    dims(1) = nlat
    dims(2) = nlon

    dimnames(1) = 'nlat'//c_null_char
    dimnames(2) = 'nlon'//c_null_char

    sa_dimnames(1) = c_loc(dimnames(1))
    sa_dimnames(2) = c_loc(dimnames(2))
    sav = 1
    lonid = -12134567
    ! note that we pass the first element of lon since lon is array object 
    ! with a decorator
    stat = nccf_def_lon_coord(ndims, dims, sa_dimnames, lon(1), sav, lonid)
    print *,'stat = ', stat, ' lonid = ', lonid
	
    stat = nccf_free_coord(lonid)
    deallocate(lon)

    print *,'stat = ', stat

end program
