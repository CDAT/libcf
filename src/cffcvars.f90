module libcf
  use iso_c_binding
  implicit none

  ! For some odd reason, netCDF f90 has NF90_GLOBAL=0, but NC_GLOBAL
  ! is really -1. WTF?
  integer, parameter:: NC_GLOBAL = -1

  ! Max len of history, institution, comment, etc.
  integer, parameter :: CF_MAX_LEN = 1024

  ! Max number of allowed coordinate variables
  integer, parameter :: CF_MAX_COORDS = 64

  ! Need a null pointer to pass to C for non-present strings and a
  ! null terminator for initializing strings.
  character, pointer :: null_ptr => NULL()
  character(len = 1), parameter:: NT = achar(0)

  ! These interfaces map to the C library calls.
  include 'libcf_api.f90'

contains

  ! These functions are the new libcf F90 API.
  function cf_def_file(ncid, convention, title, history, institution, &
       source, comment, references)
    use iso_c_binding
    use netcdf
    implicit none
    integer, intent(in) :: ncid
    character(len = *), optional, intent(in):: convention, title, history, &
         institution, source, comment, references
    character(len = CF_MAX_LEN), pointer :: p_history => NULL(), p_title => NULL()
    character(len = CF_MAX_LEN), target :: c_history, c_title
    character(len = CF_MAX_LEN), pointer :: p_institution => NULL(), p_source => NULL(), &
         p_comment => NULL(), p_references => NULL()
    character(len = CF_MAX_LEN), target :: c_institution, c_source, c_comment, &
         c_references
    integer:: cf_def_file

    cf_def_file = NF90_NOERR

    ! Handle convention.
    if (present(convention)) then
       cf_def_file = nccf_def_convention(ncid)
       if (cf_def_file .ne. 0) return
    endif

    ! Handle title and history.
    if (present(title) .or. present(history)) then
       if (present(title)) then
          p_title => c_title
          p_title = trim(title)//c_null_char
       endif
       if (present(history)) then
          p_history => c_history
          p_history = trim(history)//c_null_char
       endif
       cf_def_file = nccf_def_file(ncid, p_title, p_history)
       if (cf_def_file .ne. 0) return
    endif

    ! Handle institution, source, comment, and references.
    if (present(institution) .or. present(source) .or. present(comment) .or. &
         present(references)) then
       if (present(institution)) then
          p_institution => c_institution
          p_institution = trim(institution)//c_null_char
       endif
       if (present(source)) then
          p_source => c_source
          p_source = trim(source)//c_null_char
       endif
       if (present(comment)) then
          p_comment => c_comment
          p_comment = trim(comment)//c_null_char
       endif
       if (present(references)) then
          p_references => c_references
          p_references = trim(references)//c_null_char
       endif
       cf_def_file = nccf_def_notes(ncid, NC_GLOBAL, p_institution, p_source, &
            p_comment, p_references)
    endif
  end function cf_def_file

  function cf_def_var(ncid, varid, units, long_name, standard_name, &
       ncoord_vars, coord_varids, institution, source, comment, references)
    use iso_c_binding
    use netcdf
    implicit none
    integer(c_int), intent(in) :: ncid, varid
    integer, optional, intent(in) :: ncoord_vars
    integer, target :: c_ncoord_vars = 0
    integer, optional, dimension(:), intent(in) :: coord_varids
    integer, target, dimension(CF_MAX_COORDS) :: c_coord_varids
    character(len = *), optional, intent(in) :: institution, source, &
         comment, references, units, long_name, standard_name
    character(len = CF_MAX_LEN), target :: c_units, c_long_name, c_standard_name
    character(len = CF_MAX_LEN), target :: c_institution, c_source, &
         c_comment, c_references
    character(len = CF_MAX_LEN), pointer :: p_units => NULL(), &
         p_long_name => NULL(), p_standard_name => NULL()
    character(len = CF_MAX_LEN), pointer :: p_institution => NULL(), p_source => NULL(), &
         p_comment => NULL(), p_references => NULL()
    integer  :: cf_def_var, i

    cf_def_var = NF90_NOERR

    ! Check input consistency.
    if (present(ncoord_vars) .and. .not. present(coord_varids)) then
       cf_def_var = NF90_EINVAL
       return
    endif

    ! Handle institution, source, comment, and references.
    if (present(institution) .or. present(source) .or. present(comment) .or. &
         present(references)) then
       if (present(institution)) then
          p_institution => c_institution
          p_institution = trim(institution)//c_null_char
       endif
       if (present(source)) then
          p_source => c_source
          p_source = trim(source)//c_null_char
       endif
       if (present(comment)) then
          p_comment => c_comment
          p_comment = trim(comment)//c_null_char
       endif
       if (present(references)) then
          p_references => c_references
          p_references = trim(references)//c_null_char
       endif
       cf_def_var = nccf_def_notes(ncid, varid - 1, p_institution, p_source, &
            p_comment, p_references)
       if (cf_def_var .ne. 0) return
    endif

    ! Handle units, long_name, standard_name.
    if (present(units) .or. present(long_name) .or. present(standard_name) .or. &
         present(ncoord_vars) .or. present(coord_varids)) then
       if (present(units)) then
          p_units => c_units
          p_units = trim(units)//c_null_char
       endif
       if (present(long_name)) then
          p_long_name => c_long_name
          p_long_name = trim(long_name)//c_null_char
       endif
       if (present(standard_name)) then
          p_standard_name => c_standard_name
          p_standard_name = trim(standard_name)//c_null_char
       endif
       if (present(ncoord_vars)) c_ncoord_vars = ncoord_vars
       if (present(coord_varids)) then
          do i = 1, ncoord_vars
             c_coord_varids(i) = coord_varids(i) - 1
          end do
       endif

       cf_def_var = nccf_def_var(ncid, varid - 1, p_units, p_long_name, &
            p_standard_name, c_ncoord_vars, c_coord_varids)
    endif
  end function cf_def_var
  
  function cf_inq_var(ncid, varid, units, long_name, standard_name, &
       ncoord_vars, coord_varids, institution, source, comment, references)
    use iso_c_binding
    use netcdf
    implicit none
    integer, intent(in) :: ncid, varid
    character(len = *), optional, intent(out) :: units, long_name, standard_name
    character(len = CF_MAX_LEN), target :: c_units, c_long_name, c_standard_name
    character(len = CF_MAX_LEN), pointer :: p_units => NULL(), p_long_name => NULL(), &
         p_standard_name => NULL()
    integer(c_size_t) :: c_units_len, c_long_name_len, c_standard_name_len
    character(len = *), optional, intent(out) :: institution, source, comment, references
    character(len = CF_MAX_LEN), target :: c_institution, c_source, c_comment, c_references
    integer(c_size_t) :: c_institution_len, c_source_len, c_comment_len, c_references_len
    character(len = CF_MAX_LEN), pointer :: p_institution => NULL(), p_source => NULL(), &
         p_comment => NULL(), p_references => NULL()
    integer, intent(out), optional :: ncoord_vars
    integer, dimension(*), optional, intent(out) :: coord_varids
    integer :: c_ncoord_vars
    integer, dimension(CF_MAX_COORDS) :: c_coord_varids
    integer :: cf_inq_var, i

    if (present(units)) p_units => c_units
    if (present(long_name)) p_long_name => c_long_name
    if (present(standard_name)) p_standard_name => c_standard_name

    cf_inq_var = nccf_inq_var(ncid, varid - 1, c_units_len, p_units, c_long_name_len, &
         p_long_name, c_standard_name_len, p_standard_name, c_ncoord_vars, c_coord_varids)

    if (present(units)) then
       units = trim(c_units)
       if (index(units, C_NULL_CHAR) .gt. 0) then
          units(index(units, C_NULL_CHAR):len(units)) = ' '
       endif
    endif
    if (present(long_name)) then
       long_name = trim(c_long_name)
       if (index(long_name, C_NULL_CHAR) .gt. 0) then
          long_name(index(long_name, C_NULL_CHAR):len(long_name)) = ' '
       endif
    endif
    if (present(standard_name)) then
       standard_name = trim(c_standard_name)
       if (index(standard_name, C_NULL_CHAR) .gt. 0) then
          standard_name(index(standard_name, C_NULL_CHAR):len(standard_name)) = ' '
       endif
    endif
    if (present(ncoord_vars)) ncoord_vars = c_ncoord_vars
    if (present(coord_varids)) then
       do i = 1, c_ncoord_vars
          coord_varids(i) = c_coord_varids(i) + 1
       end do
    endif

    if (present(institution)) p_institution => c_institution
    if (present(source)) p_source => c_source 
    if (present(comment)) p_comment => c_comment
    if (present(references)) p_references => c_references

    ! Get the institution and source strings.
    cf_inq_var = nccf_inq_notes(ncid, varid - 1, c_institution_len, p_institution, &
         c_source_len, p_source, c_comment_len, p_comment, c_references_len, &
         p_references)
    if (cf_inq_var .ne. NF90_NOERR) return

    ! Copy the strings.
    if (present(institution)) then
       institution = p_institution(:len(institution))
       if (index(institution, NT) .gt. 0) then
          institution(index(institution, NT):len(institution)) = ' '
       endif
    endif
    if (present(source)) then
       source = p_source(:len(source))
       if (index(source, NT) .gt. 0) then
          source(index(source, NT):len(source)) = ' '
       endif
    endif
    if (present(comment)) then
       comment = p_comment(:len(comment))
       if (index(comment, NT) .gt. 0) then
          comment(index(comment, NT):len(comment)) = ' '
       endif
    endif
    if (present(references)) then
       references = p_references(:len(references))
       if (index(references, NT) .gt. 0) then
          references(index(references, NT):len(references)) = ' '
       endif
    endif
  end function cf_inq_var

  ! Input string MUST have length less then CF_MAX_LEN or universe
  ! might collapse into a singularity.
  function cstr(in_str)
    use iso_c_binding
    implicit none
    character(len = *), intent(in):: in_str
    character(len = CF_MAX_LEN) :: cstr

    !print *, 'in_str=', in_str
    cstr = achar(0)
    if (len(in_str) .gt. CF_MAX_LEN) return
    
    cstr = in_str
    cstr(len_trim(in_str) + 1:len_trim(in_str) + 1) = achar(0)
  end function cstr

  function cf_add_history(ncid, history)
    use iso_c_binding
    implicit none
    integer, intent(in)::ncid
    character (len = *), intent(in):: history
    integer:: cf_add_history
    integer(c_int):: c_ncid
    character(len = 256):: c_history

    c_ncid = ncid
    c_history = history
    c_history(len_trim(history) + 1:len_trim(history) + 1) = achar(0)
    cf_add_history = nccf_add_history(ncid, c_history)

  end function cf_add_history

  function cf_inq_file(ncid, title, history, institution, source, &
       comment, references)
    use iso_c_binding
    use netcdf
    implicit none
    integer, intent(in):: ncid
    character (len = *), optional, intent(out):: title, history, institution
    character (len = *), optional, intent(out):: source, comment, references
    integer:: cf_inq_file
    integer(c_size_t):: c_title_len, c_history_len, c_institution_len
    integer(c_size_t):: c_source_len, c_comment_len, c_references_len
    character(len = CF_MAX_LEN), target :: c_history, c_title, c_institution, c_source, &
         c_comment, c_references
    character(len = CF_MAX_LEN), pointer :: p_history => NULL(), p_title => NULL(), &
         p_institution => NULL(), p_source => NULL(), p_comment => NULL(), p_references => NULL()

    cf_inq_file = NF90_NOERR
    if (present(title) .or. present(history)) then
       
       if (present(title)) p_title => c_title
       if (present(history)) p_history => c_history

       cf_inq_file = nccf_inq_file(ncid, c_title_len, p_title, c_history_len, &
            p_history)
       if (cf_inq_file .ne. NF90_NOERR) return

       ! Copy the strings.
       if (present(title)) then
          title = c_title(:len(title))
          if (index(title, NT) .gt. 0) then
             title(index(title, NT):len(title)) = ' '
          endif
       endif
       if (present(history)) then
          history = c_history(:len(history))
          if (index(history, NT) .gt. 0) then
             history(index(history, NT):len(history)) = ' '
          endif
       endif

    endif ! present(title) .or. present(history)

    if (present(institution)) p_institution => c_institution
    if (present(source)) p_source => c_source 
    if (present(comment)) p_comment => c_comment
    if (present(references)) p_references => c_references

    ! Get the institution and source strings.
    cf_inq_file = nccf_inq_notes(ncid, NC_GLOBAL, c_institution_len, p_institution, &
         c_source_len, p_source, c_comment_len, p_comment, c_references_len, &
         p_references)
    if (cf_inq_file .ne. NF90_NOERR) return

    ! Copy the strings.
    if (present(institution)) then
       institution = p_institution(:len(institution))
       if (index(institution, NT) .gt. 0) then
          institution(index(institution, NT):len(institution)) = ' '
       endif
    endif
    if (present(source)) then
       source = p_source(:len(source))
       if (index(source, NT) .gt. 0) then
          source(index(source, NT):len(source)) = ' '
       endif
    endif
    if (present(comment)) then
       comment = p_comment(:len(comment))
       if (index(comment, NT) .gt. 0) then
          comment(index(comment, NT):len(comment)) = ' '
       endif
    endif
    if (present(references)) then
       references = p_references(:len(references))
       if (index(references, NT) .gt. 0) then
          references(index(references, NT):len(references)) = ' '
       endif
    endif
    
    end function cf_inq_file
end module libcf
