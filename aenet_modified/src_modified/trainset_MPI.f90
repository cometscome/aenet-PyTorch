!-----------------------------------------------------------------------
!           trainset.f90 - handling of the training set file
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!-----------------------------------------------------------------------
! 2013-05-09 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module trainset_MPI

  use aeio,    only: aeio_header,           &
                     TYPELEN, PATHLEN

  use io,      only: io_adjustl,            &
                     io_unit

  use sfsetup, only: Setup,                 &
                     save_Setup,            &
                     load_Setup,            &
                     del_Setup,             &
                     stp_init,              &
                     stp_final,             &
                     stp_normalize,         &
                     stp_assert_moduleinit, &
                     stp_nsf_max

  use trainset, only: TrnSet

  use parallel,    only: pp_init,                &
                         pp_final,               &
                         pp_bcast,               &
                         pp_barrier,             &
                         ppMaster,               &
                         ppRank,                 &
                         ppSize,&
                         pp_calc_and_bcast_footer

  implicit none
  private
  save

  public  :: new_TrnSet_MPI,          &
             close_TrnSet_MPI,        &
             save_TrnSet_info_MPI,    &
             ts_write_atom_info_MPI,  &
             ts_write_header_MPI,     &
             ts_write_sf_info_MPI,    &
             ts_write_structure_info_MPI, &
             ts_write_footer_MPI,     &
             ts_parallel_footer_MPI,   &
             ts_loadandwrite_structure_info_MPI

  private :: ts_assert_init,          &
             ts_assert_writemode,     &
             ts_assert_readmode

  !--------------------------------------------------------------------!
  ! Basis function values and derivatives may be read and written      !
  ! either using a basis function setup [type(Setup)] or directly into !
  ! double precision arrays of the correct dimensions.                 !
  !--------------------------------------------------------------------!

contains

  function new_TrnSet_MPI(nTypes, typeName, E_atom, nStrucs, file, scale, &
                      shift) result(ts)

    implicit none

    integer,                             intent(in) :: nTypes
    character(len=*), dimension(nTypes), intent(in) :: typeName
    double precision, dimension(nTypes), intent(in) :: E_atom
    integer,                             intent(in) :: nStrucs
    character(len=*),                    intent(in) :: file
    double precision, optional,          intent(in) :: scale
    double precision, optional,          intent(in) :: shift
    type(TrnSet)                                    :: ts

    logical :: fexists

!    if (ppMaster) then
!       inquire(file=trim(adjustl(file)), exist=fexists)
!       if (fexists) then
!          write(0,*) 'Error: file already exists: ', trim(adjustl(file))
!          stop
!       end if
!    end if

    allocate(ts%typeName(nTypes), ts%E_atom(nTypes))

    ts%file               = trim(adjustl(file))
    ts%nTypes             = nTypes
    ts%typeName(1:nTypes) = typeName(1:nTypes)
    ts%E_atom(1:nTypes)   = E_atom(1:nTypes)
    ts%nStrucs            = nStrucs
    ts%iStruc             = 0
    if (present(scale) .and. present(shift)) then
       ts%normalized = .true.
       ts%scale = scale
       ts%shift = shift
    else
       ts%normalized = .false.
       ts%scale = 1.0d0
       ts%shift = 0.0d0
    end if

    ts%unit   = io_unit()
!    open(ts%unit, file=trim(ts%file), status='new', action='write', &
!         form='unformatted')

    ts%nAtomsTot = 0

    ts%mode = 'write'
    ts%init = .true.

    call ts_write_header_MPI(ts)

  end function new_TrnSet_MPI

  !--------------------------------------------------------------------!

  subroutine close_TrnSet_MPI(ts, stp, status)

    implicit none

    type(TrnSet),                        intent(inout) :: ts
    type(Setup), dimension(:), optional, intent(in)    :: stp
    character(len=*),          optional, intent(in)    :: status

    if (.not. ts%init) return

    if (trim(ts%mode) == 'write') then
       if (present(stp)) then
          call ts_write_footer_MPI(ts, stp=stp)
       else
          call ts_write_footer_MPI(ts)
       end if
    end if

    if ((trim(ts%mode)=='read') .or. (trim(ts%mode)=='write')) then
       if (present(status)) then
          close(ts%unit, status=trim(status))
       else
          close(ts%unit)
       end if
    end if

    deallocate(ts%typeName, ts%E_atom)
    ts%init   = .false.

  end subroutine close_TrnSet_MPI

  !--------------------------------------------------------------------!
  !              only training set info - no actual data               !
  !--------------------------------------------------------------------!

  subroutine save_TrnSet_info_MPI(ts, file, unit)

    implicit none

    type(TrnSet),               intent(in) :: ts
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: u

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', &
            form='unformatted', action='write')
    else
       write(0,*) "Error: neither unit nor file specified in `save_TrnSet_info()'."
       stop
    end if

    write(u) ts%file
    write(u) ts%normalized
    write(u) ts%scale
    write(u) ts%shift
    write(u) ts%nTypes
    write(u) ts%typeName(1:ts%nTypes)
    write(u) ts%E_atom(1:ts%nTypes)
    write(u) ts%nAtomsTot
    write(u) ts%nStrucs
    write(u) ts%E_min, ts%E_max, ts%E_av

    if (.not. present(unit)) close(u)

  end subroutine save_TrnSet_info_MPI

  !--------------------------------------------------------------------!

  function new_TrnSet_info_MPI(nTypes) result(ts)

    implicit none

    integer,                             intent(in) :: nTypes
    type(TrnSet)                                    :: ts

    allocate(ts%typeName(nTypes), ts%E_atom(nTypes))

    ts%nTypes     = nTypes
    ts%normalized = .false.
    ts%scale      = 1.0d0
    ts%shift      = 1.0d0
    ts%file       = ''
    ts%unit       = -1
    ts%nAtomsTot  = 0
    ts%nStrucs    = 0
    ts%iStruc     = 0

    ts%init = .true.
    ts%mode = 'info'

  end function new_TrnSet_info_MPI

  !--------------------------------------------------------------------!
  !                      training set file header                      !
  !--------------------------------------------------------------------!

  subroutine ts_write_header_MPI(ts)

    implicit none

    type(TrnSet), intent(inout) :: ts
    integer :: ts_unit = 1234
    integer :: ts_unitall = 12345

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    if (ppMaster) then
       open (ts_unit, file='ts.header', form='unformatted', position='append')
       write(ts_unit) ts%nTypes
       write(ts_unit) ts%nStrucs
       write(ts_unit) ts%typeName(:)
       write(ts_unit) ts%E_atom(:)
       write(ts_unit) ts%normalized
       write(ts_unit) ts%scale
       write(ts_unit) ts%shift
       close(ts_unit)

       open (ts_unitall, file='ts.all', form='unformatted', status='replace')
       write(ts_unitall) ts%nTypes
       write(ts_unitall) ts%nStrucs
       write(ts_unitall) ts%typeName(:)
       write(ts_unitall) ts%E_atom(:)
       write(ts_unitall) ts%normalized
       write(ts_unitall) ts%scale
       write(ts_unitall) ts%shift
       close(ts_unitall)
    end if

  end subroutine ts_write_header_MPI

  !--------------------------------------------------------------------!
  !           training set file footer containing statistics           !
  !--------------------------------------------------------------------!

  subroutine ts_parallel_footer_MPI(ts, stp)

    implicit none

    type(TrnSet),                        intent(inout) :: ts
    type(Setup), dimension(:), intent(inout) :: stp

   call pp_calc_and_bcast_footer(ts,stp)

  end subroutine ts_parallel_footer_MPI

  subroutine ts_write_footer_MPI(ts, stp)

    implicit none

    type(TrnSet),                        intent(inout) :: ts
    type(Setup), dimension(:), optional, intent(in) :: stp

    integer :: itype, nTypes
    logical :: has_setups
    integer :: ts_unit = 1234
    integer :: ts_unitall = 12345

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    if (ppMaster) then

       open (ts_unit, file='ts.footer', form='unformatted', position='append')
       open (ts_unitall, file='ts.all', form='unformatted', position='append')

       write(ts_unit) ts%nAtomsTot
       write(ts_unitall) ts%nAtomsTot
       write(ts_unit) ts%E_av, ts%E_min, ts%E_max
       write(ts_unitall) ts%E_av, ts%E_min, ts%E_max
!       write(5678,*) ts%nAtomsTot
!       write(5678,*) ts%E_av, ts%E_min, ts%E_max
!
       if (present(stp)) then
          nTypes = size(stp(:))
          if (nTypes /= ts%nTypes) then
             write(0,*) "Error: wrong size of array stp in `ts_read_footer()'."
             stop
          end if
          has_setups = .true.
          write(ts_unit) has_setups
          write(ts_unitall) has_setups
!          write(5678,*) has_setups
          do itype = 1, ts%nTypes
             write(ts_unit) itype
             write(ts_unitall) itype
!             write(5678,*) itype
             call save_Setup(stp(itype), unit=ts_unit)
             call save_Setup(stp(itype), unit=ts_unitall)
          end do
       else
          has_setups = .false.
!          write(5678,*) has_setups
          write(ts_unit) has_setups
          write(ts_unitall) has_setups
       end if

       close(ts_unit)
       close(ts_unitall)

    end if

  end subroutine ts_write_footer_MPI

  subroutine ts_loadandwrite_structure_info_MPI(ifile)

   implicit none

   integer:: nAtoms
   integer:: nTypes
   double precision:: energy
   integer,          intent(in)    :: ifile
   integer :: ts_unit = 1234
   integer :: ts_unitall = 12345
   character(len=100)::filenameload
   integer::iatom
   integer::itype
   double precision, dimension(3):: cooCart
   double precision, dimension(3):: forCart

   integer:: nsf
   double precision, dimension(:),allocatable:: sfval

   double precision :: E_atom
   integer::l


   open (ts_unitall, file='ts.all', form='unformatted', position='append')
   open (ts_unit, file='ts.'//io_adjustl(ifile), form='unformatted', status='old')

   read(ts_unit) l
   filenameload = " "
   l = min(l,len(filenameload))
   read(ts_unit) filenameload(1:l)

   write(ts_unitall) len_trim(filenameload)
   write(ts_unitall) trim(filenameload)



   read(ts_unit) nAtoms, nTypes
   write(ts_unitall) nAtoms, nTypes

   read(ts_unit) energy
   write(ts_unitall) energy

   atoms : do iatom = 1, nAtoms

      read(ts_unit) itype
      write(ts_unitall) itype

      read(ts_unit) cooCart(1:3)
      write(ts_unitall) cooCart(1:3)
      read(ts_unit) forCart(1:3)
      write(ts_unitall) forCart(1:3)

      read(ts_unit) nsf
      write(ts_unitall) nsf
      allocate(sfval(1:nsf))

      read(ts_unit) sfval(1:nsf)
      write(ts_unitall) sfval(1:nsf)
      deallocate(sfval)

   end do atoms

   close(ts_unit)
   close(ts_unitall)


 end subroutine ts_loadandwrite_structure_info_MPI

  !--------------------------------------------------------------------!
  !              data from structures in the training set              !
  !--------------------------------------------------------------------!

  subroutine ts_write_structure_info_MPI(ts, filename, nAtoms, nTypes, energy, ifile)

    implicit none

    type(TrnSet),     intent(inout) :: ts
    character(len=*), intent(in)    :: filename
    integer,          intent(in)    :: nAtoms
    integer,          intent(in)    :: nTypes
    double precision, intent(in)    :: energy
    integer,          intent(in)    :: ifile
    integer :: ts_unit = 1234

    double precision :: E_atom

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    if (ts%iStruc >= ts%nStrucs) then
       write(0,*) "Error: too many files for training set."
       stop
    else
       ts%iStruc = ts%iStruc + 1
    end if

    open (ts_unit, file='ts.'//io_adjustl(ifile), form='unformatted', position='append')
    write(ts_unit) len_trim(filename)
    write(ts_unit) trim(filename)
    write(ts_unit) nAtoms, nTypes
    write(ts_unit) energy
    close(ts_unit)

    ! energy stats
    E_atom = energy/dble(nAtoms)
    if (ts%iStruc > 1) then
       ts%E_min = min(ts%E_min, E_atom)
       ts%E_max = max(ts%E_max, E_atom)
       ts%E_av  = ts%E_av + E_atom/dble(ts%nStrucs)
    else
       ts%E_min = E_atom
       ts%E_max = E_atom
       ts%E_av  = E_atom/dble(ts%nStrucs)
    end if

    ! keep track of the atoms in the training set
    ts%nAtomsTot = ts%nAtomsTot + nAtoms

  end subroutine ts_write_structure_info_MPI


  !--------------------------------------------------------------------!

  subroutine ts_write_atom_info_MPI(ts, itype, cooCart, forCart, ifile)

    implicit none

    type(TrnSet),                   intent(inout) :: ts
    integer,                        intent(in)    :: itype
    double precision, dimension(3), intent(in)    :: cooCart
    double precision, dimension(3), intent(in)    :: forCart
    integer,          intent(in)    :: ifile
    integer :: ts_unit = 1234

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    open (ts_unit, file='ts.'//io_adjustl(ifile), form='unformatted', position='append')
    write(ts_unit) itype
    write(ts_unit) cooCart(1:3)
    write(ts_unit) forCart(1:3)
    close(ts_unit)

  end subroutine ts_write_atom_info_MPI

  !--------------------------------------------------------------------!

  subroutine ts_write_sf_info_MPI(ts, nsf, sfval, ifile)

    implicit none

    type(TrnSet),                       intent(inout) :: ts
    integer,                            intent(in)    :: nsf
    double precision, dimension(nsf),   intent(in)    :: sfval
    integer,          intent(in)    :: ifile
    integer :: ts_unit = 1234

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    open (ts_unit, file='ts.'//io_adjustl(ifile), form='unformatted', position='append')
    write(ts_unit) nsf
    write(ts_unit) sfval(1:nsf)
    close(ts_unit)

  end subroutine ts_write_sf_info_MPI

  !--------------------------------------------------------------------!
  !                            state checks                            !
  !--------------------------------------------------------------------!

  subroutine ts_assert_init(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (.not. ts%init) then
       write(0,*) "Error: training set not initialized."
       stop
    end if
  end subroutine ts_assert_init

  subroutine ts_assert_writemode(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (trim(ts%mode) /= 'write') then
       write(0,*) "Error: training set not in 'write' mode."
       stop
    end if
  end subroutine ts_assert_writemode

  subroutine ts_assert_readmode(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (trim(ts%mode) /= 'read') then
       write(0,*) "Error: training set not in 'read' mode."
       stop
    end if
  end subroutine ts_assert_readmode

end module trainset_MPI
