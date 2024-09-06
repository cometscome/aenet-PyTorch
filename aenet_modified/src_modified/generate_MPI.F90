!-----------------------------------------------------------------------
!      generate.f90 - generate training sets for use with train.x
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!+ ---------------------------------------------------------------------
!+ If you make use of AENET for your publication, please cite:
!+ [1] N. Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150.
!+ [2] J. Behler and M. Parrinello, Phys. Rev. Lett. 98 (2007) 146401.
!+
!+ If you used the Chebyshev descriptor, please cite:
!+ [3] N. Artrith, A. Urban, and G. Ceder, PRB 96 (2017) 014112.
!-----------------------------------------------------------------------
! 2011-10-19 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

program generate
   use aeio,     only: aeio_readline,        &
                   aeio_header,          &
                   aeio_timestamp,       &
                   aeio_print_copyright, &
                   PATHLEN, LINELEN
   use aenet_generate_MPI,only:generate_subroutine_MPI 

  use parallel,    only: pp_init,                &
                         pp_final,               &
                         pp_bcast,               &
                         pp_barrier,             &
                         ppMaster,               &
                         ppRank,                 &
                         ppSize,                 &
                         pp_bcast_InputData,     &
                         pp_bcast_Setup   

   use aenet_mpimodule


   implicit none
   character(len=PATHLEN)                         :: inFile
   integer::ionum
   integer::ierr


   ionum= 139
   call initialize_MPI(inFile)
   call generate_subroutine_MPI(inFile,ionum)


   contains 

   subroutine initialize_MPI(inFile)

    implicit none
#ifdef PARALLEL
    include 'mpif.h'
#endif

    character(len=*), intent(out) :: inFile

    integer :: nargs
    logical :: fexists
    logical :: stopnow
    integer::ierr

#ifdef PARALLEL
    call MPI_Init(ierr)
#endif

    call pp_init()

    stopnow = .false.

    if (ppMaster) then

       !call aeio_header("generate.x - training set generation", char='=')
       !write(*,*)

       !call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')

       nargs = command_argument_count()
       if (nargs < 1) then
          write(0,*) "Error: No input file provided."
          call print_usage()
          stopnow = .true.
       end if

       call get_command_argument(1, value=inFile)
       inquire(file=trim(inFile), exist=fexists)
       if (.not. fexists) then
          write(0,*) "Error: File not found: ", trim(inFile)
          stopnow = .true.
       end if

    end if

    call pp_bcast(stopnow)
    if (stopnow) then
       stop
    end if

    call pp_bcast(inFile)

  end subroutine initialize_MPI


  subroutine print_usage()

   implicit none

   write(*,*)
   write(*,*) "generate.x -- Generate training sets for use with `train.x'"
   write(*,'(1x,70("-"))')
   write(*,*) 'Usage: generate.x <input-file>'
   write(*,*)
   write(*,*) 'See the documentation or the source code for a description of the '
   write(*,*) 'input file format.'
   write(*,*)

 end subroutine print_usage

end program