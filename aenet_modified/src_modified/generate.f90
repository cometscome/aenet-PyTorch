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
    use aenet_generate 


    implicit none
    character(len=PATHLEN)                         :: inFile
    integer::ionum


    ionum= 139
    call initialize(inFile)
    call generate_subroutine(inFile,ionum)


    contains !=============================================================!
  
  
    subroutine initialize(inFile)
  
      implicit none
  
      character(len=*), intent(out) :: inFile
  
      integer :: nargs
      logical :: fexists
  
      call aeio_header("generate.x - training set generation", char='=')
      write(ionum,*)
  
      call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')
  
      nargs = command_argument_count()
      if (nargs < 1) then
         write(0,*) "Error: No input file provided."
         call print_usage()
         call finalize()
         stop
      end if
  
      call get_command_argument(1, value=inFile)
      inquire(file=trim(inFile), exist=fexists)
      if (.not. fexists) then
         write(0,*) "Error: File not found: ", trim(inFile)
         call print_usage()
         call finalize()
         stop
      end if
  
    end subroutine initialize

    subroutine print_usage()

        implicit none

        write(ionum,*)
        write(ionum,*) "generate.x -- Generate training sets for use with `train.x'"
        write(ionum,'(1x,70("-"))')
        write(ionum,*) 'Usage: generate.x <input-file>'
        write(ionum,*)
        write(ionum,*) 'See the documentation or the source code for a description of the '
        write(ionum,*) 'input file format.'
        write(ionum,*)

  end subroutine print_usage

  
    !--------------------------------------------------------------------!
  
    subroutine finalize()
  
      implicit none
  
      integer :: itype
  
  
      call aeio_header(aeio_timestamp(), char=' ')
      call aeio_header("Training set generation done.", char='=')
  
    end subroutine finalize
end program generate