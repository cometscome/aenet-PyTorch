!-----------------------------------------------------------------------
!       predict.f90 - predict atomic energies of input structure
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
! 2011-11-17 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

program predict

  use aeio,      only: aeio_header,                    &
                       aeio_timestamp,                 &
                       aeio_print_copyright

  use aenet,     only: aenet_init,                     &
                       aenet_final,                    &
                       aenet_atomic_energy,            &
                       aenet_atomic_energy_and_forces_novirial, &
                       aenet_convert_atom_types,       &
                       aenet_free_atom_energy,         &
                       aenet_load_potential,           &
                       aenet_print_info,               &
                       aenet_Rc_min, aenet_Rc_max,     &
                       aenet_nnb_max

  use constants, only: PI

  use geometry,  only: geo_init,                       &
                       geo_final,                      &
                       pbc,                            &
                       latticeVec,                     &
                       recLattVec,                     &
                       geo_update_bounds,              &
                       origin,                         &
                       nAtoms,                         &
                       nTypes,                         &
                       atomType,                       &
                       atomTypeName,                   &
                       cooLatt

  use input,     only: InputData,                      &
                       read_InpPredict

  use io,        only: io_adjustl

  use lclist,    only: lcl_init,                       &
                       lcl_final,                      &
                       lcl_nmax_nbdist,                &
                       lcl_nbdist_cart

  use optimize,  only: opt_init,                       &
                       opt_final,                      &
                       opt_optimize_coords

  use parallel,  only: pp_init,                        &
                       pp_final,                       &
                       pp_bcast,                       &
                       pp_bcast_coo,                   &
                       pp_print_info,                  &
                       pp_bcast_InputData,             &
                       pp_bcast_latt,                  &
                       pp_sum,                         &
                       ppMaster, ppRank, ppSize
  use aenet_mpimodule
  use aenet_predict,only:predict_subroutine,print_usage

  implicit none

  !--------------------------------------------------------------------!
  ! A '*' in front of the variable name means that it is a broadcasted !
  ! variable and has the same value on each process.  A '+' means that !
  ! an array is allocated on all parallel processes, but does not      !
  ! necessarily have the same contents.                                !
  !                                                                    !
  !----------------------------- general ------------------------------!
  ! inp             structure with input data                          !
  ! inFile          name of the input file                             !
  !                                                                    !
  !---------------------------- structures ----------------------------!
  !*nFiles          number of input files/structures                   !
  ! cooFile         file name of structure file (atomic coordinates)   !
  !                                                                    !
  !------------------------------ output ------------------------------!
  ! Ecoh            cohesive energy of the current structure           !
  ! Etot            total energy                                       !
  !+forCart         cartesian atomic forces of the current structure   !
  !+atomicEnergy    energy of the individual atoms in the structure    !
  !--------------------------------------------------------------------!

  type(InputData)                               :: inp

  character(len=1024)                           :: inFile, cooFile, strucFile

  integer,          dimension(:),   allocatable :: atomType_orig

  integer                                       :: istruc, nStrucs

  double precision                              :: Ecoh, Etot, E0
  double precision, dimension(:,:), allocatable :: forCart
  double precision, dimension(:),   allocatable :: atomicEnergy

  double precision, dimension(3)                :: F_mav, F_max, F_avg
  double precision                              :: F_rms, F_rms_prev
  double precision, dimension(3)                :: dmax
  integer                                       :: imax

  integer                                       :: iter, stat
  logical                                       :: conv


  !-------------------------- initialization --------------------------!

  call initialize_MPI(inFile, strucFile, inp)
  call predict_subroutine(inFile, strucFile, inp)

   contains

   subroutine initialize_MPI(inFile, strucFile, inp)
 
     implicit none
 
     character(len=*), intent(out) :: inFile
     character(len=*), intent(out) :: strucFile
     type(InputData),  intent(out) :: inp
 
 
     logical :: fexists
     integer :: nargs
     integer :: stat
     integer :: itype
     integer::ierr

#ifdef PARALLEL
      call MPI_Init(ierr)
#endif     
 
     call pp_init()
 
     if (ppMaster) then
        nargs = command_argument_count()
        if (nargs < 1) then
           write(0,*) "Error: No input file specified."
           call print_usage()
           !call finalize()
           stop
        end if
 
        call get_command_argument(1, value=inFile)
        inquire(file=trim(inFile), exist=fexists)
        if (.not. fexists) then
           write(0,*) "Error: File not found: ", trim(inFile)
           call print_usage()
           !call finalize()
           stop
        end if
 
        ! read name of structure from command line, if present
        if (nargs > 1) then
           call get_command_argument(2, value=strucFile)
        else
           strucFile = ''
        end if
 
        ! read general input file
        inp = read_InpPredict(inFile)
     end if
     call pp_bcast(inFile)
     call pp_bcast(strucFile)
     call pp_bcast_InputData(inp)
 
     if (inp%verbosity > 0) call pp_print_info()
 
     ! initialize aenet
     call aenet_init(inp%typeName, stat)
     if (stat /= 0) then
        write(0,*) 'Error: aenet initialization failed'
        !call finalize()
        stop
     end if
 
     ! load ANN potentials
     do itype = 1, inp%nTypes
        call aenet_load_potential(itype, inp%netFile(itype), stat)
        if (stat /= 0) then
        write(0,*) 'Error: could not load ANN potentials'
           !call finalize()
           stop
        end if
     end do
 
     if (ppMaster .and. (inp%verbosity > 0)) then
        ! write header and copyright info
        call aeio_header("Atomic Energy Network Interpolation", char='=')
        call aeio_header(aeio_timestamp(), char=' ')
        write(*,*)
        call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')
     end if
 
   end subroutine initialize_MPI


end program predict
