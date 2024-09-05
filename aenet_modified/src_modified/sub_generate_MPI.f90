module aenet_generate_MPI 
   use aenet_trainbin2ascii,only:trainbin2ascii_subroutine
implicit none
contains 

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

subroutine generate_subroutine_MPI(inFile,ionum)

    use aeio,     only: aeio_readline,        &
                        aeio_header,          &
                        aeio_timestamp,       &
                        aeio_print_copyright, &
                        PATHLEN, LINELEN
  
    use geometry, only: geo_init,          &
                        geo_final,         &
                        geo_itype_of_name, &
                        geo_type_conv,     &
                        pbc,               &
                        latticeVec,        &
                        nAtoms,            &
                        nTypes,            &
                        atomType,          &
                        atomTypeName,      &
                        cooLatt,           &
                        cooCart,           &
                        forCart,           &
                        hasEnergy,         &
                        hasForces,         &
                        cohesiveEnergy,    &
                        totalEnergy
  
    use input,    only: InputData,         &
                        read_InpGenerate,  &
                        del_InputData
  
    use io,       only: io_adjustl,        &
                        io_center,         &
                        io_lower,          &
                        io_readnext,       &
                        io_unit
  
    use lclist,   only: lcl_init,          &
                        lcl_final,         &
                        lcl_print_info,    &
                        lcl_nmax_nbdist,   &
                        lcl_nbdist_cart
  
    use sfsetup,  only: Setup,                 &
                        read_Setup_parameters, &
                        save_Setup,            &
                        del_Setup,             &
                        stp_init,              &
                        stp_final,             &
                        stp_get_range,         &
                        stp_print_info,        &
                        stp_eval,              &
                        nsf_max
  
    use timing,   only: tng_init,          &
                        tng_final,         &
                        tng_timing,        &
                        tng_timing2,       &
                        tng_timing3,       &
                        tng_dump
  
    use trainset, only: TrnSet,                 &
                        new_TrnSet,             &
                        close_TrnSet,           &
                        ts_print_info,          &
                        ts_write_header,        &
                        ts_write_sf_info,       &
                        ts_write_atom_info,     &
                        ts_write_structure_info
  
    use pytorchoutput, only: pyo_write_init,               &
                             pyo_write_final,              &
                             pyo_write_header_info,        &
                             pyo_write_structure_info,     &
                             pyo_write_atom_sf_info,       &
                             pyo_select_force_structures!, &
   use pytorchoutput_MPI, only: pyo_write_structure_info_MPI,     &
                             pyo_write_atom_sf_info_MPI,&
                             pyo_loadandwrite

    use trainset_MPI, only: new_TrnSet_MPI,             &
                             close_TrnSet_MPI,           &
                             ts_write_header_MPI,        &
                             ts_write_sf_info_MPI,       &
                             ts_write_atom_info_MPI,     &
                             ts_write_structure_info_MPI,&
                             ts_parallel_footer_MPI,&
                             ts_loadandwrite_structure_info_MPI

    use parallel,    only: pp_init,                &
                         pp_final,               &
                         pp_bcast,               &
                         pp_barrier,             &
                         ppMaster,               &
                         ppRank,                 &
                         ppSize,                 &
                         pp_bcast_InputData,     &
                         pp_bcast_Setup                             
  
    implicit none
  
    !--------------------------------------------------------------------!
    ! stp(i)         structural fingerprint basis setup for atom type i  !
    ! r_min, r_max   lower and upper bound for atomic interactions       !
    ! ts             training set reference                              !
    !                                                                    !
    ! nnb_max, nnb   max. and actual number of neighboring atoms         !
    ! nbcoo(i,j)     i-th component of the coordinates of the j-th       !
    !                neighboring atom                                    !
    ! nbdist(i)      distance of the i-th neighbor                       !
    !                                                                    !
    ! sfval(i)         value of the i-th basis function                  !
    ! sfderiv_i(i,j)   i-th component of the derivative of the j-th SF   !
    !                  with respect to the central atom                  !
    !                  sfderiv_i(3,nsf_max)                              !
    ! sfderiv_j(i,j,k) i-th component of the derivative of the j-th SF   !
    !                  with respect to the coordinates of atom k         !
    !                  sfderiv_j(3,nsf_max,nnb_max)                      !
    !                                                                    !
    ! E_coh          cohesive energy                                     !
    ! nFiles_inv     = 1/inp%nStrucs                                     !
    !                                                                    !
    ! inFile         name of the input file for the generate.x program   !
    ! cooFile        name of the currently active structure file         !
    ! keyword        the last keyword read from the input file           !
    !                                                                    !
    ! do_debug       if .true., additional files containing debugging    !
    !                info will be created                                !
    !                                                                    !
    ! u_*            file units                                          !
    !--------------------------------------------------------------------!
  
    type(InputData)                                :: inp
  
    type(Setup),       dimension(:),   allocatable :: stp
    double precision                               :: r_min, r_max
    type(TrnSet)                                   :: ts
  
    integer                                        :: nnb_max, nnb
    double precision,  dimension(:,:), allocatable :: nbcoo
    double precision,  dimension(:),   allocatable :: nbdist
    integer,  dimension(:),   allocatable :: nblist
    integer,           dimension(:),   allocatable :: nbtype
  
    double precision, dimension(:),     allocatable :: sfval
    double precision, dimension(:,:),   allocatable :: sfderiv_i
    double precision, dimension(:,:,:), allocatable :: sfderiv_j
  
    double precision                               :: E_coh
    integer                                        :: ifile
    double precision                               :: nFiles_inv

    character(len=*),intent(in)                    :: inFile
    integer,intent(in)                             :: ionum
  
    !character(len=PATHLEN)                         :: inFile
    character(len=PATHLEN)                         :: cooFile
    character(len=LINELEN)                         :: keyword
  
    integer                                        :: itype1
    integer                                        :: itype, iatom
  
    integer                                        :: iline
    character(len=1024)                            :: line
  
    integer                                        :: u_in, u_tng
    logical                                        :: do_debug = .false.
    integer                                        :: u_dbg, idbg
  
    integer                                        :: u_pyo, max_nnb_trainset, pyo_forces_struc
    integer, allocatable                           :: struc_write_force(:)
  
  
    integer :: i, j, l
    logical :: lexist
    ! timing registers
    integer, parameter :: R_GEO = 1, R_NBL = 2, R_SF = 3

    character(len=1024)::outfilename_ascii
    character(len=1024)::outFileName
    logical ::to_bin,to_ascii

    logical :: stopnow
    character(len=262144)                           :: longline
  
    integer rename, status
    integer::N_do_forces
  
  !-------------------------- initialization --------------------------!

    if (ppMaster) then
      call system( 'rm -f ts.* 2> /dev/null' )
    end if
 
    call initialize_MPI()
  
    if (ppMaster) inp = read_InpGenerate(inFile)
    call pp_bcast_InputData(inp)

    allocate(stp(inp%nTypes))
    call load_symmfunc_setups_MPI(inp, stp)
  
    ! call parse_input(inFile)
  
    if (inp%do_timing .and. ppMaster) then
       u_tng = io_unit()
       call tng_init(unit=u_tng, file='generate.time', registers=3)
       write(ionum,*) 'Timing info will be written to: ', 'generate.time'
       write(ionum,*)
    end if
    if (do_debug .and. ppMaster) then
       u_dbg = io_unit()
       open(u_dbg, file='generate.debug', status='replace', action='write')
    end if
  
    ! get interaction range and max. number of atoms within range
    call stp_get_range(inp%nTypes, stp, r_min, r_max)
    nnb_max = lcl_nmax_nbdist(r_min, r_max)
    allocate(nbcoo(3,nnb_max), nbdist(nnb_max), nblist(nnb_max),nbtype(nnb_max))
  
    ! initialize workspace for structural fingerprint basis:
    call stp_init(inp%nTypes, stp, nnb_max)
    if (inp%do_timing) call tng_timing('Structural fingerprint basis initialized.')
  
    ! allocate workspace for basis function evaluation:
    allocate(sfval(nsf_max), sfderiv_i(3,nsf_max), sfderiv_j(3,nsf_max,nnb_max))
    sfval(:) = 0.0d0
    sfderiv_i(:,:) = 0.0d0
    sfderiv_j(:,:,:) = 0.0d0

    if (ppMaster) then
      call aeio_header('Generation of training set started')
      write(ionum,*)
   
      write(ionum,*) 'Number of atom types  : ', trim(io_adjustl(inp%nTypes))
      write(*,'(1x,"types                 : ")', advance='no')
      do itype = 1, inp%nTypes
         if (mod(itype,7) == 0) write(*,'(29x)')
         write(*,'(A5,1x)', advance='no') inp%typeName(itype)
      end do
      write(ionum,*)
      write(ionum,*) "Number of structures  : ", trim(io_adjustl(inp%nStrucs))
      write(ionum,*)
    end if
  
    !-------------- write basis function settings to stdout -------------!
    if (ppMaster) then
      call aeio_header("Structural fingerprint basis set-up")
      write(ionum,*)
  
      do itype1 = 1, inp%nTypes
       call stp_print_info(stp(itype1))
      end do
   end if
  
    !----------- write training set header to the output file -----------!
  
    ts = new_TrnSet_MPI(inp%nTypes, inp%typeName, inp%atomicEnergy, &
                    inp%nStrucs, trim(inp%outFileName))

    outfilename_ascii = trim(inp%outFileName)//'.ascii'                    
  
    if (inp%do_timing .and. ppMaster) call tng_timing('Training set file started.')
  
    !--------------- write pytorch force training header ----------------!

    call pp_barrier()

    if (inp%pyo_forces) then
      u_pyo = io_unit()
      if (ppMaster) then
         !call pyo_write_init(u_pyo, trim(inp%outFileName)//".forces.tmp")
         call pyo_write_init(u_pyo, "ts.forces.tmp")
         call pyo_write_header_info(u_pyo, inp%nStrucs)
         call pyo_select_force_structures(inp%nStrucs, inp%pyo_forces_percent, struc_write_force)
         N_do_forces = ubound(struc_write_force,1)
      end if
      call pp_bcast(N_do_forces)
      if (ppMaster) then
      else
         allocate(struc_write_force(N_do_forces))
      end if
      do i=1,N_do_forces
         call pp_bcast(struc_write_force(i))
      end do
      max_nnb_trainset = 0
    endif
  
    !------------------ iterate over coordinates files ------------------!
  
    if (ppMaster) then
      call aeio_header("Adding structures to the training set")
      write(ionum,*)
    end if
  
    u_in = io_unit()
    open(u_in, file=inFile, status='old', action='read')
    rewind(u_in)
  
    iline = 0
    do
       ! forward until the FILES keyword:
       call aeio_readline(u_in, iline, line)
       read(line,*) keyword
       if (trim(keyword) == 'FILES') then
          read(u_in,*)
          exit
       end if
    end do
  
    if (ppMaster) then
    ! header for stdout
      write(*,'("#",A6,2x,A6,2x,A6,2x,A15,2x,A)') &
            'N', 'nAtoms', 'nTypes', 'E/atom', 'structure file (xsf)'
   
    end if

    stopnow = .false.

    call pp_barrier()



    nFiles_inv = 1.0d0/dble(inp%nStrucs)
    structures : do ifile = 1, inp%nStrucs
  
       if (inp%do_timing .and. ppMaster) call tng_timing('Structure: '// io_adjustl(ifile))
  
       call aeio_readline(u_in, iline, line)
       cooFile = trim(line)

       if ( mod(ifile-1,ppSize) .ne. ppRank ) cycle
  
       call geo_init(cooFile, 'xsf')
       if (inp%do_timing) call tng_timing3(register=R_GEO)
       if (.not. (hasForces .and. hasEnergy)) then
          write(0,*) ">>>", hasForces, hasEnergy
          write(0,*) "Error: incomplete output data in : ", trim(cooFile)
          !call finalize()
          stop
       end if
  
       if (nTypes > inp%nTypes) then
          if (ppMaster) then
            write(ionum,*) 'Skipping ', trim(adjustl(cooFile)), &
                     ': too many atomic species'
          end if
          call geo_final()
          cycle structures
       end if
  
       if (abs(cohesiveEnergy) /= 0.0d0) then
          E_coh = cohesiveEnergy
       else
          ! if only the total energy is available, we have to calculate
          ! the cohesive energy at this point
          E_coh = totalEnergy
          do iatom = 1, nAtoms
             itype1 = atomType(iatom)
             itype1 = geo_type_conv(itype1, nTypes, atomTypeName, &
                                    inp%nTypes, inp%typeName)
             E_coh = E_coh - inp%atomicEnergy(itype1)
          end do
       end if
  
       if (ppMaster) then
         write(*,'(1x,I6,2x,I6,2x,I6,2x,ES15.8,2x,A)') &
            ifile, nAtoms, nTypes, E_coh/dble(nAtoms), &
            trim(adjustl(cooFile))
       end if
            
       call lcl_init(r_min, r_max, latticeVec, nAtoms, atomType, cooLatt, pbc)
       if (inp%do_timing .and. ppMaster) call tng_timing3(register=R_NBL)
  
       ! write structure info (atoms, types, energy) to training set file:
       call ts_write_structure_info_MPI(ts, cooFile, nAtoms, nTypes, E_coh,ifile)
       
       ! write structure info (atoms, types) to pytorch forces file:
       if (inp%pyo_forces ) then
         pyo_forces_struc = struc_write_force(ifile)
         call pyo_write_structure_info_MPI(cooFile, nAtoms, nTypes, pyo_forces_struc,ifile)
       end if
  
       atoms : do iatom = 1, nAtoms
  
          ! determine the training atom type of atom `iatom' in global
          ! index terms
          itype1 = atomType(iatom)
          itype1 = geo_type_conv(itype1, nTypes, atomTypeName, &
                                 inp%nTypes, inp%typeName)
  
          ! assert that atom type is included in the set-ups:
          if (itype1 == 0) then
             write(0,*) "Error: not a valid structure    : ", trim(cooFile)
             write(0,*) "       Additional species found."
             !call finalize()
             stop
          end if
  
          ! write atom info (species, forces) to training set file:
          call ts_write_atom_info_MPI(ts, itype1, cooCart(iatom), forCart(iatom),ifile)
  
          ! get all atoms within cut-off:
          nnb = nnb_max
          call lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, r_cut=r_max, nblist=nblist, nbtype=nbtype)
  
          if (inp%do_timing .and. ppMaster) call tng_timing3(register=R_NBL)
          !write(*,'(1x,I6,2x,A2,2x,I6)') &
          !     iatom, trim(atomTypeName(atomType(iatom))), nnb
  
          ! convert atom types to global index:
          do i = 1, nnb
             nbtype(i) = geo_type_conv(nbtype(i), nTypes, atomTypeName, &
                                      inp%nTypes, inp%typeName)
             if (nbtype(i) == 0) then
                write(0,*) "Error: atom type not found in setup."
                call finalize_MPI()
                stop
             end if
          end do
  
          ! evaluate the structural fingerprint basis function set-up:
          call stp_eval(itype1, cooCart(iatom), nnb, nbcoo, nbtype, &
                        stp(itype1), sfval=sfval, sfderiv_i=sfderiv_i, &
                        sfderiv_j=sfderiv_j)
  
          if (do_debug .and. ppMaster) then
             do idbg = 1, stp(itype1)%nsf
                write(u_dbg,'(1x,ES15.8,1x)', advance='no') sfval(idbg)
             end do
             write(u_dbg,*)
          end if
  
          if (inp%do_timing) call tng_timing3(register=R_SF)
  
          ! write basis function values and derivatives
          ! to the training set file:
          call ts_write_sf_info_MPI(ts, stp(itype1)%nsf, sfval(1:stp(itype1)%nsf),ifile)
  
          ! write basis function derivatives and neighbor list to pytorch output
          if (inp%pyo_forces .and. pyo_forces_struc==1) then
            max_nnb_trainset = max(max_nnb_trainset, nnb)
            call pyo_write_atom_sf_info_MPI(itype1, nnb, stp(itype1)%nsf, nblist(:nnb), &
                                        sfderiv_i(1:3,1:stp(itype1)%nsf), &
                                        sfderiv_j(1:3,1:stp(itype1)%nsf,1:nnb) ,ifile)
          end if
  
       end do atoms
  
       if (inp%do_timing) then
          call tng_dump(R_GEO, 'time spent reading geometries (so far)')
          call tng_dump(R_NBL, 'time spent in the neighbor list (so far)')
          call tng_dump(R_SF,  'time spent evaluating structural fingerprints (so far)')
       end if
  
       call lcl_final()
       call geo_final()
  
    end do structures
    if (ppMaster) write(ionum,*)
  
    if (inp%do_timing) then
       call tng_timing('Loop over structures done.')
       call tng_dump(R_GEO, 'total time spent reading geometries')
       call tng_dump(R_NBL, 'total time spent in the neighbor list')
       call tng_dump(R_SF,  'total time spent evaluating structural fingerprints')
    end if
  
    !--------- save basis function setups with final statistics ---------!
  
    if (ppMaster) call ts_print_info(ts)
  
    !----------------------------- finalize -----------------------------!
  
    deallocate(nbcoo, nbdist, nblist, nbtype)
    close(u_in)
    
    if (inp%pyo_forces .and. ppMaster) then
      call pyo_write_final(u_pyo, max_nnb_trainset)
    end if


    call pp_barrier()

    if (ppMaster) then
     do ifile = 1, inp%nStrucs
        call ts_loadandwrite_structure_info_MPI(ifile)
     end do
    end if
  
    call ts_parallel_footer_MPI(ts, stp)

    if (inp%pyo_forces .and. ppMaster) then
      call pyo_loadandwrite("ts.forces.tmp")
      !call pyo_loadandwrite(trim(inp%outFileName)//".forces.tmp")
    end if



    call close_TrnSet_MPI(ts, stp=stp(1:inp%nTypes))


    l = 5000
    call pp_barrier()
    if (ppMaster) then
     status = rename( "ts.all",adjustl(trim(inp%outFileName))) !rename the file
        do ifile = 1, inp%nStrucs
           open(10001, file='ts.'//io_adjustl(ifile))
           close(10001, status="delete")
        end do
    end if

    if (inp%pyo_forces .and. ppMaster) then
      status = rename( "pyo.all",trim(inp%outFileName)//".forces") !rename the file
      do ifile = 1, inp%nStrucs
         open(10001, file='ts.force.'//io_adjustl(ifile))
         close(10001, status="delete")
      end do
    end if
  
    call pp_barrier()
    call finalize_MPI()


    to_bin = .false.
    to_ascii = .true.
    outFileName = inp%outFileName
    call trainbin2ascii_subroutine(trim(outFileName), trim(outfilename_ascii),to_bin, to_ascii)
  

    return

    call finalize()


  
  contains !=============================================================!
  
  
  subroutine initialize_MPI()

    implicit none

!    character(len=*), intent(out) :: inFile

    integer :: nargs
    logical :: fexists
    logical :: stopnow

    call pp_init()

!    stopnow = .false.

    if (ppMaster) then

       call aeio_header("generate.x - training set generation", char='=')
       write(ionum,*)

       call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')

!       nargs = command_argument_count()
!       if (nargs < 1) then
!          write(0,*) "Error: No input file provided."
!          call print_usage()
!          stopnow = .true.
!       end if

!       call get_command_argument(1, value=inFile)
!       inquire(file=trim(inFile), exist=fexists)
!       if (.not. fexists) then
!          write(0,*) "Error: File not found: ", trim(inFile)
!          stopnow = .true.
!       end if

   end if
!   stopnow = .true.

!    call pp_bcast(stopnow)
!    if (stopnow) then
!       call finalize_MPI()
!       stop
!    end if

   ! call pp_bcast(inFile)

  end subroutine initialize_MPI

    !--------------------------------------------------------------------!
  
  subroutine finalize_MPI()

    implicit none

    integer :: itype

    if (allocated(sfval)) then
       deallocate(sfval, sfderiv_i, sfderiv_j)
    end if

    if (allocated(stp)) then
       call stp_final(inp%nTypes, stp)
       do itype = 1, inp%nTypes
          call del_Setup(stp(itype))
       end do
       deallocate(stp, inp%typeName, inp%atomicEnergy)
    end if

    if (ts%init) call close_TrnSet_MPI(ts)

    if (allocated(nbcoo)) deallocate(nbcoo, nbdist)

    if (inp%do_timing .and. ppMaster ) call tng_final()
    if (do_debug .and. ppMaster )  close(u_dbg)

    if (ppMaster) then
       call aeio_header(aeio_timestamp(), char=' ')
       call aeio_header("Training set generation done.", char='=')
    end if

    call pp_final()

  end subroutine finalize_MPI

  
    !--------------------------------------------------------------------!
  
    subroutine print_usage()
  
      implicit none
  
      write(ionum,*)
      write(ionum,*) "generate.x -- Generate training sets for use with `train.x'"
      write(*,'(1x,70("-"))')
      write(ionum,*) 'Usage: generate.x <input-file>'
      write(ionum,*)
      write(ionum,*) 'See the documentation or the source code for a description of the '
      write(ionum,*) 'input file format.'
      write(ionum,*)
  
    end subroutine print_usage
  
    !--------------------------------------------------------------------!
  
  subroutine load_symmfunc_setups_MPI(inp, stp)

    implicit none

    type(InputData),           intent(in)  :: inp
    type(Setup), dimension(:), intent(out) :: stp
    logical :: stopnow
    integer :: i

    stopnow = .false.
    
   
    if(ppMaster) then 
      do i = 1, inp%nTypes
         stp(i) = read_Setup_parameters(inp%setupFile(i), inp%typeName(:))
         if (.not. (trim(stp(i)%atomtype) == trim(inp%typeName(i)))) then
            write(0,*) "Error: Inconsistent atom type in setup:"
            write(0,*) "       type expected : ", trim(inp%typeName(i))
            write(0,*) "       type found    : ", trim(stp(i)%atomtype)
            stopnow = .true.
         end if
      end do
   end if

   call pp_barrier()
   

   do i = 1, inp%nTypes
      call pp_bcast_Setup(stp(i))
   end do

   
   call pp_barrier()



    call pp_bcast(stopnow)
    if (stopnow) then
       call finalize_MPI()
       stop
    end if

  end subroutine load_symmfunc_setups_MPI

 end subroutine generate_subroutine_MPI



 
end module