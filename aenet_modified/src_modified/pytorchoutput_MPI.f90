module pytorchoutput_MPI
  use io,      only: io_adjustl, io_unit

  implicit none
  private
  save

  public :: pyo_write_structure_info_MPI,     &
            pyo_write_atom_sf_info_MPI,&
            pyo_loadandwrite       

contains



  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!

  subroutine pyo_write_structure_info_MPI(filename, natoms, ntypes, pyo_forces_struc,ifile)

      implicit none
  
      integer, intent(in)          :: natoms, ntypes, pyo_forces_struc
      character(len=*), intent(in) :: filename
      integer,          intent(in)    :: ifile
      integer,parameter :: pyo_unit = 1234

      open(unit = pyo_unit, action = "write", status = "replace", file = 'ts.force.'//io_adjustl(ifile), form = "unformatted")    

      write(pyo_unit) len_trim(filename)
      write(pyo_unit) trim(filename)
      write(pyo_unit) natoms, ntypes
      write(pyo_unit) pyo_forces_struc
      
  end subroutine pyo_write_structure_info_MPI

  !--------------------------------------------------------------------!

  subroutine pyo_write_atom_sf_info_MPI(itype, nnb, nsf, nblist, sfderiv_i, sfderiv_j,ifile)

      implicit none

      integer,                            intent(in) ::  nnb, nsf, itype
      integer,          dimension(:),     intent(in) :: nblist
      double precision, dimension(:,:),   intent(in) :: sfderiv_i
      double precision, dimension(:,:,:), intent(in) :: sfderiv_j
      integer,          intent(in)    :: ifile

      double precision                               :: sfderiv_j_aux(nnb, nsf, 3), sfderiv_i_aux(nsf,3)

      integer :: ineigh, isf, icoo
      integer,parameter :: pyo_unit = 1234

      !call pyo_write_init(u_pyo, trim(inp%outFileName)//".forces.tmp")
      !call pyo_write_header_info(u_pyo, inp%nStrucs)
      !call pyo_select_force_structures(inp%nStrucs, inp%pyo_forces_percent, struc_write_force)
     

      open(unit = pyo_unit, action = "write", file = 'ts.force.'//io_adjustl(ifile), form='unformatted', position='append')  

      sfderiv_j_aux = 0.0d0
      sfderiv_i_aux = 0.0d0
      do ineigh = 1, nnb
        do isf = 1, nsf
          sfderiv_j_aux(ineigh, isf, :) = sfderiv_j(:, isf, ineigh)
        enddo
      enddo
      do isf = 1, nsf
        sfderiv_i_aux(isf,:) = sfderiv_i(:,isf)
      enddo

      write(pyo_unit) itype
      write(pyo_unit) nsf, nnb
      write(pyo_unit) nblist(1:nnb)
      write(pyo_unit) sfderiv_i_aux(1:nsf,1:3)
      write(pyo_unit) sfderiv_j_aux(1:nnb,1:nsf,1:3)
      
      
  end subroutine pyo_write_atom_sf_info_MPI

 
  subroutine pyo_loadandwrite(headfile)
    implicit none
    character(len=*),intent(in)::headfile
    integer,parameter:: pyo_head = 902
    integer::ifile
    character(len=100)::filenameload
    integer::l
    integer::nAtoms,ntypes,pyo_forces_struc

    integer,parameter :: pyo_unitall = 545
    integer::nstrucs
    integer,parameter::pyo_file = 1092
    integer::iatom
    integer::itype,nsf,nnb,max_nnb_trainset

    double precision, allocatable,dimension(:,:):: sfderiv_i_aux
    double precision, allocatable,dimension(:,:,:):: sfderiv_j_aux
    integer,          allocatable,dimension(:) :: nblist

    open (pyo_unitall, file='pyo.all', status='replace',form='unformatted')
    open(unit = pyo_head,  status = "old", file = headfile, form = "unformatted")

    

    read(pyo_head) nstrucs
    write(pyo_unitall) nstrucs

    read(pyo_head) max_nnb_trainset

    do ifile = 1, nStrucs
      ! call pyo_write_structure_info_MPI(cooFile, nAtoms, nTypes, pyo_forces_struc,ifile)
      open(unit = pyo_file,  status = "old", file = 'ts.force.'//io_adjustl(ifile), form = "unformatted")
      read(pyo_file) l
      filenameload = " "
      l = min(l,len(filenameload))
      read(pyo_file) filenameload(1:l)

      read(pyo_file) natoms, ntypes
      read(pyo_file) pyo_forces_struc

      write(pyo_unitall) len_trim(filenameload)
      write(pyo_unitall) trim(filenameload)
      write(pyo_unitall) natoms, ntypes
      write(pyo_unitall) pyo_forces_struc

      atoms : do iatom = 1, nAtoms
        if (pyo_forces_struc==1) then
          !max_nnb_trainset = max(max_nnb_trainset, nnb)
          !call pyo_write_atom_sf_info_MPI(itype1, nnb, stp(itype1)%nsf, nblist(:nnb), &
          !                            sfderiv_i(1:3,1:stp(itype1)%nsf), &
          !                            sfderiv_j(1:3,1:stp(itype1)%nsf,1:nnb) ,ifile)
          read(pyo_file) itype
          read(pyo_file) nsf, nnb
          allocate(nblist(1:nnb))

          read(pyo_file) nblist(1:nnb)
          allocate(sfderiv_i_aux(1:nsf,1:3))
          read(pyo_file) sfderiv_i_aux(1:nsf,1:3)
          allocate(sfderiv_j_aux(1:nnb,1:nsf,1:3))
          read(pyo_file) sfderiv_j_aux(1:nnb,1:nsf,1:3)

          write(pyo_unitall) itype
          write(pyo_unitall) nsf, nnb
          write(pyo_unitall) nblist(1:nnb)
          write(pyo_unitall) sfderiv_i_aux(1:nsf,1:3)
          write(pyo_unitall) sfderiv_j_aux(1:nnb,1:nsf,1:3)
          deallocate(nblist,sfderiv_i_aux,sfderiv_j_aux)
        end if
      end do atoms
      close(pyo_file)
    end do

    
    write(pyo_unitall) max_nnb_trainset

    close(pyo_unitall)
    close(pyo_head)


  end subroutine

end module pytorchoutput_MPI