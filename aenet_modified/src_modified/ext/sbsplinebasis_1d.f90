!-----------------------------------------------------------------------
! sfbasis.f90 - Basis for structural fingerprints of atomic environments
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
! 2015-12-27 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
!3rd -order Bspline
module sfbspline_1d

  use io, only: io_unit

  use chebyshev, only: chebyshev_polynomial, &
  chebyshev_polynomial_d1

  use bspline,only:make_knotsvector,bspline_basis_functions,&
        bspline_basis_functions_deriv,d => bspline_order

  implicit none
  private
  save

  public :: new_SBP1DBasis,            &
            del_SBP1DBasis,            &
            save_SBP1DBasis,           &
            load_SBP1DBasis,           &
            save_SBP1DBasis_ASCII,     &
            load_SBP1DBasis_ASCII,     &
            sbspline1d_print_info,         &
            sbspline1d_set_typeid,         &
            sbspline1d_set_typespin,       &
            sbspline1d_eval,               &
            sbspline1d_reconstruct_radial, &
            sbspline1d_reconstruct_angular

  type, public :: BsplineBasis_1d
     logical                                       :: initialized = .false.
     integer                                       :: r_points
     double precision                              :: r_Rc
     integer                                       :: r_N
     integer                                       :: a_points
     double precision                              :: a_Rc
     integer                                       :: a_N
     integer                                       :: r_i1, r_f1
     integer                                       :: r_i2, r_f2
     integer                                       :: a_i1, a_f1
     integer                                       :: a_i2, a_f2
     integer                                       :: N
     integer                                       :: num_types
     logical                                       :: multi
     character(len=2), dimension(:),   allocatable :: atom_types
     integer,          dimension(:),   allocatable :: typeid
     double precision, dimension(:),   allocatable :: typespin
     integer                                       :: num_values
     double precision, dimension(:),   allocatable :: r_knots
     double precision, dimension(:),   allocatable :: a_knots
  end type BsplineBasis_1d


  double precision, parameter, private :: PI     = 3.14159265358979d0
  double precision, parameter, private :: PI_INV = 1.0d0/PI
  double precision, parameter, private :: PI2    = 2.0d0*PI
  double precision, parameter, private :: EPS    = 1.0d-12

  

contains

  !--------------------------------------------------------------------!
  !            create a new Structural Fingerprint basis               !
  !--------------------------------------------------------------------!

  function new_SBP1DBasis(num_types, atom_types, radial_points, &
                       angular_points, radial_Rc, angular_Rc&
                       ) result(sfb)
    ! Arguments:
    !   num_types       number of atomic species
    !   atom_types(i)   i-th atomic species (2 characters)
    !   radial_points    expansion order for the radial basis
    !   angular_points   expansion order for the angular basis
    !   radial_Rc       cutoff radius for radial basis
    !   angular_Rc      cutoff radius for angular basis
    !
    ! Returns:
    !   sfb             allocated instance of BsplineBasis_1d

    implicit none

    integer,                                intent(in) :: num_types
    character(len=*), dimension(num_types), intent(in) :: atom_types
    integer,                                intent(in) :: radial_points
    integer,                                intent(in) :: angular_points
    double precision,                       intent(in) :: radial_Rc
    double precision,                       intent(in) :: angular_Rc
    type(BsplineBasis_1d)                             :: sfb

    integer :: i, s

    sfb%num_types = num_types
    sfb%r_points = radial_points
    sfb%a_points = angular_points
    sfb%r_Rc = radial_Rc
    sfb%a_Rc = angular_Rc

    sfb%r_N = sfb%r_points +2*d -d -1!+ 1

    allocate(sfb%r_knots(sfb%r_points +2*d))
    call make_knotsvector(sfb%r_knots,d,sfb%r_points ,sfb%r_Rc)
    
    allocate(sfb%a_knots(sfb%a_points +2*d))
    call make_knotsvector(sfb%a_knots,d,sfb%a_points ,2*sfb%a_Rc)
!    call make_knotsvector(sfb%a_knots,d,sfb%a_points ,sfb%a_Rc)
    
    sfb%a_N = sfb%a_points +d -1
    !sfb%a_N = sfb%a_points +1!+2*d !chebyshev for debug
    sfb%num_values = max(sfb%r_N, sfb%a_N)
    sfb%N = sfb%r_N + sfb%a_N
    sfb%r_i1 = 1
    sfb%r_f1 = sfb%r_i1 + sfb%r_N - 1
    sfb%a_i1 = sfb%r_f1 + 1
    sfb%a_f1 = sfb%a_i1 + sfb%a_N - 1
    sfb%r_i2 = sfb%a_f1 + 1
    sfb%r_f2 = sfb%r_i2 + sfb%r_N - 1
    sfb%a_i2 = sfb%r_f2 + 1
    sfb%a_f2 = sfb%a_i2 + sfb%a_N - 1
    if (sfb%num_types > 1) then
       sfb%multi = .true.
       sfb%N = 2*sfb%N
    else
       sfb%multi = .false.
    end if

    allocate(sfb%atom_types(num_types),      &
             sfb%typeid(num_types),          &
             sfb%typespin(num_types))
    sfb%atom_types = atom_types

    do i = 1, num_types
       sfb%typeid(i) = i
    end do

    s = -num_types/2
    do i = 1, num_types
       if ((s == 0) .and. (mod(num_types, 2) == 0)) s = s + 1
       sfb%typespin(i) = dble(s)
       s = s + 1
    end do

    sfb%initialized = .true.

  end function new_SBP1DBasis

  !--------------------------------------------------------------------!
  !         delete (deallocate) a Structural Fingerprint basis         !
  !--------------------------------------------------------------------!

  subroutine del_SBP1DBasis(sfb)

    implicit none

    type(BsplineBasis_1d), intent(inout) :: sfb

    call sbspline1d_assert_init(sfb)

    deallocate(sfb%atom_types, sfb%typeid, sfb%typespin)
    sfb%initialized = .false.

  end subroutine del_SBP1DBasis

  !--------------------------------------------------------------------!
  !                    print information to stdout                     !
  !--------------------------------------------------------------------!

  subroutine sbspline1d_print_info(sfb)

    implicit none

    type(BsplineBasis_1d), intent(in) :: sfb
    character(len=1024) :: frmt

    write(*,'(" Radial cutoff : ",F7.3)') sfb%r_Rc
    write(*,'(" Angular cutoff: ",F7.3)') sfb%a_Rc
    write(*,'(" Radial points  : ",I3)') sfb%r_points
    write(*,'(" Angular points : ",I3)') sfb%a_points
    
    write(*,'(" Atom types    : ")', advance='no')
    write(frmt, *) sfb%num_types
    frmt = '(' // trim(adjustl(frmt))  // '(A2,1x))'
    write(*,frmt) sfb%atom_types
    write(*,'(" Total number of basis functions: ", I3)') sfb%N

  end subroutine sbspline1d_print_info

  !============================ properties ============================!


  subroutine sbspline1d_set_typeid(sfb, typeid)

    implicit none

    type(BsplineBasis_1d), intent(inout) :: sfb
    integer, dimension(:),  intent(in)    :: typeid

    call sbspline1d_assert_init(sfb)

    if (size(typeid) < sfb%num_types) then
       write(0,*) "Error: incompatible type ID list in `sbspline1d_set_typeid'."
       stop
    end if

    sfb%typeid(1:sfb%num_types) = typeid(1:sfb%num_types)

  end subroutine sbspline1d_set_typeid

  !--------------------------------------------------------------------!

  subroutine sbspline1d_set_typespin(sfb, typespin)

    implicit none

    type(BsplineBasis_1d),         intent(inout) :: sfb
    double precision, dimension(:), intent(in)    :: typespin

    call sbspline1d_assert_init(sfb)

    if (size(typespin) < sfb%num_types) then
       write(0,*) "Error: incompatible type ID list in `sbspline1d_set_typespin'."
       stop
    end if

    sfb%typespin(1:sfb%num_types) = typespin(1:sfb%num_types)

  end subroutine sbspline1d_set_typespin


  !=============================== I/O ================================!


  !--------------------------------------------------------------------!
  !    Read/write Structural Fingerprint Basis from/to file or unit    !
  !                                                                    !
  !    The *_ASCII procedures read/write plain text files; the other   !
  !    procedures just use binary I/O.                                 !
  !--------------------------------------------------------------------!

  subroutine save_SBP1DBasis(sfb, file, unit)

    implicit none

    type(BsplineBasis_1d),     intent(in) :: sfb
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: u

    call sbspline1d_assert_init(sfb)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write', &
            form='unformatted')
    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `save_SBP1DBasis'."
       return
    end if

    write(u) sfb%r_points
    write(u) sfb%a_points
    write(u) sfb%r_Rc
    write(u) sfb%a_Rc
    write(u) sfb%r_N, sfb%a_N, sfb%N
    write(u) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    write(u) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    write(u) sfb%num_values
    write(u) sfb%num_types
    write(u) sfb%atom_types(:)
    write(u) sfb%typeid(:)
    write(u) sfb%typespin(:)

    if (.not. present(unit)) close(u)

  end subroutine save_SBP1DBasis

  !--------------------------------------------------------------------!

  function load_SBP1DBasis(file, unit) result(sfb)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(BsplineBasis_1d)                 :: sfb

    integer :: u

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), action='read', form='unformatted')
    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `load_SBP1DBasis'."
       return
    end if

    read(u) sfb%r_points
    read(u) sfb%a_points
    read(u) sfb%r_Rc
    read(u) sfb%a_Rc
    read(u) sfb%r_N, sfb%a_N, sfb%N
    read(u) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    read(u) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    read(u) sfb%num_values
    read(u) sfb%num_types
    allocate(sfb%atom_types(sfb%num_types),  &
             sfb%typeid(sfb%num_types),      &
             sfb%typespin(sfb%num_types))
    read(u) sfb%atom_types(:)
    read(u) sfb%typeid(:)
    read(u) sfb%typespin(:)

    if (sfb%num_types > 1) then
       sfb%multi = .true.
    else
       sfb%multi = .false.
    end if
    sfb%initialized = .true.

    if (.not. present(unit)) close(u)

  end function load_SBP1DBasis

  !--------------------------------------------------------------------!

  subroutine save_SBP1DBasis_ASCII(sfb, file, unit)

    implicit none

    type(BsplineBasis_1d),     intent(in) :: sfb
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    character(len=*), parameter :: DFRMT = '(4(1x,ES24.17))'
    character(len=*), parameter :: IFRMT = '(4(1x,I17))'
    character(len=*), parameter :: AFRMT = '(4(1x,A))'

    integer :: u, i

    call sbspline1d_assert_init(sfb)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write')
    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `save_SBP1DBasis_ASCII'."
       return
    end if

    write(u,*) sfb%r_points
    write(u,*) sfb%a_points
    write(u,*) sfb%r_Rc
    write(u,*) sfb%a_Rc
    write(u,*) sfb%r_N, sfb%a_N, sfb%N
    write(u,*) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    write(u,*) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    write(u,*) sfb%num_values
    write(u,*) sfb%num_types
    write(u,AFRMT) (sfb%atom_types(i), i=1,sfb%num_types)
    write(u,IFRMT) (sfb%typeid(i), i=1,sfb%num_types)
    write(u,DFRMT) (sfb%typespin(i), i=1,sfb%num_types)

    if (.not. present(unit)) close(u)

  end subroutine save_SBP1DBasis_ASCII

  !--------------------------------------------------------------------!

  function load_SBP1DBasis_ASCII(file, unit) result(sfb)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(BsplineBasis_1d)                 :: sfb

    character(len=*), parameter :: DFRMT = '(4(1x,ES24.17))'
    character(len=*), parameter :: IFRMT = '(4(1x,I17))'
    character(len=*), parameter :: AFRMT = '(4(1x,A))'

    integer :: u

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), action='read')
    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `load_SBP1DBasis_ASCII'."
       return
    end if

    read(u,*) sfb%r_points
    read(u,*) sfb%a_points
    read(u,*) sfb%r_Rc
    read(u,*) sfb%a_Rc
    read(u,*) sfb%r_N, sfb%a_N, sfb%N
    read(u,*) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    read(u,*) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    read(u,*) sfb%num_values
    read(u,*) sfb%num_types
    allocate(sfb%atom_types(sfb%num_types),  &
             sfb%typeid(sfb%num_types),      &
             sfb%typespin(sfb%num_types))
    read(u,AFRMT) sfb%atom_types(:)
    read(u,IFRMT) sfb%typeid(:)
    read(u,DFRMT) sfb%typespin(:)

    if (sfb%num_types > 1) then
       sfb%multi = .true.
    else
       sfb%multi = .false.
    end if
    sfb%initialized = .true.

    if (.not. present(unit)) close(u)

  end function load_SBP1DBasis_ASCII

  subroutine cleararray(nv,nat,deriv0,deriv1)
   implicit none
   integer,intent(in) :: nv,nat
   double precision, dimension(3,nv),     intent(out)   :: deriv0
   double precision, dimension(3,nv,nat), intent(out)   :: deriv1
   deriv0(:,:) = 0.0d0
   deriv1(:,:,:) = 0.0d0
  end subroutine

  subroutine update_deriv_r(i1,i2,j,N,sbspline_values,sbspline_deriv_i,sbspline_deriv_j,&
      values,deriv0,deriv1)
   implicit none
   integer::i1,i2,j,k,N
   double precision, dimension(:),     intent(in)   ::sbspline_values
   double precision, dimension(:,:),     intent(in)   :: sbspline_deriv_i
   double precision, dimension(:,:), intent(in)   :: sbspline_deriv_j
   double precision, dimension(:),     intent(inout)   ::values
   double precision, dimension(:,:),     intent(inout)   :: deriv0
   double precision, dimension(:,:,:), intent(inout)   :: deriv1


   values(i1:i2) = values(i1:i2) + sbspline_values(1:N)
   deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) + sbspline_deriv_i(1:3, 1:N)
   deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) + sbspline_deriv_j(1:3, 1:N)

   return
  end subroutine

  subroutine update_deriv_r_multi(i1,i2,j,N,sbspline_values,sbspline_deriv_i,sbspline_deriv_j,&
      values,deriv0,deriv1,s_j,do_deriv)
   implicit none
   integer::i1,i2,j,k,N
   double precision, dimension(:),     intent(in)   ::sbspline_values
   double precision, dimension(:,:),     intent(in)   :: sbspline_deriv_i
   double precision, dimension(:,:), intent(in)   :: sbspline_deriv_j
   double precision, dimension(:),     intent(inout)   ::values
   double precision, dimension(:,:),     intent(inout)   :: deriv0
   double precision, dimension(:,:,:), intent(inout)   :: deriv1
   double precision,intent(in):: s_j
   logical,intent(in)::do_deriv

   values(i1:i2) = values(i1:i2) + s_j*sbspline_values(1:N)
   if (do_deriv) then
      deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) &
                         + s_j*sbspline_deriv_i(1:3, 1:N)
      deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) &
                            + s_j*sbspline_deriv_j(1:3, 1:N)
   end if

   return
  end subroutine

  subroutine update_deriv(i1,i2,j,k,N,sbspline_values,sbspline_deriv_i,sbspline_deriv_j,sbspline_deriv_k,&
      values,deriv0,deriv1)
   implicit none
   integer::i1,i2,j,k,N
   double precision, dimension(:),     intent(in)   ::sbspline_values
   double precision, dimension(:,:),     intent(in)   :: sbspline_deriv_i
   double precision, dimension(:,:), intent(in)   :: sbspline_deriv_j
   double precision, dimension(:,:), intent(in)   :: sbspline_deriv_k
   double precision, dimension(:),     intent(inout)   ::values
   double precision, dimension(:,:),     intent(inout)   :: deriv0
   double precision, dimension(:,:,:), intent(inout)   :: deriv1


   values(i1:i2) = values(i1:i2) + sbspline_values(1:N)
   deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) + sbspline_deriv_i(1:3, 1:N)
   deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) + sbspline_deriv_j(1:3, 1:N)
   deriv1(1:3, i1:i2, k) = deriv1(1:3, i1:i2, k) + sbspline_deriv_k(1:3, 1:N)

   return
  end subroutine

  subroutine update_deriv_multi(i1,i2,j,k,N,sbspline_values,sbspline_deriv_i,sbspline_deriv_j,sbspline_deriv_k,&
      values,deriv0,deriv1,do_deriv,s_j,s_k)
   implicit none
   integer::i1,i2,j,k,N
   double precision, dimension(:),     intent(in)   ::sbspline_values
   double precision, dimension(:,:),     intent(in)   :: sbspline_deriv_i
   double precision, dimension(:,:), intent(in)   :: sbspline_deriv_j
   double precision, dimension(:,:), intent(in)   :: sbspline_deriv_k
   double precision, dimension(:),     intent(inout)   ::values
   double precision, dimension(:,:),     intent(inout)   :: deriv0
   double precision, dimension(:,:,:), intent(inout)   :: deriv1
   logical,intent(in)::do_deriv
   double precision,intent(in) ::s_j,s_k

   values(i1:i2) = values(i1:i2) + s_j*s_k*sbspline_values(1:N)
   if (do_deriv) then
      deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) &
                         + s_j*s_k*sbspline_deriv_i(1:3, 1:N)
      deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) &
                            + s_j*s_k*sbspline_deriv_j(1:3, 1:N)
      deriv1(1:3, i1:i2, k) = deriv1(1:3, i1:i2, k) &
                            + s_j*s_k*sbspline_deriv_k(1:3, 1:N)
   end if
   

   return
  end subroutine


  !========================= basis evaluation =========================!
  subroutine sbspline1d_eval(sfb, itype0, coo0, nat, itype1, coo1, nv, &
                      values, deriv0, deriv1)

    implicit none
    type(BsplineBasis_1d),                          intent(inout) :: sfb
    !type(FingerprintBasis),                          intent(inout) :: sfb
    integer,                                         intent(in)    :: itype0
    double precision, dimension(3),                  intent(in)    :: coo0
    integer,                                         intent(in)    :: nat
    integer,          dimension(nat),                intent(in)    :: itype1
    double precision, dimension(3,nat),              intent(in)    :: coo1
    integer,                                         intent(in)    :: nv
    double precision, dimension(nv),                 intent(out)   :: values
    double precision, dimension(3,nv),     optional, intent(out)   :: deriv0
    double precision, dimension(3,nv,nat), optional, intent(out)   :: deriv1

    double precision, dimension(sfb%num_values)   :: sfbspline1d_values
    double precision, dimension(:,:), allocatable :: sfbspline1d_deriv_i
    double precision, dimension(:,:), allocatable :: sfbspline1d_deriv_j
    double precision, dimension(:,:), allocatable :: sfbspline1d_deriv_k

    logical                        :: do_deriv
    double precision, dimension(3) :: R_ij, R_ik,R_jk
    double precision               :: d_ij, d_ik,d_jk
    double precision               :: cos_ijk
    double precision               :: s_j, s_k
    integer                        :: j, k, i1, i2, N

    !if (present(deriv0) .and. present(deriv1)) then
    !  call sfb_eval_deriv(sfb, itype0, coo0, nat, itype1, coo1, nv, &
    !!                  values, deriv0, deriv1)
    !  return
    !end if

    call sbspline1d_assert_init(sfb)

    if (nv /= sfb%N) then
       write(0,*) "Error: wrong number of basis functions in `sfb_eval'."
       stop
    end if

    if (present(deriv0) .and. present(deriv1)) then
       do_deriv = .true.
       deriv0(:,:) = 0.0d0
       deriv1(:,:,:) = 0.0d0
       !call cleararray(nv,nat,deriv0,deriv1)
       allocate(sfbspline1d_deriv_i(3, sfb%num_values), &
                sfbspline1d_deriv_j(3, sfb%num_values), &
                sfbspline1d_deriv_k(3, sfb%num_values))
    else
       do_deriv = .false.
    end if

    values(1:sfb%N) = 0.0d0
    s_j = 1.0d0

    for_j : do j = 1, nat
       R_ij = coo1(1:3, j) - coo0(1:3)
       d_ij = sqrt(dot_product(R_ij, R_ij))
       if ((d_ij <= sfb%r_Rc) .and. (d_ij > EPS)) then

          ! evaluate radial basis functions
          i1 = sfb%r_i1
          i2 = sfb%r_f1
          N = sfb%r_N
          if (do_deriv) then
             call sbspline1d_radial(sfb, R_ij, d_ij, sfbspline1d_values,sfb%r_N,sfb%r_knots,&
                             deriv_i=sfbspline1d_deriv_i, deriv_j=sfbspline1d_deriv_j)
            !call update_deriv_r(i1,i2,j,N,sfb_values,sfb_deriv_i,sfb_deriv_j,&
            !   values,deriv0,deriv1)

             values(i1:i2) = values(i1:i2) + sfbspline1d_values(1:N)
             deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) + sfbspline1d_deriv_i(1:3, 1:N)
             deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) + sfbspline1d_deriv_j(1:3, 1:N)
          else
             call sbspline1d_radial(sfb, R_ij, d_ij, sfbspline1d_values, sfb%r_N,sfb%r_knots)
             values(i1:i2) = values(i1:i2) + sfbspline1d_values(1:N)
          end if

          ! redundant radial basis in case of multi-component systems
          i1 = sfb%r_i2
          i2 = sfb%r_f2
          N = sfb%r_N
          if (sfb%multi) then
            !write(*,*) "j",j,itype1(j)
            !write(*,*) sfb%typeid(itype1(j))
            !write(*,*) sfb%typespin(sfb%typeid(itype1(j)))
             s_j = sfb%typespin(sfb%typeid(itype1(j)))
             !call update_deriv_r_multi(i1,i2,j,N,sfb_values,&
             !  sfb_deriv_i,sfb_deriv_j,&
             !  values,deriv0,deriv1,s_j,do_deriv)
             values(i1:i2) = values(i1:i2) + s_j*sfbspline1d_values(1:N)
             if (do_deriv) then
                deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) &
                                   + s_j*sfbspline1d_deriv_i(1:3, 1:N)
                deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) &
                                      + s_j*sfbspline1d_deriv_j(1:3, 1:N)
             end if
          end if

       end if  ! within radial cutoff

       if (d_ij > sfb%a_Rc) cycle for_j
       


       for_k : do k = j+1, nat
          R_ik = coo1(1:3, k) - coo0(1:3)
          d_ik = sqrt(dot_product(R_ik, R_ik))
          if ((d_ik > sfb%a_Rc) .or. (d_ik < EPS)) cycle for_k
          !cos_ijk = dot_product(R_ij, R_ik)/(d_ij*d_ik)
          R_jk = coo1(1:3, j) - coo1(1:3,k)
          d_jk = sqrt(dot_product(R_jk, R_jk))

          ! evaluate angular basis functions
          i1 = sfb%a_i1
          i2 = sfb%a_f1
          N = sfb%a_N
          if (do_deriv) then
             call sbspline1d_angular(sfb, R_ij, R_ik, d_ij, d_ik, d_jk, &
                              sfbspline1d_values, (sfb%a_points+d-1),sfb%a_knots,deriv_i=sfbspline1d_deriv_i,      &
                              deriv_j=sfbspline1d_deriv_j, deriv_k=sfbspline1d_deriv_k)
             !call update_deriv(i1,i2,j,k,N,sfb_values,sfb_deriv_i,sfb_deriv_j,sfb_deriv_k,&
             !           values,deriv0,deriv1)
             values(i1:i2) = values(i1:i2) + sfbspline1d_values(1:N)
             deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) + sfbspline1d_deriv_i(1:3, 1:N)
             deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) + sfbspline1d_deriv_j(1:3, 1:N)
             deriv1(1:3, i1:i2, k) = deriv1(1:3, i1:i2, k) + sfbspline1d_deriv_k(1:3, 1:N)
          else
             call sbspline1d_angular(sfb, R_ij, R_ik, d_ij, d_ik, d_jk, sfbspline1d_values,&
                  (sfb%a_points+d-1),sfb%a_knots)
             values(i1:i2) = values(i1:i2) + sfbspline1d_values(1:N)
          end if

          ! redundant angular basis in case of multi-component systems
          i1 = sfb%a_i2
          i2 = sfb%a_f2
          N = sfb%a_N
          if (sfb%multi) then
             s_k = sfb%typespin(sfb%typeid(itype1(k)))
             !call update_deriv_multi(i1,i2,j,k,N,sfb_values,sfb_deriv_i,&
             !     sfb_deriv_j,sfb_deriv_k,&
             !     values,deriv0,deriv1,do_deriv,s_j,s_k)

             values(i1:i2) = values(i1:i2) + s_j*s_k*sfbspline1d_values(1:N)
             if (do_deriv) then
                deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) &
                                   + s_j*s_k*sfbspline1d_deriv_i(1:3, 1:N)
                deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) &
                                      + s_j*s_k*sfbspline1d_deriv_j(1:3, 1:N)
                deriv1(1:3, i1:i2, k) = deriv1(1:3, i1:i2, k) &
                                      + s_j*s_k*sfbspline1d_deriv_k(1:3, 1:N)
             end if
          end if
       end do for_k
    end do for_j

    if (do_deriv) deallocate(sfbspline1d_deriv_i, sfbspline1d_deriv_j, sfbspline1d_deriv_k)

  end subroutine sbspline1d_eval



  !======================= Basis Set Expansion ========================!

  !--------------------------------------------------------------------!
  ! reconstruct radial distribution function from basis set expansion  !
  !--------------------------------------------------------------------!

  subroutine sbspline1d_reconstruct_radial(sfb, coeff, nx, x, y)

    implicit none

    !------------------------------------------------------------------!
    ! sfb         Instance of BsplineBasis_1d                         !
    ! coeff(i)    Coefficient of the i-th basis function               !
    ! nx          Grid points for function evaluation                  !
    ! x(i)        x value of the i-th grid point (output)              !
    ! y(i)        y (function) value of the i-th grid point (output)   !
    !------------------------------------------------------------------!

    type(BsplineBasis_1d),               intent(in)  :: sfb
    double precision, dimension(sfb%r_N), intent(in)  :: coeff
    integer,                              intent(in)  :: nx
    double precision, dimension(nx),      intent(out) :: x
    double precision, dimension(nx),      intent(out) :: y

    double precision, dimension(sfb%r_N) :: f

    double precision :: dx, r_over_Rc, w
    integer :: ix, ic

    write(*,*) "sbspline_reconstruct_radial is not implemented"
    stop

  end subroutine sbspline1d_reconstruct_radial

  !--------------------------------------------------------------------!
  ! reconstruct angular distribution function from basis set expansion !
  !--------------------------------------------------------------------!

  subroutine sbspline1d_reconstruct_angular(sfb, coeff, nx, x, y)

    implicit none

    !------------------------------------------------------------------!
    ! sfb         Instance of BsplineBasis_1d                         !
    ! coeff(i)    Coefficient of the i-th basis function               !
    ! nx          Grid points for function evaluation                  !
    ! x(i)        x value of the i-th grid point (output)              !
    ! y(i)        y (function) value of the i-th grid point (output)   !
    !------------------------------------------------------------------!

    type(BsplineBasis_1d),               intent(in)  :: sfb
    double precision, dimension(sfb%r_N), intent(in)  :: coeff
    integer,                              intent(in)  :: nx
    double precision, dimension(nx),      intent(out) :: x
    double precision, dimension(nx),      intent(out) :: y

    double precision, dimension(sfb%r_N) :: f

    double precision :: dx, r_over_PI, w
    integer :: ix, ic

    write(*,*) "sbspline_reconstruct_angular is not implemented"
    stop

  end subroutine sbspline1d_reconstruct_angular


  !======================== private/auxiliary =========================!


  !--------------------------------------------------------------------!
  !        assert that a BsplineBasis_1d has been initialized         !
  !--------------------------------------------------------------------!

  subroutine sbspline1d_assert_init(sfb)

    implicit none

    type(BsplineBasis_1d), intent(in) :: sfb

    if (.not. sfb%initialized) then
       write(0, *) "Error: BsplineBasis_1d not initialized."
       stop
    end if

  end subroutine sbspline1d_assert_init


  !====================================================================!
  !                                                                    !
  !                          cutoff function                           !
  !                                                                    !
  !====================================================================!


  pure function sbspline1d_fc(Rij, Rc) result(fc)

    implicit none

    double precision, intent(in) :: Rij, Rc
    double precision             :: fc

    if (Rij >= Rc) then
       fc  = 0.0d0
    else
       fc  =  0.5d0*(cos(PI/Rc*Rij) + 1.0d0)
    end if

  end function sbspline1d_fc

  !--------------------------------------------------------------------!

  pure function sbspline1d_fc_d1(Rij, Rc) result(dfc)

    implicit none

    double precision, intent(in) :: Rij, Rc
    double precision             :: dfc

    double precision :: a

    if (Rij >= Rc) then
       dfc = 0.0d0
    else
       a = PI/Rc
       dfc = -0.5d0*a*sin(a*Rij)
    end if

  end function sbspline1d_fc_d1


  !====================================================================!
  !                                                                    !
  !                      generic basis functions                       !
  !                                                                    !
  !====================================================================!


  subroutine sbspline1d_radial(sfb, R_ij, d_ij, values, N,knots,deriv_i, deriv_j)

    implicit none

    type(BsplineBasis_1d),                     intent(inout) :: sfb
    double precision, dimension(3),             intent(in)    :: R_ij
    double precision,                           intent(in)    :: d_ij
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j
    integer,intent(in) ::N
    double precision, dimension(:),  intent(in) ::knots(:) 

    double precision                     :: w_ij, dw_ij
    double precision, dimension(N) :: f, df
    integer                              :: i

    call sbspline1d_assert_init(sfb)

    !w_ij = sbspline_fc(d_ij, sfb%r_Rc)

    
    !f = chebyshev_polynomial(d_ij, 0.0d0, sfb%r_Rc, sfb%r_points)

    !values(1:sfb%r_N) = w_ij*f(1:sfb%r_N)

    if (present(deriv_i) .and. present(deriv_j)) then
       !dw_ij = sbspline_fc_d1(d_ij, sfb%r_Rc)
       call bspline_basis_functions_deriv(values(1:N),df(1:N),d_ij,knots,d)
       !df = chebyshev_polynomial_d1(d_ij, 0.0d0, sfb%r_Rc, sfb%r_points)
       forall (i=1:N)
          !deriv_i(:,i) = -R_ij/d_ij*(dw_ij*f(i) + w_ij*df(i))
          deriv_i(:,i) = -R_ij/d_ij*(df(i))
       end forall
       deriv_j(1:3,1:N) = -deriv_i(1:3,1:N)
    else
      call bspline_basis_functions(values(1:N),d_ij,knots,d)
    end if

  end subroutine sbspline1d_radial


  !--------------------------------------------------------------------!
  subroutine sbspline1d_angular(sfb, R_ij, R_ik, d_ij, d_ik, d_jk, values, N,knots,&
                         deriv_i, deriv_j, deriv_k)

    implicit none

    type(BsplineBasis_1d),                     intent(inout) :: sfb
    double precision, dimension(3),             intent(in)    :: R_ij, R_ik
    double precision,                           intent(in)    :: d_ij, d_ik
    double precision,                           intent(in)    :: d_jk
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j
    double precision, dimension(:,:), optional, intent(out)   :: deriv_k
    integer,intent(in) ::N
    double precision, dimension(:),  intent(in) ::knots(:) 

    double precision                     :: w_ijk
    double precision                     :: fc_j, dfc_j, fc_k, dfc_k
    double precision, dimension(N) :: f, df
    double precision                     :: id_ij2, id_ik2, id_ij_ik
    double precision, dimension(3)       :: di_R_jk, dj_R_jk, dk_R_jk
    !double precision, dimension(3)       :: di_cos_ikj, dj_cos_ikj, dk_cos_ikj
    double precision, dimension(3)       :: di_w_ijk, dj_w_ijk, dk_w_ijk
    integer                              :: i
    integer::k

    call sbspline1d_assert_init(sfb)

    fc_j = sbspline1d_fc(d_ij, sfb%a_Rc)
    fc_k = sbspline1d_fc(d_ik, sfb%a_Rc)
    w_ijk = fc_j*fc_k

    if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
       call bspline_basis_functions_deriv(f(1:N),df(1:N),d_jk,knots,d)

       dfc_j = sbspline1d_fc_d1(d_ij, sfb%a_Rc)
       dfc_k = sbspline1d_fc_d1(d_ik, sfb%a_Rc)
       !df = chebyshev_polynomial_d1(cos_ijk, -1.0d0, 1.0d0, sfb%a_order)
       !id_ij2 = 1.0d0/(d_ij*d_ij)
       !id_ik2 = 1.0d0/(d_ik*d_ik)
       !id_ij_ik = 1.0d0/(d_ij*d_ik)
       ! d/dR_i (cos_ijk)
       !di_cos_ikj = cos_ijk*(R_ij*id_ij2 + R_ik*id_ik2) - (R_ij+R_ik)*id_ij_ik
       ! d/dR_j (cos_ijk)
       !dj_cos_ikj = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
       ! d/R_j (R_jk)
       dj_R_jk = 1d0 
       ! d/dR_k (cos_ijk)
       !dk_cos_ikj = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik
       ! d/R_j (R_jk)
       dk_R_jk = -1d0
       ! d/dR_i (cos_ijk)
       !di_cos_ikj = -dj_cos_ikj - dk_cos_ikj!cos_ijk*(R_ij*id_ij2 + R_ik*id_ik2) - (R_ij+R_ik)*id_ij_ik
       ! d/dR_i (w_ijk)
       !di_w_ijk = -(dfc_j*fc_k*R_ij/d_ij + fc_j*dfc_k*R_ik/d_ik)
       ! d/dR_j (w_ijk)
       dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
       ! d/dR_k (w_ijk)
       dk_w_ijk = fc_j*dfc_k*R_ik/d_ik
       ! d/dR_i (w_ijk)
       !di_w_ijk = -dj_w_ijk  - dk_w_ijk!-(dfc_j*fc_k*R_ij/d_ij + fc_j*dfc_k*R_ik/d_ik)
       !forall (i=1:sfb%a_N)
       do concurrent(k=1:3,i=1:N)
          ! d/dR_i (w_ijk*f)
          !deriv_i(:,i) = di_w_ijk(:)*f(i) + w_ijk*df(i)*di_cos_ikj(:)
          ! d/dR_j (w_ijk*f)
          !deriv_j(:,i) = dj_w_ijk(:)*f(i) + w_ijk*df(i)*dj_cos_ikj(:)
          ! d/dR_k (w_ijk*f)
          !deriv_k(:,i) = dk_w_ijk(:)*f(i) + w_ijk*df(i)*dk_cos_ikj(:)
          ! d/dR_i (w_ijk*f)
          !deriv_i(:,i) = -deriv_j(:,i) -deriv_k(:,i)  !di_w_ijk(:)*f(i) + w_ijk*df(i)*di_cos_ikj(:)
          deriv_j(k,i) = dj_w_ijk(k)*f(i) + w_ijk*df(i)*dj_R_jk(k)
          !deriv_j(k,i) = dj_w_ijk(k)*f(i) + w_ijk*df(i)*dj_cos_ikj(k)
          ! d/dR_k (w_ijk*f)
          !deriv_k(k,i) = dk_w_ijk(k)*f(i) + w_ijk*df(i)*dk_cos_ikj(k)
          deriv_k(k,i) = dk_w_ijk(k)*f(i) + w_ijk*df(i)*dk_R_jk(k)
          ! d/dR_i (w_ijk*f)
          deriv_i(k,i) = -deriv_j(k,i) -deriv_k(k,i)  !di_w_ijk(:)*f(i) + w_ijk*df(i)*di_cos_ikj(:)
           
       end do
       !end forall
    else
      call bspline_basis_functions(f(1:N),d_jk,knots,d)
      !f = chebyshev_polynomial(cos_ijk, -1.0d0, 1.0d0, sfb%a_order)
      values(1:N) = w_ijk*f
    end if

  end subroutine sbspline1d_angular



end module sfbspline_1d
