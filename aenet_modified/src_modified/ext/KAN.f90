module KAN
    use feedforward,only: ff_activate
    implicit none
    type, public :: KAN_descriptor
        integer::npoints
        integer::nk
        real(8),allocatable   ::b_KAN(:)
        real(8),allocatable::values_r1(:,:),values_r2(:,:),values_a1(:,:),values_a2(:,:)
        real(8),allocatable::bk_r1(:),bk_r2(:),bk_a1(:),bk_a2(:)
        real(8),allocatable::points_r(:),points_a(:)
        real(8)::a_Rc,r_Rc
        integer::a_Nc,r_Nc
        real(8)::rc_min
        integer,          dimension(:),   allocatable :: typeid
        double precision, dimension(:),   allocatable :: typespin
        logical ::multi
        integer::first_fun
        integer::ntypes
    end type

    interface KAN_descriptor
        module procedure init_KAN_descriptor
        module procedure copy_KAN_descriptor
    end interface KAN_descriptor

    double precision, parameter, private :: PI     = 3.14159265358979d0
    double precision, parameter, private :: PI_INV = 1.0d0/PI
    double precision, parameter, private :: PI2    = 2.0d0*PI
    double precision, parameter, private :: EPS    = 1.0d-12


    contains

    type(KAN_descriptor) function copy_KAN_descriptor(kan) result(param)
        implicit none
        type(KAN_descriptor),intent(in)::kan
        integer::ntypes
        param%npoints = kan%npoints
        param%nk = kan%nk
        allocate(param%b_KAN(param%nk))
        allocate(param%bk_r1(param%nk))
        allocate(param%bk_a1(param%nk))

        param%multi = kan%multi

        if (param%multi) then
            allocate(param%bk_r2(param%nk))
            allocate(param%bk_a2(param%nk))
            allocate(param%values_r2(param%npoints,param%nk))
            allocate(param%values_a2(param%npoints,param%nk))
        end if

        allocate(param%values_r1(param%npoints,param%nk))
        allocate(param%values_a1(param%npoints,param%nk))

        param%b_KAN = kan%b_KAN(:)
        param%values_r1 = kan%values_r1(:,:)

        if (param%multi) then
            param%values_r2 = kan%values_r2(:,:)
        end if
        param%values_a1 = kan%values_a1(:,:)
        if (param%multi) then
            param%values_a2 = kan%values_a2(:,:)
        end if
        param%bk_r1 = kan%bk_r1(:)
        if (param%multi) then
            param%bk_r2 = kan%bk_r2(:)
        end if
        param%bk_a1 = kan%bk_a1(:)
        if (param%multi) then
            param%bk_a2 = kan%bk_a2(:)
        end if
        allocate(param%points_r(param%npoints),param%points_a(param%npoints))
        param%points_r = kan%points_r(:)
        param%points_a = kan%points_a(:)

        param%r_Rc = kan%r_Rc
        param%r_Nc = kan%r_Nc
        param%a_Rc = kan%a_Rc
        param%rc_min = kan%rc_min
        param%first_fun = kan%first_fun
        ntypes = ubound(kan%typeid,1)
        param%ntypes = ntypes
    
        allocate(param%typeid(ntypes),          &
                  param%typespin(ntypes))
        param%typeid =  kan%typeid
        param%typespin = kan%typespin

    end


    type(KAN_descriptor) function init_KAN_descriptor(ntypes,u) result(param)
        implicit none
        integer,intent(in) ::ntypes
        integer,          optional, intent(in) :: u
        integer::i
        !real(8)::x0,y0,r
        integer ::  s


        read(u) param%npoints
        read(u) param%nk
        allocate(param%b_KAN(param%nk))
        allocate(param%bk_r1(param%nk))
        allocate(param%bk_a1(param%nk))
        

        param%ntypes = ntypes

        if (ntypes > 1) then
            allocate(param%bk_r2(param%nk))
            allocate(param%bk_a2(param%nk))
            allocate(param%values_r2(param%npoints,param%nk))
            allocate(param%values_a2(param%npoints,param%nk))
        end if

        allocate(param%values_r1(param%npoints,param%nk))
        allocate(param%values_a1(param%npoints,param%nk))


        read(u) param%b_KAN
        read(u) param%values_r1
        if (ntypes > 1) then
            read(u) param%values_r2
        end if
        read(u) param%values_a1
        if (ntypes > 1) then
            read(u) param%values_a2
        end if
        read(u) param%bk_r1
        if (ntypes > 1) then
            read(u) param%bk_r2
        end if
        read(u) param%bk_a1
        if (ntypes > 1) then
            read(u) param%bk_a2
        end if
        allocate(param%points_r(param%npoints),param%points_a(param%npoints))


        read(u) param%r_Rc,param%r_Nc,param%a_Rc,param%a_Nc
        read(u) param%rc_min
        read(u) param%first_fun
        do i=1,param%npoints
            param%points_r(i) = dble(i-1)*(param%r_Rc-param%rc_min)/dble(param%npoints-1) + param%rc_min
            param%points_a(i) = dble(i-1)*2/dble(param%npoints-1) -1d0 
        end do

        if (ntypes > 1) then
            param%multi = .true.
         else
            param%multi = .false.
         end if
     
         allocate(param%typeid(ntypes),          &
                  param%typespin(ntypes))

        do i = 1, ntypes
            param%typeid(i) = i
         end do
     
         s = -ntypes/2
         do i = 1,  ntypes
            if ((s == 0) .and. (mod(ntypes, 2) == 0)) s = s + 1
            param%typespin(i) = dble(s)
            s = s + 1
         end do

        !for debug
        !open(101,file="test.txt")
        !open(102,file="test_true.txt")
        !do i=1,param%npoints
        !    call random_number(r)
        !    x0 = r*(param%r_Rc-param%rc_min) + param%rc_min
            !x0 = dble(i-1)*(param%r_Rc-param%rc_min)/dble(2*param%npoints-1) + param%rc_min
        !    call linear_interpolation(param%points_r,  param%values_r1(:,1), param%npoints, x0, y0)
        !    write(101,*) x0,y0
        !    write(102,*) param%points_r(i),param%values_r1(i,1)
        !end do
        !close(101)
        !close(102)
        !stop
    end 


    subroutine linear_interpolation(x, y, n, x0, y0)
        integer, intent(in) :: n           
        real(8), intent(in) :: x(n)        
        real(8), intent(in) :: y(n)        
        real(8), intent(in) :: x0          
        real(8), intent(out) :: y0         
        
        integer :: low, high, mid
        real(8) :: t

        low = 1
        high = n

        do while (high - low > 1)
            mid = (low + high) / 2
            if (x(mid) > x0) then
                high = mid
            else
                low = mid
            end if
        end do

        t = (x0 - x(low)) / (x(high) - x(low))
        y0 = y(low) + t * (y(high) - y(low))
    end subroutine linear_interpolation

    subroutine linear_interpolation_get_highlow(x, n, x0,high,low)
        implicit none
        integer, intent(in) :: n         
        real(8), intent(in) :: x(n)          
        real(8), intent(in) :: x0   
        integer,intent(out) ::high,low    
        integer :: mid

        low = 1
        high = n

        do while (high - low > 1)
            mid = (low + high) / 2
            if (x(mid) > x0) then
                high = mid
            else
                low = mid
            end if
        end do
    end

    subroutine linear_interpolation_multi_highlow(high,low,x, y, n, nk,x0, y0)
        implicit none
        integer, intent(in) :: n,nk           
        real(8), intent(in) :: x(n)        
        real(8), intent(in) :: y(n,nk)        
        real(8), intent(in) :: x0          
        real(8), intent(out) :: y0(nk)   
        integer,intent(in)::high,low      
        
        integer ::  mid,ik
        real(8) :: t
        t = (x0 - x(low)) / (x(high) - x(low))
        do ik = 1,nk
            y0(ik) = y(low,ik) + t * (y(high,ik) - y(low,ik))
        end do
    end


    subroutine linear_interpolation_multi_highlow_df(high,low,x, y, n, nk,x0, y0,dy0)
        implicit none
        integer, intent(in) :: n,nk           
        real(8), intent(in) :: x(n)        
        real(8), intent(in) :: y(n,nk)        
        real(8), intent(in) :: x0          
        real(8), intent(out) :: y0(nk),dy0(nk)   
        integer,intent(in)::high,low      
        
        integer ::  ik
        real(8) :: t
        t = (x0 - x(low)) / (x(high) - x(low))
        do ik = 1,nk
            y0(ik) = y(low,ik) + t * (y(high,ik) - y(low,ik))
            dy0(ik) =  (y(high,ik) - y(low,ik))/(x(high) - x(low))
        end do
    end

    subroutine linear_interpolation_multi_highlow_dfonly(high,low,x, y, n, nk,dy0)
        implicit none
        integer, intent(in) :: n,nk           
        real(8), intent(in) :: x(n)        
        real(8), intent(in) :: y(n,nk)       
        real(8), intent(out) :: dy0(nk)   
        integer,intent(in)::high,low      
        
        integer ::  ik
        do ik = 1,nk
            dy0(ik) =  (y(high,ik) - y(low,ik))/(x(high) - x(low))
        end do
    end

    subroutine linear_interpolation_multi(x, y, n, nk,x0, y0)
        implicit none
        integer, intent(in) :: n,nk           
        real(8), intent(in) :: x(n)        
        real(8), intent(in) :: y(n,nk)        
        real(8), intent(in) :: x0          
        real(8), intent(out) :: y0(nk)         
        
        integer :: low, high, mid,ik
        real(8) :: t

        low = 1
        high = n

        do while (high - low > 1)
            mid = (low + high) / 2
            if (x(mid) > x0) then
                high = mid
            else
                low = mid
            end if
        end do

        t = (x0 - x(low)) / (x(high) - x(low))
        do ik = 1,nk
            y0(ik) = y(low,ik) + t * (y(high,ik) - y(low,ik))
        end do
    end subroutine linear_interpolation_multi

  subroutine kan_radial(kan,high,low, R_ij, d_ij,points_r,values_r,values, deriv_i, deriv_j)
    implicit none

    type(KAN_descriptor),                     intent(in) :: kan
    integer,intent(in)::high,low
    double precision, dimension(3),             intent(in)    :: R_ij
    double precision,                           intent(in)    :: d_ij
    double precision, dimension(:,:),           intent(in)   :: values_r
    double precision, dimension(:),             intent(in)   ::points_r
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j

    !double precision                     :: w_ij, dw_ij
    double precision, dimension(kan%nk) :: f, df
    integer                              :: i


    call linear_interpolation_multi_highlow_df(high,low,points_r, values_r,&
            kan%npoints,kan%nk,d_ij, f,df)
    
    if (present(deriv_i) .and. present(deriv_j)) then   
        do i=1,kan%nk
            deriv_i(:,i) = -df(i)*R_ij/d_ij
        end do
        deriv_j(1:3,1:kan%nk) = - deriv_i(1:3,1:kan%nk)
    end if   

    values(1:kan%nk) = f(1:kan%nk)
  end subroutine kan_radial    

  subroutine kan_angular(kan,high,low, R_ij, R_ik, d_ij, d_ik, cos_ijk, points_a,values_a,values, &
                         deriv_i, deriv_j, deriv_k)

    implicit none

    type(KAN_descriptor),                     intent(in) :: kan
    integer,intent(in)::high,low
    double precision, dimension(3),             intent(in)    :: R_ij, R_ik
    double precision,                           intent(in)    :: d_ij, d_ik
    double precision, dimension(:,:),           intent(in)   :: values_a
    double precision, dimension(:),              intent(in)   :: points_a
    double precision,                           intent(in)    :: cos_ijk
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j
    double precision, dimension(:,:), optional, intent(out)   :: deriv_k

    double precision                     :: w_ijk
    double precision                     :: fc_j, dfc_j, fc_k, dfc_k
    double precision, dimension(kan%nk) :: f, df
    double precision                     :: id_ij2, id_ik2, id_ij_ik
    double precision, dimension(3)       :: di_cos_ikj, dj_cos_ikj, dk_cos_ikj
    double precision, dimension(3)       :: di_w_ijk, dj_w_ijk, dk_w_ijk
    integer                              :: i


    fc_j = sfb_fc(d_ij, kan%a_Rc)
    fc_k = sfb_fc(d_ik, kan%a_Rc)
    w_ijk = fc_j*fc_k

    call linear_interpolation_multi_highlow_df(high,low,points_a, values_a,&
            kan%npoints,kan%nk,cos_ijk, f,df)

    values(1:kan%nk) = w_ijk*f

    if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
       dfc_j = sfb_fc_d1(d_ij, kan%a_Rc)
       dfc_k = sfb_fc_d1(d_ik, kan%a_Rc)
       
       id_ij2 = 1.0d0/(d_ij*d_ij)
       id_ik2 = 1.0d0/(d_ik*d_ik)
       id_ij_ik = 1.0d0/(d_ij*d_ik)
       ! d/dR_i (cos_ijk)
       di_cos_ikj = cos_ijk*(R_ij*id_ij2 + R_ik*id_ik2) - (R_ij+R_ik)*id_ij_ik
       ! d/dR_j (cos_ijk)
       dj_cos_ikj = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
       ! d/dR_k (cos_ijk)
       dk_cos_ikj = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik
       ! d/dR_i (w_ijk)
       di_w_ijk = -(dfc_j*fc_k*R_ij/d_ij + fc_j*dfc_k*R_ik/d_ik)
       ! d/dR_j (w_ijk)
       dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
       ! d/dR_k (w_ijk)
       dk_w_ijk = fc_j*dfc_k*R_ik/d_ik
       forall (i=1:kan%nk)
          ! d/dR_i (w_ijk*f)
          deriv_i(:,i) = di_w_ijk(:)*f(i) + w_ijk*df(i)*di_cos_ikj(:)
          ! d/dR_j (w_ijk*f)
          deriv_j(:,i) = dj_w_ijk(:)*f(i) + w_ijk*df(i)*dj_cos_ikj(:)
          ! d/dR_k (w_ijk*f)
          deriv_k(:,i) = dk_w_ijk(:)*f(i) + w_ijk*df(i)*dk_cos_ikj(:)
       end forall
    end if

  end subroutine kan_angular  


subroutine sfb_eval_kan(kan, coo0, nat, itype1, coo1, nv,&
                      values, deriv0, deriv1)

    implicit none

    type(KAN_descriptor),                          intent(inout) :: kan
    !integer,                                         intent(in)    :: itype0
    double precision, dimension(3),                  intent(in)    :: coo0
    integer,                                         intent(in)    :: nat
    integer,          dimension(nat),                intent(in)    :: itype1
    double precision, dimension(3,nat),              intent(in)    :: coo1
    integer,                                         intent(in)    :: nv
    double precision, dimension(nv),                 intent(out)   :: values
    double precision, dimension(3,nv),     optional, intent(out)   :: deriv0
    double precision, dimension(3,nv,nat), optional, intent(out)   :: deriv1
    !logical,intent(in) ::multi
    !integer,          dimension(:),intent(in) :: typeid
    !double precision, dimension(:),intent(in) :: typespin

    double precision, dimension(kan%nk)   :: sfb_values,work1,work2
    double precision, dimension(:,:), allocatable :: sfb_deriv_i
    double precision, dimension(:,:), allocatable :: sfb_deriv_j
    double precision, dimension(:,:), allocatable :: sfb_deriv_k

    logical                        :: do_deriv
    double precision, dimension(3) :: R_ij, R_ik
    double precision               :: d_ij, d_ik
    double precision               :: cos_ijk
    double precision               :: s_j, s_k
    integer                        :: j, k
    real(8)::fc_j,fc_k,w_ijk 
    integer::high,low

    !write(*,* ) "in sfb_eval_kan"

    if (present(deriv0) .and. present(deriv1)) then
       do_deriv = .true.
       deriv0(:,:) = 0.0d0
       deriv1(:,:,:) = 0.0d0
       allocate(sfb_deriv_i(3, kan%nk), &
                sfb_deriv_j(3, kan%nk), &
                sfb_deriv_k(3, kan%nk))
    else
       do_deriv = .false.
    end if

    values= kan%b_KAN + kan%bk_r1 + kan%bk_a1
    if (kan%multi) then
        values = values + kan%bk_r2 + kan%bk_a2
    end if

    

    s_j = 1.0d0
    for_j : do j = 1, nat
       R_ij = coo1(1:3, j) - coo0(1:3)
       d_ij = sqrt(dot_product(R_ij, R_ij))
       if ((d_ij <= kan%r_Rc) .and. (d_ij > EPS)) then
            !write(*,*) j,d_ij
          ! evaluate radial basis functions
          !i1 = sfb%r_i1
          !i2 = sfb%r_f1
          !N = sfb%r_N
          call linear_interpolation_get_highlow(kan%points_r, kan%npoints, d_ij,high,low)
          if (do_deriv) then
            
            call kan_radial(kan,high,low, R_ij, d_ij,kan%points_r,kan%values_r1,sfb_values, &
                    deriv_i=sfb_deriv_i, deriv_j=sfb_deriv_j)
            !call linear_interpolation_multi_highlow(high,low,kan%points_r, kan%values_r1,kan%npoints,kan%nk,d_ij,sfb_values)

            !call linear_interpolation_multi(kan%points_r, kan%values_r1,kan%npoints,kan%nk,d_ij,sfb_values)
            values = values +  sfb_values
            !write(*,*) "v1",sfb_values
             !call sfb_radial(sfb, R_ij, d_ij, sfb_values, &
             !                deriv_i=sfb_deriv_i, deriv_j=sfb_deriv_j)
             !values(i1:i2) = values(i1:i2) + sfb_values(1:N)
             deriv0(1:3, :) = deriv0(1:3, :) + sfb_deriv_i(1:3, :)
             deriv1(1:3, :, j) = deriv1(1:3, :, j) + sfb_deriv_j(1:3, :)
          else
            call kan_radial(kan,high,low, R_ij, d_ij,kan%points_r,kan%values_r1,sfb_values)
             !call sfb_radial(sfb, R_ij, d_ij, sfb_values)
             values = values +  sfb_values
             !values(i1:i2) = values(i1:i2) + sfb_values(1:N)
          end if

          ! redundant radial basis in case of multi-component systems
          !i1 = sfb%r_i2
          !i2 = sfb%r_f2
          !N = sfb%r_N
          if (kan%multi) then
             s_j = kan%typespin(kan%typeid(itype1(j)))
             !s_j = 0 !debug
             !write(*,*) j ,s_j 
             !call linear_interpolation_multi_highlow(high,low,kan%points_r, kan%values_r2,kan%npoints,kan%nk,d_ij,sfb_values)
             !call linear_interpolation_multi(kan%points_r, kan%values_r2,kan%npoints,kan%nk,d_ij,sfb_values)
             call kan_radial(kan,high,low, R_ij, d_ij,kan%points_r,kan%values_r2,sfb_values, &
                    deriv_i=sfb_deriv_i, deriv_j=sfb_deriv_j)
             values = values + s_j*sfb_values
             !write(*,*) "v2", s_j*sfb_values
             if (do_deriv) then
                deriv0(1:3, :) = deriv0(1:3, :) &
                                   + s_j*sfb_deriv_i(1:3, :)
                deriv1(1:3, :, j) = deriv1(1:3, :, j) &
                                      + s_j*sfb_deriv_j(1:3, :)
             end if
          end if

       end if  ! within radial cutoff

       if (d_ij > kan%a_Rc) cycle for_j
       for_k : do k = j+1, nat
          R_ik = coo1(1:3, k) - coo0(1:3)
          d_ik = sqrt(dot_product(R_ik, R_ik))

          !write(*,*) j,k,d_ik
          
          if ((d_ik > kan%a_Rc) .or. (d_ik < EPS)) cycle for_k
          cos_ijk = dot_product(R_ij, R_ik)/(d_ij*d_ik)
          !if (cos_ijk > 1d0) cos_ijk = 1d0
          !if (cos_ijk < -1d0) cos_ijk = -1d0

          fc_j = sfb_fc(d_ij, kan%a_Rc)
          fc_k = sfb_fc(d_ik, kan%a_Rc)
          w_ijk = fc_j*fc_k

          ! evaluate angular basis functions
          !i1 = sfb%a_i1
          !i2 = sfb%a_f1
          !N = sfb%a_N
          call linear_interpolation_get_highlow(kan%points_a, kan%npoints, cos_ijk,high,low)
          if (do_deriv) then
            
            call kan_angular(kan,high,low, R_ij, R_ik, d_ij, d_ik, cos_ijk,kan%points_a, &
                    kan%values_a1,sfb_values, deriv_i=sfb_deriv_i,      &
                    deriv_j=sfb_deriv_j, deriv_k=sfb_deriv_k)
            !call linear_interpolation_multi_highlow(high,low,kan%points_a, kan%values_a1,kan%npoints,kan%nk,cos_ijk,sfb_values)
            !write(*,*) "--------------"
            !write(*,*) w_ijk,cos_ijk,high,low
            !write(*,*) sfb_values
            !write(*,*) "--------------"
            !stop
             !call linear_interpolation_multi(kan%points_a, kan%values_a1,kan%npoints,kan%nk,cos_ijk,sfb_values)
             !sfb_values = sfb_values*w_ijk
             values = values + sfb_values
            ! write(*,*) "a1",sfb_values*w_ijk
             !call sfb_angular(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, &
             !                 sfb_values, deriv_i=sfb_deriv_i,      &
             !                 deriv_j=sfb_deriv_j, deriv_k=sfb_deriv_k)
             !values(i1:i2) = values(i1:i2) + sfb_values(1:N)
              deriv0(1:3, :) = deriv0(1:3, :) + sfb_deriv_i(1:3, :)
              deriv1(1:3, :, j) = deriv1(1:3, :, j) + sfb_deriv_j(1:3, :)
              deriv1(1:3, :, k) = deriv1(1:3, :, k) + sfb_deriv_k(1:3, :)
          else
            call kan_angular(kan,high,low, R_ij, R_ik, d_ij, d_ik, cos_ijk,kan%points_a, &
                    kan%values_a1,sfb_values)
             !call sfb_angular(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, sfb_values)
             values = values + sfb_values
          end if

          ! redundant angular basis in case of multi-component systems
          !i1 = sfb%a_i2
          !i2 = sfb%a_f2
          !N = sfb%a_N
          if (kan%multi) then
             s_k = kan%typespin(kan%typeid(itype1(k)))
             !s_k = 0 !debug
             !s_j = 0 !debug
             !write(*,*) j ,k,s_j,s_k 
             !call linear_interpolation_multi_highlow(high,low,kan%points_a, kan%values_a2,kan%npoints,kan%nk,cos_ijk,sfb_values)
             
             !call linear_interpolation_multi(kan%points_a, kan%values_a2,kan%npoints,kan%nk,cos_ijk,sfb_values)
             !sfb_values = sfb_values*w_ijk

             

             !write(*,*) "a2",s_j*s_k*sfb_values*w_ijk
             if (do_deriv) then
                call kan_angular(kan,high,low, R_ij, R_ik, d_ij, d_ik, cos_ijk,kan%points_a, &
                    kan%values_a2,sfb_values, deriv_i=sfb_deriv_i,      &
                    deriv_j=sfb_deriv_j, deriv_k=sfb_deriv_k)
                values = values + s_j*s_k*sfb_values!*w_ijk

                deriv0(1:3, :) = deriv0(1:3, :) &
                                   + s_j*s_k*sfb_deriv_i(1:3, :)
                deriv1(1:3, :, j) = deriv1(1:3, :, j) &
                                      + s_j*s_k*sfb_deriv_j(1:3, :)
                deriv1(1:3, :, k) = deriv1(1:3, :, k) &
                                      + s_j*s_k*sfb_deriv_k(1:3, :)
             else
                call kan_angular(kan,high,low, R_ij, R_ik, d_ij, d_ik, cos_ijk,kan%points_a, &
                    kan%values_a2,sfb_values)
                values = values + s_j*s_k*sfb_values!*w_ijk
             end if
             
          end if
       end do for_k
    end do for_j

    !write(*,*) values
    call ff_activate(kan%first_fun, values, work1, work2)
    values = work1 
    do k =1,kan%nk
        deriv0(1:3,k) = deriv0(1:3,k)*work2(k)
        deriv1(1:3,k,:) = deriv1(1:3,k,:)*work2(k)
    end do

    if (do_deriv) deallocate(sfb_deriv_i, sfb_deriv_j, sfb_deriv_k)

  end subroutine sfb_eval_kan

  function sfb_fc(Rij, Rc) result(fc)

    implicit none

    double precision, intent(in) :: Rij, Rc
    double precision             :: fc

    if (Rij >= Rc) then
       fc  = 0.0d0
    else
       fc  =  0.5d0*(cos(PI/Rc*Rij) + 1.0d0)
    end if

  end function sfb_fc

  function sfb_fc_d1(Rij, Rc) result(dfc)
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

  end function sfb_fc_d1

end module 