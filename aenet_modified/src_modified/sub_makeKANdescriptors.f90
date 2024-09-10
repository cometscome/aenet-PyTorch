module makeKAN
    implicit none

    type, public :: Model_parameters
        integer              :: nlayers
        integer              :: nnodesmax
        integer              :: Wsize
        integer              :: nvalues
        integer, allocatable :: nnodes(:), fun(:), iw(:), iv(:)
        real(8), allocatable  :: W(:)

        real(8), allocatable  :: W_KAN(:,:)
        real(8),allocatable   ::b_KAN(:)
        real(8), allocatable  :: W_MLP(:)

        character(len=1024)           :: description
        character(len=2)              :: atomtype
        integer              :: nenv
        character(len=2), allocatable :: envtypes(:)

        real(8)                        :: rc_min, rc_max
        character(len=100)            :: sftype
        integer                       :: nsf, nsfparam
        integer, allocatable          :: sf(:), sfenv(:,:)
        real(8), allocatable           :: sfparam(:,:), sfval_min(:), sfval_max(:), sfval_avg(:), sfval_cov(:)
        integer::neval
        character(len=1024)           :: file
        logical                       :: normalized
        real(8)                        :: scale, shift, E_min, E_max, E_avg
        integer                       :: ntypes, natomtot, nstrucs
        character(len=2), allocatable :: type_names(:)
        real(8), allocatable           :: E_atom(:)
        real(8)::a_Rc,r_Rc
        integer::a_Nc,r_Nc
        integer::nk
        character(len=1024)           :: infile
        character(len=1024) ::model_r1,model_r2,model_a1,model_a2

        integer              :: nlayers_KAN
        integer              :: nnodesmax_KAN
        integer              :: Wsize_KAN
        integer              :: nvalues_KAN
        integer, allocatable :: nnodes_KAN(:), fun_KAN(:), iw_KAN(:), iv_KAN(:)

        integer::npoints

        real(8),allocatable::values_r1(:,:),values_r2(:,:),values_a1(:,:),values_a2(:,:)
        real(8),allocatable::bk_r1(:),bk_r2(:),bk_a1(:),bk_a2(:)
        integer::first_fun

    end type

    interface Model_parameters
        module procedure init_Model_parameters
    end interface Model_parameters


    double precision, parameter, private :: PI     = 3.14159265358979d0
    double precision, parameter, private :: PI_INV = 1.0d0/PI
    double precision, parameter, private :: PI2    = 2.0d0*PI
    double precision, parameter, private :: EPS    = 1.0d-12
    contains 

    subroutine calc_bk(bk,Wkn,nc,nk,sfval_avg,sfval_cov)
        implicit none
        integer,intent(in)::nc,nk
        real(8),intent(out)::bk(1:nk)
        real(8),intent(in)::sfval_avg(1:nc+1),sfval_cov(1:nc+1)
        real(8),intent(in)::Wkn(1:nk,1:nc+1)
        
        integer::isf
        !integer::iw1,iw2
        real(8)::scale,shift,s
        real(8)::cn(1:nc+1)



        !iw1 = 1
        !iw2 = iw1 -1 + nc +1
        do isf = 1, nc+1
            shift = sfval_avg(isf)!isf+iw1-1)
            s = sqrt(sfval_cov(isf) - shift*shift)
            !s = sqrt(sfval_cov(isf+iw1-1) - shift*shift)
            if (s <= 1.0d-10) then
                scale = 0.0d0
            else
                scale = 1.0d0/s
            end if
            !write(*,*) "scale",scale
            cn(isf) = scale*(-shift)
        end do
        bk = matmul(Wkn,cn)



    end subroutine

    subroutine write_model(param,outfile)
        implicit none
        character(len=*), intent(in)  :: outfile
        type(Model_parameters),intent(inout)::param

        param%nsf = param%nk
        deallocate(param%sf,param%sfparam,param%sfenv,param%sfval_min,param%sfval_max,param%sfval_avg,param%sfval_cov)
        allocate(param%sf(param%nsf), param%sfparam(param%nsfparam, param%nsf), &
            param%sfenv(2, param%nsf), param%sfval_min( param%nsf), &
            param%sfval_max( param%nsf), param%sfval_avg( param%nsf), param%sfval_cov( param%nsf))

        open(unit = 2, action = "write", status = "replace", file = outfile, form = "unformatted")


        write(2) param%nlayers_KAN
        write(2) param%nnodesmax_KAN
        write(2) param%Wsize_KAN
        write(2) param%nvalues_KAN
        write(2) param%nnodes_KAN(:)
        write(2) param%fun_KAN(:)
        write(2) param%iw_KAN(:)
        write(2) param%iv_KAN(:)
        write(2) param%W_MLP(:)

        write(2) param%description
        write(2) param%atomtype
        write(2) param%nenv
        write(2) param%envtypes(:)
        write(2) param%rc_min
        write(2) param%rc_max
        write(2) param%sftype

        write(2) param%nsf
        write(2) param%nsfparam
        write(2) param%sf(:)
        write(2) param%sfparam(:,:)
        write(2) param%sfenv(:,:)
        write(2) param%neval
        write(2) param%sfval_min
        write(2) param%sfval_max
        write(2) param%sfval_avg
        write(2) param%sfval_cov

        write(2) param%file
        write(2) param%normalized
        write(2) param%scale
        write(2) param%shift
        write(2) param%ntypes
        write(2) param%type_names(:)
        write(2) param%E_atom(:)
        write(2) param%natomtot
        write(2) param%nstrucs
        write(2) param%E_min, param%E_max, param%E_avg


        write(2) param%npoints
        write(2) param%nk
        write(2) param%b_KAN
        write(2) param%values_r1
        if (param%ntypes > 1) then
            write(2) param%values_r2
        end if
        write(2) param%values_a1
        if (param%ntypes > 1) then
            write(2) param%values_a2
        end if
        write(2) param%bk_r1
        if (param%ntypes > 1) then
            write(2) param%bk_r2
        end if
        write(2) param%bk_a1
        if (param%ntypes > 1) then
            write(2) param%bk_a2
        end if
        write(2) param%r_Rc,param%r_Nc,param%a_Rc,param%a_Nc
        write(2) param%rc_min
        write(2) param%first_fun

        
        !write(*,*) param%bk_r1,param%bk_r2,param%bk_a1,param%bk_a2
        
    
        close(2)
 

    end

    subroutine construct_KAN_descriptor(param,npoints)
        implicit none
        type(Model_parameters),intent(inout)::param
        integer,intent(in) ::npoints
        real(8)::r_Rc,Rmin,d_ij
        integer::i,iw1,iw2,nk
        !real(8),allocatable::values_r1(:),values_r2(:),values_a1(:),values_a2(:)
        !real(8),allocatable::bk_r1(:),bk_r2(:),bk_a1(:),bk_a2(:)
        real(8)::w_ij
        real(8)::cos_theta_ijk
        !integer::isf
        !real(8)::shift,scale,s
        real(8),allocatable::cn(:)
        real(8),allocatable::values(:)

        r_Rc = param%r_Rc
        Rmin = param%rc_min
        nk = param%nk
        param%npoints = npoints
       ! allocate(values_r1(1:nk),values_r2(1:nk),values_a1(1:nk),values_a2(1:nk))
       ! allocate(bk_r1(1:nk),bk_r2(1:nk),bk_a1(1:nk),bk_a2(1:nk))
        allocate(cn(1:param%r_Nc+1))

        allocate(param%values_r1(npoints,1:param%nk),param%values_r2(npoints,1:param%nk),&
                param%values_a1(npoints,1:param%nk),param%values_a2(npoints,1:param%nk))
        param%values_r1 = 0d0
        param%values_a1 = 0d0
        allocate(values(1:param%nk))
        values = 0d0

        param%values_r2 = 0d0
        param%values_a2 = 0d0


        allocate(param%bk_r1(1:param%nk),param%bk_r2(1:param%nk),&
                param%bk_a1(1:param%nk),param%bk_a2(1:param%nk))

        param%bk_r1 = 0d0
        param%bk_r2 = 0d0
        param%bk_a1 = 0d0
        param%bk_a2 = 0d0




 
        open(11,file= trim(param%model_r1))!trim(param%infile)//"_CKAN_r1.txt")
        if (param%ntypes > 1) then
            open(12,file= trim(param%model_r2))!trim(param%infile)//"_CKAN_r2.txt")
            open(14,file= trim(param%model_a2))!trim(param%infile)//"_CKAN_a2.txt")
        end if
        open(13,file= trim(param%model_a1))!trim(param%infile)//"_CKAN_a1.txt")
        

        param%bk_r1 = 0d0
        

        iw1 = 1
        iw2 = iw1 -1 + param%r_Nc +1
        !write(*,*) iw1,iw2
        call calc_bk(param%bk_r1,param%W_KAN(1:nk,iw1:iw2),param%r_nc,nk,param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2))

        if (param%ntypes > 1) then
            iw1 = param%r_Nc +1 + param%a_Nc +1 + 1
            !iw2 + 1
            iw2 = iw1 -1 + param%r_Nc +1
            !write(*,*) iw1,iw2
            call calc_bk(param%bk_r2,param%W_KAN(1:nk,iw1:iw2),param%r_nc,nk,param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2))
        end if

        write(11,*) "# ",nk, r_Rc, Rmin,npoints
        write(11,*) "# ",param%bk_r1

        if (param%ntypes > 1) then
            write(12,*) "# ",nk, r_Rc, Rmin,npoints
            write(12,*) "# ",param%bk_r2
        end if


        do i=1,npoints  
            d_ij = dble(i-1)*(r_Rc-Rmin)/dble(npoints-1) + Rmin
            w_ij = sfb_fc(d_ij, r_Rc)
            iw1 = 1
            iw2 = iw1 -1 + param%r_Nc +1
            call evaluate_KANfunction(param%W_KAN(1:nk,iw1:iw2),nk,d_ij,0d0,r_Rc,param%r_Nc+1,&
                param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2),values)
            param%values_r1(i,:) = values*w_ij
            write(11,*) d_ij,param%values_r1(i,:)
            
            if (param%ntypes > 1) then
                iw1 = param%r_Nc +1 + param%a_Nc +1 + 1
                iw2 = iw1 -1 + param%r_Nc +1
                !write(*,*) iw1,iw2
                call evaluate_KANfunction(param%W_KAN(1:nk,iw1:iw2),nk,d_ij,0d0,r_Rc,param%r_Nc+1,&
                    param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2),values)  
                param%values_r2(i,:) = values*w_ij
                write(12,*) d_ij,param%values_r2(i,:)
                !write(*,*) d_ij,values_r1,values_r2
            end if
            
           
        end do

        r_Rc = param%a_Rc
        Rmin = param%rc_min

        !if (param%ntypes > 1) then
        !    iw1 = 2*(param%r_Nc +1) + 1
        !else
        !    iw1 = (param%r_Nc +1) + 1
        !end if
        iw1 = (param%r_Nc +1) + 1
        iw2 = iw1 -1 + param%a_Nc +1
        !write(*,*) iw1,iw2
        call calc_bk(param%bk_a1,param%W_KAN(1:nk,iw1:iw2),param%a_nc,nk,param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2))
        write(13,*) "# ",nk, r_Rc, Rmin,npoints
        write(13,*) "# ",param%bk_r1


        if (param%ntypes > 1) then
            iw1 = param%r_Nc +1 + param%a_Nc +1 + (param%r_Nc +1) + 1
            iw2 = iw1 -1 + param%a_Nc +1
            !write(*,*) iw1,iw2
            
            call calc_bk(param%bk_a2,param%W_KAN(1:nk,iw1:iw2),param%a_nc,nk,param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2))
            write(14,*) "# ",nk, r_Rc, Rmin,npoints
            write(14,*) "# ",param%bk_r2

        end if


        do i=1,npoints  
            cos_theta_ijk = dble(i-1)*2/dble(npoints-1) -1d0
            iw1 = (param%r_Nc +1) + 1
            iw2 = iw1 -1 + param%a_Nc +1
            call evaluate_KANfunction(param%W_KAN(1:nk,iw1:iw2),nk,cos_theta_ijk,-1d0,1d0,param%a_Nc+1,&
                param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2),values)
            param%values_a1(i,:) = values
            write(13,*) cos_theta_ijk,param%values_a1(i,:)

            if (param%ntypes > 1) then
                iw1 = param%r_Nc +1 + param%a_Nc +1 + (param%r_Nc +1) + 1
                iw2 = iw1 -1 + param%a_Nc +1
               
                call evaluate_KANfunction(param%W_KAN(1:nk,iw1:iw2),nk,cos_theta_ijk,-1d0,1d0,param%a_Nc+1,&
                    param%sfval_avg(iw1:iw2),param%sfval_cov(iw1:iw2),values)
                param%values_a2(i,:) = values
                !write(*,*) d_ij,values_r1,values_r2
                
                write(14,*) cos_theta_ijk,param%values_a2(i,:)
            end if
        end do

        close(11)
        if (param%ntypes > 1) then
            close(12)
            close(14)
        end if
        close(13)
        

    end subroutine

    subroutine evaluate_KANfunction(Wkn,nk,x,xa,xb,r_N,sfval_avg,sfval_cov,values)
        implicit none
        integer,intent(in)::nk !number of KAN functions
        integer,intent(in)::r_N !number of radial basis
        real(8),intent(in)::Wkn(nk,r_N) !(nk,r_N)
        real(8),intent(in)::x,xa,xb
        real(8),intent(in)::sfval_cov(1:r_N),sfval_avg(1:r_N)
        double precision, dimension(r_N) :: f
        double precision, dimension(nk),intent(out):: values
        integer::r_order,isf
        real(8)::shift,scale,s

        r_order = r_N - 1
        f(1:r_N) = chebyshev_polynomial(x, xa, xb, r_order)
        do isf = 1, r_N
            shift = sfval_avg(isf)
            ! scale covariance to 1
            ! s = sqrt(stp%sfval_cov(isf) + shift*shift - 2.0d0*shift*stp%sfval_avg(isf))
            s = sqrt(sfval_cov(isf) - shift*shift)
            if (s <= 1.0d-10) then
                scale = 0.0d0
            else
                scale = 1.0d0/s
            end if
            !write(*,*) "scale_w",scale
            f(isf) = scale*(f(isf))!-shift)
        end do

        !f(r_N+1) = 1.0d0 
        values = matmul(Wkn,f)
    end subroutine 

    type(Model_parameters) function init_Model_parameters(infile) result(param)
        implicit none
        character(len=*), intent(in)  :: infile
        !integer,intent(in) ::npoints
        integer::iw1,iw2
        integer::i
        !logical::debugmode

        open(unit = 1, action = "read", status = "old", file = infile)
        param%infile = trim(infile)
        param%model_r1 = trim(param%infile)//"_CKAN_r1.txt"
        param%model_r2 = trim(param%infile)//"_CKAN_r2.txt"
        param%model_a1 = trim(param%infile)//"_CKAN_a1.txt"
        param%model_a2 = trim(param%infile)//"_CKAN_a2.txt"


        !write(*,*) infile
        ! Network information
        read(1,*) param%nlayers
        read(1,*) param%nnodesmax
        read(1,*) param%Wsize
        read(1,*) param%nvalues

        allocate(param%nnodes( param%nlayers), param%fun( param%nlayers-1), &
            param%iw( param%nlayers), param%iv( param%nlayers), param%W( param%Wsize))

        read(1,*) param%nnodes(:)
        read(1,*) param%fun(:)
        read(1,*) param%iw(:)
        read(1,*) param%iv(:)
        read(1,*) param%W(:)

        param%nlayers_KAN = param%nlayers -1



        param%nnodesmax_KAN = 1
        do i = 2,param%nlayers
            if (param%nnodesmax_KAN < param%nnodes(i)) then
                param%nnodesmax_KAN = param%nnodes(i)
            end if
        end do



        iw1 = param%nnodes(1)
        iw2 = param%nnodes(2)
        param%nk = iw2
        allocate(param%W_KAN(iw2,iw1))
        allocate(param%b_KAN(iw2))
        param%W_KAN = reshape(param%W(1:iw1*iw2), &
            (/ iw2,iw1 /) )
        param%b_KAN = param%W(iw1*iw2+1:iw1*iw2+iw2)

        !write(*,*) "number of layers in the network, inluding input and output layer  "
        !write(*,*) param%nlayers
        
        !write(*,*) "max. number of nodes in a layer of this network"
        !write(*,*) param%nnodesmax
        !write(*,*) "total number of weights; Wsize = size(W) "
        !write(*,*) param%Wsize
        !write(*,*) "total number of nodes/neurons in the network "
        !write(*,*) param%nvalues
        !write(*,*) "number of nodes (without bias) in the i-th layer"
        !write(*,*) param%nnodes(:)


        !write(*,*) "activation function type for the i-th layer "
        !write(*,*) param%fun(:)
        !write(*,*) "index of the last weight of the i-th layer  "
        !write(*,*) "e.g.: W(iw(2)) --> last weight in the 2nd layer "
        !write(*,*) param%iw(:)
        !write(*,*) "index of the last node in the i-th layer "
        !write(*,*) param%iv(:)
        !write(*,*) "i-th weight of the graph edges (including bias) "
        !write(*,*) param%W(:)

        !write(*,*) "weight for KAN (without bias) "
        !write(*,*) param%W_KAN
        !write(*,*) "bias for KAN  "
        !write(*,*) param%b_KAN

        !write(*,*) "max. number of nodes in a layer of KAN network"
        !write(*,*) param%nnodesmax_KAN

        allocate(param%W_MLP(1:ubound(param%W,1)-(iw1*iw2+iw2+1)+1))
        param%W_MLP(1:ubound(param%W,1)-(iw1*iw2+iw2+1)+1) = param%W(iw1*iw2+iw2+1:ubound(param%W,1))
        !write(*,*) "i-th weight of the graph edges (including bias) for new network "
        !write(*,*) param%W_MLP(:)

        !write(*,*) "total number of weights; Wsize_KAN = size(W_MLP) "
        param%Wsize_KAN = ubound(param%W_MLP,1)
        !write(*,*) param%Wsize_KAN


        allocate(param%nnodes_KAN( param%nlayers_KAN), param%fun_KAN( param%nlayers_KAN-1), &
        param%iw_KAN( param%nlayers_KAN), param%iv_KAN( param%nlayers_KAN))

        param%nnodes_KAN(1:param%nlayers_KAN) = param%nnodes(2:param%nlayers)
        param%fun_KAN(1:param%nlayers_KAN-1) = param%fun(2:param%nlayers-1)
        param%first_fun = param%fun(1)
        !write(*,*) "number of nodes (without bias) in the i-th layer for KAN"
        !write(*,*) param%nnodes_KAN(:)
        !write(*,*) "activation function type for the i-th layer for KAN"
        !write(*,*) param%fun_KAN(:)
        param%nvalues_KAN = 0
        do i=1,param%nlayers_KAN
            param%nvalues_KAN = param%nvalues_KAN +  param%nnodes_KAN(i)+1
            param%iw_KAN(i) = param%iw(i+1) - (param%nnodes(1)+1)*param%nnodes(2)
            param%iv_KAN(i) = param%iv(i+1) - (param%nnodes(1)+1)
        end do

        !write(*,*) "index of the last weight of the i-th layer for KAN "
        !write(*,*) "e.g.: W(iw_KAN(2)) --> last weight in the 2nd layer "
        !write(*,*) param%iw_KAN(:)
        !write(*,*) "index of the last node in the i-th layer for KAN"
        !write(*,*) param%iv_KAN(:)

        !write(*,*) "total number of nodes/neurons in the network "
        !write(*,*) param%nvalues_KAN
            
        ! Structural Fingerprint setup information

        read(1,*) param%description
        read(1,*) param%atomtype
        read(1,*) param%nenv

        allocate(param%envtypes(param%nenv))

        read(1,*) param%envtypes(:)
        read(1,*) param%rc_min
        read(1,*) param%rc_max
        read(1,*) param%sftype
        read(1,*) param%nsf
        read(1,*) param%nsfparam

        allocate(param%sf(param%nsf), param%sfparam(param%nsfparam, param%nsf), &
            param%sfenv(2, param%nsf), param%sfval_min( param%nsf), &
            param%sfval_max( param%nsf), param%sfval_avg( param%nsf), param%sfval_cov( param%nsf))

        read(1,*)  param%sf(:)
        read(1,*)  param%sfparam(:,:)
        read(1,*)  param%sfenv(:,:)
        read(1,*)  param%neval
        read(1,*)  param%sfval_min
        read(1,*)  param%sfval_max
        read(1,*)  param%sfval_avg
        read(1,*)  param%sfval_cov

        !write(*,*) "an optional description from the setup file "
        !write(*,*)  param%description
        !write(*,*) "species of the central atom"
        !write(*,*)  param%atomtype
        !write(*,*) "number of different surrounding species"
        !write(*,*)  param%nenv
        !write(*,*) "species of the surrounding atoms "
        !write(*,*)  param%envtypes(:)

        !write(*,*) "minimal interaction radius and max. cutoff "
        !write(*,*)  param%rc_min
        !write(*,*)  param%rc_max
        !write(*,*) "basis function type (e.g. Behler2011)    "
        !write(*,*)  param%sftype

        if (trim(param%sftype) == "Chebyshev") then
            write(*,*) "This basis function "//trim(param%sftype)//" is supported for KAN"
            param%sftype = "ChebyshevKAN"
        else
            write(*,*) "Error: This basis function "//trim(param%sftype)//" is not supported for KAN"
            stop
        end if
        
        !radial_Rc = 6.5  radial_N = 20 angular_Rc = 5.  angular_N = 6
        param%r_Rc = param%sfparam(1,1)
        param%r_nc = int(param%sfparam(2,1))
        param%a_Rc = param%sfparam(3,1)
        param%a_nc = int(param%sfparam(4,1))
        !write(*,*) "radial_Rc,radial_N,angular_Rc,angular_N"
        !write(*,*) param%r_Rc,param%r_nc,param%a_Rc,param%a_nc


        !write(*,*) "number of structural fingerprint basis functions  "
        !write(*,*)  param%nsf
        !write(*,*) "the max. number of parameters of a basis function"
        !write(*,*)  param%nsfparam
        !write(*,*) "function kind of the particular basis type   "
        !write(*,*)  param%sf(:)
        !write(*,*) "i-th parameter of the j-th basis function"
        !write(*,*) "               i <= nsfparam"
        !write(*,*)  param%sfparam(:,:)
        !write(*,*) "i-th environment species for j-th basis function"
        !write(*,*)  param%sfenv(:,:)
        !write(*,*) "number of evaluations   "
        !write(*,*)  param%neval
        !write(*,*) "lowest so far encountered value of the i-th SF "
        !write(*,*)  param%sfval_min
        !write(*,*) "largest so far encountered value of the i-th SF "
        !write(*,*)  param%sfval_max
        !write(*,*) "current average value of the i-th symm. function "
        !write(*,*)  param%sfval_avg
        !write(*,*) "current covariance of the i-th symm. function "
        !write(*,*)  param%sfval_cov
    

            ! Trainset information

        read(1,*)  param%file
        read(1,*)  param%normalized
        read(1,*)  param%scale
        read(1,*)  param%shift
        read(1,*)  param%ntypes

        allocate( param%type_names( param%ntypes),  param%E_atom( param%ntypes))

        !write(*,*) " name of the corresponding training set file "
        !write(*,*)  param%file
        !write(*,*) ".true., if the input and output values have been  normalized to the interval [-1,1] ('read' mode only)"
        !write(*,*)  param%normalized
        !write(*,*) "energy scaling factor used for the normalization "
        !write(*,*)  param%scale
        !write(*,*) "atomic energy shift used for energy normalization "
        !write(*,*)  param%shift
        !write(*,*) "number of atomic species in the training set  "
        !write(*,*)  param%ntypes

        read(1,*)  param%type_names(:)
        !write(*,*) "name of i-th atomic species    "
        !write(*,*)  param%type_names(:)

        read(1,*)  param%E_atom(:)
        !write(*,*) "E_atom"
        
        !write(*,*)  param%E_atom(:)

        read(1,*)  param%natomtot
        read(1,*)  param%nstrucs
        read(1,*)  param%E_min,  param%E_max,  param%E_avg




        !write(*,*) "natomtot"
        !write(*,*)  param%natomtot
        !write(*,*) "nstrucs"
        !write(*,*)  param%nstrucs
        !write(*,*) "E_min,E_max,E_avg"
        !write(*,*)  param%E_min,  param%E_max,  param%E_avg
   



        close(unit = 1)

    end function



    subroutine finalize()

        implicit none
    
    end subroutine finalize
    
    subroutine initialize(infile,outfile,npoints,frombinary)
        implicit none
        character(len=*), intent(out) :: infile, outfile
        integer :: iarg, nargs
        character(len=100) :: arg
        integer,intent(out)::npoints
        integer :: conversion_success
        logical,intent(out)::frombinary
    
        nargs = command_argument_count()
        if (nargs < 1) then
           write(0,*) "Error: No input file provided."
           call finalize()
           stop
        end if

        
    
        infile = ' '
        outfile = ' '
        
        iarg = 1
        call get_command_argument(iarg, value=arg)
        infile = trim(arg)
        iarg = 2
        call get_command_argument(iarg, value=arg)
        outfile = trim(arg)
        if (nargs == 3 .or. nargs == 4) then
            iarg = 3
            call get_command_argument(iarg, value=arg)
            read(arg, *, iostat=conversion_success) npoints
        else
            npoints = 1000
        end if

        if (nargs == 4) then
            iarg = 4
            call get_command_argument(iarg, value=arg)
            if (trim(arg) == "--from-binary") then
                frombinary = .true.
            else
                frombinary = .false.
            end if
        else
            frombinary = .false.
        end if



        return


    end

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

    !--------------------------------------------------------------------!
  !                   Evaluate Chebyshev polynomials                   !
  !--------------------------------------------------------------------!

  function chebyshev_polynomial(r, r0, r1, n) result(T)
    ! Arguments:
    !    r        function argument
    !    r0, r1   the Chebyshev polynomials will be rescaled from [-1,1]
    !             to the interval [r0,r1]
    !    n        maximum polynomial order
    ! Returns:
    !    T(i)  with i=1,n+1  where T(i) is the Chebyshev polynomial of
    !    order (i-1)
    !
    ! The Chebyshev polynomials obey the following recurrence relation:
    !    T[0](x) = 1
    !    T[1](x) = x
    !    T[n+1](x) = 2x T[n](x) - T[n-1](x)

    implicit none

    double precision, intent(in)     :: r, r0, r1
    integer,          intent(in)     :: n
    double precision, dimension(n+1) :: T

    integer          :: i
    double precision :: x

    x = (2.0d0*r - r0 - r1)/(r1 - r0)

    T(1) = 1.0d0
    if (n > 0) then
       T(2) = x
       do i = 3, n+1
          T(i) = 2.0d0*x*T(i-1) - T(i-2)
       end do
    end if

  end function chebyshev_polynomial

    !--------------------------------------------------------------------!
!This subroutine comes from aenet-PyTorch. 
subroutine bin2ascii(infile, outfile)
    implicit none

    character(len=*), intent(in)  :: infile, outfile

    integer              :: nlayers, nnodesmax, Wsize, nvalues
    integer, allocatable :: nnodes(:), fun(:), iw(:), iv(:)
    real*8, allocatable  :: W(:)

    character(len=1024)           :: description
    character(len=100)            :: sftype
    character(len=2)              :: atomtype
    character(len=2), allocatable :: envtypes(:)
    real*8                        :: rc_min, rc_max
    integer                       :: nsf, nsfparam, neval, nenv
    integer, allocatable          :: sf(:), sfenv(:,:)
    real*8, allocatable           :: sfparam(:,:), sfval_min(:), sfval_max(:), sfval_avg(:), sfval_cov(:)

    character(len=1024)           :: file
    logical                       :: normalized
    real*8                        :: scale, shift, E_min, E_max, E_avg
    integer                       :: ntypes, natomtot, nstrucs
    character(len=2), allocatable :: type_names(:)
    real*8, allocatable           :: E_atom(:)


    open(unit = 1, action = "read", status = "old", file = infile, form = "unformatted")
    open(unit = 2, action = "write", status = "replace", file = outfile)


    ! Network information
    read(1) nlayers
    read(1) nnodesmax
    read(1) Wsize
    read(1) nvalues

    allocate(nnodes(nlayers), fun(nlayers-1), iw(nlayers), iv(nlayers), W(Wsize))

    read(1) nnodes(:)
    read(1) fun(:)
    read(1) iw(:)
    read(1) iv(:)
    read(1) W(:)


    write(2,*) nlayers
    write(2,*) nnodesmax
    write(2,*) Wsize
    write(2,*) nvalues
    write(2,*) nnodes(:)
    write(2,*) fun(:)
    write(2,*) iw(:)
    write(2,*) iv(:)
    write(2,*) W(:)

    deallocate(nnodes, fun, iw, iv, W)



    ! Structural Fingerprint setup information
    read(1) description
    read(1) atomtype
    read(1) nenv

    allocate(envtypes(nenv))

    read(1) envtypes(:)
    read(1) rc_min
    read(1) rc_max
    read(1) sftype
    read(1) nsf
    read(1) nsfparam

    allocate(sf(nsf), sfparam(nsfparam,nsf), sfenv(2,nsf), sfval_min(nsf), &
             sfval_max(nsf), sfval_avg(nsf), sfval_cov(nsf))

    read(1) sf(:)
    read(1) sfparam(:,:)
    read(1) sfenv(:,:)
    read(1) neval
    read(1) sfval_min
    read(1) sfval_max
    read(1) sfval_avg
    read(1) sfval_cov



    write(2,*) description
    write(2,*) atomtype
    write(2,*) nenv
    write(2,"(a)") envtypes(:)
    write(2,*) rc_min
    write(2,*) rc_max
    write(2,*) sftype
    write(2,*) nsf
    write(2,*) nsfparam
    write(2,*) sf(:)
    write(2,*) sfparam(:,:)
    write(2,*) sfenv(:,:)
    write(2,*) neval
    write(2,*) sfval_min
    write(2,*) sfval_max
    write(2,*) sfval_avg
    write(2,*) sfval_cov

    deallocate(sf, sfparam, sfenv, sfval_min, sfval_max, sfval_avg, sfval_cov)



    ! Trainset information
    read(1) file
    read(1) normalized
    read(1) scale
    read(1) shift
    read(1) ntypes

    allocate(type_names(ntypes), E_atom(ntypes))

    read(1) type_names(:)
    read(1) E_atom(:)
    read(1) natomtot
    read(1) nstrucs
    read(1) E_min, E_max, E_avg


    write(2,*) file
    write(2,*) normalized
    write(2,*) scale
    write(2,*) shift
    write(2,*) ntypes
    write(2,'(A, 2X, A,2X,A,2X,A,2X,A,2X)') type_names(:)
    write(2,*) E_atom(:)
    write(2,*) natomtot
    write(2,*) nstrucs
    write(2,*) E_min, E_max, E_avg

    deallocate(type_names, E_atom)


    close(unit = 1)
    close(unit = 2)

end subroutine bin2ascii

subroutine makeKANdescriptor(infile,outfile,npoints,frombinary)
    implicit none
    character(len=*),intent(in) :: infile,outfile
    logical,intent(in)::frombinary
    integer,intent(in)::npoints
    character(len=1024) :: infile_ascii
    type(Model_parameters)::param

    write(*,*) "input file is ",trim(infile)

    if (frombinary) then
        infile_ascii = trim(infile)//".ascii"
        !write(*,*) "Binary file is converted to ascii file"
        !write(*,*) trim(infile_ascii)
        call bin2ascii(infile, infile_ascii)
    else
        infile_ascii = trim(infile)
    end if
    if (frombinary) then
        write(*,*) "Binary file "//trim(infile)//" was converted to ascii file "//trim(infile_ascii)
    end if

    param = Model_parameters(infile_ascii)
    call construct_KAN_descriptor(param,npoints)

    call write_model(param,outfile)


    write(*,*) "output file is ",trim(outfile)
    write(*,*) "npoints is ",npoints

end subroutine makeKANdescriptor 

end module makeKAN



