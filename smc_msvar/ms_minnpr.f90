module msvar_prior

  use, intrinsic :: iso_fortran_env, only: wp => real64
  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: fortress_abstract_prior, M_PI, loggampdf
  use fortress_random_t, only: fortress_random
  use fortress_linalg, only: cholesky, Kronecker, determinant, inverse
  use fortress_VAR_t, only: SimsZhaSVARPrior
  use tvt_t, only: MinnesotaPrior
  use fortress_util, only: read_array_from_file
  use logbeta, only : betaln, gamln

  implicit none

  type, public :: prior_set
     class(fortress_abstract_prior), pointer :: pr
  end type prior_set


  type, public, extends(fortress_abstract_prior) :: MSVARPrior

     character(len=144) :: prior_type
     integer :: ns_mu, ns_var, nobs, ns, nA, nF, p, constant

     
     class(prior_set), allocatable, dimension(:) :: coeff_prior

     real(wp), allocatable :: priQmu(:,:), priQvar(:,:)
     real(wp), allocatable :: priXi_mean(:), priXi_std(:)

     real(wp) :: Dnplus(6,9), Nn(9,9), M(9,6)
     
   contains

     procedure :: rvs
     procedure :: logpdf
     procedure :: to_coeff_matrices
     procedure :: AF_vec_to_phi_sigma_vec
     procedure :: phi_sigma_vec_to_AF_vec
  end type MSVARPrior

  interface MSVARPrior
     module procedure new_MSVARPrior
  end interface MSVARPrior

contains



  function logdirichletpdf(x, alpha, n)

    integer, intent(in) :: n
    real(wp), intent(in) :: x(n-1), alpha(n)

    real(wp) :: logdirichletpdf
    real(wp) :: xn
    integer :: i

    logdirichletpdf = 0.0_wp

    xn = 1.0_wp - sum(x)

    if (xn < 0.0_wp) then
       logdirichletpdf = -1000000000.0_wp
       return
    end if

    do i = 1,n-1
       logdirichletpdf = logdirichletpdf  - gamln(alpha(i)) + (alpha(i)-1.0_wp)*log(x(i))
    end do
    logdirichletpdf = logdirichletpdf - gamln(alpha(n)) + (alpha(n)-1.0_wp)*log(xn) + gamln(sum(alpha))

  end function logdirichletpdf

  function determinant_gen(mat,n) result(det)

    integer, intent(in) :: n
    real(wp), intent(in) :: mat(n,n)

    integer :: i, info, ipiv(n)

    real(wp) :: det
    real(wp) :: sgn, matcp(n,n)

    det = 1.0_wp
    sgn = 1.0_wp
    matcp = mat
    call dgetrf(n,n,matcp,n,ipiv,info)

    do i = 1,n
       det = det*matcp(i,i)
       if (ipiv(i)/=i) then
          sgn = -sgn
       end if
    end do

    det = sgn*det

  end function determinant_gen


  type(MSVARPrior) function new_MSVARPrior() result(self)

    integer :: ns_mu, ns_var

    integer :: npara, nA, nF, i
    character(len=144) :: mufile, sigmafile, qmufile, qvarfile, ximufile, xistdfile

    class(fortress_abstract_prior), allocatable, target :: prior
    class(fortress_abstract_prior), pointer :: prior_i 

    character(len=144) :: dir

    self%ns_mu = {m}
    self%ns_var = {v}
    
    self%ns = self%ns_mu * self%ns_var
    self%nobs = 3
    self%prior_type = '{prior}'
    self%p = {p}
    self%constant = {cons:d}
    self%nA = {nA}
    self%nF = {nF}

    npara = 0 
    allocate(self%coeff_prior(self%ns_mu))

    if (self%prior_type=='swz') then 
       mufile = '{mufile}'
       sigmafile = '{sigmafile}'
       do i = 1, self%ns_mu
          allocate(self%coeff_prior(i)%pr, source=SimsZhaSVARPrior(self%nA+self%nF, mufile, sigmafile))
          npara = npara + self%coeff_prior(i)%pr%npara
       end do
    elseif (self%prior_type=='rfb') then

       do i = 1,self%ns_mu
          allocate(self%coeff_prior(i)%pr, &
               source=MinnesotaPrior(self%nobs, self%p, self%constant, &
                                     {lam1}, {lam2}, {lam3}, {lam4}, {lam5}, {lamxx}, &
                                     {tau}, {ybar}, {sbar}))

          npara = npara + self%coeff_prior(i)%pr%npara
       end do

              dir = '/home/eherbst/Dropbox/var_smc_estimation/replication-code/smc_msvar/'
       open(1,file=trim(dir)//'_fortress_tmp/Dnplus.txt', status='old',action='read')
       open(2,file=trim(dir)//'_fortress_tmp/Nn.txt',status='old',action='read')
       open(3,file=trim(dir)//'_fortress_tmp/M.txt',status='old',action='read')
       do i = 1,6
          read(1,*) self%Dnplus(i,:)
       end do
       do i = 1,9
          read(2,*) self%Nn(i,:)
          read(3,*) self%M(i,:)
       end do
       close(1)
       close(2)
       close(3)

    elseif (self%prior_type=='rfb-hier') then
       ! pass 
    end if

    if (self%ns_var > 1) then
       allocate(self%priXi_mean((self%ns_var-1)*self%nobs))
       allocate(self%priXi_std((self%ns_var-1)*self%nobs))

       self%priXi_mean = 1.0d0
       self%priXi_std = 1.0d0

       !call read_array_from_file(ximufile, self%priXi_mean)
       !call read_array_from_file(xistdfile, self%priXi_std)
    end if
    npara = npara + self%nobs * (self%ns_var-1)

    allocate(self%priQmu(self%ns_mu, self%ns_mu))
    allocate(self%priQvar(self%ns_var, self%ns_var))

    qmufile = 'qmean.txt'
    qvarfile = 'qvar.txt'
    call read_array_from_file(qmufile, self%priQmu, transpose=.true.)
    call read_array_from_file(qvarfile, self%priQvar, transpose=.true.)
    npara = npara + self%ns_mu * (self%ns_mu-1) + self%ns_var * (self%ns_var-1)

    self%npara = npara
  end function new_MSVARPRior



  real(wp) function logpdf(self, para) result(lpdf)

    class(MSVARPrior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)

    integer :: i , ind0, k, sigma_ind0
    real(wp) :: a, b

    real(wp) :: phi_sigma_vec(self%coeff_prior(1)%pr%npara), A0(self%nobs,self%nobs)
    real(wp) :: log_jacobian, sigma(self%nobs, self%nobs)

    real(wp) :: temp1(self%nA, self%nobs**2), temp2(self%nA, self%nobs**2), eye(self%nobs, self%nobs)
    real(wp) :: temp3(self%nA, self%nA)
    real(wp) :: negsigmakrsigma(self%nobs**2, self%nobs**2), A0krI(self%nobs**2, self%nobs**2)

    lpdf = 0.0_wp

    ind0 = 0
    do i = 1, self%ns_mu
       if (self%prior_type=='swz') then
          lpdf = lpdf + self%coeff_prior(i)%pr%logpdf(para(ind0+1:ind0+self%coeff_prior(i)%pr%npara))
       elseif (self%prior_type=='rfb') then
          call self%AF_vec_to_phi_sigma_vec(para(ind0+1:ind0+self%coeff_prior(i)%pr%npara), phi_sigma_vec)

          sigma_ind0 = 1
          eye = 0.0_wp
          do k = 1, self%nobs
             sigma(1:k,k) = phi_sigma_vec(sigma_ind0:sigma_ind0+k-1)
             sigma_ind0 = sigma_ind0 + k
             sigma(k,1:k) = sigma(1:k,k)

             eye(k,k) = 1.0_wp
          end do

          sigma_ind0 = 1
          A0 = 0.0_wp
          do k = 1, self%nobs
             A0(1:k,k) = para(sigma_ind0:sigma_ind0+k-1)
             sigma_ind0 = sigma_ind0 + k

          end do


             
          !------------------------------------------------------------
          ! jacobian -- this can be done much faster
          !------------------------------------------------------------

          associate(ny=>self%nobs, nA=>self%nA)
          call Kronecker(ny,ny,ny,ny,ny**2,ny**2,0.0_wp,1.0_wp,-sigma,sigma,negsigmakrsigma)
          call Kronecker(ny,ny,ny,ny,ny**2,ny**2,0.0_wp,1.0_wp, A0, eye, A0krI)
          call dgemm('n','n',nA,ny**2,ny**2,1.0_wp,self%Dnplus,nA,negsigmakrsigma,ny**2,0.0_wp,temp1,nA)
          call dgemm('n','n',nA,ny**2,ny**2,2.0_wp,temp1,nA,self%Nn,ny**2,0.0_wp,temp2,nA)
          call dgemm('n','n',nA,ny**2,ny**2,1.0_wp,temp2,nA,A0krI,ny**2,0.0_wp,temp1,nA)
          call dgemm('n','n',nA,nA,ny**2,1.0_wp,temp1,nA,self%M,ny**2,0.0_wp,temp3,nA)
          log_jacobian = log(abs(determinant_gen(temp3, nA))) &
               - 1.0_wp*(ny*self%p+1)*log(abs(determinant_gen(A0,ny)))
        end associate

          lpdf= lpdf + self%coeff_prior(i)%pr%logpdf(phi_sigma_vec) + log_jacobian
       end if
       ind0 = ind0 + self%coeff_prior(i)%pr%npara
    end do

    do i = 1,self%nobs*(self%ns_var-1)
       lpdf = lpdf + loggampdf(para(ind0+i), 1.0_wp, 1.0_wp)
    end do
    ind0 = ind0 + self%nobs*(self%ns_var-1)

    do i = 1, self%ns_mu
       lpdf = lpdf + logdirichletpdf(para(ind0+1:ind0+(self%ns_mu-1)), self%priQmu(:,i), self%ns_mu)
       ind0 = ind0 + self%ns_mu-1
    end do

    do i = 1, self%ns_var
       lpdf = lpdf + logdirichletpdf(para(ind0+1:ind0+(self%ns_var-1)), self%priQvar(:,i), self%ns_var)
       ind0 = ind0 + self%ns_var-1
    end do

  end function logpdf


  function rvs(self, nsim, seed, rng) result(parasim)

    class(MSVARPrior), intent(inout) :: self

    integer, intent(in) :: nsim
    integer, optional :: seed
    type(fortress_random), optional :: rng
    real(wp) :: parasim(self%npara, nsim)
    type(fortress_random) :: use_rng


    real(wp) :: temp(nsim,1), temp_mu(nsim, self%ns_mu), temp_var(nsim, self%ns_var)
    real(wp) :: sigma_phi(self%coeff_prior(1)%pr%npara, nsim)
    integer :: i, j, ind0


    use_rng = fortress_random()

    parasim = 0.0_wp
    ind0 = 0

    do i = 1, self%ns_mu
       if (self%prior_type=='swz') then
          parasim(ind0+1:ind0+self%coeff_prior(i)%pr%npara, :) = self%coeff_prior(i)%pr%rvs(nsim)
       elseif (self%prior_type=='rfb') then

          sigma_phi = self%coeff_prior(i)%pr%rvs(nsim)
          do j = 1, nsim
             call self%phi_sigma_vec_to_AF_vec(sigma_phi(:,j), &
                  parasim(ind0+1:ind0+self%coeff_prior(i)%pr%npara, j) )
          end do
       end if
    ind0 = ind0 + self%coeff_prior(i)%pr%npara
    end do

    do i = 1,self%nobs*(self%ns_var-1)
       temp = use_rng%gamma_rvs(nsim, 1, theta=1.0_wp, k=1.0_wp)
       parasim(ind0+i, :) = temp(:,1)
    end do

    ind0 = ind0 + self%nobs*(self%ns_var-1)

    if (self%ns_mu > 1) then
       do i = 1, self%ns_mu
          do j = 1, self%ns_mu
             temp = use_rng%gamma_rvs(nsim, 1, theta=self%priQmu(j,i), k=1.0_wp)
             temp_mu(:,j) = temp(:,1)
          end do
          do j = 1, nsim
             temp_mu(j, :) = temp_mu(j, :) / sum(temp_mu(j,:))
          end do
          parasim(ind0+1:ind0+self%ns_mu-1,:) = temp_mu(:, 1:self%ns_mu-1)
          ind0 = ind0 + self%ns_mu-1
       end do
    end if

    if (self%ns_var > 1) then 
       do i = 1, self%ns_var
          do j = 1, self%ns_var
             temp = use_rng%gamma_rvs(nsim, 1, theta=self%priQvar(j,i), k=1.0_wp)
             temp_var(:,j) = temp(:,1)
          end do
          do j = 1, nsim
             temp_var(j, :) = temp_var(j, :) / sum(temp_var(j,:))
          end do
          parasim(ind0+1:ind0+self%ns_var-1,:) = temp_var(:, 1:self%ns_var-1)
          ind0 = ind0 + self%ns_var-1
       end do
    end if

  end function rvs

  subroutine AF_vec_to_phi_sigma_vec(self, AF_vec, phi_sigma_vec)
    class(MSVARPrior), intent(inout) :: self

    real(wp), intent(in) :: AF_vec(self%nA+self%nF)
    real(wp), intent(out) :: phi_sigma_vec(self%nobs*(self%nobs+1)/2+self%nF)


    real(wp) :: A(self%nobs,self%nobs), F(self%nF/self%nobs,self%nobs)
    real(wp) :: sigma(self%nobs,self%nobs), phi(self%nF/self%nobs,self%nobs)

    real(wp) :: Ai(self%nobs,self%nobs)

    integer :: info, k, j, ind0, As_ind0

    associate(n=>self%nobs, nA=>self%nA, nF=>self%nF)

      A = 0.0_wp
      As_ind0 = 1
      do j = 1,n
         A(1:j,j) = AF_vec( As_ind0 : As_ind0+j-1)
         As_ind0 = As_ind0+j 

         F(:,j) = AF_vec(nA+(nF/n)*(j-1)+1:nA+(nF/n)*j)
      end do

      
      call dgemm('n','t',n,n,n,1.0_wp,A,n,A,n,0.0_wp,sigma,n)
      call inverse(sigma,info)

      Ai = A
      call inverse(Ai,info)

      call dgemm('n','n',nF/n,n,n,1.0_wp,F,nF/n,Ai,n,0.0_wp,phi,nF/n)
      ind0 = 1
      do k = 1,n
         phi_sigma_vec(ind0:ind0+k-1) = sigma(1:k,k)
         ind0 = ind0 + k
      end do

      ind0 = n*(n+1)/2
      do k = 1,n
         phi_sigma_vec(ind0+(k-1)*nF/n+1:ind0+k*nF/n) = phi(:,k)
      end do
      
    end associate

  end subroutine AF_vec_to_phi_sigma_vec

  subroutine phi_sigma_vec_to_AF_vec(self, phi_sigma_vec, AF_vec)
    class(MSVARPrior), intent(inout) :: self

    real(wp), intent(out) :: AF_vec(self%nA+self%nF)
    real(wp), intent(in) :: phi_sigma_vec(self%nobs*(self%nobs+1)/2+self%nF)


    real(wp) :: A(self%nobs,self%nobs), F(self%nF/self%nobs,self%nobs)
    real(wp) :: sigma(self%nobs,self%nobs), phi(self%nF/self%nobs,self%nobs)

    real(wp) :: Ai(self%nobs,self%nobs)

    integer :: info, k, j, ind0, As_ind0

    associate(n=>self%nobs, nA=>self%nA, nF=>self%nF, &
         prior => self%coeff_prior(1)%pr)


      ! call phi_sigma 
      select type(prior)
      class is (MinnesotaPrior)
         call prior%para_to_sigma_phi(phi_sigma_vec, sigma, phi)
      class default
         print*,'coefficient prior misspecified'
         stop
      end select


      A = sigma

      call cholesky(A, info)
      call inverse(A, info)
      A = transpose(A)

      call dgemm('n','n',nF/n,n,n,1.0_wp, phi, nF/n, A, n, 0.0_wp, F, nF/n)
      As_ind0 = 1
      do j = 1,n
         AF_vec( As_ind0 : As_ind0+j-1) = A(1:j,j)  
         As_ind0 = As_ind0+j 

         AF_vec(nA+(nF/n)*(j-1)+1:nA+(nF/n)*j) =  F(:,j) 
      end do
      
    end associate

  end subroutine phi_sigma_vec_to_AF_vec




  subroutine to_coeff_matrices(self, para,As,Fs,Xis,Qs,RC)
    class(MSVARPrior), intent(inout) :: self
    real(wp), intent(in) :: para(self%npara)
    real(wp), intent(out) :: As(self%nobs,self%nobs,self%ns_mu), Fs(self%nF/self%nobs,self%nobs,self%ns_mu)
    real(wp), intent(out) :: Xis(self%nobs,self%nobs,self%ns_var), Qs(self%ns,self%ns)

    real(wp) :: Qmu(self%ns_mu,self%ns_mu), Qvar(self%ns_var,self%ns_var)
    integer :: RC 
    integer :: i, j, As_ind0

    associate(ns=>self%ns, &
         ns_mu=>self%ns_mu, &
         ns_var=>self%ns_var, &
         nA=>self%nA, &
         nF=>self%nF, & 
         ny=>self%nobs)


    RC = 0

    !------------------------------------------------------------
    ! conditional means
    !------------------------------------------------------------
    As(:,:,:) = 0.0_wp

    do i = 1, ns_mu
       As_ind0 = 1
       do j = 1,ny
          As(1:j,j,i) = para( (nA + nF)*(i - 1)+As_ind0 : (nA + nF)*(i - 1)+As_ind0+j-1)
          As_ind0 = As_ind0+j 

          Fs(:,j,i) = para((nA + nF)*(i-1)+nA+ (nF/ny)*(j-1) + 1:(nA + nF)*(i-1)+nA+ (nF/ny)*j)
       end do
    end do

    !------------------------------------------------------------
    ! variances
    !------------------------------------------------------------
    Xis = 0.0_wp
    do j = 1,ny
       Xis(j,j,1) = 1.0_wp
    end do

    do i = 2, ns_var
       do j = 1, ny
          Xis(j,j,i) = para(ns_mu*(nA+nF)+ny*(i-2)+j)
       end do
    end do

    !------------------------------------------------------------
    ! Qs
    !------------------------------------------------------------
    do i = 1, ns_mu
       Qmu(1:ns_mu-1,i) = para(ns_mu*(nA+nF)+(ns_var-1)*ny+(ns_mu-1)*(i-1)+1: ns_mu*(nA+nF)+(ns_var-1)*ny+(ns_mu-1)*i)
       Qmu(ns_mu,i) = 1.0_wp - sum(Qmu(1:ns_mu-1,i))

       if (any(Qmu(1:ns_mu-1,i) < 0.0_wp) .or. any(Qmu(1:ns_mu-1,i)>1.0_wp)) then
          RC = -1
       end if

       if ((Qmu(ns_mu,i) < 0.0_wp) .or. (Qmu(ns_mu,i)>1.0_wp)) then
          RC = -2
       end if

    end do

    do i = 1, ns_var
       Qvar(1:ns_var-1,i) = para(ns_mu*(nA+nF)+(ns_var-1)*ny+(ns_mu-1)*ns_mu + (ns_var-1)*(i-1)+1: &
            ns_mu*(nA+nF)+(ns_var-1)*ny+(ns_mu-1)*ns_mu + (ns_var-1)*i)
       Qvar(ns_var,i) = 1.0_wp - sum(Qvar(1:ns_var-1,i))

       if (any(Qvar(1:ns_var-1,i) < 0.0_wp) .or. any(Qvar(1:ns_var-1,i)>1.0_wp)) then
          RC = -3
       end if

       if ((Qvar(ns_var,i) < 0.0_wp) .or. (Qvar(ns_var,i)>1.0_wp)) then
          RC = -4
       end if

    end do


    call Kronecker(ns_mu,ns_mu,ns_var,ns_var,ns,ns,0.0_wp,1.0_wp,Qmu,Qvar,Qs)

    end associate
end subroutine to_coeff_matrices


end module msvar_prior


module model_t
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fortress_bayesian_model_t, only: fortress_abstract_bayesian_model
  use fortress_prior_t, only: M_PI
  use fortress_linalg, only: cholesky, Kronecker, determinant
  use msvar_prior, only: MSVARPrior

  implicit none

  type, public, extends(fortress_abstract_bayesian_model) :: model
     integer :: p, constant
     integer :: nA, nF
     integer :: ns_mu, ns_var, ns

     integer :: likT 
     real(wp), allocatable :: YY_lik(:,:), XX_lik(:,:)


   contains
     procedure :: lik
     procedure :: lik_var_ms_cond
     procedure :: lik_var_integrate
  end type model

  interface model
     module procedure new_model
  end interface model

contains

  type(model) function new_model() result(self)

    character(len=144) :: name, datafile, prior_mu_file, prior_var_file
    integer :: nobs, T, p, constant, nA, nF, i, j

    name = 'svar'
    datafile = 'data.txt'
    T = {data.shape[0]}
    nobs = {data.shape[1]}

    allocate(self%prior, source=MSVARPrior())
    
    call self%construct_model(name, datafile, self%prior%npara, nobs, T) 

    associate(prior => self%prior)
      select type(prior)
      class is (MSVARPrior)
         self%p = prior%p
         self%constant = prior%constant
         self%ns_mu = prior%ns_mu
         self%ns_var = prior%ns_var
         self%ns = prior%ns
      class default
         print*,'prior misspecified'
         stop
      end select
    end associate

    
    self%likT = self%T - self%p

    allocate(self%YY_lik(self%likT, self%nobs), &
         self%XX_lik(self%likT, self%nobs*self%p+self%constant))


    self%XX_lik = 1.0_wp
    do i = 1, self%likT
       self%YY_lik(i,:) = self%YY(:,self%p+i)
       do j = 1, self%p
          self%XX_lik(i,(j-1)*self%nobs+1:j*self%nobs) = self%YY(:,i-j+self%P)
       end do
    end do

    self%T = self%likT
  end function



  function lik(self, para, T) result(lk)
    class(model), intent(inout) :: self

  real(wp), intent(in):: para(self%npara)
  integer, intent(in), optional :: T

  real(wp) :: lk

  real(wp) :: As(self%nobs,self%nobs,self%ns_mu), Fs(self%nobs*self%p+1,self%nobs,self%ns_mu)
  real(wp) :: Xis(self%nobs,self%nobs,self%ns_var), Qs(self%ns,self%ns)
  integer :: RC

  real(wp) :: likmat(self%likT,self%ns)


    associate(prior => self%prior)
      select type(prior)
      class is (MSVARPrior)
         call prior%to_coeff_matrices(para,As,Fs,Xis,Qs,RC)
      class default
         print*,'prior misspecified'
         stop
      end select
    end associate

  if (RC < 0) then
     lk = -10000000000000.0_wp
     return
  end if

  call self%lik_var_ms_cond(As,Fs,Xis,Qs,likmat)
  call self%lik_var_integrate(likmat,Qs,lk)

end function lik


subroutine lik_var_ms_cond(self, As, Fs, Xis, Qs, likmat)

  class(model), intent(inout) :: self

  real(wp), intent(in) :: As(self%nobs,self%nobs,self%ns_mu)
  real(wp), intent(in) :: Fs(self%nobs*self%p+self%constant,self%nobs,self%ns_mu)
  real(wp), intent(in) :: Xis(self%nobs,self%nobs,self%ns_var), Qs(self%ns,self%ns)
  real(wp), intent(out) :: likmat(self%likT,self%ns)

  integer :: i,j,z
  real(wp) :: resid(self%likT,self%nobs)

  real(wp) :: AXi(self%nobs,self%nobs), FXi(self%nobs*self%p+self%constant,self%nobs)
  real(wp) :: AXiAXip(self%nobs,self%nobs)

  associate(ny=>self%nobs, p=>self%p, constant=>self%constant, T=>self%likT, &
       ns_mu=>self%ns_mu, ns_var=>self%ns_var)
  likmat = -ny/2.0*log(2.0_wp*M_PI)

  do i = 1,ns_mu
     do j = 1,ns_var

        AXi = As(:,:,i)
        FXi = Fs(:,:,i)
        ! this assumes Xi is diagonal
        do z = 1,ny
           AXi(:,z) = AXi(:,z)*Xis(z,z,j)
           FXi(:,z) = FXi(:,z)*Xis(z,z,j)
        end do

        call dgemm('n','n',T,ny,ny,1.0_wp,self%YY_lik,T,AXi,ny,0.0_wp,resid,T)
        call dgemm('n','n',T,ny,ny*p+constant,-1.0_wp,self%XX_lik,T,FXi,ny*p+constant,1.0_wp,resid,T)


        call dgemm('n','t',ny,ny,ny,1.0_wp,AXi,ny,AXi,ny,0.0_wp,AXiAXip,ny)
        likmat(:,(i-1)*ns_var+j) = likmat(:,(i-1)*ns_var+j) + 0.5_wp*log(determinant(AXiAXip,ny))
        do z = 1,ny
           likmat(:,(i-1)*ns_var+j) = likmat(:,(i-1)*ns_var+j) - 0.5_wp*resid(:,z)**2
        end do
     end do
  end do
  end associate


end subroutine lik_var_ms_cond

subroutine lik_var_integrate(self, likmat,Qs,likval)

  class(model), intent(inout) :: self

  real(wp), intent(in) :: likmat(self%T,self%ns), Qs(self%ns,self%ns)
  real(wp), intent(out) :: likval

  real(wp) :: prob(self%ns,self%likT), prob1(self%ns), prob_adj_likmat(self%ns,self%likT),prob_vec(self%ns)

  real(wp) :: Ct(self%T), sumCt

  integer :: i
  associate(ny=>self%nobs, p=>self%p, constant=>self%constant, T=>self%likT, &
       ns_mu=>self%ns_mu, ns_var=>self%ns_var, ns=>self%ns)

  Ct = maxval(likmat,dim=2)
  sumCt = sum(Ct)

  likval = 0.0_wp
  do i = 1, T

     if (i > 1) then
        prob_vec = prob(:,i-1)
     else
        prob_vec = 1.0_wp/(ns*1.0_wp)
     end if

     call dgemv('n',ns,ns,1.0_wp,Qs,ns,prob_vec,1,0.0_wp,prob1,1)

     prob_adj_likmat(:,i) = exp(likmat(i,:)-Ct(i)) * prob1
     prob(:,i) = prob_adj_likmat(:,i)/sum(prob_adj_likmat(:,i))


     likval = likval + log(sum(prob_adj_likmat(:,i)))
  end do

  likval = likval + sumCt

  ! GCC has underflow problems ... 
  if ((isnan(likval))) likval = -10000000.0d0
  
  end associate

end subroutine lik_var_integrate



end module model_t
