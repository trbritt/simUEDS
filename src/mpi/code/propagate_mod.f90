module propagate_mod

  use, intrinsic :: iso_fortran_env
  use rf_mod

  implicit none

  contains
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !+
  subroutine init_vay(tstep, ex, ey, ez, hx, hy, hz, locator, acceleration, VY_p_i1, VY_gamma_i1, ptcls, nraysp, chrgpermacro, maxrayp, &
                      xmin, ymin, zmin, dx, dy, dz, &
                      ilo, ihi, jlo, jhi, klo, khi, &
                      ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl)
    implicit none
! PREALLOCATION MUST BE DONE IN MAIN ROUTINE FOR GRID ETC

    real(dp), intent(in) :: tstep
    integer :: n
    integer :: ilo, ihi,jlo, jhi,klo, khi, ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, nraysp, maxrayp
    real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: hx,hy,hz,ex,ey,ez

    real(dp) , intent(in) :: xmin, ymin, zmin, dx, dy, dz
    real(dp), intent(in) :: chrgpermacro
    real(dp) :: dt
    real(dp) :: qm
    real(dp) :: qmt2
    real(dp), dimension(7,maxrayp), intent(in) :: ptcls !n1=7
    real(dp), dimension(1:3) :: p_prime, p, tau, t, beta_h
    real(dp) :: gamma2_prime, gamma_h, tau2, sigma, u_star
    real(dp), dimension(:,:), allocatable :: grid
    integer, dimension(:), allocatable :: locator

    real(dp), dimension(:,:), allocatable :: acceleration
    real(dp), dimension(:,:), allocatable :: VY_p_i1
    real(dp), dimension(:), allocatable :: VY_gamma_i1
    real(dp), dimension(:), allocatable :: E_h, B_h

    call generate_grid(grid, ilo, ihi, jlo, jhi, klo, khi, ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, xmin, ymin, zmin, dx, dy, dz)

    dt  = tstep
    qm = (chrgpermacro/mass) * InvRestMass
    qmt2 = qm * 0.5 * dt
    do n=1,nraysp
        call getexternalfield(n, grid, locator, maxrayp, ptcls, ex, ey, ez, hx, hy, hz, ilo, ihi, jlo, jhi, klo, khi, E_h, B_h)
        p(1) = ptcls(2,n)
        p(2) = ptcls(4,n)
        p(3) = ptcls(6,n)
        p_prime(1:3) = p(1:3) + E_h(1:3) * qmt2 / clight
        gamma2_prime = sum(p_prime**2) + 1.0d0
        tau(1:3) = B_h(1:3) * qmt2
        u_star = dot_product(p_prime, tau)
        tau2 = sum(tau**2)
        sigma = gamma2_prime - tau2
        VY_gamma_i1(n) = sqrt(0.5 * (sigma + sqrt(sigma * sigma + 4.0 * (tau2 + u_star * u_star))))
        t(1:3) = tau(1:3) / VY_gamma_i1(n)

        VY_p_i1(1:3, n) = (p_prime(1:3) + t(1:3) * dot_product(p_prime, t) + cross_product(p_prime, t)) / (1 + sum(t**2))

        gamma_h = sqrt(sum(p**2) + 1.0d0)
        beta_h(1:3) = p(1:3) / gamma_h
        acceleration(1:3, n) = qm * (cross_product(beta_h, B_h) + E_h(1:3) / clight)
    enddo
    call destroy_grid(grid)
  end subroutine init_vay
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  !+
  subroutine step_vay(time, tstep, ex, ey, ez, hx, hy, hz,  locator, acceleration, VY_p_i1, VY_gamma_i1, ptcls, nraysp, chrgpermacro, maxrayp, &
                      xmin, ymin, zmin, dx, dy, dz, &
                      ilo, ihi, jlo, jhi, klo, khi, &
                      ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl)
    ! PREALLOCATION MUST BE DONE IN MAIN ROUTINE FOR GRID ETC
    implicit none
    integer :: ilo, ihi,jlo, jhi,klo, khi,ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, nraysp, n, maxrayp
    real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: hx,hy,hz,ex,ey,ez
    real(dp) , intent(in) :: xmin, ymin, zmin, dx, dy, dz
    real(dp), dimension(7,maxrayp) :: ptcls !n1=7
    real(dp), intent(in) :: chrgpermacro
    real(dp) :: qm
    real(dp) :: qmt2
    real(dp) :: time
    real(dp), intent(in) :: tstep
    real(dp), dimension(1:3) :: p_i, beta_i, x_h, p_h, beta_h, p_prime, tau, t
    real(dp) :: t_h, gamma_h, gamma2_prime, u_star, tau2, sigma
    real(dp) :: dt
    real(dp), dimension(:,:), allocatable :: grid
    integer,  dimension(:), allocatable :: locator
    real(dp), dimension(:,:), allocatable :: acceleration
    real(dp), dimension(:,:), allocatable :: VY_p_i1
    real(dp), dimension(:), allocatable :: VY_gamma_i1
    real(dp), dimension(:), ALLOCATABLE :: E_h, B_h

    dt = tstep
    qm = (chrgpermacro/mass) * InvRestMass
    qmt2 = qm * 0.5 * dt
    call generate_grid(grid, ilo, ihi, jlo, jhi, klo, khi, ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, xmin, ymin, zmin, dx, dy, dz)

    DO n=1,nraysp
        p_i(1:3) = VY_p_i1(1:3, n)
        beta_i(1:3) = p_i(1:3) / VY_gamma_i1(n)
        t_h = time + dt
        x_h(1) = ptcls(1,n) 
        x_h(2) = ptcls(3,n) 
        x_h(3) = ptcls(5,n) 
        x_h(1:3) = x_h(1:3) + beta_i(1:3) * clight * dt
        call getexternalfield(n, grid, locator, maxrayp, ptcls, ex, ey, ez, hx, hy, hz, ilo, ihi, jlo, jhi, klo, khi, E_h, B_h)
        p_h(1:3) = VY_p_i1(1:3, n) + (E_h(1:3) / clight + cross_product(beta_i, B_h)) * qmt2
        gamma_h = sqrt(sum(p_h**2) + 1.0d0)
        beta_h(1:3) = p_h(1:3) / gamma_h
        p_prime(1:3) = p_h(1:3) + E_h(1:3) / clight * qmt2
        gamma2_prime = sum(p_prime**2) +1.0d0
        tau(1:3) = B_h(1:3) + qmt2
        u_star = dot_product(p_prime, tau)
        tau2 = sum(tau**2)
        sigma = gamma2_prime - tau2
        VY_gamma_i1(n) = sqrt(0.5 * (sigma + sqrt(sigma * sigma + 4.0 * (tau2 + u_star * u_star))))
        t(1:3) = tau(1:3) / VY_gamma_i1(n)
        VY_p_i1(1:3, n) = (p_prime(1:3) + t(1:3) * dot_product(p_prime, t) + cross_product(p_prime, t)) / (1d0+sum(t**2))

        time = t_h

        ptcls(1,n) = x_h(1)  
        ptcls(3,n) = x_h(2)  
        ptcls(5,n) = x_h(3)  
        ptcls(2,n) = p_h(1)  
        ptcls(4,n) = p_h(2)  
        ptcls(6,n) = p_h(3) 
        ptcls(7,n) = time
        acceleration(1:3, n) = (cross_product(beta_h, B_h) + E_h(1:3)/clight) * qm
    enddo
    call destroy_grid(grid)
  end subroutine step_vay
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
  subroutine destroy_propagation(acceleration, VY_p_i1, VY_gamma_i1, E_h, B_h)
    real(dp), dimension(:,:), allocatable :: acceleration
    real(dp), dimension(:,:), allocatable :: VY_p_i1
    real(dp), dimension(:), allocatable :: VY_gamma_i1
    real(dp), dimension(:), allocatable :: E_h, B_h
    deallocate(acceleration)
    deallocate(VY_p_i1)
    deallocate(VY_gamma_i1)
    deallocate(E_h, B_h)
  end subroutine
  end module