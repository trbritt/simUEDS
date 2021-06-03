module rf_mod

    use, intrinsic :: iso_fortran_env
    use fieldutils_mod

    implicit none

    real(dp), parameter :: E0_010=1.59d6
    real(dp), parameter :: B0_110=3d-3
    contains
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine TM010_gridlocalfield(time, n, maxrayp, ex, ey, ez, hx, hy, hz, ilo, ihi,jlo, jhi,klo, khi, &
                                ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl, &
                                xmin, ymin, zmin, dx, dy, dz, &
                                L, sigma_E, omega, phase, origin)
        integer, intent(in) :: ilo, ihi,jlo, jhi,klo, khi
        real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: hx,hy,hz,ex,ey,ez
        real(dp), intent(in) :: time, xmin, ymin, zmin, dx, dy, dz
        integer, intent(in) :: n, maxrayp
        integer, intent(in) :: ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl
        integer :: kk,j,i
        real(dp), intent(in) :: L, sigma_E, omega, phase !PHI IS PHASE NOT POTENTIAL HERE
        real(dp), dimension(3) :: X, origin

        !field allocations
        real(dp) :: r2, k, k2, angle, sinangle, cosangle, Eror, Bphior
        real(dp), dimension(1:N) :: Ep
        integer :: num_deriv
        num_deriv=6
        !ATTEMPT 1: use the grid locations of the spacecharge nodes only for the computation of the field,
        ! add contribution of TM010 directly to this
        do kk=lbound(ex,3),ubound(ex,3)
            do j=lbound(ex,2),ubound(ex,2)
              do i=lbound(ex,1),ubound(ex,1)
                if(j.ne.jlo_rho_gbl+(jhi_rho_gbl-jlo_rho_gbl+1)/2-1)cycle
                if(kk.ne.klo_rho_gbl+(khi_rho_gbl-klo_rho_gbl+1)/2-1)cycle
                X(1) = xmin+(i-ilo_rho_gbl)*dx
                X(2) = ymin+(j-jlo_rho_gbl)*dy
                X(3) = zmin+(kk-klo_rho_gbl)*dz

                X(1:3) = X(1:3) - origin(1:3)

                if (abs(X(3)) < L/2 + 5*sigma_E) then
                    r2 = X(1)*X(1) + X(2) * X(2)
                    k = omega / clight
                    k2 = k*k
                    angle = omega*time + phase
                    sinangle = sin(angle)
                    cosangle = cos(angle)
                    call geterfderivatives(X(3), L, sigma_E, num_deriv, Ep)
                    ez(i,j,kk) = ez(i,j,kk) + E0_010 * (64*Ep(1) + r2*(-16*k2*Ep(1) - 16*Ep(3) + r2*(k2*(k2*Ep(1) + 2*Ep(3)) + Ep(5)))) * cosangle/64.d0
                    Eror = E0_010 * cosangle * (-192*Ep(2) + r2*(24*k2*Ep(2) + 24*Ep(4) + r2*(k2*(-(k2*Ep(2)) - 2*Ep(4)) - Ep(6)))) / 384.d0
                    ex(i,j,kk) = ex(i,j,kk) + X(1)*Eror
                    ey(i,j,kk) = ey(i,j,kk) + X(2)*Eror

                    Bphior = E0_010 * sinangle * (-192*k*Ep(1) + r2*(k*(24*k2*Ep(1) + 24*Ep(3)) + k*r2*(k2*(-(k2*Ep(1)) - 2*Ep(3)) - Ep(5)))) / (384.d0*clight)
                    hx(i,j,kk) = hx(i,j,kk) - X(2)*Bphior / mu0
                    hy(i,j,kk) = hy(i,j,kk) + X(1)*Bphior / mu0
                endif
              enddo
            enddo
        enddo
        
    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine TM110_gauss_gridlocalfield(time, n, maxrayp, ex, ey, ez, hx, hy, hz, ilo, ihi,jlo, jhi,klo, khi, &
                                        ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl,&
                                        xmin, ymin, zmin, dx, dy, dz, &
                                        L, sigma_E, omega, phase, xE, origin)
    integer, intent(in) :: ilo, ihi,jlo, jhi,klo, khi
    real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: hx,hy,hz,ex,ey,ez
    real(dp), intent(in) :: time, xmin, ymin, zmin, dx, dy, dz
    integer, intent(in) :: n, maxrayp
    integer, intent(in) :: ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl
    integer :: kk,j,i
    real(dp), intent(in) :: L, sigma_E, omega, phase, xE !PHI IS PHASE NOT POTENTIAL HERE
    real(dp), dimension(3) :: X, origin


    !field allocations
    real(dp) :: A, r2, r4, k, k2, k4, x2, y2, angle, sinangle, cosangle, Eror, Bphior, z
    real(dp), dimension(1:N) :: Ep
    integer :: num_deriv
    num_deriv=6
    do kk=lbound(ex,3),ubound(ex,3)
        do j=lbound(ex,2),ubound(ex,2)
          do i=lbound(ex,1),ubound(ex,1)
            if(j.ne.jlo_rho_gbl+(jhi_rho_gbl-jlo_rho_gbl+1)/2-1)cycle
            if(kk.ne.klo_rho_gbl+(khi_rho_gbl-klo_rho_gbl+1)/2-1)cycle
            X(1) = xmin+(i-ilo_rho_gbl)*dx
            X(2) = ymin+(j-jlo_rho_gbl)*dy
            X(3) = zmin+(kk-klo_rho_gbl)*dz

            X(1:3) = X(1:3) - origin(1:3)
            ! first field
            z = X(3) - xE
            if (abs(z) < 5*sigma_E) then
                A = E0_010
                r2 = X(1)*X(1) + X(2)*X(2)
                r4 = r2*r2
                k = omega/clight
                k2 = k*k
                k4 = k2*k2
                x2 = X(1)**2
                y2 = X(2)**2

                angle = omega*time+phase
                sinangle =sin(angle)
                cosangle = cos(angle)

                call getgaussderivatives(z, sigma_E, num_deriv, Ep)
                ex(i,j,kk) = ex(i,j,kk) + A * cosangle * ((192 - 24*k2*(x2 + 3*y2) + k4*r2*(x2 + 5*y2))*Ep(1) + 6*(k2*r4 - 4*(3*x2 + y2))*Ep(3) + r2*(5*x2 + y2)*Ep(5)) / 192.d0 
                ey(i,j,kk) = ey(i,j,kk) + A * cosangle * (X(1)*X(2)*(k2*(12 - k2*r2)*Ep(1) - 12*Ep(3) + r2*Ep(5))) / 48.d0 
                ez(i,j,kk) = ez(i,j,kk) + A * cosangle * (X(1)*((192 + k2*r2*(-24 + k2*r2))*Ep(2) + r2*(2*(-12 + k2*r2)*Ep(4) + r2*Ep(6)))) / 192.d0

                hx(i,j,kk) = hx(i,j,kk) - A * sinangle * (k*X(1)*X(2)*((-12 + k2*r2)*Ep(2) + r2*Ep(4))) / (24.d0*clight*mu0)
                hy(i,j,kk) = hy(i,j,kk) + A * sinangle * (k*(X(1)-X(2))*(X(1)+X(2))*((-12 + k2*r2)*Ep(2) + r2*Ep(4))) / (48.d0*clight*mu0) 
                hz(i,j,kk) = hz(i,j,kk) - A * sinangle * (k*X(2)*((192 + k2*r2*(-24 + k2*r2))*Ep(1) + r2*(2*(-12 + k2*r2)*Ep(3) + r2*Ep(5)))) / (192.d0*clight*mu0) 
            endif

            ! second field
            z = X(3) + xE
            if (abs(z) < 5*Sigma_E) then
                A = -E0_010
                r2 = X(1)*X(1) + X(2)*X(2) 
                r4 = r2*r2 
                k  = omega/clight 
                k2 = k*k 
                k4 = k2*k2 
                x2 = X(1)*X(1) 
                y2 = X(2)*X(2) 

                angle    = omega*time+phase 
                sinangle = sin(angle)  !compiler optimizes this into sincos function
                cosangle = cos(angle) 

                !Get gaussian profile and all required derivatives
                call getgaussderivatives(z,Sigma_E, num_deriv, Ep) 

                ! Calculate fields
                ex(i,j,kk) = ex(i,j,kk) + A * cosangle * ((192 - 24*k2*(x2 + 3*y2) + k4*r2*(x2 + 5*y2))*Ep(1) + 6*(k2*r4 - 4*(3*x2 + y2))*Ep(3) + r2*(5*x2 + y2)*Ep(5)) / 192.d0 
                ey(i,j,kk) = ey(i,j,kk) + A * cosangle * (X(1)*X(2)*(k2*(12 - k2*r2)*Ep(1) - 12*Ep(3) + r2*Ep(5))) / 48.d0 
                ez(i,j,kk) = ez(i,j,kk) + A * cosangle * (X(1)*((192 + k2*r2*(-24 + k2*r2))*Ep(2) + r2*(2*(-12 + k2*r2)*Ep(4) + r2*Ep(6)))) / 192.d0 

                hx(i,j,kk) = hx(i,j,kk) - A * sinangle * (k*X(1)*X(2)*((-12 + k2*r2)*Ep(2) + r2*Ep(4))) / (24.d0*clight*mu0) 
                hy(i,j,kk) = hy(i,j,kk) + A * sinangle * (k*(X(1)-X(2))*(X(1)+X(2))*((-12 + k2*r2)*Ep(2) + r2*Ep(3))) / (48.d0*clight*mu0) 
                hz(i,j,kk) = hz(i,j,kk) - A * sinangle * (k*X(2)*((192 + k2*r2*(-24 + k2*r2))*Ep(1) + r2*(2*(-12 + k2*r2)*Ep(3) + r2*Ep(5)))) / (192.d0*clight*mu0) 
            endif
          enddo
        enddo
    enddo
    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine TM110_erf_gridlocalfield(time, n, maxrayp, ex, ey, ez, hx, hy, hz, ilo, ihi,jlo, jhi,klo, khi, &
                                        ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl, &
                                        xmin, ymin, zmin, dx, dy, dz, &
                                        L, sigma_B, omega, phase, origin)
        integer, intent(in) :: ilo, ihi,jlo, jhi,klo, khi
        real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: hx,hy,hz,ex,ey,ez
        real(dp), intent(in) :: time, xmin, ymin, zmin, dx, dy, dz
        integer, intent(in) :: n, maxrayp
        integer, intent(in) :: ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, ihi_rho_gbl, jhi_rho_gbl, khi_rho_gbl
        integer :: kk,j,i
        real(dp), intent(in) :: L, sigma_B, omega, phase !PHI IS PHASE NOT POTENTIAL HERE
        real(dp), dimension(3) :: X, origin

        !field allocations
        real(dp) :: r2, r4, k, k2, k4, x2, y2, angle, sinangle, cosangle, Eror, Bphior
        real(dp), dimension(1:N) :: Bp
        integer :: num_deriv
        num_deriv=6
        !ATTEMPT 1: use the grid locations of the spacecharge nodes only for the computation of the field,
        ! add contribution of TM010 directly to this
        do kk=lbound(ex,3),ubound(ex,3)
            do j=lbound(ex,2),ubound(ex,2)
                do i=lbound(ex,1),ubound(ex,1)
                if(j.ne.jlo_rho_gbl+(jhi_rho_gbl-jlo_rho_gbl+1)/2-1)cycle
                if(kk.ne.klo_rho_gbl+(khi_rho_gbl-klo_rho_gbl+1)/2-1)cycle
                X(1) = xmin+(i-ilo_rho_gbl)*dx
                X(2) = ymin+(j-jlo_rho_gbl)*dy
                X(3) = zmin+(kk-klo_rho_gbl)*dz

                X(1:3) = X(1:3) - origin(1:3)

                if (abs(X(3))<L/2+5*sigma_B) then
                    r2 = X(1)*X(1) + X(2)*X(2) 
                    r4 = r2*r2 
                    k  = omega/clight 
                    k2 = k*k 
                    k4 = k2*k2 
                    x2 = X(1)*X(1) 
                    y2 = X(2)*X(2) 

                    angle    = omega*time+phase 
                    sinangle = sin(angle)  ! Compiler optimizes this into sincos function
                    cosangle = cos(angle) 

                    call geterfderivatives(X(3),L,Sigma_B, num_deriv, Bp) 

                    ! // Calculate fields
                    ex(i,j,kk) = ex(i,j,kk) + B0_110 * cosangle * (clight*k*(X(1)-X(2))*(X(1)+X(2))*((-12 + k2*r2)*Bp(2) + r2*Bp(4)))/48.d0 
                    ey(i,j,kk) = ey(i,j,kk) + B0_110 * cosangle * (clight*k*X(1)*X(2)*((-12 + k2*r2)*Bp(2) + r2*Bp(4)))/24.d0 
                    ez(i,j,kk) = ez(i,j,kk) + B0_110 * cosangle * (clight*k*X(1)*((192 + k2*r2*(-24 + k2*r2))*Bp(1) + r2*(2*(-12 + k2*r2)*Bp(3) + r2*Bp(5))))/192.d0 

                    hx(i,j,kk) = hx(i,j,kk) + B0_110 * sinangle * (X(1)*X(2)*(k2*(12 - k2*r2)*Bp(1) - 12*Bp(3) + r2*Bp(5)))/(48.d0 *mu0)
                    hy(i,j,kk) = hy(i,j,kk) + B0_110 * sinangle * ((192 - 24*k2*(3*x2 + y2) + k4*r2*(5*x2 + y2))*Bp(1) + 6*(k2*r4 - 4*(x2 + 3*y2))*Bp(3) + r2*(x2 + 5*y2)*Bp(5))/(192.d0*mu0) 
                    hz(i,j,kk) = hz(i,j,kk) + B0_110 * sinangle * (X(2)*((192 + k2*r2*(-24 + k2*r2))*Bp(2) + r2*(2*(-12 + k2*r2)*Bp(4) + r2*Bp(6))))/(192.d0*mu0) 
                endif
                enddo
            enddo
        enddo

    end subroutine

end module