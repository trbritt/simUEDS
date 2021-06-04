module fieldutils_mod

    use, intrinsic :: iso_fortran_env
    use kdtree2_module
    use hdf5utils_mod
    implicit none
  
    ! Fortran 2008
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64
    integer, parameter :: qp = real128
    real(dp), parameter :: InvRestMass = 1.758820025d11
    real(dp), parameter :: clight=299792458.d0
    real(dp), parameter :: mu0=8.d0*asin(1.d0)*1.d-7
    real(dp), parameter :: mass = 0.510998910d6 / (clight * clight)
    real(dp),parameter :: pi=4.D0*datan(1.D0)
    contains
  
    function cross_product(a, b)
      real(dp), dimension(3) :: cross_product
      real(dp), dimension(3), intent(in) :: a, b
   
      cross_product(1) = a(2)*b(3) - a(3)*b(2)
      cross_product(2) = a(3)*b(1) - a(1)*b(3)
      cross_product(3) = a(1)*b(2) - b(1)*a(2)
    end function cross_product
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine getgaussderivatives(x, sigma, N, derivatives)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: sigma
      integer, intent(in)  :: N
      real(dp), dimension(1:N), intent(out) :: derivatives
      real(dp) :: scale, gauss
      integer i
      real(dp) :: acc

      if (N.lt.1.or.sigma.lt.0d0) then
        write (6,*) 'Bad result in getgaussderivatives.'
        stop
      endif
      scale = 1.0d0/(sqrt(2.0d0)*sigma)
      derivatives(1) = 1
      derivatives(2) = 2*x*scale
      do i=3,N
        derivatives(i) =  -2.0*(i-2.0d0) * derivatives(i-2) + 2.0*x*scale*derivatives(i-1)
      enddo

      do i=2, N
        acc = acc * (-scale)
        derivatives(i) = derivatives(i) * acc
      enddo
      gauss = exp(-x*x/(2*sigma*sigma))
      derivatives(1:N) = derivatives(1:N) * gauss

    end subroutine getgaussderivatives
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine geterfderivatives(x, L, sigma, N, derivatives)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: L
      real(dp), intent(in) :: sigma
      integer, intent(in)  :: N
      real(dp), dimension(1:N), intent(out) :: derivatives
      real(dp) :: scale
      integer i
      real(dp) :: acc, fac
      real(dp), dimension(1:N) :: erfp, erfm

      
      if (N.lt.1.or.sigma.lt.0d0) then
        write (6,*) 'Bad result in getgaussderivatives.'
        stop
      endif
      scale = 1.0d0/(sqrt(2.0d0)*sigma)
      if(L.eq.0d0) then
        call getgaussderivatives(x, sigma, N, derivatives)
      else
        erfp(1) = 0.5*erf(scale*(x+L/2)) 
        erfm(1) = 0.5*erf(scale*(x-L/2)) 
        erfp(2) = exp(-scale*scale * (x+L/2)*(x+L/2)) * scale / sqrt(pi) 
        erfm(2) = exp(-scale*scale * (x-L/2)*(x-L/2)) * scale / sqrt(pi) 
        do i=3,N
          erfp(i) = -2*scale*scale*((i-3.0)*erfp(i-2) + (x+L/2)*erfp(i-1))
          erfp(i) = -2*scale*scale*((i-3.0)*erfp(i-2) + (x-L/2)*erfp(i-1))
        enddo
        fac = 1/erf(scale*L/2)
        derivatives(1:N) = fac * (erfp(1:N) - erfm(1:N))
      endif

      
    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine getexternalfield(n, grid, locator, maxrayp, ptcls, ex, ey, ez, hx, hy, hz, ilo, ihi,jlo, jhi,klo, khi, &
                                E_h, B_h)
      integer, intent(in) :: ilo, ihi,jlo, jhi,klo, khi, maxrayp
      real(dp), dimension(ilo:ihi,jlo:jhi,klo:khi) :: hx,hy,hz,ex,ey,ez
      real(dp), dimension(:), allocatable :: E_h, B_h
      integer :: n,m
      real(dp), dimension(7, maxrayp) :: ptcls

      real(dp), allocatable :: data(:,:)
      type(kdtree2_result), allocatable :: results(:)
      real(dp), dimension(:,:), allocatable :: grid
      integer, dimension(:), allocatable, intent(in) :: locator


      type(kdtree2), pointer :: tree
      integer, dimension(3,4) :: indices
      integer :: tmp, i
      integer :: index_x, index_y, index_z

      allocate(results(4))
      grid(1, size(grid,2)) = ptcls(1, n)
      grid(2, size(grid,2)) = ptcls(3, n)
      grid(3, size(grid,2)) = ptcls(5, n)

      tree => kdtree2_create(grid, rearrange=.true., sort=.true.)
      call kdtree2_n_nearest_around_point(tree, idxin=size(grid,2), nn=4, correltime=5, results=results)
      write(*,*) "Nearest neighbors found at indices: ", results(:)%idx

      do m=1,4
        tmp = results(m)%idx
        index_z = tmp / (ihi*jhi)
        tmp = tmp - index_z*ihi*jhi
        index_y = tmp / ihi
        index_x = mod(tmp, ihi)
        ! this scheme returns indexes running from 0 to <i,j,k>hi-<i,j,k>lo, so add lower bound
        indices(1, m) = index_x + ilo
        indices(2, m) = index_y + jlo
        indices(3, m) = index_z + klo
      enddo
      ! write(*,*) "ilo = ", ilo, "ihi = ", ihi
      ! write(*,*) "jlo = ", jlo, "jhi = ", jhi
      ! write(*,*) "klo = ", klo, "khi = ", khi
      ! write(6,*) indices
      !compute average of 4 neighbors
      do m=1,4
        E_h(1) = E_h(1) + ex(indices(1,m), indices(2, m), indices(3,m))
        E_h(2) = E_h(2) + ey(indices(1,m), indices(2, m), indices(3,m))
        E_h(3) = E_h(3) + ez(indices(1,m), indices(2, m), indices(3,m))

        B_h(1) = B_h(1) + hx(indices(1,m), indices(2, m), indices(3,m))*mu0
        B_h(2) = B_h(2) + hy(indices(1,m), indices(2, m), indices(3,m))*mu0
        B_h(3) = B_h(3) + hz(indices(1,m), indices(2, m), indices(3,m))*mu0
      enddo
      call kdtree2_destroy(tree)
      E_h(:) = E_h(:) / 4
      B_h(:) = B_h(:) / 4
      ! write(6,*) E_h
      ! write(6,*) B_h
    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine generate_grid(grid, ilo, ihi,jlo, jhi,klo, khi, &
                            ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl, &
                            xmin, ymin, zmin, dx, dy, dz)
      real(dp), intent(in) :: xmin, ymin, zmin, dx, dy, dz
      integer, intent(in) :: ilo, ihi,jlo, jhi,klo, khi, ilo_rho_gbl, jlo_rho_gbl, klo_rho_gbl
      real(dp), dimension(:,:), allocatable :: grid
      integer :: num_points
      integer :: kk,j, i, accumulator=1
      num_points = (ihi-ilo) * (jhi-jlo) * (khi-klo)
      allocate(grid(3, 1:num_points+1)) !+1 so this grid will contain one blank point for data, set to 0 here
      
      ! grid(1,1:num_points) = xmin + ((/(i, i=ilo,ihi)/) - ilo_rho_gbl) * dx
      ! grid(2,1:num_points) = ymin + ((/(i, i=jlo,jhi)/) - jlo_rho_gbl) * dy
      ! grid(3,1:num_points) = zmin + ((/(i, i=klo,khi)/) - klo_rho_gbl) * dz
      do kk=klo,khi-1
        do j=jlo, jhi-1
          do i=ilo, ihi-1
            !https://stackoverflow.com/questions/7367770/how-to-flatten-or-index-3d-array-in-1d-array
            grid(1, accumulator) = xmin+(i-ilo_rho_gbl)*dx
            grid(2, accumulator) = ymin+(j-jlo_rho_gbl)*dy
            grid(3, accumulator) = zmin+(kk-klo_rho_gbl)*dz
            accumulator = accumulator + 1
          enddo
        enddo
      enddo
      grid(:, num_points+1) = 0.d0

    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine destroy_grid(grid)
      real(dp), dimension(:,:), allocatable:: grid
      deallocate(grid)
    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine destroy_locator(locator)
      integer, dimension(:), allocatable:: locator
      deallocate(locator)
    end subroutine
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !+
    subroutine generate_locator(locator, ilo, ihi,jlo, jhi,klo, khi)
      integer, intent(in) :: ilo, ihi,jlo, jhi,klo, khi
      integer, dimension(:), allocatable :: locator
      integer :: accumulator = 1, kk, j, i
      allocate(locator((ihi-ilo) * (jhi-jlo) * (khi-klo)))
      do kk=klo,khi-1
        do j=jlo, jhi-1
          do i=ilo, ihi-1
            locator(accumulator) = ((kk)*ihi*jhi) + ((j)*ihi) + (i)
            accumulator = accumulator + 1
          enddo
        enddo
      enddo
    end subroutine

end module