module subroutines 

    use params
    use mod_grid

    implicit none
    private

    public  :: outsea
    public  :: frise_meteo
    public  :: hxy6_optimized
    public  :: exchan
    public  :: movie
    public  :: set_weight
    public  :: vxy3_par
    public  :: ZMAX

    contains

    subroutine frise_meteo(niz,njz,pres,hc)

        implicit none
        integer, intent(in)         :: niz,njz
        real(DP), intent(in)        :: pres(niz,njz,2)
        real(DP), intent(inout)     :: hc(niz,njz,2)
        integer                     :: ii,jj

        do jj=1,njz
            do ii=1,niz
                hc(ii,jj,1) = hc(ii,jj,1) + (- pres(ii,jj,2) / (rhow*g0))
            enddo
        enddo

        end subroutine


    subroutine outsea (is,ie,js,je,ubx1, ubx2, uby1, uby2,grids,id)

        implicit none
        type(Grid), intent(inout)   :: grids(:)
        real(DP), intent(in)        :: ubx1(grids(1)%ny),ubx2(grids(1)%ny), uby1(grids(1)%nx), uby2(grids(1)%nx)
        integer, intent(in)         :: is,ie,js,je,id
        integer                     :: ib, jb
        real(DP)                    :: theta,sint,dtds

        dtds = dt / (grids(id)%dth(1)*radi)
        ib=  is
        do jb=  js, je
            theta = grids(id)%th0(1) + (jb-js) * grids(id)%dth(1)
            sint = sin(theta)
            grids(id)%hz(ib,jb,2)=grids(id)%hz(ib,jb,1)+dtds/sint*ubx1(jb)*(grids(id)%hz(ib+1,jb,1)-grids(id)%hz(ib,jb,1))
        enddo

        ib=  ie
        do jb=  js, je
            theta = grids(id)%th0(1) + (jb-js) * grids(id)%dth(1)
            sint = sin(theta)
            grids(id)%hz(ib,jb,2)=grids(id)%hz(ib,jb,1)+dtds/sint*ubx2(jb)*(grids(id)%hz(ib-1,jb,1)-grids(id)%hz(ib,jb,1))
        enddo

        jb=  js
        do ib=  is, ie
            grids(id)%hz(ib,jb,2)=grids(id)%hz(ib,jb,1)+dtds*uby1(ib)*(grids(id)%hz(ib,jb+1,1)-grids(id)%hz(ib,jb,1))
        enddo

        jb=  je
        do ib=  is, ie
            grids(id)%hz(ib,jb,2)=grids(id)%hz(ib,jb,1)+dtds*uby2(ib)*(grids(id)%hz(ib,jb-1,1)-grids(id)%hz(ib,jb,1))
        enddo

    end subroutine
    
    subroutine vxy3_par(ni, nj, vx, vy, hz, is, ie, js, je, dx, dy, dt, th0, dth, pres)

        integer, intent(in) :: ni, nj, is, ie, js, je
        real(DP), intent(in) :: dt, th0, dth
        real(DP), intent(inout) :: vx(ni, nj, 2), vy(ni, nj, 2)
        real(DP), intent(in) :: pres(ni, nj, 2), hz(ni, nj, 2)
        real(DP), intent(in) :: dx(ni, nj), dy(ni, nj)
        real(DP) :: ddx(ni, nj), ddy(ni, nj)
        real(DP) :: dtds
        real(DP) :: theta, sint, dvdx, vvy, dvdy, adv, vsqr, absv, bfc, fric, vvx
        integer :: i, j

        dtds = dt / (dth * radi)

        !$omp parallel do private(i) shared(ddx, dx, hz)
        do j= js, je
            ddx(is,j) = dx(is,j) + hz(is,j,1)
            ddx(ie+1,j) = dx(ie+1,j) + hz(ie,j,1)
            do i=is+1, ie
                if(vx(i,j,1).ge.0) then
                    ddx(i,j) = dx(i,j) + hz(i-1,j,1)
                else
                    ddx(i,j) = dx(i,j) + hz(i,j,1)
                endif
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do private(j) shared(ddy, dy, hz)
        do i=is, ie
            ddy(i,js) = dy(i,js) + hz(i,js,1)
            ddy(i,je+1) = dy(i,je+1) + hz(i,je,1)
            do j= js+1, je
                if(vy(i,j,1).ge.0) then
                    ddy(i,j) = dy(i,j) + hz(i,j-1,1)
                else
                    ddy(i,j) = dy(i,j) + hz(i,j,1)
                endif
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do private(i, theta, sint, dvdx, vvy, dvdy, adv, vsqr, absv, bfc, fric) &
        !$omp&shared(vx, vy, ddx, hz, pres)
        do j= js+1, je-1
            theta = th0 + (j - js) * dth
            sint = sin(theta)
            do i= is+1, ie
                if(dx(i,j).le.0.0 .or. ddx(i,j).le.0.10) then
                    vx(i,j,2) = 0.0
                else
                    if(vx(i,j,1).ge.0) then
                        dvdx = vx(i,j,1) - vx(i-1,j,1)
                    else
                        dvdx = vx(i+1,j,1) - vx(i,j,1)
                    endif
                    vvy = (vy(i-1,j,1) + vy(i-1,j+1,1) + vy(i,j,1) + vy(i,j+1,1)) / 4.0
                    if(vvy.ge.0) then
                        dvdy = vx(i,j,1) - vx(i,j-1,1)
                    else
                        dvdy = vx(i,j+1,1) - vx(i,j,1)
                    endif
                    adv = dtds * (vx(i,j,1) * dvdx / sint + vvy * dvdy)
                    vsqr = vx(i,j,1)**2 + vvy**2
                    if(vsqr.gt.0) then
                        absv = sqrt(vsqr)
                    else 
                        absv = 0.0
                    endif
                    if (cf.ge.0) then
                        bfc = cf
                    else
                        bfc = cf**2 * g0 * ddx(i,j)**(-1./3)
                    endif
                    fric = dt * bfc * vx(i,j,1) * absv / ddx(i,j)
                    vx(i,j,2) = vx(i,j,1) - g0*dtds/sint * (hz(i,j,1) - hz(i-1,j,1)) - adv - fric &
                                - dtds * (pres(i, j, 2) - pres(i - 1, j, 2)) / rhow / sint
                endif
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do private(i, j, theta, sint, dvdy, vvx, dvdx, adv, vsqr, absv, bfc, fric) &
        !$omp&shared(vy, vx, ddy, hz, pres)
        do j= js+1, je
            theta = th0 + (j - js) * dth
            sint = sin(theta)
            do i= is+1, ie-1
                if(dy(i,j).le.0.0 .or. ddy(i,j).le.0.10) then
                    vy(i,j,2) = 0.0
                else
                    if(vy(i,j,1).ge.0) then
                        dvdy = vy(i,j,1) - vy(i,j-1,1)
                    else
                        dvdy = vy(i,j+1,1) - vy(i,j,1)
                    endif
                    vvx = (vx(i,j-1,1) + vx(i+1,j-1,1) + vx(i,j,1) + vx(i+1,j,1)) / 4.0
                    if(vvx.ge.0) then
                        dvdx = vy(i,j,1) - vy(i-1,j,1)
                    else
                        dvdx = vy(i+1,j,1) - vy(i,j,1)
                    endif
                    adv = dtds * (vvx * dvdx / sint + vy(i,j,1) * dvdy)
                    vy(i,j,2) = vy(i,j,1) - g0*dtds * (hz(i,j,1) - hz(i,j-1,1)) - adv - fric &
                                - dtds * (pres(i, j, 2) - pres(i, j - 1, 2)) / rhow
                endif
            enddo
        enddo
        !$omp end parallel do

        end subroutine
        
    subroutine hxy6_optimized(nix, njx, vx, vy, hz, is, ie, js, je, dz, dx, dy, dtx, th0, dth)
        use omp_lib
        implicit none
        integer, intent(in)     :: nix, njx, is, ie, js, je
        real(DP), intent(in)    :: th0, dth, dtx
        real(DP), intent(in)    :: vx(nix, njx, 2), vy(nix, njx, 2), dx(nix, njx), dy(nix, njx), dz(nix, njx)
        real(DP), intent(inout) :: hz(nix, njx, 2)
        real(DP)                :: dtds, sint1, sint2, theta, dd
        real(DP)                :: ddx(nix, njx), ddy(nix, njx)
        integer                 :: i, j
    
        dtds = dtx / (dth * radi)
    
        !$omp parallel do private(i, j) shared(ddx, dx, hz, vx, is, ie, js, je)
        do j = js, je
            ddx(is, j) = dx(is, j) + hz(is, j, 1)
            ddx(ie + 1, j) = dx(ie + 1, j) + hz(ie, j, 1)
            do i = is + 1, ie
                ddx(i, j) = dx(i, j) + hz(i - merge(1, 0, vx(i, j, 2) >= 0), j, 1)
            enddo
        enddo
        !$omp end parallel do
    
        !$omp parallel do private(i, j) shared(ddy, dy, hz, vy, is, ie, js, je)
        do i = is, ie
            ddy(i, js) = dy(i, js) + hz(i, js, 1)
            ddy(i, je + 1) = dy(i, je + 1) + hz(i, je, 1)
            do j = js + 1, je
                ddy(i, j) = dy(i, j) + hz(i, j - merge(1, 0, vy(i, j, 2) >= 0), 1)
            enddo
        enddo
        !$omp end parallel do
    
        !$omp parallel do private(i, j) shared(ddx, ddy, is, ie, js, je)
        do j = js, je + 1
            do i = is, ie + 1
                if (ddx(i, j) < 0.0) ddx(i, j) = 0.0
                if (ddy(i, j) < 0.0) ddy(i, j) = 0.0
            enddo
        enddo
        !$omp end parallel do
    
        !$omp parallel do private(i, j, theta, sint1, sint2, dd) shared(hz, vx, vy, ddx, ddy, dz, is, ie, js, je, th0, dth, dtds)
        do j = js, je
            theta = th0 + (j - js) * dth
            sint1 = sin(theta)
            sint2 = sin(theta + dth)
            do i = is, ie
                hz(i, j, 2) = hz(i, j, 1) - dtds / sint1 * (vx(i + 1, j, 2) * ddx(i + 1, j) - vx(i, j, 2) * ddx(i, j) &
                    & + vy(i, j + 1, 2) * ddy(i, j + 1) * sint2 - vy(i, j, 2) * ddy(i, j) * sint1)
    
                dd = hz(i, j, 2) + dz(i, j)
                if (dd < 1.0e-6) then
                    dd = 0.0
                    hz(i, j, 2) = dd - dz(i, j)
                end if
            enddo
        enddo
        !$omp end parallel do

    end subroutine hxy6_optimized      
        
    subroutine exchan ( nix, njx, qx, qy, hz)

        integer, intent(in)     :: nix,njx
        real(DP), intent(inout) :: qx(nix,njx,2), qy(nix,njx,2), hz(nix,njx,2)
        integer                 :: i,j

        do j= 1, njx
            do i= 1, nix
                qx(i,j,1) = qx(i,j,2)
                qy(i,j,1) = qy(i,j,2)
                hz(i,j,1) = hz(i,j,2)
            enddo
        enddo

        return
        end subroutine

    SUBROUTINE movie(KK,KD,nix,njx,hz,ka,fname)

        integer, intent(in)     :: kk,kd,nix,njx,ka
        real(DP), intent(in)    :: hz(nix,njx,2)
        integer                 :: kt,i,j
        character(*), intent(in)   :: fname
        CHARACTER               :: NAME*50

        KT=KK/KD 
        WRITE(NAME,'(a, i4)') trim(fname), KT + ka*1000
        OPEN(10,FILE=trim(NAME))
        DO I=1,nix
            WRITE(10,'(10000f9.3)') (hz(I,J,2),J=1,njx)
        enddo
        CLOSE(10)

        END subroutine

    real function qqx ( q, h1, h2, d, dt, theta, dth )
        real(DP), intent(in)    :: q, h1, h2, d, theta, dth, dt
        real(DP)                :: gts, sint

        gts = g0 * dt / (dth * radi) 
        sint = sin (theta)
        if(d.le.0.0) then
            qqx=0.0
        else
            qqx = q -gts/sint * ( h1 -h2 )
        end if
        return
        end function

    real function qqy ( q, h1, h2, d, dt, dth )
        real(DP), intent(in)    :: q, h1, h2, d, dth, dt
        real(DP)                :: gts

        gts = g0 * dt / (dth * radi) 
        if(d.le.0.0) then
            qqy=0.0
        else
            qqy = q -gts  *( h1 -h2 )
        end if
        return
        end function

    subroutine set_weight(ww,nn,no,grids,ist,jst)
        
        type(Grid), intent(in)  :: grids(:)    
        integer, intent(in)     :: nn,no,ist(no),jst(no)
        real(DP), intent(out)   :: ww(nn,no)
        real(DP), allocatable   :: mu_bgo(:,:)
        real(DP), allocatable   :: mu_boo(:,:)
        real(DP)                :: mat(No,No)
        integer                 :: ii, jj, i, j, ig, nx, ny
        real(DP)                :: dist
        real(DP)                :: RHO = 1.     !< Ratio between obs and bg error
        real(DP)                :: RR  = 20000. !< Cutoff distance of error covariance (m)
        real(DP)                :: dx, dy
    
        allocate( mu_bgo(nn,no))
        allocate( mu_boo(no,no ))
        mu_bgo(:,:) = 0.0
        mu_boo(:,:) = 0.0
        nx = grids(1)%nx
        ny = grids(1)%ny
        dx = secrad * grids(1)%dx
        dy = secrad * grids(1)%dy
    
        !! Estimate background error between numerical grid and station
        do ig=1, Nx*Ny
            do i=1, No
    
            !! Gaussian correlation function
            ii = mod( ig, Nx )
            jj = (ig - ii)/Nx + 1
            call get_distance(ist(i), jst(i), ii, jj, dist, dx, dy)
            mu_bgo(ig,i) = exp( - (dist/RR)**2 )
            end do
        end do
    
    
        !! Estimate background error between stations
        do j=1, No
            do i = 1, No
            !! Gaussian correlation function
            call get_distance( ist(i), jst(i), ist(j), jst(j), dist, dx, dy )
            mu_boo(i,j) = exp( - (dist/RR)**2 )
            end do
        end do
    
    
        !! Calculate inverse matrix for obtaining weight matrix
        do j=1, No
            do i=1, No
                mat(i,j) = mu_boo(i,j)
                if(i==j) mat(i,j) = mat(i,j) + rho**2
            end do
        end do
    
        !! invert weight vector
        do ig=1, 3*nx*ny
            call matrix__gs(no, mat, mu_bgo(ig,:), ww(ig,:))
        end do
        do ig=1, 3*nx*ny
            do i = 1, no
            if (ww(ig,i) < 1.0d-7) then
                ww(ig,i) = 0.0d0
            endif
            enddo
        enddo
    
        deallocate( mu_bgo, mu_boo )
    
    end subroutine set_weight

    subroutine get_distance( i0, j0, i1, j1, dist, dx, dy )

        implicit none
        integer, intent(in)     :: i0, j0, i1, j1
        real(DP), intent(out)   :: dist
        real(DP), intent(in)    :: dx, dy
        real(DP)                :: lat0, lon0, lat1, lon1
        real(DP)                :: dlat, dlon, a, c

        ! Convert grid indices to radians
        lat0 = i0 * dy
        lon0 = j0 * dx
        lat1 = i1 * dy
        lon1 = j1 * dx

        ! Differences in latitude and longitude
        dlat = lat1 - lat0
        dlon = lon1 - lon0

        ! Haversine formula 
        a = sin(dlat / 2.0_DP)**2 + cos(lat0) * cos(lat1) * sin(dlon / 2.0_DP)**2
        c = 2.0_DP * atan2(sqrt(a), sqrt(1.0_DP - a))

        dist = radi * c  ! Distance in meters

    end subroutine get_distance

    subroutine matrix__gs( n, m, v, a )

        integer,  intent(in)  :: n       !< model size
        real(DP), intent(in)  :: m(n,n)  !< coefficient matrix
        real(DP), intent(in)  :: v(n)    !< values
        real(DP), intent(out) :: a(n)    !< answers
        real(DP), parameter   :: TOL = 1e-6
        real(DP) :: a0(n)
        integer :: iter
        integer :: i, j
        real(DP) :: m_max, wk

        a(:) = 0
        a0(:) = 0
        iter = 0
        m_max = maxval(abs(m))

        do
            do i=1, n
                wk = 0
                do j=1, i-1
                    wk = wk + m(i,j) * a(j)
                end do
                do j=i+1, n
                    wk = wk + m(i,j) * a(j)
                end do
                a(i) = (v(i) - wk) / m(i,i)
            end do
            if( maxval(abs(a0 - a))/m_max < TOL ) exit
            a0(:) = a(:)
            iter = iter + 1
            if(iter > 10000) then
            write(*,*) "does not converge"
            exit
            end if
        end do

    end subroutine matrix__gs
    
    SUBROUTINE ZMAX(ni,nj,hz,ZMX)
        
        integer, intent(in)     :: ni,nj
        real(DP), intent(in)    :: hz(ni,nj,2)
        real(DP), intent(inout)   :: zmx(ni,nj)
        integer                 :: i,j

        DO J=1,nj
            DO I=1,ni
                IF(hz(I,J,2).GT.ZMX(I,J)) ZMX(I,J)=hz(I,J,2)
            enddo
        enddo
  
        RETURN
        END
end module

    