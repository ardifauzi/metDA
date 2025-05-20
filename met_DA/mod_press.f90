module pressure

    use params
    use mod_grid

    implicit none
    private

    public  :: pre2D_nest

    contains


    subroutine pre2D_nest(grids,num_grids,itt,dt)

        implicit none

        type(Grid), intent(inout) :: grids(:)        
        integer, intent(in)     :: itt,num_grids
        real(DP), intent(in)    :: dt
    
        ! Constants
        real(DP), parameter :: P0 = 1.7 * 100.0D0           
        real(DP), parameter :: L0 = 150000          
        real(DP), parameter :: V0 = 55.42D0        
        real(DP), parameter :: phi_p = 44.0D0      
        real(DP), parameter :: omega = 0.00005D0        

        ! Variables
        real(DP), allocatable   :: xx(:,:), yy(:,:)
        real(DP)    :: lon_origin, lat_origin
        real(DP)    :: displacement, time, L_half,a,b,c
        real(DP),allocatable  :: P(:),sigmax(:),sigmay(:),x(:),y(:)
        integer     :: i, j, g, nx, ny,k,nga
        real(DP), allocatable :: pres(:,:), lon(:), lat(:)
        ! coordinate of starting point pressure
        real(DP), parameter :: Y0 = 47.0000 
        real(DP), parameter :: X0 = -12.0000

        ! Field parameters
        nga = 5
        allocate(x(nga),y(nga),P(nga),sigmax(nga),sigmay(nga))

        P(1) = 200
        x(1) = X0 + (-350000) * (cos(phi_p*degrad)/radi) * (1/degrad)
        y(1) = Y0 + (-350000) * (sin(phi_p*degrad)/radi*cos(phi_p*degrad)) * (1/degrad)
        sigmay(1) = 170000
        P(2) = -100
        x(2) = X0 + (-65000) * (cos(phi_p*degrad)/radi) * (1/degrad)
        y(2) = Y0 + (-65000) * (sin(phi_p*degrad)/radi*cos(phi_p*degrad)) * (1/degrad)
        sigmay(2) = 200000
        P(3) = 80
        x(3) = X0 + (400000) * (cos(phi_p*degrad)/radi) * (1/degrad)
        y(3) = Y0 + (400000) * (sin(phi_p*degrad)/radi*cos(phi_p*degrad)) * (1/degrad)
        sigmay(3) = 150000
        P(4) = 50
        x(4) = X0 + (-750000) * (cos(phi_p*degrad)/radi) * (1/degrad)
        y(4) = Y0 + (-750000) * (sin(phi_p*degrad)/radi*cos(phi_p*degrad)) * (1/degrad)
        sigmay(4) = 150000
        P(5) = 80
        x(5) = X0 + (-970000) * (cos(phi_p*degrad)/radi) * (1/degrad)
        y(5) = Y0 + (-970000) * (sin(phi_p*degrad)/radi*cos(phi_p*degrad)) * (1/degrad)
        sigmay(5) = 90000

        do g = 1, num_grids
            nx = grids(g)%nx
            ny = grids(g)%ny
            
            allocate(pres(nx,ny))
            allocate(xx(nx, ny), yy(nx, ny))
            allocate(lon(nx),lat(ny))
 
            ! Define longitude and latitude grids
            do i = 1, nx
                grids(g)%lon(i) = (grids(g)%x0 + (i - 1) * grids(g)%dx/3600) * degrad  ! Longitude: top to bottom
            enddo
            do j = 1, ny
                grids(g)%lat(j) = (grids(g)%y0 - (j - 1) *  grids(g)%dy/3600) * degrad  ! Latitude: left to right
            enddo
            lon = grids(g)%lon
            lat = grids(g)%lat

            ! Time loop
            time = itt * dt

            ! Compute the pressure field
            if (mso == 0) then
                pres = 0.0d0
                do k = 1, nga
                    lon_origin = x(k) * degrad
                    lat_origin = y(k) * degrad
                    do i = 1, nx
                        do j = 1, ny
                            xx(i, j) = radi * cos(lat_origin) * &
                                    (lon(i) - lon_origin) 
                            yy(i, j) = radi * (lat(j) - lat_origin)
                            displacement = xx(i,j) * sin(phi_p * degrad) + &
                                        yy(i,j) * cos(phi_p * degrad) - &
                                        V0 * time
                            pres(i, j) = pres(i, j) + P(k) * exp(-((displacement / (sigmay(k) / 2.0D0))**2))
                        enddo
                    enddo
                enddo
            endif
            grids(g)%pres(:,:,2) = pres
            deallocate(xx,yy)
            deallocate(pres)
            deallocate(lon,lat)
        
        enddo
        
    end subroutine

end module