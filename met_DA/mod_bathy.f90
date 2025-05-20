module bathy

    use params
    use mod_grid

    implicit none
    private

    public :: bath_new

    contains

    subroutine bath_new(grids,p)
        implicit none
        type(Grid), intent(inout) :: grids(:)

        integer, intent(out) :: p(8, 35)
        integer :: g, i, j
        character(len=1000)    :: fname
    
        do g = 1, size(grids)-1
            write(fname,'(a)')grids(g)%bathdir
            open(21, file=trim(fname), status='old')
            do i = 1, grids(g)%nx
                read(21, *) (grids(g)%ddz(i, j), j = 1, grids(g)%ny)
            end do
            close(21)

            do i = 1, grids(g)%nx
                grids(g)%lon(i) = (grids(g)%x0 + (i - 1) * grids(g)%dx/3600) * degrad  ! Longitude: top to bottom
            end do
            do j = 1, grids(g)%ny
                grids(g)%lat(j) = (grids(g)%y0 - (j - 1) * grids(g)%dy/3600) * degrad  ! Latitude: left to right
            end do
            write(*, *) "Bathymetry data loaded for grid ID:", grids(g)%id
    
            ! Process the coarse grid (g == 1)
            if (g == 1) then
                grids(g)%th0 = (90.0 - grids(g)%y0) * 60.0 * minrad
                grids(g)%dth = secrad * grids(g)%dx
                call ztoxy(grids(g)%nx, grids(g)%ny, grids(g)%ddz, grids(g)%ddx, grids(g)%ddy)
                write(*, *) "Coarse grid bathymetry processed for grid ID:", grids(g)%id
            else
                ! Process finer grids (g > 1)
                p(1, g-1) = nint(((grids(g)%x0 - grids(grids(g)%parent_id)%x0) * 60.0) / (grids(grids(g)%parent_id)%dx / 60.0))
                p(2, g-1) = nint(p(1, g-1) + (grids(g)%nx - 1) / 3.0 + 1)
                p(3, g-1) = nint(((grids(grids(g)%parent_id)%y0 - grids(g)%y0) * 60.0) / (grids(grids(g)%parent_id)%dx / 60.0))
                p(4, g-1) = nint(p(3, g-1) + (grids(g)%ny - 1) / 3.0 + 1)
                p(5, g-1) = 1
                p(6, g-1) = grids(g)%nx - 1
                p(7, g-1) = 1
                p(8, g-1) = grids(g)%ny - 1
    
                grids(g)%th0 = (90.0 - grids(g)%y0) * 60.0 * minrad
                grids(g)%dth = secrad* grids(g)%dx 
                call ztoxy(grids(g)%nx, grids(g)%ny, grids(g)%ddz, grids(g)%ddx, grids(g)%ddy)
                write(*, *) "Finer grid bathymetry processed for grid ID:", grids(g)%id
            end if
        end do
    end subroutine
         

end module