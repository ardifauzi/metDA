module boundary

    use params
    use mod_grid

    implicit none

    public :: bndry

    contains

    subroutine bndry(is,ie,js,je,ubx1,ubx2,uby1,uby2,grids)

        type(Grid), intent(inout)   :: grids(:)
        integer, intent(in)         :: is,ie,js,je
        real(DP), intent(out)       :: ubx1(grids(1)%ny), ubx2(grids(1)%ny), uby1(grids(1)%nx), uby2(grids(1)%nx)
        integer                     :: ib, jb

        ! west boundary
        ib=is
        do jb=  js, je
            if(grids(1)%ddx(ib,jb) <= 0.0) then
                ubx1(jb) = 0.0
            else
                ubx1(jb) = sqrt(g0*grids(1)%ddx(ib,jb))
            end if
        enddo
        ! east boundary
        ib=ie
        do jb=  js, je
            if(grids(1)%ddx(ib,jb) <= 0.0) then
                ubx2(jb) = 0.0
            else
                ubx2(jb) = sqrt(g0*grids(1)%ddx(ib,jb))
            end if
        enddo
        ! north boundary
        jb=js
        do ib=  is, ie
            if(grids(1)%ddy(ib,jb) <= 0.0) then
                uby1(ib) = 0.0
            else
                uby1(ib) = sqrt(g0*grids(1)%ddy(ib,jb))
            end if
        enddo
        ! south boundary
        jb=je
        do ib=  is, ie
            if(grids(1)%ddy(ib,jb) <= 0.0) then
                uby2(ib) = 0.0
            else
                uby2(ib) = sqrt(g0*grids(1)%ddy(ib,jb))
            end if
        enddo

    end subroutine
    
    end module