module mod_grid
    use params
    
    implicit none
    public  :: initialize_grids 
    public  :: display_grid_info

    type :: Grid
        integer :: id            ! Unique identifier for the grid
        integer :: level         ! Level in the hierarchy (0 for parent, 1+ for children)
        integer :: nx, ny        ! Grid dimensions
        real(DP) :: x0, y0       ! Upper-left coordinates of the grid
        real(DP) :: dx, dy       ! Grid resolution (spacing in x and y directions)
        integer :: parent_id     ! ID of the parent grid (-1 for root grids)
        integer, allocatable :: children_ids(:) ! IDs of child grids
        character(len=256) :: bathdir    ! Bathymetry file for the grid
        real(DP), allocatable :: vx(:,:,:), vy(:,:,:)   ! Averaged bathymetry (dx, dy)
        real(DP), allocatable :: hz(:,:,:)           
        real(DP), allocatable :: ddx(:,:)          
        real(DP), allocatable :: ddy(:,:)           
        real(DP), allocatable :: ddz(:,:)           
        real(DP), allocatable :: dth(:)             
        real(DP), allocatable :: th0(:)            
        real(DP), allocatable :: pres(:,:,:)             
        real(DP), allocatable :: lon(:)             
        real(DP), allocatable :: lat(:)             
    end type Grid

    type(Grid), allocatable :: grids(:) 

    contains

    subroutine initialize_grids(config_file, num_grids)
        implicit none
        character(len=*), intent(in)    :: config_file
        integer, intent(out)            :: num_grids
        integer :: id, level, nx, ny, parent_id, ios, g
        real(DP) :: x0, y0, dx, dy
        character(len=256) :: line, bathdir
        integer, allocatable :: child_counts(:)
        
        open(unit=10, file=config_file, status='old', iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error: Unable to open configuration file: ", config_file
            stop
        endif
        
        read(10, '(i3)', iostat=ios) num_grids
        allocate(grids(1:2))
        allocate(child_counts(1:num_grids))
        child_counts = 0
        
        do g = 1, num_grids
            read(10, '(A)', iostat=ios) line
            write(*,*) trim(line)
            read(line, '(I3,I3,I7,I7,I3,f12.6,f12.6,f9.5,f9.5,A100)') id, level, nx, ny, parent_id, x0, y0, dx, dy, bathdir
        
            grids(g)%id = id
            grids(g)%level = level
            grids(g)%nx = nx
            grids(g)%ny = ny
            grids(g)%parent_id = parent_id
            grids(g)%x0 = x0
            grids(g)%y0 = y0
            grids(g)%dx = dx
            grids(g)%dy = dy
            grids(g)%bathdir = trim(bathdir)
        
            allocate(grids(g)%vx(nx, ny, 2))
            allocate(grids(g)%vy(nx, ny, 2))
            allocate(grids(g)%hz(nx, ny, 2))
            allocate(grids(g)%ddx(nx, ny))
            allocate(grids(g)%ddy(nx, ny))
            allocate(grids(g)%ddz(nx, ny))
            allocate(grids(g)%dth(1))
            allocate(grids(g)%th0(1))
            allocate(grids(g)%pres(nx, ny, 2))
            allocate(grids(g)%lon(nx))
            allocate(grids(g)%lat(ny))     
   
            grids(g)%vx = 0.0
            grids(g)%vy = 0.0
            grids(g)%hz = 0.0
        
            if (parent_id > 0) then
                child_counts(parent_id) = child_counts(parent_id) + 1
            endif
        end do
        
        do g = 1, num_grids
            if (child_counts(grids(g)%id) > 0) then
                allocate(grids(g)%children_ids(child_counts(grids(g)%id)))
            endif
        end do
        
        child_counts = 0
        do g = 1, num_grids
            parent_id = grids(g)%parent_id
            if (parent_id > 0) then
                child_counts(parent_id) = child_counts(parent_id) + 1
                grids(parent_id)%children_ids(child_counts(parent_id)) = grids(g)%id
            endif
        end do
        
        close(10)
    end subroutine    
    
    subroutine display_grid_info()
        implicit none
        integer :: g
    
        write(*,*) "Grid Configuration:"
        do g = 1, size(grids)
            write(*,*) "Grid ID:", grids(g)%id
            write(*,*) "  Level:", grids(g)%level
            write(*,*) "  Dimensions (nx, ny):", grids(g)%nx, grids(g)%ny
            write(*,*) "  Upper-left (x0, y0):", grids(g)%x0, grids(g)%y0
            write(*,*) "  Resolution (dx, dy):", grids(g)%dx, grids(g)%dy
            write(*,*) "  Parent ID:", grids(g)%parent_id
            write(*,*) "  Bathymetry File:", trim(grids(g)%bathdir)
            if (allocated(grids(g)%children_ids)) then
                write(*,*) "  Children ID:", grids(g)%children_ids
            else
                write(*,*) "  Children: None"
            end if
        end do
    end subroutine    
    
end module 
