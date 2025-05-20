program meteo 

    use params
    use mod_grid
    use bathy
    use gauge
    use boundary
    use pressure
    use subroutines

    implicit none

    real(DP), allocatable       :: ubx1i(:), ubx2i(:), uby1i(:), uby2i(:)
    real(DP), allocatable       :: htgi(:,:),htgi2(:,:),ptgi(:,:),htgi_da(:,:)
    integer, allocatable        :: isti(:), jsti(:), kfinei(:), kfinei2(:)
    integer, allocatable        :: isti_da(:), jsti_da(:), kfinei_da(:)
    integer, allocatable        :: ipti(:), jpti(:), kpfinei(:)
    real(DP), allocatable       :: d_obs(:)
    real(DP), allocatable       :: ww(:,:), mf(:), ma(:)
    real(DP), allocatable       :: yt(:,:), yf(:) 
    real(DP), allocatable       :: eta(:,:), vx(:), vy(:),pres(:,:,:)
    real(DP), allocatable       :: zmx_syn(:,:),zmx_da(:,:)
    
    integer         :: it,isi,iei,jsi,jei,itgout,itg,itda
    real(DP)        :: itmax,ttt,prog,elapsed,remaining
    integer         :: pp(8,35),itmout,num_grids
    integer         :: start,clock_rate,clock_max,finish
    integer         :: g, NN, ni, nj,i,j
    character(12)   :: config_file = "gridfile.cfg"

    call system_clock(start,clock_rate,clock_max)
    call initialize_grids(config_file,num_grids)
    call display_grid_info()
    
    ni = grids(1)%nx
    nj = grids(1)%ny
    NN  = 3 * grids(1)%nx * grids(1)%ny
    allocate (ubx1i(grids(1)%ny), ubx2i(grids(1)%ny), uby1i(grids(1)%nx), uby2i(grids(1)%nx))
    allocate (htgi(72000,ntg),ptgi(72000,ntw),htgi_da(72000,ntd),htgi2(72000,ntd))
    allocate (isti(ntg), jsti(ntg), kfinei(ntg), kfinei2(ntg))
    allocate (isti_da(ntd), jsti_da(ntd), kfinei_da(ntd))
    allocate (ipti(ntw), jpti(ntw), kpfinei(ntw))
    allocate( ww(NN,ntd) )
    allocate( d_obs(ntd) )
    allocate( mf(NN), ma(NN) )
    allocate( yt(32400,ntd), yf(ntd) )
    allocate( eta(grids(1)%nx,grids(1)%ny), vx(grids(1)%nx*grids(1)%ny), vy(grids(1)%nx*grids(1)%ny) )
    allocate( pres(grids(1)%nx,grids(1)%ny,2) )
    allocate( zmx_syn(grids(1)%nx,grids(1)%ny),zmx_da(grids(1)%nx,grids(1)%ny) )

    itmax   = tmax/dt + 0.5
    itgout  = gaugeout / dt + 0.5
    itmout  = movout / dt + 0.5
    itg     = 0
    itda    = 0

    isi = 1
    iei = grids(1)%nx - 1
    jsi = 1
    jei = grids(1)%ny - 1

    call bath_new(grids,pp)
    call obsgauge(isti, jsti, kfinei, grids)
    call dagauge(isti_da, jsti_da, kfinei_da, grids)
    call presgauge(ipti, jpti, kpfinei, grids)
    call bndry(isi,iei,jsi,jei,ubx1i,ubx2i,uby1i,uby2i,grids)
    call set_weight(ww,nn,ntd,grids,isti_da,jsti_da)
    
    grids(2)%nx = grids(1)%nx
    grids(2)%ny = grids(1)%ny
    grids(2)%vx = grids(1)%vx
    grids(2)%vy = grids(1)%vy
    grids(2)%hz = grids(1)%hz
    grids(2)%ddx = grids(1)%ddx
    grids(2)%ddy = grids(1)%ddy
    grids(2)%ddz = grids(1)%ddz
    grids(2)%th0 = grids(1)%th0
    grids(2)%dth = grids(1)%dth
    pres(:,:,:) = 0.0d0
    kfinei2 = kfinei + 1
    
    if (grids(1)%dx*30.87/dt < sqrt(2*g0*maxval(grids(1)%ddz))) then
      write(*,*) "ERROR: stability condition is violated: ",  &
      grids(1)%dx*30.87/dt/sqrt(2*g0*maxval(grids(1)%ddz))
      stop
    else
      write(*,*) "CFL condition OK: ",  &
      grids(1)%dx*30.87/dt/sqrt(2*g0*maxval(grids(1)%ddz))
    endif

    open(21, file='tsunami_obs.dat', status='old')
    do i = 1, 10800
        read(21, *) (yt(i, j), j = 1, ntd)
    end do
    close(21)

    do it= 1, int(itmax)
     
      !! DA
      call vxy3_par(grids(1)%nx,grids(1)%ny,grids(2)%vx,grids(2)%vy,grids(2)%hz,&
                          isi,iei,jsi,jei,grids(1)%ddx,grids(1)%ddy,dt,grids(1)%th0(1),grids(1)%dth(1),pres)
      call hxy6_optimized(grids(1)%nx,grids(1)%ny,grids(2)%vx,grids(2)%vy,grids(2)%hz,&
                          isi,iei,jsi,jei,grids(1)%ddz,grids(1)%ddx,grids(2)%ddy,dt,grids(2)%th0(1),grids(2)%dth(1))
      call outsea(isi,iei,jsi,jei,ubx1i, ubx2i, uby1i, uby2i,grids,2)
 
      call dagage(grids,yf,isti_da, jsti_da, 2)

      !! Residual
      d_obs(:) = yt(it,:) - yf(:)

      !! Assimilation
      mf(     1:  grids(1)%nx * grids(1)%ny)                        = reshape( grids(2)%hz(:,:,2), [ni * nj] )
      mf(  grids(1)%nx * grids(1)%ny+1:2*grids(1)%nx * grids(1)%ny) = reshape( grids(2)%vx(:,:,2), [ni * nj])
      mf(2*grids(1)%nx * grids(1)%ny+1:3*grids(1)%nx * grids(1)%ny) = reshape( grids(2)%vy(:,:,2), [ni * nj])
     
      ma = mf + matmul( ww, d_obs )

      grids(2)%hz(:,:,2)  = reshape(ma(     1:  grids(1)%nx * grids(1)%ny), [grids(1)%nx, grids(1)%ny] ) 
      grids(2)%vx(:,:,2)  = reshape(ma(  grids(1)%nx* grids(1)%ny+1:2*grids(1)%nx * grids(1)%ny), [grids(1)%nx, grids(1)%ny] ) 
      grids(2)%vy(:,:,2)  = reshape(ma(2*grids(1)%nx* grids(1)%ny+1:3*grids(1)%nx * grids(1)%ny), [grids(1)%nx, grids(1)%ny] ) 
       
      ! Output snapshots
      if(mod(it,itmout).eq.0)  then
        call movie(it,itmout,grids(1)%nx,grids(1)%ny,grids(2)%hz,1,"zda")
      endif
      
      ! Save tsunami waveforms
      if(mod(it,itgout).eq.0)  then
        itg = itg + 1
        call tgage(itg,grids,htgi,isti, jsti, kfinei,ntg)
        call tgage(itg,grids,htgi2,isti, jsti, kfinei2,ntg)
      endif
      call tgwrt(itg,htgi,ntg,"tsunami_syn_ga.dat")
      call tgwrt(itg,htgi2,ntg,"tsunami_da_ga.dat")
      itda = itda + 1
      call tgage(itda,grids,htgi_da,isti_da, jsti_da, kfinei_da,ntd)
      call tgwrt(itda,htgi_da,ntd,"tsunami_da_obs.dat")

      call zmax(grids(1)%nx, grids(2)%ny,grids(2)%hz,zmx_da)

      ! Exchange flux and height for next step
      do g = 1, 2
        call exchan(grids(g)%nx,grids(g)%ny,grids(g)%vx,grids(g)%vy,grids(g)%hz) 
      enddo
    enddo

    open(12,file='zmax_da.dat', status= "unknown")
    do i = 1, grids(1)%nx
      write(12,'(10000F9.3)') (zmx_da(i,j), j=1,grids(1)%ny)
    enddo
    close(12)

end program 