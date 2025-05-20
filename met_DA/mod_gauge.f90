module gauge

    use params
    use mod_grid

    implicit none
    private

    public  :: obsgauge
    public  :: dagauge
    public  :: presgauge
    public  :: tgwrt
    public  :: tgage
    public  :: pgwrt
    public  :: pgage
    public  :: dagage

    contains
    
    subroutine obsgauge(ist, jst, kfine,grids)
        
        type(Grid), intent(in)      :: grids(:)
        integer, intent(out)        :: ist(ntg), jst(ntg), kfine(ntg)
        integer                     :: i
        character(len=10)           :: sname(ntg)

        open(21,file='tg.list')
        read(21,*)
        do i = 1, ntg
            read(21,'(a10,i6,i6,i6)') sname(i),ist(i),jst(i),kfine(i)
            if (ist(i).gt.1.or.ist(i).lt.grids(kfine(i)+1)%nx.and. &
                jst(i).gt.1.or.jst(i).lt.grids(kfine(i)+1)%ny) then
                write(*,*) sname(i),ist(i),jst(i), kfine(i)
            else
                write(*,'(a,a,a,i6,a,i6,a,i3)') "Warning!! check ", sname(i),"ist= ",ist(i),"jst= ", jst(i),"region= ", kfine(i)
            end if
        enddo
        close(21)
        write(*,*)  'tide stations are read!!'
        
        end subroutine

    subroutine presgauge(ist, jst, kfine,grids)

        type(Grid), intent(in)      :: grids(:)
        integer, intent(out)        :: ist(ntw), jst(ntw), kfine(ntw)
        integer                     :: i
        character(len=10)           :: sname(ntw)

        open(21,file='pg.list')
        read(21,*)
        do i = 1, ntw
            read(21,'(a10,i6,i6,i6)') sname(i),ist(i),jst(i),kfine(i)
            if (ist(i).gt.1.or.ist(i).lt.grids(kfine(i)+1)%nx.and. &
                jst(i).gt.1.or.jst(i).lt.grids(kfine(i)+1)%ny) then
                write(*,*) sname(i),ist(i),jst(i), kfine(i)
            else
                write(*,'(a,a,a,i6,a,i6,a,i3)') "Warning!! check ", sname(i),"ist= ",ist(i),"jst= ", jst(i),"region= ", kfine(i)
            end if
        enddo
        close(21)
        write(*,*)  'weather stations are read!!'
        
        return
        end subroutine
    
    subroutine dagauge(ist, jst, kfine,grids)
    
        type(Grid), intent(in)      :: grids(:)
        integer, intent(out)        :: ist(ntd), jst(ntd), kfine(ntd)
        integer                     :: i
        character(len=10)           :: sname(ntd)

        open(21,file='tg_da.list')
        read(21,*)
        do i = 1, ntd
            read(21,'(a10,i6,i6,i6)') sname(i),ist(i),jst(i),kfine(i)
            if (ist(i).gt.1.or.ist(i).lt.grids(kfine(i)+1)%nx.and. &
                jst(i).gt.1.or.jst(i).lt.grids(kfine(i)+1)%ny) then
                write(*,*) sname(i),ist(i),jst(i), kfine(i)
            else
                write(*,'(a,a,a,i6,a,i6,a,i3)') "Warning!! check ", sname(i),"ist= ",ist(i),"jst= ", jst(i),"region= ", kfine(i)
            end if
        enddo
        close(21)
        write(*,*)  'da stations are read!!'
        
        end subroutine

    subroutine tgage(it,grids,htg,ist,jst,kfine,ntg)

        type(Grid), intent(in)      :: grids(:)
        real(DP), intent(out)       :: htg(72000,500)
        integer, intent(in)         :: ist(ntg), jst(ntg), kfine(ntg),ntg
        integer                     :: i,it

        do i=1,ntg
            htg(it,i) = grids(kfine(i))%hz(ist(i),jst(i),2)
        enddo

    end subroutine

    subroutine dagage(grids,htg,ist,jst,id)
        
        implicit none

        type(Grid), intent(in)      :: grids(:)
        real(DP), intent(out)       :: htg(ntd)
        integer, intent(in)         :: ist(ntd), jst(ntd),id
        integer                     :: i

        do i=1,ntd
            htg(i) = grids(id)%hz(ist(i),jst(i),2)
        enddo
        

    end subroutine

    subroutine pgage(it,grids,htg,ist,jst,kfine)

        type(Grid), intent(in)      :: grids(:)
        real(DP), intent(out)       :: htg(72000,500)
        integer, intent(in)         :: ist(ntw), jst(ntw), kfine(ntw)
        integer                     :: i,it

    do i=1,ntw
        htg(it,i) = grids(kfine(i))%pres(ist(i),jst(i),2)
    enddo

    end subroutine

    subroutine tgwrt(npnt,htg,ntg,char)

        real(DP), intent(in)    :: htg(72000,500)
        integer, intent(in)     :: npnt,ntg
        character(*), intent(in):: char
        integer                 :: jt,i
        character(len=20)       :: fname

        WRITE(fname,'(a)') trim(char)
        open(10,file=trim(fname),status='unknown')
        do jt=1,npnt
            write(10,'(500f14.6)') (htg(jt,i),i=1,ntg)
        enddo
        close(10)

        end subroutine
    
    subroutine pgwrt(npnt,htg)

        real(DP), intent(in)    :: htg(24400,500)
        integer, intent(in)     :: npnt
        integer                 :: jt,i

        open(10,file='pressure.dat',form='formatted')

        do jt=1,npnt
            write(10,'(500f14.6)') (htg(jt,i),i=1,ntw)
        enddo
        close(10)

        end subroutine

end module