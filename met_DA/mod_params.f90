module params

    implicit none 
    public

    integer, parameter, public :: DP = selected_real_kind(9)

    real(DP), parameter :: g0 = 9.80665                 ! Gravity
    real(DP), parameter :: pi = atan(1.0) * 4           ! Pi 
    real(DP), parameter :: rhow = 1020.0D0              ! Seawater density
    real(DP), parameter :: degrad = pi / 180.0          ! Degree to radian
    real(DP), parameter :: secrad = pi / (180.0*60*60)  ! Second to radian
    real(DP), parameter :: minrad = pi / (180.0 * 60.0) ! Minute to radian
    real(DP), parameter :: radi = 6371000.0d0           ! Earth's radius
    real(DP), parameter :: cf = -0.025                  ! Manning's friction coef

    integer, parameter  :: icor     = 0                 ! Coriolis effect, 0 for with Coriolis
    real(DP), parameter :: dt       = 1.00d0            ! Computational time interval (sec)
    integer, parameter  :: tmax     = 3600*9            ! Maximum computational time (sec)
    real(DP), parameter :: gaugeout = 60.0d0            ! Output interval at gauges
    real(DP), parameter :: movout   = 60.0d0*5          ! Output interval of movies
    integer, parameter  :: ntg      = 6                 ! No. of observational gauges
    integer, parameter  :: ntd      = 11                ! No. of observational gauges
    integer, parameter  :: ntw      = 4                 ! No. of observational gauges

end module params 
