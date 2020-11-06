	include "noises/NOISES.F"
	include "sub2s/detector_position.F"
	include "sub2s/source_direction.F"
	include "sub2s/HF2_GR.F"
	include "sub2s/response.F"
	include "sub2s/di.F"
	include "sub2s/INVERSE_FISHER.F"
	include "sub2s/Fisher_LIGO.F"
	include "sub2s/Delta_Omegas_cal.F"
	include "sub2s/dL.F"
	include "sub2s/t_c.F"

	program main
	    implicit none
	    double precision, parameter :: PI = 3.141592653589793
	    integer          :: i, j, k
	    double precision :: para(1:9)
	    double precision :: DO1, DO2, DO3, DO4
	    double precision :: DOI1, DOI2, DOI3, DOI4
	    integer          ::  N_detector
	    double precision :: m1, m2
	    double precision :: DlndL, DlndLI
	    integer :: imax, n = 9
	    double precision :: Fij_s(1:7 ,1:9, 1:9)
	    double precision :: Fij(1:9, 1:9), cov(1:9, 1:9)
	    double precision :: de(1:7 ,1:4)
	    double precision, allocatable :: fbase(:, :), Nbase(:, :)
	    integer :: counter, stat
	    double precision :: SNR_s(7), SNR, SNRI
	    double precision :: z0, dL, dL_c
	    double precision :: omega1, omega2
	    integer :: qt0(1:6)
	    double precision :: t_c


	    z0 = 0.01D0
	    call dL_Func(dL, z0)
	    dL_c = 1.0D0 * dL / 1000.0D0

!!!!!!!     The time of the GW event.
	    qt0 = (/2017, 8, 17, 12, 41, 4/)
	    call t_c_Func(qt0, t_c)



!	    detector data (/Latitude,Longitude,Orientation,The angle between two arms/)

	    de(1, :) = (/30.56, -90.77, 243.0, 90.0/)     !LIGO Livingston
	    de(2, :) = (/46.45, -119.41, 171.8, 90.0/)    !LIGO Handford
	    de(3, :) = (/43.63, 10.5, 116.5, 90.0/)       !VIRGO
	    de(4, :) = (/36.25, 137.18, 0.0, 90.0/)       !KAGRA
	    de(5, :) = (/19.09, 74.05, 0.0, 90.0/)        !LIGO India

	    imax = 95841
	    allocate(fbase(1:3, 1:imax))
	    allocate(Nbase(1:3, 1:imax))
	    open(unit=10, file="noises/L.dat", Status="OLD")
	    open(unit=20, file="noises/H.dat", Status="OLD")
	    open(unit=30, file="noises/V.dat", Status="OLD")
	    do counter = 1, imax, 1
	        read(10, *) fbase(1, counter), Nbase(1, counter)
	        read(20, *) fbase(2, counter), Nbase(2, counter)
	        read(30, *) fbase(3, counter), Nbase(3, counter)
	    end do
	    close(10)
	    close(20)
	    close(30)


	    N_detector=3
	    para(1)  =  197.45  !!! alpha, degree
	    para(2)  =  -23.381   !!! delta, degree
	    para(3)  =  0.1   !!! varphi, degree
	    para(4)  =  28.0   !!! iota, degree
	    m1 = (2.26 + 1.36) / 2.0
	    m2 = (1.36 + 0.86) / 2.0
	    para(8)  =  (m1 + m2) * (1 + z0)                !!! M, M_sun
	    para(9)  =  m1 * m2 / ((m1 + m2)**2)   !!! eta
	    para(5)  =  t_c           !!! t_c: second
	    para(6)  =  0.0          !!! psi_c
	    para(7)  =  log(1000.0) !!! ln(d_L/Mpc)
!!!!! calculate Fisher matrix in stationary case
	    Fij = 0.0
	    Fij_s = 0.0
	    SNR = 0.0
	    SNR_s = 0.0
	    do k = 1, 3
	        call Fisher_LIGO_stationary(de(k,1),de(k,2),de(k,3),de(k,4),  ! input: (degree, GW detector)
     $	                         para,                   ! nine parameters
     $                           fbase(k,:),Nbase(k,:),imax,   ! input for noise
     $                           Fij_s(k,:,:),SNR_s(k))                   ! output
	        do i = 1,9
	            do j = 1,9
	                Fij(i,j) = Fij(i,j)+Fij_s(k,i,j)
	            end do
	        end do
	        SNR = SNR + SNR_s(k) * SNR_s(k)
	    end do
	    call Delta_Omega(Fij,para(2),DO1,DO2,DO3,DlndL)
!	    write(*, *) SNR_s
	    SNR = sqrt(SNR)
!!!!!     calculate Fisher matrix in general case
	    Fij = 0.0
	    Fij_s = 0.0
	    SNRI = 0.0
	    SNR_s = 0.0
	    do k = 1, 3
	        call Fisher_LIGO_general(de(k,1),de(k,2),de(k,3),de(k,4),  ! input: (degree, GW detector)
     $	                          para,                   ! nine parameters
     $                            fbase(k,:),Nbase(k,:),imax,   ! input for noise
     $                            Fij_s(k,:,:),SNR_s(k))                   ! output
	        do i = 1,9
	            do j = 1,9
	                Fij(i,j) = Fij(i,j)+Fij_s(k,i,j)
	            end do
	        end do
	        SNRI = SNRI + SNR_s(k) * SNR_s(k)
	    end do
!	    write(*, *) SNR_s
	    call Delta_Omega(Fij,para(2),DOI1,DOI2,DOI3,DlndLI)
	    SNRI = sqrt(SNRI)
	    DO1 = DO1 * dL_c * dL_c
	    DO2 = DO1 * dL_c * dL_c
	    DO3 = DO3 * dL_c * dL_c
	    DOI1 = DOI1 * dL_c * dL_c
	    DOI2 = DOI2 * dL_c * dL_c
	    DOI3 = DOI3 * dL_c * dL_c
	    SNR = 1.0D0 * SNR / dL_c
	    SNRI = 1.0D0 * SNRI / dL_c
	    DlndL = dlndL * dL_c
	    DlndLI = DlndLI * dL_c
	    omega1 = DO1 * -1.0D0 * log(0.1D0)
	    omega2 = DOI1 * -1.0D0 * log(0.1D0)
	    write(*, *) omega2, DlndLI
!!!!        Calculate the covariance matrix
	    call INVERSEMATRIX(Cov, Fij, n)
	    do i = 1, 9
	        do j = 1, 9
	            Cov(i, j) = Cov(i, j) * dL_c * dL_c
	        end do
	    end do
	    write(*, *) Cov(1, 1), Cov(1, 2), Cov(2, 2), Cov(7, 7)
	    write(*, *) SNR, SNRI
	end
