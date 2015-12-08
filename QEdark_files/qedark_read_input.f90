!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Adrian Soto
! September 2015
! Stony Brook University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read QEdark input data from infile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE qedark_read_input(&
     infile, calculation_mode, restart, &
     nbndval, nbndcond, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     num_mx, mx_NU, & 
     er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     num_q_bins, deltaq, &
     scissor_correction, scissorgap)
  
  USE kinds,                       ONLY: DP
  USE klist,                       ONLY: nelec
  USE wvfct,                       ONLY: nbnd
  USE noncollin_module,            ONLY: noncolin

  IMPLICIT NONE

  CHARACTER(len=*) :: infile

  CHARACTER(len=20) :: calculation_mode
  LOGICAL, INTENT(OUT) :: restart
  
  INTEGER, INTENT(OUT) :: nbndval
  INTEGER, INTENT(OUT) :: nbndcond

  REAL(DP), INTENT(OUT) :: vearth_SI     ! Earth's average velocity around the galactic center in SI units  
  REAL(DP), INTENT(OUT) :: vesc_SI       ! Earth's escape velocity in SI units                              
  REAL(DP), INTENT(OUT) :: v0_SI         ! DM typical velocity in SI                                        
  REAL(DP), INTENT(OUT) :: deltav_SI(4)  ! velocity change due to Earth's rotation. 4 points: spring, summer, fall, winter

  
  INTEGER, PARAMETER :: max_num_mx=99
  INTEGER, INTENT(OUT) :: num_mx     
  REAL(DP), INTENT(OUT) :: mx_NU(max_num_mx)  
  

  INTEGER, INTENT(OUT) :: er_bin_type              ! 0 for no bins. 1 for linear. 2 for exponential.
  INTEGER, INTENT(OUT) :: num_er_bins              ! Number of recoil energy bins          
  REAL(DP), INTENT(OUT) :: ermax_NU
  REAL(DP), INTENT(OUT) :: er_binsize 

  INTEGER, INTENT(OUT) :: num_q_bins
  REAL(DP), INTENT(OUT) :: deltaq

  LOGICAL, INTENT(OUT)  :: scissor_correction 
  REAL(DP), INTENT(OUT) :: scissorgap              ! True (experimental) gap in eV


  INTEGER :: open_status

  NAMELIST / dm_parameters / restart, &
       calculation_mode, &
       nbndval, nbndcond, &
       vearth_SI, vesc_SI, v0_SI, deltav_SI, &
       num_mx, mx_NU, &
       er_bin_type, num_er_bins, ermax_NU, er_binsize, &
       num_q_bins, deltaq, &
       scissor_correction, scissorgap



  ! Default values
  restart = .false.
  calculation_mode = 'f2'
  
  IF ( noncolin .eqv. .false.) THEN
     nbndval = nelec/2
  ELSE
     nbndval = nelec
  ENDIF
  nbndcond = nbnd - nbndval
  
  vearth_SI = 0.0_DP
  vesc_SI = 0.0_DP
  v0_SI = 0.0_DP
  
  deltav_SI(:) = 0.0_DP
  mx_NU(:) = 0.510998910 ! == electron rest mass
  
  er_bin_type = 0
  num_er_bins = 1
  ermax_NU = 0.0_DP
  er_binsize = 1.0_DP
  
  num_q_bins = 1
  deltaq=0.1
  
  scissor_correction = .false.
  scissorgap = 0.0_DP
  


!  OPEN (7, FILE=infile, status='old', iostat=open_status, ACTION='read')

  WRITE (*,*), "Reading DM data from ", infile

  OPEN (7, FILE = infile, iostat=open_status, ACTION='read')
  IF ( open_status /= 0 ) then
     PRINT *, 'Could not open. EXIT', infile 
     STOP
  ENDIF

  READ (UNIT=7, NML=dm_parameters)

!  WRITE(*, NML=dm_parameters)
    
  CLOSE(7)

END SUBROUTINE qedark_read_input
