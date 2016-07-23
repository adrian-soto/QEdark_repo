!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Adrian Soto
!  06-10-2015
!  Stony Brook University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of the code QEdark v1.1.0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Routine(s) to evaluate the form factor as in eq (3.16)
! in Essig et.al. --  arXiV:1509.01598 and the integral of this 
! quantity multiplied by the weight eta(vmin(Ei, Ei')), needed
! to calculate the DM-electron scattering rate.
!
!
!
! This code is an extension to the pwscf code of the 
! Quantum ESPRESSO package. It has been written to work with
! QE version v.5.1.2. It may also work with other versions but this 
! is not guranteed. 
!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!           QUICK INSTALLATION GUIDE
!
! 1) Download QuantumESPRESSO v.5.1.2 from
! www.qe-forge.org/gf/download/frsrelease/185/753/espresso-5.1.2.tar.gz
!
!
! 2) Modify pwscf.f90 by including the code line
!       call formfact()
!    just above the line
!       CALL stop_run( exit_status )
!
!
! 
! 3) Modify PW/src/Makefile to include dependencies on two
!    new files. Simply add
!
!      qedark.o \
!      qedark_routines.o \
!      qedark_f2.o \
!      qedark_multimass.o \
!      qedark_read_input.o \
!
!    directly below 
!       wyckoff.o \
!
!
! 4) Copy the files qedark.f90, qedark_routines.f90, 
!    qedark_f2.f90 and qedark_multimass.f90
!    into PW/src/ directory
!
!
! 5) Run configure script in QE main directory as
!       ./configure --disable-parallel --enable-openmp
!
!
! 6) Compile by typing the following command in QE main 
!    directory
!       make pw
!
!
!

SUBROUTINE qedark()  
  !
  ! Main driver of QEdark code.
  !
  !
  USE wavefunctions_module,           ONLY: evc                  ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: igk, nbnd, npwx, et,g2kin  !, btype
  USE klist,                          ONLY: nks, nkstot, ngk, wk, xk, nelec
  USE lsda_mod,                       ONLY: nspin
  USE io_files,                       ONLY: nwordwfc, iunwfc, iunigk
  USE buffers,                        ONLY: get_buffer
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: bg, tpiba, tpiba2, omega
  

  use omp_lib

  
  IMPLICIT NONE

  CHARACTER(20), PARAMETER :: infile="dm.in"                      ! Name of qedark input file


  CHARACTER(20) :: calculation_mode
  LOGICAL :: restart_mode
  LOGICAL :: f1exists, f2exists, f3exists

  REAL(DP) :: vearth_SI                                           ! Earth's average velocity around the galactic center in SI units 
  REAL(DP) :: vesc_SI                                             ! Earth's escape velocity in SI units 
  REAL(DP) :: deltav_SI(4)                                        ! Change in Earth's velocity over 4 seasons (spring, summer, fall winter) in SI
  REAL(DP) :: v0_SI                                               ! DM typical velocity in SI


  INTEGER, PARAMETER :: max_num_mx=99
  INTEGER :: num_mx
  REAL(DP) :: mx_NU(max_num_mx)                                   ! Dark matter particle mass in NU 
  
  INTEGER :: er_bin_type                                          ! Type of bins for integrated form factors
  INTEGER :: num_er_bins                                          ! Number of bins for integrated form factors

  REAL(DP) :: ermax_NU                                            ! Max recoil energy cutoff in eV
  REAL(DP) :: er_binsize                                          ! Bin size (Ry) for linear bins with user-provided size 
  
  INTEGER :: num_q_bins
  REAL(DP) :: deltaq

  INTEGER :: nksf                                                 ! Number of k-points in formfactor calculation.
  
  INTEGER :: nbndval, nbndcond                                    ! Number of occupied and unoccupied Kohn-Sham orbitals. TO BE SET BY USER!
  
  LOGICAL :: scissor_correction 
  REAL(DP) :: scissorgap                                          ! Band gap (in eV) for scissor operator corrected band energies 


  INTEGER :: ierr


  CALL start_clock( ' QEdark ')
  print *, " "  
  print *, " "
  print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
  print *, "--- Entering QEdark v1.1.0"
  print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
  print *, " "
  print *, " "
  

  IF (nspin .ne. 1) THEN
     PRINT *, " "
     PRINT *, " Spin-polarized calculation == TRUE"
     PRINT *, " "
     !CALL errore ('QEdark', 'Form factor calculation works only for spin-unpolarized systems!', 1)
  ELSE
     PRINT *, " "
     PRINT *, " Spin-polarized calculation == FALSE"
     PRINT *, " "
  ENDIF


  ! Read parameters from file
  WRITE(*,*), " "
  CALL qedark_read_input(&
     infile, calculation_mode, restart_mode, &
     nksf, nbndval, nbndcond, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     num_mx, mx_NU, &
     Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     num_q_bins, deltaq, &
     scissor_correction, scissorgap)


  IF (num_mx > 100) &
       CALL errore('QEdark', 'Number of mass points exceeds max_num_mx.', ABS(ierr) )

  
  IF (calculation_mode == "f2") THEN
     CALL qedark_f2( &
          restart_mode, &
          nksf, nbndval, nbndcond, &
          vearth_SI, vesc_SI, v0_SI, deltav_SI, &
          Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
          num_q_bins, deltaq, &
          scissor_correction, scissorgap)

  ELSEIF (calculation_mode == "f2_3d") THEN
     CALL qedark_f2_3d( &
          restart_mode, &
          nksf, nbndval, nbndcond, &
          vearth_SI, vesc_SI, v0_SI, deltav_SI, &
          Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
          num_q_bins, deltaq, &
          scissor_correction, scissorgap)
     
  ELSEIF (calculation_mode == "multimass") THEN
     CALL qedark_multimass( &
          restart_mode, &
          nksf, nbndval, nbndcond, &
          vearth_SI, vesc_SI, v0_SI, deltav_SI, &
          num_mx, mx_NU, &
          Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
          scissor_correction, scissorgap)
     
  ELSE
     CALL errore('QEdark', 'calculation_mode not recognized. EXITING.', ABS(ierr) )

  END IF


  print *, " "
  print *, " "
  print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
  print *, "-- Exiting QEdark v1.1.0"
  print *, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
  print *, " "
  print *, " "


  CALL stop_clock( ' QEdark ')
  
END SUBROUTINE QEdark
