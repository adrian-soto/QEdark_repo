!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Adrian Soto
!  17-09-2015
!  Stony Brook University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of the code QEdark v.0.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutines for QEdark code
!
!
SUBROUTINE SaveRestartParams( &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, mx_NU, &
     Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     do_scissor_correction, scissorgap, do_vacuum_level, vacuum_level)

  USE kinds,                         ONLY: DP  
  USE wvfct,                         ONLY: ecutwfc, nbnd
  USE klist,                         ONLY: nks
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN):: do_scissor_correction, do_vacuum_level
  INTEGER, INTENT(IN):: Er_bin_type, num_er_bins
  REAL(DP), INTENT(IN):: vearth_SI, vesc_SI, v0_SI, deltav_SI(3), mx_NU, &
       ermax_NU, er_binsize, scissorgap, vacuum_level

  
  CHARACTER(70) :: FMTparams

  FMTparams = '(F12.6, 2I6.2, 7ES15.7, 2I6.2, 2F12.6, L5, F12.6, L5, ES15.7)'

  OPEN(UNIT=910, FILE='params.restart', STATUS='replace', ACCESS='sequential', FORM='formatted')
  WRITE(910,FMTparams), ecutwfc, nbnd, nks, vearth_SI, vesc_SI, v0_SI, deltav_SI(1:3), & 
       mx_NU, Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     do_scissor_correction, scissorgap, do_vacuum_level, vacuum_level

  CLOSE(910)

END SUBROUTINE SaveRestartParams



SUBROUTINE SaveRestartData(ik1, ik2, num_er_bins, nmonths, ctot, cbinned)

  USE kinds,                         ONLY: DP  
  USE wvfct,                         ONLY: ecutwfc, nbnd
  USE klist,                         ONLY: nks
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ik1, ik2
  INTEGER, INTENT(IN) :: num_er_bins, nmonths
  REAL(DP), INTENT(IN) :: ctot(3, nmonths)
  REAL(DP), INTENT(IN) :: cbinned(3, nmonths, num_er_bins+1)


  
  CHARACTER(26) :: FMTkpts
  
  FMTkpts = '(2I6.2)'

  ! Save last pair of kpoints (checkpoint)
  OPEN (UNIT=911, FILE='kpts.restart', STATUS='replace', ACCESS='sequential', FORM='formatted')

  WRITE(911, FMTkpts), ik1, ik2
  CLOSE (911)
  
  ! Save last integrated formfactor
  CALL C2file("C.restart", num_er_bins, nmonths, ctot, cbinned)

END SUBROUTINE SaveRestartData




SUBROUTINE ReadRestart(num_er_bins, nmonths, &
     restartmode, do_scissor_correction, do_vacuum_level, Er_bin_type, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     mx_NU, ermax_NU, er_binsize, scissorgap, vacuum_level, &
     ik1init, ik2init, ctot, cbinned)

  USE kinds,                         ONLY: DP  
  USE wvfct,                         ONLY: ecutwfc, nbnd
  USE klist,                         ONLY: nks
   
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_er_bins, nmonths
  
  LOGICAL, INTENT(IN) :: restartmode, do_scissor_correction, do_vacuum_level
  INTEGER, INTENT(IN) :: Er_bin_type
  REAL(DP), INTENT(IN) :: vearth_SI, vesc_SI, v0_SI, deltav_SI(3), &
       mx_NU, ermax_NU, er_binsize, scissorgap, vacuum_level
  INTEGER, INTENT(OUT) :: ik1init, ik2init
  REAL(DP), INTENT(OUT) :: ctot(3, nmonths)
  REAL(DP), INTENT(OUT) :: cbinned(3, nmonths, num_er_bins+1)

  CHARACTER(26) :: FMTkpts      ! Format specifier for kpoint indices
  CHARACTER(70) :: FMTparams      ! Format specifier for restart parameters
  CHARACTER(26) :: FMTC          ! Format specifier for integrated formfactors
  INTEGER :: idmff, imonth
  CHARACTER(2) :: numbins
  

  LOGICAL :: restartmode_restart, do_scissor_correction_restart, do_vacuum_level_restart
  INTEGER :: nbnd_restart, nks_restart, Er_bin_type_restart, num_er_bins_restart
  REAL(DP) :: ecutwfc_restart, vearth_SI_restart, vesc_SI_restart, v0_SI_restart, deltav_SI_restart(3), &
       mx_NU_restart, ermax_NU_restart, er_binsize_restart, scissorgap_restart, vacuum_level_restart 
  

  !print *, "ASC-- DEBUG-- From ReadRestart(): reading integrated formfactors"
     
  ! Read restart parameters
  
  
  ! Formats
  FMTparams = '(F12.6, 2I6.2, 7ES15.7, 2I6.2, 2F12.6, L5, F12.6, L5, ES15.7)'
  FMTkpts = '(2I6.2)'
  
  WRITE (numbins, '(I3)') num_er_bins + 1 ! Write the number of bins to string variable numbins 
  FMTC = "(" //'2I6.2,' // 'E19.9E3,' // numbins // 'E19.9E3' // ")"
  
  
  ! Read
  OPEN(UNIT=901, FILE='params.restart', STATUS='old', ACTION='read', ACCESS='sequential', FORM='formatted')
  READ(901, FMTparams), &
       ecutwfc_restart, nbnd_restart, nks_restart, &
       vearth_SI_restart, vesc_SI_restart, v0_SI_restart, deltav_SI_restart, mx_NU_restart, &
     Er_bin_type_restart, num_er_bins_restart, ermax_NU_restart, er_binsize_restart, &
     do_scissor_correction_restart, scissorgap_restart, do_vacuum_level_restart, vacuum_level_restart
  CLOSE(901)
  
 

  ! Check if parameters from previous and current run are the same
   IF(do_scissor_correction_restart .neqv. do_scissor_correction) THEN
     CALL errore('ReadRestart', 'do_scissor_correction has different value than in previous run', 1)

  ELSEIF(do_vacuum_level_restart .neqv. do_vacuum_level) THEN
     CALL errore('ReadRestart', 'do_scissor_correction has different value than in previous run', 1)

  ELSEIF (nbnd_restart /= nbnd) THEN
     CALL errore ('ReadRestart', 'nbnd_restart is different than nbnd in input file', 1)

  ELSEIF (num_er_bins_restart /= num_er_bins) THEN
     CALL errore ('ReadRestart', 'num_er_bins_restart is different than num_er_bins in dm.in', 1)

  ELSEIF (nbnd_restart /= nbnd) THEN
     CALL errore ('ReadRestart', 'nbnd_restart is different than nbnd in input file', 1)

  ELSEIF (abs(vearth_SI_restart-vearth_SI) > 0.001 ) THEN
     CALL errore ('ReadRestart', 'vearth_SI is different than in previous run', 1)

  ELSEIF (abs(vesc_SI_restart-vesc_SI)/vesc_SI > 0.001 ) THEN
     CALL errore ('ReadRestart', 'vesc_SI is different than in previous run', 1)

  ELSEIF (abs(v0_SI_restart-v0_SI)/v0_SI > 0.001 ) THEN
     CALL errore ('ReadRestart', 'v0_SI is different than in previous run', 1)

  ELSEIF (abs((deltav_SI_restart(1)-deltav_SI(1))/deltav_SI(1)) > 0.001 ) THEN
     CALL errore ('ReadRestart', 'deltav_SI(1) is different than in previous run', 1)

  ELSEIF (abs((deltav_SI_restart(2)-deltav_SI(2))/deltav_SI(2)) > 0.001 ) THEN
     CALL errore ('ReadRestart', 'deltav_SI(2) is different than in previous run', 1)

  ELSEIF (abs((deltav_SI_restart(3)-deltav_SI(3))/deltav_SI(3)) > 0.001 ) THEN
     CALL errore ('ReadRestart', 'deltav_SI(3) is different than in previous run', 1)

  ELSEIF (abs(ermax_NU_restart-ermax_NU)/ermax_NU > 0.001 ) THEN
     CALL errore ('ReadRestart', 'deltav_SI is different than in previous run', 1)

  ELSEIF (abs(scissorgap_restart-scissorgap)/scissorgap > 0.001 ) THEN
     CALL errore ('ReadRestart', 'scissorgap is different than in previous run', 1)

  ELSEIF (abs(vacuum_level_restart-vacuum_level)/vacuum_level > 0.001 ) THEN
          print *, "ASC-- DEBUG-- vacuum_level_restart, vacuum_level", vacuum_level_restart,vacuum_level
   
     CALL errore ('ReadRestart', 'vacuum_level is different than in previous run', 1)
     
  ELSEIF (abs(ecutwfc_restart-ecutwfc)/ecutwfc > 0.001 ) THEN
     CALL errore ('ReadRestart', 'ecutwfc_restart is different than ecutfwc in input file', 1)
  
  ELSEIF (nks_restart /= nks) THEN
     CALL errore ('ReadRestart', 'nks is different than nks in input file', 1)
  
  ELSE
     
     
     OPEN(UNIT=902, FILE='kpts.restart', STATUS='old', ACTION='read', ACCESS='sequential', FORM='formatted')
     READ(902, FMTkpts), ik1init, ik2init
     CLOSE(902)
     
     
     OPEN(UNIT=903, FILE='C.restart', STATUS='old', ACTION='read', ACCESS='sequential', FORM='formatted')
     
     imonth=1
     DO WHILE(imonth < nmonths)
        idmff=1
        DO WHILE(idmff < 3)
           !
           ! Third column is the total integrated C. 
           ! Forth column onwards are the integrated C, classified
           ! in bins depending on the electron recoil energy.
           !
           READ(903, FMTC), &
                idmff,&
                imonth, &
                ctot(idmff, imonth), &
                cbinned(idmff, imonth, 1:num_er_bins+1)
        ENDDO
     ENDDO
     
     CLOSE(903)
     

  ENDIF

END SUBROUTINE ReadRestart


SUBROUTINE C2file(filename, num_er_bins, nmonths, ctot, cbinned)
  !
  ! Adrian Soto
  ! 12-01-2015
  ! Stony Brook University
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print output of f2 calculation
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  USE kinds,                         ONLY: DP

  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: num_er_bins, nmonths
  REAL(DP), INTENT(IN) :: ctot(3, nmonths)
  REAL(DP), INTENT(IN) :: cbinned(3, nmonths, num_er_bins+1)
  
  INTEGER :: imonth, idmff
  CHARACTER(3) :: numbins ! WARNING: max numbins is 999
  CHARACTER(26) :: FMT              ! Format specifier  



  ! Print to file
  OPEN (UNIT=188, FILE=filename, STATUS='replace', ACCESS='sequential', FORM='formatted')
  
  ! Write the number of bins to string variable numbins 
  WRITE (numbins, '(I3)') num_er_bins + 1 
  
  ! Create format descriptor
  FMT = "(" //'2I6.2,' // 'E19.9E3,' // numbins // 'E19.9E3' // ")"
  
  DO imonth=1, nmonths
     DO idmff=1, 3
        !
        ! Third column is the total integrated C. 
        ! Forth column onwards are the integrated C, classified
        ! in bins depending on the electron recoil energy.
        !
        WRITE(188, FMT), &
             idmff,&
             imonth, &
             ctot(idmff, imonth), &
             cbinned(idmff, imonth, 1:num_er_bins+1)
        
        
     ENDDO
     
  ENDDO

  CLOSE(188)

END SUBROUTINE C2file
  


SUBROUTINE C2file_onlyf2(filename, numqbins, num_er_bins, ctot)
  !
  ! Adrian Soto
  ! 09-06-2015
  ! Stony Brook University
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print output of f2 calculation
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE kinds,                         ONLY: DP

  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: num_er_bins, numqbins
  REAL(DP), INTENT(IN) :: ctot(numqbins+1, num_er_bins+1)

  
  INTEGER :: iq, iE
  CHARACTER(3) :: numbinsE ! WARNING: max numbins is 999
  CHARACTER(3) :: numbinsq ! WARNING: max numbins is 999
  CHARACTER(26) :: FMT              ! Format specifier  


  ! Print to file
  OPEN (UNIT=188, FILE=filename, STATUS='replace', ACCESS='sequential', FORM='formatted')
  
  ! Write the number of bins to string variable numbins 
  WRITE (numbinse, '(I3)') num_er_bins + 1 
  WRITE (numbinsq, '(I3)') numqbins + 1 
  
  ! Create format descriptor
  FMT = 'E19.9E3,'
  

  WRITE(188,*) numqbins, num_er_bins

  DO iE=1, num_er_bins
     DO iq=1, numqbins
        !
        !WRITE(188, FMT), ctot(iq, iE)
        WRITE(188, *) ctot(iq, iE)
     ENDDO
  ENDDO

  CLOSE(188)

  print *, " "
  print *, " Output written to ", filename
  print *, " "

END SUBROUTINE C2file_onlyf2



SUBROUTINE C2file_f2_3d(filename, dq, numqbins, dEr, num_Er_bins, ctot)
  !
  ! Adrian Soto
  ! 24-03-2016
  ! Stony Brook University
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Print output of f2 calculation
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE kinds,                         ONLY: DP

  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: num_Er_bins, numqbins
  REAL(DP), INTENT(IN):: dq, dEr
  REAL(DP), INTENT(IN):: ctot(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1)

  
  INTEGER :: iqx, iqy, iqz, iE
  CHARACTER(3) :: numbinsE ! WARNING: max numbins is 999
  CHARACTER(3) :: numbinsq ! WARNING: max numbins is 999
  CHARACTER(26) :: FMT              ! Format specifier  


  ! Print to file
  OPEN (UNIT=188, FILE=filename, STATUS='replace', ACCESS='sequential', FORM='formatted')
  
  ! Write the number of bins to string variable numbins 
  

  ! Create format descriptor for real numbers
  FMT = "(ES12.6)"  

  WRITE(188,*) numqbins, dq
  WRITE(188,*) numqbins, dq
  WRITE(188,*) numqbins, dq
  WRITE(188,*) num_Er_bins, dEr

  DO iE=1, num_er_bins
     DO iqz=1, numqbins
        DO iqy=1, numqbins
           DO iqx=1, numqbins
              WRITE(188, FMT) ctot(iqx, iqy, iqz, iE)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  CLOSE(188)

  print *, " "
  print *, " Output written to ", filename
  print *, " "

END SUBROUTINE C2file_f2_3d

  



SUBROUTINE scissor(nocc, nunocc, gap, energies)
  !
  ! Adrian Soto
  ! 23-06-2014
  ! Stony Brook University
  ! 
  ! Apply scissor operator to correct band energies:
  ! shift occupied bands down and unoccupied bands up so
  ! that the band gap matches a user-specified band gap 
  ! energy.
  !
  ! Internally, QE keeps the energies in Ry. Even though the band
  ! gap input is taken in eV.
  
  
  USE kinds,                         ONLY: DP
  USE klist,                         ONLY: nks

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nocc
  INTEGER, INTENT(IN) :: nunocc
  REAL(DP), INTENT(IN) :: gap
  REAL(DP), INTENT(INOUT) :: energies(nocc+nunocc, nks)

  REAL(DP) :: maxocc, minunocc                         ! Minimum energy gap found
 ! INTEGER :: ioccgap, iunoccgap, ikoccgap, ikunoccgap  ! Indices for gap finding. Uncomment the lines involving these variables if interested.
  INTEGER :: ik, ibnd                                  ! Loop indices
  REAL(DP) :: shift                                    ! (wanted_gap - DFT_gap)/2.0
  
  REAL(DP), PARAMETER :: Ry2eV = 13.60569253_DP        ! ==2/alpha 



  WRITE (*,*), "Performing scissor correction to the band energies"
 

!  ioccgap=0
!  iunoccgap=0
!  ikoccgap=0
!  ikunoccgap=0


  ! (1) Find DFT gap

  ! Find largest occupied energy. Brute-force search algorithm.
  maxocc = energies(1,1)
  
  DO ik=1, nks
     DO ibnd=1, nocc
        IF (energies(ibnd, ik) > maxocc) THEN
           
!           ikoccgap = ik
!           ioccgap = ibnd
           maxocc = energies(ibnd, ik)
           
        ENDIF
     ENDDO
  ENDDO


  ! Find smallest unoccupied energy. Brute-force search algorithm.
  minunocc = energies(nocc+1, 1)

  DO ik=1, nks
     DO ibnd=nocc+1, nocc+nunocc
        IF (energies(ibnd, ik) < minunocc) THEN           

!           ikunoccgap = ik
!           iunoccgap = ibnd
           minunocc = energies(ibnd, ik)

        ENDIF
     ENDDO
  ENDDO

  
  ! (2) Correct band energies (in Ry)
  shift = (gap/Ry2eV - (minunocc - maxocc))/2.0_DP
  
  

  ! Shift down occupied band energies

  DO ik=1, nks
     DO ibnd=1, nocc
        energies(ibnd, ik) = energies(ibnd, ik) - shift
     ENDDO
  ENDDO
  

  ! Shift up unoccupied band energies by gap/2.0

  DO ik=1, nks
     DO ibnd=nocc+1, nocc+nunocc
        energies(ibnd, ik) = energies(ibnd, ik) + shift
     ENDDO
  ENDDO


    
   WRITE (*,*), "Band gap has been set to ", gap, "eV"

END SUBROUTINE scissor



FUNCTION eucl_norm(v, a)
  !
  ! Euclidean norm of vector v in R^3.
  ! The basis a(3,3) needs to be provided.
  ! 
  
  USE kinds,                         ONLY: DP
  
  IMPLICIT NONE
  
  REAL(DP), INTENT(IN), DIMENSION(3) :: v
  REAL(DP), INTENT(IN), DIMENSION(3,3) :: a
  
  REAL(DP) :: eucl_norm
  
  REAL(DP) :: aa11, aa22, aa33, aa12, aa13, aa23 ! Basis vector dot products

  ! Digonal elements
 aa11 = a(1,1)**2 + a(2,1)**2 + a(3,1)**2                 ! a1.a1
 aa22 = a(1,2)**2 + a(2,2)**2 + a(3,2)**2                 ! a2.a2
 aa33 = a(1,3)**2 + a(2,3)**2 + a(3,3)**2                 ! a3.a3
        
 ! Cross off diagonal elements
 aa12 = a(1,1)*a(1,2) + a(2,1)*a(2,2) + a(3,1)*a(3,2)    ! a1.a2
 aa13 = a(1,1)*a(1,3) + a(2,1)*a(2,3) + a(3,1)*a(3,3)    ! a1.a3
 aa23 = a(1,3)*a(1,2) + a(2,3)*a(2,2) + a(3,3)*a(3,2)    ! a2.a3

 eucl_norm = DSQRT(  &
      v(1)**2   * aa11 + &
      v(2)**2   * aa22 + &
      v(3)**2   * aa33 + &
      2.0_DP * v(1)*v(2) * aa12 + &
      2.0_DP * v(1)*v(3) * aa13 + &
      2.0_DP * v(3)*v(2) * aa23 )
 
END FUNCTION eucl_norm





FUNCTION eucl_norm_fast(v, aa)
  !
  ! Euclidean norm of vector v in R^3.
  ! Providing the dot products of the basis vectors aa
  ! we skip a few operations
  ! 
  
  USE kinds,                         ONLY: DP
  
  IMPLICIT NONE
  
  REAL(DP), INTENT(IN), DIMENSION(3) :: v  ! Vector in R^3
  REAL(DP), INTENT(IN), DIMENSION(6) :: aa ! Basis vector dot products
  
  REAL(DP) :: eucl_norm_fast
  

 eucl_norm_fast = DSQRT(  &
      v(1)**2   * aa(1) + &
      v(2)**2   * aa(2) + &
      v(3)**2   * aa(3) + &
      2.0_DP * v(1)*v(2) * aa(4) + &
      2.0_DP * v(1)*v(3) * aa(5) + &
      2.0_DP * v(3)*v(2) * aa(6) )
 
END FUNCTION eucl_norm_fast





SUBROUTINE create_dk_table(numk, kcoord, dk)
  !
  ! Create all possible dk= k1-k2
  ! This quantity is computed in the basis provided by input
  ! and the units of k are preserved
  !
  USE kinds,                          ONLY: DP
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: numk                   ! Number of k-points in mesh
  REAL(DP), INTENT(IN) :: kcoord(3,numk)        ! k-point coordinates
  REAL(DP), INTENT(OUT) :: dk(3,numk,numk)      ! Table storing all dk's
  
  INTEGER :: ik1, ik2                           ! k-vector indices

  
  DO ik2=1, numk
     DO ik1=1, numk
        dk(1,ik1,ik2) = kcoord(1,ik1) - kcoord(1,ik2)
        dk(2,ik1,ik2) = kcoord(2,ik1) - kcoord(2,ik2)
        dk(3,ik1,ik2) = kcoord(3,ik1) - kcoord(3,ik2)
     ENDDO
  ENDDO


END SUBROUTINE create_dk_table




SUBROUTINE find_ik(kvec, numk, kcoord, tolerance, ik)
!
! Given a k-vector, find its corresponding index ik in the k-point list
!

  USE kinds,                          ONLY: DP
  
  IMPLICIT NONE
  
  REAL(DP), INTENT(IN) :: kvec(3)
  INTEGER, INTENT(IN) :: numk
  REAL(DP), INTENT(IN) :: kcoord(3,numk)
  REAL(DP), INTENT(IN) :: tolerance
  INTEGER, INTENT(OUT) :: ik
  
  INTEGER :: ikaux
  REAL(DP) :: delta, dk(3)
  

  ik=0 ! If not found, 0 will be returned
  
  DO ikaux=1, numk
     dk(:) = kcoord(:,ikaux) - kvec(:)
     delta = SQRT( dk(1)**2 + dk(2)**2 + dk(3)**2 )
  
     IF (delta < tolerance) THEN
        ik = ikaux
        EXIT
     ENDIF
  ENDDO

  IF (ik==0) print *, "-- ERROR in find_ik subroutine: ik==0"

END SUBROUTINE find_ik




SUBROUTINE all_igk(numk, numpw, alligk)
  !
  ! Load igk in an array from file for all k-points.
  !
  ! Remember that igk(:) gives the index for G-vector list
  ! for a given k-point.
  !

  USE wvfct,                          ONLY: igk 
  USE io_files,                       ONLY: iunigk
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: numk, numpw               ! Number of k-points and number of plane waves (==num of G-vectors)
  INTEGER, INTENT(OUT) :: alligk(numpw, numk)       ! List containing all igk(:)

  INTEGER :: ik
  
  IF ( numk > 1 ) THEN
     REWIND( iunigk )
     DO ik=1, numk
        READ( iunigk ) igk
        alligk(:, ik) = igk(:)
     ENDDO
  ENDIF

END SUBROUTINE all_igk




SUBROUTINE sum_or_not(ik1, ik2, alligk, tolerance, gsumtable)
  !
  ! Create truth table to discriminate which terms contribute
  ! to the form factor sum eq. (3.16) and which do not for a given 
  ! couple of k-points indexed by ik1 and ik2. Thus this routine 
  ! needs to be called everytime k1 and k2 change.
  !
  !
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: npwx
  USE klist,                          ONLY: nks
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: tpiba

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ik1, ik2                          ! k-vector indices
  INTEGER, INTENT(IN) :: alligk(npwx, nks)
  REAL(DP), INTENT(IN) :: tolerance                        ! Tolerance for G-vector distance in RAU
  LOGICAL, INTENT(OUT) :: gsumtable(npwx, npwx, npwx)      ! Truth table for G-vector summation

  REAL(DP), DIMENSION(3) :: g1, g2, gaux                   ! G-vector coordinates
  INTEGER :: ig1, ig2, igaux                               ! G-vector indices
  REAL(DP) :: distance


  gsumtable(:,:,:) = .false.


  DO ig1=1, npwx
     IF (alligk(ig1, ik1) == 0) CYCLE
     g1(:) = g(:, alligk(ig1, ik1))

     DO ig2=1, npwx
        IF (alligk(ig2, ik2) == 0) CYCLE
        g2(:) = g(:, alligk(ig2, ik2))

        DO igaux=1, npwx               
           IF(alligk(igaux,ik2) == 0) CYCLE

           gaux(:) = g(:, alligk(igaux, ik2)) 

           ! Check if distance:=|g2-g1-gaux| < tolerance
           ! Do this operation on array gaux(:) to avoid declaring a new one

           gaux(:) = g2(:) - g1(:) - gaux(:)
           distance = tpiba * SQRT(SUM(gaux(:)**2))
           
           IF (distance < tolerance) THEN
                gsumtable(igaux,ig2,ig1) = .true.
                
             ENDIF

           ENDDO
        ENDDO
     ENDDO

     
END SUBROUTINE sum_or_not







SUBROUTINE which_sums_old(ik1, ik2, alligk, tolerance, gsumindex)
  !
  ! Index the terms in the sum that contribute
  ! to the form factor sum eq. (3.16) for a given 
  ! couple of k-points indexed by ik1 and ik2. Thus this routine 
  ! needs to be called everytime k1 and k2 change.
  !
  !
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: npwx
  USE klist,                          ONLY: nks
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: tpiba

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ik1, ik2                          ! k-vector indices
  INTEGER, INTENT(IN) :: alligk(npwx, nks)
  REAL(DP), INTENT(IN) :: tolerance                        ! Tolerance for G-vector distance in RAU
  INTEGER, INTENT(OUT) :: gsumindex(npwx, npwx)            ! Index table for G-vector summation

  REAL(DP), DIMENSION(3) :: g1, g2, gaux                   ! G-vector coordinates
  INTEGER :: ig1, ig2, igaux                               ! G-vector indices
  REAL(DP) :: distance


  ! 0 is an exit value: should be discarded
  gsumindex(:,:) = 0


  DO ig1=1, npwx
     IF (alligk(ig1, ik1) == 0) CYCLE
     g1(:) = g(:, alligk(ig1, ik1))

     DO ig2=1, npwx
        IF (alligk(ig2, ik2) == 0) CYCLE
        g2(:) = g(:, alligk(ig2, ik2))

        DO igaux=1, npwx               
           IF(alligk(igaux,ik2) == 0) CYCLE

           gaux(:) = g(:, alligk(igaux, ik2)) 
           
           ! In order to get Gaux = G2 - G1,
           ! check if distance:=|g2-g1-gaux| < tolerance
           ! Do this operation on array gaux(:) to avoid declaring a new one

           ! In the band2band case G2 = G + G'' = G1 + Gaux
           gaux(:) = g2(:) - (g1(:) + gaux(:))
           distance = tpiba * sqrt(sum(gaux(:)**2))                      
           
           IF (distance < tolerance) THEN
                gsumindex(ig2,ig1) = igaux
                
             ENDIF

           ENDDO
        ENDDO
     ENDDO

     
END SUBROUTINE which_sums_old










SUBROUTINE which_sums(ik1, ik2, alligk, tolerance, gsumindex)
  !
  ! Index the terms in the sum that contribute
  ! to the form factor sum eq. (3.16) for a given 
  ! couple of k-points indexed by ik1 and ik2. Thus this routine 
  ! needs to be called everytime k1 and k2 change.
  !
  !
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: npwx
  USE klist,                          ONLY: nks
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: tpiba

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ik1, ik2                          ! k-vector indices
  INTEGER, INTENT(IN) :: alligk(npwx, nks)
  REAL(DP), INTENT(IN) :: tolerance                        ! Tolerance for G-vector distance in RAU
  INTEGER, INTENT(OUT) :: gsumindex(npwx, npwx)            ! Index table for G-vector summation

  REAL(DP), DIMENSION(3) :: g1, g2, gaux                   ! G-vector coordinates
  INTEGER :: ig1, ig2, igaux                               ! G-vector indices
  REAL(DP) :: distance


  ! 0 is an exit value: should be discarded
  gsumindex(:,:) = 0


  DO ig1=1, npwx
     IF (alligk(ig1, ik1) == 0) CYCLE
     g1(:) = g(:, alligk(ig1, ik1))

     DO ig2=1, npwx
        IF (alligk(ig2, ik2) == 0) CYCLE
        g2(:) = g(:, alligk(ig2, ik2))

        DO igaux=1, npwx               
           IF(alligk(igaux,ik2) == 0) CYCLE

           gaux(:) = g(:, alligk(igaux, ik2)) 
           
           ! In order to get Gaux = G1 + G2,
           ! check if distance:=|g1+g2-gaux| < tolerance
           !
           ! Do this operation in place to avoid declaring a new array

           gaux(:) = g1(:) + g2(:) - gaux(:)
           distance = tpiba * sqrt(sum(gaux(:)**2))                      
           
           IF (distance < tolerance) THEN
              gsumindex(ig2,ig1) = igaux
             ENDIF

           ENDDO
        ENDDO
     ENDDO

     
END SUBROUTINE which_sums





SUBROUTINE which_sums_vacuum(ik1, ik2, alligk, tolerance, gsumindex)
  !
  ! Index the terms in the sum that contributes
  ! to the band2vacuum_formfactor sum for a given 
  ! couple of k-points indexed by ik1 and ik2. Thus this routine 
  ! needs to be called everytime k1 and k2 change.
  !
  !
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: npwx
  USE klist,                          ONLY: nks
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: tpiba

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ik1, ik2                          ! k-vector indices
  INTEGER, INTENT(IN) :: alligk(npwx, nks)
  REAL(DP), INTENT(IN) :: tolerance                        ! Tolerance for G-vector distance in RAU
  INTEGER, INTENT(OUT) :: gsumindex(npwx, npwx)            ! Index table for G-vector summation

  REAL(DP), DIMENSION(3) :: g1, g2, gaux                   ! G-vector coordinates
  INTEGER :: ig1, ig2, igaux                               ! G-vector indices
  REAL(DP) :: distance


  ! 0 is an exit value: should be discarded
  gsumindex(:,:) = 0


  DO ig1=1, npwx
     IF (alligk(ig1, ik1) == 0) CYCLE
     g1(:) = g(:, alligk(ig1, ik1))

     DO ig2=1, npwx
        IF (alligk(ig2, ik2) == 0) CYCLE
        g2(:) = g(:, alligk(ig2, ik2))

        DO igaux=1, npwx               
           IF(alligk(igaux,ik2) == 0) CYCLE

           gaux(:) = g(:, alligk(igaux, ik2)) 

           ! Check if distance:=|g2-g1-gaux| < tolerance
           ! Do this operation on array gaux(:) to avoid declaring a new one
           
           ! In the band2vacuum case we want G2 = G - G'' = G1 + Gaux
           gaux(:) = g2(:) - (g1(:) - gaux(:))
           distance = tpiba * SQRT(SUM(gaux(:)**2))                      
           
           IF (distance < tolerance) THEN
                gsumindex(ig2,ig1) = igaux
                
             ENDIF

           ENDDO
        ENDDO
     ENDDO

     
END SUBROUTINE which_sums_vacuum







SUBROUTINE count_k1plusk2(numk, k1coord, k2coord, qcoord, tolerance, count)
  !
  ! Routine for debugging purposes
  ! 
  ! Count how many k-points in the mesh can be written as q=k1+k2' 
  ! We also include the possibility of q=k1+k2
  !
  ! TODO: Include the case q=k1+k2+-bi, where bi is 
  ! a reciprocal lattice basis vector (in case the sum falls outside
  ! 1BZ).
  !
  
  USE kinds,                          ONLY: DP

  IMPLICIT NONE
  

  INTEGER, INTENT(IN) :: numk                                          ! number of k-points
  REAL(DP), DIMENSION(3,numk), INTENT(IN) :: k1coord, k2coord, qcoord  ! k-point coordinates
  REAL(DP), INTENT(IN) :: tolerance                                    ! allowed error in k-point distance to discriminate counts
  INTEGER, INTENT(OUT) :: count                                        ! number of succesive findings
  
  INTEGER :: ik1, ik2, iq
  REAL(DP) :: dq(3)
  REAL(DP) :: delta
  LOGICAL :: found

  count = 0
  DO iq=1, numk

     found = .false.
     DO ik1=1, numk
        DO ik2=1, numk

           dq(:) = qcoord(:,iq) - k1coord(:,ik1) - k2coord(:,ik2)
           delta = SQRT( dq(1)**2 + dq(2)**2 + dq(3)**2 )

           IF (delta < tolerance) found=.true.              

           IF (found) EXIT
        ENDDO
        IF (found) EXIT
     ENDDO
     IF (found) count = count + 1
  ENDDO

END SUBROUTINE count_k1plusk2




SUBROUTINE check_wf_normalization()

USE kinds,                          ONLY: DP
USE electrons_base,                 ONLY: nelt
USE wavefunctions_module,           ONLY: evc
USE wvfct,                          ONLY: igk, nbnd 
USE klist,                          ONLY: nks, ngk, wk
USE io_files,                       ONLY: nwordwfc, iunwfc, iunigk
USE buffers,                        ONLY: get_buffer


IMPLICIT NONE


INTEGER :: ik             ! k-point index
INTEGER :: ig             ! G-vector index
INTEGER :: iband          ! band index

REAL(DP) :: wfnorm        ! wavefunction norm
REAL(DP) :: acc           ! accumulate sum here

INTEGER :: numval         ! number of occupied (valence) orbitals

LOGICAL :: odd=.false.    ! Odd number of electrons?

numval = nelt/2
IF (MOD(nelt, 2) .ne. 0) THEN
   odd = .true.
   numval = numval + 1 ! If number of electrons is odd
ENDIF


! Check wavefunction normalization
     wfnorm = 0.0_DP

     IF (odd) THEN        
        DO ik = 1, nks
           call get_buffer (evc, nwordwfc, iunwfc, ik)
           DO iband=1, numval
              acc = 0.0_DP
              DO ig=1, SIZE(evc(:,1))
                 ! 2 electrons per band offer 4 possibilities for multiplying equivalent wfs (up-up, up-down, down-up and down-down) 
                 ! Multiply by 4.0 outside the G-vector loop for fully occupied bands
                 ! Multiply by 2,0 outside the G-vector loop for half occupied band
                 acc = acc + CONJG(evc(ig, iband)) * evc(ig, iband) 
              ENDDO
              ! acc, upon exiting the loop, is equal to 1.0 (wfs are normalized on each k-point)
              IF (iband .ne. numval) wfnorm = wfnorm + 4.0 * REAL(acc) * wk(ik) * wk(ik)
              IF (iband .eq. numval) wfnorm = wfnorm + 2.0 * REAL(acc) * wk(ik) * wk(ik)
           ENDDO
        ENDDO        
        
     ELSE

        DO ik = 1, nks
           call get_buffer (evc, nwordwfc, iunwfc, ik)
           DO iband=1, numval
              acc = 0.0_DP
              DO ig=1, SIZE(evc(:,1))
                 ! 2 electrons per band offer 4 possibilities for multiplying equivalent wfs (up-up, up-down, down-up and down-down)
                 ! Multiply by 4.0 outside the G-vector loop
                 acc = acc + CONJG(evc(ig, iband)) * evc(ig, iband) 
              ENDDO
              ! acc, upon exiting the loop, is equal to 1.0 (wfs are normalized on each k-point)
              wfnorm = wfnorm + 4.0 * REAL(acc) * wk(ik) * wk(ik)
           ENDDO
        ENDDO        

     ENDIF
        
END SUBROUTINE check_wf_normalization




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stuff for Dark Matter integration                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION vmin(deltaE, q, mx)
  !
  ! Evaluate vmin = deltaE / q + q/2*mx
  ! for cross section calculation.
  ! 
  !
  ! For QE compatibility, the input values should
  ! all be passed in RAU. Output is also in RAU
  ! 
  ! 
  ! -1.0 is the error value
  !
  
  USE kinds,                         ONLY: DP
  
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: deltaE               ! Energy difference between valence and conduction electrons (Ry)
  REAL(DP), INTENT(IN) :: q                    ! Norm of momentum transfer (1/bohr)
  REAL(DP), INTENT(IN) :: mx                   ! Mass of dark matter particle (RAU)
  REAL(DP) :: vmin

  REAL(DP), PARAMETER :: twobyalpha = 274.07199814110116      ! ==2/alpha, which is the conversion factor for velocities from N.U. to R.A.U.  
!  REAL(DP), PARAMETER :: speedoflight_kmps = 2.98E5           ! In km/s
  REAL(DP), PARAMETER :: zero_tol = 1.0E-8

  vmin = -1.0_DP

  IF (deltaE < 0.0) THEN
     CALL errore ('formfact', 'vmin -- deltaE must be a positive number!', 1)
  ELSEIF (q < zero_tol) THEN
     CALL errore ('formfact', 'vmin -- q negative or zero', 1)
  ELSEIF (mx < zero_tol) THEN
     CALL errore ('formfact', 'vmin -- mx is negative or zero', 1)
  ELSE
     vmin = deltaE/q + q/(2.0_DP*mx) 
  ENDIF
  
END FUNCTION vmin







FUNCTION eta(vesc, vearth, v0, vmin)
  !
  ! Evaluation of eta(v_min).
  !
  ! WARNING: All velocities are in this formula are in the 
  ! same units. eta contains the correct prefactors 
  ! in contrast with eta(v_min)

  USE kinds,                         ONLY: DP

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: vesc, vearth, v0, vmin
  REAL(DP) :: eta

  REAL(DP) :: derf ! for Portland compiler (Hopper)

!!!!  REAL(DP), PARAMETER :: twobyalpha = 274.07199814110116 ! ==2/alpha, which is the conversion factor for velocities from N.U. to R.A.U.  
!!!!  REAL(DP), PARAMETER :: K = 2.504E-09

  REAL(DP), PARAMETER :: pi = 3.141592653589793238462643383_DP
!!!!  REAL(DP), PARAMETER :: v0 = 230.0                              ! v0 in km/s

! These paramters can be used if vearth=240, v0=230 and vesc=600
! to make the evaluation of eta faster.

!  REAL(DP), PARAMETER :: A = 0.0069611336747931645
!  REAL(DP), PARAMETER :: B = 2.6681527182273292
!  REAL(DP), PARAMETER :: C = 0.000045560513361222344 
!  REAL(DP), PARAMETER :: D = 2.6681571651485263      
!  REAL(DP), PARAMETER :: E = 0.9997750917632876      


    
!!!!!!  vmin_rau = vmin/twobyalpha    ! Convert to natural units to evaluate eta


  ! Velocities below are all in km/s
  REAL(DP) :: kk         ! tien -- normalization factor
  kk = 6.75E7            ! tien -- (km/s)^3
  IF (vmin < 0.0 .or. vesc < 0.0 .or. vearth < 0.0) THEN
     CALL errore ('formfact', 'One if the input velocities in function eta() is negative', 1)

  ELSEIF (vmin <= vesc - vearth) THEN ! evaluate eta_1

     eta = -4.0_DP*DEXP(-(vesc**2/v0**2))*vearth + &
          pi**(0.5)*v0*(DERF((vmin + vearth)/v0) - DERF((vmin - vearth)/v0))
     eta = eta*v0**2*pi/(2.0_DP*vearth*kk)

  ELSEIF (vmin > vesc - vearth .and. vmin <= vesc + vearth) THEN ! evaluate eta_2

     eta = -2.0_DP*DEXP(-(vesc/v0)**2)*(vesc + vearth - vmin) + &
          pi**(0.5)*v0*(DERF(vesc/v0) - DERF((vmin - vearth)/v0))
     eta = eta*v0**2*pi/(2.0_DP*vearth*kk)
  ELSEIF (vmin > vesc + vearth) THEN
     ! vmin exceeds allowed limit-- assign zero weight
     eta = 0.0_DP
  ELSE
     print *, "vmin=", vmin
     print *, "vesc, vearth", vesc, vearth
     CALL errore ('formfact', 'vmin has a non-allowed value!', 1)
  ENDIF

END FUNCTION eta





FUNCTION bzvol(b)
  
  USE kinds,                           ONLY: DP
  USE cell_base,                       ONLY: tpiba

  IMPLICIT NONE

  Real(DP), INTENT(IN) :: B(3,3)        ! Reciprocal lattice basis vectors
  REAL(DP) :: bzvol                     ! 1BZ volume in atomic units

  ! Component expansion of the triple product formula  A.(B x C) = epsilon_{i,j,k} A^i B^j C^k
  

  bzvol =   b(1,1) * b(2,2) * b(3,3) &
          + b(2,1) * b(3,2) * b(1,3) &
          + b(3,1) * b(1,2) * b(2,3) &
          - b(2,1) * b(1,2) * b(3,3) &
          - b(1,1) * b(3,2) * b(2,3) &
          - b(3,1) * b(2,2) * b(1,3)
  


     ! convert BZ volume to atomic units
  bzvol = bzvol * tpiba**3  
    
  
END FUNCTION bzvol




SUBROUTINE k_integral(bzvol, vec, integral)
  ! 
  ! Integrate vec(ik) over 1BZ  
  !
  !
  USE kinds,                         ONLY: DP
  USE klist,                         ONLY: nks, wk

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: bzvol          ! 1BZ volume
  REAL(DP), INTENT(IN) :: vec(nks)       ! k-point dependent quantity to be integrated
  REAL(DP), INTENT(OUT) :: integral      ! integral
  
  INTEGER :: ik

  IF ( SIZE(SHAPE(vec)) .ne. 1 .and. SIZE(vec) .ne. nks) THEN
     CALL errore('formfact', 'Cannot integrate over 1BZ - array dimensions are wrong!', 1) 
  ELSE
  
     integral = 0.0_DP
     DO ik=1, nks
        integral = integral + vec(ik) * wk(ik)
     ENDDO
  
     integral = integral * bzvol
     
  ENDIF

END SUBROUTINE k_integral



SUBROUTINE bdotb(bb)
  !
  ! Calculate all 6 possible dot products between reciprocal lattice 
  ! basis vectors and store them and store them in bb(6). This is done 
  ! in units of tpiba, so bb will need to be multiplied by tpiba2 to 
  ! convert to A.U.
  ! 
  !
  USE kinds,                      ONLY: DP
  USE cell_base,                  ONLY: bg
  
  IMPLICIT NONE
  
  REAL(DP), INTENT(OUT) :: bb(6)           ! 6 terms of bi.bj
  
  
  ! Symmetric terms
  bb(1) = bg(1,1)**2 + bg(2,1)**2 + bg(3,1)**2                 ! b1.b1
  bb(2) = bg(1,2)**2 + bg(2,2)**2 + bg(3,2)**2                 ! b2.b2
  bb(3) = bg(1,3)**2 + bg(2,3)**2 + bg(3,3)**2                 ! b3.b3
  
  ! Cross terms
  bb(4) = bg(1,1)*bg(1,2) + bg(2,1)*bg(2,2) + bg(3,1)*bg(3,2)  ! b1.b2
  bb(5) = bg(1,1)*bg(1,3) + bg(2,1)*bg(2,3) + bg(3,1)*bg(3,3)  ! b1.b3
  bb(6) = bg(1,3)*bg(1,2) + bg(2,3)*bg(2,2) + bg(3,3)*bg(3,2)  ! b2.b3

END SUBROUTINE bdotb
    
    

SUBROUTINE create_bins(bintype, min, max, numbins, binsize, binedges)
  !
  !
  ! Set bin edges for binning data. In the picture 
  ! below, ei are the edge values and bi the bins
  !
  ! e1    e2    e3    e4    e5    e6
  ! |     |     |     |     |     |
  ! |_____|_____|_____|_____|_____| ...
  !   b1    b2    b3    b4    b5
  !
  USE kinds,                         ONLY: DP
 
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: bintype                                    ! Type of bin edges (linear, exponential,...)
  REAL(DP), INTENT(IN) :: min, max                                  ! min and max values  
  INTEGER, INTENT(IN) :: numbins                                    ! number of bins
  REAL(DP), INTENT(IN) :: binsize                                   ! bin size -- only for case (2)
  REAL(DP), INTENT(OUT) :: binedges(numbins+1)                      ! array storing the values of the bin edges


  INTEGER :: ierr
  REAL(DP) :: maxaux
  SELECT CASE(bintype)
  CASE (1)
     CALL linear_bins(min, max, numbins, binedges)
  CASE (2)
     maxaux = min + float(numbins) * binsize
     CALL linear_bins(min, maxaux, numbins, binedges)
  CASE (3)
     CALL exponential_bins(min, max, numbins, binedges)
  CASE default
     WRITE(*,*), "ERROR: wrong bin type selected. Creating linear bins."
     CALL linear_bins(min, max, numbins, binedges)
  END SELECT
  
END SUBROUTINE create_bins



subroutine exponential_bins(min, max, nbins, binedge)
    
  USE kinds,                         ONLY: DP
  
  !
  implicit none
  
  REAL(DP), INTENT(IN) :: min, max
  INTEGER, INTENT(IN) :: nbins
  REAL(DP), INTENT(OUT) :: binedge(nbins+1)

  REAL(DP) :: x, y, deltax
  INTEGER :: i
  
  deltax = (max-min)/float(nbins)
  
  binedge(1) = min
  DO i=1, nbins
     x = float(i)*deltax
     y = min + max*(2**(x/max) - 1)     
     binedge(i+1) = y
     
  ENDDO
  
end subroutine exponential_bins



subroutine linear_bins(min, max, nbins, binedge)
    
  USE kinds,                         ONLY: DP

  implicit none
    
  REAL(DP), INTENT(IN) :: min, max
  INTEGER, INTENT(IN) :: nbins
  REAL(DP), INTENT(OUT) :: binedge(nbins+1)

  REAL(DP) :: x, y, deltax
  INTEGER :: i
      
  deltax = (max-min)/float(nbins)
  
  DO i=0, nbins
     x = float(i)*deltax
     y = min + x
     binedge(i+1) = y
  ENDDO
  
end subroutine linear_bins




FUNCTION find_bin(numbins, binedges, value) 
  !
  ! Adrian Soto
  ! 26-06-2014
  ! Stony Brook University
  !
  ! Given a value, and the values of the bin edges,
  ! return the bin index. In the picture below,
  ! ei are the edge values and bi the bins
  !
  !e1=emin   e2    e3    e4    e5    e6  emax-1  emax  infty
  !   |      |     |     |     |     |     |     |     |
  !   |______|_____|_____|_____|_____| ... |_____|_____|
  !     b1    b2     b3    b4    b5        bmax-1 bmax
  !
  !
  ! WARNING: The edge values are assumed to be
  ! sorted in increasing order.
  !
  !
  USE kinds,                         ONLY: DP
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: numbins
  REAL(DP), INTENT(IN) :: binedges(numbins+1)  
  REAL(DP), INTENT(IN) :: value
  INTEGER :: find_bin  

  INTEGER :: i

  DO i=1, numbins

     IF (value > binedges(i) .and. value < binedges(i+1)) THEN
        find_bin = i
        EXIT
     ENDIF
     
     IF (value > binedges(numbins+1)) THEN
        find_bin = numbins+1
     ENDIF

  ENDDO

END FUNCTION



FUNCTION find_ehomo(nocc, nunocc, energies)
  !
  ! Adrian Soto
  ! 19-09-2014
  ! Stony Brook University
  !
  ! 
  ! Find the HOMO energy in Ry.
  !
  !
  ! NOTE: in the subroutine print_ks_energies()
  !       there HOMO energy is found and stored 
  !       in the variable ehomo using a different 
  !       searching method (probably better than 
  !       this one).
  !
  USE kinds,                         ONLY: DP
  USE klist,                         ONLY: nks

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nocc
  INTEGER, INTENT(IN) :: nunocc
  REAL(DP), INTENT(IN) :: energies(nocc+nunocc, nks)

  INTEGER :: ik, ibnd                                  ! Loop indices
  REAL(DP) :: find_ehomo ! Highest Valence Energy (Ry)
  
  ! Scan all occupied bands to find the largest value
  find_ehomo=energies(1,1)
  DO ik=1, nks
     DO ibnd=1, nocc 
        IF (energies(ibnd, ik) > find_ehomo) THEN
           find_ehomo=energies(ibnd,ik)
        ENDIF

     ENDDO
 ENDDO
    
END FUNCTION find_ehomo


FUNCTION find_elumo(nocc, nunocc, energies)
  !
  ! Adrian Soto
  ! 19-09-2014
  ! Stony Brook University
  !
  ! 
  ! Find the LUMO energy in Ry.
  !
  !
  ! NOTE: in the subroutine print_ks_energies()
  !       there HOMO energy is found and stored 
  !       in the variable ehomo using a different 
  !       searching method (probably better than 
  !       this one).
  !
  USE kinds,                         ONLY: DP
  USE klist,                         ONLY: nks

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nocc
  INTEGER, INTENT(IN) :: nunocc
  REAL(DP), INTENT(IN) :: energies(nocc+nunocc, nks)

  INTEGER :: ik, ibnd                                  ! Loop indices
  REAL(DP) :: find_elumo ! Lowest unoccupied energy (Ry)
  
  ! Scan all occupied bands to find the smallest value
  find_elumo=energies(nocc+nunocc,1)
  DO ik=1, nks
     DO ibnd=1, nunocc 
        IF (energies(nocc+ibnd, ik) < find_elumo) THEN
           find_elumo=energies(nocc+ibnd,ik)
        ENDIF

     ENDDO
 ENDDO
    
END FUNCTION find_elumo




SUBROUTINE print_DM_data(vesc_kmps, vearth_kmps, v0_kmps, mx_NU, ermax_NU)

  USE kinds,                         ONLY: DP
  
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: vesc_kmps, vearth_kmps, v0_kmps, mx_NU, ermax_NU

  WRITE(*,*), " "
  WRITE(*,*), "escape velocity = ", vesc_kmps, "km/s"
  WRITE(*,*), "Earth velocity = ", vearth_kmps, "km/s"
  WRITE(*,*), "DM mean velocity = ", v0_kmps, "km/s"
  WRITE(*,*), "DM mass = ", mx_NU/1.0E6_DP, "MeV/c^2" 
  WRITE(*,*), "Max recoil energy cutoff = ", ermax_NU, "eV" 
  WRITE(*,*), " "

END SUBROUTINE print_DM_data




SUBROUTINE qspace(bzvol, print2file)
  !
  ! Adrian Soto
  ! 25-10-2014
  !
  ! Calculates qspace information.
  !
  USE kinds,                         ONLY: DP
  USE klist,                         ONLY: nks, ngk, wk, xk
  USE gvect,                         ONLY: g
  USE cell_base,                     ONLY: tpiba, bg
  USE wvfct,                         ONLY: igk, ecutwfc
  USE io_files,                      ONLY: iunigk


  IMPLICIT NONE
  
  REAL(DP), INTENT(IN) :: bzvol
  LOGICAL, INTENT(IN) :: print2file
 
  REAL(DP), PARAMETER :: PI=3.1415926535897932385 
    
  INTEGER :: ik, ig
  REAL(DP) :: q(3), qnorm
  REAL(DP) :: qvol, qmax, spherevol
  REAL(DP) :: bginv(3,3)  ! inverse of bg: transforms from RLV to cartesian
  REAL(DP) :: bzvnew

  qvol = 0.0_DP
  qmax = 0.0_DP
  q=0.0_DP
  
  IF (print2file) THEN
     OPEN (UNIT=314, FILE='rlbasis.debug', STATUS='replace', ACCESS='sequential', FORM='formatted')
     
     ! Print RL vectors 
     WRITE (314, *), bg(1,:)
     WRITE (314, *), bg(2,:)   
     WRITE (314, *), bg(3,:)
     WRITE (314, *), " "
     
     bginv(:,:) = bg(:,:)
     CALL gauss(bginv, 3)
     
     ! Print inverse RL vectors, i.e. direct lattice vectors.
     WRITE (314, *), bginv(1,:)
     WRITE (314, *), bginv(2,:)   
     WRITE (314, *), bginv(3,:)
     
     CLOSE(314)
     
  ENDIF
 
  ! check bzv by calculating triple product
  bzvnew = tpiba**3 * ( &
         bg(1,1)*bg(2,2)*bg(3,3) &
       + bg(2,1)*bg(3,2)*bg(1,3) &
       + bg(3,1)*bg(1,2)*bg(2,3) &
       - bg(1,3)*bg(2,2)*bg(3,1) &
       - bg(1,2)*bg(2,1)*bg(3,3) &
       - bg(1,1)*bg(2,3)*bg(3,2) )
  
  IF (abs(bzvnew - bzvol) > 1.0E-06) &
       print *, "-- WARNING!!!!!! there is a discrepancy in the Brillouin zone volume. Check RL vectors"


  IF(print2file) THEN
     print *, " "
     print *, "-- Printing q-vectors to file qvecs.debug"
     print *, " "
     OPEN (UNIT=314, FILE='qvecs.debug', STATUS='replace', ACCESS='sequential', FORM='formatted')
  ENDIF

  
  REWIND( iunigk )
  qvol=0.0_DP
  DO ik=1, nks
     DO ig=1, ngk(ik)
        
        ! Add volume element for volume calculation
        qvol = qvol + 0.5 * wk(ik) * bzvol
        
        ! q vector in cartesian coordinates
        IF (igk(ig) == 0) CYCLE
        q(:) = xk(:,ik) + g(:,igk(ig))
        qnorm = tpiba * SQRT(SUM(q(:)**2))
        
        IF (qnorm**2 > ecutwfc) CYCLE
        
        
        IF (qnorm > qmax) qmax = qnorm
        IF (print2file) WRITE (314, *), q(1:3) 
        
     ENDDO
  ENDDO

  IF(print2file) CLOSE(314)

  spherevol=(4.0/3.0) * PI * qmax**3 

  print *, " "
  print *, "q-space volume = ", qvol, " bohr^-3"
  print *, "qmax = ", qmax, " bohr^-1"
  print *, "spherevol =", spherevol, " bohr^-3"
  print *, "spherevol/qvol =", spherevol/qvol
  print *, " "
  
END SUBROUTINE qspace




! --------------------------------------------------------------------
SUBROUTINE Gauss (a,n)       ! Invert matrix by Gauss method
! --------------------------------------------------------------------
! Serial routine from 
! http://computer-programming-forum.com/49-fortran/96e6be410fdc1511.htm
!  

  USE kinds, only:DP
  
  IMPLICIT NONE

  INTEGER :: n
  REAL(DP) :: a(n,n)
  
  ! - - - Local Variables - - -
  REAL(DP) :: b(n,n), c, d, temp(n)
  INTEGER :: i, j, k, m, imax(1), ipvt(n)
  ! - - - - - - - - - - - - - -
  
  b = a
  ipvt = (/ (i, i = 1, n) /)
  
  DO k = 1,n
     imax = MAXLOC(ABS(b(k:n,k)))
     m = k-1+imax(1)
     
     IF (m /= k) THEN
        ipvt( (/m,k/) ) = ipvt( (/k,m/) )
        b((/m,k/),:) = b((/k,m/),:)
     END IF
     d = 1/b(k,k)
     
     temp = b(:,k)
     DO j = 1, n
        c = b(k,j)*d
        b(:,j) = b(:,j)-temp*c
        b(k,j) = c
     END DO
     b(:,k) = temp*(-d)
     b(k,k) = d
  END DO
  
  a(:,ipvt) = b
  
END SUBROUTINE Gauss
