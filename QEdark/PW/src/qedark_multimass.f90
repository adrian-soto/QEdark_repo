!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Adrian Soto
!  06-10-2015
!  Stony Brook University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of the code QEdark v1.1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
SUBROUTINE qedark_multimass(restartmode, &
     nksf, numval, numcond, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     num_mx, mx_NU, &
     Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     do_scissor_correction, scissorgap)

  !
  ! Evaluate and print the INTEGRATED form factors (scattering rate up to
  ! constant) squared weighted by the eta function for multiple DM masses. 
  !
  !
  USE wavefunctions_module,           ONLY: evc    ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: igk, nbnd, npwx, et, g2kin  !, btype
  USE klist,                          ONLY: nks, ngk, wk, xk, nelec
  USE lsda_mod,                       ONLY: nspin
  USE io_files,                       ONLY: nwordwfc, iunwfc, iunigk
  USE buffers,                        ONLY: get_buffer
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: bg, tpiba, tpiba2, omega
  USE noncollin_module,               ONLY: noncolin

  use omp_lib

  
  IMPLICIT NONE

  INTEGER, EXTERNAL :: find_bin
  REAL(DP), EXTERNAL :: vmin, eta, bzvol, eucl_norm_fast, find_ehomo, find_elumo
  

  ! Variables specific to int_formfact(..)
  REAL(DP), PARAMETER :: Ry2eV = 13.60569253_DP                   ! ==2/alpha 
  REAL(DP), PARAMETER :: twobyalpha=274.07199814110116            ! ==2/alpha, which is the conversion factor for velocities from N.U. to Rydberg A.U. 
  REAL(DP), PARAMETER :: speedoflight=3.0E08_DP                   ! Speed of light in SI units 
  REAL(DP), PARAMETER :: RAU2kmps=(speedoflight/1000.0_DP)/twobyalpha 
  
  
  LOGICAL :: restartmode
  LOGICAL :: f1exists, f2exists, f3exists
  INTEGER :: ik1init, ik2init
  
  LOGICAL :: do_scissor_correction != .false.                      ! 
  REAL(DP) :: scissorgap                                          ! Band gap (in eV) for scissor operator corrected band energies 

  REAL(DP) :: ehomo                                               ! Highest energy of HOMO (Ry)
  REAL(DP) :: elumo                                               ! Highest energy of LUMO (Ry)
  REAL(DP) :: pwq(3)                                              ! PW momentum vector matching momentum conservation

  !REAL(DP), PARAMETER :: vearth_SI = 2.40E5_DP                    ! Earth's average velocity around the galactic center in SI units 
  REAL(DP) :: vearth_SI                                           ! Earth's average velocity around the galactic center in SI units 
  REAL(DP) :: vearth_RAU !=(twobyalpha/speedoflight)*vearth_SI      ! Earth's average velocity in RAU
  REAL(DP) :: vearth_kmps !=vearth_SI/1000.0                        ! Earth's average velocity in km/s

  REAL(DP) :: vesc_SI                                             ! Earth's escape velocity in SI units 
  REAL(DP) :: vesc_kmps != vesc_SI/1000.0                            ! Earth's escape velocity in km/s
  REAL(DP) :: vesc_RAU !=(twobyalpha/speedoflight)*vesc_SI          ! Earth's escape velocity in RAU

  REAL(DP) :: v0_SI                                                 ! DM typical velocity in SI
  REAL(DP) :: v0_kmps !=v0_SI/1000.0                                ! DM typical velocity in km/s
  REAL(DP) :: v0_RAU !=(twobyalpha/speedoflight)*v0_SI              ! DM typical velocity in RAU

  INTEGER :: idmff, imonth
  INTEGER :: nmonths = 3  
  REAL(DP) :: deltav_SI(3)
  REAL(DP) :: deltav_kmps(3) !=v0_SI/1000.0                                
  REAL(DP) :: deltav_RAU(3) !=(twobyalpha/speedoflight)*v0_SI              
  
  REAL(DP) :: me_NU = 0.510998910E6_DP                             ! Electron mass in NU 
  INTEGER, PARAMETER :: max_num_mx=99
  INTEGER :: imx, num_mx
  REAL(DP) :: mx_NU(max_num_mx) != 100.0E6_DP                     ! Dark matter particle mass in NU 
  REAL(DP) :: mx_RAU(max_num_mx) != 0.5_DP*(mx_NU/me_NU)          ! Dark matter particle mass in RAU (me=0.5)

  REAL(DP) :: vmin_aux(max_num_mx)                                ! auxiliary vmin for looping

  REAL(DP) :: bb(6)                                               ! Dot products of reciprocal lattice basis vectors
  REAL(DP) :: bzv                                                 ! 1BZ volume


  REAL(DP) :: q(3)                                                ! Coordinates of q vector
  REAL(DP) :: qnorm                                               ! Norm of q vector
  REAL(DP) :: wq                                                  ! weight of q-point for q-space integration
  
  
  INTEGER :: er_bin_type                                          ! Type of bins for integrated form factors
  INTEGER :: num_er_bins                                          ! Number of bins for integrated form factors
  REAL(DP) :: er_binsize                                          ! Bin size (Ry) for linear bins with user-provided size 

  INTEGER :: ibin  
  REAL(DP) :: deltaE
  REAL(DP) :: ermax_NU !=0.0_DP                                   ! Max recoil energy cutoff in eV
  REAL(DP) :: ermax_RAU ! = ercut_NU / Ry2eV                      ! Max recoil energy cutoff in Ry
  REAL(DP), ALLOCATABLE :: binedges(:)
  
  CHARACTER(2) :: numbins                                         ! String containing number of bins. It can't be larger than 99
  CHARACTER(26) :: FMT                                            ! Format specifier
  
  
  REAL(DP) :: aux(max_num_mx)
  REAL(DP) :: aux_dmff(max_num_mx)
  
  COMPLEX(DP) :: f1(4)                                               ! Current |f|  used for iteration, 4 for all spin pairs (uu,ud,du,dd)
  COMPLEX(DP) :: f2                                                  ! Current |f|^2  used for iteration


  REAL(DP), ALLOCATABLE :: ctot(:,:,:), cbinned(:,:,:,:)             ! Binned integral. This is the final program output.


!  CHARACTER(20) :: formfactorfile = "formfactor.dat"              ! Output file name
 CHARACTER(20) :: outfile = "C.dat"                               ! Output file name 

  INTEGER :: ik1, ik2                                             ! Indices for k-point loops 
  INTEGER :: ig1, ig2                                             ! Indices for G-vector loops
  INTEGER :: iband1, iband2                                       ! Indices for band loops 

  INTEGER :: numpwk1, numpwk2                                     ! Number of plane waves for a given k-vector  
  COMPLEX(DP), ALLOCATABLE :: evcouter(:,:)                       ! Auxiliar wavefunction array for looping. Will be used in the OUTERmost loop

  INTEGER, ALLOCATABLE :: alligk(:,:)                             ! Here we load from file the igk(:) for all k-points, i.e. alligk(ig, ik) = igk(ig)
  REAL(DP) :: dk(3,nks,nks)                                       ! Table storing all k1-k2 values  TODO: this table is antisymmetric and can be reduced

  INTEGER :: gsi(npwx, npwx)                                      ! G-vector sum index table for one pair of k-points

  INTEGER :: nksf                                                 ! Number of k-points for formfactor calculation  

  INTEGER :: numval, numcond                                      ! Number of occupied and unoccupied Kohn-Sham orbitals. TO BE SET BY USER!
  INTEGER :: numvaltot, numcondtot                                ! Total # bands in DFT run. numvaltot==nelec/2 and numcondtot=nbnd-numvaltot
  INTEGER :: ivalbottom, ivaltop                                  ! Minimum and maximum values for valence band index
  INTEGER :: icondbottom, icondtop                                  ! Minimum and maximum values for conduction band index  

  INTEGER :: ierr                                                 ! Error index

  REAL(DP) :: tol = 1.0E-6                                        ! Tolerance for G-vector distances

  
  !DEBUGGING VARIABLES
  LOGICAL :: runff                                                ! For debugging purposes: set to false to skip the form factor evaluation




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE INSTRUCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  runff = .true. ! Set to false to not skip ff calculation --for debugging purposes


  CALL start_clock( ' qedark_multimass ')


  print *, "           -------             "
  print *, " calculation_mode set to multimass"
  print *, " Calculating momentum-integrated, energy-dependent "
  print *, " DM-electron scattering rates for multiple DM masses."
  print *, "           -------             "




  
!  IF (nspin .ne. 1) THEN
!     CALL errore ('qedark_multimass', 'Form factor calculation works only for spin-unpolarized systems!', 1)
!  ENDIF


  IF ( nksf > nks .or. nksf < 0 ) &
       CALL errore( 'qedark_f2 ',' nksf has a non-allowed value. Check input. ', ABS(ierr) )


  ! Band indices
  IF ( noncolin .eqv. .false.) THEN
     numvaltot = nelec/2
  ELSE
     numvaltot = nelec
  ENDIF
  numcondtot= nbnd-numvaltot
  ivalbottom = numvaltot-numval+1 
  ivaltop = numvaltot
  icondbottom = numvaltot+1
  icondtop = numvaltot + numcond


  IF( numval>numvaltot .or. numcond>numcondtot .or. numval<1 .or. numcond<1 ) &
       CALL errore( 'formfactor ',' Check numval and numcond values in input ', ABS(ierr) )  
  
  IF (numval /= numvaltot) THEN
     PRINT *, " " 
     PRINT *, "  WARNING: this calculation will NOT iterate over all valence bands" 
     PRINT *, " ivalbottom, ivaltop = ", ivalbottom, ivaltop
     PRINT *, " "
  ENDIF


  IF (numcond /= numcondtot) THEN
     PRINT *, " " 
     PRINT *, "  WARNING: this calculation will NOT iterate over all conduction bands" 
     PRINT *, " icondbottom, icondtop = ", ivalbottom, ivaltop
     PRINT *, " icondbottom, icondtop = ", icondbottom, icondtop
     PRINT *, " "
  ENDIF


  !CALL volume (tpiba, bg(:,1), bg(:,2), bg(:,3), bzv) ! could also do 
  bzv = bzvol(bg) 


  print *, "  npwx      ==", npwx
  print *, "  nbnd      ==", nbnd
  print *, "  nksf      ==", nksf
  print *, "  numval    ==", numval
  print *, "  numvaltot ==", numvaltot
  print *, "  numcond   ==", numcond
  print *, "  numcondtot==", numcondtot
  PRINT *, "  cell vol  == ", omega, "bohr^(3)"
  PRINT *, "  1BZ vol   == ", bzv, "bohr^(-3)"
  


  !CALL qspace(bzv, .false.)


  
  IF (noncolin .eqv. .false.) THEN
     ALLOCATE ( evcouter(npwx, nbnd) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_multimass ',' error allocating evcouter ', ABS(ierr) )

  ELSE ! spin up and spin down wf coefficients both in the same array of size 2*npwx
     ALLOCATE ( evcouter(2*npwx, nbnd) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_multimass ',' error allocating evcouter ', ABS(ierr) )


  ENDIF





  ALLOCATE ( alligk(npwx, nks) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_multimass ',' error allocating alligk ', ABS(ierr) )
  alligk(:,:)=0.0_DP ! when alligk=0, ig doesn't correspond to a G-vector in the set and should be disregarded 
  


  IF (num_er_bins > 999) CALL errore( 'qedark_multimass ','Number of energy recoil bins cannot exceed 999 ', ABS(ierr) )


  vearth_RAU = (twobyalpha/speedoflight)*vearth_SI
  vearth_kmps = vearth_SI/1000.0
  vesc_kmps = vesc_SI/1000.0
  vesc_RAU = (twobyalpha/speedoflight)*vesc_SI
  v0_kmps = v0_SI/1000.00                      
  v0_RAU = (twobyalpha/speedoflight)*v0_SI    
  deltav_kmps(:) = deltav_SI(:)/1000.0
  deltav_RAU(:) = (twobyalpha/speedoflight)*deltav_SI(:)
  mx_RAU(:) = 0.5_DP*(mx_NU(:)/me_NU)
  ermax_RAU = ermax_NU / Ry2eV 


!!!!!!!!!  CALL print_DM_data(vesc_kmps, vearth_kmps, v0_kmps, mx_NU, ermax_NU)


  print *, " "
  IF (do_scissor_correction) &
     CALL scissor(numvaltot, numcondtot, scissorgap, et)



  ! Stuff for binning integral. 
  ! Needs to go after ermax_RAU is calculated
  IF (num_er_bins > 0) THEN
     
     ALLOCATE (binedges(num_er_bins+1), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_multimass ',' error allocating alligk ', ABS(ierr) )
     
     CALL create_bins(er_bin_type, 0.0_DP, ermax_RAU, &
          num_er_bins, er_binsize, binedges)
     WRITE (*,*), " "
     WRITE (*,*) "Creating bins for C integral ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Energy bin size=", Ry2eV*er_binsize, "eV"
     
     

     ! This array contains the integrated formfactor. Second index labels month. Third array labels functional dependence of F_DM(q)     
     ALLOCATE( ctot(3, nmonths, num_mx) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_multimass',' cannot allocate ctot ', ABS(ierr) )
     ctot(:,:,:) = 0.0_DP
     

     ALLOCATE( cbinned(3, nmonths, num_er_bins+1, num_mx) , STAT=ierr )  ! +1 because last bin contains contributions with deltaE up to infinity
     IF( ierr /= 0 ) &
          CALL errore( 'formfactor ',' cannot allocate cbinned ', ABS(ierr) )
     cbinned(:,:,:,:) = 0.0_DP
     
  ENDIF


  
  ! Initialize tables
  CALL bdotb(bb)
  CALL all_igk(nks, npwx, alligk)
  CALL create_dk_table(nks, xk, dk)

  
  print *, " "
  print *, " "
  
  IF (runff) THEN
 
     !  WRITE(18,*) "            ik          iq         i          i'          |f_\{i-->i'}(ik,iq)|^2"     




     ik1init=1
     ik2init=1
     


     ! Loop over k-vector index of outer wavefunction evcouter
     !$omp parallel &
     !$omp private(ig1, ig2, iband1, iband2, q, qnorm, deltaE, ibin,  &
     !$omp ik1, ik2, evcouter, evc, gsi, wq, &
     !$omp f1, f2, vmin_aux, aux, aux_dmff, idmff, imonth) &
     !$omp reduction(+ : ctot, cbinned) 
     !$omp do

     DO ik2=ik2init, nksf 
        !numpwk2 = ngk(ik2)     ! This can be used instead of npwx in the G-vector loops to avoid extra checks
        

        !print *, "Iterating... @ ik2=", ik2, "from thread", omp_get_thread_num() 
        print *, "Iterating... @ ik2=", ik2
        

        ! Load wavefunctions from file         
        CALL get_buffer (evcouter, nwordwfc, iunwfc, ik2)                  

        ! Loop over k-vector index of inner wavefunction evc
        DO ik1=ik1init, nksf
           !numpwk1 = ngk(ik1)    ! This can be used instead of npwx in the G-vector loops to avoid extra checks
           

           !!!!!!!! CALL SaveRestartData(ik1, ik2, num_er_bins, nmonths, ctot, cbinned)
                      
           ! Update G-vector table for current (ik1, ik2)

           CALL which_sums(ik1, ik2, alligk, tol, gsi)           

                      
           ! Load wavefunctions from file
           CALL get_buffer (evc, nwordwfc, iunwfc, ik1) 
           
           ! Loop over band index of outer wavefunction (conduction)
           DO iband2=icondbottom, icondtop

           ! Loop over band index of inner wavefunction (valence)
              DO iband1=ivalbottom, ivaltop

                 ! Loop over G-vector index of outer wavefunction

                 DO ig2=1, ngk(ik2) !numpwk2

                    ! Make sure that G-vector exists
                    IF (alligk(ig2,ik2) < 1) CYCLE                          
                    
                    ! Components of q in cartesian coordinates
                    q(:) = xk(:,ik2) - xk(:,ik1) + g(:, alligk(ig2,ik2)) 
                    
                    ! Calculate |q| and convert to RAU multiplying by tpiba
                    qnorm = tpiba * sqrt(sum(q(:)**2))
                                        
                    IF (qnorm < tol) THEN
                       ! |q|==0.0 and cannot divide by it --> leave this term out of the integral
                       CYCLE 
                    ENDIF

                    
                    ! Initialize f1 for this (ik1, ik2, ig2)
                    f1(:)=0.0

                    ! Energy difference of the interband transition
                    deltaE = et(iband2, ik2) - et(iband1, ik1)
                    
                    
                    ! Find corresponding energy bin
                    IF (num_er_bins > 0) &
                         ibin = find_bin(num_er_bins, binedges, deltaE)    
                    
                    

                    ! Sum over G to calculate formfactor                       
                    
                    DO ig1=1, ngk(ik1) 
                       ! Make sure that G-vector exists
                       IF (alligk(ig1,ik1) < 1) CYCLE                          
                       IF (gsi(ig2,ig1) < 1) CYCLE                          
                       
                       ! Calculate form factor term.
                       ! Not storing the full form factor in the memory since 
                       ! we're dumping each form factor term into the integral.
                       
                       IF (noncolin .eqv. .false.) THEN
                          f1(1) = f1(1) +  CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1, iband1) 
                       
                       ELSE

                          f1(1) = f1(1) + CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1, iband1)            ! u-->u
                          f1(2) = f1(2) + CONJG( evcouter(gsi(ig2,ig1)+npwx, iband2) ) * evc(ig1, iband1)       ! u-->d
                          f1(3) = f1(3) + CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1+npwx, iband1)       ! d-->u
                          f1(4) = f1(4) + CONJG( evcouter(gsi(ig2,ig1)+npwx, iband2) ) * evc(ig1+npwx, iband1)  ! d-->d

                       ENDIF
                       
                    ENDDO !G-vector sum in (12)
                    
                    
                    ! Take norm squared and sum over spins if needed
                    IF (noncolin .eqv. .false.) THEN
                       f2=DBLE( CONJG(f1(1))*f1(1) )
                    ELSE
                       f1(1)=SUM(f1(:))
                       f2=DBLE( CONJG(f1(1))*f1(1) )
                    ENDIF

                    ! Calculate vmin and convert to km/s
                    DO imx=1, num_mx
                       vmin_aux(imx) = vmin(deltaE, qnorm, mx_RAU(imx) ) * RAU2kmps
                    ENDDO
                    
                    
                    
                    ! In general the q-vector does not sit on the k-mesh. 
                    ! Set weight to average of the weights of the two 
                    ! k-points involved in the evaluation of q
                    wq = 0.5*(wk(ik2) + wk(ik1))
                    
                    ! Loops over F_DM and month
                    DO imonth=1, nmonths
                       
                       ! Current contribution to integral
                       ! Weight by (eta/|q|) and integrate over k and q
                       DO imx=1, num_mx
                          aux(imx) = f2 * &
                               eta(vesc_kmps, vearth_kmps+deltav_kmps(imonth), v0_kmps, vmin_aux(imx)) * &
                               (0.25 * wk(ik1)* wq * bzv**2) / qnorm
                       ENDDO
                       
                       DO idmff=1, 3
                          ! Dark matter formfactor F_DM(q) ~ 1, 1/q or 1/q**2     
                          ! Multiply by F_DM^2
                          DO imx=1, num_mx
                             aux_dmff(imx) = aux(imx) / (qnorm **( 2*(idmff-1)) )
                          ENDDO
                          
                          ! Add contribution to corresponding bin
                          IF (num_er_bins > 0) THEN
                             ! Add contribution to corresponding bin
                             DO imx=1, num_mx
                                cbinned(idmff, imonth, ibin, imx) = cbinned(idmff, imonth, ibin, imx) + aux_dmff(imx)
                             ENDDO
                          ENDIF
                          
                          ! Add contribution to total
                          DO imx=1, num_mx
                             ctot(idmff, imonth, imx) = ctot(idmff, imonth, imx) + aux_dmff(imx)
                          ENDDO
                       
                       ENDDO

                    ENDDO

                 ENDDO

              ENDDO

           ENDDO !iband2 loop

        ENDDO
        
        !IF(restartmode) ik1init=1 ! Because ik1init needs to be set back to 1 for the next ik2 iteration!!!
     ENDDO
     !$omp end do
     !$omp end parallel
        


     print *, " " 
     print *, " Calculation finished "
     print *, " "                 


     ! Print to files
     DO imx=1, num_mx
        write (outfile, "(A2,I2.2,A4)") "C.", imx, ".dat"        
        WRITE (*,*), "Writing output to ", outfile
        CALL C2file(outfile, num_er_bins, nmonths, ctot(:,:,imx), cbinned(:,:,:,imx) )
     ENDDO
     
   
  ENDIF  ! runff
  

  DEALLOCATE (evcouter)     
  

  CALL stop_clock( ' qedark_multimass ' )
  
END SUBROUTINE qedark_multimass
