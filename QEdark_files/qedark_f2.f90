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
SUBROUTINE qedark_f2( restartmode, &
     nksf, numval, numcond, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     numqbins, dq, &
     do_scissor_correction, scissorgap)
  !
  ! Main driver to evaluate and print the spherically form factor
  ! squared as a function of |q| and E.
  ! 
  ! 
  !
  !
  USE constants, ONLY: rytoev
  USE wavefunctions_module,           ONLY: evc     ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: igk, nbnd, npwx, et,g2kin  !, btype
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
  REAL(DP), EXTERNAL :: vmin, bzvol, eucl_norm_fast
  

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

  INTEGER :: iE ! for Energy bin  
  REAL(DP) :: deltaE
  REAL(DP) :: ermax_NU !=0.0_DP                                   ! Max recoil energy cutoff in eV
  REAL(DP) :: ermax_RAU ! = ercut_NU / Ry2eV                      ! Max recoil energy cutoff in Ry
  REAL(DP), ALLOCATABLE :: binedgesE(:)
  
  CHARACTER(2) :: numbins                                         ! String containing number of bins. It can't be larger than 99
  CHARACTER(26) :: FMT                                            ! Format specifier
  
  
  REAL(DP) :: aux(max_num_mx)
  REAL(DP) :: aux_dmff(max_num_mx)
  
  COMPLEX(DP) :: f1(4)                                               ! Current |f|  used for iteration, 4 possible spin pairs
  COMPLEX(DP) :: f2                                                  ! Current |f|^2  used for iteration
  !REAL(DP), ALLOCATABLE :: ctot(:,:), cbinned(:,:,:)              ! Integrated form factor C^{i-->i'}


  REAL(DP), ALLOCATABLE :: ctot(:,:)                              ! Integrated form factor. No band indices


  CHARACTER(20) :: outfile = "C.dat"                               ! Output file name 

  INTEGER :: ik1, ik2                                             ! Indices for k-point loops 
  INTEGER :: ig1, ig2                                             ! Indices for G-vector loops
  INTEGER :: iband1, iband2                                       ! Indices for band loops 

  COMPLEX(DP), ALLOCATABLE :: evcouter(:,:)                       ! Auxiliar wavefunction array for looping. Will be used in the OUTERmost loop

  INTEGER, ALLOCATABLE :: alligk(:,:)                             ! Here we load from file the igk(:) for all k-points, i.e. alligk(ig, ik) = igk(ig)
  REAL(DP) :: dk(3,nks,nks)                                       ! Table storing all k1-k2 values  TODO: this table is antisymmetric and can be reduced
  INTEGER :: gsi(npwx, npwx)                                      ! Sum index table for band2band

  INTEGER :: nksf                                                 ! Number of k-points for formfactor calculation
  
  INTEGER :: numval, numcond                                      ! Number of occupied and unoccupied Kohn-Sham orbitals. TO BE SET BY USER!
  INTEGER :: numvaltot, numcondtot                                ! Total # bands in DFT run. numvaltot==nelec/2 and numcondtot=nbnd-numvaltot 
  INTEGER :: ivalbottom, ivaltop                                  ! Minimum and maximum values for valence band index
  INTEGER :: icondbottom, icondtop                                ! Minimum and maximum values for conduction index
  
  INTEGER :: ierr                                                 ! Error index

  REAL(DP) :: tol = 1.0E-6                                        ! Tolerance for G-vector distances


  INTEGER :: numqbins
  INTEGER :: iq                                                   ! for |q| bin
  REAL(DP) :: dq                                                  ! size of |q| bin
  REAL(DP), ALLOCATABLE :: binedgesq(:)  
  
  !DEBUGGING VARIABLES
  LOGICAL :: runff                                                ! For debugging purposes: set to false to skip the form factor evaluation                       

  runff = .true. ! Set to false to not skip ff calculation --for debugging purposes


  CALL start_clock( ' qedark_f2 ')

  print *, "           -------             "
  print *, " calculation_mode == f2"
  print *, " Calculating formfactor squared binned in E and q"
  print *, "           -------             "



!  IF (nspin .ne. 1) THEN
!     CALL errore ('qedark_f2', 'Form factor calculation works only for spin-unpolarized systems!', 1)
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
       CALL errore( 'qedark_f2 ',' Check numval and numcond values in input ', ABS(ierr) )

  IF (numval /= numvaltot) THEN
     PRINT *, " "
     PRINT *, "  WARNING: this calculation will NOT iterate over all valence bands"
     PRINT *, " ivalbottom, ivaltop = ", ivalbottom, ivaltop
     PRINT *, " "
  ENDIF


  IF (numcond /= numcondtot) THEN
     PRINT *, " "
     PRINT *, "  WARNING: this calculation will NOT iterate over all conduction bands"
     PRINT *, " icondbottom, icondtop = ", icondbottom, icondtop
     PRINT *, " "
  ENDIF

    

  !CALL volume (tpiba, bg(:,1), bg(:,2), bg(:,3), bzv) ! could also do 
  bzv = bzvol(bg) 


  print *, "  npwx      ==", npwx
  print *, "  nbnd      ==", nbnd
  print *, "  numval    ==", numval
  print *, "  numvaltot ==", numvaltot
  print *, "  numcond   ==", numcond
  print *, "  numcondtot==", numcondtot
  PRINT *, "  cell vol  == ", omega, "bohr^(3)"
  PRINT *, "  1BZ vol   == ", bzv, "bohr^(-3)"


  !CALL qspace(bzv, .false.)



  IF (nspin .ne. 1) THEN
     CALL errore ('qedark_f2', 'Form factor calculation works only for spin-unpolarized systems!', 1)
  ENDIF

  IF (noncolin .eqv. .false.) THEN
     ALLOCATE ( evcouter(npwx, nbnd) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2',' error allocating evcouter ', ABS(ierr) )

  ELSE
     ALLOCATE ( evcouter(2*npwx, nbnd) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2',' error allocating evcouter ', ABS(ierr) )
  ENDIF





  ALLOCATE ( alligk(npwx, nks) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2',' error allocating alligk ', ABS(ierr) )
  alligk(:,:)=0.0_DP ! when alligk=0, ig doesn't correspond to a G-vector in the set and should be disregarded 

  



  WRITE(*,*), " "


  IF (num_er_bins > 999) CALL errore( 'qedark_f2','Number of energy recoil bins cannot exceed 999 ', ABS(ierr) )
  IF (numqbins > 999) CALL errore( 'qedark_f2','Number of momentum transfer bins cannot exceed 999 ', ABS(ierr) )
  

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
  IF (do_scissor_correction) THEN
     CALL scissor(numvaltot, numcondtot, scissorgap, et)
  ENDIF


  ! Stuff for binning integrated form factors. 
  ! Needs to go after ermax_RAU is calculated
  IF (num_er_bins > 0) THEN
     
     ALLOCATE (binedgese(num_er_bins+1), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2',' error allocating binedgesE ', ABS(ierr) )
     
     CALL create_bins(er_bin_type, 0.0_DP, ermax_RAU, &
          num_er_bins, er_binsize, binedgese)
     WRITE (*,*), " "
     WRITE (*,*) "Creating E bins for formfactor sum ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Energy bin size=", Ry2eV*er_binsize, "eV"
  ENDIF



  ! Stuff for binning integrated form factors. 
  ! Needs to go after ermax_RAU is calculated
  IF (num_er_bins > 0) THEN
     ALLOCATE (binedgesq(num_er_bins+1), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2 ',' error allocating binedgesq ', ABS(ierr) )

          
     CALL create_bins(er_bin_type, 0.0_DP, dq*numqbins, &
          numqbins, dq, binedgesq)
     WRITE (*,*), " "
     WRITE (*,*) "Creating q bins for formfactor sum ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Momentum transfer bin size (tpiba)=", dq
  ENDIF



  ALLOCATE( ctot(numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2 ',' cannot allocate ctot ', ABS(ierr) )
  ctot(:,:) = 0.0_DP
  


  
  ! Initialize tables
  CALL bdotb(bb)
  CALL all_igk(nks, npwx, alligk)
  CALL create_dk_table(nks, xk, dk)

  

  print *, " "
  print *, " "
  print *, " "
  
  IF (runff) THEN
 
     !  WRITE(18,*) "            ik          iq         i          i'          |f_\{i-->i'}(ik,iq)|^2"     


     ik1init=1
     ik2init=1
     


     ! Loop over k-vector index of outer wavefunction evcouter
     !$omp parallel &
     !$omp private(ig1, ig2, iband1, iband2, q, qnorm, deltaE, iE, iq,  &
     !$omp ik1, ik2, evcouter, evc, gsi, wq, f1, f2) &
     !$omp reduction(+ : ctot) 
     !$omp do

     DO ik2=ik2init, nksf 

        !print *, "Iterating... @ ik2=", ik2, "from thread", omp_get_thread_num() 
        print *, "Iterating... @ ik2=", ik2, "from thread"

        ! Load wavefunctions from file         
        CALL get_buffer (evcouter, nwordwfc, iunwfc, ik2)                  

                
        ! Loop over k-vector index of inner wavefunction evc
        DO ik1=ik1init, nksf

           ! Load wavefunctions from file
           CALL get_buffer (evc, nwordwfc, iunwfc, ik1)  
           
           ! Update G-vector table for current (ik1, ik2)
           CALL which_sums(ik1, ik2, alligk, tol, gsi)
           
           ! weight of current point in q-vector space
           wq = 0.5*(wk(ik2) + wk(ik1))                     

           ! Loop over band index of outer wavefunction
           DO iband2=icondbottom, icondtop

           ! Loop over band index of inner wavefunction
              DO iband1=ivalbottom, ivaltop
                
                 deltaE = et(iband2, ik2) - et(iband1, ik1)
                 iE = find_bin(num_er_bins, binedgesE, deltaE)

                 ! Loop over G-vector index of outer wavefunction
                 DO ig2=1, ngk(ik2) 

                    ! Make sure that G-vector exists
                    IF (alligk(ig2,ik2) < 1) CYCLE                          
                    
                    ! Components of q in cartesian coordinates
                    q(:) = xk(:,ik2) - xk(:,ik1) + g(:, alligk(ig2,ik2)) 
                    
                    ! Calculate |q| and convert to RAU multiplying by tpiba
                    qnorm = tpiba * dsqrt(sum(q(:)**2))
                    
                    IF (qnorm < tol) THEN
                       ! |q|==0.0 and cannot divide by it --> leave this term out of the integral
                       CYCLE 
                    ENDIF
                            

                    iq = find_bin(numqbins, binedgesQ, qnorm)    
                    
                    ! Initialize f1 for this (ik1, ik2, ig2)
                    f1(:)=0.0
                    
                    ! Sum over G to calculate formfactor                       
                    DO ig1=1, ngk(ik1) 
                       ! Make sure that G-vector exists
                       IF (alligk(ig1,ik1) < 1) CYCLE                          
                       IF (gsi(ig2,ig1) < 1) CYCLE                          
                       

                       IF (noncolin .eqv. .false.) THEN
                          f1 = f1 + CONJG( evcouter(gsi(ig2,ig1) , iband2) ) * evc(ig1, iband1) 
                       
                       ELSE
                          f1(1) = f1(1) + CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1, iband1)            ! u-->u
                          f1(2) = f1(2) + CONJG( evcouter(gsi(ig2,ig1)+npwx, iband2) ) * evc(ig1, iband1)       ! u-->d
                          f1(3) = f1(3) + CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1+npwx, iband1)       ! d-->u
                          f1(4) = f1(4) + CONJG( evcouter(gsi(ig2,ig1)+npwx, iband2) ) * evc(ig1+npwx, iband1)  ! d-->d
                       ENDIF
                          

                    ENDDO !G-vector sum for formfactor                    
                    
                    ! Take norm squared and sum over spins (when non-collinear)
                    IF(noncolin .eqv. .false.) THEN
                       f2= CONJG(f1(1))*f1(1) 
                    ELSE
                       f1(1)= SUM(f1(:))
                       f2   = DBLE( CONJG(f1(1))*f1(1) )
                    ENDIF

                    ! Add contribution to corresponding bin
                    ctot(iq, iE) = ctot(iq, iE) + f2 * wq

                 ENDDO
                 
              ENDDO
              
           ENDDO !iband2 loop
           
        ENDDO
        
     ENDDO
     !$omp end do
     !$omp end parallel
     
     
     
     
     ! Print to file
     CALL C2file_onlyf2(outfile, numqbins, num_er_bins, ctot(:,:) )


  ENDIF  ! runff
  
  DEALLOCATE (evcouter)     
  
  
  print *, " " 
  
  
  CALL stop_clock( ' qedark_f2 ' )
  
END SUBROUTINE qedark_f2
