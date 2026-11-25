!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dfunct
!
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE newd( ) 
  !--------------------------------------------------------------------------
  !! Wrapper routine for Norm-Conserving and Ultrasoft/PAW cases
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : ityp, ntyp => nsp
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : nh
  IMPLICIT NONE
  !
  !! do not do anything if there are no projectors at all, in order to
  !! prevent trouble with zero-size allocation / copy to gpu
  !
  IF ( ALL( nh(1:ntyp) == 0 ) ) RETURN
  !
  !$acc enter data create(ityp)
  !$acc update device(ityp)
  IF ( .NOT. okvan ) THEN
     !
     CALL newd_nous ( )
     !
  ELSE
     !
     CALL newd_us ( )
     !
  END IF
  !$acc exit data delete(ityp)
  !
END SUBROUTINE newd
!
!----------------------------------------------------------------------------
SUBROUTINE newd_nous( ) 
  !----------------------------------------------------------------------------
  !
  !! no ultrasoft potentials: use bare coefficients for projectors
  !
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : deeq, deeq_nc, dvan, dvan_so
  USE uspp_param,           ONLY : nh
  USE noncollin_module,     ONLY : noncolin, domag, nspin_mag, lspinorb
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht
  ! counters and auxiliary variables
  !
  !$acc data present (deeq, deeq_nc, dvan, dvan_so)
  DO nt = 1, ntyp
     !
     nht = nh(nt)
     if ( nht <= 0 ) CYCLE
     !
     IF ( lspinorb ) THEN
        !
        !$acc parallel loop collapse(4)
        DO is =  1, nspin
           DO na = 1, nat
              DO jh = 1, nht
                 DO ih = 1, nht
                    IF ( ityp(na) == nt ) deeq_nc(ih,jh,na,is) = dvan_so(ih,jh,is,nt)
                 END DO
              END DO
           END DO
        END DO
        !
     ELSE IF ( noncolin ) THEN
        !
        !$acc parallel loop collapse(3)
        DO na = 1, nat
           DO jh = 1, nht
              DO ih = 1, nht
                 IF ( ityp(na) == nt ) THEN
                    deeq_nc(ih,jh,na,1) = dvan(ih,jh,nt)
                    deeq_nc(ih,jh,na,2) = ( 0.D0, 0.D0 )
                    deeq_nc(ih,jh,na,3) = ( 0.D0, 0.D0 )
                    deeq_nc(ih,jh,na,4) = dvan(ih,jh,nt)
                 END IF
              END DO
           END DO
        END DO
        !
     ELSE
        !
        !$acc parallel loop collapse(4)
        DO is = 1, nspin
           DO na = 1, nat
              DO jh = 1, nht
                 DO ih = 1, nht
                    !
                    IF ( ityp(na) == nt ) deeq(ih,jh,na,is) = dvan(ih,jh,nt)
                    !
                 END DO
              END DO
           END DO
        END DO
        !
     END IF
     !
     ! ... sync with CPU (not sure this is needed)
     if (noncolin) then
        !$acc update self(deeq_nc)
     else
        !$acc update self(deeq)
     endif
     !
  END DO
  !$acc end data
  !
END SUBROUTINE newd_nous
!
!-------------------------------------------------------------------------
SUBROUTINE newq_acc(vr,deeq,skip_vltot)
  !----------------------------------------------------------------------
  !! This routine computes the integral of the perturbed potential with
  !! the Q function
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega, tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : g, gg, ngm, gstart, mill, eigts1, eigts2, eigts3
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vltot
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : nspin_mag
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! Input: potential , output: contribution to integral
  REAL(kind=dp), intent(in)  :: vr(dfftp%nnr,nspin)
  REAL(kind=dp), intent(out) :: deeq( nhm, nhm, nat, nspin )
  LOGICAL, intent(in) :: skip_vltot !If .false. vltot is added to vr when necessary
  ! INTERNAL
  INTEGER :: ngm_s, ngm_e, ngm_l
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ig, nt, ih, jh, na, is, ijh, nij, nb, nab, nhnt
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:), aux(:,:), qgm(:,:)
    ! work space
  REAL(DP), ALLOCATABLE :: ylmk0(:,:), qmod(:), deeaux(:,:)
    ! spherical harmonics, modulus of G
  REAL(DP) :: fact
  !
  ! variable to map index of atoms of the same type
  INTEGER, ALLOCATABLE :: na_to_nab(:)
  !
  IF ( gamma_only ) THEN
     fact = 2.0_dp
  ELSE
     fact = 1.0_dp
  ENDIF
  !
  !$acc kernels
  deeq(:,:,:,:) = 0.D0
  !$acc end kernels
  !
  ! With k-point parallelization, distribute G-vectors across processors
  ! ngm_s = index of first G-vector for this processor
  ! ngm_e = index of last  G-vector for this processor
  ! ngm_l = local number of G-vectors 
  !
  CALL divide (inter_pool_comm, ngm, ngm_s, ngm_e)
  ngm_l = ngm_e-ngm_s+1
  IF ( ngm_l < 1 ) CALL errore('newq_acc','no G-vectors?!?',1)
  !
  ALLOCATE(na_to_nab(nat))
#if defined(_OPENACC)
  ALLOCATE( vaux(ngm_l,nspin_mag), qmod(ngm_l), ylmk0( ngm_l, lmaxq*lmaxq ) )
#else
  ! ... CPU: qmod and ylmk0 are computed per G block in newq_blocked
  ALLOCATE( vaux(ngm_l,nspin_mag), qmod(1), ylmk0(1,1) )
#endif
  !$acc data create( na_to_nab, vaux, qmod, ylmk0 )
  !
#if defined(_OPENACC)
  CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
  !
  !$acc parallel loop
  DO ig = 1, ngm_l
     qmod (ig) = SQRT(gg(ngm_s+ig-1))*tpiba
  ENDDO
#endif
  ! ... Fourier transform of the total effective potential
  !
  DO is = 1, nspin_mag
     IF ( (nspin_mag==4 .AND. is/=1) .OR. skip_vltot ) THEN
        CALL rho_r2g( dfftp, vr(:,is:is), vaux(:,is:is), igs=ngm_s )
     ELSE
        CALL rho_r2g( dfftp, vr(:,is:is), vaux(:,is:is), v=vltot, igs=ngm_s )
     END IF
  END DO
  !
#if !defined(_OPENACC)
  ! ... CPU: G-blocked algorithm, no ngm-sized Q(G) or V(G)S(G) arrays
  CALL newq_blocked( ngm_s, ngm_l, vaux, fact, deeq )
#else
  DO nt = 1, ntyp
     !
     IF ( upf(nt)%tvanp ) THEN
        !
        ! count max number of atoms of type nt, create mapping table
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
           IF ( ityp(na) == nt ) THEN
              na_to_nab(na) = nab
           ELSE
              na_to_nab(na) = -1
           END IF
        END DO
        IF ( nab == 0 ) CYCLE ! No atoms for this type (?!?)
        !$acc update device(na_to_nab)
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nhnt = nh(nt)
        nij = nh(nt)*(nh(nt)+1)/2
        ALLOCATE ( qgm(ngm_l,nij) )
        ALLOCATE ( aux (ngm_l, nab ), deeaux(nij, nab) )
        !
        !$acc data create( qgm, aux, deeaux ) present(eigts1, eigts2, eigts3, mill)
        !
        ! ... Compute and store Q(G) for this atomic species 
        ! ... (without structure factor)
        !
        ijh = 0
        DO ih = 1, nhnt
           DO jh = ih, nhnt
              ijh = ijh + 1
              CALL qvan2 ( ngm_l, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
           END DO
        END DO
        !
        ! ... Compute and store V(G) times the structure factor e^(-iG*tau)
        !
        DO is = 1, nspin_mag
           !$acc parallel loop collapse(2)
           DO na = 1, nat
              DO ig = 1, ngm_l
                 nb = na_to_nab(na)
                 IF (nb > 0) &
                    aux(ig,nb) = vaux(ig,is) * CONJG ( &
                      eigts1(mill(1,ngm_s+ig-1),na) * &
                      eigts2(mill(2,ngm_s+ig-1),na) * &
                      eigts3(mill(3,ngm_s+ig-1),na) )
              END DO
           END DO
           !
           ! ... here we compute the integral Q*V for all atoms of this kind
           !
           !$acc host_data use_device(qgm,aux,deeaux)
           CALL MYDGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm, 2*ngm_l, aux, &
                    2*ngm_l, 0.0_dp, deeaux, nij )
           IF ( gamma_only .AND. gstart == 2 ) &
                CALL MYDGER(nij, nab,-1.0_dp, qgm, 2*ngm_l,aux,2*ngm_l,deeaux,nij)
           !$acc end host_data
           !
           nhnt = nh(nt)
           !$acc parallel loop collapse(3)
           DO na = 1, nat
              DO ih = 1, nhnt
                 DO jh = 1, nhnt
                    nb = na_to_nab(na)
                    IF (nb > 0) THEN
                       ijh = jh + ((ih-1)*(2*nhnt-ih))/2
                       IF (jh >= ih) deeq(ih,jh,na,is) = omega * deeaux(ijh,nb)
                       IF (jh > ih) deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
                    END IF
                 END DO
              END DO
           END DO
           !
        END DO
        !$acc end data
        !
        DEALLOCATE ( deeaux, aux )
        DEALLOCATE ( qgm )
        !
     END IF
     !
  END DO
#endif
  !
  !$acc host_data use_device(deeq)
  CALL mp_sum( deeq( :, :, :, 1:nspin_mag ), inter_pool_comm )
  CALL mp_sum( deeq( :, :, :, 1:nspin_mag ), intra_bgrp_comm )
  !$acc end host_data
  !$acc end data
  DEALLOCATE( qmod, ylmk0, vaux )
  DEALLOCATE(na_to_nab)
  !
END SUBROUTINE newq_acc
!
#if !defined(_OPENACC)
!----------------------------------------------------------------------
SUBROUTINE newq_blocked( ngm_s, ngm_l, vaux, fact, deeq )
  !----------------------------------------------------------------------
  !! CPU version of the integral \(\int V_{eff}(r) Q_{ij}(r-\tau_a) dr\),
  !! computed in G space on blocks of \(\text{gblk}\) G-vectors owned by one
  !! OpenMP thread each. Per block and species: spherical harmonics, |G|,
  !! \(Q_{ij}(G)\) (qvan2, once for all spin components), \(V(G)S^*_a(G)\)
  !! and a DGEMM accumulating \(\sum_G Q^*_{ij}(G) V(G) S^*_a(G)\) into a
  !! per-thread partial sum. Partial sums are added in thread order, so the
  !! result is deterministic for a given number of threads.
  !
  !! Compared with the original algorithm this removes the work arrays
  !! qgm(ngm_l,nij), aux(ngm_l,nat) and ylmk0(ngm_l,lmaxq^2); extra memory is
  !! O(gblk*(nij+nat+lmaxq^2)) per thread plus nij*nat*nspin per thread for
  !! the partial sums. Results differ from the original only by the order of
  !! the floating-point summation over G (rounding level).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega, tpiba
  USE gvect,                ONLY : g, gg, gstart, mill, eigts1, eigts2, eigts3
  USE lsda_mod,             ONLY : nspin
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : nspin_mag
!$ USE omp_lib,             ONLY : omp_get_max_threads, omp_get_thread_num
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: ngm_s, ngm_l
  !! global index of the first local G-vector, number of local G-vectors
  COMPLEX(DP), INTENT(IN)    :: vaux(ngm_l,nspin_mag)
  !! V_eff(G) on the local G-vectors
  REAL(DP),    INTENT(IN)    :: fact
  !! 2 with gamma tricks, 1 otherwise
  REAL(DP),    INTENT(INOUT) :: deeq(nhm,nhm,nat,nspin)
  !! output (must be zeroed by the caller)
  !
  INTEGER, PARAMETER :: gblk = 512
  !
  INTEGER :: nt, na, nb, nab, nij, ig0, n, ig, igg, ih, jh, ijh, is, nth, ith, it
  INTEGER,     ALLOCATABLE :: atoms(:)
  REAL(DP),    ALLOCATABLE :: deeaux(:,:,:,:), qmod(:), ylm(:)
  COMPLEX(DP), ALLOCATABLE :: qgm(:,:), aux(:,:)
  !
  nth = 1
!$ nth = omp_get_max_threads()
  !
  DO nt = 1, ntyp
     !
     IF ( .NOT. upf(nt)%tvanp ) CYCLE
     nab = COUNT( ityp(1:nat) == nt )
     IF ( nab == 0 ) CYCLE
     nij = nh(nt)*(nh(nt)+1)/2
     !
     ALLOCATE( atoms(nab) )
     nb = 0
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           nb = nb + 1
           atoms(nb) = na
        ENDIF
     ENDDO
     ALLOCATE( deeaux(nij,nab,nspin_mag,0:nth-1) )
     deeaux(:,:,:,:) = 0.0_DP
     !
     !$omp parallel default(shared) &
     !$omp   private(ith,ig0,n,igg,ig,nb,na,is,ih,jh,ijh,qgm,aux,qmod,ylm)
     ith = 0
!$   ith = omp_get_thread_num()
     ALLOCATE( qgm(gblk,nij), aux(gblk,nab), qmod(gblk) )
     ! ... ylm is 1D: passed (sequence association) as ylm(n,lmaxq**2)
     ALLOCATE( ylm(gblk*lmaxq*lmaxq) )
     !
     !$omp do schedule(static)
     DO ig0 = 1, ngm_l, gblk
        !
        n   = MIN( gblk, ngm_l - ig0 + 1 )
        igg = ngm_s + ig0 - 1           ! global index of the first G of the block
        !
        ! ... Q_ij(G) for this block (no structure factor)
        CALL ylmr2( lmaxq*lmaxq, n, g(1,igg), gg(igg), ylm )
        DO ig = 1, n
           qmod(ig) = SQRT( gg(igg+ig-1) ) * tpiba
        ENDDO
        ijh = 0
        DO ih = 1, nh(nt)
           DO jh = ih, nh(nt)
              ijh = ijh + 1
              CALL qvan2( n, ih, jh, nt, qmod, qgm(1,ijh), ylm )
           ENDDO
        ENDDO
        !
        DO is = 1, nspin_mag
           ! ... V(G) times the structure factor e^(-iG*tau)
           DO nb = 1, nab
              na = atoms(nb)
              DO ig = 1, n
                 aux(ig,nb) = vaux(ig0+ig-1,is) * CONJG( &
                              eigts1(mill(1,igg+ig-1),na) * &
                              eigts2(mill(2,igg+ig-1),na) * &
                              eigts3(mill(3,igg+ig-1),na) )
              ENDDO
           ENDDO
           ! ... partial integral Q*V over this block, for all atoms of this kind
           CALL DGEMM( 'C', 'N', nij, nab, 2*n, fact, qgm, 2*gblk, aux, 2*gblk, &
                       1.0_DP, deeaux(1,1,is,ith), nij )
           ! ... G=0 is counted once with gamma tricks (as in newq_acc)
           IF ( gamma_only .AND. gstart == 2 .AND. ig0 == 1 ) &
              CALL DGER( nij, nab, -1.0_DP, qgm, 2*gblk, aux, 2*gblk, &
                         deeaux(1,1,is,ith), nij )
        ENDDO
        !
     ENDDO
     !$omp end do
     !
     DEALLOCATE( qgm, aux, qmod, ylm )
     !$omp end parallel
     !
     ! ... deterministic reduction of the per-thread partial sums
     DO it = 1, nth-1
        deeaux(:,:,:,0) = deeaux(:,:,:,0) + deeaux(:,:,:,it)
     ENDDO
     !
     DO is = 1, nspin_mag
        DO nb = 1, nab
           na = atoms(nb)
           ijh = 0
           DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                 ijh = ijh + 1
                 deeq(ih,jh,na,is) = omega * deeaux(ijh,nb,is,0)
                 deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
     DEALLOCATE( deeaux, atoms )
     !
  ENDDO
  !
END SUBROUTINE newq_blocked
#endif
  !
!----------------------------------------------------------------------------
SUBROUTINE newd_us( ) 
  !----------------------------------------------------------------------------
  !! This routine computes the integral of the effective potential with
  !! the Q function and adds it to the bare ionic D term which is used
  !! to compute the non-local term in the US scheme.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : okvan, deeq, deeq_nc, dvan, dvan_so
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE noncollin_module,     ONLY : noncolin, domag, nspin_mag, lspinorb
  USE uspp,                 ONLY : nhtol, nhtolm
  USE scf,                  ONLY : v
  USE realus,               ONLY : newq_r
  USE control_flags,        ONLY : tqr
  USE ldaU,                 ONLY : lda_plus_U, Hubbard_projectors
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht
  ! counters and auxiliary variables
  !
  CALL start_clock( 'newd' )
  !
  IF (tqr) THEN
     CALL newq_r(v%of_r,deeq,.false.)
     !$acc update device(deeq)
  ELSE
     CALL newq_acc(v%of_r,deeq,.false.)
  END IF
  !
  call add_paw_to_deeq_acc(deeq)
  !
  types : &
  DO nt = 1, ntyp
     !
     if_noncolin:&
     IF ( noncolin ) THEN
        !
        IF (upf(nt)%has_so) THEN
           !
           CALL newd_so_acc(nt)
           !
        ELSE
           !
           CALL newd_nc_acc(nt)
           !
        END IF
        !
     ELSE if_noncolin
        !
        nht = nh(nt)
        !$acc parallel loop collapse(4) present(ityp,dvan,deeq)
        DO is = 1, nspin
           DO na = 1, nat
              DO ih = 1, nht
                 DO jh = 1, nht
                    IF ( ityp(na) == nt ) THEN
                       deeq(ih,jh,na,is) = deeq(ih,jh,na,is) + dvan(ih,jh,nt)
                    END IF
                 END DO
              END DO
           END DO
        END DO
        !
     END IF if_noncolin
     !
  END DO types
  !
  IF (lda_plus_U .AND. (Hubbard_projectors == 'pseudo')) THEN
    CALL add_vhub_to_deeq_acc(deeq)
  ENDIF
  !
  CALL stop_clock( 'newd' )
  !
  ! ... sync with CPU (not sure this is needed)
  if (noncolin) then
     !$acc update self(deeq_nc)
  else
     !$acc update self(deeq)
  endif
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_so_acc(nt)
      !------------------------------------------------------------------------
      !
      USE upf_spinorb,    ONLY : fcoef
      USE ions_base,      ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER :: nt

      INTEGER :: ijs, is1, is2, kh, lh, nhnt, ih, jh, na
      !
      nhnt = nh(nt)
      ijs = 0
      !
      DO is1 = 1, 2
         !
         DO is2 =1, 2
            !
            ijs = ijs + 1
            !
            IF (domag) THEN
               !$acc parallel loop collapse(3) present(deeq_nc,deeq,fcoef,ityp)
               DO na = 1, nat
                  !
                  DO ih = 1, nhnt
                     !
                     DO jh = 1, nhnt
                        !
                        IF ( ityp(na) == nt ) THEN
                           !
                           deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                           !
                           DO kh = 1, nhnt
                              !
                              DO lh = 1, nhnt
                                 !
                                 deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                      deeq(kh,lh,na,1)*         &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt) + &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt)) + &
                                   deeq(kh,lh,na,2)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt) + &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                   (0.D0,-1.D0)*deeq(kh,lh,na,3)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt) - &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                   deeq(kh,lh,na,4)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt) - &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt))   
                                 !
                              END DO
                              !
                           END DO
                           !
                        END IF
                        !
                     END DO
                  END DO
                  !
               END DO
               !
            ELSE
               !
               !$acc parallel loop collapse(3) present(deeq_nc,deeq,fcoef,ityp)
               DO na = 1, nat
                  !
                  DO ih = 1, nhnt
                     !
                     DO jh = 1, nhnt
                        !
                        IF ( ityp(na) == nt ) THEN
                           !
                           deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                           !
                           DO kh = 1, nhnt
                              !
                              DO lh = 1, nhnt
                                 !
                                 deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) + &
                                      deeq(kh,lh,na,1)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt) + &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt) ) 
                                 !
                              END DO
                              !
                           END DO
                           !
                        END IF
                        !
                     END DO
                     !
                  END DO
                  !
               END DO
               !
            END IF
            !
         END DO
         !
      END DO
      !
    RETURN
      !
    END SUBROUTINE newd_so_acc
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_nc_acc(nt)
      !------------------------------------------------------------------------
      !
      USE ions_base,     ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nt
      INTEGER :: nhnt, na, ih, jh
      !
      nhnt = nh(nt)
      !
      !$acc parallel loop collapse(3) present(deeq_nc,deeq,ityp)
      DO na = 1, nat
         DO ih = 1, nhnt
            DO jh = 1, nhnt
               !
               IF ( ityp(na) == nt ) THEN
                  !
                  IF (lspinorb) THEN
                     deeq_nc(ih,jh,na,1) = dvan_so(ih,jh,1,nt) + &
                                           deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
                     !
                     deeq_nc(ih,jh,na,4) = dvan_so(ih,jh,4,nt) + &
                                           deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
                     !
                  ELSE
                     deeq_nc(ih,jh,na,1) = dvan(ih,jh,nt) + &
                                           deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
                     !
                     deeq_nc(ih,jh,na,4) = dvan(ih,jh,nt) + &
                                           deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
                     !
                  END IF
                  deeq_nc(ih,jh,na,2) = deeq(ih,jh,na,2) - &
                                        ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
                  !
                  deeq_nc(ih,jh,na,3) = deeq(ih,jh,na,2) + &
                                        ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
                  !
               END IF
               !
            END DO
         END DO
      END DO
      !
    RETURN
    END SUBROUTINE newd_nc_acc
    !
  END SUBROUTINE newd_us

END MODULE dfunct
