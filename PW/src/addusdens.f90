!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusdens(rho)
  !----------------------------------------------------------------------
  !! Add US contribution to the charge density to \(\text{rho}(G)\).
  !
  USE realus,               ONLY : addusdens_r
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho(dfftp%ngm,nspin_mag)
  !! Charge density in G space
  !
  IF ( tqr ) THEN
     CALL addusdens_r( rho )
  ELSE
     CALL addusdens_g( rho )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE addusdens
!
!
!----------------------------------------------------------------------
SUBROUTINE addusdens_g(rho)
  !----------------------------------------------------------------------
  !! This routine adds to the charge density \(\text{rho}(G)\) in reciprocal space
  !! the part which is due to the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, gg, g, &
                                   eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : okvan, becsum
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho(dfftp%ngm,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ngm_s, ngm_e, ngm_l, ngm_s_tmp, ngm_e_tmp, ngm_l_tmp
  ! local number of G-vectors, starting/ending indices
  INTEGER :: ig, na, nt, ih, jh, ijh, is, nab, nb, nij, ij, im
  ! counters
  REAL(DP), ALLOCATABLE :: tbecsum(:,:,:)
  ! \sum_kv <\psi_kv|\beta_l><beta_m|\psi_kv> for each species of atoms
  REAL(DP), ALLOCATABLE :: qmod(:), ylmk0(:,:)
  ! modulus of G, spherical harmonics
  COMPLEX(DP), ALLOCATABLE :: skk(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux2(:,:)
  ! structure factors, US contribution to rho
  COMPLEX(DP), ALLOCATABLE :: aux(:,:), qgm(:)
  ! work space for rho(G,nspin), Fourier transform of q
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock( 'addusdens' )
  !
  ALLOCATE( aux(ngm,nspin_mag) )
  !$acc data create(aux)
  !
  !$acc kernels
  aux = (0.d0,0.d0)
  !$acc end kernels
  !
  ! ... With k-point/bgrp parallelization, distribute G-vectors across all processors
  ! ... ngm_s = index of first G-vector for this processor (in the k-point x bgrp pool)
  ! ... ngm_e = index of last  G-vector for this processor (in the k-point x bgrp pool)
  ! ... ngm_l = local number of G-vectors 
  !
  CALL divide( inter_pool_comm, ngm, ngm_s_tmp, ngm_e_tmp )
  ngm_l_tmp = ngm_e_tmp - ngm_s_tmp + 1
  CALL divide( inter_bgrp_comm, ngm_l_tmp, ngm_s, ngm_e )
  ngm_l = ngm_e - ngm_s + 1 
  ngm_s = ngm_s + ngm_s_tmp - 1
  ngm_e = ngm_e + ngm_s_tmp - 1
  !
  ! ... for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GO TO 10
  !
#if !defined(_OPENACC)
  ! ... CPU: G-blocked algorithm, no ngm-sized work arrays (see below)
  CALL addusdens_g_blocked( ngm_s, ngm_l, aux )
  GO TO 10
#endif
  !
  ALLOCATE( qmod(ngm_l), qgm(ngm_l) )
  ALLOCATE( ylmk0(ngm_l,lmaxq*lmaxq) )
  !$acc data create(qmod,qgm,ylmk0)
  !
  CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
  !
#if defined(_OPENACC)
  !$acc parallel loop
#else
  !$omp parallel do
#endif
  DO ig = 1, ngm_l
     qmod(ig) = SQRT(gg(ngm_s+ig-1))*tpiba
  ENDDO
  !
  ! ... [largest size for buffer: nij = nhm*(nhm+1)/2]
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! ... nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ! ... count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na)==nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE( skk(ngm_l,nab), tbecsum(nij,nab,nspin_mag), aux2(ngm_l,nij) )
        !$acc data create(skk,tbecsum,aux2)
        !
        CALL start_clock( 'addusd:skk' )
        nb = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) THEN
              nb = nb + 1
              !tbecsum(:,nb,:) = becsum(1:nij,na,1:nspin_mag)
              !
              !$acc parallel loop collapse(2) present(becsum)
              DO im = 1, nspin_mag
                 DO ij = 1, nij
                   tbecsum(ij,nb,im) = becsum(ij,na,im)
                 ENDDO
              ENDDO
              !
#if defined(_OPENACC)
              !$acc parallel loop present(eigts1,eigts2,eigts3,mill)
#else
              !$omp parallel do
#endif
              DO ig = 1, ngm_l
                 skk(ig,nb) = eigts1(mill(1,ngm_s+ig-1),na) * &
                              eigts2(mill(2,ngm_s+ig-1),na) * &
                              eigts3(mill(3,ngm_s+ig-1),na)
              ENDDO
           ENDIF
        ENDDO
        CALL stop_clock( 'addusd:skk' )
        !
        DO is = 1, nspin_mag
           ! ... sum over atoms
           !
           !$acc host_data use_device(skk,tbecsum,aux2)
           CALL MYDGEMM( 'N', 'T', 2*ngm_l, nij, nab, 1.0_dp, skk, 2*ngm_l, &
                         tbecsum(1,1,is), nij, 0.0_dp, aux2, 2*ngm_l )
           !$acc end host_data
           !
           ! ... sum over lm indices of Q_{lm}
           ijh = 0
           DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                 ijh = ijh + 1
                 CALL qvan2( ngm_l, ih, jh, nt, qmod, qgm, ylmk0 )
#if defined(_OPENACC)
                 !$acc parallel loop
#else
                 !$omp parallel do
#endif
                 DO ig = 1, ngm_l
                    aux(ngm_s+ig-1,is) = aux(ngm_s+ig-1,is) + aux2(ig,ijh)*qgm(ig)
                 ENDDO
                 !
             ENDDO
           ENDDO
        ENDDO
        !
        !$acc end data
        DEALLOCATE( tbecsum, skk, aux2 )
        !
     ENDIF
  ENDDO
  !
  !$acc end data
  DEALLOCATE( ylmk0 )
  DEALLOCATE( qgm, qmod )
  !
  10 CONTINUE
  !
  !$acc update self(aux)
  !$acc end data
  !
  CALL mp_sum( aux, inter_bgrp_comm )
  CALL mp_sum( aux, inter_pool_comm )
  !
  ! ... add aux to the charge density in reciprocal space
  !
  rho(:,:) = rho(:,:) + aux(:,:)
  !
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'addusdens' )
  !
  RETURN
  !
END SUBROUTINE addusdens_g
!
#if !defined(_OPENACC)
!----------------------------------------------------------------------
SUBROUTINE addusdens_g_blocked( ngm_s, ngm_l, aux )
  !----------------------------------------------------------------------
  !! CPU version of the G-space augmentation charge, processing the local
  !! G-vectors ngm_s ... ngm_s+ngm_l-1 in blocks of \(\text{gblk}\) owned by
  !! one OpenMP thread each. For every block and species:
  !
  !! - structure factors \(S_a(G)\) of the atoms of that species;
  !! - \(B_{ij}(G) = \sum_a S_a(G)\,\text{becsum}_{ij,a}\) (one small DGEMM
  !!   per spin);
  !! - spherical harmonics, \(|G|\) and \(Q_{ij}(G)\) (qvan2) for the block,
  !!   computed once and used for all spin components;
  !! - \(\rho_{aug}(G) \mathrel{+}= \sum_{ij} B_{ij}(G)\,Q_{ij}(G)\).
  !
  !! Compared with the original algorithm this removes the work arrays of
  !! size ngm_l*nij (aux2), ngm_l*nat (skk) and ngm_l*lmaxq^2 (ylmk0): the
  !! extra memory is O(gblk*(nij*nspin+nat+lmaxq^2)) per thread. All
  !! G-dependent work is threaded, and qvan2 is no longer repeated for each
  !! spin component. For each G the ij terms are accumulated in the original
  !! order.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : tpiba
  USE gvect,                ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : nspin_mag
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : upf, lmaxq, nh
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngm_s
  !! global index of the first local G-vector
  INTEGER, INTENT(IN) :: ngm_l
  !! number of local G-vectors
  COMPLEX(DP), INTENT(INOUT) :: aux(ngm,nspin_mag)
  !! augmentation charge, accumulated
  !
  INTEGER, PARAMETER :: gblk = 512
  !
  INTEGER :: nt, na, nb, nab, nij, ig0, n, ig, igg, ih, jh, ijh, is
  INTEGER,     ALLOCATABLE :: atoms(:)
  REAL(DP),    ALLOCATABLE :: tbecsum(:,:,:), qmod(:), ylm(:)
  COMPLEX(DP), ALLOCATABLE :: skk(:,:), aux2(:,:,:), qg(:)
  !
  DO nt = 1, ntyp
     !
     IF ( .NOT. upf(nt)%tvanp ) CYCLE
     nab = COUNT( ityp(1:nat) == nt )
     IF ( nab == 0 ) CYCLE
     nij = nh(nt)*(nh(nt)+1)/2
     !
     ALLOCATE( atoms(nab), tbecsum(nij,nab,nspin_mag) )
     nb = 0
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           nb = nb + 1
           atoms(nb) = na
           tbecsum(:,nb,:) = becsum(1:nij,na,1:nspin_mag)
        ENDIF
     ENDDO
     !
     !$omp parallel default(shared) &
     !$omp   private(ig0,n,igg,ig,nb,na,is,ih,jh,ijh,skk,aux2,qg,qmod,ylm)
     ALLOCATE( skk(gblk,nab), aux2(gblk,nij,nspin_mag), qg(gblk), qmod(gblk) )
     ! ... ylm is 1D: it is passed (sequence association) as ylm(n,lmaxq**2)
     ! ... with leading dimension equal to the actual block size n
     ALLOCATE( ylm(gblk*lmaxq*lmaxq) )
     !
     !$omp do schedule(static)
     DO ig0 = 1, ngm_l, gblk
        !
        n   = MIN( gblk, ngm_l - ig0 + 1 )
        igg = ngm_s + ig0 - 1           ! global index of the first G of the block
        !
        ! ... structure factors of all atoms of this species
        DO nb = 1, nab
           na = atoms(nb)
           DO ig = 1, n
              skk(ig,nb) = eigts1(mill(1,igg+ig-1),na) * &
                           eigts2(mill(2,igg+ig-1),na) * &
                           eigts3(mill(3,igg+ig-1),na)
           ENDDO
        ENDDO
        !
        ! ... B_ij(G) = sum over atoms (real becsum: complex as 2 reals)
        DO is = 1, nspin_mag
           CALL DGEMM( 'N', 'T', 2*n, nij, nab, 1.0_DP, skk, 2*gblk, &
                       tbecsum(1,1,is), nij, 0.0_DP, aux2(1,1,is), 2*gblk )
        ENDDO
        !
        ! ... |G| and Y_lm(G) for this block
        CALL ylmr2( lmaxq*lmaxq, n, g(1,igg), gg(igg), ylm )
        DO ig = 1, n
           qmod(ig) = SQRT( gg(igg+ig-1) ) * tpiba
        ENDDO
        !
        ! ... sum over ij of B_ij(G) Q_ij(G); Q_ij computed once for all spins
        ijh = 0
        DO ih = 1, nh(nt)
           DO jh = ih, nh(nt)
              ijh = ijh + 1
              CALL qvan2( n, ih, jh, nt, qmod, qg, ylm )
              DO is = 1, nspin_mag
                 DO ig = 1, n
                    aux(igg+ig-1,is) = aux(igg+ig-1,is) + aux2(ig,ijh,is)*qg(ig)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        !
     ENDDO
     !$omp end do
     !
     DEALLOCATE( skk, aux2, qg, qmod, ylm )
     !$omp end parallel
     !
     DEALLOCATE( atoms, tbecsum )
     !
  ENDDO
  !
END SUBROUTINE addusdens_g_blocked
#endif
