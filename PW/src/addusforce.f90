!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusforce( forcenl )
  !----------------------------------------------------------------------
  !! Wrapper to \(\texttt{addusforce_g}\) or \(\texttt{addusforce_r}\).
  !
  USE kinds,         ONLY : dp
  USE ions_base,     ONLY : nat
  USE control_flags, ONLY : tqr
  USE realus,        ONLY : addusforce_r
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3,nat)
  !! the non-local contribution to the force
  !
  IF ( tqr ) THEN
     CALL addusforce_r( forcenl )
  ELSE
     CALL addusforce_g( forcenl )
  ENDIF
  !
END SUBROUTINE addusforce
!
!----------------------------------------------------------------------
SUBROUTINE addusforce_g( forcenl )
  !----------------------------------------------------------------------
  !! This routine computes the contribution to atomic forces due
  !! to the dependence of the Q function on the atomic position.
  !! \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
  !!    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
  !! where:
  !! \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
  !!    w_i\langle \beta_m|\psi_i\rangle \]
  !! On output: the contribution is added to \(\text{forcenl}\).
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE cell_base,          ONLY : omega, tpiba
  USE fft_base,           ONLY : dfftp
  USE fft_rho,            ONLY : rho_r2g
  USE gvect,              ONLY : ngm, gg, g, eigts1, eigts2, eigts3,&
                                 mill
  USE noncollin_module,   ONLY : nspin_mag
  USE scf,                ONLY : v, vltot
  USE uspp,               ONLY : becsum, okvan
  USE uspp_param,         ONLY : upf, lmaxq, nh, nhm
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE mp_pools,           ONLY : inter_pool_comm
  USE mp,                 ONLY : mp_sum
  USE control_flags,      ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: forcenl(3,nat)
  !! the non-local contribution to the force
  !
  ! ... local variables
  !
  INTEGER :: ngm_s, ngm_e, ngm_l
  INTEGER :: ig, nt, ih, jh, ijh, nij, ipol, is, na, nb, nab, ir
  REAL(DP), ALLOCATABLE :: forceq(:,:)
  REAL(DP), ALLOCATABLE :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:)
  REAL(DP) :: fact, forceqx, forceqy, forceqz
  COMPLEX(DP) :: cfac
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:,:), vg(:,:), qgm(:,:)
  !
  IF (.NOT.okvan) RETURN
  !
  ALLOCATE( forceq(3,nat) )
  forceq(:,:) = 0._DP
  IF ( gamma_only ) THEN
     fact = 2._DP*omega
  ELSE
     fact = 1._DP*omega
  ENDIF
  !
  ! ... Fourier transform of the total effective potential
  !
  ALLOCATE( vg(ngm,nspin_mag) )
  !$acc data copyin( v, vltot ) create( vg ) 
  !$acc data copyin( v%of_r )
  !
  DO is = 1, nspin_mag
    IF ( nspin_mag == 4 .and. is /= 1 ) THEN
       CALL rho_r2g( dfftp, v%of_r(:,is), vg(:,is:is) )
    ELSE
       CALL rho_r2g( dfftp, v%of_r(:,is), vg(:,is:is), v=vltot )
    ENDIF
    ! ... Note the factors -i and 2pi/a *units of G) here in V(G) !
    !$acc kernels
    vg(:,is) = vg(:,is) * tpiba * (0.d0,-1.d0)
    !$acc end kernels
  ENDDO
  !
  ! ... With k-point parallelization, distribute G-vectors across processors
  ! ... ngm_s = index of first G-vector for this processor
  ! ... ngm_e = index of last  G-vector for this processor
  ! ... ngm_l = local number of G-vectors
  !
  CALL divide( inter_pool_comm, ngm, ngm_s, ngm_e )
  ngm_l = ngm_e-ngm_s+1
  !
  ! ... for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GO TO 10
  !
#if !defined(_OPENACC)
  ! ... CPU: G-blocked algorithm, no ngm-sized Q(G) or iG*V(G)*S(G) arrays
  CALL addusforce_g_blocked( ngm_s, ngm_l, vg, fact, forceq )
  GO TO 10
#endif
  !
  ALLOCATE( ylmk0(ngm_l,lmaxq*lmaxq), qmod(ngm_l) )
  !$acc data create(ylmk0,qmod)
  !
  CALL ylmr2( lmaxq*lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0 )
  !
  !$acc parallel loop
  DO ig = 1, ngm_l
     qmod(ig) = SQRT( gg(ngm_s+ig-1) )*tpiba
  ENDDO
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! ... nij = max number of (ih,jh) pairs per atom type nt
        ! ... qgm contains the Q functions in G space
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ALLOCATE( qgm(ngm_l,nij) )
        !$acc data create(qgm)
        !
        ijh = 0
        DO ih = 1, nh(nt)
           DO jh = ih, nh(nt)
              ijh = ijh + 1
              CALL qvan2( ngm_l, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
           ENDDO
        ENDDO
        !
        ! ... nab = number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na)==nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE( aux1(ngm_l,nab,3), ddeeq(nij,nab,3,nspin_mag) )
        !$acc data create(ddeeq)
        !$acc data create(aux1)
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 nb = nb + 1
                 !
                 ! ... aux1 = product of potential, structure factor and iG
                 !
#if defined(_OPENACC)
!$acc parallel loop present(eigts1,eigts2,eigts3,mill,g)
#else
!$omp parallel do default(shared) private(ig,cfac)
#endif
                 DO ig = 1, ngm_l
                    cfac = vg(ngm_s+ig-1,is) * &
                           CONJG(eigts1(mill(1,ngm_s+ig-1),na) * &
                                 eigts2(mill(2,ngm_s+ig-1),na) * &
                                 eigts3(mill(3,ngm_s+ig-1),na) )
                    aux1(ig,nb,1) = g(1,ngm_s+ig-1) * cfac
                    aux1(ig,nb,2) = g(2,ngm_s+ig-1) * cfac
                    aux1(ig,nb,3) = g(3,ngm_s+ig-1) * cfac
                 ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
                 !
              ENDIF
           ENDDO
           !
           ! ... ddeeq = dot product of aux1 with the Q functions
           ! ... No need for special treatment of the G=0 term (is zero)
           !
           !$acc host_data use_device(qgm,aux1,ddeeq)
           DO ipol = 1, 3
              CALL MYDGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm, 2*ngm_l, &
                            aux1(1,1,ipol), 2*ngm_l, 0.0_DP, ddeeq(1,1,ipol,is), nij )
           ENDDO
           !$acc end host_data
           !
        ENDDO
        !
        !$acc end data
        DEALLOCATE( aux1 )
        !
        DO is = 1, nspin_mag
           nb = 0
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 nb = nb + 1
                 forceqx = 0
                 forceqy = 0
                 forceqz = 0
                 !$acc parallel loop reduction(+:forceqx,forceqy,forceqz) present(becsum)
                 DO ijh = 1, nij
                   forceqx = forceqx + ddeeq(ijh,nb,1,is) * becsum(ijh,na,is)
                   forceqy = forceqy + ddeeq(ijh,nb,2,is) * becsum(ijh,na,is)
                   forceqz = forceqz + ddeeq(ijh,nb,3,is) * becsum(ijh,na,is)
                 ENDDO
                 forceq(1,na) = forceq(1,na) + forceqx
                 forceq(2,na) = forceq(2,na) + forceqy
                 forceq(3,na) = forceq(3,na) + forceqz
              ENDIF
           ENDDO
        ENDDO
        !$acc end data
        DEALLOCATE( ddeeq )
        !$acc end data
        DEALLOCATE( qgm )
     ENDIF
  ENDDO
  !
  !$acc end data
  DEALLOCATE( ylmk0, qmod )
  !
  10 CONTINUE
  CALL mp_sum( forceq, inter_pool_comm )
  CALL mp_sum( forceq, intra_bgrp_comm )
  !
  forcenl(:,:) = forcenl(:,:) + forceq(:,:)
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( vg )
  DEALLOCATE( forceq )
  !
  RETURN
  !
END SUBROUTINE addusforce_g
!
#if !defined(_OPENACC)
!----------------------------------------------------------------------
SUBROUTINE addusforce_g_blocked( ngm_s, ngm_l, vg, fact, forceq )
  !----------------------------------------------------------------------
  !! CPU version of the augmentation contribution to the forces,
  !! \[ F_{j,a} = \sum_{ij} \text{becsum}_{ij,a} \sum_G Q^*_{ij}(G)\,
  !!    iG_j V(G) S^*_a(G) \]
  !! computed on blocks of \(\text{gblk}\) local G-vectors, each owned by one
  !! OpenMP thread. Per block and species: Y_lm(G), |G|, Q_ij(G) (qvan2,
  !! once for all spin components and directions), \(G_j V(G) S^*_a(G)\)
  !! and one DGEMM per direction accumulating into per-thread partial sums
  !! ddeeq(ij,atom,direction,spin), which are reduced in thread order
  !! (deterministic for a given number of threads).
  !
  !! Compared with the original algorithm this removes the work arrays
  !! qgm(ngm_l,nij), aux1(ngm_l,nat,3) and ylmk0(ngm_l,lmaxq^2); extra memory
  !! is O(gblk*(nij+3*nat+lmaxq^2)) per thread plus 3*nij*nat*nspin per
  !! thread. Results differ from the original only by the order of the
  !! floating-point summation over G (rounding level).
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE cell_base,          ONLY : tpiba
  USE gvect,              ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,   ONLY : nspin_mag
  USE uspp,               ONLY : becsum
  USE uspp_param,         ONLY : upf, lmaxq, nh
!$ USE omp_lib,           ONLY : omp_get_max_threads, omp_get_thread_num
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: ngm_s, ngm_l
  !! global index of the first local G-vector, number of local G-vectors
  COMPLEX(DP), INTENT(IN)    :: vg(ngm,nspin_mag)
  !! -i*tpiba*V_eff(G) on all G-vectors of this process
  REAL(DP),    INTENT(IN)    :: fact
  !! omega (2*omega with gamma tricks)
  REAL(DP),    INTENT(INOUT) :: forceq(3,nat)
  !! augmentation force, accumulated
  !
  INTEGER, PARAMETER :: gblk = 512
  !
  INTEGER :: nt, na, nb, nab, nij, ig0, n, ig, igg, ih, jh, ijh, is, ipol, nth, ith, it
  INTEGER,     ALLOCATABLE :: atoms(:)
  REAL(DP),    ALLOCATABLE :: ddeeq(:,:,:,:,:), qmod(:), ylm(:)
  COMPLEX(DP), ALLOCATABLE :: qgm(:,:), aux1(:,:,:)
  COMPLEX(DP) :: cfac
  REAL(DP)    :: fx, fy, fz
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
     ALLOCATE( ddeeq(nij,nab,3,nspin_mag,0:nth-1) )
     ddeeq(:,:,:,:,:) = 0.0_DP
     !
     !$omp parallel default(shared) &
     !$omp   private(ith,ig0,n,igg,ig,nb,na,is,ipol,ih,jh,ijh,cfac,qgm,aux1,qmod,ylm)
     ith = 0
!$   ith = omp_get_thread_num()
     ALLOCATE( qgm(gblk,nij), aux1(gblk,nab,3), qmod(gblk) )
     ! ... ylm is 1D: passed (sequence association) as ylm(n,lmaxq**2)
     ALLOCATE( ylm(gblk*lmaxq*lmaxq) )
     !
     !$omp do schedule(static)
     DO ig0 = 1, ngm_l, gblk
        !
        n   = MIN( gblk, ngm_l - ig0 + 1 )
        igg = ngm_s + ig0 - 1           ! global index of the first G of the block
        !
        ! ... Q_ij(G) for this block
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
           ! ... aux1 = product of potential, structure factor and iG
           DO nb = 1, nab
              na = atoms(nb)
              DO ig = 1, n
                 cfac = vg(igg+ig-1,is) * &
                        CONJG( eigts1(mill(1,igg+ig-1),na) * &
                               eigts2(mill(2,igg+ig-1),na) * &
                               eigts3(mill(3,igg+ig-1),na) )
                 aux1(ig,nb,1) = g(1,igg+ig-1) * cfac
                 aux1(ig,nb,2) = g(2,igg+ig-1) * cfac
                 aux1(ig,nb,3) = g(3,igg+ig-1) * cfac
              ENDDO
           ENDDO
           ! ... partial dot products with Q_ij(G) over this block
           ! ... (no special treatment of G=0: the iG factor makes it zero)
           DO ipol = 1, 3
              CALL DGEMM( 'C', 'N', nij, nab, 2*n, fact, qgm, 2*gblk, &
                          aux1(1,1,ipol), 2*gblk, 1.0_DP, ddeeq(1,1,ipol,is,ith), nij )
           ENDDO
        ENDDO
        !
     ENDDO
     !$omp end do
     !
     DEALLOCATE( qgm, aux1, qmod, ylm )
     !$omp end parallel
     !
     ! ... deterministic reduction of the per-thread partial sums
     DO it = 1, nth-1
        ddeeq(:,:,:,:,0) = ddeeq(:,:,:,:,0) + ddeeq(:,:,:,:,it)
     ENDDO
     !
     ! ... contract with becsum (same order as addusforce_g)
     DO is = 1, nspin_mag
        DO nb = 1, nab
           na = atoms(nb)
           fx = 0.0_DP ; fy = 0.0_DP ; fz = 0.0_DP
           DO ijh = 1, nij
              fx = fx + ddeeq(ijh,nb,1,is,0) * becsum(ijh,na,is)
              fy = fy + ddeeq(ijh,nb,2,is,0) * becsum(ijh,na,is)
              fz = fz + ddeeq(ijh,nb,3,is,0) * becsum(ijh,na,is)
           ENDDO
           forceq(1,na) = forceq(1,na) + fx
           forceq(2,na) = forceq(2,na) + fy
           forceq(3,na) = forceq(3,na) + fz
        ENDDO
     ENDDO
     !
     DEALLOCATE( ddeeq, atoms )
     !
  ENDDO
  !
END SUBROUTINE addusforce_g_blocked
#endif
