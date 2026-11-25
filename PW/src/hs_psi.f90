!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE hs_psi_fused
  !----------------------------------------------------------------------------
  !! Control data for the fused application of the nonlocal part of H and of
  !! S-1 (ultrasoft/PAW, k-point version, CPU, G space): see hs_psi.
  !
  IMPLICIT NONE
  SAVE
  !
  LOGICAL :: defer_vnl = .FALSE.
  !! if .TRUE. h_psi_ computes becp but does NOT add V_NL|psi> to hpsi:
  !! the caller (hs_psi) does it together with S|psi> in add_vuspsi_spsi_k
  !
CONTAINS
  !
  LOGICAL FUNCTION hs_fused_ok()
    !! .TRUE. if the fused path is implemented for the current setup
    USE uspp,             ONLY : nkb, okvan
    USE control_flags,    ONLY : gamma_only, use_gpu, scissor
    USE noncollin_module, ONLY : noncolin
    USE realus,           ONLY : real_space
    USE bp,               ONLY : lelfield
    USE xc_lib,           ONLY : exx_is_active
    !
    hs_fused_ok = okvan .AND. nkb > 0 .AND. .NOT. ( gamma_only .OR. noncolin &
                  .OR. real_space .OR. use_gpu .OR. scissor .OR. lelfield    &
                  .OR. exx_is_active() )
    !
  END FUNCTION hs_fused_ok
  !
END MODULE hs_psi_fused
!
!----------------------------------------------------------------------------
SUBROUTINE hs_psi( lda, n, m, psi, hpsi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian and the S matrix
  ! ... to a vector psi and puts the result in hpsi and spsi
  ! ... Wrapper routine - calls h_psi_ and s_psi_
  !
  ! ... No bgrp parallelization here !
  !
  ! ... For ultrasoft/PAW with k-points (CPU, G-space projectors) the nonlocal
  ! ... part of H|psi> and (S-1)|psi> are applied in a single sweep over vkb
  !
  USE kinds,  ONLY: DP
  USE noncollin_module, ONLY: npol
  USE hs_psi_fused,     ONLY: defer_vnl, hs_fused_ok
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX (DP) :: psi(lda*npol, m), hpsi(lda*npol,m), spsi(lda*npol,m)
  !
  CALL start_clock( 'hs_psi' )
  !
  IF ( hs_fused_ok() ) THEN
     !
     defer_vnl = .TRUE.
     CALL h_psi_ ( lda, n, m, psi, hpsi ) ! kinetic + local + becp; V_NL deferred
     defer_vnl = .FALSE.
     CALL add_vuspsi_spsi_k ( lda, n, m, psi, hpsi, spsi )
     !
  ELSE
     !
     CALL h_psi_ ( lda, n, m, psi, hpsi ) ! apply H to m wfcs (no bgrp parallelization here)
     CALL s_psi_ ( lda, n, m, psi, spsi ) ! apply S to m wfcs (no bgrp parallelization here)
     !
  ENDIF
  !
  CALL stop_clock( 'hs_psi' )
  !
  RETURN
  !
END SUBROUTINE hs_psi
!
!----------------------------------------------------------------------------
SUBROUTINE add_vuspsi_spsi_k( lda, n, m, psi, hpsi, spsi )
  !----------------------------------------------------------------------------
  !! Fused version of add_vuspsi_k + s_psi_k (collinear, k-points, CPU):
  !!   hpsi <- hpsi + sum_ij |beta_i> D_ij <beta_j|psi>
  !!   spsi <- psi  + sum_ij |beta_i> q_ij <beta_j|psi>
  !! Requires becp%k = <beta|psi> (computed by h_psi_).
  !! The D- and q-weighted coefficients are stored side by side in ps(nkb,2m)
  !! and vkb is swept once, block by block over G, so that each block of vkb
  !! is read from memory once and reused from cache for the two products.
  !! No work array scaling with npw is needed.
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nat, ntyp => nsp, ityp
  USE lsda_mod,         ONLY: current_spin
  USE uspp,             ONLY: vkb, nkb, deeq, qq_at, ofsbeta
  USE uspp_param,       ONLY: upf, nh
  USE becmod,           ONLY: becp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN)    :: psi(lda,m)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  COMPLEX(DP), INTENT(OUT)   :: spsi(lda,m)
  !
  COMPLEX(DP), PARAMETER :: one = (1.0_DP,0.0_DP), zero = (0.0_DP,0.0_DP)
  COMPLEX(DP), ALLOCATABLE :: ps(:,:), aux(:,:)
  INTEGER :: na, nt, ig0, nb, nstep, ibnd, ierr
  INTEGER, SAVE :: nblock = -1
  CHARACTER(LEN=16) :: envval
  !
  CALL start_clock( 'vus_s_psi' )
  !
  IF ( nblock < 0 ) THEN
     ! rows of vkb per block (default 256: smaller blocks are slower with MKL ZGEMM)
     nblock = 0
     CALL get_environment_variable( 'QE_HS_BLOCK', envval, STATUS=ierr )
     IF ( ierr == 0 ) READ( envval, *, IOSTAT=ierr ) nblock
     IF ( ierr /= 0 ) nblock = 0
  ENDIF
  IF ( nblock > 0 ) THEN
     nstep = nblock
  ELSE
     nstep = 256
  ENDIF
  !
  ALLOCATE( ps(nkb,2*m), STAT=ierr )
  IF ( ierr /= 0 ) CALL errore( 'add_vuspsi_spsi_k', 'cannot allocate ps', ABS(ierr) )
  !
  DO nt = 1, ntyp
     IF ( nh(nt) == 0 ) CYCLE
     ALLOCATE( aux(nh(nt),nh(nt)) )
     DO na = 1, nat
        IF ( ityp(na) /= nt ) CYCLE
        aux(:,:) = CMPLX( deeq(1:nh(nt),1:nh(nt),na,current_spin), 0.0_DP, KIND=DP )
        CALL ZGEMM( 'N', 'N', nh(nt), m, nh(nt), one, aux, nh(nt), &
                    becp%k(ofsbeta(na)+1,1), nkb, zero, ps(ofsbeta(na)+1,1), nkb )
        IF ( upf(nt)%tvanp ) THEN
           aux(:,:) = CMPLX( qq_at(1:nh(nt),1:nh(nt),na), 0.0_DP, KIND=DP )
           CALL ZGEMM( 'N', 'N', nh(nt), m, nh(nt), one, aux, nh(nt), &
                       becp%k(ofsbeta(na)+1,1), nkb, zero, ps(ofsbeta(na)+1,m+1), nkb )
        ELSE
           ps(ofsbeta(na)+1:ofsbeta(na)+nh(nt),m+1:2*m) = zero
        ENDIF
     ENDDO
     DEALLOCATE( aux )
  ENDDO
  !
  !$omp parallel do private(nb, ibnd) schedule(static)
  DO ig0 = 1, n, nstep
     nb = MIN( nstep, n - ig0 + 1 )
     DO ibnd = 1, m
        spsi(ig0:ig0+nb-1,ibnd) = psi(ig0:ig0+nb-1,ibnd)
     ENDDO
     CALL ZGEMM( 'N', 'N', nb, m, nkb, one, vkb(ig0,1), lda, ps(1,1), nkb, &
                 one, hpsi(ig0,1), lda )
     CALL ZGEMM( 'N', 'N', nb, m, nkb, one, vkb(ig0,1), lda, ps(1,m+1), nkb, &
                 one, spsi(ig0,1), lda )
  ENDDO
  !$omp end parallel do
  IF ( lda > n ) spsi(n+1:lda,1:m) = psi(n+1:lda,1:m)
  !
  DEALLOCATE( ps )
  !
  CALL stop_clock( 'vus_s_psi' )
  !
END SUBROUTINE add_vuspsi_spsi_k
