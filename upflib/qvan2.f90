!
! Copyright (C) 2023 Quantum ESPRESSO  Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE qvan2( ngy, ih, jh, np, qmod, qg, ylmk0 )
  !-----------------------------------------------------------------------
  !! This routine computes the Fourier transform of the Q functions.
  !
  !! The interpolation table for the radial Fourier transform is stored 
  !! in qrad.
  !
  !! The formula implemented here is:
  !! \[ q(g,i,j) = \sum_\text{lm} (-i)^l \text{ap}(\text{lm},i,j) 
  !! \text{yr}_\text{lm}(g) \text{qrad}(g,l,i,j) \]
  !
  USE upf_kinds,   ONLY: DP
  USE qrad_mod,    ONLY: dq, tab_qrad
  USE uspp_param,  ONLY: lmaxq, nbetam
  USE uspp,        ONLY: nlx, lpl, lpx, ap, indv, nhtolm
!$ USE omp_lib,    ONLY: omp_in_parallel, omp_get_num_threads, omp_get_thread_num
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngy
  !! number of G vectors to compute
  INTEGER, INTENT(IN) :: ih
  !! first index of Q
  INTEGER, INTENT(IN) :: jh
  !! second index of Q
  INTEGER, INTENT(IN) :: np
  !! index of pseudopotentials
  REAL(DP), INTENT(IN) :: ylmk0(ngy,lmaxq*lmaxq)
  !! spherical harmonics
  REAL(DP), INTENT(IN) :: qmod(ngy)
  !! moduli of the q+g vectors 
  REAL(DP), INTENT(OUT) :: qg(2,ngy)
  !! the Fourier transform of interest
  !
  ! ... local variables
  !
  REAL(DP) :: sig
  ! the nonzero real or imaginary part of (-i)^L
  !
  REAL(DP), PARAMETER :: sixth = 1.0_DP / 6.0_DP
  !
  INTEGER :: nb, mb, ijv, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3, ind
  ! nb,mb  : atomic index corresponding to ih,jh
  ! ijv    : combined index (nb,mb)
  ! ivl,jvl: combined LM index corresponding to ih,jh
  ! ig     : counter on g vectors
  ! lp     : combined LM index
  ! l-1    is the angular momentum L
  ! lm     : all possible LM's compatible with ih,jh
  ! i0-i3  : counters for interpolation table
  ! ind    : ind=1 if the results is real (l even), ind=2 if complex (l odd)
  !
  REAL(DP) :: dqi, qm, px, ux, vx, wx, uvx, pwx, work
#if !defined(_OPENACC)
  INTEGER, PARAMETER :: maxlm = 64
  INTEGER :: nlm, lp_v(maxlm), l_v(maxlm), ind_v(maxlm)
  INTEGER :: nth, ith, chunk, ig1, ig2
  REAL(DP) :: sig_v(maxlm)
  LOGICAL :: use_omp
#endif
  ! 1 divided dq
  ! qmod/dq
  ! measures for interpolation table
  ! auxiliary variables for intepolation
  ! auxiliary variables
  !
!$acc data present_or_copyin(qmod,ylmk0) present_or_copyout(qg) present(tab_qrad)
  !
  ! This should not happen, but better to check
  ! FIXME: why is the following not working?
  !IF ( INT(qmod(ngy)/dq)+4 > size(qrad,1) ) CALL upf_error &
  !     ('qvan2', 'internal error: dimension of interpolation table', 1 )
  !
  ! ... computes the indices which correspond to ih,jh
  !
  nb = indv(ih,np)
  mb = indv(jh,np)
  !
  IF (nb >= mb) THEN
     ijv = nb * (nb - 1) / 2 + mb
  ELSE
     ijv = mb * (mb - 1) / 2 + nb
  ENDIF
  !
  ivl = nhtolm(ih,np)
  jvl = nhtolm(jh,np)
  !
  IF (nb > nbetam .OR. mb > nbetam) &
       CALL upf_error( ' qvan2 ', ' wrong dimensions (1)', MAX(nb,mb) )
  IF (ivl > nlx .OR. jvl > nlx) &
       CALL upf_error( ' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl) )
  !
  ! ... and makes the sum over the non zero LM
  !
  dqi = 1.0_DP / dq
#if !defined(_OPENACC)
  !
  ! ... CPU path: loop over G outermost (OpenMP), LM innermost. The cubic
  ! ... interpolation weights depend only on |G| and are computed once per G
  ! ... instead of once per (G,LM). The LM sum for each G is accumulated in
  ! ... the same order as the original code, so results are bit-identical.
  !
  nlm = lpx(ivl,jvl)
  IF ( nlm > maxlm ) CALL upf_error( 'qvan2', ' too many LM components ', nlm )
  DO lm = 1, nlm
     lp = lpl(ivl,jvl,lm)
     IF ( lp < 1 .OR. lp > 49 ) CALL upf_error( 'qvan2', ' lp wrong ', MAX(lp,1) )
     lp_v(lm) = lp
     CALL qvan2_lsig( lp, l_v(lm), sig, ind_v(lm) )
     sig_v(lm) = sig * ap(lp, ivl, jvl)
  ENDDO
  !
  ! ... Threads are spawned only for large G sets and only when not already
  ! ... inside a parallel region: blocked callers (addusdens, newq) call
  ! ... qvan2 from within their own OpenMP region on small G blocks, where a
  ! ... nested (or even an inactive) parallel region per call would cost more
  ! ... than the work itself.
  !
  use_omp = .FALSE.
!$ use_omp = ( ngy >= 2048 ) .AND. ( .NOT. omp_in_parallel() )
  IF ( use_omp ) THEN
!$omp parallel private(nth,ith,chunk,ig1,ig2)
     nth = 1 ; ith = 0
!$   nth = omp_get_num_threads()
!$   ith = omp_get_thread_num()
     chunk = ( ngy + nth - 1 ) / nth
     ig1 = ith*chunk + 1
     ig2 = MIN( ngy, ig1 + chunk - 1 )
     IF ( ig1 <= ig2 ) CALL qvan2_kernel( ig1, ig2 )
!$omp end parallel
  ELSE
     CALL qvan2_kernel( 1, ngy )
  ENDIF
  !
#else
  !$acc kernels
  qg = 0.0_DP
  !$acc end kernels
  !
  DO lm = 1, lpx(ivl,jvl)
     lp = lpl(ivl,jvl,lm)
     IF ( lp < 1 .OR. lp > 49 ) CALL upf_error( 'qvan2', ' lp wrong ', MAX(lp,1) )
     !
     ! ... finds angular momentum l corresponding to combined index lp (l is 
     !     actually l+1 because this is the way qrad is stored, check init_us_1)
     !
     IF ( lp == 1 ) THEN
        l = 1
        sig = 1.0_DP
        ind = 1
     ELSEIF ( lp <=  4 ) THEN
        l = 2
        sig =-1.0_DP
        ind = 2
     ELSEIF ( lp <=  9 ) THEN
        l = 3
        sig =-1.0_DP
        ind = 1
     ELSEIF ( lp <= 16 ) THEN
        l = 4
        sig = 1.0_DP
        ind = 2
     ELSEIF ( lp <= 25 ) THEN
        l = 5
        sig = 1.0_DP
        ind = 1
     ELSEIF ( lp <= 36 ) THEN
        l = 6
        sig =-1.0_DP
        ind = 2
     ELSE
        l = 7
        sig =-1.0_DP
        ind = 1
     ENDIF
     !
     sig = sig * ap(lp, ivl, jvl)
     !
!$acc parallel loop
     DO ig = 1, ngy
        !
        qm = qmod (ig) * dqi
        px = qm - INT(qm)
        ux = 1.0_DP - px
        vx = 2.0_DP - px
        wx = 3.0_DP - px
        i0 = INT(qm) + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        uvx = ux * vx * sixth
        pwx = px * wx * 0.5_DP
        work = tab_qrad(i0,ijv,l,np) * uvx * wx + &
               tab_qrad(i1,ijv,l,np) * pwx * vx - &
               tab_qrad(i2,ijv,l,np) * pwx * ux + &
               tab_qrad(i3,ijv,l,np) * px * uvx
        qg(ind,ig) = qg(ind,ig) + sig * ylmk0(ig,lp) * work
        !
     ENDDO
  ENDDO
#endif
!$acc end data
  !
  !
  RETURN
  !
CONTAINS
  !
#if !defined(_OPENACC)
  SUBROUTINE qvan2_kernel( ig1, ig2 )
    !! Q_ij(G) for G-vectors ig1..ig2: G outer, LM inner. The cubic
    !! interpolation weights depend only on |G| and are computed once per G.
    !! The LM sum per G is accumulated in the original order (bit-identical).
    INTEGER, INTENT(IN) :: ig1, ig2
    INTEGER  :: ig, lm, i0
    REAL(DP) :: qm, px, ux, vx, wx, uvx, pwx, work, acc(2)
    DO ig = ig1, ig2
       qm = qmod (ig) * dqi
       px = qm - INT(qm)
       ux = 1.0_DP - px
       vx = 2.0_DP - px
       wx = 3.0_DP - px
       i0 = INT(qm) + 1
       uvx = ux * vx * sixth
       pwx = px * wx * 0.5_DP
       acc(1) = 0.0_DP
       acc(2) = 0.0_DP
       DO lm = 1, nlm
          work = tab_qrad(i0  ,ijv,l_v(lm),np) * uvx * wx + &
                 tab_qrad(i0+1,ijv,l_v(lm),np) * pwx * vx - &
                 tab_qrad(i0+2,ijv,l_v(lm),np) * pwx * ux + &
                 tab_qrad(i0+3,ijv,l_v(lm),np) * px * uvx
          acc(ind_v(lm)) = acc(ind_v(lm)) + sig_v(lm) * ylmk0(ig,lp_v(lm)) * work
       ENDDO
       qg(1,ig) = acc(1)
       qg(2,ig) = acc(2)
    ENDDO
  END SUBROUTINE qvan2_kernel
#endif
  !
  SUBROUTINE qvan2_lsig( lp, l, sig, ind )
    !! Angular momentum (l+1), sign of (-i)^L and real/imag slot for combined index lp
    INTEGER, INTENT(IN) :: lp
    INTEGER, INTENT(OUT) :: l, ind
    REAL(DP), INTENT(OUT) :: sig
    IF ( lp == 1 ) THEN
       l = 1 ; sig = 1.0_DP ; ind = 1
    ELSEIF ( lp <=  4 ) THEN
       l = 2 ; sig =-1.0_DP ; ind = 2
    ELSEIF ( lp <=  9 ) THEN
       l = 3 ; sig =-1.0_DP ; ind = 1
    ELSEIF ( lp <= 16 ) THEN
       l = 4 ; sig = 1.0_DP ; ind = 2
    ELSEIF ( lp <= 25 ) THEN
       l = 5 ; sig = 1.0_DP ; ind = 1
    ELSEIF ( lp <= 36 ) THEN
       l = 6 ; sig =-1.0_DP ; ind = 2
    ELSE
       l = 7 ; sig =-1.0_DP ; ind = 1
    ENDIF
  END SUBROUTINE qvan2_lsig
  !
END SUBROUTINE qvan2

!----------------------------------------------------------------------
subroutine compute_qqr ( tpiba, q, omega, qq_nt )
  !----------------------------------------------------------------------
  !
  !   The qq are the G=0 components of Q or, in general, Q(q)
  !
  USE upf_kinds,    ONLY : DP
  USE uspp_param,   ONLY : upf, lmaxq, nh, nhm, nsp
  !
  real(DP), intent(in)  :: q(3)
  !! input wave-vector, can be q = 0 (scf) or q != 0 (Berry's phase)
  real(DP), intent(in)  :: omega
  !! cell size
  real(DP), intent(in) :: tpiba
  !! 2pi/a (can be set to anything if q=0)
  real(DP), intent(out) :: qq_nt(nhm,nhm,nsp)
  !! ouput Q(q)
  real(DP) :: ylmk0 (lmaxq*lmaxq,1)
  real(DP) :: q_(3,1)
  real(DP) :: qmod(1)
  real(DP) :: qgm(2,1)
  !! local variables - the presence of an index "1" prevents mismatch errors
  !! when calling ylmr2 and qvan2 (both expect arrays with that dimension)
  !
  q_(:,1) = q(:)
  qmod(1) = q(1)**2 + q(2)**2 + q(1)**2
  call ylmr2 (lmaxq * lmaxq, 1, q_, qmod, ylmk0)
  !
  qmod(1) = sqrt ( qmod(1) ) * tpiba
  do nt = 1, nsp
     if ( upf(nt)%tvanp ) then
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (1, ih, jh, nt, qmod, qgm, ylmk0)
              qq_nt(ih,jh,nt) = omega * qgm (1,1)
              qq_nt(jh,ih,nt) = omega * qgm (1,1)
           enddo
        enddo
     endif
  enddo
  !
end subroutine compute_qqr
!
!----------------------------------------------------------------------
subroutine compute_qqc ( tpiba, q, omega, qq_nt )
  !----------------------------------------------------------------------
  !
  !   Compute complex Q(q) for q /= 0 
  !
  USE upf_kinds,    ONLY : DP
  USE uspp_param,   ONLY : upf, lmaxq, nh, nhm, nsp
  !
  real(DP), intent(in)  :: q(3)
  !! input wave-vector, can be q = 0 (scf) or q != 0 (Berry's phase)
  real(DP), intent(in)  :: omega
  !! cell size
  real(DP), intent(in) :: tpiba
  !! 2pi/a (can be set to anything if q=0)
  complex(DP), intent(out) :: qq_nt(nhm,nhm,nsp)
  !! output Q(q)
  real(DP) :: ylmk0 (lmaxq*lmaxq,1)
  real(DP) :: q_(3,1)
  real(DP) :: qmod(1)
  real(DP) :: qgm(2,1)
  !! local variables - the presence of an index "1" prevents mismatch errors
  !! when calling ylmr2 and qvan2 (both expect arrays with that dimension)
  !
  q_(:,1) = q(:)
  qmod(1) = q(1)**2 + q(2)**2 + q(1)**2
  call ylmr2 (lmaxq * lmaxq, 1, q_, qmod, ylmk0)
  !
  qmod(1) = sqrt ( qmod(1) ) * tpiba
  do nt = 1, nsp
     if ( upf(nt)%tvanp ) then
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (1, ih, jh, nt, qmod, qgm, ylmk0)
              qq_nt(ih,jh,nt) = omega * CMPLX(qgm (1,1),qgm(2,1),kind=dp )
              qq_nt(jh,ih,nt) = omega * CMPLX(qgm (1,1),qgm(2,1),kind=dp )
           enddo
        enddo
     endif
  enddo
  !
end subroutine compute_qqc
