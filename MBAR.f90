module MBAR_m
  use precision_m
  use constant_m
  use io_m
  implicit none
  
  contains
    subroutine MBAR(maxIteration,criterion,nSimulations,totalNumSnapshots,reducedEnergies, &
         & nSnapshotsInSimulation,freeEnergies,weights,covMat,solver)
      integer(kind=4), intent(in) :: maxIteration 
      real(kind=fp_kind), intent(in) :: criterion
      integer(kind=4), intent(in) :: nSimulations
      integer(kind=4), intent(in) :: totalNumSnapshots
      real(kind=fp_kind), intent(in) :: reducedEnergies(totalNumSnapshots,nSimulations)
      integer(kind=4), intent(in) :: nSnapshotsInSimulation(nSimulations)
      real(kind=fp_kind), intent(out) :: freeEnergies(nSimulations)
      real(kind=fp_kind), intent(out) :: weights(totalNumSnapshots,nSimulations)
      real(kind=fp_kind), intent(out) :: covMat(nSimulations,nSimulations)
      character(len=3), intent(in), optional :: solver

      real(kind=fp_kind) :: minReducedEnergy
      real(kind=fp_kind), allocatable :: shiftedReducedEnergies(:,:)

      real(kind=fp_kind), allocatable ::oldFreeEnergies(:)
      real(kind=fp_kind), allocatable :: oldExpFreeEnergies(:)
      real(kind=fp_kind), allocatable :: expShiftedReducedEnergies(:,:)
     
      character(len=3) :: solv = 'NRM'

      real(kind=fp_kind), allocatable :: numerator(:,:), denominator(:)
      integer(kind=4) :: Iteration
      integer(kind=4) :: IndexW, IndexS
      integer(kind=4) :: JndexW, JndexS
      integer(kind=4) :: KndexW

      real(kind=fp_kind) :: freeEnergyRmsd
      real(kind=fp_kind), allocatable :: jacobian(:,:), residual(:)
      real(kind=fp_kind) :: rmsd

      real(kind=fp_kind), allocatable :: jacobianInv(:,:)
      real(kind=fp_kind), allocatable :: deltaFreeEnergies(:)

      real(kind=fp_kind) :: freeEnergyShift

      real(kind=fp_kind) :: maxRelativeDelta

      logical :: file_exists

      if(idebug == 1) write(6,*) 'Entering MBAR'

      if(present(solver))solv = solver

      allocate(oldFreeEnergies(nSimulations))
      allocate(shiftedReducedEnergies(totalNumSnapshots,nSimulations))
      allocate(oldExpFreeEnergies(nSimulations))
      allocate(expShiftedReducedEnergies(totalNumSnapshots,nSimulations))
 
      minReducedEnergy = minval(reducedEnergies)   
      shiftedReducedEnergies = reducedEnergies - minReducedEnergy
      expShiftedReducedEnergies = exp(-shiftedReducedEnergies)
 
      forall(IndexW = 1:nSimulations) freeEnergies(IndexW) = minval(shiftedReducedEnergies(:,IndexW))

      allocate(numerator(totalNumSnapshots,nSimulations))
      allocate(denominator(totalNumSnapshots))
      if(solv == 'NRM')then
        allocate(jacobian(nSimulations,nSimulations))
        allocate(jacobianInv(nSimulations,nSimulations))
        allocate(residual(nSimulations))
        allocate(deltaFreeEnergies(nSimulations))
 
        write(6,'(1X,A)')'Solving non-linear equations using Newton-Raphson method:'
        Iteration = 0
        do while (Iteration < maxIteration)
          oldFreeEnergies = freeEnergies
          oldExpFreeEnergies = exp(-oldFreeEnergies)
 
          forall (JndexS = 1: totalNumSnapshots, IndexW = 1 : nSimulations) numerator(JndexS,IndexW) = &
              &  expShiftedReducedEnergies(JndexS,IndexW) / oldExpFreeEnergies(IndexW)
          forall (JndexS = 1: totalNumSnapshots) denominator(JndexS) = &
              &  sum(nSnapshotsInSimulation(:)*numerator(JndexS,:))

          forall(IndexW = 1 : nSimulations, JndexS = 1 : totalNumSnapshots) &
             & weights(JndexS,IndexW) = numerator(JndexS,IndexW)/denominator(JndexS)
         
          forall (IndexW=1:nSimulations) residual(IndexW) = 1.0d0 - sum(weights(:,IndexW))
  
!          if(Iteration==0)write(6,'(1X,A,I6,A,ES12.5)')'Iteration: ', Iteration, & 
!               & ', sum of |residual|: ', sum(abs(residual))
 
          Iteration = Iteration + 1
 
          forall(IndexW = 1 : nSimulations, JndexW = 1 : nSimulations) jacobian(IndexW, JndexW) = &
              & sum(numerator(:,IndexW) * nSnapshotsInSimulation(JndexW) * &
              &     numerator(:,JndexW)/denominator(:)**2) 
          forall(IndexW = 1 : nSimulations) jacobian(IndexW,IndexW) = &
              & jacobian(IndexW,IndexW) - sum(numerator(:,IndexW)/denominator(:))
 
          call squareMatrixInv(nSimulations,jacobian,jacobianInv)
 
          forall(IndexW = 1 : nSimulations) deltaFreeEnergies(IndexW) = &
             & sum(jacobianInv(IndexW,:)*residual(:))

          freeEnergies = oldFreeEnergies - deltaFreeEnergies
          freeEnergies = freeEnergies - freeEnergies(1)
          freeEnergyRmsd = rmsd(nSimulations,oldFreeEnergies,freeEnergies)
          write(6,'(1X,A,I6,A,ES12.5,A,ES12.5)')'Iteration: ', Iteration, ', sum of |residual|: ', &
                & sum(abs(residual)), ', RMSD of free energy: ', freeEnergyRmsd
          if(maxval(abs((freeEnergies(2:)-oldFreeEnergies(2:))/oldFreeEnergies(2:)))<criterion)exit
        end do
      else if(solv=='DI0' .or. solv=='DIR') then
        write(6,'(1X,A)')'Solving non-linear equations using Direct-Iteration:'
        Iteration = 0
        if(solv=='DIR')then
          INQUIRE(FILE="fort.66", EXIST=file_exists)
          if(file_exists)read(66)freeEnergies
        end if
        do while (Iteration < maxIteration)
          Iteration = Iteration + 1
          oldFreeEnergies = freeEnergies
          oldExpFreeEnergies = exp(-oldFreeEnergies)
 
          forall (JndexS = 1: totalNumSnapshots, IndexW = 1 : nSimulations) numerator(JndexS,IndexW) = &
              &  expShiftedReducedEnergies(JndexS,IndexW)
          forall (JndexS = 1: totalNumSnapshots) denominator(JndexS) = &
              &  sum(nSnapshotsInSimulation(:)*numerator(JndexS,:)/oldExpFreeEnergies(:))

          forall(IndexW = 1 : nSimulations, JndexS = 1 : totalNumSnapshots) &
             & weights(JndexS,IndexW) = numerator(JndexS,IndexW)/denominator(JndexS)

          forall(IndexW = 1 : nSimulations) freeEnergies(IndexW) = sum(weights(:,IndexW))

          freeEnergies = -log(freeEnergies)
          freeEnergies = freeEnergies - freeEnergies(1)
          freeEnergyRmsd = rmsd(nSimulations,oldFreeEnergies,freeEnergies)
          maxRelativeDelta = maxval(abs((freeEnergies(2:)-oldFreeEnergies(2:))/oldFreeEnergies(2:)))
          write(6,'(1X,A,I6,A,ES12.5,A,ES12.5)')'Iteration: ', Iteration, ', RMSD of reduced free energy: ', &
             & freeEnergyRmsd, ' Max relative delta: ', maxRelativeDelta
          if(iteration>1.and.maxRelativeDelta<criterion)then
            write(6,'(A)')'Convergence criterion met'
            exit
          end if
        end do
      end if
      if(Iteration >= maxIteration)then
         write(6,*)'Maximum iteration arrived'
         write(66)freeEnergies
         stop
      end if

      forall (JndexS = 1 : totalNumSnapshots, IndexW = 1 : nSimulations) &
          & numerator(JndexS,IndexW) = expShiftedReducedEnergies(JndexS,IndexW)/exp(-freeEnergies(IndexW))
      forall (JndexS = 1 : totalNumSnapshots) denominator(JndexS) = &
          & sum(nSnapshotsInSimulation(:) * numerator(JndexS,:))
      forall(IndexW = 1 : nSimulations, JndexS = 1 : totalNumSnapshots) &
          & weights(JndexS,IndexW) = numerator(JndexS,IndexW)/denominator(JndexS)
     
      call ComputCovMatFromWeights(totalNumSnapshots,nSimulations,nSnapshotsInSimulation,weights,covMat)

      if(allocated(oldFreeEnergies))deallocate(oldFreeEnergies)
      if(allocated(shiftedReducedEnergies))deallocate(shiftedReducedEnergies)
      if(allocated(jacobian))deallocate(jacobian)
      if(allocated(jacobianInv))deallocate(jacobianInv)
      if(allocated(residual))deallocate(residual)
      if(allocated(oldExpFreeEnergies))deallocate(oldExpFreeEnergies)
      if(allocated(expShiftedReducedEnergies))deallocate(expShiftedReducedEnergies)
      if(allocated(deltaFreeEnergies))deallocate(deltaFreeEnergies)
      if(allocated(numerator))deallocate(numerator)
      if(allocated(denominator))deallocate(denominator)
    end subroutine MBAR
  
    subroutine MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation, &
         & freeEnergies,targetReducedEnergies,weights,freeEnergy,freeEnergySD)
      implicit none
      integer(kind=4), intent(in) :: nSimulations
      integer(kind=4), intent(in) :: totalNumSnapshots
      real(kind=fp_kind), intent(in) :: reducedEnergies(totalNumSnapshots,nSimulations)
      integer(kind=4), intent(in) :: nSnapshotsInSimulation(nSimulations)
      real(kind=fp_kind), intent(in) :: freeEnergies(nSimulations)
      real(kind=fp_kind), intent(in) :: targetReducedEnergies(totalNumSnapshots)
      real(kind=fp_kind), intent(in out) :: weights(totalNumSnapshots,nSimulations+1)
      real(kind=fp_kind), intent(out) :: freeEnergy
      real(kind=fp_kind), intent(out) :: freeEnergySD

      integer(kind=4) :: extNSnapshotsInSimulation(nSimulations+1)
      real(kind=fp_kind) :: extCovFreeEnergies(nSimulations+1,nSimulations+1)
      integer(kind=4) :: IndexS, IndexW
      integer(kind=4) :: JndexS, JndexW
      integer(kind=4) :: KndexS, KndexW

      if(idebug == 1) write(6,*) 'Entering MBAR_weight'
      forall (JndexS = 1: totalNumSnapshots) weights(JndexS,nSimulations+1) = &
          & exp(-targetReducedEnergies(JndexS))/ &
          & sum(nSnapshotsInSimulation(:)*exp(freeEnergies(:)-reducedEnergies(JndexS,:)))
      freeenergy = -log(sum(weights(:,nSimulations+1)))

      extNSnapshotsInSimulation = 0
      extNSnapshotsInSimulation(1:nSimulations) = nSnapshotsInSimulation(1:nSimulations)
      weights(:,nSimulations+1) = weights(:,nSimulations+1)/sum(weights(:,nSimulations+1))
      call ComputCovMatFromWeights(totalNumSnapshots,nSimulations+1,extNSnapshotsInSimulation,weights,extCovFreeEnergies)
      freeEnergySD = sqrt(extCovFreeEnergies(nSimulations+1,nSimulations+1) - &
                     &  2*extCovFreeEnergies(nSimulations+1,1) + extCovFreeEnergies(1,1) )
    end subroutine MBAR_weight

    subroutine squareMatrixInv(m,a,ainv)
      implicit none
      integer(kind=4), intent(in) :: m
      real(kind=fp_kind), intent(in) :: a(m,m)
      real(kind=fp_kind), intent(out) :: ainv(m,m)

      real(kind=fp_kind), allocatable :: acopy(:,:)
      real(kind=fp_kind), allocatable :: sigma(:), u(:,:), vt(:,:)
      real(kind=fp_kind), allocatable :: work(:)
      integer(kind=4) :: lwork, info

      real(kind=fp_kind), allocatable :: vs(:,:)

      integer(kind=4) :: i, j, k

      lwork = 10*m
      allocate(acopy(m,m))
      allocate(work(lwork))
      allocate(sigma(m))
      allocate(u(m,m))
      allocate(vt(m,m))
      allocate(vs(m,m))
      acopy=a
      call dgesvd('A','A',m,m,acopy,m,sigma,u,m,vt,m,work,lwork,info)
      if(info<0)then
        write(6,*)'dgesvd failed with info=',info
        stop
      else if(info>0) then
        write(6,*)'dgesvd failed with info=',info
        stop
      else
        do i = 1, m
          do j =1, m
            vs(i,j)=vt(j,i)/sigma(j)
          end do
        end do
        ainv=0.d0
        do i = 1, m
          do j = 1, m
            do k = 1, m
              ainv(i,j)=ainv(i,j)+vs(i,k)*u(j,k)
            end do
          end do
        end do
      end if
      deallocate(acopy)
      deallocate(work)
      deallocate(sigma)
      deallocate(u)
      deallocate(vt)
      deallocate(vs)
    end subroutine squareMatrixInv

    subroutine squareMatrixInv2(m,a,ainv)
      implicit none
      integer(kind=4), intent(in) :: m
      real(kind=fp_kind), intent(in) :: a(m,m)
      real(kind=fp_kind), intent(out) :: ainv(m,m)

      real(kind=fp_kind), allocatable :: awork(:,:)
      real(kind=fp_kind), allocatable :: work(:)
      real(kind=fp_kind), allocatable :: ipiv(:)
      integer(kind=4) :: lwork
      integer(kind=4) :: info

      integer(kind=4) :: i, j, k

      allocate(ipiv(m))    
      ainv = a
      call dgetrf( m, m, ainv, m, ipiv, info )
      if(info < 0) then
        write(6,'('' squareMatrixInv2: LU decomp has an illegal value'')')
        stop
      elseif(info > 0) then
        write(6,'('' squareMatrixInv2: LU decomp gave almost-singular U'')')
        stop
      endif
      
      lwork = 10*m
      allocate(work(lwork))
      call dgetri( m, ainv, m, ipiv, work, lwork, info )
      if(info < 0) then
        write(6,*) 'parameter ',-info, ' has an illegal value' 
        stop
      elseif(info > 0) then
        write(6,*) 'diagonal element ', info, ' of the factor U is zero'
        stop
      endif
      deallocate(ipiv)
      deallocate(work)
    end subroutine squareMatrixInv2

    subroutine ComputCovMatFromWeights(N,K,nvec,weights,covMat)
      implicit none
      integer(kind=4), intent(in) :: N, K
      integer(kind=4), intent(in) :: nvec(K)
      real(kind=fp_kind), intent(in) :: weights(N,K)
      real(kind=fp_kind), intent(out) :: covMat(K,K)
      
      real(kind=fp_kind) :: nmat(K,K)

      real(kind=fp_kind), allocatable :: w(:,:)
      real(kind=fp_kind), allocatable :: sigma(:), u(:,:), vt(:,:)
      real(kind=fp_kind), allocatable :: work(:)
      integer(kind=4) :: lwork, info

      real(kind=fp_kind), allocatable :: sigmaKK(:,:)
      real(kind=fp_kind), allocatable :: svtnvs(:,:)
      real(kind=fp_kind), allocatable :: tmp(:,:)
      real(kind=fp_kind), allocatable :: v(:,:)

      real(kind=fp_kind), allocatable :: sigma2(:), u2(:,:), vt2(:,:)
      real(kind=fp_kind), allocatable :: v2(:,:), sigmaKK2(:,:), ut2(:,:)
      real(kind=fp_kind), allocatable :: innerinverse(:,:)

      integer(kind=4) :: indexN, indexK
      integer(kind=4) :: jndexN, jndexK
      integer(kind=4) :: kndexN, kndexK

      lwork = 10*N
      allocate(w(N,K))
      allocate(work(lwork))
     
      w = weights

      allocate(sigma(N))
      allocate(u(N,K))
      allocate(vt(K,K))

      call dgesvd('S','S',N,K,w,N,sigma,u,N,vt,K,work,lwork,info)
      if(allocated(w))deallocate(w)

      if(info < 0)then
        write(6,*)'dgesvd of weight matrix failed with info=',info
        stop
      else if(info > 0) then
        write(6,*)'dgesvd of weight matrix failed with info=',info
        stop
      else
        allocate(sigmaKK(K,K))
        allocate(svtnvs(K,K))
        allocate(v(K,K))
        allocate(tmp(K,K))
        v=transpose(vt)

        sigmaKK = 0.d0
        forall (indexK = 1:K) sigmaKK(indexK,indexK) = sigma(indexK)

        nmat = 0.d0
        forall (indexK = 1:K) nmat(indexK,indexK) = dble(nvec(indexK))

        svtnvs = matmul(sigmaKK,vt)
        svtnvs = matmul(svtnvs,nmat)
        svtnvs = matmul(svtnvs,v)
        svtnvs = -matmul(svtnvs,sigmaKK)

        forall (indexK = 1:K)  svtnvs(indexK,indexK) = 1.0d0 + svtnvs(indexK,indexK)
        if(allocated(work))deallocate(work)
        lwork = 10*K
        allocate(work(lwork))
        allocate(sigma2(K))
        allocate(u2(K,K))
        allocate(vt2(K,K))
        call dgesvd('A','A',K,K,svtnvs,K,sigma2,u2,K,vt2,K,work,lwork,info)
        if(info<0)then
          write(6,*)'dgesvd of svtnvs matrix failed with info=',info
          stop
        else if(info>0) then
          write(6,*)'dgesvd of svtnvs matrix failed with info=',info
          stop
        else
          allocate(innerinverse(K,K))
          allocate(v2(K,K))
          v2 = transpose(vt2)
          allocate(sigmaKK2(K,K))
          sigmaKK2 = 0.d0
          do indexK = 1, K
            if(abs(sigma2(indexK))>1.E-10)sigmaKK2(indexK,indexK)=1.0d0/sigma2(indexK)
          end do
          allocate(ut2(K,K))
          ut2 = transpose(u2)

          innerinverse = matmul(v2,sigmaKK2)
          innerinverse = matmul(innerinverse,ut2)

          if(allocated(v2))deallocate(v2)
          if(allocated(sigmaKK2))deallocate(sigmaKK2)
          if(allocated(ut2))deallocate(ut2)
          if(allocated(sigma2))deallocate(sigma2)
          if(allocated(u2))deallocate(u2)
          if(allocated(vt2))deallocate(vt2)

          covMat = matmul(sigmaKK,innerinverse)
          covMat = matmul(covMat,sigmaKK)
          covMat = matmul(v,covMat)
          covMat = matmul(covMat,vt)

          if(allocated(innerinverse))deallocate(innerinverse)
        end if
        if(allocated(sigmaKK))deallocate(sigmaKK)
        if(allocated(svtnvs))deallocate(svtnvs)
        if(allocated(v))deallocate(v)
        if(allocated(tmp))deallocate(tmp)
      end if
      if(allocated(work))deallocate(work)
      if(allocated(sigma))deallocate(sigma)
      if(allocated(u))deallocate(u)
      if(allocated(vt))deallocate(vt)
    end subroutine ComputCovMatFromWeights
end module MBAR_m
