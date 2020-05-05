module MBAR_m
  use precision_m
  use constant_m
  use io_m
  use snapshot_m
  use simulation_m
  use reducedHamiltonian_m
  implicit none
  
  contains
    subroutine MBAR(maxIteration,criterion,nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies)
      integer(kind=4), intent(in) :: maxIteration 
      real(kind=fp_kind), intent(in) :: criterion
      integer(kind=4), intent(in) :: nSimulations
      integer(kind=4), intent(in) :: totalNumSnapshots
      real(kind=fp_kind), intent(in) :: reducedEnergies(totalNumSnapshots,nSimulations)
      integer(kind=4), intent(in) :: nSnapshotsInSimulation(nSimulations)
      real(kind=fp_kind), intent(out) :: freeEnergies(nSimulations)

      real(kind=fp_kind) :: minReducedEnergy
      real(kind=fp_kind), allocatable :: shiftedReducedEnergies(:,:)

      real(kind=fp_kind), allocatable ::oldFreeEnergies(:)
      real(kind=fp_kind), allocatable :: oldExpFreeEnergies(:)
      real(kind=fp_kind), allocatable :: expShiftedReducedEnergies(:,:)


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

      if(idebug == 1) write(*,*) 'Entering MBAR'

      do IndexW = 1, nSimulations
        write(77,'(10(1X,G12.5))')reducedEnergies(:,IndexW)
      end do
      allocate(oldFreeEnergies(nSimulations))
      allocate(shiftedReducedEnergies(totalNumSnapshots,nSimulations))
      allocate(oldExpFreeEnergies(nSimulations))
      allocate(expShiftedReducedEnergies(totalNumSnapshots,nSimulations))

      minReducedEnergy = minval(reducedEnergies)   
      shiftedReducedEnergies = reducedEnergies - minReducedEnergy
      do IndexW = 1, nSimulations
        write(88,'(10(1X,G12.5))')shiftedReducedEnergies(:,IndexW)
      end do
      expShiftedReducedEnergies = exp(-shiftedReducedEnergies)

      do IndexW = 1, nSimulations
        freeEnergies(IndexW) = minval(shiftedReducedEnergies(:,IndexW))
        write(99,*)minval(shiftedReducedEnergies(:,IndexW))
      end do

      allocate(jacobian(nSimulations,nSimulations))
      allocate(jacobianInv(nSimulations,nSimulations))
      allocate(residual(nSimulations))
      allocate(deltaFreeEnergies(nSimulations))
      allocate(numerator(totalNumSnapshots,nSimulations))
      allocate(denominator(totalNumSnapshots))
      Iteration = 0
      do while (Iteration < maxIteration)
        oldFreeEnergies = freeEnergies
        oldExpFreeEnergies = exp(-oldFreeEnergies)

        denominator = 0.d0
        do JndexS = 1, totalNumSnapshots
          do IndexW = 1, nSimulations
            numerator(JndexS,IndexW) = expShiftedReducedEnergies(JndexS,IndexW)/oldExpFreeEnergies(IndexW)
          end do
          do KndexW = 1, nSimulations
            denominator(JndexS) = denominator(JndexS) + &
              & nSnapshotsInSimulation(KndexW) * numerator(JndexS,KndexW)
          end do
        end do

        residual = 1.0d0
        do IndexW = 1, nSimulations
          do JndexS = 1, totalNumSnapshots
            residual(IndexW) = residual(IndexW) - numerator(JndexS,IndexW)/denominator(JndexS)
          end do
        end do
 
        if(Iteration==0)write(*,'(A,I6,A,G12.5)')'iteration: ', Iteration, ', sum of |residual|: ', sum(abs(residual))

        Iteration = Iteration + 1

        jacobian = 0.d0
        do IndexW = 1, nSimulations
          do JndexW = 1, nSimulations
            do JndexS = 1, totalNumSnapshots
              jacobian(IndexW, JndexW) = jacobian(IndexW, JndexW) + &
                 & numerator(JndexS,IndexW) * nSnapshotsInSimulation(JndexW) * &
                 & numerator(JndexS,JndexW)/denominator(JndexS)**2
              if(IndexW == JndexW) jacobian(IndexW, JndexW) = jacobian(IndexW,JndexW) - &
                 & numerator(JndexS,IndexW)/denominator(JndexS)
            end do
          end do
        end do

        call squareMatrixInv(nSimulations,jacobian,jacobianInv)

        deltaFreeEnergies = 0.d0
        do IndexW = 1, nSimulations
          do JndexW = 1, nSimulations
            deltaFreeEnergies(IndexW) = deltaFreeEnergies(IndexW) + &
              & jacobianInv(IndexW,JndexW)*residual(JndexW)
          end do
        end do
        freeEnergies = oldFreeEnergies - deltaFreeEnergies
        freeEnergies = freeEnergies - freeEnergies(1)
        freeEnergyRmsd = rmsd(nSimulations,oldFreeEnergies,freeEnergies)
        write(*,'(A,I6,A,G12.5,A,G12.5)')'iteration: ', Iteration,', RMSD of free energy: ', freeEnergyRmsd, ', sum of |residual|: ', sum(abs(residual))
        if(freeEnergyRmsd<criterion)exit
      end do
      if(Iteration >= maxIteration)then
         write(*,*)'Maximum iteration arrived'
         stop
      else 
        write(*,'(1X,G12.5)')FreeEnergies
        do IndexW = 1, nSimulations
!          simulatedReducedHamiltonian(IndexW)%freeEnergy=freeEnergies(IndexW)
        end do
      end if
      deallocate(oldFreeEnergies)
      deallocate(shiftedReducedEnergies)
      deallocate(jacobian)
      deallocate(jacobianInv)
      deallocate(residual)
      deallocate(oldExpFreeEnergies)
      deallocate(expShiftedReducedEnergies)
      deallocate(deltaFreeEnergies)
      deallocate(numerator)
      deallocate(denominator)
    end subroutine MBAR
  
    subroutine MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies,targetReducedEnergies,weights)
      implicit none
      integer(kind=4), intent(in) :: nSimulations
      integer(kind=4), intent(in) :: totalNumSnapshots
      real(kind=fp_kind), intent(in) :: reducedEnergies(totalNumSnapshots,nSimulations)
      integer(kind=4), intent(in) :: nSnapshotsInSimulation(nSimulations)
      real(kind=fp_kind), intent(in) :: freeEnergies(nSimulations)
      real(kind=fp_kind), intent(in) :: targetReducedEnergies(totalNumSnapshots)
      real(kind=fp_kind), intent(out) :: weights(totalNumSnapshots)


      real(kind=fp_kind) :: numerator, denominator
      integer(kind=4) :: IndexS, IndexW
      integer(kind=4) :: JndexS, JndexW

      if(idebug == 1) write(*,*) 'Entering MBAR_weight'
      do JndexS = 1, totalNumSnapshots
        numerator = exp(-targetReducedEnergies(JndexS))
        denominator = 0.d0
        do JndexW = 1, nSimulations
          denominator = denominator + nSnapshotsInSimulation(JndexW)*exp(-(reducedEnergies(JndexS,JndexW)-freeEnergies(JndexW)))
        end do
        weights(JndexS) = numerator/denominator
      end do
      weights = weights/sum(weights)
    end subroutine MBAR_weight

    subroutine squareMatrixInv(m,a,ainv)
      implicit none
      integer(kind=4), intent(in) :: m
      real(kind=fp_kind), intent(in) :: a(m,m)
      real(kind=fp_kind), intent(out) :: ainv(m,m)

      real(kind=fp_kind), allocatable :: sigma(:), u(:,:), vt(:,:)
      real(kind=fp_kind), allocatable :: work(:)
      integer(kind=4) :: lwork, info

      real(kind=fp_kind), allocatable :: vs(:,:)

      integer(kind=4) :: i, j, k

      lwork = 10*m
      allocate(work(lwork))
      allocate(sigma(m))
      allocate(u(m,m))
      allocate(vt(m,m))
      allocate(vs(m,m))
      call dgesvd('A','A',m,m,a,m,sigma,u,m,vt,m,work,lwork,info)
      if(info<0)then
        write(*,*)'dgesvd failed with info=',info
        stop
      else if(info>0) then
        write(*,*)'dgesvd failed with info=',info
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
      if(info .lt. 0) then
        write(6,'('' squareMatrixInv2: LU decomp has an illegal value'')')
        stop
      elseif(info.gt. 0) then
        write(6,'('' squareMatrixInv2: LU decomp gave almost-singular U'')')
        stop
      endif
      
      lwork = 10*m
      allocate(work(lwork))
      call dgetri( m, ainv, m, ipiv, work, lwork, info )
      if(info .lt. 0) then
        write(6,*) 'parameter ',-info, ' has an illegal value' 
        stop
      elseif(info.gt. 0) then
        write(6,*) 'diagonal element ', info, ' of the factor U is zero'
        stop
      endif
      deallocate(ipiv)
      deallocate(work)
    end subroutine squareMatrixInv2
end module MBAR_m
