module MBAR_m
  use precision_m
  use constant_m
  use io_m
  use snapshot_m
  use simulation_m
  use reducedHamiltonian_m
  implicit none
  
  contains
    subroutine MBAR(maxIteration,criterion,nSimulations,TotalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies)
      integer(kind=4), intent(in) :: maxIteration 
      real(kind=fp_kind), intent(in) :: criterion
      integer(kind=4), intent(in) :: nSimulations
      integer(kind=4), intent(in) :: TotalNumSnapshots
      real(kind=fp_kind), intent(in) :: reducedEnergies(TotalNumSnapshots,nSimulations)
      integer(kind=4), intent(in) :: nSnapshotsInSimulation(nSimulations)
      real(kind=fp_kind), intent(out) :: freeEnergies(nSimulations)

      real(kind=fp_kind) :: minReducedEnergy
      real(kind=fp_kind), allocatable :: oldFreeEnergies(:)
      real(kind=fp_kind), allocatable :: shiftedReducedEnergies(:,:)

      real(kind=fp_kind) :: numerator, denominator
      integer(kind=4) :: Iteration
      integer(kind=4) :: IndexW, IndexS
      integer(kind=4) :: JndexW, JndexS

      real(kind=fp_kind) :: freeEnergyRmsd
      real(kind=fp_kind) :: rmsd
 
      if(idebug == 1) write(*,*) 'Entering MBAR'

      allocate(oldFreeEnergies(nSimulations))
      allocate(shiftedReducedEnergies(TotalNumSnapshots,nSimulations))
  
      freeEnergies = 0.d0
      minReducedEnergy = minval(reducedEnergies)   
      shiftedReducedEnergies = reducedEnergies - minReducedEnergy

      Iteration = 0
      do while (Iteration <maxIteration)
        Iteration = Iteration + 1
        oldFreeEnergies = freeEnergies - freeEnergies(1)
        do IndexW = 1, nSimulations
          freeEnergies(IndexW) = 0.d0
          do JndexS = 1, TotalNumSnapshots
            numerator = exp(-shiftedReducedEnergies(JndexS,IndexW))
            denominator = 0.d0
            do JndexW = 1, nSimulations
              denominator = denominator + nSnapshotsInSimulation(JndexW)*exp(-(shiftedReducedEnergies(JndexS,JndexW)-oldFreeEnergies(JndexW)))
            end do
            freeEnergies(IndexW) = freeEnergies(IndexW) + &
               & numerator/denominator
          end do

          freeEnergies(IndexW) = log(freeEnergies(IndexW))
        end do
        freeEnergies = freeEnergies - freeEnergies(1)
        freeEnergyRmsd = rmsd(nSimulations,oldFreeEnergies,freeEnergies)
        write(*,*)'RMSD of free energy: ',freeEnergyRmsd
        if(freeEnergyRmsd<criterion)exit
      end do
      if(Iteration >= maxIteration)then
         write(*,*)'Maximum iteration arrived'
         stop
      else 
        do IndexW = 1, nSimulations
!          simulatedReducedHamiltonian(IndexW)%freeEnergy=freeEnergies(IndexW)
        end do
      end if
      deallocate(oldFreeEnergies)
      deallocate(shiftedReducedEnergies)
    end subroutine MBAR
  
    subroutine MBAR_weight(nSimulations,totalNumSnapshots,reducedEnergies,nSnapshotsInSimulation,freeEnergies,targetReducedEnergies,weights)
      implicit none
      integer(kind=4), intent(in) :: nSimulations
      integer(kind=4), intent(in) :: TotalNumSnapshots
      real(kind=fp_kind), intent(in) :: reducedEnergies(TotalNumSnapshots,nSimulations)
      integer(kind=4), intent(in) :: nSnapshotsInSimulation(nSimulations)
      real(kind=fp_kind), intent(in) :: freeEnergies(nSimulations)
      real(kind=fp_kind), intent(in) :: targetReducedEnergies(TotalNumSnapshots)
      real(kind=fp_kind), intent(out) :: weights(TotalNumSnapshots)


      real(kind=fp_kind) :: numerator, denominator
      integer(kind=4) :: IndexS, IndexW
      integer(kind=4) :: JndexS, JndexW

      if(idebug == 1) write(*,*) 'Entering MBAR_weight'
      do JndexS = 1, TotalNumSnapshots
        numerator = exp(-targetReducedEnergies(JndexS))
        denominator = 0.d0
        do JndexW = 1, nSimulations
          denominator = denominator + nSnapshotsInSimulation(JndexW)*exp(-(reducedEnergies(JndexS,JndexW)-freeEnergies(JndexW)))
        end do
        weights(JndexS) = numerator/denominator
      end do
      weights = weights/sum(weights)
    end subroutine MBAR_weight
end module MBAR_m
