module control_m
  use precision_m
  use constant_m
  use io_m 
  implicit none
  public
  integer(kind=4), parameter :: idebug = 0
  integer(Kind=4) :: nSimulations
  real(kind=fp_kind) :: targetBeta=1.d0/(298*kB)
  integer(kind=4) :: iperiodic
  real(kind=fp_kind) :: period=0.d0
  integer(kind=4) :: nbins
  real(kind=fp_kind) :: binmin, binmax
  integer(kind=4) :: iGaussSmooth

  contains
    subroutine controlInitialize()
      implicit none
      write(6,'(1X,A)')'Number of simulations:'
      read*,nSimulations
      write(6,'(1X,A)')'Name of the meta file:'
      read*,metaFile
      write(6,'(1X,A)')'Is the coordinate periodic? 0: No. 1: Yes.'
      read*,iperiodic
      if(iperiodic>0)iperiodic=1
      if(iperiodic==1) then
        write(6,'(1X,A)')'The period is:'
        read*,period
      end if
      write(6,'(1X,A)')'Number of bins for output:'
      read*,nbins
      write(6,'(1X,A)')'Minimum and maximum of the bins:'
      read*,binmin,binmax
      write(6,'(1X,A)')'Target temperature:'
      read*,targetBeta
      targetBeta = 1.d0/ (targetBeta*kB)
      write(6,'(1X,A)')'Name of the file containing the energies of the target Hamiltonian:'
      read*,targetHamiltonianEnergyFile
      write(6,'(1X,A)')'Name of the out file containing the PMF of the target Hamiltonian:'
      read*,targetHamiltonianPmfFile
      write(6,'(1X,A)')'Name of the out file containing the weights of the target Hamiltonian:'
      read*, targetWeightsFile
      write(6,'(1X,A)')'Perform Gaussian smoothing on the data? 0: No. 1: Yes.'
      read*,iGaussSmooth

    end subroutine controlInitialize

    subroutine controlFinalize()
      implicit none
      
    end subroutine controlFinalize
end module control_m
