module bin_m
  use precision_m
  implicit none
  private
  integer(kind=4), public :: nbins
  real(kind=fp_kind), public :: binmin, binmax, binwidth
  type, public :: bin_t
    real(kind=fp_kind) :: bincenter
    real(kind=fp_kind) :: binwidth
    integer(kind=4) :: nSnapshotsInBin
    real(kind=fp_kind) :: sumOfWeightsInBin
    real(kind=fp_kind) :: reweightingEntropy
    real(kind=fp_kind) :: pmf
    real(kind=fp_kind)  :: pmfSE
    contains

  end type bin_t
  type (bin_t), allocatable, public :: bins(:)
  public :: deleteBinInfo
  public :: initbins
contains
    subroutine initbins()
      implicit none
      integer(kind=4) :: IndexB
      binwidth = (binmax - binmin)/nbins
      allocate(bins(nbins))
      forall(IndexB = 1:nbins) 
         bins(IndexB)%bincenter = binmin + (IndexB-0.5)*binwidth
         bins(IndexB)%binwidth = binwidth
         bins(IndexB)%nSnapshotsInBin = 0
         bins(IndexB)%sumOfWeightsInBin = 0.d0
         bins(IndexB)%reweightingEntropy = 0.d0
         bins(IndexB)%pmf = 0.d0
         bins(IndexB)%pmfSE = 0.d0
      end forall
    end subroutine initbins

    subroutine deleteBinInfo()
      if(allocated(bins))deallocate(bins)
    end subroutine deleteBinInfo


end module bin_m
