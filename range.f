! version %I% of %G%
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
  implicit none
  integer n1, n2, nprocs, irank, ista, iend, iwork1, iwork2
  iwork1=(n2-n1+1) / nprocs
  iwork2=MOD(n2-n1+1, nprocs)
  ista = irank*iwork1 + n1 + MIN(irank,iwork2)
  iend = ista + iwork1 - 1
  IF(iwork2>irank) iend = iend + 1
end subroutine para_range
