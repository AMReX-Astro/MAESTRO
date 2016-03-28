!     You can made CUDA Fortran programs like this with 
!     $ pgf90 -Mcuda deadbeef.F
      Use CudaFor
      Implicit None
C
      Integer*8, Device, Allocatable :: A(:)
      Integer*8 N,FreeMem,TotGlbMem
      Integer*8 iNaN,iCero,iSet
      Integer*4 iStat
      Logical IsKrnOK
      Data iNaN/z'DEADBEEFDEADBEEF'/
      Data iCero/0/
c--------------------------------------------------------------
c
      iStat = cudaMemGetInfo(FreeMem, TotGlbMem)
      N = FreeMem/8 - 100000
      Allocate(A(N))
      iSet = iCero
      iSet = iNan
      iStat = cudaMemSet(A,iSet,N)
      If(IsKrnOK(' Problem with cudaMemSet')) then
        If(iSet.eq.iCero) then
          write(6,*) ' Device Memory is "0000000"'
        Else
          write(6,*) ' Device Memory is "DEADBEEF"'
        EndIf
      EndIf
      End

c
c==============================
*Deck IsKrnOK
      Logical Function IsKrnOK(Txt)
c
      Use CUDAfor
c
      Implicit None
C
C     Checks for (A)Syncronous errors from CUDA kernel launches
C     .True.    Everything is OK
c     .False.   At least one error
C     Resets the Error variable
C
      Character*(*) Txt
      Integer*4 iErrS,iErrA
c
      IsKrnOK = .True.
      iErrS = cudaGetLastError()
      iErrA = cudaDeviceSynchronize()
      If(iErrS.ne.cudaSuccess)  then
        IsKrnOK = .False.
        write(6,'(3(a,1x))') ' Sync Kernel Error',Txt,
     $            cudaGetErrorString(iErrS)
      EndIf
      If(iErrA.ne.cudaSuccess)  then
        IsKrnOK = .False.
        write(6,'(3(a,1x))') ' ASync Kernel Error',Txt,
     $            cudaGetErrorString(cudaGetLastError())
      EndIf

      Return
      End

