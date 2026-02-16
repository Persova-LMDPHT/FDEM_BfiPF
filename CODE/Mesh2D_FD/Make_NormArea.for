!=======================================================================
!   GENERAL REMARKS
!
!   This code is freely available under the following conditions:
!
!   1) The code is to be used only for non-commercial purposes.
!   2) No changes and modifications to the code without prior permission of the developer.
!   3) No forwarding the code to a third party without prior permission of the developer.
!
!		FDEM_BfiPF
!   This file contains an auxiliary procedure. 
!
!   Written by Dr. Marina Persova
!   Novosibirsk State Technical University   
!   Novosibirsk, 630073, Russia
!   mpersova@mail.ru
!   Version 2.0 December 10, 2024
!  
!======================================================================! 

      subroutine MakeStructMesh()
      implicit real*8(a-h,o-z)  
	parameter (npar1 = 1000)  	      
      common/Rlayer/Hlay(npar1),SLay(npar1),SLayZ(npar1)
      common/RlayerSN/SLayNorm(npar1),SLayZNorm(npar1)
      common/dim_Rlayer/Klay
      do i=1,Klay+1
      SLayNorm(i)=SLay(i)
      SLayZNorm(i)=SLayZ(i)
      enddo
      return   
      end