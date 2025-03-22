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
!   This file contains procedures for correction of the model 
!   in case the layer sizes are larger than the computational domain size. 
!
!   Written by Dr. Marina Persova
!   Novosibirsk State Technical University   
!   Novosibirsk, 630073, Russia
!   mpersova@mail.ru
!   Version 2.0 December 10, 2024
!  
!======================================================================! 
!
! correction of the model in case the layer sizes are larger than the computational domain size 
      subroutine ModelCorrection()              
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)   	    
      common/dim_Rlayer/Klay
      common/dim_RlayerM/KlayM
      common/Rlayer/Hlay(npar1),SLay(npar1),SLayZ(npar1)
      common/RlayerM/HlayM(npar1),SLayM(npar1),SLayZM(npar1)
      common/RlayerSN/SLayNorm(npar1),SLayZNorm(npar1)      
      common/RlayerMSN/SLayNormM(npar1),SLayZNormM(npar1)      
	common/DprLay/dpr(npar1)  
	common/DprLayM/dprM(npar1)  
      common/settings2d/bak2d,rmin2d,zmin2d,
     *ukoef1_2d,ukoef2_2d,ukoef3_2d
      common/bak2d/ZB1_2D,ZB2_2D
      common/eps_set2d/eps
!
      ZB1_2D=-bak2d-eps
      if (ZB1_2D.gt.Hlay(1)) ZB1_2D=Hlay(1)-100d0
      ZB2_2D=bak2d+eps
!     
! exclusion from the model of layers located outside the computational domain    
      call CorrectLay2D(KLay,Hlay,SLay,SLayZ,SLayNorm,SLayZNorm,dpr,
     *-bak2d-eps,bak2d+eps) 
      ZB1=-baks
      ZB2=baks     
      KlayM=Klay
      do i=1,KlayM-1
      HlayM(i)=Hlay(i)
      SLayM(i)=SLay(i)
      SLayZM(i)=SLayZ(i)
      dprM(i)=dpr(i)
      SLayNormM(i)=SLayNorm(i)     
      SLayZNormM(i)=SLayZNorm(i)          
      enddo 
      SLayM(KlayM)=SLay(KlayM)
      SLayZM(KlayM)=SLayZ(KlayM)
      dprM(KlayM)=dpr(KlayM)          
      SLayNormM(KlayM)=SLayNorm(KlayM)     
      SLayZNormM(KlayM)=SLayZNorm(KlayM)          
! exclusion from the model of layers located outside the computational domain    
      call CorrectLay2D(KLayM,HlayM,SLayM,SLayZM,SLayNormM,SLayZNormM,
     *dprM,-baks-eps,baks+eps)                        
      return
      end
!
!
! exclusion from the model of layers located outside the computational domain    
      subroutine CorrectLay2D(KLay,Hlay,SLay,SLayZ,
     *SLayNorm,SLayZNorm,dpr,ZNE,ZVE) 
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      dimension Hlayt(npar1),Slayt(npar1),SLayZt(npar1),dprt(npar1) 
      dimension Hlay(npar1),Slay(npar1),SLayZ(npar1),dpr(npar1) 
      dimension SLayNormt(npar1),SLayZNormt(npar1)
      dimension SLayNorm(npar1),SLayZNorm(npar1)
      Klayt=0
      iend=Klay
      do i=1,Klay-1
      if (Hlay(i).lt.ZNE) goto 1
      if (Hlay(i).gt.ZVE) then
      iend=i
      goto 2
      endif
      Klayt=Klayt+1
      Hlayt(Klayt)=Hlay(i)
      Slayt(Klayt)=SLay(i)
      SLayZt(Klayt)=SLayZ(i)
      dprt(Klayt)=dpr(i)
      SLayNormt(Klayt)=SLayNorm(i)     
      SLayZNormt(Klayt)=SLayZNorm(i)    
1     enddo
2     Klayt=Klayt+1
      SLayt(Klayt)=SLay(iend)
      SLayZt(Klayt)=SLayZ(iend)
      dprt(Klayt)=dpr(iend)
      SLayNormt(Klayt)=SLayNorm(iend)     
      SLayZNormt(Klayt)=SLayZNorm(iend)    
! Reassignment of values
      Klay=Klayt
      do i=1,Klay-1
      Hlay(i)=Hlayt(i)
      Slay(i)=SLayt(i)
      SLayZ(i)=SLayZt(i)
      dpr(i)=dprt(i)
      SLayNorm(i)=SLayNormt(i)     
      SLayZNorm(i)=SLayZNormt(i)          
      enddo  
      Slay(Klay)=SLayt(Klay)
      SLayZ(Klay)=SLayZt(Klay)
      dpr(Klay)=dprt(Klay)
      SLayNorm(Klay)=SLayNormt(Klay)     
      SLayZNorm(Klay)=SLayZNormt(Klay)          
      return
      end    
              