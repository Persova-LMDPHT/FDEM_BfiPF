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
!   This file contains functions for reading parameters for building a 2D mesh
!
!   Written by Dr. Marina Persova
!   Novosibirsk State Technical University   
!   Novosibirsk, 630073, Russia
!   mpersova@mail.ru
!   Version 2.0 December 10, 2024
!  
!======================================================================! 
!
! reading the settings from file RegularMeshBuilderSettings2D.cfg
      subroutine ReadSettings()      
      implicit real*8(a-h,o-z)
! size of the calculation domain (distance from receivers)
      common/settings2d/bak2d
! mesh construction parameters  
! meshes parameters for Bfi-formulation
      common/settings2d_ArH/rmin2d_ArH,zmin2d_ArH,
     *ukoef1_2d_ArH,ukoef2_2d_ArH,ukoef3_2d_ArH
! meshes parameters for AV-formulation
      common/settings2d_AxH/rmin2d_AxH,zmin2d_AxH,
     *ukoef1_2d_AxH,ukoef2_2d_AxH,ukoef3_2d_AxH                
      open(2,err=4,file='RegularMeshBuilderSettings2D.cfg',status='old')    
! computational domain size
      read(2,*) bak2d
! meshes parameters for Bfi-formulation
      read(2,*) rmin2d_ArH,zmin2d_ArH,ukoef1_2d_ArH,
     *ukoef2_2d_ArH,ukoef3_2d_ArH
! meshes parameters for AV-formulation
      read(2,*) rmin2d_AxH,zmin2d_AxH,ukoef1_2d_AxH,
     *ukoef2_2d_AxH,ukoef3_2d_AxH   
      close(2)      
      return
4     open(2,file='error_convert')
      write(2,*) 'no file RegularMeshBuilderSettings2D.cfg'
      stop
      close(2)            
      end  
!
! reading file with layer parameters
      subroutine ReadZSig2d() 
      implicit real*8(a-h,o-z)       
 	parameter (npar = 100000)      
	parameter (npar1 = 1000)  
	common/DprLay/dpr(npar1)    
! Hlay - top coordinates of the layers;  SLay - lateral conductivity;  SLayZ - vertical conductivity
      common/Rlayer/Hlay(npar1),SLay(npar1),SLayZ(npar1)
! Klay - the number of layers
      common/dim_Rlayer/Klay
      Klay=0
      open(2,err=4,file='z_sig_2d',status='old')
      read(2,*) Klay      
      do i=1,Klay
      read(2,*) Hlay(i),SLay(i)
      SLayZ(i)=SLay(i)
      dpr(i)=0      
      enddo
      Klay=Klay+1
      SLay(Klay)=0d0
      SLayZ(Klay)=0d0
      dpr(Klay)=0d0
      close(2)  
      
      open(2,err=4,file='mlayersAdditional',status='old')
      read(2,*)
      do i=1,Klay-1
      read(2,*) ro1,ro2,ro3,omut,dpr(Klay-i)
      SLayZ(Klay-i)=1./ro3
      enddo  
      close(2)
     
      return
4     open(2,file='message_convert')
      write(2,*) 'no file z_sig_2d or mlayersAdditional'
c      stop
      close(2)            
      end  
!
! reading the file with installation parameters
      subroutine ReadSetPar()  
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)           
      common/kolGenG/KOL_SOR,KOL_HEL      
      common/ZG_/ZG_(2,npar)                          
      common/TransmitZ/ZminGen,ZmaxGen,ZminGenH,ZmaxGenH   
      common/TransmitallZ/ZG_ALL(2*npar),ZG_ALLH(2*npar)    
      real*4 xp,yp,zp
      pi=3.1415926535897932384626433832795
      KOL_SOR=0
      open(2,file='group')
      read(2,*) kol_gr
      do i=1,kol_gr
      read(2,*) n1,nb,ne
      KOL_SORL=ne-nb+1
      KOL_SOR=KOL_SOR+KOL_SORL
      enddo
      close(2)  
      KOL_HEL=KOL_SOR   
c
c
      open(2,file='sours')
      ZmaxGenH=-1e10
      ZminGenH=1e10           
      do i=1,KOL_SOR
      read(2,*) ax,ay,az,bx,by,bz
      ZG_(1,i)=az
      ZG_(2,i)=az
      ZG_ALLH(i)=az
      if (az.lt.ZminGenH) ZminGenH=az
      if (az.gt.ZmaxGenH) ZmaxGenH=az      
      enddo
      close(2)
      return
      end


