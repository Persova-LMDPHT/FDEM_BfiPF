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
!   This file contains the procedures for ñonstruction and writing 
!   of meshes for calculating the primary field.
!
!   Written by Dr. Marina Persova
!   Novosibirsk State Technical University   
!   Novosibirsk, 630073, Russia
!   mpersova@mail.ru
!   Version 2.0 December 10, 2024
!  
!======================================================================! 
!
!
! ñonstruction and writing of meshes for calculating the primary field
      subroutine MakeMesh1D()
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)     
! current meshes parameters
      common/settings2d/bak2d,rmin2d,zmin2d,
     *ukoef1_2d,ukoef2_2d,ukoef3_2d	
! meshes parameters for Bfi-formulation	
      common/settings2d_ArH/rmin2d_ArH,zmin2d_ArH,
     *ukoef1_2d_ArH,ukoef2_2d_ArH,ukoef3_2d_ArH
! meshes parameters for AV-formulation
      common/settings2d_AxH/rmin2d_AxH,zmin2d_AxH,
     *ukoef1_2d_AxH,ukoef2_2d_AxH,ukoef3_2d_AxH	 
      common/TransmitZ/ZminGen,ZmaxGen,ZminGenH,ZmaxGenH
      common/kolGenG/KOL_SOR,KOL_HEL            
      common/TransmitallZ/ZG_ALL(2*npar),ZG_ALLH(2*npar)        
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1)             
      common/eps_set2d/eps
      dimension ZG_ALLT(2*npar)
! mesh construction for Bfi-formulation *.dat
      do i=1,KLay-1
      ZG_ALLT(i)=Hlay(i)
      enddo
      ii=i
      do i=1,KOL_HEL
      ZG_ALLT(ii)=ZG_ALLH(i)
      ii=ii+1
      enddo
      KOL=KLay-1+KOL_HEL
! assigning the mesh parameters for the Bfi-formulation to the current parameters
      rmin2d=rmin2d_ArH
      zmin2d=zmin2d_ArH
      ukoef1_2d=ukoef1_2d_ArH
      ukoef2_2d=ukoef2_2d_ArH
      ukoef3_2d=ukoef3_2d_ArH
! building a mesh vertically with condensation to the source position and layer boundaries
      call Set2DzNew(Hlay(1),Hlay(KLay-1),ZG_ALLT,0,KOL)          
! mesh construction along the r-axis     
      call Set2Dr1()
! writing mesh files
      call WriteMesh(0)   
      call WriteHELs(0)
      
! mesh construction for AV-formulation *.dat1
! assigning the mesh parameters for the AV-formulation to the current parameters
      rmin2d=rmin2d_AxH
      zmin2d=zmin2d_AxH
      ukoef1_2d=ukoef1_2d_AxH
      ukoef2_2d=ukoef2_2d_AxH
      ukoef3_2d=ukoef3_2d_AxH
! building a mesh vertically with condensation only to the source position
      call Set2DzNew(ZminGenH,ZmaxGenH,ZG_ALLH,0,KOL_HEL)          
! mesh construction along the r-axis     
      call Set2Dr1()
! writing mesh files
      call WriteMesh(1)   
      call WriteHELs(1)
          
      return
      end
!
!
!
! Construction of a 2D mesh with condensation to each surface from ZG_ALL array
      subroutine Set2DzNew(ZminGen,ZmaxGen,ZG_ALL,IFL,kpp)            
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)   
      common/kolGenG/kPOL      	    	   
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1)
      common/settings2d/bak2d,rmin2d,zmin2d,
     *ukoef1_2d,ukoef2_2d,ukoef3_2d
      common/bak2d/ZB1_2D,ZB2_2D     
      common/r0_/r0_
      common/DSizeT/DSizeT         
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
      common/eps_set2d/eps
      common/P1/ZH1(npar1),ZK1(npar1),LZ1(npar1),ZC1(2,npar1),KZ1(npar1)
      common/P2/ZH2(npar1),ZK2(npar1),LZ2(npar1),ZC2(2,npar1),KZ2(npar1)
      common/P3/ZH3(npar1),ZK3(npar1),LZ3(npar1),ZC3(2,npar1),KZ3(npar1)
      dimension ZG_ALL(2*npar)  
      dimension ZG_ALL2(2*npar),ZG_ALL3(2*npar),nZG_ALL3(2*npar) 
! Sorting the z-coordinates, leaving only the distinct z-coordinates    
      kPOLt=0
      do i=1,kpp      
      call InsertinPlist(kPOLt,ZG_ALL(i),ZG_ALL2)
      enddo
! Max and min coordinates      
      ZminGen=ZG_ALL2(1)
      ZmaxGen=ZG_ALL2(kPOLt)
! Insert intermediate z-coordinates from which condensation is performed
      kPOLtn=1
      ZG_ALL3(1)=ZG_ALL2(1)
      nZG_ALL3(1)=1
      do i=2,kPOLt
      if ((ZG_ALL2(i)-ZG_ALL2(i-1)).gt.(zmin2d+eps)) then
      kPOLtn=kPOLtn+1
      ZG_ALL3(kPOLtn)=(ZG_ALL2(i)+ZG_ALL2(i-1))*0.5d0
      nZG_ALL3(kPOLtn)=0
      kPOLtn=kPOLtn+1
      ZG_ALL3(kPOLtn)=ZG_ALL2(i)
      nZG_ALL3(kPOLtn)=1      
      else
      kPOLtn=kPOLtn+1
      ZG_ALL3(kPOLtn)=ZG_ALL2(i)
      nZG_ALL3(kPOLtn)=1            
      endif
      enddo
! build a mesh in Z
      KZMC=1
      ZMC(KZMC)=ZB1_2D
      call SetZN(ZminGen,ukoef3_2d,zmin2d,ZB1_2D)  
      nFQ=1
      do i=2,kPOLtn
      if (nZG_ALL3(i).eq.0) then
      call SetZV(ZG_ALL3(i-1),ukoef2_2d,zmin2d,ZG_ALL3(i))
      nFQ=0
      goto 5
      endif
      if ((nZG_ALL3(i).eq.1).and.(nFQ.eq.0)) then
      call SetZN(ZG_ALL3(i),ukoef2_2d,zmin2d,ZG_ALL3(i-1))
      nFQ=1
      goto 5
      endif
      if ((nZG_ALL3(i).eq.1).and.(nFQ.eq.1)) then
      call SetZM(ZG_ALL3(i-1),ukoef2_2d,zmin2d,ZG_ALL3(i))
      nFQ=1
      goto 5
      endif
5     enddo 
      call SetZV(ZmaxGen,ukoef3_2d,zmin2d,ZB2_2D)
50    return
      end
!
!
! mesh construction along the r-axis     
      subroutine Set2Dr1()   
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/settings2d/bak2d,rmin2d,zmin2d,
     *ukoef1_2d,ukoef2_2d,ukoef3_2d
      common/r0_/r0_
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
! mesh construction with spacing to the remote boundary
      HR=rmin2d
      call makeset(r0_,bak2d,HR,ukoef1_2d,1,KOLR)
      KRMC=1
      RMC(KRMC)=r0_
      do j=1,KOLR
      KRMC=KRMC+1
      RMC(KRMC)=RMC(KRMC-1)+HR
      HR=HR*ukoef1_2d
      enddo
      return
      end
!      
! initialization of mesh parameters on a subinterval
      subroutine initset(K1,ZC11,ZC12,ZH1,ZK1,LZ1,Z1,Z2,H,UK,LN)
      implicit real*8(a-h,o-z)
      ZC11=Z1
      ZC12=Z2
      ZH1=H
      ZK1=UK
      LZ1=LN
      K1=K1+1
      return
      end
!
! Build a mesh down from the condensation surfaces
      subroutine SetZN(ZminGen,ukoef1_2d,zmin2d,bak2d)
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1)
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
      common/eps_set2d/eps
      common/P1/ZH1(npar1),ZK1(npar1),LZ1(npar1),ZC1(2,npar1),KZ1(npar1)
      Hbeg=zmin2d   
      K1=1
      ifl1=0
      do i=Klay-1,1,-1
      if (Hlay(i).gt.(ZminGen+eps)) goto 1
      if (Hlay(i).lt.(bak2d-eps)) goto 1
      if (abs(Hlay(i)-ZminGen).lt.eps) then
      if (i.eq.1) then
! Build a mesh to a remote boundary
      call initset(K1,ZC1(1,K1),ZC1(2,K1),ZH1(K1),ZK1(K1),LZ1(K1),
     *bak2d,Hlay(1),zmin2d,ukoef1_2d,-1)
      call makeset(ZC1(1,K1-1),ZC1(2,K1-1),ZH1(K1-1),
     *ZK1(K1-1),LZ1(K1-1),KZ1(K1-1))
      ifl1=1
      endif
      else      
      call initset(K1,ZC1(1,K1),ZC1(2,K1),ZH1(K1),ZK1(K1),LZ1(K1),
     *Hlay(i),ZminGen,zmin2d,ukoef1_2d,-1)
      call makeset(ZC1(1,K1-1),ZC1(2,K1-1),ZH1(K1-1),
     *ZK1(K1-1),LZ1(K1-1),KZ1(K1-1))
      Hbeg=ZH1(K1-1)
      endif
      ibeg=i
      goto 2     
1     enddo
! All layers above ZminGen
! Build a grid to a remote boundary
      call initset(K1,ZC1(1,K1),ZC1(2,K1),ZH1(K1),ZK1(K1),LZ1(K1),
     *bak2d,ZminGen,zmin2d,ukoef1_2d,-1)
      call makeset(ZC1(1,K1-1),ZC1(2,K1-1),ZH1(K1-1),
     *ZK1(K1-1),LZ1(K1-1),KZ1(K1-1))
      ifl1=1            
2     if (ifl1.eq.0) then
3     if (ibeg.gt.1) then
      if (Hlay(ibeg-1).gt.(bak2d+eps)) then 
      call initset(K1,ZC1(1,K1),ZC1(2,K1),ZH1(K1),ZK1(K1),LZ1(K1),
     *Hlay(ibeg-1),Hlay(ibeg),Hbeg,ukoef1_2d,-1)
      call makeset(ZC1(1,K1-1),ZC1(2,K1-1),ZH1(K1-1),
     *ZK1(K1-1),LZ1(K1-1),KZ1(K1-1))
      Hbeg=ZH1(K1-1)
      ibeg=ibeg-1
      goto 3
      endif
      endif
      call initset(K1,ZC1(1,K1),ZC1(2,K1),ZH1(K1),ZK1(K1),LZ1(K1),
     *bak2d,Hlay(ibeg),Hbeg,ukoef1_2d,-1)
      call makeset(ZC1(1,K1-1),ZC1(2,K1-1),ZH1(K1-1),
     *ZK1(K1-1),LZ1(K1-1),KZ1(K1-1))
      ifl1=1                  
      endif
!      
      do i=K1-1,1,-1      
      do j=1,KZ1(i)
      KZMC=KZMC+1
      ZMC(KZMC)=ZMC(KZMC-1)+ZH1(i)
      ZH1(i)=ZH1(i)*ZK1(i)
      enddo
      enddo
!      
      return
      end
!
! Construct a mesh between the condensation surfaces
      subroutine SetZM(ZminGen,ukoef1_2d,zmin2d,ZmaxGen)
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1)
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
      common/eps_set2d/eps
      common/P2/ZH2(npar1),ZK2(npar1),LZ2(npar1),ZC2(2,npar1),KZ2(npar1)
      K2=1
      ZZ1=ZminGen
      Hbeg=zmin2d
      ukoef0=1d0
      do i=1,KLay-1
      if ((Hlay(i).gt.(ZminGen+eps)).
     *and.(Hlay(i).lt.(ZmaxGen-eps))) then
      call initset(K2,ZC2(1,K2),ZC2(2,K2),ZH2(K2),ZK2(K2),LZ2(K2),
     *ZZ1,Hlay(i),Hbeg,ukoef0,1)
      call makeset(ZC2(1,K2-1),ZC2(2,K2-1),ZH2(K2-1),
     *ZK2(K2-1),LZ2(K2-1),KZ2(K2-1))
      Hbeg=ZH2(K2-1)*(ZK2(K2-1)**(KZ2(K2-1)-1))
      ZZ1=Hlay(i)
      endif
      enddo      
      call initset(K2,ZC2(1,K2),ZC2(2,K2),ZH2(K2),ZK2(K2),LZ2(K2),
     *ZZ1,ZmaxGen,Hbeg,ukoef0,1)
      call makeset(ZC2(1,K2-1),ZC2(2,K2-1),ZH2(K2-1),
     *ZK2(K2-1),LZ2(K2-1),KZ2(K2-1))
!
!
      do i=1,K2-1
      do j=1,KZ2(i)
      KZMC=KZMC+1
      ZMC(KZMC)=ZMC(KZMC-1)+ZH2(i)
      ZH2(i)=ZH2(i)*ZK2(i)
      enddo
      enddo
      return
      end
!
! Build a mesh up from the condensation surfaces
      subroutine SetZV(ZmaxGen,ukoef1_2d,zmin2d,bak2d)
      implicit real*8(a-h,o-z)
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1)
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
      common/eps_set2d/eps
      common/P3/ZH3(npar1),ZK3(npar1),LZ3(npar1),ZC3(2,npar1),KZ3(npar1)
      Hbeg=zmin2d      
      K3=1
      ifl3=0
      
      do i=1,KLay-1
      if (Hlay(i).lt.(ZmaxGen-eps)) goto 11
      if (Hlay(i).gt.(bak2d+eps)) goto 11
      if (abs(Hlay(i)-ZmaxGen).lt.eps) then
      if (i.eq.KLay-1) then
! Build a mesh to a remote boundary
      call initset(K3,ZC3(1,K3),ZC3(2,K3),ZH3(K3),ZK3(K3),LZ3(K3),
     *Hlay(KLay-1),bak2d,zmin2d,ukoef1_2d,1)
      call makeset(ZC3(1,K3-1),ZC3(2,K3-1),ZH3(K3-1),
     *ZK3(K3-1),LZ3(K3-1),KZ3(K3-1))
      ifl3=1
      endif
      else
      call initset(K3,ZC3(1,K3),ZC3(2,K3),ZH3(K3),ZK3(K3),LZ3(K3),
     *ZmaxGen,Hlay(i),zmin2d,ukoef1_2d,1)
      call makeset(ZC3(1,K3-1),ZC3(2,K3-1),ZH3(K3-1),
     *ZK3(K3-1),LZ3(K3-1),KZ3(K3-1))
      Hbeg=ZH3(K3-1)*(ZK3(K3-1)**(KZ3(K3-1)-1))
      endif
      ibeg=i
      goto 12     
11    enddo
! All layers below ZmaxGen
! Build a mesh to a remote boundary
      call initset(K3,ZC3(1,K3),ZC3(2,K3),ZH3(K3),ZK3(K3),LZ3(K3),
     *ZmaxGen,bak2d,zmin2d,ukoef1_2d,1)
      call makeset(ZC3(1,K3-1),ZC3(2,K3-1),ZH3(K3-1),
     *ZK3(K3-1),LZ3(K3-1),KZ3(K3-1))
      ifl3=1            
12    if (ifl3.eq.0) then
13    if (ibeg.lt.(KLay-1)) then
      if (HLay(ibeg+1).lt.(bak2d-eps)) then 
      call initset(K3,ZC3(1,K3),ZC3(2,K3),ZH3(K3),ZK3(K3),LZ3(K3),
     *HLay(ibeg),HLay(ibeg+1),Hbeg,ukoef1_2d,1)
      call makeset(ZC3(1,K3-1),ZC3(2,K3-1),ZH3(K3-1),
     *ZK3(K3-1),LZ3(K3-1),KZ3(K3-1))
      Hbeg=ZH3(K3-1)*(ZK3(K3-1)**(KZ3(K3-1)-1))
      ibeg=ibeg+1
      goto 13
      endif
      endif      
      call initset(K3,ZC3(1,K3),ZC3(2,K3),ZH3(K3),ZK3(K3),LZ3(K3),
     *HLay(ibeg),bak2d,Hbeg,ukoef1_2d,1)
      call makeset(ZC3(1,K3-1),ZC3(2,K3-1),ZH3(K3-1),
     *ZK3(K3-1),LZ3(K3-1),KZ3(K3-1))
      ifl3=1                  
      endif
c
      do i=1,K3-1
      do j=1,KZ3(i)
      KZMC=KZMC+1
      ZMC(KZMC)=ZMC(KZMC-1)+ZH3(i)
      ZH3(i)=ZH3(i)*ZK3(i)
      enddo
      enddo
      return
      end
!
!
! mesh construction in the interval tbeg - ts, ht - step from tbeg, tkoef - coefficient,
! lnt - direction of condensation, 
! lnt=1: the mesh is built using the initial parameters,
! lnt=-1: the obtained parameters are inverted, 
! kt - number of steps
	subroutine makeset (tbeg,ts,ht,tkoef,lnt,kt)
      implicit double precision (a-h,o-z)
      common/eps_set2d/eps
      epsminstep=eps
      epsset=eps
! equality check
	if (ts-tbeg.le.epsminstep) then
	ht=0
	kt=0
	return
	endif
! initial step determination procedure
	if (ts-tbeg.le.ht*(1.+epsset)) then
	ht=ts-tbeg
	kt=1
	return
	endif
	if (abs(tkoef-1d0).le.epsset) then
	kt=(ts-tbeg)/ht
	if ((kt+1)*ht-(ts-tbeg).lt.(2d0/3d0)*ht) kt=kt+1
	if (kt.gt.0) then
	htsloe=(ts-tbeg)/kt
	else
	kt=1
	htsloe=ht
	endif
	else
	kt=0
	htt=ht
	t=tbeg
	dt=ts-tbeg
	i=0
2	if(t.ge.ts) then
	goto 3
	else
	t=t+htt*(tkoef**i)
	i=i+1
	goto 2
	endif
3     hlast=htt*(tkoef**(i-1))
	if ((t-ts).lt.(2d0/3d0)*hlast) then
	if (lnt.gt.0) then
	htsloe=dt*(1-tkoef)/(1-tkoef**i)	
	else
	htsloe=dt*(1-tkoef)/(1-tkoef**i)	
	htsloe=htsloe*(tkoef**(i-1))	
	tkoef=1./tkoef
	endif
	kt=i
	else
      if (lnt.gt.0) then
	htsloe=dt*(1-tkoef)/(1-tkoef**(i-1))	
	else
	htsloe=dt*(1-tkoef)/(1-tkoef**(i-1))	
	htsloe=htsloe*(tkoef**(i-2))	
	tkoef=1./tkoef
	endif
	kt=i-1
	endif
	endif
	ht=htsloe
	return
	end
!
!
!   Writing current value to files
      subroutine WriteHELs(NB)      
	parameter (npar = 100000)      
	parameter (npar1 = 1000) 
	implicit real*8(a-h,o-z)     
      common/kolGenG/KOL_SOR,KOL_HEL      
      character*1 Ext
      if (NB.eq.0) then
      call WriteHELs0()
      return
      endif
      write(Ext,'(i1.1)') NB
      TokBase=1d0
      open(3,file='currentval'//Ext)
      do i=1,KOL_HEL
      write(3,*) TokBase
      enddo
      close(3)
      return
      end
!
      subroutine WriteHELs0()      
	parameter (npar = 100000)      
	parameter (npar1 = 1000) 
	implicit real*8(a-h,o-z)     
      common/kolGenG/KOL_SOR,KOL_HEL      
      TokBase=1d0
      open(3,file='currentval')
      do i=1,KOL_HEL
      write(3,*) TokBase
      enddo
      close(3)
      return
      end            
!
! Writing a mesh to files
      subroutine WriteMesh(N)
      implicit real*8(a-h,o-z)   
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1),SLayZ(npar1)
      common/RlayerSN/SLayNorm(npar1),SLayZNorm(npar1)            
      common/dpr/dpr0      
      common/umu/umu0 
      character*1 Ext
!
!
      if (N.eq.0) then
      call WriteMesh0()
      return
      endif
!
!
      open(2,file='sigma')
      nm=0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,SLayNorm(i)
      enddo
      close(2) 
      open(2,file='sigmaZ')
      nm=0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,SLayZNorm(i)
      enddo
      close(2) 
      open(2,file='dpr')
      nm=0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,dpr0
      enddo
      close(2) 
      open(2,file='mu')
      nm=0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,umu0
      enddo
      close(2) 
!
!      
      write(Ext,'(i1.1)') N
	open(2,file='rz.dat'//Ext,access='direct',recl=16)
	do i=1,KZMC
	do j=1,KRMC
	write(2,rec=(i-1)*KRMC+j) RMC(j),ZMC(i) 
	enddo
	enddo
	close(2)
!
	open(2,file='nvtr.dat'//Ext,access='direct',recl=24)
	do i=1,KZMC-1
	do j=1,KRMC-1
	write(2,rec=(i-1)*(KRMC-1)+j) i*KRMC+j,i*KRMC+j+1,(i-1)*KRMC+j,
     *(i-1)*KRMC+j+1,0,1
	enddo
	enddo
	close(2)
!	
	open(3,file='nvkat2d.dat'//Ext,access='direct',recl=4)
	do i=1,KZMC-1
	zz=(ZMC(i)+ZMC(i+1))/2.
	call MakeMat(zz,nmat)
	do j=1,KRMC-1
	write(3,rec=(i-1)*(KRMC-1)+j) nmat
	enddo
	enddo
	close(3)
!
      open(2,file='rz.txt'//Ext)
      write(2,*) KRMC,KZMC
      close(2)
      open(2,file='r.dat'//Ext,access='direct',recl=8)
      do i=1,KRMC
      write(2,rec=i) RMC(i)
      enddo
      close(2)
      open(2,file='z.dat'//Ext,access='direct',recl=8)
      do i=1,KZMC
      write(2,rec=i) ZMC(i)
      enddo
      close(2)
!
!     
!
      open(2,file='tsize.dat'//Ext)
      write(2,*) '0'
      write(2,*)
      close(2)
!
      open(2,file='l1.dat'//Ext,access='direct',recl=4)
      KL1=0
      do i=1,KRMC
      KL1=KL1+1
      write(2,rec=KL1) i
      KL1=KL1+1
      write(2,rec=KL1) (KZMC-1)*KRMC+i
      enddo
      do i=2,KZMC-1
      KL1=KL1+1
      write(2,rec=KL1) i*KRMC
      if (N.ne.1) then
      KL1=KL1+1
      write(2,rec=KL1) (i-1)*KRMC+1
      endif      
      enddo
      close(2)
!      
!
	kuzlov=KRMC*KZMC
	kpram=(KRMC-1)*(KZMC-1)
!
      nn4=0
!
       open (1,file='inf2tr.dat'//Ext,status='unknown')
      write(1,105)0,0,1
105    format(' islau=',i8,' indku1=',i8,' indfro=',i8)
       write(1,101)kuzlov,kpram,KL1,0,0
101    format('kuzlov=',i8,'    ktr=',i8,'    kt1=',i8,'  kreb2=',i8,  
     *'  kreb3=',i8)
       write(1,102)2,2,2,8
102    format('kisrr1=',i8,' kisrr2=',i8,' kisrr3=',i8,'  kbrsr=',i8)
       write(1,103)nn4
103    format(' kreb4=',I8)
 
100   format(7x,i8,4(8x,i8))
       close(1)
!      
!
      return
      end
!
!
      subroutine WriteMesh0()
      implicit real*8(a-h,o-z)   
	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/DIMMesh1D/KRMC,KZMC
      common/Mesh1D/RMC(npar),ZMC(npar)
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1),SLayZ(npar1)
      common/RlayerSN/SLayNorm(npar1),SLayZNorm(npar1)            
      common/dpr/dpr0      
      common/umu/umu0 
c
c
      open(2,file='sigma')
      nm=0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,SLayNorm(i)
      enddo
      close(2) 
      open(2,file='sigmaZ')
      nm=0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,SLayZNorm(i)
      enddo
      close(2) 
      open(2,file='dpr')
      nm=0
c      write(2,*) nm,dpr0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,dpr0
      enddo
      close(2) 
      open(2,file='mu')
      nm=0
c      write(2,*) nm,umu0
      do i=KLay,1,-1
      nm=nm+1
      write(2,*) nm,umu0
      enddo
      close(2) 
c
c      
	open(2,file='rz.dat',access='direct',recl=16)
	do i=1,KZMC
	do j=1,KRMC
	write(2,rec=(i-1)*KRMC+j) RMC(j),ZMC(i) 
	enddo
	enddo
	close(2)
c
	open(2,file='nvtr.dat',access='direct',recl=24)
	do i=1,KZMC-1
	do j=1,KRMC-1
	write(2,rec=(i-1)*(KRMC-1)+j) i*KRMC+j,i*KRMC+j+1,(i-1)*KRMC+j,
     *(i-1)*KRMC+j+1,0,1
	enddo
	enddo
	close(2)
c	
	open(3,file='nvkat2d.dat',access='direct',recl=4)
	do i=1,KZMC-1
	zz=(ZMC(i)+ZMC(i+1))/2.
	call MakeMat(zz,nmat)
	do j=1,KRMC-1
	write(3,rec=(i-1)*(KRMC-1)+j) nmat
	enddo
	enddo
	close(3)
c
      open(2,file='rz.txt')
      write(2,*) KRMC,KZMC
      close(2)
      open(2,file='r.dat',access='direct',recl=8)
      do i=1,KRMC
      write(2,rec=i) RMC(i)
      enddo
      close(2)
      open(2,file='z.dat',access='direct',recl=8)
      do i=1,KZMC
      write(2,rec=i) ZMC(i)
      enddo
      close(2)
c

c
      open(2,file='tsize.dat')
      write(2,*) '0'
      write(2,*)
      close(2)
c
      open(2,file='l1.dat',access='direct',recl=4)
      KL1=0
      do i=1,KRMC
      KL1=KL1+1
      write(2,rec=KL1) i
      KL1=KL1+1
      write(2,rec=KL1) (KZMC-1)*KRMC+i
      enddo
      do i=2,KZMC-1
      KL1=KL1+1
      write(2,rec=KL1) i*KRMC
c      if (N.ne.1) then
c      KL1=KL1+1
c      write(2,rec=KL1) (i-1)*KRMC+1
c      endif      
      enddo
      close(2)
      
c
	kuzlov=KRMC*KZMC
	kpram=(KRMC-1)*(KZMC-1)
c
      nn4=0
c
       open (1,file='inf2tr.dat',status='unknown')
      write(1,105)0,0,1
105    format(' islau=',i8,' indku1=',i8,' indfro=',i8)
       write(1,101)kuzlov,kpram,KL1,0,0
101    format('kuzlov=',i8,'    ktr=',i8,'    kt1=',i8,'  kreb2=',i8,  
     *'  kreb3=',i8)
       write(1,102)2,2,2,8
102    format('kisrr1=',i8,' kisrr2=',i8,' kisrr3=',i8,'  kbrsr=',i8)
       write(1,103)nn4
103    format(' kreb4=',I8)
 
100   format(7x,i8,4(8x,i8))
       close(1)
      

      return
      end
!
! material number definition
	subroutine MakeMat(zz,nmat)
      implicit real*8(a-h,o-z)   
 	parameter (npar = 100000)      
	parameter (npar1 = 1000)      
      common/dim_Rlayer/Klay
      common/Rlayer/Hlay(npar1),SLay(npar1),SLayZ(npar1)
      matbeg=Klay
      if (abs(zz).lt.1) then
      eps=1e-5
      else 
      eps=zz*1e-5    
      endif
      do i=1,KLay-1
      if (zz.lt.(Hlay(i)+eps)) then
      nmat=matbeg
      goto 1
      endif
      matbeg=matbeg-1
      enddo
      nmat=1
1     return
      end
      
! Sorting Plist array    
       subroutine InsertinPlist(Kol,Q,Plist)
       implicit real*8(a-h,o-z)
       common/eps_set2d/eps
       dimension Plist(*)  
       j=1 
2      if (j.le.(Kol-1)) then
       if (abs(Q-Plist(j)).lt.eps) goto 1
!       
       if ((Q.gt.(Plist(j))).and.(Q.lt.(Plist(j+1)))) then    
       if (abs(Q-Plist(j+1)).lt.eps) goto 1
       do m=Kol,j+1,-1
       Plist(m+1)=Plist(m)
       enddo
       Plist(j+1)=Q
       Kol=Kol+1
       goto 1
       endif
!       
       j=j+1
       goto 2
       endif 
!            
       if (Kol.ne.0) then
       if (abs(Q-Plist(Kol)).lt.eps) goto 1
       endif
!
       if (Kol.ne.0) then
       if (Q.lt.(Plist(1))) then
       do m=Kol,1,-1
       Plist(m+1)=Plist(m)
       enddo
       Plist(1)=Q
       Kol=Kol+1
       goto 1
       endif
       endif
!       
       Kol=Kol+1
       Plist(Kol)=Q
1      return
       end

