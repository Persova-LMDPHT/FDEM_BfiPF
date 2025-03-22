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
!   This file contains function calls for constructing a 2D mesh. 
!
!   Written by Dr. Marina Persova
!   Novosibirsk State Technical University   
!   Novosibirsk, 630073, Russia
!   mpersova@mail.ru
!   Version 2.0 December 10, 2024
!  
!======================================================================! 

      implicit real*8(a-h,o-z)
! maximum number of steps for each coordinate
	parameter (npar = 100000)      
! maximum number of layers      
	parameter (npar1 = 1000)     
! relative magnetic and dielectric permeabilities
      common/dpr/dpr0      
      common/umu/umu0
! value for checking for matching values      
      common/eps_set2d/eps
! initial coordinate along the r-axis
      common/r0_/r0_
!
      r0_=1d-2
      eps=1d-6 
      dpr0=0d0
      umu0=1d0                
! reading the settings from file RegularMeshBuilderSettings2D.cfg
      call ReadSettings()
! reading file with layer parameters
      call ReadZSig2d()  
! reading the file with installation parameters
      call ReadSetPar()  
! auxiliary procedure
      call MakeStructMesh()           
! correction of the model in case the layer sizes are larger than the computational domain size
      call ModelCorrection()  
! construction and writing of meshes for calculating the primary field
      call MakeMesh1D()                     
      end      


         