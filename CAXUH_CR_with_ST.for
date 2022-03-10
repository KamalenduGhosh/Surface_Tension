c **********************************************************************
c 2D Axisymmetric finite deformation formulation
c 7-node Crouzeix-Raviart element
c incompressible non-Gaussian material + Iron/ferrofluid and compressible Air
c hybrid formulation
c W=3^(1-a1)/2a1*mu1*[I1^a1-N^a1]+3^(1-a2)/2a2*mu2*[I1^a2-N^a2]-(mu1+mu2)lnJ+mup/2(J-1)^2+(AMKE-AEPS)/2I4-AMKE/2I5  for PDMS
c W=mu/2[I1-3]-S(I5)   for iron/ferrofluid core
c W=mu/2[I1-3]+lambda/2(J-1)^2-muln(J)
c **********************************************************************

c **********************************************************************
c User element statement in the input file  for incompressible materials:
c
*User Element, Nodes=7, Type=U1, IProperties=1, Properties=11, Coordinates=2, Variables=1
c 1,2,11
c 2,1,2,11
c 3,1,2,11
c 4,1,2,11
c 5,1,2,11
c 6,1,2,11
c 7,1,2,11,12,13,14
**
**
c User element statement in the input file  for compressible materials: 
c
*User Element, Nodes=7, Type=U2, IProperties=1, Properties=11, Coordinates=2, Variables=1
c  1,2,11
c  2,1,2,11
c 3,1,2,11
c 4,1,2,11
c 5,1,2,11
c 6,1,2,11
c 7,1,2,11
**
**
c
c Uvarm Variables (used for visualization)
c
c  1) S11
c  2) S12
c  3) S21
c  4) S22
c  5) D1
c  6) D2
c  7) D3
c
c Material Properties Vector for PDMS
c
c  mu 1= props(1)       ! Initial Shear modulus 1
c  lambda  = props(2)   ! Lame first parameter
c  AEPS= props(3)       ! permeability of air=permb of matrix
c  AMKE=props(4)        ! permeability of air=permb of matrix
c  c=props(5)           ! concentration of particles c=0 for only PDMS
c  ms=props(6)          ! doesn't take part in the calculation for PDMS
c  medium=props(7)      ! medium=1 for PDMS
c  props(8)=N           ! N=3 for 3D domain
c  alpha1=props(9)      ! a1
c  alpha2=props(10)     ! a2
c  mu2=props(11)        ! Initial Shear modulus 2
c  radius=props(12)     ! Radius of the Inclusion
c  gamma=props(13)      ! Surface Tension at the inclusion interface
c  nsri = 0             ! reduced integration flag
c
c Material Properties Vector for iron/ferrofluid core 
c
c  mu 1= props(1)       ! Initial Shear modulus 1=shear modulus of iron/ferrofluid
c  lambda  = props(2)   ! Lame first parameter of iron/ferrofluid
c  AEPS= props(3)       ! =permeability of iron/ferrofluid--- change if AMKE neq AEPS
c  AMKE=props(4)        ! permeability of iron/ferrofluid
c  c=props(5)           ! any value, doesnt take part in constitutive model
c  ms=props(6)          ! magnetic saturation ms
c  medium=props(7)      ! medium=2 for iron/ferrofluid
c  props(8)=N           ! N=3 for 3D domain
c  alpha1=props(9)      ! a1=1 for gaussian Neo Hookean
c  alpha2=props(10)     ! a2=1 for gaussian Neo Hookean
c  mu2=props(11)        ! Initial Shear modulus 2=0 for gaussian Neo Hookean
c  nsri = 0             ! reduced integration flag
c
cc Material Properties Vector for compressible Air
c
c  mu 1= props(1)       ! Initial Shear modulus 1=shear modulus of Air
c  lambda  = props(2)   ! Lame first parameter of Air
c  AEPS= props(3)       ! =permeability of Air--- change if AMKE neq AEPS
c  AMKE=props(4)        ! permeability of Air
c  c=props(5)           ! any value, doesnt take part in constitutive model
c  ms=props(6)          ! magnetic saturation ms (Any value does not matter)
c  medium=props(7)      ! medium=0 for compressible air
c  props(8)=N           ! N=3 for 3D domain
c  alpha1=props(9)      ! a1=1 for gaussian Neo Hookean
c  alpha2=props(10)     ! a2=1 for gaussian Neo Hookean
c  mu2=props(11)        ! Initial Shear modulus 2=0 for gaussian Neo Hookean
c  nsri = 0             ! reduced integration flag
c
c **********************************************************************


      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and
      !    that value must be set here the same.
      !
      !  nUvrm
      !   Total number of Uvarm variables

      integer numElem,ElemOffset,err,nUvrm

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=232786)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=232801)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter(nUvrm=13)

      real*8, allocatable :: globalSdv(:,:,:)

      end module global
      
***********************************************************************
      subroutine uel (rhs, amatrx, svars, energy, ndofel, nrhs,
     1 nsvars, props, nprops, coords, mcrd, nnode, u, du, v, a, jtype,
     2 time, dtime, kstep, kinc, jelem, params, ndload, jdltyp,
     3 adlmag, predef, npredf, lflags, mlvarx, ddlmag, mdload, pnewdt,
     4 jprops, njprop, period)

      ! This subroutine is called by ABAQUS
      ! Subroutine for 10-node tetrahedral element

      use global
      include 'aba_param.inc'
c
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),
     1 svars(*),energy(8),coords(mcrd,nnode),u(ndofel),
     2 du(mlvarx,*),v(ndofel),a(ndofel),time(2),params(*),
     3 jdltyp(mdload,*),adlmag(mdload,*),ddlmag(mdload,*),
     4 predef(2,npredf,nnode),lflags(*),jprops(*)
c
      dimension xi(mcrd),nshp(nnode),F(mcrd+1,mcrd+1)
      dimension xjac(mcrd,mcrd),xjacinv(mcrd,mcrd)
      dimension E(mcrd+1),Ft(mcrd+1,mcrd+1)
      dimension PV(mcrd+1,mcrd+1)
      dimension T(mcrd+1,mcrd+1)
      dimension pgauss(ndofel),sgauss(ndofel,ndofel)
      dimension psgauss(ndofel),ssgauss(ndofel,ndofel)
      dimension wsgp(5), xisgp(5)
      real*8, allocatable :: wgp(:),xigp(:,:),gmtx(:,:),gmtxt(:,:),
     + xfe(:),xke(:,:),auxm1(:,:),AI(:),R(:),R_s(:),
     + xkesurften(:,:),xfesurften(:),auxm2(:,:),gmtxs(:,:),gmtxst(:,:)
      real*8 RXi,RXis
      dimension Fsurf(mcrd+1,mcrd+1),refn(mcrd+1),currn(mcrd+1)
      dimension Po(mcrd+1,mcrd+1),P(mcrd+1,mcrd+1)
c      dimension xFinvtN(mcrd)
      real*8 detjacsurf,sarea ,gam_st, Jhat


      ! iwr if set to 1, element matrices (stiffness, load, etc...) are written in
      ! the .msg file LOT OF MEMORY AND COMPUTATIONAL TIME REQUIRED
      ! ngp number of Gauss' points
      ! imat material number ???
      ! ninv number of invariants
      ! ndf number of dofs per node
      ! ivarm uvarm variable counter
      ! nsri selective reduced integration flag 0 normal 1 reduced
      iwr=0
      imat=props(7)
c      write(*,*),'imat',imat,'nnode',nnode,'ndofel',ndofel
      ninv=5
      ndpress=3
      ndf=3
      nsri=0
      if (imat .eq. 0) then
         nsize=ndf*mcrd    ! 6 partial derivatives
      else 
         nsize=ndf*mcrd+1   ! 6 partial derivatives and 1 pressure
      end if

c      if (imat .neq. props(7)) then
c         write(*,*),'imat not equal to jprops props(7)'
c         write(7,*),'imat not equal to props(7)'
c         stop
c      end if
c
c Allocate memory for element matrices
c
      allocate(gmtx(nsize+1,ndofel),gmtxt(ndofel,nsize+1))
      allocate(xfe(nsize+1),xke(nsize+1,nsize+1),auxm1(ndofel,nsize+1))
      allocate(gmtxs(nsize+1,ndofel),gmtxst(ndofel,nsize+1))
      allocate(xfesurften(nsize+1),xkesurften(nsize+1,nsize+1))
      allocate(auxm2(ndofel,nsize+1))
      allocate(AI(ninv))
c
c Allocate memory for globalSdv allocatable array for the uvarm variables
c     
      ngp=6
      allocate(wgp(ngp))
      allocate(xigp(mcrd,ngp))
      allocate(R(ngp))
      allocate(R_s(5))
      allocate(globalSdv(numElem,ngp,nUvrm),stat=err)
c      call kinitia(globalSdv,numElem*ngp*nUvrm)
c
c Gauss Integration
c
      call kinitia(rhs,mlvarx*1)
      call kinitia(amatrx,ndofel*ndofel)
      call kinitia(energy,8)
      call kinitia(R,ngp)
      call kinitia(R_s,5)

      if (lflags(3).eq.1.or.lflags(3).eq.5) then !stiff and res or res
      if (lflags(1).eq.1.or.lflags(1).eq.71) then
      do ksri=1,(nsri+1) !only integrate with full integration the deviatoric part for ksri=1
                         !then reduce integration for the volumetric part for ksri=2
c         if (ksri.eq.2) then !activate only 1 Gauss Point
c            deallocate(wgp,xigp)
c            ngp=1
c            allocate(wgp(ngp))
c            allocate(xigp(mcrd,ngp))
c         end if

         call kgausspoints(wgp,xigp,ngp,mcrd)
c Calculate R
         do igauss=1,ngp        !Loop over Gauss points
            do i=1,mcrd         !Get location of current Gauss' point
               xi(i)=xigp(i,igauss)
            end do
            call kshp(RXi,coords,xi,mcrd,nnode,iwr)
            R(igauss)=R(igauss)+RXi
         end do

         do igauss=1,ngp        !Loop over Gauss points
            ivarm = 0
            do i=1,mcrd         !Get location of current Gauss' point
               xi(i)=xigp(i,igauss)
            end do
c     
c     Compute at current Gauss point, gradient matrix (gmtx)
c     and jacobians (xjac,xjacinv)
c     
            if (imat.gt.0) then
         
               call kshpquadpressNI(gmtx,xjac,xjacinv,detjac,coords,xi,
     +              R(igauss),ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr) !
            
c
c     Compute deformation gradient F at integration points
c
               call kfbgradEpressNI(F,E,press,u,gmtx,ndf,mcrd,ndofel,nsize,kinc,
     +         kstep,props,nprops)

c     
c     Evaluate force vector xfe & stifness matrix xke for the current element
c     
               call kumat12NI(xfe,xke,AI,F,E,press,props,imat,ndf,mcrd,
     +              ndofel,nsize,nprops,ninv,iwr,lflags(3),ndpress)
            else
               call kshpquad(gmtx,xjac,xjacinv,detjac,coords,xi,
     +              R(igauss),ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr) !

               call kfbgradE(F,E,u,gmtx,ndf,mcrd,ndofel,nsize)
               call kumat12(xfe,xke,AI,F,E,press,props,imat,ndf,mcrd,
     +              ndofel,nsize,nprops,ninv,iwr,lflags(3),ndpress)
            end if

c     
c     Computation of Residual = [xfe] . [G]
c     
            call kinitia(pgauss,ndofel)
            call ktranspose(gmtx,gmtxt,nsize+1,ndofel)
            call kmult(gmtxt,xfe,pgauss,ndofel,nsize+1,1)
c
c     Computation of Element Stiffness Matrix = [G]^T . xke . [G]
c
            if (lflags(3).eq.1) then
               call kinitia(auxm1,ndofel*(nsize+1))
               call ktranspose(gmtx,gmtxt,nsize+1,ndofel)
               call kmult(gmtxt,xke,auxm1,ndofel,nsize+1,nsize+1)
               call kmult(auxm1,gmtx,sgauss,ndofel,nsize+1,ndofel)
            end if
c     
c     Gauss integration (sum) of residual & stiffness for the element
c     
            volume=wgp(igauss)*detjac*2*3.1428571428571430d0*R(igauss)
            do i=1,ndofel
               rhs(i,1)=rhs(i,1)-pgauss(i)*volume
            end do
            if (lflags(3).eq.1) then
               do i=1,ndofel
                  do j=1,ndofel
                     amatrx(i,j)=amatrx(i,j)+sgauss(i,j)*volume
                  end do
               end do
            end if
c
c     Gauss integration (sum) of strain energy for the element
            
!             if (imat.eq.1) then
! c
!                call kwisoNI(W,AI,press,props,ninv,nprops)
!             else 
!                call kwiso(W,AI,props,ninv,nprops)
!             end if
            call kwisoNI(W,AI,press,props,ninv,nprops,time(1))

            energy(2)=energy(2)+W*volume
c            call ktranspose(F,Ft,mcrd+1,mcrd+1)
c            call kmult(F,Ft,PV,mcrd+1,mcrd+1,mcrd+1)  
c            call cauchystr(Ft,xfe(1),xfe(2),
c     +           xfe(3),xfe(4),xfe(5),AI(3),T)
c     
c     UVARM output
c
            globalSdv(jelem,igauss,1)=AI(1)
            globalSdv(jelem,igauss,2)=AI(2)
            globalSdv(jelem,igauss,3)=AI(4)
            globalSdv(jelem,igauss,4)=AI(5)
		globalSdv(jelem,igauss,6)=F(1,1)
		globalSdv(jelem,igauss,7)=F(1,2)
		globalSdv(jelem,igauss,8)=F(2,1)
		globalSdv(jelem,igauss,9)=F(2,2)
		globalSdv(jelem,igauss,10)=F(3,3)
		globalSdv(jelem,igauss,11)=AI(3)    !xfe(7)/props(3)
		globalSdv(jelem,igauss,12)=press   !p
         end do !end Loop over Gauss points
      end do !end loop for reduced integration
c     *****************************************************************************  
c     Surface Tension Calculations
      if(ndload .eq. 1.0) then
            if(ndload .gt. 1.0) then
               write(*,*)'TERMINATE and REMESH'
               write(*,*)'Element Number=',jelem,'has more than 1 bdry faces'
            end if
            gam_st = props(13)
            ! if (kstep.eq.1) then
            !       gam_st = adlmag(1,1)
            !   else
            !       gam_st = props(13)
            ! end if
            ! write(7,*) 'gamma_st=', gam_st,'time',time(1),'step=',kstep
            if (jdltyp(1,1) .eq. 3) then  ! U3 Face 1-2
                  call ksurgausspoints(wsgp,xisgp) 
                  do isgauss=1,5              ! Start Surface Gauss Loop
                        xi(1) = xisgp(isgauss)
                        xi(2) = 1.0d0 - xi(1)
                        call kshp(RXis,coords,xi,mcrd,nnode,iwr)
                        R_s(isgauss) = R_s(isgauss) + RXis
                  end do
                  do isgauss=1,5              
                        xi(1) = xisgp(isgauss)
                        xi(2) = 1.0d0 - xi(1)

                        call kshpquadpressNI(gmtx,xjac,xjacinv,detjac,coords,xi,
     +                     R_s(isgauss),ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr)
                        call ksurfF(xi,Fsurf,refn,Po,gmtxs,F,gmtx,u,ndf,mcrd,coords
     +                             ,ndofel,nsize,props,nprops,nnode)
                        call ksurnormal(currn,P,refn,F,mcrd)
                        call kumatsurften(xkesurften,xfesurften,F,Fsurf,refn,Po,
     +                              currn,P,gam_st,nsize,mcrd,Jhat)
                        call ktranspose(gmtxs,gmtxst,nsize+1,ndofel) 
                        call kmult(gmtxst,xkesurften,auxm2,ndofel,nsize+1,nsize+1)
                        call kmult(auxm2,gmtxs,ssgauss,ndofel,nsize+1,ndofel)
                        call kmult(gmtxst,xfesurften,psgauss,ndofel,nsize+1,1)
                        call ksurfjac(detjacsurf,coords,
     +                               xi,ndf,mcrd,nnode,ndofel,nsize,iwr,1)

                        ! if (W .lt. 0.0d0) then
                        ! write(7,*) 'Jhat=',Jhat,'jelem=',jelem,'kinc=',kinc
                        ! end if 
                        sarea = 2*3.1428571428571430d0*R_s(isgauss)*wsgp(isgauss)*detjacsurf      
                          
                        
                        do icount=1,ndofel
                              rhs(icount,1)=rhs(icount,1)
     +                     -psgauss(icount)*sarea
                        end do  ! Final RHS vector
          
                            do icount=1,ndofel
                              do jcount=1,ndofel
                                  amatrx(icount,jcount)=
     +                            amatrx(icount,jcount)
     +                            +ssgauss(icount,jcount)*sarea
                              end do
                            end do  ! Final amatrx
                            call kwisosurf(Wsurf,Jhat,gam_st)
                            energy(2)=energy(2)+Wsurf*sarea
                        
                  end do                 

            else if (jdltyp(1,1) .eq. 4) then  ! U4 Face 3-1
                  call ksurgausspoints(wsgp,xisgp) 
                  do isgauss=1,5              ! Start Surface Gauss Loop
                        xi(1) = xisgp(isgauss)
                        xi(2) = 0.0d0
                        call kshp(RXis,coords,xi,mcrd,nnode,iwr)
                        R_s(isgauss) = R_s(isgauss) + RXis
                  end do
                  do isgauss=1,5              
                        xi(1) = xisgp(isgauss)
                        xi(2) = 0.0d0

                        call kshpquadpressNI(gmtx,xjac,xjacinv,detjac,coords,xi,
     +                     R_s(isgauss),ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr)
                        call ksurfF(xi,Fsurf,refn,Po,gmtxs,F,gmtx,u,ndf,mcrd,coords
     +                             ,ndofel,nsize,props,nprops,nnode)
                        call ksurnormal(currn,P,refn,F,mcrd)
                        call kumatsurften(xkesurften,xfesurften,F,Fsurf,refn,Po,
     +                              currn,P,gam_st,nsize,mcrd,Jhat)
                        call ktranspose(gmtxs,gmtxst,nsize+1,ndofel) 
                        call kmult(gmtxst,xkesurften,auxm2,ndofel,nsize+1,nsize+1)
                        call kmult(auxm2,gmtxs,ssgauss,ndofel,nsize+1,ndofel)
                        call kmult(gmtxst,xfesurften,psgauss,ndofel,nsize+1,1)
                        call ksurfjac(detjacsurf,coords,
     +                               xi,ndf,mcrd,nnode,ndofel,nsize,iwr,2)

                        sarea = 2*3.1428571428571430d0*R_s(isgauss)*wsgp(isgauss)*detjacsurf      
                          
                        
                        do icount=1,ndofel
                              rhs(icount,1)=rhs(icount,1)
     +                     -psgauss(icount)*sarea
                            end do  ! Final RHS vector
          
                            do icount=1,ndofel
                              do jcount=1,ndofel
                                  amatrx(icount,jcount)=
     +                      amatrx(icount,jcount)
     +                      +ssgauss(icount,jcount)*sarea
                              end do
                            end do  ! Final amatrx     
                            call kwisosurf(Wsurf,Jhat,gam_st)
                            energy(2)=energy(2)+Wsurf*sarea                  
                        
                  end do

            else if (jdltyp(1,1) .eq. 5) then  ! U5 Face 2-3
                  call ksurgausspoints(wsgp,xisgp) 
                  do isgauss=1,5              ! Start Surface Gauss Loop
                        xi(1) = 0.0d0
                        xi(2) = xisgp(isgauss)
                        call kshp(RXis,coords,xi,mcrd,nnode,iwr)
                        R_s(isgauss) = R_s(isgauss) + RXis
                  end do
                  do isgauss=1,5              
                        xi(1) = 0.0d0
                        xi(2) = xisgp(isgauss)

                        call kshpquadpressNI(gmtx,xjac,xjacinv,detjac,coords,xi,
     +                     R_s(isgauss),ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr)
                        call ksurfF(xi,Fsurf,refn,Po,gmtxs,F,gmtx,u,ndf,mcrd,coords
     +                             ,ndofel,nsize,props,nprops,nnode)
                        call ksurnormal(currn,P,refn,F,mcrd)
                        call kumatsurften(xkesurften,xfesurften,F,Fsurf,refn,Po,
     +                              currn,P,gam_st,nsize,mcrd,Jhat)
                        call ktranspose(gmtxs,gmtxst,nsize+1,ndofel) 
                        call kmult(gmtxst,xkesurften,auxm2,ndofel,nsize+1,nsize+1)
                        call kmult(auxm2,gmtxs,ssgauss,ndofel,nsize+1,ndofel)
                        call kmult(gmtxst,xfesurften,psgauss,ndofel,nsize+1,1)
                        call ksurfjac(detjacsurf,coords,
     +                               xi,ndf,mcrd,nnode,ndofel,nsize,iwr,3)

                        sarea = 2*3.1428571428571430d0*R_s(isgauss)*wsgp(isgauss)*detjacsurf      
                          
                        
                        do icount=1,ndofel
                              rhs(icount,1)=rhs(icount,1)
     +                     -psgauss(icount)*sarea
                            end do  ! Final RHS vector
          
                            do icount=1,ndofel
                              do jcount=1,ndofel
                                  amatrx(icount,jcount)=
     +                      amatrx(icount,jcount)
     +                      +ssgauss(icount,jcount)*sarea
                              end do
                            end do  ! Final amatrx   
                            call kwisosurf(Wsurf,Jhat,gam_st)
                            energy(2)=energy(2)+Wsurf*sarea                   
                        
                  end do
            end if
      end if
      do igauss=1,ngp
            globalSdv(jelem,igauss,13) = energy(2)
      end do
c     *****************************************************************************  
c
c     Written Output
c
      if (iwr.ge.1) then
            write(7,*) '-------------------'
            write(7,*) 'element number=',jelem
            write(7,*) 'time=',time
            write(7,*) 'ndf,mcrd,nnode,ndofel'
            write(7,2001) ndf,mcrd,nnode,ndofel
            write(7,*) 'nmatp,ninv,imat,iwr'
            write(7,2001) nprops,ninv,imat,iwr
            write(7,*) 'G, K, eps'
            write(7,2003) (props(i),i=1,nprops)
            write(7,*) 'imat,ngaus,iwr'
            write(7,2001) imat,ngp,iwr
            write(7,*) 'xl'
            do i=1,nnode
               write(7,2002) (coords(j,i),j=1,mcrd)
            enddo
            write(7,*) 'volume, detjac'
            write(7,2002) volume,detjac
            write(7,*) 'gauss points'
            write(7,*) 'weight'
            write(7,2002) (wgp(i),i=1,ngp)
            write(7,*) 'xi'
            write(7,2002) (xi(i),i=1,mcrd)
            write(7,*) 'xigp'
            do i=1,ngp
               write(7,2002) (xigp(j,i),j=1,3)
            end do
            write(7,*) 'u'
            write(7,2002) (u(j),j=1,ndofel)
c            write(7,*) 'gmtx'
c            do i=1,nsize
c               write(7,2002) (gmtx(i,j),j=1,ndofel)
c            enddo
            write(7,*) 'F'
            do i=1,mcrd
               write(7,2002) (F(i,j),j=1,mcrd)
            enddo
            write(7,*) 'E'
            do i=1,mcrd
               write(7,2002) E(i)
            enddo
            write(7,*) 'pressure'
            write(7,*) press            
            write(7,*) 'xfe'
            write(7,2002) (xfe(j),j=1,nsize)
c            write(7,*) 'xke - Stiffness matrix'
c            do i=1,nsize
c               write(7,2003) (xke(i,j),j=1,nsize)
c            end do
c            write(7,*) 'Residual p'
c            do i=1,ndofel
c               write(7,2002) rhs(i,1)
c            end do
c            write(7,*) 'Stiffness K'
c            do i=1,ndofel
c               write(7,2003) (amatrx(i,j),j=1,ndofel)
c            end do
            write(7,*) 'Energy W'
            write(7,*) energy(2)
            write(7,*) '-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'
            write(7,*)
            write(7,*)
            write(7,*)
            write(7,*)
      
        end if
      end if
      end if

      return
 2001 format(500(I6))
 2002 format(500(F14.8))
 2003 format(500(E13.5))
      end

c***********************************************************************
c******************************subroutines******************************
c***********************************************************************
c**********************************************************************
      subroutine uvarm(uvar,direct,t,time,dtime,cmname,orname,
     1 nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     2 jmac,jmatyp,matlayo,laccfla)

      ! this subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.

      use global

      include 'aba_param.inc'

      character*80 cmname,orname
      character*3 flgray(15)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
      do i=1,nUvrm
         uvar(i) = globalSdv(noel-ElemOffset,npt,i)
      end do

c      for example
c      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
c      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
c      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      return
      end subroutine uvarm

****************************************************************************
c***********************************************************************
c
c     Energy function as in Notes_xfeM.pdf
c
      subroutine kumat12(xfe,xke,AI,F,E,press,props,imat,ndf,mcrd,
     + ndofel,nsize,nprops,ninv,iwr,irsflag,ndpress)
c
      implicit double precision(a-h,o-z)
c
      dimension AI(ninv)
      dimension xfe(nsize+1),
     +  xke(nsize+1,nsize+1),delta(mcrd+1,mcrd+1)
      dimension F(mcrd+1,mcrd+1),Finv(mcrd+1,mcrd+1),Cg(mcrd+1,mcrd+1)
      dimension Ft(mcrd+1,mcrd+1)
      dimension props(nprops)
      dimension DWDI(ninv),DWDI2(ninv,ninv)
      dimension DIDF(ninv,mcrd+1,mcrd+1),DIDF2(ninv,mcrd+1,mcrd+1,mcrd+1
     +                  ,mcrd+1)
      dimension DWDF(mcrd+1,mcrd+1),DWDF2(mcrd+1,mcrd+1,mcrd+1,mcrd+1)
      dimension DPDF(mcrd+1,mcrd+1),DPDF2(mcrd+1,mcrd+1,mcrd+1,mcrd+1)
      dimension E(mcrd+1)
      dimension DIDE(ninv,mcrd+1),DIDE2(ninv,mcrd+1,mcrd+1)
      dimension DIDFE(ninv,mcrd+1,mcrd+1,mcrd+1)
      dimension DWDE(mcrd+1),DWDE2(mcrd+1,mcrd+1),DWDFE(mcrd+1,mcrd+1,mcrd+1)
      dimension DPDE(mcrd+1),DPDE2(mcrd+1,mcrd+1),DPDFE(mcrd+1,mcrd+1,mcrd+1)
c      dimension propsnew(nprops)
c
c.... Evaluate Kronecker delta
      call kidtens(delta,mcrd+1)
c
c     Evaluate invariants, and derivatives of the energy function
c
      
      call kinvariso(AI,F,E,Cg,mcrd,ninv)
c      call NRnupsolve(AI(4), AI(5),props,nprops,pnup)
c      call effprops(props,nprops,pnup,propsnew)
c         write(7,*),'----------------------------------------'
c         write(7,*),'I4=',AI(4),'I5=',AI(5)
c         write(7,*),'nup=',pnup
c         write(7,*),'eff shear mod=',propsnew(1)/10.0d0**6.0d0,'lambda=',propsnew(2),'AEPS=',propsnew(3),'AMKE=',propsnew(4)
c         write(7,*),'----------------------------------------'
      call kdwdiiso(DWDI,AI,props,ninv,nprops)
      call kdwdi2iso(DWDI2,AI,props,ninv,nprops)
      call kdidfdfdeiso(DIDF,DIDF2,DIDE,DIDE2,DIDFE,
     +                     F,E,detF,ninv,mcrd)

c
c     Evaluate 1st & 2nd order derivatives
c
c     DWDF=dW/dFij
c     DWDF2=d2W/dFijdFkl
c     DWDE=dW/dEi
c     DWDE2=d2W/dEidEj
c     DWDFE=d2W/DFijdEk
c     DWDp=dW/dp
c     DWDp2=d2W/dp2
c     DWDFp=d2W/dFijdp
c     DWDEp=d2W/dEidp
c
c
      call kinitia(DWDF,(mcrd+1)*(mcrd+1))
      call kinitia(DWDF2,(mcrd+1)*(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DWDE,(mcrd+1))
      call kinitia(DWDE2,(mcrd+1)*(mcrd+1))
      call kinitia(DWDFE,(mcrd+1)*(mcrd+1)*(mcrd+1))
c
      do 10 i=1,mcrd+1
      do 10 j=1,mcrd+1
      do 10 k=1,mcrd+1
      do 10 l=1,mcrd+1
       do 20 ip=1,ninv
        if((j+k+l).eq.3) DWDE(i)=DWDE(i)+DWDI(ip)*DIDE(ip,i)
        if((k+l).eq.2) DWDF(i,j)=DWDF(i,j)+DWDI(ip)*DIDF(ip,i,j)
        if((k+l).eq.2) DWDE2(i,j)=DWDE2(i,j)+DWDI(ip)*DIDE2(ip,i,j)
        if(l.eq.1) DWDFE(i,j,k)=DWDFE(i,j,k)+DWDI(ip)*DIDFE(ip,i,j,k)
c.......
        DWDF2(i,j,k,l)=DWDF2(i,j,k,l)+DWDI(ip)*DIDF2(ip,i,j,k,l)

        do 30 iq=1,ninv          
c......... 
        if((k+l).eq.2) DWDE2(i,j)=DWDE2(i,j)+DWDI2(ip,iq)*DIDE(ip,i)*
     +                                       DIDE(iq,j)   
        if(l.eq.1) DWDFE(i,j,k)=DWDFE(i,j,k)+DWDI2(ip,iq)*DIDF(ip,i,j)*
     +                                       DIDE(iq,k)
        DWDF2(i,j,k,l)=DWDF2(i,j,k,l)+DWDI2(ip,iq)*DIDF(ip,i,j)*
     +          DIDF(iq,k,l)
 30     continue
 20    continue
 10   continue
c
c     Evaluate 1st & 2nd order derivatives
c     of potential energy P
c
c     DPDF=dP/dFij
c     DPDp=dP/dp  pressure
c     DPDF2=d2P/dFijdFkl
c     DPDE=dP/dEi
c     DPDE2=d2P/dEidEj
c     DPDFE=d2P/DFijdEk  
c     DPDFp=d2P/dFijdp
c     DPDEp=d2P/dEidp
c     DPDp2=d2P/dp2
c
c
      call kinitia(DPDF,(mcrd+1)*(mcrd+1))
      call kinitia(DPDF2,(mcrd+1)*(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DPDE,mcrd+1)
      call kinitia(DPDE2,(mcrd+1)*(mcrd+1))
      call kinitia(DPDFE,(mcrd+1)*(mcrd+1)*(mcrd+1))
      
c
c     Subsidiary calculations for later use in the
c     computation of the derivatives of P
c

c...  Compute inverse of F and auxiliary scalars
      call kinverse(F,Finv,detF,mcrd+1)

c.....
      do 11 i=1,mcrd+1
         DPDE(i)=DWDE(i)
c........
         do 21 j=1,mcrd+1
            DPDE2(i,j)=DWDE2(i,j)
            DPDF(i,j)=DWDF(i,j)
c...........
            do 31 k=1,mcrd+1
               DPDFE(i,j,k)=DWDFE(i,j,k)
c..............
               do 41 l=1,mcrd+1
                  DPDF2(i,j,k,l)=DWDF2(i,j,k,l)
 41            continue
c..............
 31         continue
c...........
 21      continue
c........
 11   continue
c.....
c
c
c
c     Evaluation of the force vector xfe
c-----------------------------------------
c
c     xfe={dP/dFij,dPdBi}
c
      call kinitia(xfe,nsize)
      xfe(1)=DPDF(1,1)
      xfe(2)=DPDF(1,2)
      xfe(3)=DPDF(2,1)
      xfe(4)=DPDF(2,2)
      xfe(5)=DPDF(3,3)
      xfe(6)=DPDE(1)
      xfe(7)=DPDE(2)
c
c
c
c     Evaluation of the stiffness matrix xke
c-------------------------------------------
c
c     xke={DPDF2...} - {DP2FF...}
c
      call kinitia(xke,(nsize+1)*(nsize+1))
c      if (irsflag.eq.1) then
         xke(1,1)=DPDF2(1,1,1,1)
         xke(1,2)=DPDF2(1,1,1,2)
         xke(1,3)=DPDF2(1,1,2,1)
         xke(1,4)=DPDF2(1,1,2,2)
         xke(1,5)=DPDF2(1,1,3,3)
         xke(1,6)=DPDFE(1,1,1)
         xke(1,7)=DPDFE(1,1,2)

         xke(2,1)=xke(1,2)
         xke(2,2)=DPDF2(1,2,1,2)
         xke(2,3)=DPDF2(1,2,2,1)
         xke(2,4)=DPDF2(1,2,2,2)
         xke(2,5)=DPDF2(1,2,3,3)
         xke(2,6)=DPDFE(1,2,1)
         xke(2,7)=DPDFE(1,2,2)

         xke(3,1)=xke(1,3)
         xke(3,2)=xke(2,3)
         xke(3,3)=DPDF2(2,1,2,1)
         xke(3,4)=DPDF2(2,1,2,2)
         xke(3,5)=DPDF2(2,1,3,3)
         xke(3,6)=DPDFE(2,1,1)
         xke(3,7)=DPDFE(2,1,2)

         xke(4,1)=xke(1,4)
         xke(4,2)=xke(2,4)
         xke(4,3)=xke(3,4)
         xke(4,4)=DPDF2(2,2,2,2)
         xke(4,5)=DPDF2(2,2,3,3)
         xke(4,6)=DPDFE(2,2,1)
         xke(4,7)=DPDFE(2,2,2)

         xke(5,1)=xke(1,5)
         xke(5,2)=xke(2,5)
         xke(5,3)=xke(3,5)
         xke(5,4)=xke(4,5)
         xke(5,5)=DPDF2(3,3,3,3)
         xke(5,6)=DPDFE(3,3,1)
         xke(5,7)=DPDFE(3,3,2)

         xke(6,1)=xke(1,6)
         xke(6,2)=xke(2,6)
         xke(6,3)=xke(3,6)
         xke(6,4)=xke(4,6)
         xke(6,5)=xke(5,6)
         xke(6,6)=DPDE2(1,1)
         xke(6,7)=DPDE2(1,2)

         xke(nsize+1,1)=xke(1,nsize+1)
         xke(nsize+1,2)=xke(2,nsize+1)
         xke(nsize+1,3)=xke(3,nsize+1)
         xke(nsize+1,4)=xke(4,nsize+1)
         xke(nsize+1,5)=xke(5,nsize+1)
         xke(nsize+1,6)=xke(6,nsize+1)
         xke(nsize+1,7)=DPDE2(2,2)

c      end if



c
c     Write in 'debug.dat'
c-------------------------------------------
c

      if (iwr.gt.1) then
         write(7,*) 'Invariants'
         write(7,*) (AI(i),i=1,ninv)
c         write(7,*) 'delta'
c         do i=1,mcrd
c            write(7,*) (delta(i,j),j=1,mcrd)
c         enddo
         write(7,*) 'F'
         do i=1,mcrd
            write(7,*) (F(i,j),j=1,mcrd)
         enddo
         write(7,*) 'E'
         write(7,*) (E(i),i=1,mcrd)
c         write(7,*) 'Finv'
c         do i=1,mcrd
c            write(7,*) (Finv(i,j),j=1,mcrd)
c         enddo
c         write(7,*) 'Cg'
c         do i=1,mcrd
c            write(7,*) (Cg(i,j),j=1,mcrd)
c         enddo
c         write(7,*) 'DWDF'
c         do i=1,mcrd
c            write(7,*) (DWDF(i,j),j=1,mcrd)
c         enddo
c         write(7,*) 'DPDF'
c         do i=1,mcrd
c            write(7,*) (DPDF(i,j),j=1,mcrd)
c         enddo
         write(7,*) 'DWDI'
         write(7,*) (DWDI(j),j=1,ninv)
c         write(7,*) 'DIDF'
c         do i=1,ninv
c            write(7,*) '..........................................'
c            do j=1,mcrd
c               write(7,*) (DIDF(i,j,k),k=1,mcrd)
c            end do   
c         end do
c         
c
c
c         write(7,*) 'DWDF2'
c         do i=1,mcrd
c            write(7,*) '..........................................'
c            do j=1,mcrd
c               write(7,*) '..........................................'
c               do k=1,mcrd
c                  write(7,*) (DWDF2(i,j,k,l),l=1,mcrd)
c               enddo
c            enddo
c         enddo
c
c         write(7,*) 'DPDF2'
c         do i=1,mcrd
c            write(7,*) '..........................................'
c            do j=1,mcrd
c               write(7,*) '..........................................'
c               do k=1,mcrd
c                  write(7,*) (DPDF2(i,j,k,l),l=1,mcrd)
c               enddo
c            enddo
c         enddo
c
c         write(7,*) '-------------------------------------------------'

         write(7,*) 'xfe'
         write(7,*) (xfe(j),j=1,nsize)
         write(7,*) 'xke'
         do i=1,nsize
           write(7,*) (xke(i,j),j=1,nsize)
         enddo
         write(7,*) '-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'
         write(7,*)
         write(7,*)
      endif

c      close(7)

      return
      end
c
c**********************************************************************
c***********************************************************************
c
c     Energy function as in Notes_xfeM.pdf
c
      subroutine kumat12NI(xfe,xke,AI,F,E,press,props,imat,ndf,mcrd,
     + ndofel,nsize,nprops,ninv,iwr,irsflag,ndpress)
c
      implicit double precision(a-h,o-z)
c
      dimension AI(ninv)
      dimension xfe(nsize+1),
     +  xke(nsize+1,nsize+1),delta(mcrd+1,mcrd+1)
      dimension F(mcrd+1,mcrd+1),Finv(mcrd+1,mcrd+1),Cg(mcrd+1,mcrd+1)
      dimension Ft(mcrd+1,mcrd+1)
      dimension props(nprops)
      dimension DWDI(ninv),DWDI2(ninv,ninv)
      dimension DIDF(ninv,mcrd+1,mcrd+1),DIDF2(ninv,mcrd+1,mcrd+1,mcrd+1
     +              ,mcrd+1)
      dimension DWDF(mcrd+1,mcrd+1),DWDF2(mcrd+1,mcrd+1,mcrd+1,mcrd+1)
      dimension DWDFp(mcrd+1,mcrd+1),DWDIp(ninv)
      dimension DPDF(mcrd+1,mcrd+1),DPDF2(mcrd+1,mcrd+1,mcrd+1,mcrd+1)
      dimension DPDFp(mcrd+1,mcrd+1)
      dimension E(mcrd+1)
      dimension DIDE(ninv,mcrd+1),DIDE2(ninv,mcrd+1,mcrd+1)
      dimension DIDFE(ninv,mcrd+1,mcrd+1,mcrd+1)
      dimension DWDE(mcrd+1),DWDE2(mcrd+1,mcrd+1),DWDFE(mcrd+1,mcrd+1,mcrd+1)
      dimension DPDE(mcrd+1),DPDE2(mcrd+1,mcrd+1),DPDFE(mcrd+1,mcrd+1,mcrd+1)
      dimension DWDEp(mcrd+1),DPDEp(mcrd+1)
      dimension propsnew(nprops)
c
c.... Evaluate Kronecker delta
      call kidtens(delta,mcrd+1)
c
c     Evaluate invariants, and derivatives of the energy function
c
      if (imat.gt.0) then
         call kinvariso(AI,F,E,Cg,mcrd,ninv)
         call kdwdiisoNI(DWDI,AI,press,props,ninv,nprops)
         call kdwdi2isoNI(DWDI2,AI,press,props,ninv,nprops)
         call kdwdpisoNI(DWDp,AI,press,props,ninv,nprops)
         call kdwdp2isoNI(DWDp2,AI,press,props,ninv,nprops)
         call kdwdipisoNI(DWDIp,AI,press,props,ninv,nprops)
         call kdidfdfdeiso(DIDF,DIDF2,DIDE,DIDE2,DIDFE,
     +                     F,E,detF,ninv,mcrd)
      endif
c
c     Evaluate 1st & 2nd order derivatives
c
c     DWDF=dW/dFij
c     DWDF2=d2W/dFijdFkl
c     DWDE=dW/dEi
c     DWDE2=d2W/dEidEj
c     DWDFE=d2W/DFijdEk
c     DWDp=dW/dp
c     DWDp2=d2W/dp2
c     DWDFp=d2W/dFijdp
c     DWDEp=d2W/dEidp
c
c
      call kinitia(DWDF,(mcrd+1)*(mcrd+1))
      call kinitia(DWDF2,(mcrd+1)*(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DWDFp,(mcrd+1)*(mcrd+1))
      call kinitia(DWDE,(mcrd+1))
      call kinitia(DWDE2,(mcrd+1)*(mcrd+1))
      call kinitia(DWDFE,(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DWDEp,mcrd+1)
c
      do 10 i=1,mcrd+1
      do 10 j=1,mcrd+1
      do 10 k=1,mcrd+1
      do 10 l=1,mcrd+1
       do 20 ip=1,ninv
        if((j+k+l).eq.3) DWDE(i)=DWDE(i)+DWDI(ip)*DIDE(ip,i)
        if((k+l).eq.2) DWDF(i,j)=DWDF(i,j)+DWDI(ip)*DIDF(ip,i,j)
        if((k+l).eq.2) DWDE2(i,j)=DWDE2(i,j)+DWDI(ip)*DIDE2(ip,i,j)
        if(l.eq.1) DWDFE(i,j,k)=DWDFE(i,j,k)+DWDI(ip)*DIDFE(ip,i,j,k)
        if((j+k+l).eq.3) DWDEp(i)=DWDEp(i)+DWDIp(ip)*DIDE(ip,i)
        if((k+l).eq.2) DWDFp(i,j)=DWDFp(i,j)+DWDIp(ip)*DIDF(ip,i,j)
c.......
        DWDF2(i,j,k,l)=DWDF2(i,j,k,l)+DWDI(ip)*DIDF2(ip,i,j,k,l)

        do 30 iq=1,ninv          
c......... 
        if((k+l).eq.2) DWDE2(i,j)=DWDE2(i,j)+DWDI2(ip,iq)*DIDE(ip,i)*
     +                                       DIDE(iq,j)   
        if(l.eq.1) DWDFE(i,j,k)=DWDFE(i,j,k)+DWDI2(ip,iq)*DIDF(ip,i,j)*
     +                                       DIDE(iq,k)
        DWDF2(i,j,k,l)=DWDF2(i,j,k,l)+DWDI2(ip,iq)*DIDF(ip,i,j)*
     +          DIDF(iq,k,l)
 30     continue
 20    continue
 10   continue
c
c     Evaluate 1st & 2nd order derivatives
c     of potential energy P
c
c     DPDF=dP/dFij
c     DPDp=dP/dp  pressure
c     DPDF2=d2P/dFijdFkl
c     DPDE=dP/dEi
c     DPDE2=d2P/dEidEj
c     DPDFE=d2P/DFijdEk  
c     DPDFp=d2P/dFijdp
c     DPDEp=d2P/dEidp
c     DPDp2=d2P/dp2
c
c
      call kinitia(DPDF,(mcrd+1)*(mcrd+1))
      call kinitia(DPDFp,(mcrd+1)*(mcrd+1))
      call kinitia(DPDF2,(mcrd+1)*(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DPDE,mcrd+1)
      call kinitia(DPDE2,(mcrd+1)*(mcrd+1))
      call kinitia(DPDFE,(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DPDEp,mcrd+1)
      
c
c     Subsidiary calculations for later use in the
c     computation of the derivatives of P
c

c...  Compute inverse of F and auxiliary scalars
      call kinverse(F,Finv,detF,mcrd+1)

c.....
      do 11 i=1,mcrd+1
         DPDE(i)=DWDE(i)
         DPDEp(i)=DWDEp(i)
c........
         do 21 j=1,mcrd+1
            DPDE2(i,j)=DWDE2(i,j)
            DPDF(i,j)=DWDF(i,j)
            DPDFp(i,j)=DWDFp(i,j)
c...........
            do 31 k=1,mcrd+1
               DPDFE(i,j,k)=DWDFE(i,j,k)
c..............
               do 41 l=1,mcrd+1
                  DPDF2(i,j,k,l)=DWDF2(i,j,k,l)
 41            continue
c..............
 31         continue
c...........
 21      continue
c........
 11   continue
c.....
      DPDp=DWDp
      DPDp2=DWDp2
c
c
c
c     Evaluation of the force vector xfe
c-----------------------------------------
c
c     xfe={dP/dFij,dPdBi}
c
      call kinitia(xfe,nsize)
      xfe(1)=DPDF(1,1)
      xfe(2)=DPDF(1,2)
      xfe(3)=DPDF(2,1)
      xfe(4)=DPDF(2,2)
      xfe(5)=DPDF(3,3)
      xfe(6)=DPDE(1)
      xfe(7)=DPDE(2)
      xfe(nsize+1)=DPDp
c
c
c
c     Evaluation of the stiffness matrix xke
c-------------------------------------------
c
c     xke={DPDF2...} - {DP2FF...}
c
      call kinitia(xke,(nsize+1)*(nsize+1))
c      if (irsflag.eq.1) then
         xke(1,1)=DPDF2(1,1,1,1)
         xke(1,2)=DPDF2(1,1,1,2)
         xke(1,3)=DPDF2(1,1,2,1)
         xke(1,4)=DPDF2(1,1,2,2)
         xke(1,5)=DPDF2(1,1,3,3)
         xke(1,6)=DPDFE(1,1,1)
         xke(1,7)=DPDFE(1,1,2)
         xke(1,nsize+1)=DPDFp(1,1)

         xke(2,1)=xke(1,2)
         xke(2,2)=DPDF2(1,2,1,2)
         xke(2,3)=DPDF2(1,2,2,1)
         xke(2,4)=DPDF2(1,2,2,2)
         xke(2,5)=DPDF2(1,2,3,3)
         xke(2,6)=DPDFE(1,2,1)
         xke(2,7)=DPDFE(1,2,2)
         xke(2,nsize+1)=DPDFp(1,2)

         xke(3,1)=xke(1,3)
         xke(3,2)=xke(2,3)
         xke(3,3)=DPDF2(2,1,2,1)
         xke(3,4)=DPDF2(2,1,2,2)
         xke(3,5)=DPDF2(2,1,3,3)
         xke(3,6)=DPDFE(2,1,1)
         xke(3,7)=DPDFE(2,1,2)
         xke(3,nsize+1)=DPDFp(2,1)

         xke(4,1)=xke(1,4)
         xke(4,2)=xke(2,4)
         xke(4,3)=xke(3,4)
         xke(4,4)=DPDF2(2,2,2,2)
         xke(4,5)=DPDF2(2,2,3,3)
         xke(4,6)=DPDFE(2,2,1)
         xke(4,7)=DPDFE(2,2,2)
         xke(4,nsize+1)=DPDFp(2,2)

         xke(5,1)=xke(1,5)
         xke(5,2)=xke(2,5)
         xke(5,3)=xke(3,5)
         xke(5,4)=xke(4,5)
         xke(5,5)=DPDF2(3,3,3,3)
         xke(5,6)=DPDFE(3,3,1)
         xke(5,7)=DPDFE(3,3,2)
         xke(5,nsize+1)=DPDFp(3,3)

         xke(6,1)=xke(1,6)
         xke(6,2)=xke(2,6)
         xke(6,3)=xke(3,6)
         xke(6,4)=xke(4,6)
         xke(6,5)=xke(5,6)
         xke(6,6)=DPDE2(1,1)
         xke(6,7)=DPDE2(1,2)
         xke(6,nsize+1)=DPDEp(1)

         xke(nsize,1)=xke(1,nsize)
         xke(nsize,2)=xke(2,nsize)
         xke(nsize,3)=xke(3,nsize)
         xke(nsize,4)=xke(4,nsize)
         xke(nsize,5)=xke(5,nsize)
         xke(nsize,6)=xke(6,nsize)
         xke(nsize,7)=DPDE2(2,2)
         xke(nsize,nsize+1)=DPDEp(2)

         xke(nsize+1,1)=xke(1,nsize+1)
         xke(nsize+1,2)=xke(2,nsize+1)
         xke(nsize+1,3)=xke(3,nsize+1)
         xke(nsize+1,4)=xke(4,nsize+1)
         xke(nsize+1,5)=xke(5,nsize+1)
         xke(nsize+1,6)=xke(6,nsize+1)
         xke(nsize+1,7)=xke(7,nsize+1)
         xke(nsize+1,nsize+1)=DPDp2

c      end if



c
c     Write in 'debug.dat'
c-------------------------------------------
c

      if (iwr.gt.1) then
         write(7,*) 'Invariants'
         write(7,*) (AI(i),i=1,ninv)
c         write(7,*) 'delta'
c         do i=1,mcrd
c            write(7,*) (delta(i,j),j=1,mcrd)
c         enddo
         write(7,*) 'F'
         do i=1,mcrd
            write(7,*) (F(i,j),j=1,mcrd)
         enddo
         write(7,*) 'E'
         write(7,*) (E(i),i=1,mcrd)
c         write(7,*) 'Finv'
c         do i=1,mcrd
c            write(7,*) (Finv(i,j),j=1,mcrd)
c         enddo
c         write(7,*) 'Cg'
c         do i=1,mcrd
c            write(7,*) (Cg(i,j),j=1,mcrd)
c         enddo
c         write(7,*) 'DWDF'
c         do i=1,mcrd
c            write(7,*) (DWDF(i,j),j=1,mcrd)
c         enddo
c         write(7,*) 'DPDF'
c         do i=1,mcrd
c            write(7,*) (DPDF(i,j),j=1,mcrd)
c         enddo
         write(7,*) 'DWDI'
         write(7,*) (DWDI(j),j=1,ninv)
c         write(7,*) 'DIDF'
c         do i=1,ninv
c            write(7,*) '..........................................'
c            do j=1,mcrd
c               write(7,*) (DIDF(i,j,k),k=1,mcrd)
c            end do   
c         end do
c         
c
         write(7,*) 'DWDI2'
         do i=1,ninv
            write(7,*) (DWDI2(i,j),j=1,ninv)
         end do
         write(7,*) DWDI2(4,1)
         write(7,*) 'DWDp'
         write(7,*) DWDp
         write(7,*) 'DWDp2'
         write(7,*) DWDp2
         write(7,*) 'DWDIp'
         write(7,*) (DWDIp(j),j=1,ninv)
c
c         write(7,*) 'DWDF2'
c         do i=1,mcrd
c            write(7,*) '..........................................'
c            do j=1,mcrd
c               write(7,*) '..........................................'
c               do k=1,mcrd
c                  write(7,*) (DWDF2(i,j,k,l),l=1,mcrd)
c               enddo
c            enddo
c         enddo
c
c         write(7,*) 'DPDF2'
c         do i=1,mcrd
c            write(7,*) '..........................................'
c            do j=1,mcrd
c               write(7,*) '..........................................'
c               do k=1,mcrd
c                  write(7,*) (DPDF2(i,j,k,l),l=1,mcrd)
c               enddo
c            enddo
c         enddo
c
c         write(7,*) '-------------------------------------------------'

         write(7,*) 'xfe'
         write(7,*) (xfe(j),j=1,nsize)
         write(7,*) 'xke'
         do i=1,nsize
           write(7,*) (xke(i,j),j=1,nsize)
         enddo
         write(7,*) '-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'
         write(7,*)
         write(7,*)
      endif

c      close(7)

      return
      end
c
c**********************************************************************
c**********************************************************************
c
      subroutine kwiso(W,AI,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension props(nprops)
c
c     Asign material parameters
c
      AG=props(1)  !Shear modulus
      AEPS=props(3) !Dielectric permittivity
      AMKE=props(4) !Dielectric permittivity
      
      W=AG/2.0d0*(AI(1)-3.0d0-2.0d0*DLog(AI(3)))+
     +        AMKE*AI(3)/2.0d0*AI(5)+(AMKE-AEPS)*AI(4)/2.0d0       

      return
      end
c
c**********************************************************************
c
      subroutine kwisoNI(W,AI,press,props,ninv,nprops,time)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension props(nprops)
      real*8 time
c
c     Asign material parameters
c
      AG1 = props(1)  !G1  
      AG2 = props(11) !G2
      AL = props(2)  !Lame modulus
      AEPS = props(3) !Dielectric permittivity
      AMKE = props(4) !Dielectric permittivity
      medium = props(7)
      AG = AG1+AG2    !G1+G2 
      ms = props(6)
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0
c
c     
      AJ = (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press +
     +  SQRT(16.0d0*AG*AL + 
     +  (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)

      if (ms .eq. 2.0d0) then
            AG = 0.0d0
      end if


      W=AG/2.0d0*(AI(1)-3.0d0-2.0d0*DLog(AJ))+
     +        AL/2.0d0*(AJ-1.0d0)**2.0d0+press*(AI(3)-1.0d0)-
     +        press*(AJ-1.0d0) + pgam*AJ

      if (ms .eq. 2.0d0) then
            if (time .eq. 0.0d0) then
                  W = pgam
            end if
      end if
            
      ! W = 0.0d0

      ! w = AG*(AI(3) - 3.0d0)/2.0d0
     

      return
      end
c
c**********************************************************************
      subroutine kwisosurf(W,Jhat,gam_st)
c
      implicit double precision(a-h,o-z)
      real*8 gam_st,Jhat
                
      W = gam_st*Jhat 
      ! W = 0.0
            
      if (W .lt. 0.0d0) then
            write(7,*) 'SURFACE ENERGY IS NEGATIVE'
      end if 
                  
c      write(7,*)'gamma=',gam_st
                  
                      
      return
      end
c**********************************************************************
c
      subroutine kdwdiiso(DWDI,AI,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension DWDI(ninv),props(nprops)
c
c     Asign material parameters
c
      AG=props(1)   !Shear modulus
      AL=props(2)  !Lame modulus
      AEPS=props(3) !Dielectric permittivity
      AMKE=props(4) !Dielectric permittivity
c
c     
      call kinitia(DWDI,ninv)
      DWDI(1)=AG/2.0d0
      DWDI(3)=-AG/AI(3)-AI(5)*AMKE/2.0d0+AL*(AI(3)-1)
      DWDI(4)=(AMKE-AEPS)/2.0d0
      DWDI(5)=-AMKE/2.0d0*AI(3)
      

      return
      end
c
c**********************************************************************
c
      subroutine kdwdiisoNI(DWDI,AI,press,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension DWDI(ninv),props(11)
      real*8 AG1,AG2,a1,a2,AG,AN,AIsph
      real*8 c,mu0,ms,mbeta,arg,coth
      real*8 cosech,medium
      mu0=1.256*10.0d0**(-6.0d0)
c
c     Asign material parameters
c
      AG1 = props(1)  !G1  
      AG2 = props(11) !G2
      AL = props(2)  !Lame modulus
      AEPS = props(3) !Dielectric permittivity
      AMKE = props(4) !Dielectric permittivity
      medium = props(7)
      AG = AG1+AG2    !G1+G2 
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0
c
c     
      AJ = (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press +
     +  SQRT(16.0d0*AG*AL + 
     +  (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)

      call kinitia(DWDI,ninv)

      DWDI(1) = AG/2.0d0
      DWDI(3) = press
      DWDI(4) = (AMKE-AEPS)/2.0d0
      DWDI(5) = -(AEPS*(2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press + 
     +         SQRT(16.0d0*AG*AL + 
     +     (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0)))/(8.0d0*AL)
     

      return
      end
c
c
c**********************************************************************
c**********************************************************************
c
      subroutine kdwdi2iso(DWDI2,AI,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension DWDI2(ninv,ninv),props(nprops)

c
c     Asign material parameters
c
      AG=props(1)   !Shear modulus
      AL=props(2)  !Lame modulus
      AEPS=props(3) !Dielectric permittivity
      AMKE=props(4) !Dielectric permittivity
c
      call kinitia(DWDI2,ninv*ninv)
      DWDI2(3,3)=AG/(AI(3)**2.0d0)+AL
      DWDI2(3,5)=-AMKE/2.0d0
      DWDI2(5,3)=-AMKE/2.0d0

c    
c-------------------------------------------------------    
      return
      end
c
c**********************************************************************
c
      subroutine kdwdi2isoNI(DWDI2,AI,press,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension DWDI2(ninv,ninv),props(11)
      real*8 AG1,AG2,a1,a2,AN,AG,AIsph
      real*8 c,mu0,ms,mbeta,arg,coth
      real*8 cosech,medium,cosh
      real*8 denom1,denom2,num1,num2
      real*8 num3,num4
      mu0=1.256*10.0d0**(-6.0d0)
c
c     Asign material parameters
c
      AG1 = props(1)  !G1  
      AG2 = props(11) !G2
      AL = props(2)  !Lame modulus
      AEPS = props(3) !Dielectric permittivity
      AMKE = props(4) !Dielectric permittivity
      medium = props(7)
      AG = AG1+AG2    !G1+G2 
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0
c
c     
      AJ = (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press +
     +  SQRT(16.0d0*AG*AL + 
     +  (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)

      call kinitia(DWDI2,ninv*ninv)

      DWDI2(5,5) = -(AEPS*(AEPS + (AEPS*(2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press))
     + /SQRT(16.0d0*AG*AL + (2.0d0*AL + 
     +  AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0)))/(8.0d0*AL)

    
c-------------------------------------------------------    
      return
      end
c
c**********************************************************************
c**********************************************************************
c
      subroutine kdwdpisoNI(DWDp,AI,press,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension props(11)
      real*8 AG1,AG2,a1,a2,AN
c
c     Asign material parameters
c
      AG1 = props(1)  !G1  
      AG2 = props(11) !G2
      AL = props(2)  !Lame modulus
      AEPS = props(3) !Dielectric permittivity
      AMKE = props(4) !Dielectric permittivity
      medium = props(7)
      AG = AG1+AG2    !G1+G2 
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0

c
c     
      AJ = (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press +
     +  SQRT(16.0d0*AG*AL + 
     +  (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)
c
      DWDp = -(AEPS*AI(5) + AL*(2.0d0 - 4.0d0*AI(3)) -2.0d0*pgam+ 
     +       2.0d0*press + SQRT(16.0d0*AG*AL + 
     +       (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)
c
      return
      end
c
c**********************************************************************
c**********************************************************************
c
      subroutine kdwdp2isoNI(DWDp2,AI,press,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension props(nprops)
      real*8 AG1,AG2,a1,a2,AG,AN
c
c     Asign material parameters
c
      AG1 = props(1)  !G1  
      AG2 = props(11) !G2
      AL = props(2)  !Lame modulus
      AEPS = props(3) !Dielectric permittivity
      AMKE = props(4) !Dielectric permittivity
      medium = props(7)
      AG = AG1+AG2    !G1+G2 
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0
c
c     
      AJ = (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press +
     +  SQRT(16.0d0*AG*AL + 
     +  (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)
c
      DWDp2 = -(2.0d0 + (2.0d0*(2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press))
     +  /SQRT(16.0d0*AG*AL + (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 
     +  2.0d0*press)**2.0d0))/(4.0d0*AL)
c    
      return
      end
c
c**********************************************************************
c**********************************************************************
c
c**********************************************************************
c
      subroutine kdwdipisoNI(DWDIp,AI,press,props,ninv,nprops)
c
      implicit double precision(a-h,o-z)
      dimension AI(ninv)
      dimension DWDIp(ninv),props(nprops)
      real*8 AG1,AG2,a1,a2,AG,AN
c
c     Asign material parameters
c
      AG1 = props(1)  !G1  
      AG2 = props(11) !G2
      AL = props(2)  !Lame modulus
      AEPS = props(3) !Dielectric permittivity
      AMKE = props(4) !Dielectric permittivity
      medium = props(7)
      AG = AG1+AG2    !G1+G2 
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0
c
c           
      AJ = (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press +
     +  SQRT(16.0d0*AG*AL + 
     +  (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press)**2.0d0))/(4.0d0*AL)

      call kinitia(DWDIp,ninv)
c    
      DWDIp(3) = 1.0d0
      DWDIp(5) = -(AEPS + (AEPS*(2*AL + AEPS*AI(5) -2.0d0*pgam+ 2.0d0*press))
     +    /SQRT(16.0d0*AG*AL + (2.0d0*AL + AEPS*AI(5) -2.0d0*pgam
     +    + 2.0d0*press)**2.0d0))/(4.0d0*AL)
c 
      return
      end
c
c**********************************************************************
c**********************************************************************
c**********************************************************************
c
      subroutine kinvariso(AI,F,E,Cg,mcrd,ninv)
c
      implicit double precision(a-h,o-z)
c      
      dimension AI(ninv)
      dimension F(mcrd+1,mcrd+1),Ft(mcrd+1,mcrd+1),Cg(mcrd+1,mcrd+1)
      dimension E(mcrd+1),FiTE(mcrd+1),Finv(mcrd+1,mcrd+1),FinvT(mcrd+1,
     +          mcrd+1)
c
c     Compute Isotropic invariants I1 - I6
c     (For the moment only isotropic models are included)
c
      call kinitia(AI,ninv)
c
c
c     Compute tensorial quantities Ft, C,
c
      call ktranspose(F,Ft,mcrd+1,mcrd+1)
      call kmult(Ft,F,Cg,mcrd+1,mcrd+1,mcrd+1)
      call kinverse(F,Finv,detF,mcrd+1)
      call ktranspose(Finv,FinvT,mcrd+1,mcrd+1)
      call kmult(FinvT,E,FiTE,mcrd+1,mcrd+1,1)
c
c  Mechanical 
c  invariants for isotropic materials
c--------------------------------------
c
c.... I3=det(F)
      AI(3)=detF

      do i=1,mcrd+1
c.... I1=Trace(C)
         AI(1)=AI(1)+Cg(i,i)
c.... I5=FiTE.FiTE  and I4=E.E
         AI(4)=AI(4)+E(i)**2.0d0
         AI(5)=AI(5)+FiTE(i)**2.0d0
      end do
c
      return
      end
c
c**********************************************************************
c
     
      subroutine kdidfdfdeiso(DIDF,DIDF2,DIDE,DIDE2,DIDFE,
     +                        F,E,detF,ninv,mcrd)
c
      implicit double precision(a-h,o-z)
      
      dimension F(mcrd+1,mcrd+1),Finv(mcrd+1,mcrd+1),Ft(mcrd+1,mcrd+1)
      dimension DIDF(ninv,mcrd+1,mcrd+1)
      dimension DIDF2(ninv,mcrd+1,mcrd+1,mcrd+1,mcrd+1)
      dimension delta(mcrd+1,mcrd+1)
      dimension E(mcrd+1),FiTE(mcrd+1),FinvT(mcrd+1,mcrd+1)
      dimension Cg(mcrd+1,mcrd+1), Cginv(mcrd+1,mcrd+1)
      dimension FiFiTE(mcrd+1)
      dimension DIDE(ninv,mcrd+1),DIDE2(ninv,mcrd+1,mcrd+1)
      dimension DIDFE(ninv,mcrd+1,mcrd+1,mcrd+1)
c
c
c    Initialization of vectors and matrices
c
      call kinitia(DIDF,ninv*(mcrd+1)*(mcrd+1))
      call kinitia(DIDF2,ninv*(mcrd+1)*(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(DIDE,ninv*(mcrd+1))
      call kinitia(DIDE2,ninv*(mcrd+1)*(mcrd+1))
      call kinitia(DIDFE,ninv*(mcrd+1)*(mcrd+1)*(mcrd+1))
      call kinitia(Cg,(mcrd+1)*(mcrd+1))
      call kinitia(Ft,(mcrd+1)*(mcrd+1))
c
      call kidtens(delta,mcrd+1)
      call kinverse(F,Finv,detF,mcrd+1)
      call ktranspose(Finv,FinvT,mcrd+1,mcrd+1)
c 
c--------------------------------------------------------
c    Computation of auxiliary vectors & tensors
c
      call ktranspose(F,Ft,mcrd+1,mcrd+1)     ! Ft = F^T  (F transpose)
      call kmult(Ft,F,Cg,mcrd+1,mcrd+1,mcrd+1) ! C = F^T . F   (Right Green tensor)
      call kinverse(Cg,Cginv,detC,mcrd+1)
      call kmult(FinvT,E,FiTE,mcrd+1,mcrd+1,1)
      call kmult(Finv,FiTE,FiFiTE,mcrd+1,mcrd+1,1)
c
c    Computation of 1st and 2nd derivatives of the invariants
c
      do 20 i=1,mcrd+1
         DIDE(4,i)=2.0d0*E(i)
         DIDE(5,i)=2.0d0*FiFiTE(i)
         do 21 j=1,mcrd+1
            DIDF(1,i,j)=2.d0*F(i,j)
            DIDF(3,i,j)=detF*Finv(j,i)
            DIDF(5,i,j)=-2.0d0*FiTE(i)*FiFiTE(j)
c
            DIDE2(4,i,j)=2.0d0*delta(i,j)
            DIDE2(5,i,j)=2.0d0*Cginv(i,j)
            do 22 k=1,mcrd+1               
               DIDFE(5,i,j,k)=-2.0d0*Cginv(j,k)*FiTE(i)
     +                        -2.0d0*Finv(k,i)*FiFiTE(j)
               do 23 l=1,mcrd+1
                  DIDF2(1,i,j,k,l)=2.d0*delta(i,k)*delta(j,l)
                  DIDF2(3,i,j,k,l)=detF*(Finv(j,i)*Finv(l,k)-Finv(j,k)*
     +                 Finv(l,i))
                  DIDF2(5,i,j,k,l)=2.0d0*Cginv(j,l)*FiTE(i)*FiTE(k)+
     +                             2.0d0*Finv(l,i)*FiTE(k)*FiFiTE(j)+
     +                             2.0d0*Finv(j,k)*FiTE(i)*FiFiTE(l)
 23            continue    !loop l
 22         continue    !loop k
 21      continue    !loop j)
 20   continue    !loop i)
c
c
      return
      end  
c
c***********************************************************************

c***********************************************************************
c
c     Compute R

      subroutine kshp(RXi,coords,xi,mcrd,nnode,iwr)

      implicit double precision(a-h,o-z)

      dimension  coords(mcrd,nnode)
      dimension xi(mcrd),xnshp(nnode)
      real*8 RXi
c     Shape Functions
      xnshp(1)=xi(1)*(-1.0+2.0*xi(1)-
     + 3.0*xi(2)*(-1.0+xi(2)+xi(1)))
      xnshp(2)=xi(2)*(-1.0+2.0*xi(2)-
     + 3.0*xi(1)*(-1.0+xi(2)+xi(1)))
      xnshp(3)=-((-1.0+xi(2)+xi(1))*(1.0-
     + 2.0*xi(1)+xi(2)*(-2.0+3.0*xi(1))))
      xnshp(4)=4.0*xi(2)*xi(1)*(-2.0+3.0*xi(2)+3.0*xi(1))
      xnshp(5)=4.0*xi(2)*(-1.0+xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xnshp(6)=4.0*(-1.0+3.0*xi(2))*xi(1)*(-1.0+xi(2)+xi(1))
      xnshp(7)=-27.0*xi(2)*xi(1)*(-1.0+xi(2)+xi(1))

      RXi=0.0d0
      
      do i=1,nnode
        RXi=RXi+xnshp(i)*coords(1,i)
      end do
c
      return
      end
c
c***********************************************************************
c***********************************************************************
c
      subroutine kshpquad(gmtx,xjac,xjacinv,detjac,coordsabq,xi,
     +     R,ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr)

       implicit double precision(a-h,o-z)

      dimension coordsabq(mcrd,nnode)
      real*8 R
      dimension coords(mcrd,nnode),xnshp(nnode)
      dimension xi(mcrd),xnshppress(ndpress)
      dimension xngshpquad(nnode,mcrd)
      dimension xjac(mcrd,mcrd),xjacinv(mcrd,mcrd)
      dimension gmtx(nsize+1,ndofel)
      dimension xjimtx(nsize+1,nsize+1),xngmtx(nsize+1,ndofel)
c
c     Extract coordinates of mesh nodes only
c
      do i=1,mcrd
         do j=1,7
            coords(i,j)=coordsabq(i,j) 
         end do
      end do

c     Shape Functions
      xnshp(1)=xi(1)*(-1.0+2.0*xi(1)-
     + 3.0*xi(2)*(-1.0+xi(2)+xi(1)))
      xnshp(2)=xi(2)*(-1.0+2.0*xi(2)-
     + 3.0*xi(1)*(-1.0+xi(2)+xi(1)))
      xnshp(3)=-((-1.0+xi(2)+xi(1))*(1.0-
     + 2.0*xi(1)+xi(2)*(-2.0+3.0*xi(1))))
      xnshp(4)=4.0*xi(2)*xi(1)*(-2.0+3.0*xi(2)+3.0*xi(1))
      xnshp(5)=4.0*xi(2)*(-1.0+xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xnshp(6)=4.0*(-1.0+3.0*xi(2))*xi(1)*(-1.0+xi(2)+xi(1))
      xnshp(7)=-27.0*xi(2)*xi(1)*(-1.0+xi(2)+xi(1))

c
c     Quadratic Shape functions derivatives
c
c     ngshp(j,i)=dN(j)/dX(i)
      xngshpquad(1,1)=-1.0-3.0*xi(2)**2.0+xi(2)*(3.0-6.0*xi(1))+4.0*xi(1)
      xngshpquad(1,2)=-3.0*xi(1)*(-1+2.0*xi(2)+xi(1))
      xngshpquad(2,1)=-3.0*xi(2)*(-1.0+xi(2)+2.0*xi(1))
      xngshpquad(2,2)=-1.0 + xi(2)*(4.0-6.0*xi(1))+3.0*xi(1)-3.0*xi(1)**2.0
      xngshpquad(3,1)=-3.0-3.0*xi(2)**2.0+xi(2)*(7.0-6.0*xi(1))+4.0*xi(1)
      xngshpquad(3,2)=-3.0+xi(2)*(4.0-6.0*xi(1))+7.0*xi(1)-3.0*xi(1)**2.0
      xngshpquad(4,1)=4.0*xi(2)*(-2.0+3.0*xi(2)+6.0*xi(1))
      xngshpquad(4,2)=4.0*xi(1)*(-2.0+6.0*xi(2)+3.0*xi(1))
      xngshpquad(5,1)=4.0*xi(2)*(-4.0+3.0*xi(2)+6.0*xi(1))
      xngshpquad(5,2)=4.0*(-1.0+2.0*xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xngshpquad(6,1)=4.0*(-1.0+3.0*xi(2))*(-1.0+xi(2)+2.0*xi(1))
      xngshpquad(6,2)=4.0*xi(1)*(-4.0+6.0*xi(2)+3.0*xi(1))
      xngshpquad(7,1)=-27.0*xi(2)*(-1+xi(2)+2.0*xi(1))
      xngshpquad(7,2)=-27.0*xi(1)*(-1+2.0*xi(2)+xi(1))
	  
c
c     xjac(i,j) = sum_I {ngshp(I,k)*X(l,I)}
c
      do k=1,mcrd
      do l=1,mcrd
         xjac(k,l)=0.d0
         do i=1,nnode
            xjac(k,l)=xjac(k,l)+xngshpquad(i,k)*coords(l,i)
         end do
      end do
      end do
      call kinverse(xjac,xjacinv,detjac,mcrd)
c
c     Form the Jacobian Matrix [JI] -> xjimtx
c
      call kinitia(xjimtx,(nsize+1)*(nsize+1))
c      do 10 k=1,ndf-1 !Magnetic case
      do 10 k=1,mcrd
      do 10 i=1,mcrd
      do 10 j=1,mcrd
         xjimtx(mcrd*(k-1)+i,mcrd*(k-1)+j)=xjacinv(i,j)
 10   continue
c
c     For the electric case the last 3x3 entries
c     of [JI] must be modified according to
c     E1=-dphi/dX1 E2=-dphi/dX2 E3=-dphi/dX3
c
      xjimtx(5,5)=1.0d0        ! Axisymmetric u/R
      xjimtx(6,6)=-xjacinv(1,1)
      xjimtx(6,7)=-xjacinv(1,2)
      xjimtx(7,6)=-xjacinv(2,1)
      xjimtx(7,7)=-xjacinv(2,2)
c
c
c     Form the NG matrix [NG] -> ngmtxc
c
c
c
c      
      call kinitia(xngmtx,(nsize+1)*ndofel)
	  do 11 k=1,ndf-1
	  do 11 i=1,mcrd
	  do 11 j=1,nnode
	    xngmtx(mcrd*(k-1)+i,ndf*(j-1)+k)=xngshpquad(j,i)
 11   continue
      do j=1,nnode
         xngmtx(5,ndf*(j-1)+1)=xnshp(j)/R
      end do
      do j=1,nnode
         xngmtx(6,ndf*(j-1)+3)=xngshpquad(j,1)
         xngmtx(7,ndf*(j-1)+3)=xngshpquad(j,2)
      end do
c
c
c     Compute [G]=[JI] x [NG] -> gmtx(i,j)=xjimtx(i,k)*xngmtx(k,j)
c
      call kinitia(gmtx,(nsize+1)*ndofel)
c
      call kmult(xjimtx,xngmtx,gmtx,nsize+1,nsize+1,ndofel)
      return
      end
c
c***********************************************************************
c
      subroutine kshpquadpressNI(gmtx,xjac,xjacinv,detjac,coordsabq,xi,
     +     R,ndf,mcrd,nnode,ndofel,ndpress,nsize,iwr)

       implicit double precision(a-h,o-z)

      dimension coordsabq(mcrd,nnode)
      real*8 R
      dimension coords(mcrd,nnode),xnshp(nnode)
      dimension xi(mcrd),xnshppress(ndpress)
      dimension xngshpquad(nnode,mcrd)
      dimension xjac(mcrd,mcrd),xjacinv(mcrd,mcrd)
      dimension gmtx(nsize+1,ndofel)
      dimension xjimtx(nsize+1,nsize+1),xngmtx(nsize+1,ndofel)
c
c     Extract coordinates of mesh nodes only
c
      do i=1,mcrd
         do j=1,7
            coords(i,j)=coordsabq(i,j) 
         end do
      end do

c     Shape Functions
      xnshp(1)=xi(1)*(-1.0+2.0*xi(1)-
     + 3.0*xi(2)*(-1.0+xi(2)+xi(1)))
      xnshp(2)=xi(2)*(-1.0+2.0*xi(2)-
     + 3.0*xi(1)*(-1.0+xi(2)+xi(1)))
      xnshp(3)=-((-1.0+xi(2)+xi(1))*(1.0-
     + 2.0*xi(1)+xi(2)*(-2.0+3.0*xi(1))))
      xnshp(4)=4.0*xi(2)*xi(1)*(-2.0+3.0*xi(2)+3.0*xi(1))
      xnshp(5)=4.0*xi(2)*(-1.0+xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xnshp(6)=4.0*(-1.0+3.0*xi(2))*xi(1)*(-1.0+xi(2)+xi(1))
      xnshp(7)=-27.0*xi(2)*xi(1)*(-1.0+xi(2)+xi(1))

c
c     Quadratic Shape functions derivatives
c
c     ngshp(j,i)=dN(j)/dX(i)
      xngshpquad(1,1)=-1.0-3.0*xi(2)**2.0+xi(2)*(3.0-6.0*xi(1))+4.0*xi(1)
      xngshpquad(1,2)=-3.0*xi(1)*(-1+2.0*xi(2)+xi(1))
      xngshpquad(2,1)=-3.0*xi(2)*(-1.0+xi(2)+2.0*xi(1))
      xngshpquad(2,2)=-1.0 + xi(2)*(4.0-6.0*xi(1))+3.0*xi(1)-3.0*xi(1)**2.0
      xngshpquad(3,1)=-3.0-3.0*xi(2)**2.0+xi(2)*(7.0-6.0*xi(1))+4.0*xi(1)
      xngshpquad(3,2)=-3.0+xi(2)*(4.0-6.0*xi(1))+7.0*xi(1)-3.0*xi(1)**2.0
      xngshpquad(4,1)=4.0*xi(2)*(-2.0+3.0*xi(2)+6.0*xi(1))
      xngshpquad(4,2)=4.0*xi(1)*(-2.0+6.0*xi(2)+3.0*xi(1))
      xngshpquad(5,1)=4.0*xi(2)*(-4.0+3.0*xi(2)+6.0*xi(1))
      xngshpquad(5,2)=4.0*(-1.0+2.0*xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xngshpquad(6,1)=4.0*(-1.0+3.0*xi(2))*(-1.0+xi(2)+2.0*xi(1))
      xngshpquad(6,2)=4.0*xi(1)*(-4.0+6.0*xi(2)+3.0*xi(1))
      xngshpquad(7,1)=-27.0*xi(2)*(-1+xi(2)+2.0*xi(1))
      xngshpquad(7,2)=-27.0*xi(1)*(-1+2.0*xi(2)+xi(1))
	  
c
c     Pressure shape functions
c
      xnshppress(1)=1.0
	xnshppress(2)=xi(1)
	xnshppress(3)=xi(2)
c
c     xjac(i,j) = sum_I {ngshp(I,k)*X(l,I)}
c
      do k=1,mcrd
      do l=1,mcrd
         xjac(k,l)=0.d0
         do i=1,nnode
            xjac(k,l)=xjac(k,l)+xngshpquad(i,k)*coords(l,i)
         end do
      end do
      end do
      call kinverse(xjac,xjacinv,detjac,mcrd)
c
c     Form the Jacobian Matrix [JI] -> xjimtx
c
      call kinitia(xjimtx,(nsize+1)*(nsize+1))
c      do 10 k=1,ndf-1 !Magnetic case
      do 10 k=1,mcrd
      do 10 i=1,mcrd
      do 10 j=1,mcrd
         xjimtx(mcrd*(k-1)+i,mcrd*(k-1)+j)=xjacinv(i,j)
 10   continue
c
c     For the electric case the last 3x3 entries
c     of [JI] must be modified according to
c     E1=-dphi/dX1 E2=-dphi/dX2 E3=-dphi/dX3
c
      xjimtx(5,5)=1.0d0        ! Axisymmetric u/R
      xjimtx(6,6)=-xjacinv(1,1)
      xjimtx(6,7)=-xjacinv(1,2)
      xjimtx(7,6)=-xjacinv(2,1)
      xjimtx(7,7)=-xjacinv(2,2)
c
c     Pressure dof
c
      xjimtx(nsize+1,nsize+1)=1.0d0
c
c     Form the NG matrix [NG] -> ngmtxc
c
c
c
c      
      call kinitia(xngmtx,(nsize+1)*ndofel)
	  do 11 k=1,ndf-1
	  do 11 i=1,mcrd
	  do 11 j=1,nnode
	    xngmtx(mcrd*(k-1)+i,ndf*(j-1)+k)=xngshpquad(j,i)
 11   continue
      do j=1,nnode
         xngmtx(5,ndf*(j-1)+1)=xnshp(j)/R
      end do
      do j=1,nnode
         xngmtx(6,ndf*(j-1)+3)=xngshpquad(j,1)
         xngmtx(7,ndf*(j-1)+3)=xngshpquad(j,2)
      end do
      do i=1,ndpress
         xngmtx(nsize+1,ndofel-ndpress+i)=xnshppress(i)
      end do
c
c
c     Compute [G]=[JI] x [NG] -> gmtx(i,j)=xjimtx(i,k)*xngmtx(k,j)
c
      call kinitia(gmtx,(nsize+1)*ndofel)
c
      call kmult(xjimtx,xngmtx,gmtx,nsize+1,nsize+1,ndofel)
      return
      end

c
c***********************************************************************
c***********************************************************************
c
c     Compute deformation gradient at the nodes given the displacement
c     and gradient matrix G

      subroutine kfbgradE(F,E,uvec,gmtx,ndf,mcrd,ndofel
     +           ,nsize)

      implicit double precision(a-h,o-z)

      dimension  uvec(ndofel),gmtx(nsize+1,ndofel),
     + gradu(nsize+1),F(mcrd+1,mcrd+1),E(mcrd+1)

      call kinitia(gradu,nsize+1)
      call kmult(gmtx,uvec,gradu,nsize+1,ndofel,1)
c      write(*,*) gradu
c
c     Following lines should be modified for 3d or generalized pl-e
c     or axisymmetry
c
      call kinitia(F,(mcrd+1)*(mcrd+1))
      call kinitia(E,mcrd+1)
      F(1,1)=1.d0+gradu(1)
      F(1,2)=gradu(2)
      F(2,1)=gradu(3)
      F(2,2)=1.d0+gradu(4)
      F(3,3)=1.0+gradu(5)
      E(1)=gradu(6)
      E(2)=gradu(7)
c
      return
      end
c
c***********************************************************************
c
c     Compute deformation gradient at the nodes given the displacement
c     and gradient matrix G

      subroutine kfbgradEpressNI(F,E,press,uvec,gmtx,ndf,mcrd,ndofel
     +           ,nsize,kinc,kstep,props,nprops)

      implicit double precision(a-h,o-z)

      dimension  uvec(ndofel),gmtx(nsize+1,ndofel),
     + gradupress(nsize+1),F(mcrd+1,mcrd+1),E(mcrd+1)
      dimension props(nprops)

      call kinitia(gradupress,nsize+1)
      call kmult(gmtx,uvec,gradupress,nsize+1,ndofel,1)
c      write(*,*) gradupress
c
c     Following lines should be modified for 3d or generalized pl-e
c     or axisymmetry
c
      call kinitia(F,(mcrd+1)*(mcrd+1))
      call kinitia(E,mcrd+1)
      F(1,1)=1.d0+gradupress(1)
      F(1,2)=gradupress(2)
      F(2,1)=gradupress(3)
      F(2,2)=1.d0+gradupress(4)
      F(3,3)=1.0+gradupress(5)
      E(1)=gradupress(6)
      E(2)=gradupress(7)
      gam_st = props(13)
      pgam = -2.0d0*gam_st/props(12)
      ! pgam = 0.0d0


      if((kinc.lt.1).and.(kstep.eq.1)) then
            press=-props(1) -props(11) + pgam
      else 
            press=gradupress(nsize+1)
      end if
c
      return
      end
c
c***********************************************************************
c

      subroutine kgausspoints(wgp,xigp,ngp,mcrd)

      implicit double precision (a-h,o-z)
      integer ngp,mcrd
      dimension wgp(ngp),xigp(mcrd,ngp)

      if(ngp.eq.1) then
         wgp(1)=0.5d0
         xigp(1,1)=0.3333333333333333d0
         xigp(2,1)=0.3333333333333333d0
      else if(ngp.eq.3) then
         wgp(1)=0.166666666666667d0
         wgp(2)=0.166666666666667d0
         wgp(3)=0.166666666666667d0
         xigp(1,1)=0.1666666666666667d0
         xigp(2,1)=0.1666666666666667d0
         xigp(1,2)=0.6666666666666667d0
         xigp(2,2)=0.6666666666666667d0
         xigp(1,3)=0.1666666666666667d0
         xigp(1,3)=0.1666666666666667d0
         xigp(2,3)=0.6666666666666667d0
      else if(ngp.eq.4) then
	     wgp(1)=0.281250000000000
         wgp(2)=0.260416666666667
         wgp(3)=0.260416666666667
		 wgp(4)=0.260416666666667
		 
         xigp(1,1)=0.333333333333333d0
         xigp(2,1)=0.333333333333333d0
         xigp(1,2)=0.2d0
         xigp(2,2)=0.2d0
         xigp(1,3)=0.6d0
         xigp(2,3)=0.2d0
		 xigp(1,4)=0.2d0
		 xigp(2,4)=0.6d0
      else if(ngp.eq.6) then
         wgp(1)=0.111690794839006
		 wgp(2)=0.111690794839006
		 wgp(3)=0.111690794839006
		 wgp(4)=0.054975871827661
		 wgp(5)=0.054975871827661
		 wgp(6)=0.054975871827661
c
         xigp(1,1)=0.108103018168070d0
		 xigp(2,1)=0.445948490915965d0
		 xigp(1,2)=0.445948490915965d0
		 xigp(2,2)=0.108103018168070d0
c
         xigp(1,3)=0.445948490915965d0
		 xigp(2,3)=0.445948490915965d0
		 xigp(1,4)=0.816847572980459d0
		 xigp(2,4)=0.091576213509771d0
		 xigp(1,5)=0.091576213509771d0
		 xigp(2,5)=0.816847572980459d0
		 xigp(1,6)=0.091576213509771d0
		 xigp(2,6)=0.091576213509771d0
	  else
	     write(*,*)'Unknown number of Gauss points'
		 write(7,*)'Unknown number of Gauss points'
      end if

      return
      end

c
c***********************************************************************
c
      subroutine kidtens(A,mcrd)
c
      implicit double precision (a-h, o-z)
      dimension A(mcrd,mcrd)
c
      call kinitia(A,mcrd*mcrd)
      do i=1,mcrd
        A(i,i)=1.d0
      enddo
c
      return
      end
c
c***********************************************************************
c
      subroutine kdet(deta,a,mcrd)

      implicit double precision(a-h,o-z)
c
      dimension a(mcrd,mcrd)

      if (mcrd.eq.1) then
         deta=a(1,1)
         return
      endif

      if (mcrd.eq.2) then
         deta=a(1,1)*a(2,2)-a(1,2)*a(2,1)
         return
      endif

      if (mcrd.eq.3) then
         deta = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2))-
     +       a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1))+
     +       a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
         return
      end if

      return
      end
c
c***********************************************************************
c
      subroutine kinverse(a,ainv,deta,mcrd)

      implicit double precision(a-h,o-z)
c
      dimension a(mcrd,mcrd),ainv(mcrd,mcrd)


      call kinitia(ainv,mcrd*mcrd)
      call kdet(deta,a,mcrd)
      if (dabs(deta).lt.1.d-30) then
         write(*,*) 'Singular Det found - ** Program STOP **'
         stop
      endif

      if (mcrd.eq.1) then
         ainv(1,1)=1.d0/a(1,1)
         return
      endif

      if (mcrd.eq.2) then
         ainv(1,1)=a(2,2)/deta
         ainv(1,2)=-a(1,2)/deta
         ainv(2,1)=-a(2,1)/deta
         ainv(2,2)=a(1,1)/deta
         return
      endif

      if (mcrd.eq.3) then
         ainv(1,1) = (-a(2,3)*a(3,2) + a(2,2)*a(3,3))/deta
         ainv(1,2) = ( a(1,3)*a(3,2) - a(1,2)*a(3,3))/deta
         ainv(1,3) = (-a(1,3)*a(2,2) + a(1,2)*a(2,3))/deta
         ainv(2,1) = ( a(2,3)*a(3,1) - a(2,1)*a(3,3))/deta
         ainv(2,2) = (-a(1,3)*a(3,1) + a(1,1)*a(3,3))/deta
         ainv(2,3) = ( a(1,3)*a(2,1) - a(1,1)*a(2,3))/deta
         ainv(3,1) = (-a(2,2)*a(3,1) + a(2,1)*a(3,2))/deta
         ainv(3,2) = ( a(1,2)*a(3,1) - a(1,1)*a(3,2))/deta
         ainv(3,3) = (-a(1,2)*a(2,1) + a(1,1)*a(2,2))/deta
         return
      endif

      return
      end
c
c***********************************************************************
c
      subroutine kinitia(a,n)
c
      implicit double precision(a-h,o-z)
c
      dimension a(n)
c
      do i=1,n
        a(i)=0.d0
      end do
c
      return
      end
c
c
C***********************************************************************
c
      subroutine kmult(a,b,c,l,m,n)
c
      implicit double precision(a-h,o-z)

      dimension a(l,m),b(m,n),c(l,n)
c
      do 10 i=1,l
      DO 10 j=1,n
      aux=0.d0
      do 20 k=1,m
 20   aux=aux+a(i,k)*b(k,j)
 10   c(i,j)=aux
c
      return
      end
c
c***********************************************************************
c
      subroutine ktranspose(a,at,m,n)
c
      implicit double precision (a-h,o-z)
c
      dimension a(m,n),at(n,m)
c
      do 10 j=1,n
      do 10 i=1,m
 10   at(j,i)=a(i,j)
c
      return
      end
c
c**********************************************************************
C***********************************************************************
c
      subroutine effprops(props,nprops,nup,propsnew1)
c
      implicit double precision (a-h,o-z)
      integer nprops,medium
      dimension props(nprops),propsnew1(nprops)
      real*8 gm,mu0,mup,c,ms,denom,nup

      call kinitia(propsnew1,nprops)

      gm=props(1)     !shear modulus of matrix
      mu0=props(3)    ! permeability of air
      mup=props(4)    ! permeability of particle
      c=props(5)      ! concentration of particle
      ms=props(6)     ! magnetic saturation
      medium=props(7) ! medium=0 air medium=1 homogenized medium

      if (medium .eq. 0) then
        propsnew1(1)=gm
        propsnew1(2)=propsnew1(1)*10.0d0     ! change to 10000.0d0 for incompressible air
        propsnew1(3)=mu0
        propsnew1(4)=mu0
        propsnew1(8)=3.0
        propsnew1(9)=1.00
        propsnew1(10)=1.00
        propsnew1(11)=0.00
c		write(7,*),'THIS IS AIR'

      else 
c	    write(7,*),'THIS IS MAT'
        denom=5.0d0*((2.0d0 + c)*mu0+(1.0d0-c)*nup)**2.0d0


        propsnew1(1)=gm                                    ! effective shear modulus
        propsnew1(3)=mu0+(3.0d0*c*mu0*(nup-mu0))/((2.0d0+c)*mu0+(1.0d0-c)*nup)        !AEPS
        propsnew1(4)=mu0+(3.0d0*c*(10.0d0+2.0d0*c
     +               +3.0d0*(c**2.0d0))*(nup-mu0)*(mu0**2.0d0))/denom
     +               +(3.0d0*c*(1.0d0-c)*(5.0d0+3.0d0*c)*(nup-mu0)*mu0*nup)/denom     !AMKE
        propsnew1(5)=c                                                                !c
        propsnew1(6)=props(6)                                                         ! ms
        propsnew1(7)=props(7)                                                         ! medium
        propsnew1(8)=props(8)                                                              ! N
        propsnew1(9)=props(9)                                                  ! alpha1
        propsnew1(10)=props(10)                                                     ! alpha2
        propsnew1(11)=props(11)                                        ! AG2
        propsnew1(2)=(gm+propsnew1(11))*10000.0d0                                           ! lambda =10^4*effective shear modulus

      end if

c      write(7,*),'denom=',denom,'medium',medium

      return
      end
c
c***********************************************************************
C***********************************************************************
c
      subroutine NRnupsolve(Ifour,Ifive,props,nprops,nup1)
c
      implicit double precision (a-h,o-z)
      real*8 Ifour,Ifive,diff,nup0,mu0,nuptemp,EQNprime1,EQN1,nup1
      dimension props(nprops)
	  integer nprops,medium
	  medium=props(7) ! medium=0 air medium=1 homogenized medium

      nup0=props(4)
      mu0=props(3)

      nup1=nup0;
      diff=1.0d0
	  
	  if (medium .eq. 0) then
	      nup1=nup0
	  else
          if (Ifive .lt. 10.0d0**3.0d0) then
             nup1=nup0
          else 

             do while(ABS(diff)/nup1 .gt. 0.0001)
                call Equation(Ifour,Ifive,nup1,props,nprops,EQN1,EQNprime1)
                nuptemp=nup1-EQN1/EQNprime1;
                diff=nup1-nuptemp
                nup1=nuptemp
             end do
          end if
      end if

      return
      end
c
c***********************************************************************
C***********************************************************************
c
      subroutine Equation(I4,I5,nup,props,nprops,EQN,EQNprime)
c
      implicit double precision (a-h,o-z)
      real*8 I4,I5,nup,EQN,EQNprime,k1
      real*8 coth,cosech,arg,denom,k2
      real*8 ms,mu0,c,mup,denom2
      dimension props(nprops)

      c=props(5)
      mu0=props(3)
      mup=props(4)
      ms=props(6)

      denom=5.0*((2.0 + c)*mu0+(1.0-c)*nup)**3.0d0
      denom2=5.0*((2.0 + c)*mu0+(1.0-c)*nup)**4.0d0
      
      k1=-54.0*(1.0-c)*c*I4*(mu0**2.0d0)*(-mu0+nup)/denom
     +     +9*I5*(mu0**2.0d0)*((10.0-c+6*c**2.0)*mu0
     +     +(5.0+c-6.0*c**2.0)*nup)/denom

      k2=162.0*((1-c)**2.0d0)*c*I4*(mu0**2.0d0)*(-mu0+nup)/denom2
     +   -54.0*(1.0-c)*c*I4*(mu0**2.0d0)/denom
     +   +9*(5.0+c-6*(c**2.0d0))*I5*(mu0**2.0d0)/denom
     +   -(27.0*(1.0-c)*I5*(mu0**2.0d0)*((10.0-c+6*(c**2.0d0))*mu0
     +   +(5.0+c-6*(c**2.0d0))*nup))/denom2

      arg=3.0*SQRT(k1)*(-mu0+mup)/(ms*mu0)
      call COTHYP(arg,coth)
      call CSCHYP(arg,cosech)
c      write(7,*),'arg=',arg,'coth=',coth,'cosech=',cosech

      EQN=-nup+(1.0/3.0)*mu0*(3.0
     +    +(ms**2.0d0)*mu0/(k1*mu0-k1*mup)
     +    +3.0*ms*coth/SQRT(k1))


      
      EQNprime=-1+(1.0/3.0)*mu0*(-3.0*k2*coth*ms/(2.0d0*(k1**1.50))
     +         -9.0*k2*(cosech**2.0d0)*(-mu0-mup)/(2.0*k1*mu0)
     +         -(ms**2.0d0)*mu0*(k2*mu0-k2*mup)/((k1*mu0
     +         -k1*mup)**2.0d0))

c
      return
      end
c
c***********************************************************************
C***********************************************************************
c
      subroutine COTHYP(ARG,cothans)
c
      implicit double precision (a-h,o-z)
      real*8 ARG,cothans
      
      cothans=1/TANH(ARG)
      
c
      return
      end
c
c***********************************************************************
C***********************************************************************
c
      subroutine CSCHYP(ARG,cosechans)
c
      implicit double precision (a-h,o-z)
      real*8 ARG,cosechans
      
      cosechans=1/SINH(ARG)
      
c
      return
      end
c
c***********************************************************************
C***********************************************************************
c
      subroutine COSHYP(ARG,cosh)
c
      implicit double precision (a-h,o-z)
      real*8 ARG,cosh
      
      cosh=SINH(ARG)/TANH(ARG)
      
c
      return
      end
c
c***********************************************************************
      subroutine cauchystr(Ft,xfe1,xfe2,
     +           xfe3,xfe4,xfe5,J,T)
c
      implicit double precision (a-h,o-z)
      real*8 xfe1,xfe2,xfe3,xfe4,xfe5,J
      dimension Ft(3,3),T(3,3)
      dimension S(3,3)
      
      call kinitia(S,3*3)
      S(1,1)=xfe1/J
      S(1,2)=xfe2/J
      S(2,1)=xfe3/J
      S(2,2)=xfe4/J
      S(3,3)=xfe5/J
      call kmult(S,Ft,T,3,3,3)
      
c
      return
      end
c**********************************************************************
c**********************************end*********************************
c**********************************************************************
      subroutine ksurgausspoints(wgp,xigp)

      implicit double precision (a-h,o-z)
      dimension wgp(5),xigp(5)
c      dimension wgpn(5),xigpn(5)

c      wgp(1) = 0.5688888888888890d0
c      wgp(2) = 0.4786286704993660d0
c	wgp(3) = 0.4786286704993660d0
c	wgp(4) = 0.2369268850561890d0
c	wgp(5) = 0.2369268850561890d0

      wgp(1) = 0.2844444444444440d0
      wgp(2) = 0.2393143352496830d0
	wgp(3) = 0.2393143352496830d0
	wgp(4) = 0.1184634425280950d0
	wgp(5) = 0.1184634425280950d0
c
c      xigp(1) = 0.0d0
c	xigp(2) = 0.5384693101056830d0
c	xigp(3) = -0.5384693101056830d0
c	xigp(4) = 0.9061798459386640d0
c      xigp(5) = -0.9061798459386640d0

      xigp(1) = 0.500000000000000d0 
	xigp(2) = 0.7692346550528410d0
	xigp(3) = 0.2307653449471590d0  
	xigp(4) = 0.9530899229693320d0
      xigp(5) = 0.0469100770306680d0
	

      return
      end
c**********************************************************************
c**********************************************************************
c     Compute the surface deformation gradient evaluated at surface gauss poins 

      subroutine ksurfF(xi,Fsurf,refN,Po,gmtxs,F,gmtx,uvec,ndf,mcrd,
     +                coords,ndofel,nsize,props,nprops,nnode)

      implicit double precision(a-h,o-z)

      dimension  uvec(ndofel),gmtx(nsize+1,ndofel),
     + gradupress(nsize+1),F(mcrd+1,mcrd+1),Fsurf(mcrd+1,mcrd+1),
     + props(nprops),refN(mcrd+1),Po(mcrd+1,mcrd+1),
     + gmtxs(nsize+1,ndofel),PP(nsize+1,nsize+1),
     + coords(mcrd,nnode),xnshp(nnode),xi(mcrd),
     + xip(mcrd)

      call kinitia(gradupress,nsize+1)
      call kinitia(refN,mcrd+1)
      call kinitia(Po,mcrd+1,mcrd+1)
      call kinitia(xip,mcrd)
      call kmult(gmtx,uvec,gradupress,nsize+1,ndofel,1)
c
c   F is the bulk deformation gradient evaluated at the current surface gauss point
c**********************************************************************
c**********************************************************************
      call kinitia(F,(mcrd+1)*(mcrd+1))
      F(1,1)=1.d0+gradupress(1)
      F(1,2)=gradupress(2)
      F(2,1)=gradupress(3)
      F(2,2)=1.d0+gradupress(4)
      F(3,3)=1.0+gradupress(5)

c     CR shape functions
    
      xnshp(1)=xi(1)*(-1.0+2.0*xi(1)-
     + 3.0*xi(2)*(-1.0+xi(2)+xi(1)))
      xnshp(2)=xi(2)*(-1.0+2.0*xi(2)-
     + 3.0*xi(1)*(-1.0+xi(2)+xi(1)))
      xnshp(3)=-((-1.0+xi(2)+xi(1))*(1.0-
     + 2.0*xi(1)+xi(2)*(-2.0+3.0*xi(1))))
      xnshp(4)=4.0*xi(2)*xi(1)*(-2.0+3.0*xi(2)+3.0*xi(1))
      xnshp(5)=4.0*xi(2)*(-1.0+xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xnshp(6)=4.0*(-1.0+3.0*xi(2))*xi(1)*(-1.0+xi(2)+xi(1))
      xnshp(7)=-27.0*xi(2)*xi(1)*(-1.0+xi(2)+xi(1))


      do i=1,mcrd
        do j=1,nnode
            xip(i)=xip(i) + xnshp(j)*coords(i,j)
        enddo
      enddo

c   Normal to the undeformed surface N=refN

      do icount=1,mcrd
        refN(icount) = xip(icount)/props(12)
      end do 

c   Compute the projection to the undeformed surface Po=I-N \otimes N
      Po(1,1) = 1.0d0 - refN(1)*refN(1)
      Po(1,2) = -refN(1)*refN(2)

      Po(2,1) = -refN(2)*refN(1)
      Po(2,2) = 1.0d0 - refN(2)*refN(2)

      Po(3,3) = 1.0d0

      ! write(7,*)'N1=',xnshp(1),'N2=',xnshp(2),'N3=',xnshp(3)
      ! write(7,*)'N4=',xnshp(4),'N5=',xnshp(5),'N6=',xnshp(6)
      ! write(7,*)'N7=',xnshp(7)
      ! write(7,*)'Nr=',refN(1),'Nz=',refN(2),'Nt=',refN(3)
      ! write(7,*)'|N|=',sqrt(refN(1)**2.0d0 + refN(2)**2.0d0)
      ! write(*,*)'nrpops=',nprops,'nsize=',nsize

c   Compute the surface deformation gradient at the current surface Gauss point
c   Fsurf=FPo
      call kmult(F,Po,Fsurf,mcrd+1,mcrd+1,mcrd+1)

c   Find gmtxs=PP*gmtx
      call kinitia(PP,(nsize+1)*(nsize+1))

      do 10 k=1,2
      do 10 i=1,mcrd
      do 10 j=1,mcrd
         PP(mcrd*(k-1)+i,mcrd*(k-1)+j)=Po(i,j)
 10   continue

      PP(5,5)=1.0d0
      PP(6,6)=1.0d0
      PP(7,7)=1.0d0
      PP(8,8)=1.0d0

      call kmult(PP,gmtx,gmtxs,nsize+1,nsize+1,ndofel)

c      write(7,*)'********PPPPPPPP********'
c      do i=1,nsize
c            write(7,*) 'PP(',i,':)',(PP(i,j),j=1,nsize+1)
c      end do

      return
      end
c**********************************************************************
c**********************************************************************
c     Compute Deformed Surface Normal "n" and the projection tensor
c     P=I-n\otimes n
      subroutine ksurnormal(currN,P,refN,F,mcrd)
c
      implicit double precision(a-h,o-z)
      
      dimension currN(mcrd+1),refN(mcrd+1)
      dimension F(mcrd+1,mcrd+1)
      dimension Finv(mcrd+1,mcrd+1),Finvt(mcrd+1,mcrd+1)
      dimension FinvtN(mcrd+1),P(mcrd+1,mcrd+1)
      real*8    detF,detFinvtN
      
      call kinverse(F,Finv,detF,mcrd+1)
      call ktranspose(Finv,Finvt,mcrd+1,mcrd+1)
      call kmult(Finvt,refN,FinvtN,mcrd+1,mcrd+1,1)
      
      detFinvtN=SQRT(FinvtN(1)**2.0d0+
     +     FinvtN(2)**2.0d0+FinvtN(3)**2.0d0)
      
c   n=F^(-T)N/|F^(-T)N|
      
      currN(1)=FinvtN(1)/detFinvtN
      currN(2)=FinvtN(2)/detFinvtN
      currN(3)=FinvtN(3)/detFinvtN
      
c   Compute the projection to the deformed surface P=I-n \otimes n
      P(1,1)=1.0d0-currN(1)*currN(1)
      P(1,2)=-currN(1)*currN(2)
      P(1,3)=-currN(1)*currN(3)
      
      P(2,1)=-currN(2)*currN(1)
      P(2,2)=1.0d0-currN(2)*currN(2)
      P(2,3)=-currN(2)*currN(3)
      
      P(3,1)=-currN(3)*currN(1)
      P(3,2)=-currN(3)*currN(2)
      P(3,3)=1.0d0-currN(3)*currN(3)

      ! write(7,*)'Nr=',currN(1),'Nz=',currN(2),'Nt=',currN(3)
      ! write(7,*)'|N|=',sqrt(currN(1)**2.0d0 + currN(2)**2.0d0 + currN(3)**2.0d0)
      
      return
      end
c
c***********************************************************************
      subroutine kumatsurften(xke,xfe,F,Fsurf,refn,Po,
     +           currn,P,gamma,nsize,mcrd,Jhat)
 
      implicit double precision(a-h,o-z)
 
      dimension F(mcrd+1,mcrd+1),Finv(mcrd+1,mcrd+1)
      dimension Finvt(mcrd+1,mcrd+1),Fsurf(mcrd+1,mcrd+1)
      dimension Po(mcrd+1,mcrd+1),P(mcrd+1,mcrd+1)
      dimension fhat(mcrd+1,mcrd+1),fhatinv(mcrd+1,mcrd+1)
      dimension fhatt(mcrd+1,mcrd+1),fhatfhatt(mcrd+1,mcrd+1)
      dimension perp(mcrd+1,mcrd+1)
      dimension currn(mcrd+1),refn(mcrd+1)
      dimension auxN(mcrd+1)
      dimension xke(nsize+1,nsize+1),xkeA(nsize+1,nsize+1)
      dimension xkeD1(nsize+1,nsize+1),xkeD2(nsize+1,nsize+1)
      dimension xfe(nsize+1)
      real*8 Jhat,mag,detF,gamma
 
      call kinitia(Finv,(mcrd+1)*(mcrd+1))
      call kinitia(Finvt,(mcrd+1)*(mcrd+1))
      call kinitia(fhat,(mcrd+1)*(mcrd+1))
      call kinitia(fhatinv,(mcrd+1)*(mcrd+1))
      call kinitia(fhatt,(mcrd+1)*(mcrd+1))
      call kinitia(fhatfhatt,(mcrd+1)*(mcrd+1))
      call kinitia(iperp,(mcrd+1)*(mcrd+1))
     
 
      call kinitia(xke,(nsize+1)*(nsize+1))
      call kinitia(xkeA,(nsize+1)*(nsize+1))
      call kinitia(xkeD1,(nsize+1)*(nsize+1))
      call kinitia(xkeD2,(nsize+1)*(nsize+1))
      call kinitia(xfe,nsize+1)
 
c     Surface Det[Fhat]=Jhat=J|F^(-T)N|  
 
      call kinverse(F,Finv,detF,mcrd+1)
      call ktranspose(Finv,Finvt,mcrd+1,mcrd+1)
c     auxN=F^(-T)N
      call kmult(Finvt,refn,auxN,mcrd+1,mcrd+1,1)
 
      mag=SQRT(auxN(1)**2.0d0+auxN(2)**2.0d0+auxN(3)**2.0d0)
      Jhat=detF*mag
 
c     Find fhat (inverse of the surface deformation gradient) fhat=Fhat^(-1)P 
      call kmult(Finv,P,fhat,mcrd+1,mcrd+1,mcrd+1)
      call ktranspose(fhat,fhatt,mcrd+1,mcrd+1)
 
c     perp=I-P=n \otimes n
      perp(1,1) = currn(1)*currn(1)     
      perp(1,2) = currn(1)*currn(2)
      perp(1,3) = currn(1)*currn(3)
 
      perp(2,1) = currn(2)*currn(1)
      perp(2,2) = currn(2)*currn(2)
      perp(2,3) = currn(2)*currn(3)
 
      perp(3,1) = currn(3)*currn(1)
      perp(3,2) = currn(3)*currn(2)
      perp(3,3) = currn(3)*currn(3)
 
c     Fill xkeA= f^T \otimes f^T    
 
      xkeA(1,1) = fhatt(1,1)*fhatt(1,1)
      xkeA(1,2) = fhatt(1,1)*fhatt(1,2)
      xkeA(1,3) = fhatt(1,1)*fhatt(2,1)
      xkeA(1,4) = fhatt(1,1)*fhatt(2,2)
      xkeA(1,5) = fhatt(1,1)*fhatt(3,3)
 
      xkeA(2,1) = xkeA(1,2)
      xkeA(2,2) = fhatt(1,2)*fhatt(1,2)
      xkeA(2,3) = fhatt(1,2)*fhatt(2,1)
      xkeA(2,4) = fhatt(1,2)*fhatt(2,2)
      xkeA(2,5) = fhatt(1,2)*fhatt(3,3)
 
      xkeA(3,1)=xkeA(1,3)
      xkeA(3,2)=xkeA(2,3)
      xkeA(3,3)=fhatt(2,1)*fhatt(2,1)
      xkeA(3,4)=fhatt(2,1)*fhatt(2,2)
      xkeA(3,5)=fhatt(2,1)*fhatt(3,3)
 
      xkeA(4,1)=xkeA(1,4)
      xkeA(4,2)=xkeA(2,4)
      xkeA(4,3)=xkeA(3,4)
      xkeA(4,4)=fhatt(2,2)*fhatt(2,2)
      xkeA(4,5)=fhatt(2,2)*fhatt(3,3)

 
      xkeA(5,1)=xkeA(1,5)
      xkeA(5,2)=xkeA(2,5)
      xkeA(5,3)=xkeA(3,5)
      xkeA(5,4)=xkeA(4,5)
      xkeA(5,5)=fhatt(3,3)*fhatt(3,3)
      
 
 
c     xkeD1 (-fhatt(il)Fhat(jk))
 
      xkeD1(1,1)=-fhatt(1,1)*fhat(1,1)
      xkeD1(1,2)=-fhatt(1,2)*fhat(1,1)
      xkeD1(1,3)=-fhatt(1,1)*fhat(1,2)
      xkeD1(1,4)=-fhatt(1,2)*fhat(1,2)
      xkeD1(1,5)=-fhatt(1,3)*fhat(1,3)

 
 
      xkeD1(2,1)=-fhatt(1,1)*fhat(2,1)
      xkeD1(2,2)=-fhatt(1,2)*fhat(2,1)
      xkeD1(2,3)=-fhatt(1,1)*fhat(2,2)
      xkeD1(2,4)=-fhatt(1,2)*fhat(2,2)
      xkeD1(2,5)=-fhatt(1,3)*fhat(2,3)
 
 
      xkeD1(3,1)=-fhatt(2,1)*fhat(1,1)
      xkeD1(3,2)=-fhatt(2,2)*fhat(1,1)
      xkeD1(3,3)=-fhatt(2,1)*fhat(1,2)
      xkeD1(3,4)=-fhatt(2,2)*fhat(1,2)
      xkeD1(3,5)=-fhatt(2,3)*fhat(1,3)
 
 
      xkeD1(4,1)=-fhatt(2,1)*fhat(2,1)
      xkeD1(4,2)=-fhatt(2,2)*fhat(2,1)
      xkeD1(4,3)=-fhatt(2,1)*fhat(2,2)
      xkeD1(4,4)=-fhatt(2,2)*fhat(2,2)
      xkeD1(4,5)=-fhatt(2,3)*fhat(2,3)
 
 
      xkeD1(5,1)=-fhatt(3,1)*fhat(3,1)
      xkeD1(5,2)=-fhatt(3,2)*fhat(3,1)
      xkeD1(5,3)=-fhatt(3,1)*fhat(3,2)
      xkeD1(5,4)=-fhatt(3,2)*fhat(3,2)
      xkeD1(5,5)=-fhatt(3,3)*fhat(3,3)
  
 
c     xkeD2= (n \otimes n)(ik)fhatfhatt(jl)
 
      call kmult(fhat,fhatt,fhatfhatt,mcrd+1,mcrd+1,mcrd+1)
 
      xkeD2(1,1)=perp(1,1)*fhatfhatt(1,1)
      xkeD2(1,2)=perp(1,1)*fhatfhatt(1,2)
      xkeD2(1,3)=perp(1,2)*fhatfhatt(1,1)
      xkeD2(1,4)=perp(1,2)*fhatfhatt(1,2)
      xkeD2(1,5)=perp(1,3)*fhatfhatt(1,3)
 
 
      xkeD2(2,1)=perp(1,1)*fhatfhatt(2,1)
      xkeD2(2,2)=perp(1,1)*fhatfhatt(2,2)
      xkeD2(2,3)=perp(1,2)*fhatfhatt(2,1)
      xkeD2(2,4)=perp(1,2)*fhatfhatt(2,2)
      xkeD2(2,5)=perp(1,3)*fhatfhatt(2,3)
 
 
      xkeD2(3,1)=perp(2,1)*fhatfhatt(1,1)
      xkeD2(3,2)=perp(2,1)*fhatfhatt(1,2)
      xkeD2(3,3)=perp(2,2)*fhatfhatt(1,1)
      xkeD2(3,4)=perp(2,2)*fhatfhatt(1,2)
      xkeD2(3,5)=perp(2,3)*fhatfhatt(1,3)

 
      xkeD2(4,1)=perp(2,1)*fhatfhatt(2,1)
      xkeD2(4,2)=perp(2,1)*fhatfhatt(2,2)
      xkeD2(4,3)=perp(2,2)*fhatfhatt(2,1)
      xkeD2(4,4)=perp(2,2)*fhatfhatt(2,2)
      xkeD2(4,5)=perp(2,3)*fhatfhatt(2,3)
 
 
      xkeD2(5,1)=perp(3,1)*fhatfhatt(3,1)
      xkeD2(5,2)=perp(3,1)*fhatfhatt(3,2)
      xkeD2(5,3)=perp(3,2)*fhatfhatt(3,1)
      xkeD2(5,4)=perp(3,2)*fhatfhatt(3,2)
      xkeD2(5,5)=perp(3,3)*fhatfhatt(3,3)

 
c     xke=gamma*Jhat*(xkeA+xkeD1+xkeD2)
 
      do i=1,nsize+1
          do j=1,nsize+1
            xke(i,j)=Jhat*gamma*(xkeA(i,j)+xkeD1(i,j)+
     +           xkeD2(i,j))
          end do
      end do
 
c     xfe=gamma*Jhat*fhatt
 
      xfe(1)=gamma*Jhat*fhatt(1,1)
      xfe(2)=gamma*Jhat*fhatt(1,2)
      xfe(3)=gamma*Jhat*fhatt(2,1)
      xfe(4)=gamma*Jhat*fhatt(2,2)
      xfe(5)=gamma*Jhat*fhatt(3,3)

 
      return
      end
c***********************************************************************
c***********************************************************************
c
      subroutine ksurfjac(detjacsurf,coordsabq,
     +     xi,ndf,mcrd,nnode,ndofel,nsize,iwr,mtype)

      implicit double precision(a-h,o-z)

      dimension coordsabq(mcrd,nnode)
      dimension coords(mcrd,nnode)
      dimension xi(mcrd)
      dimension xngshpquad(nnode,mcrd)
      dimension xjac(mcrd,mcrd),xjact(mcrd,mcrd)
      dimension xdX(mcrd),xjactn(mcrd)
      real*8 detjacsurf

      call kinitia(xjact,mcrd)
      call kinitia(xdX,mcrd)
      call kinitia(xjactn,mcrd)
      if (mtype .eq. 1) then
c      write(7,*)'this is surface=',mtype
c     dX Face 1 is (-1,1)
            xdX(1) = -1.0d0
            xdX(2) = 1.0d0
      else if (mtype .eq. 2) then
c     Normal of Face 2 is (1, 0) 
            xdX(1) = 1.0d0
            xdX(2) = 0.0d0
      else if (mtype .eq. 3) then
c     Normal of Face 3 is (-1,0,0)       
            xdX(1) = 0.0d0
            xdX(2) = -1.0d0

      end if
      
c
c     Extract coordinates of mesh nodes only
c
      do i=1,mcrd
         do j=1,nnode
            coords(i,j)=coordsabq(i,j)
         end do
      end do

c
c     Quadratic Shape functions derivatives
c
c     ngshp(j,i)=dN(j)/dX(i)


      xngshpquad(1,1)=-1.0-3.0*xi(2)**2.0+xi(2)*(3.0-6.0*xi(1))+4.0*xi(1)
      xngshpquad(1,2)=-3.0*xi(1)*(-1+2.0*xi(2)+xi(1))
      xngshpquad(2,1)=-3.0*xi(2)*(-1.0+xi(2)+2.0*xi(1))
      xngshpquad(2,2)=-1.0 + xi(2)*(4.0-6.0*xi(1))+3.0*xi(1)-3.0*xi(1)**2.0
      xngshpquad(3,1)=-3.0-3.0*xi(2)**2.0+xi(2)*(7.0-6.0*xi(1))+4.0*xi(1)
      xngshpquad(3,2)=-3.0+xi(2)*(4.0-6.0*xi(1))+7.0*xi(1)-3.0*xi(1)**2.0
      xngshpquad(4,1)=4.0*xi(2)*(-2.0+3.0*xi(2)+6.0*xi(1))
      xngshpquad(4,2)=4.0*xi(1)*(-2.0+6.0*xi(2)+3.0*xi(1))
      xngshpquad(5,1)=4.0*xi(2)*(-4.0+3.0*xi(2)+6.0*xi(1))
      xngshpquad(5,2)=4.0*(-1.0+2.0*xi(2)+xi(1))*(-1.0+3.0*xi(1))
      xngshpquad(6,1)=4.0*(-1.0+3.0*xi(2))*(-1.0+xi(2)+2.0*xi(1))
      xngshpquad(6,2)=4.0*xi(1)*(-4.0+6.0*xi(2)+3.0*xi(1))
      xngshpquad(7,1)=-27.0*xi(2)*(-1+xi(2)+2.0*xi(1))
      xngshpquad(7,2)=-27.0*xi(1)*(-1+2.0*xi(2)+xi(1))

c
c     xjac=F^(T), xjact=F^T
c
      do k=1,mcrd
      do l=1,mcrd
         xjac(k,l)=0.d0
         do i=1,nnode
            xjac(k,l)=xjac(k,l)+xngshpquad(i,k)*coords(l,i)
         end do
      end do
      end do
      
      call ktranspose(xjac,xjact,mcrd,mcrd)

c      do i=1,2
c            write(7,*) (xjact(j,i),j=1,2)
c      end do

      call kmult(xjact,xdX,xjactn,mcrd,mcrd,1)
c      write(7,*)'FdX(1)=',xjactn(1),'FdX(2)=',xjactn(2)
c     detjacsurf=|FdX| =dx    
      detjacsurf = SQRT(xjactn(1)**2.0d0
     +     +xjactn(2)**2.0d0)

c      write(7,*)'dx=',detjacsurf

c
c
      return
      end

c
c***********************************************************************