!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2015 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  *NEW* Analysis routine for discs by CJN
!  ***PLEASE EMAIL chris.nixon@jila.colorado.edu IF YOU WANT TO USE IT***
!  (c) 2014 Chris Nixon
!
! Method:
! "primary(secondary)" implies tag=1(2) not anything to do with mass or disc hosting etc.
! Circumprimary   disc: only particles inside the primary   Roche-lobe & bound
! Circumsecondary disc: only particles inside the secondary Roche-lobe & bound
! Circumbinary    disc: only particles at R > a (semi-major axis, not separation) & bound
!
! Creates: binary.dat      for binary params
!          circumprimary_*.dat   for circumprimary disc params
!          circumsecondary_*.dat for circumsecondary disc params
!          circumbinary_*.dat    for circumbinary disc params
!
!  REFERENCES: None
!
!  OWNER: Chris Nixon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, infile_utils, io, options, part, physcon, setbinary
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'binaryanalysis'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)

 use io,       only:fatal
 use part,     only:xyzmh_ptmass,vxyz_ptmass,ihacc
 use physcon,  only:pi
 use setbinary,only:Rochelobe_estimate
 use eos,      only:gamma
 use options,  only:ieos

 integer, parameter :: ngrid = 100
 integer, parameter :: iout  = 149
 integer, parameter :: iparams = 150
 integer, parameter :: itime = 151
 character(len=*), intent(in) :: dumpfile
 real,dimension(:,:), intent(in) :: xyzh,vxyzu
 real, intent(in) :: pmass,time
 integer, intent(in) :: npart,iunit,numfile

 character(len=50) :: output, timeoutput
 character(len=15),dimension(3) :: disc_type

 logical :: use_a = .true. 
 logical :: exists

 integer :: nsinks,i,j,ipri,isec,ibin,nunbound,icell,ierr
 integer,dimension(3) :: n
 integer,dimension(npart) :: imysink
 integer :: check
 integer :: i1,i2,i3
 integer :: irel
 integer :: ianalysis

 real,parameter :: G = 1.0
 real :: ui,ri,Ei,ai,a,ecc,dr,csi,ecci,xprir,xsecr,Limag,Hi,Ltot
 real :: psi_x,psi_y,psi_z,psi,rmax,rmin,tilt,twist,Lmagi,ai_cell,omegai
 real,dimension(3) :: xi,vi,Li,mptmass,xcom,xpricom,xseccom,rtest
 real,dimension(3) :: eccveci,vhi,nodei
 real,dimension(3) :: eccvecb,nodeb,vh,Lbin
 real,dimension(ngrid) :: rad,sigma,h_smooth,ecc_cell,ecc2_cell,Prec_cell
 real,dimension(ngrid) :: ninbin,cs_cell,omega_cell,H_cell,area,E_cell,Ephi_cell
 real,dimension(3,ngrid) :: L_cell,unitL_cell
 real,dimension(3,3) :: xptmass,vptmass
 real :: ecc1,ecc2,ecc3
 real :: ephi1,ephi2,ephi3
 real :: inc1,inc2,inc3
 real :: ie1,ie2,ie3
 real :: ibd1,ibd2,ibd3
 real :: prec1,prec2,prec3
 real :: o1,o2,o3
 real :: r1,r2, r3
 real :: masst
 real :: precveci,precvecb,eccphii,eccphib

! Use two variables solely to remove compiler warnings...
 write(*,'("Performing analysis on ",a,"... which is unit ",i5,"...")') trim(dumpfile),iunit

! read the .in file
 call read_dotin('disc.in',ieos,iparams,ierr)
! call read_dotin('planet.in',ieos,iparams,ierr)
 if (ierr /= 0) call fatal('analysis','could not open/read disc.in')

 if (gamma < tiny(gamma)) call fatal(analysistype,'gamma not set...',var='gamma',val=gamma)

! Check nsinks is >=2
 nsinks = size(xyzmh_ptmass(1,:))
 if (nsinks < 2) call fatal(analysistype,'Not enough sinks...')

! Labels for binary - in principle allows more than 2 sinks
 ipri = 1
 isec = 2
 ibin = 3

! print*,xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,2),xyzmh_ptmass(1:3,3)

! choose which disc to calulate eccentricities as a function ot time etc
 ianalysis=ipri  !ibin 

 ecc1 = -1.0
 ecc2 = -1.0
 ecc3 = -1.0
 inc1 = -1.0
 inc2 = -1.0
 inc3 = -1.0
 prec1=-1.0
 prec2=-1.0
 prec3=-1.0
 o1=-1.0
 o2=-1.0
 o3=-1.0
 ie1=-1.0
 ie2=-1.0
 ie3=-1.0
 ephi1=-1.0
 ephi2=-1.0
 ephi3=-1.0
 ibd1=-1.0
 ibd2=-1.0
 ibd3=-1.0

 write(disc_type(1),'("circumprimary")')
 write(disc_type(2),'("circumsecondary")')
 write(disc_type(3),'("circumbinary")')

! Calculate binary params (and write to binary.dat):
 call get_binary_params(ipri,isec,ibin,irel,xyzmh_ptmass,vxyz_ptmass,time,a,ecc,G,eccvecb,nodeb,precvecb,eccphib,vh,Lbin)
 mptmass(ipri) = xyzmh_ptmass(4,ipri)
 mptmass(isec) = xyzmh_ptmass(4,isec)
 mptmass(ibin) = mptmass(ipri) + mptmass(isec)

 xptmass(:,ipri) = xyzmh_ptmass(1:3,ipri)
 xptmass(:,isec) = xyzmh_ptmass(1:3,isec)
 xptmass(:,ibin) = (mptmass(ipri)*xptmass(:,ipri) + mptmass(isec)*xptmass(:,isec))/mptmass(ibin)
 vptmass(:,ipri) = vxyz_ptmass(1:3,ipri)
 vptmass(:,isec) = vxyz_ptmass(1:3,isec)
 vptmass(:,ibin) = (mptmass(ipri)*vptmass(:,ipri) + mptmass(isec)*vptmass(:,isec))/mptmass(ibin)
! location of centre of mass of binary
 xcom(:) = xptmass(:,ibin)

! Calculate Roche lobe size of primary - based on semi-major axis, not instantaneous separation
 xpricom(:) = xptmass(1:3,ipri) - xcom(:)
 xseccom(:) = xptmass(1:3,isec) - xcom(:)
 xprir = sqrt(dot_product(xpricom,xpricom)) ! distance from primary to c. o. m.
 xsecr = sqrt(dot_product(xseccom,xseccom)) ! distance from secondary to c. o. m.
! rtest(ipri) = Rochelobe_estimate(mptmass(isec),mptmass(ipri),a) ! max radius of circumprimary disc
 rtest(ipri) = 100.00D0
 rtest(isec) = Rochelobe_estimate(mptmass(ipri),mptmass(isec),a) ! max radius of circumsecondary disc
 rtest(ibin) = max(a,xprir+rtest(ipri),xsecr+rtest(isec)) ! this is an underestimate of the circumbinary inner edge...

 write(*,'("Roche-lobe of primary is: ",es17.10)') Rochelobe_estimate(mptmass(isec),mptmass(ipri),a)
 write(*,'("Roche-lobe of secondary is: ",es17.10)') Rochelobe_estimate(mptmass(ipri),mptmass(isec),a)

 write(*,'("Rejecting particles that are a distance ",es17.10, " from the primary")') rtest(ipri)
 write(*,'("Rejecting particles that are a distance ",es17.10, " from the secondary")') rtest(isec)

 n(:) = 0 !number of particles in each 'disc'
 nunbound = 0
 imysink(:) = 0

! FIRST LOOP: tag particles for each disc
 do j=1,3 ! SINK LOOP

    write(*,'("Analysing ",a15)') disc_type(j)

    sigma(:)      = 0.0
    L_cell(:,:)   = 0.0
    h_smooth(:)   = 0.0
    cs_cell(:)    = 0.0
    omega_cell(:) = 0.0
    H_cell(:)     = 0.0
    ecc_cell(:)   = 0.0
    E_cell(:)     = 0.0
    Prec_cell(:)  = 0.0

    ninbin(:) = 0

    rmin = huge(rmin)
    rmax = tiny(rmax)

    do i=1,npart ! PARTICLE LOOP
! First skip dead particles
       if (xyzh(4,i) < tiny(xyzh)) cycle

! Second check if particles in right place
       xi(1:3) = xyzh(1:3,i) - xptmass(1:3,j)

       ri = sqrt(dot_product(xi,xi))
       if (j /= ibin .and. ri > rtest(j)) cycle ! If pri/sec then check particles are inside RL
       if (j == ibin .and. ri < rtest(j)) cycle ! If bin then check particles are outside binary

! Local storage of particle velocity and thermal energy
       vi(1:3) = vxyzu(1:3,i) - vptmass(1:3,j)
       if (size(vxyzu(:,1)) == 4) then
          ui = vxyzu(4,i)
       else
          call get_utherm(ieos,xi(1),xi(2),xi(3),gamma,ui,csi)
! Use the following line if you are using eos centred on sink particles 
!          call get_utherm(ieos,xyzh(1,i),xyzh(2,i),xyzh(3,i),gamma,ui,csi)
       endif

! By this point, particles are interesting (i.e. in the right area).
! Now check they are also bound to their sink (i.e. specific energy negative):
       Ei = get_particle_energy(G,mptmass(j),ri,vi,ui)
       if (Ei > 0.0) nunbound = nunbound + 1
       if (Ei > 0.0) cycle !unbound

! Calculate particle eccentricity and semi-major axis
       call cross(xi,vi,Li)
       Limag = sqrt(dot_product(Li,Li))
       if (abs(Ei) < tiny(Ei)) stop 'initial particle energy problem'
       call get_ae(Limag,Ei,mptmass(j),pmass,ai,ecci)

       imysink(i) = j
       n(j) = n(j) + 1

! min and max of semi-major axis for grid
       if (use_a) then
          rmin = min(rmin,ai)
          rmax = max(rmax,ai)
       else
          rmin = min(rmin,ri)
          rmax = max(rmax,ri)
       endif

    enddo ! END PARTICLE LOOP
!!!!! rgm !!!!
    if(rmax.gt.100.0D0)then
       rmax=100.0D0
    endif

    rmin = 0.99*rmin
    rmax = 1.01*rmax

    if(j == ianalysis)then
       r1=40.0D0
       r2=50.0D0
       r3=70.0D0

       i1=huge(i1)
       i2=huge(i2)
       i3=huge(i3)
    endif

! Set up the radius array
    dr = (rmax-rmin)/real(ngrid-1)
    do i=1,ngrid
       rad(i)=rmin + real(i-1)*dr
       area(i) = (pi*((rad(i)+dr/2.0)**2-(rad(i)- dr/2.0)**2))

       if (j == ianalysis .and. i /= 1) then
          if (rad(i-1) < r1 .and. rad(i) > r1) i1=i
          if (rad(i-1) < r2 .and. rad(i) > r2) i2=i
          if (rad(i-1) < r3 .and. rad(i) > r3) i3=i
       endif

    enddo

!    write(*,'("Setting up grid from ",es17.10," to ",es17.10," with ",i3," points...")') rmin,rmax,ngrid

! SECOND LOOP: calculate particle properties and place into bins
    do i=1,npart ! SECOND PARTICLE LOOP
       if (imysink(i) /= j) cycle
       
       xi(1:3) = xyzh(1:3,i) - xptmass(1:3,j)
       vi(1:3) = vxyzu(1:3,i) - vptmass(1:3,j)
       ri = sqrt(dot_product(xi,xi))
       if (ri < rmin .or. ri > rmax) cycle
       if (size(vxyzu(:,1)) == 4) then
          ui = vxyzu(4,i)
          csi = sqrt(gamma*(gamma-1)*ui)
       else
!          call get_utherm(ieos,xi(1),xi(2),xi(3),gamma,ui,csi)
! Use the following line if you are using eos centred on sink particles 
           call get_utherm(ieos,xyzh(1,i),xyzh(2,i),xyzh(3,i),gamma,ui,csi)
       endif
       Ei = get_particle_energy(G,mptmass(j),ri,vi,ui)
       if (abs(Ei) < tiny(Ei)) stop 'Ei is zero...'

! Calculate particle eccentricity and semi-major axis
       call cross(xi,vi,Li)
       Limag = sqrt(dot_product(Li,Li))
!!       call get_ae(Limag,Ei,mptmass(j),pmass,ai,ecci)
       irel=ianalysis
       call get_aergm(irel,G,xi,vi,ri,Li,Limag,Ei,mptmass(j),pmass,vhi,ai,ecci,eccveci,nodei,precveci,eccphii)
!        print*,ecci/sqrt(dot_product(eccveci,eccveci)),eccveci(1),eccveci(2),eccveci(3)
!       if(precveci.gt.3.1415)then
!          print*,ecci,sqrt(dot_product(eccveci,eccveci),precveci
!       endif
       if (use_a) then
          icell = int((ai-rmin)/dr) + 1
       else
          icell = int((ri-rad(1))/dr) + 1
       endif
       if (icell > ngrid .or. icell < 1) print*, "icell = ", icell
!call fatal(analysistype,'particles where they shouldn''t be...')

! Angular momentum with pmass for later
       Li(:) = pmass*Li(:)

! Omegai and hence Hi
       if (use_a) then
          omegai = sqrt(G*mptmass(j)/ai**3)
       else
          omegai = sqrt(G*mptmass(j)/ri**3)
       endif
       Hi     = csi/omegai
		 
! Add particle to icell
       sigma(icell)      = sigma(icell)      + pmass/area(icell)
       L_cell(:,icell)   = L_cell(:,icell)   + Li(:)
       h_smooth(icell)   = h_smooth(icell)   + xyzh(4,i)
       cs_cell(icell)    = cs_cell(icell)    + csi
       omega_cell(icell) = omega_cell(icell) + omegai
       H_cell(icell)     = H_cell(icell)     + Hi
!       ecc_cell(icell)   = ecc_cell(icell)   + ecci
       ecc_cell(icell)   = ecc_cell(icell)   + sqrt(dot_product(eccveci,eccveci))
       E_cell(icell)     = E_cell(icell)     + Ei
       Ephi_cell(icell)  = Ephi_cell(icell)  + eccphii
       Prec_cell(icell)  = Prec_cell(icell)  + precveci
       ninbin(icell) = ninbin(icell) + 1
    enddo
! Average arrays that need it
    do i=1,ngrid
       if (ninbin(i) /= 0) then
          h_smooth(i)    = h_smooth(i)/ninbin(i)
          cs_cell(i)     = cs_cell(i)/ninbin(i)
          omega_cell(i)  = omega_cell(i)/ninbin(i)
          H_cell(i)      = H_cell(i)/ninbin(i)
          ecc_cell(i)    = ecc_cell(i)/ninbin(i)
          E_cell(i)      = E_cell(i)/ninbin(i)
          Prec_cell(i)   = Prec_cell(i)/ninbin(i)
          Ephi_cell(i)   = Ephi_cell(i)/ninbin(i)
       endif
    enddo
    if (j == ipri) then
       write(output,'(a13,"_",i5.5,".dat")') disc_type(j),numfile
    elseif (j == isec) then
       write(output,'(a15,"_",i5.5,".dat")') disc_type(j),numfile
    elseif (j == ibin) then
       write(output,'(a12,"_",i5.5,".dat")') disc_type(j),numfile
    else
       call fatal(analysistype,'j is some funny number...')
    endif

    open(iout,file=trim(output))
    write(iout,'("# Analysis data at t = ",es20.12)') time
    write(iout,"('#',14(1x,'[',i2.2,1x,a11,']',2x))") &
         1,'radius', &
         2,'sigma', &
         3,'<h>/H', &
         4,'lx', &
         5,'ly', &
         6,'lz', &
         7,'tilt', &
         8,'twist', &
         9,'psi', &
         10,'ecc', &
         11,'ecc2', &
         12,'H', &
!         12,'preccvec',&
!         13,'ninbin'
         13, 'i_bd',&
         14,'preccvec'
    

! Work out unitL array
    do i=1,ngrid
       if (ninbin(i) /= 0) then
          Ltot = sqrt(dot_product(L_cell(:,i),L_cell(:,i)))
          unitL_cell(:,i) = L_cell(:,i)/Ltot
       endif
    enddo

! Loop, calculate final quantities and then print to file
    do i=1,ngrid
       if (ninbin(i) /= 0) then
          if (i /= 1 .and. i /= ngrid) then
             if (ninbin(i-1) /= 0 .and. ninbin(i+1) /= 0) then
                psi_x=(unitL_cell(1,i+1)-unitL_cell(1,i-1))/(rad(i+1)-rad(i-1))
                psi_y=(unitL_cell(2,i+1)-unitL_cell(2,i-1))/(rad(i+1)-rad(i-1))
                psi_z=(unitL_cell(3,i+1)-unitL_cell(3,i-1))/(rad(i+1)-rad(i-1))
                psi=sqrt(psi_x**2 + psi_y**2 + psi_z**2)*rad(i)
             else
                psi=0.
             endif
          else
             psi=0.
          endif
          if (psi /= psi) then
             print*, 'Error with psi... at i = ',i, (rad(i+1)-rad(i-1)), (unitL_cell(3,i+1)-unitL_cell(3,i-1))
             print*, psi_x,psi_y,psi_z
             print*, rad(i-1),rad(i),rad(i+1)
             print*, unitL_cell(1,i-1),unitL_cell(2,i-1),unitL_cell(3,i-1)
             print*, unitL_cell(1,i+1),unitL_cell(2,i+1),unitL_cell(3,i+1)
             stop
          endif
       else
          unitL_cell(:,i) = 0.0
          psi = 0.0
       endif

       if (ninbin(i) > 0) then
          tilt  = acos(unitL_cell(3,i))
          twist = atan2(unitL_cell(2,i),unitL_cell(1,i))
       else
          tilt = 0.0
          twist = 0.0
       endif

       if (ninbin(i) > 0) then
          Lmagi = sqrt(dot_product(L_cell(:,i),L_cell(:,i)))/(pmass*real(ninbin(i)))
          call get_ae(Lmagi,E_cell(i),mptmass(j),pmass,ai_cell,ecc2_cell(i))
       else
          ecc2_cell(i) = 0.0
       endif

! PUT A CHECK HERE THAT ai_cell is consistent with rad(i)...

       if (ninbin(i) > 0) then
          write(iout,'(14(es18.10,1X))') rad(i),sigma(i),h_smooth(i)/H_cell(i), &
               unitL_cell(1,i),unitL_cell(2,i),unitL_cell(3,i), &
               tilt,twist,psi,ecc_cell(i),ecc2_cell(i),H_cell(i),cs_cell(i)**(2), Prec_cell(i) !real(ninbin(i))
       endif

       if (j == ianalysis .and. i1 /= huge(i1)) then
          ecc1=ecc_cell(i1)
          inc1=acos(unitL_cell(3,i1))
          prec1=Prec_cell(i1)
          o1=atan2(unitL_cell(2,i1),unitL_cell(1,i1))
!          ie1=acos(dot_product(eccvecb(:),unitL_cell(:,i1))/sqrt(dot_product(eccvecb,eccvecb)))
          ie1=acos(dot_product(eccvecb(:),unitL_cell(:,i1))/sqrt(dot_product(eccvecb,eccvecb)))
          ibd1=acos(dot_product(Lbin(:),unitL_cell(:,i1))/sqrt(dot_product(Lbin,Lbin)))
          ephi1=Ephi_cell(i1)
       endif
       if (j == ianalysis .and. i2 /= huge(i2)) then
          ecc2=ecc_cell(i2)
          inc2=acos(unitL_cell(3,i2))
          prec2=Prec_cell(i2)
          o2=atan2(unitL_cell(2,i2),unitL_cell(1,i2))
          ie2=acos(dot_product(eccvecb(:),unitL_cell(:,i2))/sqrt(dot_product(eccvecb,eccvecb)))
          ibd2=acos(dot_product(Lbin(:),unitL_cell(:,i2))/sqrt(dot_product(Lbin,Lbin)))
          ephi2=Ephi_cell(i2)
        endif
       if (j == ianalysis .and. i3 /= huge(i3)) then
          ecc3=ecc_cell(i3)
          inc3=acos(unitL_cell(3,i3))
          prec3=Prec_cell(i3)
          o3=atan2(unitL_cell(2,i3),unitL_cell(1,i3))
          ie3=acos(dot_product(eccvecb(:),unitL_cell(:,i3))/sqrt(dot_product(eccvecb,eccvecb)))
          ibd3=acos(dot_product(Lbin(:),unitL_cell(:,i3))/sqrt(dot_product(Lbin,Lbin)))
          ephi3=Ephi_cell(i3)
       endif

    enddo
    close(iout)
 enddo ! END SINK LOOP

 write(timeoutput,'("time.dat")')
 if (time <= tiny(time)) then
    open(itime,file=timeoutput,status='replace',action='write',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open time.dat file at t=0.0')
    write(itime,"('#',20(1x,'[',i2.2,1x,a11,']',2x))") &
         1,'time', &
         2,'e_1(r=40)', &
         3,'e_2(r=50)', &
         4,'e_3(r=70)', &
         5,'i_1', &
         6,'i_2', &
         7,'i_3', &
         8,'o_1', &
         9,'o_2', &
         10,'o_3', &
         11,'p_1', &
         12,'p_2',&
         13,'p_3',&
         14,'ephi_1', &
         15,'ephi_2',&
         16,'ephi_3',&
         17,'ibd_1',&
         18,'ibd_2',&
         19,'ibd_3',&
         14,'mass'
 else
    inquire(file=timeoutput,exist=exists)
    if (.not. exists) call fatal(analysistype,'t /= 0.0, but the time analysis output file does not exist...')
    open(itime,file=timeoutput,status='old',action='write',position='append',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open time.dat file during run')
 endif
!rgm calculate total disc mass 
 masst=0.0D0
 do i=1,ngrid
    masst=masst+2.0D0*3.1415D0*rad(i)*sigma(i)*dr
 enddo
!
 write(itime,'(20(ES18.10,1X))') time,ecc1,ecc2,ecc3,inc1,inc2,inc3,o1,o2,o3,prec1,prec2,prec3,ephi1,ephi2,ephi3,&
      ibd1,ibd2,ibd3,masst
! write(itime,'(14(ES18.10,1X))') time,ecc1,ecc2,ecc3,inc1,inc2,inc3,o1,o2,o3,prec1,prec2,prec3,masst
 close(itime)

 write(*,'("Number of particles in primary,secondary and binary discs: ",i9,", ",i9," and ",i9)') n(1),n(2),n(3)
 write(*,'("Number of unbound particles: ",i6)') nunbound
 write(*,'("Finished analysis.")')
 write(*,'(" ")')
 write(*,'(" ")')

end subroutine do_analysis
!----------------------------------------------------------------
!+
!  Read information from .in file
!+
!----------------------------------------------------------------
subroutine read_dotin(filename,ieos,iunit,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 implicit none
 character(len=*), intent(in)  :: filename
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ieos,ierr
 type(inopts), dimension(:), allocatable :: db

! Read in parameters from the .in file
 open(unit=iunit,file=filename,status='old',iostat=ierr,form='formatted')
 if (ierr /= 0) return

 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(ieos,'ieos',db,ierr)
 call close_db(db)
 close(iunit)

end subroutine read_dotin
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Calculate the binary paramters
! (c) 2014 Chris Nixon
!-----------------------------------------------------------------------
subroutine get_binary_params(ipri,isec,ibin,irel,xyzmh_ptmass,vxyz_ptmass,time,a,ecc,G,eccvecb,nodeb,precvecb,eccphib,vh,L)
!-----------------------------------------------------------------------
 use io, only:fatal

 implicit none

 integer,intent(in) :: ipri,isec,ibin
 real,intent(in) :: time,G
 real,dimension(:,:),intent(in) :: xyzmh_ptmass,vxyz_ptmass
 real,dimension(3),intent(out) :: L
 real,intent(out) :: a,ecc,precvecb,eccphib
 real,intent(out) :: eccvecb(3),nodeb(3),vh(3)
 integer,intent(out) :: irel

 logical :: exists
 character(len=10) :: output
 integer,parameter :: iunit = 150
 integer :: check
 real :: rbin,mpri,msec,E,Lmag
 real,dimension(3) :: xpri,vpri,xsec,vsec,dr,dv


 write(output,"(a10)") 'binary.dat'

 mpri = xyzmh_ptmass(4,ipri)
 msec = xyzmh_ptmass(4,isec)

 xpri(:) = xyzmh_ptmass(1:3,ipri)
 vpri(:) = vxyz_ptmass(1:3,ipri)
 xsec(:) = xyzmh_ptmass(1:3,isec)
 vsec(:) = vxyz_ptmass(1:3,isec)
 dr(:) = xpri(:) - xsec(:)
 dv(:) = vpri(:) - vsec(:)
 rbin  = sqrt(dot_product(dr,dr))

! Calculate the binary specific relative ang. mom and energy
 call cross(dr,dv,L)
 Lmag = sqrt(dot_product(L,L))
 E = 0.5*dot_product(dv,dv) - G*(mpri+msec)/rbin

 if (abs(E) < tiny(E)) stop 'binary energy problem'
!!! call get_ae(Lmag,E,mpri,msec,a,ecc)
 irel=ibin
 call get_aergm(irel,G,dr,dv,rbin,L,Lmag,E,mpri,msec,vh,a,ecc,eccvecb,nodeb,precvecb,eccphib)

 if (time <= tiny(time)) then
    open(iunit,file=output,status='replace',action='write',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open binary.dat file at t=0.0')
    write(iunit,"('#',8(1x,'[',i2.2,1x,a11,']',2x))") &
         1,'time', &
         2,'a', &
         3,'ecc', &
         4,'ecc_2',&
         5,'e_phi',&
         6,'e_lz',&
         7,'lb_phi',&
         8,'lb_inc'
 else
    inquire(file=output,exist=exists)
    if (.not. exists) call fatal(analysistype,'t /= 0.0, but the analysis output file does not exist...')
    open(iunit,file=output,status='old',action='write',position='append',iostat=check)
    if (check /= 0) call fatal(analysistype,'unable to open binary.dat file during run')
 endif
 write(iunit,'(8(ES18.10,1X))') time,a,ecc,sqrt(dot_product(eccvecb,eccvecb)),&
 &atan2(eccvecb(2),eccvecb(1)),acos(eccvecb(3)/sqrt(dot_product(eccvecb,eccvecb))),&
 atan2(L(2),L(1)),acos(L(3)/sqrt(dot_product(L,L)))
 close(iunit)
!-----------------------------------------------------------------------
end subroutine get_binary_params
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
subroutine get_ae(Lmag,E,m1,m2,a,ecc)
!-----------------------------------------------------------------------
! Return the semi-major axis and eccentricity between two objects
! Chris Nixon 02/10/10
!-----------------------------------------------------------------------
 use part,     only:xyzmh_ptmass,vxyz_ptmass
 implicit none
 real,intent(out) :: a,ecc
 real,intent(in) :: Lmag,E
 real :: m1,m2
 
 if (Lmag < tiny(Lmag)) stop 'Lmag is zero in get_ae'
 if (abs(E) < tiny(E)) stop 'E is zero in get_ae'

! Hence obtain the binary eccentricity
 ecc = sqrt(1.0 + (2.0*E*Lmag**2)/((m1+m2)**2))

! and semi-major axis
 a = Lmag*Lmag/((m1+m2)*(1.0-ecc*ecc))
!-----------------------------------------------------------------------
end subroutine get_ae
!-----------------------------------------------------------------------
subroutine get_aergm(irel,G,xi,vi,ri,Li,Lmag,E,m1,m2,vh,a,ecc,eccvec,node,precvec,eccphi)
!-----------------------------------------------------------------------
! Return the semi-major axis and eccentricity between two objects
! Rebecca Martin 
!-----------------------------------------------------------------------
 use part,     only:xyzmh_ptmass,vxyz_ptmass
 implicit none
 real,intent(out) :: a,ecc,eccvec(3),vh(3),node(3),precvec,eccphi
 real,intent(in) :: Lmag,E
 real,intent(in) :: Li(3),xi(3),vi(3),ri,G
 integer,intent(in) :: irel


 real :: mu,m1,m2

 if (Lmag < tiny(Lmag)) stop 'Lmag is zero in get_ae'
 if (abs(E) < tiny(E)) stop 'E is zero in get_ae'

! Hence obtain the binary eccentricity
 ecc = sqrt(1.0 + (2.0*E*Lmag**2)/((m1+m2)**2))

! and semi-major axis
 a = Lmag*Lmag/((m1+m2)*(1.0-ecc*ecc))

 if(irel.eq.3)then
    mu = G*(m1+xyzmh_ptmass(4,2))
 elseif(irel.eq.1)then
    mu=G*m1
 elseif(irel.eq.2)then
    mu=G*xyzmh_ptmass(4,2)
 else
    print*,'confused about which star to calculate eccentricity'
    stop
 endif
 call cross(vi,Li,vh)
 eccvec=vh/mu - xi/ri
! print*,ecc,sqrt(dot_product(eccvec,eccvec)),vh,mu,xi,ri
! print*,vh,ri,xi,eccvec,mu
 node(1)=-Li(2)
 node(2)=Li(1)
 node(3)=0.0

 precvec=acos(dot_product(node,eccvec)/sqrt(dot_product(node,node))/sqrt(dot_product(eccvec,eccvec)))
 eccphi=atan2(eccvec(2),eccvec(1))
! print*,ecc,sqrt(dot_product(eccvec,eccvec)),precvec
!-----------------------------------------------------------------------
end subroutine get_aergm
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
subroutine cross(a,b,c)
!-----------------------------------------------------------------------
! Return the vector cross product of two 3d vectors
! Chris Nixon 02/10/10
!-----------------------------------------------------------------------
 implicit none
 real,intent(in),dimension(3)  :: a,b
 real,intent(out),dimension(3) :: c

 c(1) = a(2)*b(3)-b(2)*a(3)
 c(2) = a(3)*b(1)-b(3)*a(1)
 c(3) = a(1)*b(2)-b(1)*a(2)

!-----------------------------------------------------------------------
end subroutine cross
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Function to calculate particle energy
!-----------------------------------------------------------------------
real(kind=8) function get_particle_energy(G,msink,ri,vi,ui)
 implicit none
 real,intent(in) :: G,msink,ri,ui
 real,dimension(3),intent(in) :: vi

 get_particle_energy = -G*msink/ri + 0.5*dot_product(vi,vi) + ui

!-----------------------------------------------------------------------
end function get_particle_energy
!-----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-----------------------------------------------------------------------
! Subroutine to return particle internal energy and sound speed
!-----------------------------------------------------------------------
subroutine get_utherm(ieos,xi,yi,zi,gamma,ui,csi)
 use eos, only:equationofstate
 implicit none

 integer,intent(in) :: ieos
 real,intent(in) :: xi,yi,zi,gamma
 real,intent(out) :: ui,csi

 real :: rhoi = 1.0 ! this is essentially a dummy variable, not needed here
!(only needed if adiabatic, but this routine is not called in that case...)
 real :: ponrhoi, tempi

 call equationofstate(ieos,ponrhoi,csi,rhoi,xi,yi,zi,tempi)

 if (gamma == 1.0) then
    ui = ponrhoi
 else
    ui = ponrhoi/(gamma-1.0)
 endif

!-----------------------------------------------------------------------
end subroutine get_utherm
!-----------------------------------------------------------------------
end module analysis

