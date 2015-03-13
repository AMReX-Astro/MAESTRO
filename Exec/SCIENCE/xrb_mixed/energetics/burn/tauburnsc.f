c   screening corrections included; xrb time scale
      program tauburn
c
c     for a given composition, temperature and density, this
c     routine calculates the inductive burning time at constant 
c     pressure
c     adjust initial comp, rho, and carbon cutoff before running

      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'

      dimension avenuc(10000),timeav(10000)

      iredo = 0


      hydrogen = 0.723d0
      helium = 0.243d0
      oxy14 = 0.004d0
      oxy15 = 0.030d0

      iprint = 6
 
 
9     snucint = 0.0d0
      timezero = 0.0d0
      dth = 1.0d-20
      time = 0.0d0
      temp = 9.5d8
      t9 = temp/1.d9
      rhozero = 5.93d5
c      rhozero = 2.0d6
      rho = rhozero

      call initburn (hydrogen,helium,oxy14,oxy15,dth)
    
c     call fxt-eos to get cold fuel pressure and energy

      call press(rhozero,t9,zbar,abar,p,dpdr,etot,cp,cv,cs)
      pzero=p
      time=timezero
      maxeset = 0 
      enucmax = 0.0d0

c     tstop is how long after max egeneration (aka the flame)
c      to stop the run. Approximately the time scale for the 
c      current burning process to be truncated by expansion

      tstop = 100.0d0
      tburn = 0.0d0

      do 100 i = 1,10000

      call dtnuc(nucleu)

      call step(xmassconv)

c  average energy generation over last nav time steps for stability
c    in nse

      nav = 5
      timeav(i) = dth
      avenuc(i) = eb*dth
      if (i.gt.nav) go to 2
      enuc = eb
      go to 5

 2    xx1 = 0.0d0
      xx2 = 0.0d0
      do 4 kk = 1,nav
      xx1 = xx1 + timeav(i-kk+1)
 4    xx2 = xx2 + avenuc(i-kk+1)
      ynav = nav
      dtq = xx1/ynav  
      enuc = xx2/xx1    

 5    enucmax = max(enucmax,enuc)
      dq=enuc*dtq

      if (iredo.eq.0) time = time + dth
      if (iredo.eq.1) time = time - dth

      temp = 1.00d9*t9
      call press(rho,t9,zbar,abar,p,dpdr,etot,cp,cv,cs)
      temp = temp+dq/cp
      t9=temp/1.d9

c     iterate to get new density at new temperature with old pressure

      r=rho
      do 30 j =1,20
      call press(r,t9,zbar,abar,p,dpdr,etot1,cp,cv,cs)
      r=r+0.5*(pzero-p)/dpdr
 30   continue
      rho=r

c      if (iredo.eq.1) snucint = snucint+enuc*dth
      snucint = snucint+enuc*dtq

      if (12.0d0*y(nc12).lt.0.6d0*x120) time1 = time1+dth

      bea = 0.0d0
      do 705 ns = 1,i2
 705  bea = bea+y(ns)*q(ns)
      bea = bea +28.296d0*aaa

      ixmax = 0
      xmax =0.0d0
      do 706 ns = 1,i2
      xa=na(ns)
      if (xa*y(ns).lt.xmax) go to 706
      xmax = xa*y(ns)
      ixmax=ns
 706  continue

      xa=na(ixmax)
      write (6,50) time,dth,t9,rho,aap,aaa*4.0d0,14.0d0*y(no14),
     1     15.0d0*y(no15),enuc,nz(ixmax),na(ixmax),xa*y(ixmax)
 50    format(1pe12.3,8e12.3,2i4,e12.3)

      tburn = tburn + dth

      if (tburn.gt.tstop) call exit(1)
      if(abs(xmassconv).gt.1.d-3) go to 801
c      if((12.0d0*y(nc12)+16.0d0*y(no16)+28.0d0*y(nsi28)).
c     1        lt.1.d-1.and.dth.gt.1.0d-2) go to 801
c      if (enuc.lt.0.8d0*enucmax.and.maxeset.eq.0
c     1     .and.enuc.gt.1.d15) go to 99 
      go to 100
 99   maxeset = 1
      tburn = 0.0d0
 100  continue

 801  continue

      abaredit = 0.0d0
      bea = 0.0d0
      do 805 ns = 1,i2
      bea = bea+y(ns)*q(ns)
 805  abaredit = abaredit+y(ns)
      abaredit = aan + aap + aaa/4.0d0 + abaredit
      abaredit = 1.0d0/abaredit
      bea = bea + aaa*28.296d0
      xxx = (bea - 7.82828d0)*1.602d-6*6.023d23
      write (6,806) abaredit,bea,snucint,xxx
 806  format (1pe12.3,3e12.3)
c      call bwrit(nucleu,time,iprint)
       if (iredo.eq.1)  call exit(1)
       timezero = time
       maxeset = 0 
       enucmax = 0.0d0
       iredo = 1
       go to 9
       
 101  call exit(1)
      end

c************************************************************


      subroutine press (rho,t9,zbar,abar,p,dpdr,etot,cp,cv,cs)

      implicit none
      save
      include 'vector_eos.dek'

c..
c..tests the helmholtz eos routine
c..
c..ionmax  = number of isotopes in the network
c..xmass   = mass fractions
c..ymass   = molar fractions
c..aion    = number of nucleons
c..zion    = number of protons

      integer          ionmax
      parameter        (ionmax=9)
      double precision t9,rho,abar,zbar,p,dpdr,etot,cp,cv,cs


c..set the input vector
      temp_row(1) = 1.0d9*t9
      den_row(1)  = rho
      abar_row(1) = abar
      zbar_row(1) = zbar

c..here the pipeline is only 1 element long
      jlo_eos = 1
      jhi_eos = 1


c..call the eos

      call helmeos


c..get p and its derivative and the energy

      p = ptot_row(1)
      dpdr = dpd_row(1)
      etot = etot_row(1)
      cp = cp_row(1)
      cv = cv_row(1)
      cs = cs_row(1)

      return
      end   




      subroutine helmeos
      implicit none
      save
      include 'vector_eos.dek'


c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999


c..declare
      double precision pi,amu,kerg,clight,avo,qe,h,ssol,asol
      parameter       (pi      = 3.1415926535897932384d0,
     1                  amu    = 1.6605402d-24,
     2                  kerg   = 1.380658d-16,
     3                  clight = 2.99792458d10, 
     4                  avo    = 6.0221367d23,
     5                  qe     = 4.8032068d-10,  
     6                  h      = 6.6260755d-27,
     7                  ssol   = 5.67051d-5,
     8                  asol   = 4.0d0 * ssol / clight)

      integer          i,j
      double precision x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida,
     1                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     2                 dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt,
     3                 deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt,
     4                 dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion,
     5                 sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd,
     6                 dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,
     7                 gam1,gam2,gam3,chit,chid,nabad,sound,etaele,
     8                 detadt,detadd,xnefer,dxnedt,dxnedd,s,
     9                 temp,den,abar,zbar,ytot1,ye,
     &                 sioncon,forth,forpi,kergavo,ikavo,asoli3,light2

      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)

c..for the abar derivatives
      double precision dpradda,deradda,dsradda,
     1                 dpionda,deionda,dsionda,
     2                 dpepda,deepda,dsepda,
     3                 dpresda,denerda,dentrda,
     4                 detada,dxneda


c..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz,
     1                 dpiondz,deiondz,dsiondz,
     2                 dpepdz,deepdz,dsepdz,
     3                 dpresdz,denerdz,dentrdz,
     4                 detadz,dxnedz


c..for the tables, in general
      integer          imax,jmax
      parameter        (imax = 211, jmax = 71)
      double precision d(imax),t(jmax)

c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)

c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),dpdfdd(imax,jmax),
     2                 dpdftt(imax,jmax),dpdfdt(imax,jmax)

c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),efdd(imax,jmax),eftt(imax,jmax),
     2                 efdt(imax,jmax)

c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),xfdd(imax,jmax),xftt(imax,jmax),
     2                 xfdt(imax,jmax)

c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     6                 ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md

c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)


c..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     1                 a1,b1,c1,d1,e1,a2,b2,c2,
     3                 ecoul,decouldd,decouldt,decoulda,decouldz,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     6                 tmelt,tfermi,rhocond,z2,x1,x2,third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third = 1.0d0/3.0d0,
     9                  esqu  = qe * qe)


c..for initialization
      integer          ifirst
      data             ifirst/0/ 


c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)

c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1

c..open the table
      open(unit=2,file='helm_table.dat',status='old')

c..read the helmholtz free energy table
       tlo   = 4.0d0
       thi   = 11.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -10.0d0
       dhi   = 11.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo

c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo

       close(unit=2)
       write(6,*)
       write(6,*) 'finished reading eos table'
       write(6,04) 'imax=',imax,' jmax=',jmax
       write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
       write(6,*)

      end if



c..start of vectorization loop, normal executaion starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = ytot1 * zbar

c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


c..radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0


       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0

c..ion section:
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnida  = -xni * ytot1

       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
       dpionda = dxnida * kt 
       dpiondz = 0.0d0

       eion    = 1.5d0 * pion*deni
       deiondd = (1.5d0 * dpiondd - eion)*deni
       deiondt = 1.5d0 * dpiondt*deni
       deionda = 1.5d0 * dpionda*deni
       deiondz = 0.0d0
    

       x       = abar*abar*sqrt(abar) * deni/avo
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)
       sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi
     1            - kergavo * deni * ytot1
       dsiondt = (dpiondt*deni + deiondt)*tempi - 
     1           (pion*deni + eion) * tempi*tempi 
     2           + 1.5d0 * kergavo * tempi*ytot1
       x       = avo*kerg/abar
       dsionda = (dpionda*deni + deionda)*tempi 
     1           + kergavo*ytot1*ytot1* (2.5d0 - y)
       dsiondz = 0.0d0


c..electron-positron section:
c..assume complete ionization 
       xnem    = xni * zbar

c..enter the table with ye*den
       din = ye*den

c..bomb proof the input
       if (temp .gt. t(jmax)) then
        write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
        write(6,*) 'temp too hot, off grid'       
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (temp .lt. t(1)) then
        write(6,01) 'temp=',temp,' t(1)=',t(1)
        write(6,*) 'temp too cold, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .gt. d(imax)) then
        write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
        write(6,*) 'ye*den too big, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .lt. d(1)) then
        write(6,01) 'ye*den=',din,' d(1)=',d(1)
        write(6,*) 'ye*den too small, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if

c..hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))


c..access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)


c..various differences
       xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

c..the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)

       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)

c..derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)

       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)

       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)

       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)

c..second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)

c       ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c       ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c       ddsi2d =   ddpsi2(xd)

c       ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c       ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c       ddsi2md =  ddpsi2(mxd)


c..the free energy
       free  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density
       df_d  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..derivative with respect to temperature
       df_t = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density**2
c       df_dd = h5(iat,jat,
c     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2         ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

c..derivative with respect to temperature**2
       df_tt = h5(iat,jat,
     1       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to temperature and density
       df_dt = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



c..now get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)

       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)

       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)

       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)


c..derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)


c..look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

c..pressure derivative with density
       dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(ye * dpepdd,0.0d0)



c..look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)


c..electron chemical potential etaele
       etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
       x       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = ye * x

c..derivative with respect to temperature
       detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      detada = -x * din * ytot1
      detadz =  x * den * ytot1



c..look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

c..electron + positron number densities
      xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
      x        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
      x = max(x,0.0d0)
      dxnedd   = ye * x

c..derivative with respect to temperature
      dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
      dxneda = -x * din * ytot1
      dxnedz =  x  * den * ytot1
       


c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of this one
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt
c       dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/ye - 2.0d0 * din * df_d
       dpepda  = -ytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x
       dsepda  = ytot1 * (ye * df_dt * din - sele)
       dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd
       deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
       deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




c..coulomb section:
c..initialize


        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        dpcoulda = 0.0d0
        dpcouldz = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        decoulda = 0.0d0
        decouldz = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0
        dscoulda = 0.0d0
        dscouldz = 0.0d0


c..uniform background corrections only 
c..from yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter
        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami 
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
         if (plasg .ge. 1.0) then
          x        = plasg**(0.25d0) 
          y        = avo * ytot1 * kerg 
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = third * den * ecoul
          scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1              + d1 * (log(plasg) - 1.0d0) - e1)

          y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
          decouldd = y * plasgdd 
          decouldt = y * plasgdt + ecoul/temp
          decoulda = y * plasgda - ecoul/abar
          decouldz = y * plasgdz

          y        = third * den
          dpcouldd = third * ecoul + y*decouldd
          dpcouldt = y * decouldt
          dpcoulda = y * decoulda
          dpcouldz = y * decouldz


          y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
          dscouldd = y * plasgdd
          dscouldt = y * plasgdt
          dscoulda = y * plasgda - scoul/abar
          dscouldz = y * plasgdz


c..yakovlev & shalybkov 1989 equations 102, 103, 104
         else if (plasg .lt. 1.0) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - third * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0d0 * pcoul/den
          scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

          s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
          dpcouldd = -dpiondd*z - pion*s*plasgdd
          dpcouldt = -dpiondt*z - pion*s*plasgdt
          dpcoulda = -dpionda*z - pion*s*plasgda
          dpcouldz = -dpiondz*z - pion*s*plasgdz

          s        = 3.0d0/den
          decouldd = s * dpcouldd - ecoul/den
          decouldt = s * dpcouldt
          decoulda = s * dpcoulda
          decouldz = s * dpcouldz

          s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
          dscouldd = s * plasgdd
          dscouldt = s * plasgdt
          dscoulda = s * plasgda - scoul/abar
          dscouldz = s * plasgdz
         end if



c..bomb proof
        x   = prad + pion + pele + pcoul
        if (x .le. 0.0) then

c         write(6,*) 
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*) 

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if





c..sum all the components
       pres    = prad + pion + pele + pcoul
       ener    = erad + eion + eele + ecoul
       entr    = srad + sion + sele + scoul

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd 
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz

       denerdd = deraddd + deiondd + deepdd + decouldd
       denerdt = deraddt + deiondt + deepdt + decouldt
       denerda = deradda + deionda + deepda + decoulda
       denerdz = deraddz + deiondz + deepdz + decouldz

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
       dentrda = dsradda + dsionda + dsepda + dscoulda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz


c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)


c..maxwell relations; each is zero if the consistency is perfect
       x   = den * den
       dse = temp*dentrdt/denerdt - 1.0d0
       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
       dsp = -dentrdd*x/dpresdt - 1.0d0


c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda   
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda   
        dez_row(j)    = denerdz

        stot_row(j)   = entr 
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda        
        dsz_row(j)    = dentrdz

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad 

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion 
        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda  
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda   
        deepz_row(j)  = deepdz

        sele_row(j)   = sele 
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt 
        dsepd_row(j)  = dsepdd 
        dsepa_row(j)  = dsepda        
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul 
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound

c..end of vectorization loop
      enddo
      return
      end

c************** start burn routines ***********************

      subroutine initburn (hydrogen,helium,oxy14,oxy15,dtime)

      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'

c..
c..   this routine initializes the burn run
c..   call once at beginning of problem
c..

      common/ntwk/ iz(100),inmin(100),inmax(100),iipass
c..
      common x1(nburn),nx1(nburn),xz(286)
      dimension sol(286),over(286),xmax(286)
      integer izsol(286),iasol(286),jcode(286),
     1       iprogn(286), iprogz(286), iproga(286)
c..
c..
c..data statements begin here
      data ifirst/0/, initsc/0/

      character*2 isymb(0:98)

      data isymb/
     &     ' n',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     1     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc',
     2     'Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
     3     'As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc',
     4     'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe',
     5     'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb',
     6     'Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os',
     7     'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr',
     8     'Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf'/

c.. mass fractions anders & grevese 1989

        data (sol(i),i=1,45)/
     1 7.0573E-01, 4.8010E-05, 2.9291E-05, 2.7521E-01, 6.4957E-10,
     2 9.3490E-09, 1.6619E-10, 1.0674E-09, 4.7301E-09, 3.0324E-03,
     3 3.6501E-05, 1.1049E-03, 4.3634E-06, 9.5918E-03, 3.8873E-06,
     4 2.1673E-05, 4.0515E-07, 1.6189E-03, 4.1274E-06, 1.3022E-04,
     5 3.3394E-05, 5.1480E-04, 6.7664E-05, 7.7605E-05, 5.8052E-05,
     6 6.5301E-04, 3.4257E-05, 2.3524E-05, 8.1551E-06, 3.9581E-04,
     7 3.2221E-06, 1.8663E-05, 9.3793E-08, 2.5320E-06, 8.5449E-07,
     8 7.7402E-05, 1.5379E-05, 2.6307E-08, 3.4725E-06, 4.4519E-10,
     9 2.6342E-07, 5.9898E-05, 4.1964E-07, 8.9734E-07, 1.4135E-06/

        data (sol(i),i=46,90)/
     1 2.7926E-09, 1.3841E-07, 3.8929E-08, 2.2340E-07, 2.0805E-07,
     2 2.1491E-06, 1.6361E-07, 1.6442E-07, 9.2579E-10, 3.7669E-07,
     3 7.4240E-07, 1.4863E-05, 1.7160E-06, 4.3573E-07, 1.3286E-05,
     4 7.1301E-05, 1.1686E-03, 2.8548E-05, 3.6971E-06, 3.3579E-06,
     5 4.9441E-05, 1.9578E-05, 8.5944E-07, 2.7759E-06, 7.2687E-07,
     6 5.7528E-07, 2.6471E-07, 9.9237E-07, 5.8765E-07, 8.7619E-08,
     7 4.0593E-07, 1.3811E-08, 3.9619E-08, 2.7119E-08, 4.3204E-08,
     8 5.9372E-08, 1.7136E-08, 8.1237E-08, 1.7840E-08, 1.2445E-08,
     9 1.0295E-09, 1.0766E-08, 9.1542E-09, 2.9003E-08, 6.2529E-08/

        data (sol(i),i=91,135)/
     1 1.1823E-08, 1.1950E-08, 1.2006E-08, 3.0187E-10, 2.0216E-09,
     2 1.0682E-08, 1.0833E-08, 5.4607E-08, 1.7055E-08, 1.1008E-08,
     3 4.3353E-09, 2.8047E-10, 5.0468E-09, 3.6091E-09, 4.3183E-08,
     4 1.0446E-08, 1.3363E-08, 2.9463E-09, 4.5612E-09, 4.7079E-09,
     5 7.7706E-10, 1.6420E-09, 8.7966E-10, 5.6114E-10, 9.7562E-10,
     6 1.0320E-09, 5.9868E-10, 1.5245E-09, 6.2225E-10, 2.5012E-10,
     7 8.6761E-11, 5.9099E-10, 5.9190E-10, 8.0731E-10, 1.5171E-09,
     8 9.1547E-10, 8.9625E-10, 3.6637E-11, 4.0775E-10, 8.2335E-10,
     9 1.0189E-09, 1.0053E-09, 4.5354E-10, 6.8205E-10, 6.4517E-10/

        data (sol(i),i=136,180)/
     1 5.3893E-11, 3.9065E-11, 5.5927E-10, 5.7839E-10, 1.0992E-09,
     2 5.6309E-10, 1.3351E-09, 3.5504E-10, 2.2581E-11, 5.1197E-10,
     3 1.0539E-10, 7.1802E-11, 3.9852E-11, 1.6285E-09, 8.6713E-10,
     4 2.7609E-09, 9.8731E-10, 3.7639E-09, 5.4622E-10, 6.9318E-10,
     5 5.4174E-10, 4.1069E-10, 1.3052E-11, 3.8266E-10, 1.3316E-10,
     6 7.1827E-10, 1.0814E-09, 3.1553E-09, 4.9538E-09, 5.3600E-09,
     7 2.8912E-09, 1.7910E-11, 1.6223E-11, 3.3349E-10, 4.1767E-09,
     8 6.7411E-10, 3.3799E-09, 4.1403E-09, 1.5558E-09, 1.2832E-09,
     9 1.2515E-09, 1.5652E-11, 1.5125E-11, 3.6946E-10, 1.0108E-09/

        data (sol(i),i=181,225)/
     1 1.2144E-09, 1.7466E-09, 1.1240E-08, 1.3858E-12, 1.5681E-09,
     2 7.4306E-12, 9.9136E-12, 3.5767E-09, 4.5258E-10, 5.9562E-10,
     3 8.0817E-10, 3.6533E-10, 7.1757E-10, 2.5198E-10, 5.2441E-10,
     4 1.7857E-10, 1.7719E-10, 2.9140E-11, 1.4390E-10, 1.0931E-10,
     5 1.3417E-10, 7.2470E-11, 2.6491E-10, 2.2827E-10, 1.7761E-10,
     6 1.9660E-10, 2.5376E-12, 2.8008E-11, 1.9133E-10, 2.6675E-10,
     7 2.0492E-10, 3.2772E-10, 2.9180E-10, 2.8274E-10, 8.6812E-13,
     8 1.4787E-12, 3.7315E-11, 3.0340E-10, 4.1387E-10, 4.0489E-10,
     9 4.6047E-10, 3.7104E-10, 1.4342E-12, 1.6759E-11, 3.5397E-10/

        data (sol(i),i=226,270)/
     1 2.4332E-10, 2.8557E-10, 1.6082E-10, 1.6159E-10, 1.3599E-12,
     2 3.2509E-11, 1.5312E-10, 2.3624E-10, 1.7504E-10, 3.4682E-10,
     3 1.4023E-10, 1.5803E-10, 4.2293E-12, 1.0783E-12, 3.4992E-11,
     4 1.2581E-10, 1.8550E-10, 9.3272E-11, 2.4131E-10, 1.1292E-14,
     5 9.4772E-11, 7.8768E-13, 1.6113E-10, 8.7950E-11, 1.8989E-10,
     6 1.7878E-10, 9.0315E-11, 1.5326E-10, 5.6782E-13, 5.0342E-11,
     7 5.1086E-11, 4.2704E-10, 5.2110E-10, 8.5547E-10, 1.3453E-09,
     8 1.1933E-09, 2.0211E-09, 8.1702E-13, 5.0994E-11, 2.1641E-09,
     9 2.2344E-09, 1.6757E-09, 4.8231E-10, 9.3184E-10, 2.3797E-12/

        data (sol(i),i=271,286)/
     1 1.7079E-10, 2.8843E-10, 3.9764E-10, 2.2828E-10, 5.1607E-10,
     2 1.2023E-10, 2.7882E-10, 6.7411E-10, 3.1529E-10, 3.1369E-09,
     3 3.4034E-09, 9.6809E-09, 7.6127E-10, 1.9659E-10, 3.8519E-13,
     4 5.3760E-11/


c.. izsol and iasol are arrays of z,a of stable isotopes H-Pu238

        data (izsol(i),i=1,117)/
     1   1,   1,   2,   2,   3,   3,   4,   5,   5,   6,   6,   7,   7,
     2   8,   8,   8,   9,  10,  10,  10,  11,  12,  12,  12,  13,  14,
     3  14,  14,  15,  16,  16,  16,  16,  17,  17,  18,  18,  18,  19,
     4  19,  19,  20,  20,  20,  20,  20,  20,  21,  22,  22,  22,  22,
     5  22,  23,  23,  24,  24,  24,  24,  25,  26,  26,  26,  26,  27,
     6  28,  28,  28,  28,  28,  29,  29,  30,  30,  30,  30,  30,  31,
     7  31,  32,  32,  32,  32,  32,  33,  34,  34,  34,  34,  34,  34,
     8  35,  35,  36,  36,  36,  36,  36,  36,  37,  37,  38,  38,  38,
     9  38,  39,  40,  40,  40,  40,  40,  41,  42,  42,  42,  42,  42/
 
        data (izsol(i),i=118,234)/
     1  42,  42,  44,  44,  44,  44,  44,  44,  44,  45,  46,  46,  46,
     2  46,  46,  46,  47,  47,  48,  48,  48,  48,  48,  48,  48,  48,
     3  49,  49,  50,  50,  50,  50,  50,  50,  50,  50,  50,  50,  51,
     4  51,  52,  52,  52,  52,  52,  52,  52,  52,  53,  54,  54,  54,
     5  54,  54,  54,  54,  54,  54,  55,  56,  56,  56,  56,  56,  56,
     6  56,  57,  57,  58,  58,  58,  58,  59,  60,  60,  60,  60,  60,
     7  60,  60,  62,  62,  62,  62,  62,  62,  62,  63,  63,  64,  64,
     8  64,  64,  64,  64,  64,  65,  66,  66,  66,  66,  66,  66,  66,
     9  67,  68,  68,  68,  68,  68,  68,  69,  70,  70,  70,  70,  70/

        data (izsol(i),i=235,286)/
     1  70,  70,  71,  71,  72,  72,  72,  72,  72,  72,  73,  73,  74,
     2  74,  74,  74,  74,  75,  75,  76,  76,  76,  76,  76,  76,  76,
     3  77,  77,  78,  78,  78,  78,  78,  78,  79,  80,  80,  80,  80,
     4  80,  80,  80,  81,  81,  82,  82,  82,  82,  83,  90,  92,  92/

c.. iasol

        data (iasol(i),i=1,117)/
     1   1,   2,   3,   4,   6,   7,   9,  10,  11,  12,  13,  14,  15,
     2  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
     3  29,  30,  31,  32,  33,  34,  36,  35,  37,  36,  38,  40,  39,
     4  40,  41,  40,  42,  43,  44,  46,  48,  45,  46,  47,  48,  49,
     5  50,  50,  51,  50,  52,  53,  54,  55,  54,  56,  57,  58,  59,
     6  58,  60,  61,  62,  64,  63,  65,  64,  66,  67,  68,  70,  69,
     7  71,  70,  72,  73,  74,  76,  75,  74,  76,  77,  78,  80,  82,
     8  79,  81,  78,  80,  82,  83,  84,  86,  85,  87,  84,  86,  87,
     9  88,  89,  90,  91,  92,  94,  96,  93,  92,  94,  95,  96,  97/

        data (iasol(i),i=118,234)/
     1  98, 100,  96,  98,  99, 100, 101, 102, 104, 103, 102, 104, 105,
     2 106, 108, 110, 107, 109, 106, 108, 110, 111, 112, 113, 114, 116,
     3 113, 115, 112, 114, 115, 116, 117, 118, 119, 120, 122, 124, 121,
     4 123, 120, 122, 123, 124, 125, 126, 128, 130, 127, 124, 126, 128,
     5 129, 130, 131, 132, 134, 136, 133, 130, 132, 134, 135, 136, 137,
     6 138, 138, 139, 136, 138, 140, 142, 141, 142, 143, 144, 145, 146,
     7 148, 150, 144, 147, 148, 149, 150, 152, 154, 151, 153, 152, 154,
     8 155, 156, 157, 158, 160, 159, 156, 158, 160, 161, 162, 163, 164,
     9 165, 162, 164, 166, 167, 168, 170, 169, 168, 170, 171, 172, 173/

        data (iasol(i),i=235,286)/
     1 174, 176, 175, 176, 174, 176, 177, 178, 179, 180, 180, 181, 180,
     2 182, 183, 184, 186, 185, 187, 184, 186, 187, 188, 189, 190, 192,
     3 191, 193, 190, 192, 194, 195, 196, 198, 197, 196, 198, 199, 200,
     4 201, 202, 204, 203, 205, 204, 206, 207, 208, 209, 232, 235, 238/

c.. jcode is an array that tells the type of radioactive progenitors
c.. each stable species can have
c.. jcode = 0 of the species is the only stable one of that a
c..       = 1 if the species can have only proton-rich progenitors
c..       = 2 if the species can have only neutron-rich progenitors
c..       = 3 if the species can only be made as itself (eg k40)

        data (jcode(i),i=1,117)/
     1   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     2   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     3   0,   0,   0,   0,   0,   0,   2,   0,   0,   1,   0,   2,   0,
     4   3,   0,   1,   0,   0,   0,   2,   2,   0,   1,   0,   1,   0,
     5   2,   3,   0,   1,   0,   0,   2,   0,   1,   0,   0,   2,   0,
     6   1,   0,   0,   0,   2,   0,   0,   1,   0,   0,   0,   2,   0,
     7   0,   1,   0,   0,   2,   2,   0,   1,   1,   0,   2,   2,   2,
     8   0,   0,   1,   1,   1,   0,   2,   2,   0,   2,   1,   1,   1,
     9   0,   0,   0,   0,   2,   2,   2,   0,   1,   1,   0,   3,   0/
 
        data (jcode(i),i=118,234)/
     1   2,   2,   1,   1,   0,   1,   0,   2,   2,   0,   1,   1,   0,
     2   2,   2,   2,   0,   0,   1,   1,   1,   0,   2,   2,   2,   2,
     3   1,   2,   1,   1,   1,   1,   0,   0,   0,   2,   2,   2,   0,
     4   2,   1,   1,   1,   3,   0,   2,   2,   2,   0,   1,   1,   1,
     5   0,   3,   0,   2,   2,   2,   0,   1,   1,   1,   0,   3,   0,
     6   2,   3,   0,   1,   1,   0,   2,   0,   1,   0,   2,   0,   0,
     7   2,   2,   1,   0,   1,   0,   1,   2,   2,   0,   0,   1,   1,
     8   0,   2,   0,   2,   2,   0,   1,   1,   1,   0,   2,   0,   2,
     9   0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0,   0/
 
        data (jcode(i),i=235,286)/
     1   2,   2,   0,   3,   1,   1,   0,   0,   0,   2,   3,   0,   1,
     2   0,   0,   2,   2,   0,   2,   1,   1,   1,   0,   0,   2,   2,
     3   0,   0,   1,   1,   0,   0,   2,   2,   0,   1,   1,   0,   0,
     4   0,   0,   2,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0/

c..  
c..setrat reads rate parameters, z, n, and binding energies
c..

      open(unit=5,file='rath00_8.1.bdat',
     1     status='old',err=1000)
      go to 1001

1000  write(6,*) 'bdat file not found:','rath00_8.1.bdat'
      call exit(1)

1001  open(unit=4,file='network',status='old',err=1002)
      go to 1003

1002  write(6,*) 'network file not found:','network'
      call exit(1)

1003  call setrat

      imax2=nburn-3
      if (i2.le.imax2) go to 5
      write (6,7)
    7 format (//,'error in dimension. halt.')
      call exit(1)

    5 if (kmax-ncsize) 30,30,10
   10 write (6,20) kmax
   20 format (1x,' error in dimension of c. kmax=',i5)
      call exit(1)

 30   write (6,*) 'network initialized',nppf

      itot=i2+3
      do 55 j=1,itot
 55   y(j)=0.0d0
      aaa=helium/4.0d0
      aap=hydrogen

      y(no14) = oxy14/14.0d0
      y(no15) = oxy15/15.0d0
      write(6,*) no14,no15,oxy14
      fdtn = 2.0d0
      chimin = 1.0d-7
      delchi = 0.15d0
      mazful = 0
c..  
c..   generate average zbar, z2bar, abar, ye for screen
c..  
      i3=i2
      abm1=0.

      do 13 il=1,i2
 13   abm1=abm1+y(il)
      abm1=abm1+aap+aan+aaa
      abar=1./abm1
      ye=0.
      do 14 il=1,i2
      xz(il) = nz(il)
 14   ye=ye+xz(il)*y(il)
      ye=ye+aap+2.0d0*aaa
      zbar=abar*ye
      z2=0.
      do 16 il=1,i2
      xz(il)=xz(il)*xz(il)*y(il)
 16   z2=z2+xz(il)
      z2=z2+aap+4.0d0*aaa
      z2bar=z2*abar

      return 
      end



      subroutine step(xtotm1) 

      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'
      dimension xz(nburn), x1(nburn), nx1(nburn)

      itot = i2+3
c
c..   begin evolutionary loop
c

c..  
c..   generate average zbar, z2bar, abar, ye for screen
c..  
      i3=i2+3
      abm1=0.

      do 13 il=1,i3
 13   abm1=abm1+y(il)
      abar=1./abm1
      ye=0.
      do 14 il=1,i3
      xz(il) = nz(il)
 14   ye=ye+xz(il)*y(il)
      zbar=abar*ye
      z2=0.
      do 16 il=1,i3
      xz(il)=xz(il)*xz(il)*y(il)
 16   z2=z2+xz(il)
      z2bar=z2*abar

c..   screen computes screening corrections to nuclear reaction rates
c..   store arrays of correction factors for p and alpha interactions
c..   according to the z of the target nucleus in the mass
c..   increasing reaction.  screening for 3a, 12c +12c, 12c + 16o,
c..   and 16o + 16o are stored in sc3a, sc1212, sc1216,
c..   and sc1616 respectively.
c..  
125   tt=1.0e+09*t9
      iscrn = 0
      do 101 i=1,i2
       z1=nz(i)
       a1=na(i)
       iscrn = iscrn + 1
       call screen (tt,rho,z1,1.d0,a1,1.d0,scfp,iscrn,initsc)
       iscrn = iscrn + 1
       call screen (tt,rho,z1,2.d0,a1,4.d0,scfa,iscrn,initsc)
       scfacp(i)=scfp
       scfaca(i)=scfa
c       scfacp(i) = 1.0d0
c       scfaca(i) = 1.0d0
101   continue
      iscrn = iscrn + 1
      call screen (tt,rho,2.d0,2.d0,4.d0,4.d0,sca1,iscrn,initsc)
      iscrn = iscrn + 1
      call screen (tt,rho,4.d0,2.d0,8.d0,4.d0,sca2,iscrn,initsc)
      sc3a=sca1*sca2
c      sc3a = 1.0d0
      iscrn = iscrn + 1
      call screen (tt,rho,6.d0,6.d0,12.d0,12.d0,sc1212,iscrn,initsc)
c      sc1212 = 1.0d0
      iscrn = iscrn + 1
      call screen (tt,rho,6.d0,8.d0,12.d0,16.d0,sc1216,iscrn,initsc)
c      sc1216 = 1.0d0
      iscrn = iscrn + 1
      call screen (tt,rho,8.d0,8.d0,16.d0,16.d0,sc1616,iscrn,initsc)
c      sc1616 = 1.0d0
c..
c..reset the screen initialization flag
      initsc = 1
c..
c..   rate determines forward and inverse rates for general reactions
c..   spcrat and specl do the same for special reactions

      if (c12agm.eq.0.0d0) c12agm = 1.0d0
      call rate
      call spcrat
      call specl(c12agm)

c..
c.. no neutrino processing
      inu = 0
c..
c..tag all the nonzero locations; copy them; write some information

9800  if (ifirst .eq. 0) then
       ifirst = 1
       call builda(inu)
       do 128 j=1,nzo
        ivect(j) = irow(j)
        jvect(j) = icol(j)
128    continue
       open(unit=16,file='jac.dat',status='unknown')
       write(16,*) ' '
       write(16,*) nzo,' nonzero matrix locations'
       write(16,*) nterm,' jacobian contributions'
       write(16,*) ' '
       close(unit=16)
c..
c..force the diagonal to be the pivot; the symbolic decomposition
       do 129 j=1,nzo
        a(j) = 1.0e-10
        if (ivect(j) .eq. jvect(j)) a(j) = 1.0
129    continue
       call ma28ad(itot,nzo,a,naij,irow,naij,icol,
     1             upiv,ikeep,iw,wrk,iflag)
       if (iflag .lt. 0) then
        write(6,*) 'error flag from ma28ad ',iflag
        stop 'error from ma28ad'
       end if
      end if
c..
c..compt forms the right hand side and the negative jacobian of
c..the reaction equations for both general and special reactions

      call compt(inu)
c..
c..form the a matrix, a = 1/dt - jacobian
      rdt=1.0/dth
      do 130 j=1,nzo
       if (ivect(j) .eq. jvect(j)) a(j) = a(j) + rdt
130   continue
c..
c..the numerical decomposition
      call ma28bd(itot,nzo,a,naij,ivect,jvect,icol,ikeep,iw,wrk,iflag)
      if (flag .lt. 0) then
       write(6,*) 'error flag from ma28bd ',flag
       stop 'error from ma28bd'
      end if
c..
c..backsubstitution 

      call ma28cd(itot,a,naij,icol,ikeep,b,wrk,1)

c
c..   update particle densities
c
      do 140  j=1,i2
      y(j)=y(j)+b(j)
  140 if (y(j).lt.1.0d-99) y(j)=0.0
      aan=aan+b(i2+1)
      aap=aap+b(i2+2)
      aaa=aaa+b(itot)
      y(i2+1)=aan
      y(i2+2)=aap
      y(i2+3)=aaa

c..  check to see if any of these went negative, if so, stop
c      if (y(i2+1).lt.0.0.or.y(i2+2).lt.0.0.or.y(i2+3).lt.0.0) then
c         call bwrit(nucleu,tim,n6)
c         write(6,141) y(i2+1),y(i2+2),y(i2+3)*4.0
c141      format(2x,'stop! x(n)= ',1pe12.5,' x(p)= ',e12.5,' x(a)= ',
c     1          e12.5)
c         stop
c      end if
c

c..   calculate energy generation rate
c
      fak=b(itot)
      eb=28.295*fak
      dn=0.
      do 160  j=1,i2
      eb=eb+q(j)*b(j)
  160 dn=dn+b(j)
      eb=9.65e+17*eb/dth
      s=eb
c..  
c..   calculate etan and mass fractions
c..  
      etan=0.
      call qjsub0 (nx1,nn,nz,i2)
      call qvflot (x1,nx1,i2)
      call qvdot (etan,x1,y,i2)
      call qvflot (x1,na,i2)
      call qvmpy0 (x1,x1,y,i2)
      etan=etan+y(i2+1)-y(i2+2)

c.. sum the mass fractions
      xtot=0.
      do i=1,i2
         xtot=xtot+x1(i)
      end do
      xtotm1=1.-(xtot+y(i2+1)+y(i2+2)+4.0*y(i2+3))
 
      return
      end


      subroutine naray
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'
c..  
c..   subroutine builds the integer array nrr(7,i) which specifies the
c..   network location of elements coupled to i by various reactions.
c..  
c..   1 = ng          4 = ap
c..   2 = pn          5 = an
c..   3 = pg          6 = ag
c..                   7 = b-
      common jz(7),jn(7)
      jz(1)=0
      jz(2)=1
      jz(3)=1
      jz(4)=1
      jz(5)=2
      jz(6)=2
      jz(7)=-1
      jn(1)=1
      jn(2)=-1
      jn(3)=0
      jn(4)=2
      jn(5)=1
      jn(6)=2
      jn(7)=1
      idex=7*nburn
      call qvseti (nz0,nrreq,idex)
      do 100 i=1,i2
      do 100 n=1,7
      kz=nz(i)+jz(n)
      kn=nn(i)+jn(n)
      do 100 k=1,i2
      if(kz.ne.nz(k))go to 100
      if(kn.ne.nn(k))go to 100
      nrr(n,i)=k
  100 continue

c..   modified ix/28/88 to build location array for neutrino
c..   interactions as well
c..  
c..   1 = nu, e-,n               4 = nu  e+,n (like gp)
c..   2 = nu, e-   (like pn)     5 = nu, e+   (like np)
c..   3 = nu  e-,p (like gn)     6 = nu, e+,p
c..  

      jz(1)=1
      jz(2)=1
      jz(3)=0
      jz(4)=-1
      jz(5)=-1
      jz(6)=-2
      jz(7)=-2
      jn(1)=-2
      jn(2)=-1
      jn(3)=-1
      jn(4)=0
      jn(5)=1
      jn(6)=1
      jn(7)=-2
      do 200 i=1,i2
      do 200 n=1,7
      nrrneut(n,i)=0.
      kz=nz(i)+jz(n)
      kn=nn(i)+jn(n)
      do 200 k=1,i2
      if(kz.ne.nz(k))go to 200
      if(kn.ne.nn(k))go to 200
      nrrneut(n,i)=k
  200 continue

      return
      end
c..
c..
c..
      subroutine setrat
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'

c..
c..communicate
      common/ntwk/ iz(100),inmin(100),inmax(100),iipass
      dimension dum(100)
  
c..  
c..  
      nh2 = 0
      nh3 = 0
      nhe3 = 0
      nli6 = 0
      nli7 = 0
      nbe7 = 0
      nbe9=0
      nb8 = 0
      nb10 = 0
      nb11 = 0
      nc11 = 0
      nc12 = 0
      nc13 = 0
      nc14 = 0
      n13 = 0
      n14 = 0
      n15 = 0
      no14 = 0
      no15 = 0
      no16 = 0
      no17 = 0
      no18 = 0
      nf17 = 0
      nf18 = 0
      nf19 = 0
      ne19 = 0
      ne20 = 0
      ne21 = 0
      ne22 = 0
      na21 = 0
      na22 = 0
      na23 = 0
      mg23 = 0
      mg24 = 0
      mg25 = 0
      mg26 = 0
      nal25 = 0
      nal26 = 0
      nal27 = 0
      nsi27 = 0
      nsi28 = 0
      nsi29 = 0
      nsi30 = 0
      np30 = 0
      np31 = 0
      ns31 = 0
      nful=0
      nfulnot=0
c..  
c..   read in network to be employed
c..  
      ii=1
c..   protection for overwriting arrays
      if (i.gt.nburn-3) then
         write(6,*) ' i.ge.(nburn-3), stop'
         call exit(1)
      end if
    1 read (4,2,end=3)  iz(ii),inmin(ii),inmax(ii)
    2 format (3i5)
      if (iz(ii).eq.99) go to 3
      ii=ii+1
      go to 1
    3 ii=ii-1
      iipass=ii
c..  
c..   begin loop to read in cross-section data
c..  

      k=1
      i=0
      istart = 0
   10 i=i+1
c..  
c..   read first card of element parameter deck
c..  
c..   i is the code number of element z(i),n(i).
c..  
c..   j = 1 for ng              j = 6 for an
c..   j = 2 for pn              j = 7 for ag
c..   j = 3 for ground state    j = 8 for semi-empirical electron captur
c..         beta decay          j = 9 for semi-empirical positron decay
c..   j = 4 for pg              j = 10 for semi-empirical beta decay
c..   j = 5 for ap
c..  
c..   ic1(j,i) = type formula to be used to calculate rate
c..  
c..   ic2(j,i) = number of constants in fitting formula
c..            for reaction j on species i
c..  

      dum(1)=0.d0
      read (5,20,err=1321)  nz(i),na(i),(ic1(j,i),ic2(j,i), j=1,10)
   20 format (2i6,20i3)
      go to 1323
1321  write(43,1322) nz(i),na(i),(ic1(j,i),ic2(j,i), j=1,10),i,j
1322  format(i7)
      call exit(1)
1323  continue


c.... now first line in bdat indicates either 4 or 5 prm fit to pf
c.... if nz(i).and.na(i).both.lt.0 then number of prms is 5
c.... read next line to start the parameter read

      if (i.eq.1.and.istart.eq.0) then
       if (nz(i).gt.0.and.na(i).gt.0) then
        nppf=5
        istart = 1
       else
        nppf=6
        istart = 1
        read (5,20)  nz(i),na(i),(ic1(j,i),ic2(j,i), j=1,10)
       end if
      end if

      if (nz(i).eq.99) go to 60
c..  
c..   read temperature dependent partition function information
c..  
      llmin=nppf*(i-1)+1
      llmax=llmin+(nppf-1)

      if (nppf.eq.5) read (5,777)  nz(i),nn(i),q(i),
     1       (as(ll),ll=llmin,llmax),ist(i)
      if (nppf.eq.6) read (5,21)  nz(i),nn(i),q(i),
     1        (as(ll),ll=llmin,llmax),ist(i)
   21 format(2i3,f11.4,f5.1,5e12.3,i2)
 777  format(2i3,f11.4,f5.1,4e12.3,i2)

      if (ist(i).ne.0.) read (5,22) (gs(ll),ll=6*i-5,6*i-6+2*ist(i))
   22 format(f10.4,f10.3,f10.4,f10.3,f10.4,f10.3,f10.4,f10.3)
c
      do 25  jj=1,ii
      if (nz(i).ne.iz(jj)) go to 25
      if (nn(i).ge.inmin(jj).and.nn(i).le.inmax(jj)) go to 29
      go to 2600
   25 continue
 2600 continue
      do 28  jj=1,10
      if (ic1(jj,i))28,28,27
   27 kmx=ic2(jj,i)

c.. jan 2000
c.. format 40 remains, but now read with format 240

c      read (5,40)  (dum(j1),j1=1,kmx)
      read (5,240)  (dum(j1),j1=1,kmx)
   28 continue

c    problem if i = 1, can end up reinitializing nppf

      i=i-1
      go to 10
c
c..   read parameters for reaction j on species i
c
   29 do 50  j=1,8
      if (ic1(j,i)) 50,50,30
   30 kmax=k+ic2(j,i)-1

c.. jan 2000 read (5,40) now read(5,240)

c      read (5,40)  (c(j1),j1=k,kmax)
      read (5,240)  (c(j1),j1=k,kmax)
   40 format (7e10.3)
  240 format (7e13.6)


      k=kmax+1
   50 continue

      go to 10
c
 60   i2=i-1
      nn(i2+1)=1
      nz(i2+1)=0
      nn(i2+2)=0
      nz(i2+2)=1
      nn(i2+3)=2
      nz(i2+3)=2
      na(i2+1)=1
      na(i2+2)=1
      na(i2+3)=4
c..  
c..   read data for electron capture on proton and positron capture
c..   on neutron.
c..  
      read (5,40) ((rrpen(j,i),i=1,7),j=1,6)
      read (5,40) ((rrnep(j,i),i=1,7),j=1,6)

      do 360 i = 1,i2
  360 nn(i) = na(i) - nz(i)
c..  
c..   naray designates linkages among various isotopes in the network
c..  
      call naray

c..   read data for fuller weak rates.  data is found at end of main
c..   raection deck and is ordered in sequence of decreasing q-value
c..   for electron capture (in electron rest masses). icode keeps the
c..   matrix location of the isotope that is beta-decaying.
c..   data is tabular with 6 values of density and 7 of temperature.
c..   five quantities are tabulated: positron decay rate, effective
c..   electron capture ft value, beta decay rate, neutrino loss rate,
c..   and anti-neutrino loss rate. see also weak1

  190 read (5,200) mz,mn,qful
  200 format (2i5,f10.4)
      if (mz.eq.99) go to 230
      do 210  i=1,i2
      if (nz(i).ne.mz) go to 210
      if (nn(i).eq.mn) go to 220
  210 continue

c..   element not in network. skip data.

  213 do 217  mm=1,35
      read (5,215) idumful
  215 format (a2)
  217 continue
      go to 190

  220 if (nrr(2,i).eq.0) go to 213
      nful=nful+1
      if (qful.gt.-1.) nfulnot=nfulnot+1
      read (5,225) ((datful(nful,j,k),j=1,6),k=1,7)
  225 format (6f10.3)
      read (5,225) ((datful(nfuldim+nful,j,k),j=1,6),k=1,7)
      read (5,225)  ((datful(2*nfuldim+nful,j,k),j=1,6),k=1,7)
      read (5,225)  ((datful(3*nfuldim+nful,j,k),j=1,6),k=1,7)
      read (5,225)  ((datful(4*nfuldim+nful,j,k),j=1,6),k=1,7)
      icode(nful)=i
      qn(nful)=qful
      go to 190

c
c..   search for the locations of key species in heavy ion reactions
c..        -c 12,o 16,ne20,na23,mg23, mg24,al27,si27,si28,p 30,p 31,s 31
c


  230 do 120  i=1,i2
      if (nz(i).eq.1.and.nn(i).eq.1) nh2=i
      if (nz(i).eq.1.and.nn(i).eq.2) nh3=i
      if (nz(i).eq.2.and.nn(i).eq.1) nhe3=i
      if (nz(i).eq.3.and.nn(i).eq.3) nli6=i
      if (nz(i).eq.3.and.nn(i).eq.4) nli7=i
      if (nz(i).eq.4.and.nn(i).eq.3) nbe7=i
      if (nz(i).eq.4.and.nn(i).eq.5) nbe9=i
      if (nz(i).eq.5.and.nn(i).eq.3) nb8=i
      if (nz(i).eq.5.and.nn(i).eq.5) nb10=i
      if (nz(i).eq.5.and.nn(i).eq.6) nb11=i
      if (nz(i).eq.6.and.nn(i).eq.5) nc11=i
      if (nz(i).eq.6.and.nn(i).eq.6) nc12=i
      if (nz(i).eq.6.and.nn(i).eq.7) nc13=i
      if (nz(i).eq.6.and.nn(i).eq.8) nc14=i
      if (nz(i).eq.7.and.nn(i).eq.6) n13=i
      if (nz(i).eq.7.and.nn(i).eq.7) n14=i
      if (nz(i).eq.7.and.nn(i).eq.8) n15=i
      if (nz(i).eq.8.and.nn(i).eq.6) no14=i
      if (nz(i).eq.8.and.nn(i).eq.7) no15=i
      if (nz(i).eq.8.and.nn(i).eq.8) no16=i
      if (nz(i).eq.8.and.nn(i).eq.9) no17=i
      if (nz(i).eq.8.and.nn(i).eq.10) no18=i
      if (nz(i).eq.9.and.nn(i).eq.8) nf17=i
      if (nz(i).eq.9.and.nn(i).eq.9) nf18=i
      if (nz(i).eq.9.and.nn(i).eq.10) nf19=i
      if (nz(i).eq.10.and.nn(i).eq.9) ne19=i
      if (nz(i).eq.10.and.nn(i).eq.10) ne20=i
      if (nz(i).eq.10.and.nn(i).eq.11) ne21=i
      if (nz(i).eq.10.and.nn(i).eq.12) ne22=i
      if (nz(i).eq.11.and.nn(i).eq.9) na20=i
      if (nz(i).eq.11.and.nn(i).eq.10) na21=i
      if (nz(i).eq.11.and.nn(i).eq.11) na22=i
      if (nz(i).eq.11.and.nn(i).eq.12) na23=i
      if (nz(i).eq.12.and.nn(i).eq.10) mg22=i
      if (nz(i).eq.12.and.nn(i).eq.11) mg23=i
      if (nz(i).eq.12.and.nn(i).eq.12) mg24=i
      if (nz(i).eq.12.and.nn(i).eq.13) mg25=i
      if (nz(i).eq.12.and.nn(i).eq.14) mg26=i
      if (nz(i).eq.13.and.nn(i).eq.10) nal23=i
      if (nz(i).eq.13.and.nn(i).eq.11) nal24=i
      if (nz(i).eq.13.and.nn(i).eq.12) nal25=i
      if (nz(i).eq.13.and.nn(i).eq.13) nal26=i
      if (nz(i).eq.13.and.nn(i).eq.14) nal27=i
      if (nz(i).eq.14.and.nn(i).eq.12) nsi26=i
      if (nz(i).eq.14.and.nn(i).eq.13) nsi27=i
      if (nz(i).eq.14.and.nn(i).eq.14) nsi28=i
      if (nz(i).eq.14.and.nn(i).eq.15) nsi29=i
      if (nz(i).eq.14.and.nn(i).eq.16) nsi30=i
      if (nz(i).eq.15.and.nn(i).eq.12) np27=i
      if (nz(i).eq.15.and.nn(i).eq.13) np28=i
      if (nz(i).eq.15.and.nn(i).eq.14) np29=i
      if (nz(i).eq.15.and.nn(i).eq.15) np30=i
      if (nz(i).eq.15.and.nn(i).eq.16) np31=i
      if (nz(i).eq.16.and.nn(i).eq.15) ns31=i
      if (ns31.ne.0) go to 130
  120 continue
c
  130 continue
c      write (6,140)
c  140 format (//,20x,'matrix locations for key elements',/)
c      write (6,150)
c  150 format (2x,'c12 ',2x,'c13 ',2x,'c14 ',2x,'n13 ',2x,'n14 ',
c     1     2x,'n15 ',2x,'o14 ',2x,'o15 ',2x,'o16 ',2x,'o17 ',2x,'o18 ',
c     2     2x,'nf17',2x,'f18 ',2x,'f19 ',2x,'ne19',2x,'ne20')
c      write (6,160) nc12,nc13,nc14,n13,n14,n15,no14,no15,no16,no17,no18,
c     1     nf17,nf18,nf19,ne19,ne20
c  160 format (15i6)
c      write (6,170)
c  170 format (/,2x,'ne21',2x,'ne22',2x,'na23',2x,'mg23',2x,'mg24',
c     1     2x,'mg25',2x,'mg26',2x,'al27',2x,'si27',2x,'si28',2x,'p30 ',
c     2     2x,'p31 ',2x,'s31 ')
c      write (6,180) ne21,ne22,na23,mg23,mg24,mg25,mg26,nal27,nsi27,
c     1     nsi28,np30,np31,ns31
c  180 format (13i6,/)
c      write(6,*) ' cf88 rates used'

      close(unit=5)

c..
c..bullet proof the heavy ion reactions
      iistop = 0
      if (nc12 .ne. 0 .or. no16 .ne. 0) then
       if (ne20  .eq. 0)  iistop = 1
       if (na23  .eq. 0)  iistop = 1
       if (mg23  .eq. 0)  iistop = 1
       if (nsi28 .eq. 0)  iistop = 1
       if (np30  .eq. 0)  iistop = 1
       if (np31  .eq. 0)  iistop = 1
       if (ns31  .eq. 0)  iistop = 1
       if (nal27 .eq. 0)  iistop = 1
       if (nsi27 .eq. 0)  iistop = 1
       if (mg24  .eq. 0)  iistop = 1
      end if
c..
      if (iistop .eq. 1) then
       write(6,*) 'missing some heavy ion matrix pointers'
       stop 'error in setting pointers'
      end if


      return
      end


c***************************************************************

      subroutine rate
      implicit real*8 (a-h,o-z)
c....
c
c may 14, 2001  - new subr for RATE derived from rtbl00_6.f
c
c....
c....
c.... generates general cross sections according to
c.... data read by setrat.
c....
c.... updated by rdh 1/22/00 to accept basis function fits
c.... to experimental data
c....

c.. may 14, 2001  must use burnf include deck - check array sizes
c      use zcom
c      include 'rtbl.dek'

      include 'burnf.dek'

      dimension recp(nburn)
      integer lk0,lk1,lk2,lk3,lk4,lk5,lk6,lk7,lk8,nrate
      character*2 cj

      real*8    fivsix,third
      parameter (fivsix = 5./6., third = 1./3.)

c.. iliadis t9 grid - t9ga(31)
      integer nt9
      parameter (nt9=31)
      dimension t9ga(nt9)
      data t9ga/0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,
     & 0.09,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
     & 1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10./

c....
c.... create various necessary variables
c....
      t9m1=1./t9
      rhom1=1./rho
      rt9=11.6048*t9m1
      t329=t9*sqrt(t9)
      t139=t9**third
      t9m13=1./t139
      t9m23=t9m13*t9m13
      t9m32=1./t329
      t92=t9*t9
      t93=t92*t9
      t539=t9*t139*t139

c.. interpolation constant (in t9) needed for experimental fits

      call hunt(t9ga,nt9,t9,jlo)

c.. protect rate by using lowest fitted value if t9.le.0.01,
c.. or highest fitted value if t9.ge.10.0

      if (jlo.eq.0) then
       jlo = 1
       t9g = t9
       ct9 = 0.
      else if (jlo.eq.nt9) then
       jlo = nt9-1
       t9g = t9ga(nt9)
       ct9 = 1.
      else
       t9g = t9
       ct9 = (t9g-t9ga(jlo))/(t9ga(jlo+1)-t9ga(jlo))
      end if
c.. necessary variables

      t92g=t9g*t9g
      t93g=t92g*t9g
      t913g=t9g**(1./3.)
      t923g=t913g*t913g
      t9m13g=1./t913g
      t9m23g=1./t923g

c....
c.... generate electron abundance
c....
      smu=ye
      rmu = 1./smu
c....
c.... zero sig array (lcm) and other rates
c....
      k=1
      call qvset (0.d0,recp,nburn)
      call blockcopy(recp,sigeq,nburn)
      call blockcopy(recp,sigeq(nburn+1),nburn)
      call blockcopy(recp,sigeq(2*nburn+1),nburn)
      call blockcopy(recp,sigeq(3*nburn+1),nburn)
      call blockcopy(recp,sigeq(4*nburn+1),nburn)
      call blockcopy(recp,sigeq(5*nburn+1),nburn)
      call blockcopy(recp,sigeq(6*nburn+1),nburn)
      call blockcopy(recp,sigeq(7*nburn+1),nburn)
      call blockcopy(recp,sigeq(8*nburn+1),nburn)
      call blockcopy(recp,sigeq(9*nburn+1),nburn)
      call blockcopy(recp,sigeq(10*nburn+1),nburn)
      call blockcopy(recp,sigeq(11*nburn+1),nburn)
      call blockcopy(recp,sigeq(12*nburn+1),nburn)
      call blockcopy(recp,sigeq(13*nburn+1),nburn)

      rpen=0.d0
      rnep=0.d0
      spen=0.d0
      snep=0.d0
      rectot=0.d0
      redtot=0.d0
      rpdtot=0.d0
      eectot=0.d0
      eedtot=0.d0
      epdtot=0.d0
      wrate=0.d0
      snuw=0.d0
c....
c....  calculate the temperature dependent partition functions w(i)
c.... g(i) is the ratio of the actual temperature dependent
c.... partition function to the ground state partition function.
c.... ground state partition functions have already been folded
c.... into the reverse rate factors used in rate.
c....
      do 12 i = 1,i2
      jd1=nppf*(i-1)+1
      jd2=jd1+1
      jd3=jd2+1
      jd4=jd3+1
      jd5=jd4+1
      jd6=jd5+1
      g(i)=0.0
      g0=1.0
      if (as(jd2).eq.0.) go to 11
      if (nppf.eq.5) then
       g(i)=exp(as(jd2)*t9m1+as(jd3)+as(jd4)*t9+as(jd5)*t92)
      else
       g(i)=exp(as(jd2)*t9m1+as(jd3)+as(jd4)*t9+as(jd5)*t92+as(jd6)*t93)
      end if
      if(ist(i).eq.0) go to 11
      do 21  jxx=6*(i-1)+1,6*(i-1)+2*ist(i)-1,2
   21 g0=g0+gs(jxx+1)*exp(-gs(jxx)*t9m1)
   11 g(i) = g0+g(i)
   12 w(i)=as(jd1)*g(i)

c....
c.... begin loop to generate sig for strong and electromagnetic
c.... reactions and for ground state beta decays.
c....
      do 190  i=1,i2
c....
c.... j = 1 for (ng)               j = 4 for (pg)
c.... j = 2 for (pn)               j = 5 for (ap)
c.... j = 3 for ground             j = 6 for (an)
c....     state beta decay         j = 7 for (ag)
c..                                j = 8 g.s. alpha decay
c....
c...
c... 2-28-00
c..
c.. add alpha decay, j now goes from 1,8
c.. nrr(jn,i) must be modified
c..
      do 180  j=1,8
c      write(6,780) i,j,nz(i),na(i),ic1(j,i),ic2(j,i)
780   format(2i5,2i4,2i5)
      jn=j
      if (jn-3) 26,25,25
   25 jn=j-1
      if (j.eq.8) jn=jn-1
   26 if (nrr(jn,i)) 27,180,27
   27 if (ic1(j,i)) 30,180,30
   30 j2=2*(j-1)+1
      j1=j2+1
c....
c.... must use at most 10 fitting constants for strong and
c.... electromagnetic rates
c....
      k0=k
      k1=k+1
      k2=k+2
      k3=k+3
      k4=k+4
      k5=k+5
      k6=k+6
      k7=k+7
      k8=k+8
      k9=k+9
c....
c.... ic1(j,i) tells which formula to use for reaction j on species i
c....
      jdex=ic1(j,i)
      go to (40,50,60,70,80,85,182,183,185,187,191,201,301,
     1 310,310,310,310,340,350,360,370,372),jdex
   40 continue
c....
c.... jdex = 1
c.... calculate rates for (ng) reactions and their inverses
c....
      sig(j2,i) = rho*c(k0)*(t9/0.348)**c(k3)*
     1     exp(c(k1)*(t9-0.348)+c(k2)*(t9-0.348)**2)
      sig(j1,i) = sig(j2,i)*rhom1*t329*c(k7)*exp(-c(k9)*t9m1)*
     1     g(i)/g(nrr(jn,i))
      go to 180
c....
c.... jdex = 2
c.... exoergic (np) and (na) reactions
c....
   50 sig(j1,i) = rho*c(k0)*exp(c(k1)*t9+c(k2)*t92
     1     +c(k3)*t93)
      sig(j2,i) = sig(j1,i)*c(k7)*g(nrr(jn,i))/g(i)*
     1     exp(-c(k9)*t9m1)
      go to 180
c....
c.... jdex = 3
c.... charged particle induced reactions
c....
   60 sig(j2,i) = rho*t9m23*exp(c(k0)-c(k1)*t9m13*
     1     (1.+c(k2)*t9+c(k3)*t92+c(k4)*t93))
      if(j.eq.2.or.j.eq.5.or.j.eq.6) go to 61
      sig(j1,i) = sig(j2,i)*rhom1*t329*g(i)/g(nrr(jn,i))*
     1     c(k6)*exp(-rt9*c(k5))
      go to 180
   61 sig(j1,i) = sig(j2,i)*c(k6)*exp(-rt9*c(k5))*g(i)/g(nrr(jn,i))
      go to 180
c....
c.... jdex = 4
c.... (ap) reaction given as (pa) reaction fit
c....
   70 sig(j1,i) = rho*t9m23*exp(c(k0)-c(k1)*t9m13*
     1     (1.+c(k2)*t9+c(k3)*t92+c(k4)*t93))
      sig(j2,i) = sig(j1,i)*c(k6)*exp(-rt9*c(k5))*g(nrr(jn,i))/g(i)
      go to 180
c....
c.... jdex = 5
c.... rates for (np) and (na) reactions with low positive q-values
c....
   80 sig(j1,i) = rho*c(k0)*exp(-c(k8)/(t9+c(k4))**third*
     1     ((1.+c(k1)*t9+c(k2)*t92+c(k3)*t93))
     2     +c(k9)/(t9+c(k5)))
      sig(j2,i) = sig(j1,i)*c(k7)*exp(-c(k9)*t9m1)*
     1     g(nrr(jn,i))/g(i)
      go to 180

c.... jdex = 6
c.... ground state plus first excited state positron decay.
c.... excited state assumed in thermal equilibrium.

   85 sig(6,i)=(c(k0)*c(k1)+c(k2)*c(k3)*exp(-rt9*c(k4)))
     1     /(c(k0)+c(k2)*exp(-rt9*c(k4)))
      go to 180
c....
c.... jdex = 7
c.... ground state beta decay
c....
  182 sig(5,i) = c(k0)
      go to 180
c....
c.... jdex = 8
c.... ground state positron decay or electron capture
c....
  183 sig(6,i) = c(k0)
      go to 180
c....
c....  jdex = 9
c.... new rate formula (9 6) for hot hydrogen burning rates (rkw)
c.... rates from oap eqn. 30 (t9.ne.t9alpha). esp. (p,g),(a,p).
c.... this is used for all * (pg) and (ap) reactions calculated by
c.... woosley (trs80) 1981 (private communication 1-2-86)
c....
  185 sig(j2,i)=rho*c(k0)*t9m23/(1.+c(k1)*t9)**fivsix*
     1     exp(c(k2)-c(k3)*((1.+c(k1)*t9)/t9)**third)
      if (j.eq.2.or.j.eq.5.or.j.eq.6) go to 91
      sig(j1,i)=c(k4)*rhom1*sig(j2,i)*t329*g(i)/g(nrr(jn,i))*
     1     exp(c(k5)*t9m1)
      go to 180
   91 sig(j1,i)=c(k4)*sig(j2,i)*g(i)/g(nrr(jn,i))*exp(c(k5)*t9m1)
      go to 180
c....
c....  jdex = 10
c.... fit for general reaction rate from wiesher et.al 1985, eqn. (9)
c.... new mg-al(p,g) resonance reactions with variable # of constants
c.... given as sum on i{ {a_i/t932*exp(-b_i/t9)}+c/t923*exp(-d/t913)}
c.... also called for wallace and woosley 81 rates with c=d=0.
c....
  187 nres = (ic2(j,i)-4)/2
      sigrs = 0.
      do 1875 nr=1,nres
      ir=k+2*(nr-1)
 1875 sigrs= sigrs + c(ir)*exp(-c(ir+1)*t9m1)/t329
      ir=k+ic2(j,i)-4
      sig(j2,i)=rho*(sigrs+c(ir)*exp(-c(ir+1)*t9m13)*t9m23)
      if(j.eq.2.or.j.eq.5.or.j.eq.6) go to 188
      sig(j1,i)=c(ir+2)*rhom1*sig(j2,i)*t329*g(i)/g(nrr(jn,i))
     1         *exp(c(ir+3)*t9m1)
      go to 180
  188 sig(j1,i)=c(ir+2)*sig(j2,i)*g(i)/g(nrr(jn,i))*exp(c(ir+3)*t9m1)
      go to 180
c....
c.... jdex = 11
c.... charged particle induced reactions t9**4 poly fit
c.... 9 fit constants read by setrat
c....
  191 sig(j2,i) = rho*t9m23*exp(c(k0)-c(k1)*t9m13*
     1     (1.+c(k2)*t9+c(k3)*t92+c(k4)*t93+c(k8)*t93*t9))
      if(j.eq.2.or.j.eq.5.or.j.eq.6) go to 192
      sig(j1,i) = sig(j2,i)*rhom1*t329*g(i)/g(nrr(jn,i))*
     1     c(k6)*exp(-rt9*c(k5))
      go to 180
  192 sig(j1,i) = sig(j2,i)*c(k6)*exp(-rt9*c(k5))*g(i)/g(nrr(jn,i))
      go to 180
c....
c.... jdex = 12
c.... (ap) reaction given as (pa) reaction fit
c.... t9**4 poly fit
c....
  201 sig(j1,i) = rho*t9m23*exp(c(k0)-c(k1)*t9m13*
     1     (1.+c(k2)*t9+c(k3)*t92+c(k4)*t93+c(k8)*t93*t9))
      sig(j2,i) = sig(j1,i)*c(k6)*exp(-rt9*c(k5))*g(nrr(jn,i))/g(i)
      go to 180
c....
c.... jdex = 13
c.... (pn) or (an) approximation to charged particle reaction from oap-422
c.... eq (30), fit constants evaluated by program ncross.for r.d.h. 1990
c....
c.... special for ncross rates
301   t9a=t9/(1.+ c(k3)*t9)
      t9a56=t9a**fivsix
      t9a13=t9a**third
      t9am13=1./t9a13
      facf=c(k0)*t9a56*t9m32
      xfacf=c(k1)-c(k2)*t9am13
      sig(j2,i)=rho*facf*exp(xfacf)
      xfacr=xfacf-c(k5)*t9m1
      sig(j1,i)=rho*facf*c(k4)*exp(xfacr)*g(nrr(jn,i))/g(i)
      go to 180
c...
c... jdex = 14, 15, 16, and 17
c... rates from f.k.thielemann reaclib deck
c... assigned for exoergic frwd (14), rev (15),
c... beta- (16), and beta+ or ec (17) rate fits
c... do not evaluate if t9.lt.0.01
c...
310   if (t9.lt.0.01) go to 180
      nrate = (ic2(j,i)-2)/7
      lk0 = k
      do 320 jk=1,nrate
       lk1 = lk0+1
       lk2 = lk1+1
       lk3 = lk2+1
       lk4 = lk3+1
       lk5 = lk4+1
       lk6 = lk5+1
       lk7 = lk6+1
       lk8 = lk7+1
       frate = 0.
       frate = exp(c(lk0)+c(lk1)*t9m1+c(lk2)*t9m13+c(lk3)*t139
     1  +c(lk4)*t9+c(lk5)*t539+c(lk6)*log(t9))
       if (ic1(j,i).eq.14) sig(j2,i) = sig(j2,i)+frate
       if (ic1(j,i).eq.15) sig(j1,i) = sig(j1,i)+frate
       if (ic1(j,i).eq.16) sig(5,i) = sig(5,i)+frate
       if (ic1(j,i).eq.17) sig(6,i) = sig(6,i)+frate
       lk0 = lk0+7
320   continue
c... no rev rate for weak
      if (ic1(j,i).ge.16.or.ic1(j,i).eq.17) go to 180
c... eval rev rate
      if (ic1(j,i).eq.14) then
       sig(j2,i) = sig(j2,i)*rho
      else
       sig(j1,i) = sig(j1,i)*rho
      end if
      if (j.eq.2.or.j.eq.5.or.j.eq.6) then
       if (ic1(j,i).eq.14) then
        sig(j1,i) = sig(j2,i)*g(i)/g(nrr(jn,i))*c(lk7)*exp(-c(lk8)*t9m1)
       else
        sig(j2,i) = sig(j1,i)*g(nrr(jn,i))/g(i)*c(lk7)*exp(-c(lk8)*t9m1)
       end if
      else
       if (ic1(j,i).eq.14) then
        sig(j1,i) = sig(j2,i)*rhom1*g(i)/g(nrr(jn,i))*c(lk7)*t329
     1     *exp(-c(lk8)*t9m1)
       else
        sig(j2,i) = sig(j1,i)*rhom1*g(nrr(jn,i))/g(i)*c(lk7)*t329
     1     *exp(-c(lk8)*t9m1)
       end if
      end if
      go to 180

c.. jdex = 18
c.. particle capture rates (n,g), (p,g), (a,g) from
c.. rath 2000, these are treated specially because
c.. they have negative Q values in the forward direction.
c.. these only exist in the forward (x,g) direction

340   if (t9.lt.0.01) go to 180
      lk0 = k
      lk1 = lk0+1
      lk2 = lk1+1
      lk3 = lk2+1
      lk4 = lk3+1
      lk5 = lk4+1
      lk6 = lk5+1
      lk7 = lk6+1
      lk8 = lk7+1

c.. the forward (x,g) rate, multiplied by rho

      sig(j2,i) = exp(c(lk0)+c(lk1)*t9m1
     1  +c(lk2)*t9m13+ c(lk3)*t139
     2  +c(lk4)*t9+c(lk5)*t539
     3  +c(lk6)*log(t9))
      sig(j2,i) = sig(j2,i)*rho

c.. the reverse (g,x) rate
c.. since we are NOT multiplying by the forward rate,
c.. the density is NOT taken out with rhom1
c.. tommy's fit with the three rev prmtrs calc the
c.. reverse rate, which must be multiplied by the
c.. ratio of the partition functions

      sig(j1,i) = exp(c(lk7)+c(lk8)*t9m1
     1  +c(lk2)*t9m13+ c(lk3)*t139
     2  +c(lk4)*t9+c(lk5)*t539
     3  +(c(lk6)+1.5)*log(t9))
      sig(j1,i) = sig(j1,i)*g(i)/g(nrr(jn,i))

      go to 180

c.. jdex = 19
c.. g.s. 1 parameter alpha decay from nwc
c.. apply by adding to (g,a) reaction channel
c.. The rate for alpha-decay is stored in nucleus
c.. card deck (z-2,a-4) in bdat. Ex: te106ad found
c.. in sn102 card deck just after sn102(a,g)

350   sig(14,i)=sig(14,i)+c(k0)
      go to 180

c.... jdex = 20
c.... ground state plus first excited state beta- decay.
c.... excited state assumed in thermal equilibrium.

360   sig(5,i)=(c(k0)*c(k1)+c(k2)*c(k3)*exp(-rt9*c(k4)))
     1     /(c(k0)+c(k2)*exp(-rt9*c(k4)))
      go to 180

c....
c.... jdex = 21
c.... basis function fits to experimental data
c.... charged particle induced reactions
c....
370   sig(j2,i) = t9m23g*exp(c(k0)-c(k1)*t9m13g*
     1     (1.+c(k2)*t9g+c(k3)*t92g+c(k4)*t93g))
c.. recall that the nt9 correction factor parameters start
c.. at (k+ic2(j,i)-nt9)
      klo = jlo+(k+ic2(j,i)-nt9-1)
      cfi = c(klo) + (c(klo+1)-c(klo))*ct9
      sig(j2,i) = sig(j2,i)*rho*exp(cfi)
      if(j.eq.2.or.j.eq.5.or.j.eq.6) go to 371
      sig(j1,i) = sig(j2,i)*rhom1*t329*g(i)/g(nrr(jn,i))*
     1     c(k6)*exp(-rt9*c(k5))
      go to 180
371   sig(j1,i) = sig(j2,i)*c(k6)*exp(-rt9*c(k5))*g(i)/g(nrr(jn,i))
      go to 180
c....
c.... jdex = 22
c.... (ap) reaction given as (pa) reaction fit
c....
372   sig(j1,i) = t9m23g*exp(c(k0)-c(k1)*t9m13g*
     1     (1.+c(k2)*t9g+c(k3)*t92g+c(k4)*t93g))
c.. recall that the nt9 correction factor parameters start at ic2(j,i)-nt9
      klo = jlo+(k+ic2(j,i)-nt9-1)
      cfi = c(klo) + (c(klo+1)-c(klo))*ct9
      sig(j1,i) = sig(j1,i)*rho*exp(cfi)
      sig(j2,i) = sig(j1,i)*c(k6)*exp(-rt9*c(k5))*g(nrr(jn,i))/g(i)
      go to 180

  180 k=k+ic2(j,i)
  190 continue
c....
c.... test to see if counting came out right
c....
      k=k-1
      if (kmax.lt.0) then
          write(3,199) kmax
199   format(' kmax .le.0! ',i7)
      end if
1999  format(2x,'kmax=',i6)

      if (k-kmax) 200,220,200
  200 write (6,210)  k,kmax
  210 format (//,17herror in rate  k=,i5,5hkmax=,i5)
      call exit(1)

c... drop hnsn and mzrk rates for version 92.1,
c....
c.... compute semi-empirical weak interaction rates according
c.... to prescriptions of mazurek and hansen. use ground state
c.... numbers if t9.lt.1.0 or rho.lt.10**5. for any rho and t9
c.... the ground state decay rate will be assumed to be a lower limit
c.... to the weak interaction rate.
c....

c... only call georges rates if t9.ge.1 and rho.ge.1.e6

  220 if (t9.lt.1.0) return
      if (rho.lt.1.d+06) return
      
c.... if mazful=1 then overwrite weak interaction rates with rates
c.... from fuller as available.

      if (mazful.ne.1) return
      tt=1.d+9*t9
      yee=ye
      rho1=rho
      call es0 (rho1,tt,yee,eta)
      call ecapnuc(eta,tt)
      call weak2 (rho1,tt,eta)
      do 280  j=1,nful
      ied=icode(j)
      iec=nrr(2,ied)
      if (iec.ne.0) sig(6,ied)=recful(j)+rpdful(j)
  280 sig(5,ied)=redful(j)
 900  continue
      return
      end

      SUBROUTINE HUNT(XX,N,X,JLO)
      implicit none

      integer jlo,n
      real*8 x,XX(N)

c.. given a array xx(1:n), and given a value x, returns a
c.. value jlo such that x is between xx(jlo) and xx(jlo+1).
c.. xx(1:n) must be monotonic, either increasing or
c.. decreasing. jlo=0 or jlo=n is returned to indicate x
c.. is out of range. jlo on input is taken as the initial
c.. guess for jlo on output.

      integer inc,jhi,jm
      LOGICAL ASCND
c..
      ASCND=XX(N).GT.XX(1)
c.. true if ascending order of table, false otherwise

      IF(JLO.LE.0.OR.JLO.GT.N)THEN

c.. input guess not useful. go immediately to bisection.

        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF

c.. set the hunting increment

      INC=1

c.. hunt up:

      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN

c.. done hunting, since off end of table.

          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN

c.. not done hunting

          JLO=JHI

c.. so double the increment and try again

          INC=INC+INC
          GO TO 1

c.. done hunting, value bracketed

        ENDIF

      ELSE

c.. hunt down

        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN

c.. done hunting, since off end of table

          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN

c.. not done hunting

          JHI=JLO

c.. so double the increment and try again

          INC=INC+INC
          GO TO 2
        ENDIF

c.. done hunting, value bracketed

      ENDIF

c.. hunt is done, so begin final bisection phase

3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END

c..
c..
c..
      subroutine neutrino (tnumu,tnue,enutot,tim,psi,taunu,
     1     timf,r09,v9nu,iret)

      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'

c..  initial coding by sew september 28, 1988
c..       last modification november 23, 1988

c..  this subroutine modifies the cross section array to take into
c..  account the effects of a neutrino flux of e, mu, and tau
c..  type neutrinos and antineutrinos. neutral current interactions
c.   are added to (g,n) and (g,p) reactions for the n and p ejecting
c..  inelastic channels. the complex interactions of 12c are treated
c..  separately bith here and in the matrix constructing routine
c..  compt. for charged current reactions three channels of de-excitatio
c..  of the daughter nucleus are considered. these are stored in a
c..  new array of cross sections called signuc(4,nburn). signuc stands
c..  for "cross section - neutrino - charged". it is a part of the
c..  zburn cliche.
c
c..        signuc(1,i) is the cross section for the nucleus i
c..                    to interact with a neutrino and eject
c..                    an electron and a neutron. it is like a (p,2n)
c..                    reaction (except that no proton is absorbed
c..                    and only one neutron is ejected) and has no
c..                    present analogue in the code.
c..        signuc(2,i) is the cross section for the nucleus i
c..                    to interact with a neutrino and eject
c..                    an electron and a proton. it makes the same produ
c..                    nucleus as a (g,n) reaction but ejects a proton
c..                    rather than a neutron.
c..        signuc(3,i) is the cross section for the nucleus i
c..                    to interact with an antineutrino and eject
c..                    a positron and a neutron. it goes to the same
c..                    product nucleus as a (g,p) reaction but makes
c..                    a neutron rather than a proton.
c..        signuc(4,i) is the cross section for the nucleus i
c..                    to interact with an antineutrino and eject
c..                    a positron and a proton. it is like a (n,2p)
c..                    reaction (except that no neutron is absorbed and
c..                    only one proton is ejected). it has no present
c..                    analogue in the code.


c..  the gamma-ray channel of the charged current interaction is
c..  a modification to the beta or positron decay rate (sig(5 or 6),
c..  remember to modifify positron rate of nucleus that flows into i
c..  not i itself)

c..  array signuc is only given a finite value if both parent and daught
c..  are carried in the network. branching ratios come from a special
c..  subroutine called branch provided by haxton which reads a data deck
c..  called branch.dat (assigned to for060 in sub branch).

c-----------------------------------------------------------------------

c..   calling arguments:

c..       tnumu = mu,tau neutrino temperature (in mev)
c..       tnue  = e-neutrino temperature (in mev)
c..       enutot = total neutrino energy (binding energy of neutron star
c..                in units of 10**53 ergs
c..       tim   = elapsed time in present run
c..       psi   = 0. if not an explosive run (i.e. psi=0 if t,rho const.
c..       taunu = e-folding time for neutrino burst
c..       timf  = time after which hyrostatic pre-processing stops
c..       r09   = inner boundary of the mass zone in question (in 10^9 c
c..       v9nu  = homologous expansion velocity (in 10^9 cm/s)

c--------------------------------------------------------------------

c..   output:

c..     modifications to sig array and generation of signuc array

c---------------------------------------------------------------------

c..    use common block information defined in main program

c..  data statement to give t-dependent branching ratios for 12c
c..  given by haxton as percentages. divide by 100 before using.

      dimension bnc12(14),bpc12(14),bac12(14),bnpc12(14),
     1     bpac12(14),bnac12(14),bhe3c12(14),bh3pc12(14),
     2     bhe3nc12(14),br(3)

c..  data from haxton ix/18/88; this is the "large pn set".

      data bnc12 /  0.008, 0.249, 1.039, 2.366, 4.042, 5.844, 7.597,
     1              9.199,10.608,11.818,12.846,13.711,14.438,15.048/

      data bpc12 /  0.452, 2.941, 7.152,12.513,18.323,23.957,29.025,
     1             33.371,36.991,39.964,42.393,44.379,46.009,47.356/

      data bac12 /  4.046, 4.142, 4.147, 4.108, 4.031, 3.931, 3.821,
     1              3.711, 3.608, 3.516, 3.436, 3.366, 3.308, 3.259/

      data bnpc12 / 0.000, 0.001, 0.022, 0.111, 0.300, 0.580, 0.918,
     1              1.283, 1.648, 1.998, 2.325, 2.624, 2.894, 3.136/

      data bpac12 / 0.000, 0.001, 0.013, 0.054, 0.135, 0.248, 0.383,
     1              0.526, 0.669, 0.805, 0.932, 1.049, 1.155, 1.251/

      data bnac12 / 0.000, 0.001, 0.013, 0.058, 0.142, 0.261, 0.400,
     1              0.545, 0.688, 0.822, 0.946, 1.058, 1.159, 1.248/

      data bhe3c12 /0.000, 0.000, 0.001, 0.004, 0.010, 0.020, 0.031,
     1              0.043, 0.055, 0.066, 0.077, 0.087, 0.096, 0.104/

      data bh3pc12 /0.000, 0.000, 0.002, 0.015, 0.047, 0.102, 0.175,
     1              0.261, 0.354, 0.449, 0.543, 0.633, 0.719, 0.799/

      data bhe3nc12 /0.000, 0.000, 0.002, 0.011, 0.034, 0.074, 0.128,
     1               0.193, 0.263, 0.336, 0.408, 0.478, 0.545, 0.608/


c..    zero signuc array

      do 20  j=1,nburn
      do 10  i=1,4
   10 signuc(i,j)=0.
   20 continue

c..    calculate flux of mu, tau, and e neutrinos in units
c..    of 10**40 cm**(-2) normalized to 10**9 cm.
c..    flux(nu) assumed equal to flux(nubar); mu's equal tau's
c..    neutrino energy assumed to be 3.15 times the temperature
c..    (for purposes of obtaining flux only). flxmutau is for each
c..    separate flavour and type of mu and tau neutrino.
c..    for example if enutot is 3 x 10**53 ergs and tmutau is 10,
c..    flxmutau is 0.0315. 10**(-5) is 10**53/(10**(40)*10**(18))

      emutau = (2./3.) * enutot

      if (tnumu.gt.0.) flxmutau=1.e-05*emutau
     1     /(3.15*tnumu*1.602e-06*4.*3.14159)
      if (tnue.gt.0.)  flxe=0.25*1.e-05*emutau
     1     /(3.15*tnue*1.602e-06*4.*3.14159)

c..   tnue=0 is an indicator to turn off the charged current reactions

      if (tnue.eq.0.) flxe=0.

c..   multiply by appropriate time dependent factor

      if (psi.eq.0) then
c.. pre processing
         fluxnu = (r09**(-2.0))*exp(-tim/taunu)/taunu
      else
c.. post processing
         tprime = exp(-(tim+timf)/taunu)
         fluxnu = tprime/taunu/((r09+v9nu*tim)**2)
      end if
      flxmutau=flxmutau*fluxnu
      flxe=flxe*fluxnu
      t10=tnumu/10.
      t10e=tnue/10.
      itnu=tnumu+0.01


c..   now calculate general cross sections and adjustments to rates
c..   for nuclei heavier than c12. assume network ordered by increasing
c..   z and a, so nuclei heavier than c12 will have network numbers
c..   greater than nc12

c..   neutral current cross sections begin with x
c..   charged current cross sections begin with y
c..   all cross sections are in units of 10**(-40) cm**2

   35 do 1000 j=nc12,i2
      ajr=na(j)
      inj=nn(j)
      izj=nz(j)

c..   nrrneut array index for loctions:
c..   1 = nu, e-,n               4 = nu  e+,n (like gp)
c..   2 = nu, e-   (like pn)     5 = nu, e+   (like np)
c..   3 = nu  e-,p (like gn)     6 = nu, e+,p
c..                              7 = nu, alpha

c..   signuc(1,2,3,4;j) goes to nrrneut(1,3,4,6;j) respectively
      if (ajr.lt.12.) go to 9999

c..   get branching ratios by calling subroutine branch
c..   first for the neutral particle channels then for the charged
c..   particle channels (provided that the nuclei are present in
c..   the network.)

      bnbeta=0.
      bpbeta=0.
      bgbeta=0.
      bnpos=0.
      bppos=0.
      bgpos=0.
      bn=0.
      bp=0.
      ba=0.
      call branch (izj,inj,0,tnumu,br,iret)
      bn=br(1)
      bp=br(2)
      ba=br(3)

      if (tnue.eq.0.) go to 45
      call branch (izj,inj,1,tnue,br,iret)
      bnbeta=br(1)
      bpbeta=br(2)

c..   assume that gamma branch is 1 minus the sum of the n,
c..   p, and a branches

      sub=br(1)+br(2)+br(3)
      if (sub.gt.0.) bgbeta=1.-sub

      call branch (izj,inj,-1,tnue,br,iret)
      bnpos=br(1)
      bppos=br(2)
      sub=br(1)+br(2)+br(3)
      if (sub.gt.0.) bgpos=1.-sub

c..  now get total cross sections by calling function signeut
c..  from haxton (x/30/88) for neutral currents and
c..  sigchrg from haxton (xi/5/88) for charged currents. these
c..  functions are good for all atomic masses to 80
c..  and tnu to 15 mev. use tnue for charged current part
c..  tnumu for neutral current part.

  45  xaj=0.
      yajb=0.
      yajp=0.
      if (tnumu.gt.0.) xaj=(signeut(inj,izj,tnumu,iret)/100.)
     1     *flxmutau*ajr
      if (tnue.gt.0.) yajb=(sigchrg(inj,izj,tnue,0,iret)/100.)
     1     *flxe*ajr
      if (tnue.gt.0.) yajp=(sigchrg(inj,izj,tnue,1,iret)/100.)
     1     *flxe*ajr

c..   only do charged currents for 12c here. do neutral
c..   currents separately after main loop.
c..   neutrino particle ejection is added to photodisintegration
c..  rates [(g,n),(g,p), and (g,a)] going into appropriate nuclei.

  800 if (j.eq.nc12) go to 900
      if (nrrneut(3,j).gt.0) sig(2,nrrneut(3,j))=
     1     sig(2,nrrneut(3,j))+xaj*bn
      if (nrrneut(4,j).gt.0) sig(8,nrrneut(4,j))=
     1     sig(8,nrrneut(4,j))+xaj*bp
      if (nrrneut(7,j).gt.0) sig(14,nrrneut(7,j))=
     1     sig(14,nrrneut(7,j))+xaj*ba
  900 signuc(1,j)=yajb*bnbeta
      signuc(2,j)=yajb*bpbeta
      signuc(3,j)=yajp*bnpos
      signuc(4,j)=yajp*bppos
      if (nrrneut(2,j).gt.0) sig(5,j)=sig(5,j)+yajb*bgbeta
      if (nrrneut(5,j).gt.0) sig(6,nrrneut(5,j))=
     1    sig(6,nrrneut(5,j))+yajp*bgpos

 1000 continue

c..   apply individual adjustments to h1, he4, and c12

c..   first the capture of electron antineutrinos on the proton
c..    use only the antineutrino, not the full flux of electron
c..    flavored neutrinos. flxe is correct. units of cross section are
c..    10**(-40) cm**2.

      sigpnubn=1.16*t10e**2
C... !!!!!!! cancel out for test 
c      rpen=rpen+sigpnubn*flxe

c..   now helium, neutral currents only.

      sighe4=(signeut(2,2,tnumu,iret)/100.)*flxmutau*4.
c..he4(g,n)
      if (nhe3.gt.0) sig(2,nhe3)=sig(2,nhe3)+0.484*sighe4
c..he4(g,p)
      if (nh3.gt.0)  sig(8,nh3)=sig(8,nh3)+0.516*sighe4

c..  charged current reactions for 12c done above. do only neutral
c..  currents here

      sigc12=(signeut(6,6,tnumu,iret)/100.)*flxmutau*12.
c..12c(g,n)
      if (nc11.gt.0) sig(2,nc11)=sig(2,nc11)+
     1     (bnc12(itnu)/100.)*sigc12
c..12c(g,p)
      if (nb11.gt.0) sig(8,nb11)=sig(8,nb11)+
     1     (bpc12(itnu)/100.)*sigc12

      ral=ral+(bac12(itnu)/100.)*sigc12
c..12c(nu,nup 3a)
      rc12np=(bnpc12(itnu)/100.)*sigc12
c..12c(nu,nup pn)
      rc12pa=(bpac12(itnu)/100.)*sigc12
c..12c(nu,nup pa)
      rc12na=(bnac12(itnu)/100.)*sigc12
c..12c(nu nup na)
      rc12he3=(bhe3c12(itnu)/100.)*sigc12
c..12c(nu,nup he3)
      rc12h3p=(bh3pc12(itnu)/100.)*sigc12
c..12c(nu nup h3 p)
      rc12he3n=(bhe3nc12(itnu)/100.)*sigc12
c..12c(nu nup he3 n)

      return

9999  write (6,999) j
 999  format (1x,' error for ion ', i2, ' in sub neutrino')
      call exit(1)

      end
c..
c..
c..
      function signeut(nn,nnz,t,iret)

      implicit real*8 (a-h,o-z)
      save

c..   calculates the quantity
c..   0.5/a * <sig(nu,nu)+sig(nubar,nubar)>
c..   in units of 10**-42 cm**2, given certain restrictions on the
c..   neutron number n, the proton number nz, and the neutrino
c..   temperature t
c..   the exact fits are use for n=nz even-even nuclei from a=4 to
c..   a=60, and for 56fe and 14n
c..   other n=nz nuclei are obtained by simple interpolations in a
c..   between the "alpha-stable" nuclei
c..   for n.ne.nz nuclei, one interpolates in a for the negative
c..   parity and positive parity responses separately.  a blocking
c..   factor is then evaluated that corrects the positive parity
c..   result for the associated changes in available initial
c..   particles or final holes.

c..   first we fix a problem associated with a=4,16,40,and 80 nuclei
c..   with n.ne.z: then interpolation in a fails
c..   because the gt strength in 4he,16o,40ca, and 80zr vanishes.
c..   it is then better to treat these nuclei as below

      signeut=0.

      if (nn+nnz.eq.4.and.nn.ne.nnz) go to 45
      if (nn+nnz.eq.16.and.nn.ne.nnz) go to 45
      if (nn+nnz.eq.40.and.nn.ne.nnz) go to 45
      if (nn+nnz.eq.80.and.nn.ne.nnz) go to 46
      n=nn
      nz=nnz
      go to 48

   46 n=nn-1
      nz=nnz-1
      go to 48
   45 n=nn+1
      nz=nnz+1
   48 continue
      a=n+nz

      if (t.lt.4.or.t.gt.15) go to 100
      if (n.eq.7.and.nz.eq.7) go to 1
      if (n.eq.30.and.nz.eq.26) go to 2
      if (a.ge.4..and.a.le.16.) go to 3
      if (a.ge.16..and.a.le.40.) go to 4
      if (a.ge.40..and.a.le.80.) go to 5

c..    restrictions on nuclei:
c..            a.ge.4
c..            a.le.80

      go to 110

c..   do 14n as a special case

    1 al=2.82
      bl=.2018
      cl=2.18
      go to 50

c..   do 56fe as a special case

    2 al=3.48
      bl=.2580
      cl=1.57
      go to 50

c..   p-shell nuclei

    3 if (n.ne.nz.or.(n/2)*2.ne.n) go to 13
      ia=a/4
      go to (70,71,72,73),ia

c..   4he

   70 al=2.11
      bl=.2697
      cl=3.61
      go to 50

c..   8be

   71 al=1.33
      bl=.2198
      cl=3.11
      go to 50

c..   12c

   72 al=1.96
      bl=.216
      cl=2.68
      go to 50

c..   16o

   73 al=1.63
      bl=.2572
      cl=3.12
      go to 50
   13 i1=a-10
      i2=i1*i1
      if (n.ne.nz) go to 15

c..   interpolate in a

      al=1.61687505245-.02025*i1+.00703124888*i2
      bl=.21220622957-.0010325*i1+.00142343855*i2
      cl=2.83624958992-.0475*i1+.01468751207*i2
      go to 50

c..   n.ne.z: interpolate in a, correct pos. par. for blocking etc.

   15 al=.70081251860+.00592499971*i1-.01257812511*i2
      bl=.2688062489-.0094625*i1-.00145156216*i2
      cl=1.78162562847-.00785002671*i1+.0775937289*i2
      am=1.30249989033-.03135000542*i1+.01012500562*i2
      bm=.29544377327+.00060749950*i1-.00024843775*i2
      cm=3.02937507629-.05420003086*i1+.00228124484*i2
      al=al*block(n,nz)
      go to 60

c..   sd shell

    4 if (n.ne.nz.or.(n/2)*2.ne.n) go to 113
      ia=a/4-4
      go to (80,81,82,83,84,85),ia

c..   20ne

   80 al=2.36
      bl=.2935
      cl=2.40
      go to 50

c..   24mg

   81 al=4.16
      bl=.2764
      cl=1.99
      go to 50

c..   28si

   82 al=4.09
      bl=.2650
      cl=1.92
      go to 50

c..   32s

   83 al=3.70
      bl=.2606
      cl=1.88
      go to 50

c..   36ar

   84 al=1.85
      bl=.2751
      cl=2.09
      go to 50

c..   40ca

   85 al=1.334
      bl=.3516
      cl=2.08
      go to 50

  113 i1=a-28
      i2=i1*i1
      if (n.ne.nz) go to 115

c..   interpolate in a for n.eq.z non-alpha-like nuclei

      al=3.92818975449-.02114285715*i1-.01869047061*i2
      bl=.26500955224+.0020589286*i1+.00027752953*i2
      cl=1.89190471172-.03437500447*i1+.00499256188*i2
      go to 50

c..   n.ne.z: interpolate in a pos. and neg. par., plus blocking

  115 al=2.32561922073-.00150892884*i1-.01760788821*i2
      bl=.31139522791+.00431696465*i1-.00073608616*i2
      cl=1.29461872578-.09465179592*i1+.01153051015*i2
      am=2.00266671181-.01641071774*i1-.00404166616*i2
      bm=.33077144623+.00213839277*i1-.00005558062*i2
      cm=2.27480983734-.02359822765*i1+.00095609826*i2
      al=al*block(n,nz)
      go to 60

c..   pf shell

    5 if (n.ne.nz.or.(n/2)*2.ne.n) go to 213
      if (n.gt.30) go to 213
      ia=a/4-10
      go to (92,93,94,95,96),ia

c..   44ti

   92 al=1.92
      bl=.2176
      cl=1.98
      go to 50

c..   48cr

   93 al=2.50
      bl=.2491
      cl=1.71
      go to 50

c..   52fe

   94 al=3.134
      bl=.2763
      cl=1.663
      go to 50

c..   56ni

   95 al=3.291
      bl=.2614
      cl=1.56
      go to 50

c..   60zn

   96 al=4.338
      bl=.2931
      cl=1.494
      go to 50

  213 i1=a-60
      i2=i1*i1
      if (n.ne.nz) go to 215

c..   n.eq.nz but not alpha-like: interpolate in a

      al=3.73080277443+.00646111881*i1-.00631654356*i2
      bl=.26686736941+.00182040210*i1+.00017159728*i2
      cl=1.50438034534-.00862230826*i1+.00107378582*i2
      go to 50

c..   n.ne.z: interpolate pos. and neg. separately, using pos.
c..   parity blocking function

  215 al=2.83276081085+.00624889275*i1-.00721747614*i2
      bl=.31738245487+.00116630865*i1-.00002931619*i2
      cl=1.2693990469-.00007593031*i1-.00016705488*i2
      am=1.23050379753-.00457822578*i1+.00025808148*i2
      bm=.35802352428+.00046168378*i1-.00000249797*i2
      cm=1.98840069771-.00905555021*i1-.00012111371*i2
      al=al*block(n,nz)
      go to 60

   50 signeut=al*(.1*t-bl)**cl
      return
   60 signeut=al*(.1*t-bl)**cl+am*(.1*t-bm)**cm
      return

  100 if (iret.eq.0) write(6,200) t
  200 format(2x,'t out of allowed region: t = ',f10.5)
      return
  110 if (iret.eq.0) write(6,201) n,nz
  201 format(2x,'n = ',i3,' nz = ',i3,' out of range')
      return

      end
c..
c..
c..
      function block(n,nz)

      implicit real*8 (a-h,o-z)
      save

      xn=n
      xnz=nz
      a2=0.5*(n+nz)
      block=(blocka(xn)+blocka(xnz))/(2.*blocka(a2))
      return
      end
c..
c..
c..
      function blocka(xn)

      implicit real*8 (a-h,o-z)
      save

c..   we've calculated the inelastic probabilities by hand for
c..   the operator sigma*tao3, using the shell order:
c..       1s,1p,1d5/2,2s1/2,1d3/2,1f7/2,2p3/2,1f5/2,2p1/2,1g9/2
c

      dimension val(101)

c..   the 1st,third, etc entries correspond to integers

      data val/0.,.25,.5,.25,0.,.4167,.8333,1.0833,1.333,1.4167,
     1 1.5,1.4167,1.3333,1.0833,.8333,.4167,0.,.45,.9,1.25,1.6,
     2 1.85,2.1,2.25,2.4,2.45,2.5,2.45,2.4,2.65,2.9,2.65,2.4,
     3 2.25,2.1,1.85,1.6,1.25,.9,.45,0.,.4643,.9286,1.3214,1.7143,
     4 2.0357,2.3571,2.6071,2.8571,3.0357,3.2143,3.3214,3.4286,
     5 3.4643,3.5,3.4643,3.4286,3.8452,4.2619,4.5119,4.7619,
     6 4.8452,4.9286,4.8452,4.7619,4.6548,4.5476,4.369,4.1905,
     7 3.9405,3.6905,3.369,3.0476,2.6548,2.2619,1.7976,1.3333,
     8 1.0833,.8333,.4167,0.,.4722,.9444,1.3611,1.7778,2.1389,
     9 2.5000,2.8056,3.1111,3.3611,3.6111,3.8056,4.,4.1389,
     9 4.2778,4.3611,4.4444,4.4722,4.5,4.4722,4.4444/

      if (xn.gt.50..or.xn.lt.0.) stop 411
      i=2*xn+1.01
      blocka=val(i)
      return
      end
c..
c..
c..
      function sigchrg(nn,nnz,t,iso,iret)

      implicit real*8 (a-h,o-z)
      save

c..   calculates the quantity
c..   1/a <sig(nu,e-)> (iso.eq.0) or 1/a <sig(nubar,e+)> (iso.ne.0)
c..   in units of 10**-42 cm**2, given certain restrictions on the
c..   neutron number n, the proton number nz, and the neutrino
c..   temperature t
c..   the exact fits are use for n=nz even-even nuclei from a=4 to
c..   a=60, and for 56fe and 14n
c..   other n=nz nuclei are obtained by simple interpolations in a
c..   between the "alpha-stable" nuclei
c..   for n.ne.nz nuclei, one interpolates in a for the negative
c..   parity and positive parity responses separately.  a blocking
c..   factor is then evaluated that corrects the positive parity
c..   result for the associated changes in available initial
c..   particles or final holes.
c
c..   first we fix a problem associated with a=4,16,40,and 80 nuclei
c..   with n.ne.z: then interpolation in a fails
c..   because the gt strength in 4he,16o,40ca, and 80zr vanishes.
c..   it is then better to treat these nuclei as below
c

      sigchrg=0.

      if (nn+nnz.eq.4.and.nn.ne.nnz) go to 45
      if (nn+nnz.eq.16.and.nn.ne.nnz) go to 45
      if (nn+nnz.eq.40.and.nn.ne.nnz) go to 45
      if (nn+nnz.eq.80.and.nn.ne.nnz) go to 46
      n=nn
      nz=nnz
      go to 48
   46 n=nn-1
      nz=nnz-1
      go to 48
   45 n=nn+1
      nz=nnz+1
   48 continue
      a=n+nz

      if (t.lt.2.5.or.t.gt.11) go to 100
      if (n.eq.7.and.nz.eq.7) go to 1
      if (n.eq.30.and.nz.eq.26) go to 2
      if (a.ge.4..and.a.le.16.) go to 3
      if (a.ge.16..and.a.le.40.) go to 4
      if (a.ge.40..and.a.le.80.) go to 5

c..    restrictions on nuclei:
c..            a.ge.4
c..            a.le.80

      go to 110

c..   do 14n as a special case

    1 if (iso.ne.0) go to 51
      al=5.52
      bl=.131
      cl=2.98
      go to 50
   51 al=2.52
      bl=.1189
      cl=2.40
      go to 50

c..   do 56fe as a special case

    2 if (iso.ne.0) go to 52
      al=11.86
      bl=.1598
      cl=2.06
      go to 50
   52 al=3.03
      bl=.1388
      cl=2.18
      go to 50

c..   p-shell nuclei

    3 if (n.ne.nz.or.(n/2)*2.ne.n) go to 13
      ia=a/4
      if (iso.ne.0) go to 53
      go to (70,71,72,73),ia

c..   4he

   70 al=5.52
      bl=.1790
      cl=4.38
      go to 50

c..   8be

   71 al=5.06
      bl=.1511
      cl=4.27
      go to 50

c..   12c

   72 al=4.26
      bl=.144
      cl=3.59
      go to 50

c..   16o

   73 al=3.94
      bl=.1657
      cl=4.1
      go to 50
   53 go to (74,75,76,77),ia

c..   4he

   74 al=2.56
      bl=.1936
      cl=3.87
      go to 50

c..   8be

   75 al=2.35
      bl=.1666
      cl=3.75
      go to 50

c..   12c

   76 al=1.94
      bl=.147
      cl=3.01
      go to 50

c..   16o

   77 al=1.64
      bl=.1632
      cl=3.52
      go to 50
   13 i1=a-10
      i2=i1*i1
      if (n.ne.nz) go to 15

c..   interpolate in a

      if (iso.ne.0) go to 14
      al=4.6681242-.13774998*i1+.00171878189*i2
      bl=.14332503-.001225*i1+.0008062489*i2
      cl=3.9081256-.03725*i1+.00921873*i2
      go to 50
   14 al=2.15624976-.079*i1-.001562491*i2
      bl=.152412498-.002845*i1+.000721875*i2
      cl=3.36312556-.04374997*i1+.00921873*i2
      go to 50

c..   n.ne.z: interpolate in a, correct pos. par. for blocking etc.

   15 if (iso.ne.0) go to 16
      al=1.1251875-.008125*i1-.016546875*i2
      bl=.20015-.006605*i1-.0016*i2
      cl=2.341-.0129*i1+.09525*i2
      am=4.4437508-.1665*i1-.004687518*i2
      bm=.18454376-.0008725*i1-.00016093776*i2
      cm=4.19625-.0325*i1-.00281249*i2
      al=al*block1(iso,n,nz)
      go to 60
   16 al=.52931255+.001625*i1-.008453126*i2
      bl=.2135-.007535*i1-.0017125*i2
      cl=1.895625-.03775*i1+.086093746*i2
      am=1.98-.0905*i1-.00125*i2
      bm=.193137538-.00223*i1-.0002718765*i2
      cm=3.67937517-.03875*i1-.00234375*i2
      al=al*block1(iso,n,nz)
      go to 60

c..   sd shell

    4 if (n.ne.nz.or.(n/2)*2.ne.n) go to 113
      ia=a/4-4
      if (iso.ne.0) go to 153
      go to (80,81,82,83,84,85),ia

c..   20ne

   80 al=5.43
      bl=.1775
      cl=3.52
      go to 50

c..   24mg

   81 al=7.66
      bl=.1868
      cl=3.09
      go to 50

c..   28si

   82 al=8.02
      bl=.1845
      cl=2.94
      go to 50

c..   32s

   83 al=7.21
      bl=.1833
      cl=2.93
      go to 50

c..   36ar

   84 al=4.01
      bl=.1709
      cl=3.29
      go to 50

c..   40ca

   85 al=3.21
      bl=.2179
      cl=3.39
      go to 50
  153 go to (86,87,88,89,90,91),ia

c..   20ne

   86 al=3.11
      bl=.1664
      cl=2.99
      go to 50

c..   24mg

   87 al=4.03
      bl=.1607
      cl=2.54
      go to 50

c..   28si

   88 al=4.02
      bl=.1369
      cl=2.40
      go to 50

c..   32s

   89 al=3.55
      bl=.1166
      cl=2.41
      go to 50

c..   36ar

   90 al=2.18
      bl=.0883
      cl=2.87
      go to 50

c..   40ca

   91 al=2.13
      bl=.196
      cl=2.65
      go to 50
  113 i1=a-28
      i2=i1*i1
      if (n.ne.nz) go to 115

c..   interpolate in a for n.eq.z non-alpha-like nuclei

      if (iso.ne.0) go to 114
      al=7.5895247-.048928577*i1-.0304613169*i2
      bl=.180481+.0012491068*i1+.00005185959*i2
      cl=2.9595236778-.0245535709*i1+.00567708677*i2
      go to 50
  114 al=3.9023811817-.00776786*i1-.01485863142*i2
      bl=.127038091-.00090982235*i1+.0003098959*i2
      cl=2.463809967-.026607159*i1+.00476190029*i2
      go to 50

c..   n.ne.z: interpolate in a pos. and neg. par., plus blocking

  115 if (iso.ne.0) go to 116
      al=3.60890508-.0008214287*i1-.026659230*i2
      bl=.216109514+.0024982139*i1-.000551488*i2
      cl=2.040476322-.056428585*i1+.01800595038*i2
      am=4.67952442-.03964287043*i1-.00863095652*i2
      bm=.20543333+.00156607106*i1-.00005877976*i2
      cm=3.47666645-.01553572062*i1+.00110119302*i2
      al=al*block2(iso,n,nz)
      go to 60

  116 al=1.555761933-.00157142908*i1-.01161458343*i2
      bl=.19789047539+.0009357144*i1-.00053244049*i2
      cl=1.44333302975-.106250025*i1+.01072917134*i2
      am=2.98857140541+.00562500115*i1-.00756696286*i2
      bm=.188547611+.00069553562*i1-.00003534201*i2
      cm=2.857619047-.02383929*i1+.001063988*i2
      al=al*block2(iso,n,nz)
      go to 60

c..   pf shell

    5 if (n.ne.nz.or.(n/2)*2.ne.n) go to 213
      if (n.gt.30) go to 213
      ia=a/4-10
      if (iso.ne.0) go to 253
      go to (92,93,94,95,96),ia

c..   44ti

   92 al=4.39
      bl=.1606
      cl=3.12
      go to 50

c..   48cr

   93 al=5.39
      bl=.1576
      cl=2.66
      go to 50

c..   52fe

   94 al=6.69
      bl=.1923
      cl=2.8
      go to 50

c..   56ni

   95 al=7.87
      bl=.1906
      cl=2.59
      go to 50

c..   60zn

   96 al=7.8
      bl=.1978
      cl=2.81
      go to 50
  253 go to (65,66,67,68,69),ia

c..   44ti

   65 al=1.90
      bl=.0234
      cl=2.79
      go to 50

c..   48cr

   66 al=2.4
      bl=.0748
      cl=2.43
      go to 50

c..   52fe

   67 al=2.78
      bl=.0718
      cl=2.26
      go to 50

c..   56ni

   68 al=2.78
      bl=.0325
      cl=2.10
      go to 50

c..   60zn

   69 al=3.25
      bl=.115
      cl=2.07
      go to 50
  213 i1=a-60
      i2=i1*i1
      if (n.ne.nz) go to 215

c..   n.eq.nz but not alpha-like: interpolate in a

      if (iso.ne.0) go to 214
      al=7.70684767+.01913354173*i1-.01113483496*i2
      bl=.186798751+.00099851389*i1+.00005964799*i2
      cl=2.63336658478-.00336403586*i1+.00152490078*i2
      go to 50
  214 al=2.9913826+.01504343*i1-.00209995429*i2
      bl=.06176359579+.00157862774*i1+.00026039057*i2
      cl=2.081406116-.016317559*i1+.00095894939*i2
      go to 50

c..   n.ne.z: interpolate pos. and neg. separately, using pos.
c..   parity blocking function

  215 if (iso.ne.0) go to 216
      al=4.95927858+.00332073285*i1-.01252509281*i2
      bl=.2114610433+.000002319509*i1-.000029167542*i2
      cl=2.33968067+.00088740385*i1-.00084552233*i2
      am=3.48863196+.004576978*i1-.00004479993*i2
      bm=.22204378+.000271869*i1-.00000033726*i2
      cm=3.3287153244-.00730226515*i1-.00009617128*i2
      al=al*block3(iso,n,nz)
      go to 60

  216 al=1.344244838-.000560901*i1-.003343085*i2
      bl=.183059737+.0009808325*i1-.00008514435*i2
      cl=1.19753778+.001738597*i1+.00041935878*i2
      am=2.50919485+.0058666477*i1-.00054663862*i2
      bm=.187255308-.00013076293*i1+.00001131436*i2
      cm=2.4000082016-.01373144*i1+.00004741201*i2
      al=al*block3(iso,n,nz)
      go to 60

   50 sigchrg=al*(.1*t-bl)**cl
      return
   60 sigchrg=al*(.1*t-bl)**cl+am*(.1*t-bm)**cm
      return
  100 if (iret.eq.0) write(6,200) t
  200 format(2x,'t out of allowed region: t = ',f10.5)
      return
  110 if (iret.eq.0) write(6,201) n,nz
  201 format(2x,'n = ',i3,' nz = ',i3,' out of range')
      return

      end
c..
c..
c..
      function block1(iso,n,nz)

      implicit real*8 (a-h,o-z)
      save

c..   blocking factor for p-shell

      a=n+nz

      if (n.gt.8.or.nz.lt.2) go to 10
      if (nz.gt.8.or.n.lt.2) go to 20

c..   n,nz in p-shell

      if (iso.ne.0) go to 1

c..   (nu,e-)

      block1=4.*(n-2.)*(8.-nz)/((a-4.)*(16.-a))
      return

c..   (nubar,e+)

    1 block1=4.*(nz-2.)*(8.-n)/((a-4.)*(16.-a))
      return
   10 if (iso.ne.0) go to 11
      block1=24*(n-nz)/((a-4.)*(16.-a))
      return
   11 block1=0
      return
   20 if (iso.ne.0) go to 12
      block1=0.
      return
   12 block1=24*(nz-n)/((a-4.)*(16.-a))
      return

      end
c..
c..
c..
      function block2(iso,n,nz)

      implicit real*8 (a-h,o-z)
      save

c..   blocking factor for sd-shell

      a=n+nz
      xn=n
      xnz=nz
      if (n.gt.20.or.nz.lt.8) go to 10
      if (nz.gt.20.or.n.lt.8) go to 20

c..   in sd shell entirely

      if (iso.ne.0) go to 1

c..   (nu,e-)

      xn1=xns(xn)
      xph1=2.-xns(xnz)
      xn2=xnd(xn)
      xph2=10.-xnd(xnz)
      xna1=xns(a/2.)
      xna2=xnd(a/2.)
      block2=(xn1*xph1/2.+xn2*xph2/10.)/(xna1*(2.-xna1)/2.+xna2*
     1 (10.-xna2)/10.)
      return

c..   (nubar,e+)

    1 xn1=xns(xnz)
      xph1=2.-xns(xn)
      xn2=xnd(xnz)
      xph2=10.-xnd(xn)
      xna1=xns(a/2.)
      xna2=xnd(a/2.)
      block2=(xn1*xph1/2.+xn2*xph2/10.)/(xna1*(2.-xna1)/2.+xna2*
     1 (10.-xna2)/10.)
      return

   10 if (iso.ne.0) go to 11
      xna1=xns(a/2.)
      xna2=xnd(a/2.)
      block2=(n-nz)/(xna1*(2-xna1)/2.+xna2*(10-xna2)/10.)
      return

   11 block2=0.
      return

   20 if (iso.ne.0.) go to 12
      block2=0.
      return

   12 xna1=xns(a/2.)
      xna2=xnd(a/2.)
      block2=(nz-n)/(xna1*(2-xna1)/2.+xna2*(10-xna2)/10.)
      return

      end
c..
c..
c..
      function xnd(x)

      implicit real*8 (a-h,o-z)
      save

      if (x.le.14.) go to 1
      if (x.le.16.) go to 2
      xnd=x-10.
      return
    1 xnd=x-8.
      return
    2 xnd=6.
      return
      end
c..
c..
c..
      function xns(x)

      implicit real*8 (a-h,o-z)
      save

      if (x.le.14.) go to 1
      if (x.le.16.) go to 2
      xns=2.
      return
    1 xns=0.
      return
    2 xns=x-14.
      return
      end
c..
c..
c..
      function block3(iso,n,nz)

      implicit real*8 (a-h,o-z)
      save

c..   blocking factor for pf-shell

      a=n+nz
      xn=n
      xnz=nz
      if (n.gt.40.or.nz.lt.20) go to 10
      if (nz.gt.40.or.n.lt.20) go to 20
      if (iso.ne.0) go to 1

c..   (nu,e-)

      xn1=xnp(xn)
      xph1=6.-xnp(xnz)
      xn2=xnf(xn)
      xph2=14.-xnf(xnz)
      xna1=xnp(a/2.)
      xna2=xnf(a/2.)
      block3=(xn1*xph1/6.+xn2*xph2/14.)/(xna1*(6.-xna1)/6.+xna2*
     1 (14.-xna2)/14.)
      return

    1 xn1=xnp(xnz)
      xph1=6.-xnp(xn)
      xn2=xnf(xnz)
      xph2=14.-xnf(xn)
      xna1=xnp(a/2.)
      xna2=xnf(a/2.)
      block3=(xn1*xph1/6.+xn2*xph2/14.)/(xna1*(6.-xna1)/6.+xna2*
     1 (14.-xna2)/14.)
      return

   10 if (iso.ne.0) go to 11
      xna1=xnp(a/2.)
      xna2=xnf(a/2.)
      block3=(n-nz)/(xna1*(6.-xna1)/6.+xna2*(14.-xna2)/14.)
      return

   11 block3=0.
      return

   20 if (iso.ne.0) go to 12
      block3=0.
      return

   12 xna1=xnp(a/2.)
      xna2=xnf(a/2.)
      block3=(nz-n)/(xna1*(6.-xna1)/6.+xna2*(14.-xna2)/14.)
      return

      end
c..
c..
c..
      function xnp(x)

      implicit real*8 (a-h,o-z)
      save

      if (x.lt.28.) go to 1
      if (x.lt.32.) go to 2
      if (x.lt.38.) go to 3
      xnp=x-34.
      return
    1 xnp=0.
      return
    2 xnp=x-28.
      return
    3 xnp=4.
      return
      end
c..
c..
c..
      function xnf(x)

      implicit real*8 (a-h,o-z)
      save

      if (x.lt.28.) go to 1
      if (x.lt.32.) go to 2
      if (x.lt.38.) go to 3
      xnf=14.
      return
    1 xnf=x-20.
      return
    2 xnf=8.
      return
    3 xnf=x-24.
      return
      end
c..
c..
c..
      subroutine branch(nz,n,icode,t,v,iret)

      implicit real*8 (a-h,o-z)
      save
      character*40 filename2

c..   returns the neutron, proton, and alpha branches in v(3),
c..   respectively
c..   nz is the target (not daughter) proton number, n the neutron,
c..   t the temperature with allowed values
c..       neutral:  between 4 and 12 mev
c..       charged:  between 3 and 7 mev
c..   simple linear interpolation on 1 mev nodes is used
c..   icode is the reaction switch:
c..       icode=-1  antineutrino charged current
c..       icode=0   neutral current
c..       icode=1   neutrino charged current
c..   if there is no table entry, zeros are returned

      dimension v(3),stor1(9,36,50),x(9),y(9),z(9)
      dimension xstor1(5,36,50),xstor2(5,36,50),xstor3(5,36,50)
      dimension ystor1(5,36,50),ystor2(5,36,50),ystor3(5,36,50)
      dimension x1(5),y1(5),z1(5)
      dimension stor2(9,36,50),stor3(9,36,50)

      data icodd/0/

      v(1)=0.
      v(2)=0.
      v(3)=0.
      if (nz.gt.36.or.n.gt.50) go to 51
      if (icodd.ne.0) go to 10

      filename2='branch.dat'
      open (unit=60, file=filename2(1:ifile), status='old',
     1     form='formatted',err=999)
      do 1 i=1,50
      do 1 j=1,36
      do 1 k=1,9
      stor2(k,j,i)=0.
      stor3(k,j,i)=0.
    1 stor1(k,j,i)=0.
      do 31 i=1,50
      do 31 j=1,36
      do 31 k=1,5
      xstor1(k,j,i)=0.
      xstor2(k,j,i)=0.
      xstor3(k,j,i)=0.
      ystor1(k,j,i)=0.
      ystor2(k,j,i)=0.
   31 ystor3(k,j,i)=0.
    2 read(60,100) i1,i2,x,y,z
      if (i1.eq.-1) go to 3
  100 format(1x,2i2,1x,9f7.4,/,6x,9f7.4,/,6x,9f7.4)
      do 4 i=1,9
      stor1(i,i1,i2)=x(i)
      stor2(i,i1,i2)=y(i)
    4 stor3(i,i1,i2)=z(i)
      go to 2
    3 read(60,101) i1,i2,x1,y1,z1
      if (i1.eq.-1) go to 5
  101 format(1x,2i2,1x,5f7.4,/,6x,5f7.4,/,6x,5f7.4)
      do 6 i=1,5
      xstor1(i,i1,i2)=x1(i)
      xstor2(i,i1,i2)=y1(i)
    6 xstor3(i,i1,i2)=z1(i)
      go to 3
    5 read(60,101) i1,i2,x1,y1,z1
      if (i1.eq.-1) go to 8
      do 7 i=1,5
      ystor1(i,i1,i2)=x1(i)
      ystor2(i,i1,i2)=y1(i)
    7 ystor3(i,i1,i2)=z1(i)
      go to 5
    8 icodd=1
      close (unit=60,status='keep',err=999)


   10 if (icode.eq.1) go to 40
      if (icode.eq.-1) go to 41

      i1=t-3
      if (t.eq.12.) i1=8
      if (i1.gt.8.or.i1.lt.1) go to 50
      i2=i1+1
      t1=i1+3
      t2=i2+3
      v1a=stor1(i1,nz,n)
      v2a=stor1(i2,nz,n)
      v1b=stor2(i1,nz,n)
      v2b=stor2(i2,nz,n)
      v1c=stor3(i1,nz,n)
      v2c=stor3(i2,nz,n)
      go to 45

   40 i1=t-2
      if (t.eq.7.) i1=4
      if (i1.gt.4.or.i1.lt.1) go to 50
      i2=i1+1
      t1=i1+2
      t2=i2+2
      v1a=xstor1(i1,nz,n)
      v2a=xstor1(i2,nz,n)
      v1b=xstor2(i1,nz,n)
      v2b=xstor2(i2,nz,n)
      v1c=xstor3(i1,nz,n)
      v2c=xstor3(i2,nz,n)
      go to 45

   41 i1=t-2
      if (t.eq.7.) i1=4
      if (i1.gt.4.or.i1.lt.1) go to 50
      i2=i1+1
      t1=i1+2
      t2=i2+2
      v1a=ystor1(i1,nz,n)
      v2a=ystor1(i2,nz,n)
      v1b=ystor2(i1,nz,n)
      v2b=ystor2(i2,nz,n)
      v1c=ystor3(i1,nz,n)
      v2c=ystor3(i2,nz,n)

   45 v(1)=(v2a-v1a)*t-v2a*t1+v1a*t2
      v(2)=(v2b-v1b)*t-v2b*t1+v1b*t2
      v(3)=(v2c-v1c)*t-v2c*t1+v1c*t2
      return

   50 if (iret.eq.0) write(6,200) t
  200 format(2x,'t = ',f5.1,' mev',2x,'out of range in branch')
      return
   51 if (iret.eq.0) write(6,201) nz,n
  201 format(2x,'nz = ',i3,1x,'n = ',i3,1x,'out of range in branch')
      return
  999 write (6,1000)
 1000 format (1x,'error opening or closing branch.dat file')
      call exit(1)
      end
c..
c..
c.. subroutine weak 1 removed from here forever
c..
      subroutine screen (t,d,z1,z2,a1,a2,scfac,iscrn,init)
      implicit real*8 (a-h,o-z)
      save

c..   this subroutine calculates screening factors for nuclear reaction
c..   rates in the weak, intermediate , and strong regimes
c..   given the temperature (t--degk), the density (d--g/cc), the
c..   atomic numbers and weights of the elements in the reaction
c..   channel with the largest coulomb barrier (z1,z2,a1,a2),
c..   and the mean plasma parameters
c..   calculated in main and passed over in common aver:
c..   (mean atomic number--zbar, mean square of the atomic number--z2bar
c..   mean atomic weight--abar, and total number of moles of
c..   nuclei per gram--ytot1).
c..   the unscreened rate is to be multiplied by the dimensionless
c..   output parameter scfac to get the corrected rate.
c..   the treatment is based on graboske, dewit, grossman, and cooper
c..   ap j. 181,457 (1973) for weak screening and on
c..   alastuey and jancovici, ap.j. 226, 1034, 1978, with plasma
c..   parameters from itoh, totsuji, setsuo, and dewitt, ap.j. 234,
c..   1079,1979,  for strong screening (rkw modification).
c..
c..space for the temp and den independent parts
      integer      nmax,iscrn,init
      parameter    (nmax=13000)
      real*8       zhat(nmax),zhat2(nmax),zfac1(nmax),zfac2(nmax),
     1             zfac3(nmax),zfac4(nmax),a3,a4,a5,a6,
     2             x13,x14,x53,x512,theta
      parameter    (x13   = 1./3., 
     1              x14   = 1./4.,
     2              x53   = 5./3.,
     3              x512  = 5./12.,
     4              theta = 1.0)
c..
c..plasma parameter communication
      common/aver/zbar,abar,z2bar,ye
c..
c..calculate average plasma parameters
      ytot1 = 1./abar
      if ( (z1*z2) .le. 0.) go to 10
      qlam0  = 1.88e+8*sqrt(d/(abar*t*t*t))
      ztilda = sqrt(z2bar+zbar*theta)
      qlam0z = qlam0*ztilda
      gamp   = 2.27493e+5*(d*zbar*ytot1)**x13/t
      taufac = 4.248719e+3/t**x13
c..
c..calculate temperature and density independent screen factors only once
      if (init .eq. 0) then
       if (iscrn .gt. nmax) then
        write(6,*) 'error: iscrn > nmax in routine screen',iscrn,nmax
        write(6,*) 'check iscrn and/or increase parameter nmax'
        stop 'error in routine screen'
       end if
       zhat(iscrn)  = (z1+z2)**x53 - z1**x53 - z2**x53
       zhat2(iscrn) = (z1+z2)**x512 - z1**x512 - z2**x512
       zfac1(iscrn) = 2.**x13 * z1*z2/(z1+z2)**x13
       zfac2(iscrn) = (z1**2*z2**2*a1*a2/(a1+a2))**x13 
       zfac3(iscrn) = (z1+z2)**x13/(2.**x13*z1*z2)
       zfac4(iscrn) = log(z1*z2/(z1+z2))
      end if
c..
c..screening factors
      gamef  = gamp * zfac1(iscrn) 
      tau12  = taufac * zfac2(iscrn)
      alph12 = 3.*gamef/tau12
c..  
c..   limit alph12 to 1.6 to prevent unphysical behavior
c..   (h dec. as rho inc.) at high rho.  this should really
c..   be replaced by a pycnonuclear reaction rate formula.

      if (alph12.le.1.6) go to 20
      alph12 = 1.6
      gamef  = 1.6*tau12/3.
      gamp   = gamef * zfac3(iscrn)

c..label weak removed from beginning of this line
20    h12w=z1*z2*qlam0z
      h12=h12w
      if(gamef.gt.0.3) go to 30
      go to 50

30    c   = 0.896434*gamp*zhat(iscrn) - 3.44740*gamp**x14*zhat2(iscrn)
     1      - 0.5551 * (log(gamp) + x53 * zfac4(iscrn)) - 2.996
      a3  = alph12 * alph12 * alph12
      a4  = a3 * alph12
      a5  = a4 * alph12
      a6  = a5 * alph12
      h12 = c - (tau12/3.) * (0.15625*a3 - 0.014*a4 - 0.0128*a5) 
     1        - gamef * (0.0055*a4 - 0.0098*a5 + 0.0048*a6)
      h12 = log(max(1.-0.0562*a3,0.77d0)) + h12
      if (gamef.gt.0.8) go to 50
      h12 = h12w * ((0.8-gamef)/0.5)+h12*((gamef-0.3)/0.5)

50    if (h12 .gt. 300.) h12=300.
      if (h12 .lt. 0.) h12=0.
      scfac = exp(h12)
      return

10    scfac=1.
      return
      end
c..
c..
c..
c..
c..
      subroutine weak2 (dd,tt,eta)
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'

c..   this subroutine calculates the net nuclear weak interaction rate
c..   (wrate, moles/g/sec) and the absolute value of the neutrino energy
c..   loss rate (snuw, erg/g/sec) associated with weak interactions
c..   involving nuclei. detailed weak reaction rate data by isotope
c..   is also calculated and stored in weakcom as specified below.

c..   the rates for electron capture, positron decay, and their
c..   associated neutrino energy losses are calculated by interpolating
c..   the data of fuller, fowler, and newman (ap.j.suppl. 1982), as
c..   revised, fit, and tabulated in may, 1982.
c..  
c..   for each isotope the base 10 logarithms of the following
c..    1) the positron decay rate (reactions per second per nucleus),
c..    2) the effective ft value for electron capture (dimensionless)
c..    3) the beta decay rate (reactions per second per nucleus)
c..    4) the mean energy of neutrinos emitted by positron decay
c..       and electron capture combined (mev per second per nucleus)
c..    5) the mean energy of anti-neutrinos emitted by beta decay.
c..  
c..   each table consists of a grid of six density*ye points
c..   (1.e+6,1.e+7,1.e+8,1.e+9,1.e+10,1.e+11) arrayed horizontally
c..   and seven temperature points (t9=1.,1.5,2.,3.,5.,10,.30.)
c..   arrayed vertically.
c..   the electron capture q-value in units of electron masses
c..   is stored in qn.
c..   the data for individual isotopes is stored in datful in the
c..   reverse order of their q-values.

c..   required input

c..   dd..  density (g/cc).
c..   tt..  temperature (degrees kelvin).
c..   ye..  electron abundance (moles/g).
c..   eta.. electron degeneracy parameter (mue/kt, dimensionless)
c..         chiu's convention (contains no electron rest mass).
c..   yk(j).abundances of the isotopes in the qnse network (moles/g)
c..   rpen..electron capture rate on a free proton (capture/s/proton
c..   rnep..positron capture rate on a free neutron (captures/s/neut
c..   spen,snep.. associated neutrino loss rates, respectively,
c..                for the above processes (erg/s/proton or neutron).
c..   datful..  weak interaction data stored in lcm


c..   detailed output

c..   recful(j).electron capture rate on element qkec(j)  (captures/nucl
c..   eec(j)..  electron capture neutrino loss rate for element qkec(j).
c..              (mec**2/nucleus/s).
c..   rpdful(j)..positron decay rate for element pkec(j).  (decays/nucle
c..   redful(j)..electron decay rate for element ekec(j).  (decays/nucle
c..   rectot..  total electron capture rate.  (moles/g/sec).
c..   eectot..  total electron capture and positron decay neutrino
c..              loss rate (erg/g/s).
c..   rpdtot..  total positron decay rate. (mole/g/sec).
c..   epdtot..  total positron decay neutrino loss rate. (erg/g/s).
c..             (now combined into eectot and epdtot set =0.)
c..   redtot..  total electron decay rate. (mole/g/sec).
c..   eedtot..  total electron decay neutrino loss rate. (erg/g/s).

c..   misc. parameters and data for fuller et al.'s rates

      parameter (pie=3.1415926536, const1o3=1./3.)
      parameter (constpi2o2=pie*pie/2.,cons2pi2o3=2.*pie*pie/3.)
      dimension  rvf(6),tvf(7)
      dimension rffdm(4),rffd0(4),rffd1(4),rffd2(4)
      dimension tffdm(5),tffd0(5),tffd1(5),tffd2(5)
      data idweak/0/
      data rvf/6.,7.,8.,9.,10.,11./
      data tvf/1.,1.5,2.,3.,5.,10.,30./

c..   calculate coefficients and load secondary data arrays
c..   if this is first pass

      if(idweak.ge.1) go to 50
      idweak=1
      ref(1)=0.d0
      rtf(1)=0.d0

c..   calculate the constant cubic interpolation parameters
c..   for fuller et al.'s rates

      do 10 k=2,4
      rffdm(k)=1./((rvf(k-1)-rvf(k))*(rvf(k-1)-rvf(k+1))*
     1  (rvf(k-1)-rvf(k+2)))
      rffd0(k)=1./((rvf(k)-rvf(k-1))*(rvf(k)-rvf(k+1))*
     1  (rvf(k)-rvf(k+2)))
      rffd1(k)=1./((rvf(k+1)-rvf(k-1))*(rvf(k+1)-rvf(k))*
     1  (rvf(k+1)-rvf(k+2)))
      rffd2(k)=1./((rvf(k+2)-rvf(k-1))*(rvf(k+2)-rvf(k))*
     1  (rvf(k+2)-rvf(k+1)))
 10   continue

      do 20 j=2,5
      tffdm(j)=1./((tvf(j-1)-tvf(j))*(tvf(j-1)-tvf(j+1))*
     1  (tvf(j-1)-tvf(j+2)))
      tffd0(j)=1./((tvf(j)-tvf(j-1))*(tvf(j)-tvf(j+1))*
     1  (tvf(j)-tvf(j+2)))
      tffd1(j)=1./((tvf(j+1)-tvf(j-1))*(tvf(j+1)-tvf(j))*
     1  (tvf(j+1)-tvf(j+2)))
      tffd2(j)=1./((tvf(j+2)-tvf(j-1))*(tvf(j+2)-tvf(j))*
     1  (tvf(j+2)-tvf(j+1)))
 20   continue

c..   calculate coefficients for fuller's phase space fitting factors

      do 30 j=1,nfulnot
      qc0(j)=(qn(j)+1.)**2
      qc1(j)=2.*(qn(j)+1.)*(qn(j)+2.)
      qc2(j)=(qn(j)+1.)**2+4.*(qn(j)+1.)+1.
      qc3(j)=2.*(qn(j)+2.)
 30   continue

      do 40 j=nfulnot+1,nful
      qc0(j)=0.d0
      qc1(j)=0.d0
      qc2(j)=qn(j)*qn(j)
      qc3(j)=-2.*qn(j)
 40   continue
      call blockcopy(qc0,qfac,iblkcpy)

c..   calculate parameters

 50   t=tt
      d=dd

      call qvset (0.0d0,recful,nfuldim)
      call qvset (0.0d0,rpdful,nfuldim)
      call qvset (0.0d0,redful,nfuldim)
      call qvset (0.0d0,eec,nfuldim)
      call qvset (0.0d0,eed,nfuldim)
      call qvset (0.0d0,ratf,nif)
      call qvset (0.0d0,ref,4*nif)
      call qvset (0.0d0,rtf,nif)
      call qvset (0.0d0,qc0,nfuldim)
      call qvset (0.0d0,qc1,nfuldim)
      call qvset (0.0d0,qc2,nfuldim)
      call qvset (0.0d0,qc3,nfuldim)
      call qvset (0.0d0,etaf,nfuldim)
      call qvset (0.0d0,ff0,nfuldim)
      call qvset (0.0d0,f1,nfuldim)
      call qvset (0.0d0,f2,nfuldim)
      call qvset (0.0d0,f3,nfuldim)
      call qvset (0.0d0,f4,nfuldim)
      call blockcopy(qfac,qc0,iblkcpy)
c..
c..      call qvset (0.d0,recful,nset45)
c..      call blockcopy(qfac,qc0,iblkcpy)

      tf=max(1.d+9,t)
      tf=min(3.d+10,tf)
      t9f=tf*1.d-9
      tme=tf/5.93014d+9
      tme2=tme*tme
      tme3=tme2*tme
      tme4=tme3*tme
      tme5=tme4*tme
      tme6=tme5*tme
      tmei=1./tme
      xlog2=1./log10(2.d0)
      r=log10(d*ye)
      rf=min(11.d0,r)
      rf=max(6.d0,rf)

c..   calculate cubic interpolation parameters
c..   for fuller et al.'s weak interaction rates

      if(t9f.lt.2.) jpf=2
      if((t9f.ge.2.).and.(t9f.lt.3.)) jpf=3
      if((t9f.ge.3.).and.(t9f.lt.5.)) jpf=4
      if(t9f.ge.5.) jpf=5
      kpf=min(max(2,int(rf)-5),4)
      rffm=rf-rvf(kpf-1)
      rff0=rf-rvf(kpf)
      rff1=rf-rvf(kpf+1)
      rff2=rf-rvf(kpf+2)
      dfacfm=rff0*rff1*rff2*rffdm(kpf)
      dfacf0=rffm*rff1*rff2*rffd0(kpf)
      dfacf1=rffm*rff0*rff2*rffd1(kpf)
      dfacf2=rffm*rff0*rff1*rffd2(kpf)
      tffm=t9f-tvf(jpf-1)
      tff0=t9f-tvf(jpf)
      tff1=t9f-tvf(jpf+1)
      tff2=t9f-tvf(jpf+2)
      tfacfm=tff0*tff1*tff2*tffdm(jpf)
      tfacf0=tffm*tff1*tff2*tffd0(jpf)
      tfacf1=tffm*tff0*tff2*tffd1(jpf)
      tfacf2=tffm*tff0*tff1*tffd2(jpf)

c..   do a cubic interpolation in fuller et al.'s data tables
c..   to find the base 10 logarithms of the:
c..    1) positron decay rates (per nucleus per second),
c..    2) effective ft values for electron capture (dimensionless),
c..    3) beta decay rate (per nucleus per second)
c..    4) combined average energy of the neutrinos emitted
c..       in these reactions (mev per reaction)
c..    5) mean antineutrino energy from beta decay (mev per sec per nuc)
c..   and store them sequentially in the ref array.
c..   ref(i,j)=dfacfm*ratf(i,kpf-1,j)+dfacf0*ratf(i,kpf,j)
c..   +dfacf1*ratf(i,kpf+1,j)+dfacf2*ratf(i,kpf+2,j)    for j=jpf-1 to j
c..   rtf(i)=tfacfm*ref(i,jpf-1)+tfacf0*ref(i,jpf)+tfacf1*ref(i,jpf+1)
c..   +tfacf2+ref(i,jpf+2)

      lrf=1
      lf=ndif*(jpf-2)+nif*(kpf-2)+1
      do 60 j=jpf-1,jpf+2

      call blockcopy (datfuleq(lf),ratf(1),nif)
      call qvmpy1(ref(lrf),ratf(1),dfacfm,nif)
      call blockcopy (datfuleq(lf+nif),ratf(1),nif)
      call qma2(ref(lrf),ratf(1),dfacf0,ref(lrf),nif)
      call blockcopy (datfuleq(lf+nif2),ratf(1),nif)
      call qma2(ref(lrf),ratf(1),dfacf1,ref(lrf),nif)
      call blockcopy (datfuleq(lf+nif3),ratf(1),nif)
      call qma2(ref(lrf),ratf(1),dfacf2,ref(lrf),nif)

      lrf=lrf+nif
      lf=lf+ndif
 60   continue

      call qvmpy1(rtf,ref(1),tfacfm,nif)
      call qma2(rtf,ref(nif+1),tfacf0,rtf,nif)
      call qma2(rtf,ref(nif2+1),tfacf1,rtf,nif)
      call qma2(rtf,ref(nif3+1),tfacf2,rtf,nif)

c..   find threshold for degenerate phase space factor calculation
c..   etaf(j)=eta+(qn(j)+1.)/tme

      etaf(nfulnot)=eta
      nfult=nfulnot+1
      call qma3(etaf(nfult),qn(nfult),tmei,eta+tmei,nful-nfulnot)
      nfuln=nfulnot
      if(etaf(nfulnot).lt.0.) go to 100

      do 70 j=nfulnot,nful
      if(etaf(j).lt.0.) go to 80
      f1(j)=exp(-etaf(j))
 70   continue
      j=nful+1

c..   calculate fermi integral approximations in the degenerate
c..   case (etaf(j).ge.0.)
c..   'f1(j)'=expf(-etaf(j))
c..   f2(j)=etaf(j)**3/3.+4.*etaf(j)+2.*f1(j)
c..   f3(j)=etaf(j)**4/4.+0.5*pie**2*etaf(j)**2+12.-6.*f1(j)
c..   f4(j)=etaf(j)**5/5.+2.*pie**2*etaf(j)**3/3.+48.*etaf(j)+24.*f1(j)

 80   nfuln=j
      nfuldeg=j-nfulnot
      if(nfuldeg.le.0) go to 100

      call qvmpy0(ff0(nfulnot),etaf(nfulnot),etaf(nfulnot),nfuldeg)
      call qma3(f4(nfulnot),ff0(nfulnot),0.2d0,cons2pi2o3,nfuldeg)
      call qma1(f4(nfulnot),f4(nfulnot),ff0(nfulnot),48.d0,nfuldeg)
      call qvmpy0(f4(nfulnot),f4(nfulnot),etaf(nfulnot),nfuldeg)
      call qma2(f4(nfulnot),f1(nfulnot),24.d0,f4(nfulnot),nfuldeg)
      call qma3(f3(nfulnot),ff0(nfulnot),0.25d0,constpi2o2,nfuldeg)
      call qma1(f3(nfulnot),f3(nfulnot),ff0(nfulnot),12.d0,nfuldeg)
      call qma2(f3(nfulnot),f1(nfulnot),-6.d0,f3(nfulnot),nfuldeg)
      call qma3(f2(nfulnot),ff0(nfulnot),const1o3,4.d0,nfuldeg)
      call qvmpy0(f2(nfulnot),f2(nfulnot),etaf(nfulnot),nfuldeg)
      call qma2(f2(nfulnot),f1(nfulnot),2.d0,f2(nfulnot),nfuldeg)

      if(nfuln.gt.nful) go to 90

c..   calculate fermi integral approximations in the non-degenerate
c..   case (etaf(j).lt.0.)
c..   f1(j)=expf(etaf(j)), f2(j)=2.*f1(j), f3(j)=6.*f1(j), f4(j)=24.*f1(

 100  nfulnodeg=nful-nfuln+1
      do 95 j=nfuln,nful
95    f1(j)=exp(etaf(j))

      call qvmpy1(f2(nfuln),f1(nfuln),2.d0,nfulnodeg)
      call qvmpy1(f3(nfuln),f1(nfuln),6.d0,nfulnodeg)
      call qvmpy1(f4(nfuln),f1(nfuln),24.d0,nfulnodeg)

c..   calculate the additional fermi integrals needed for the
c..   no threshold case

 90   if(eta.lt.-15.) ff0(nfulnot)=exp(eta)
      if(eta.gt.100.) ff0(nfulnot)=eta
      if((eta.ge.-15.).and.(eta.le.100.)) ff0(nfulnot)=
     1  log(1.+exp(eta))
      if(eta.ge.0.) f1(nfulnot)=0.5*eta*eta+2.-f1(nfulnot)

c..   calculate values of f/(alpha*gbar) for electron captures with thre
c..   f/(alpha*gbar)=ff0(j)=tme5*f4(etaf(j))+qc3(j)*tme4*f3(etaf(j))
c..    +qc2(j)*tme3*f2(etaf(j))

      nfulthres=nful-nfulnot

      call qvmpy1(f4(nfult),f4(nfult),tme5,nfulthres)
      call qmm2(f3(nfult),f3(nfult),tme4,qc3(nfult),nfulthres)
      call qmm2(f2(nfult),f2(nfult),tme3,qc2(nfult),nfulthres)
      call qaa0(ff0(nfult),f2(nfult),f3(nfult),f4(nfult),nfulthres)

c..   calculate values of f/(alpha*gbar) for electron captures
c..   without thresholds
c..   f/(alpha*gbar)=ff0(j)=tme5*f4(eta)+qc3(j)*tme4*f3(eta)+qc2(j)*tme3
c..     +qc1(j)*tme2*f1(eta)+qc0(j)*tme*ff0(eta)

      call qma3(ff0,qc0,tme*ff0(nfulnot),tme5*f4(nfulnot),nfulnot)
      call qma2(ff0,qc1,tme2*f1(nfulnot),ff0,nfulnot)
      call qma2(ff0,qc2,tme3*f2(nfulnot),ff0,nfulnot)
      call qma2(ff0,qc3,tme4*f3(nfulnot),ff0,nfulnot)

c..   complete the calculation of the electron capture
c..   rates using the interpolated ft values and the f/(alpha*gbar) phas
c..   space fitting factors.
c..   also finish the calculation of electron and positron decay rates a
c..   associated neutrino energy loss rates
c..   recful(j)=ff0(j)*10**(-rtf(nfuldim+j)-0.1591745)
c..   rpdful(j)=10.**(rtf(j))
c..   redful(j)=10.**(2*nfuldim+j))
c..   eec(j)=10.**(rtf(3*nfuldim+j))*1.9569259*(recful(j)+rpdful(j))
c..   eed(j)=10.**(rtf(4*nfuldim+j))*1.9569259*redful(j))

      call qvsub2(recful,-0.1591745d0,rtf(nfuldim+1),nful)
      call qvcopy(rtf,rpdful,nful)
      call qvcopy(rtf(2*nfuldim+1),redful,nful)
      call qvcopy(rtf(3*nfuldim+1),eec,nful)
      call qvcopy(rtf(4*nfuldim+1),eed,nful)

      do 96 j=1,nful
      recful(j)=10.**(recful(j))
      rpdful(j)=10.**(rpdful(j))
      redful(j)=10.**(redful(j))
      eec(j)=10.**(eec(j))
      eed(j)=10.**(eed(j))
 96   continue

      call qvmpy0(recful,ff0,recful,nful)
      call qam0(eec,recful,rpdful,eec,nful)
      call qvmpy1(eec,eec,1.9569259d0,nful)
      call qvmpy0(eed,redful,eed,nful)
      call qvmpy1(eed,eed,1.9569259d0,nful)

c..   include free nucleon rates

      recful(nful+1)=-rnep
      recful(nful+2)=rpen
      eec(nful+1)=snep/8.18683d-7
      eec(nful+2)=spen/8.18683d-7

c..   calculate total rates
c..   note that the neutrino energy loss from
c..   positron decay has been included in eec and eectot so
c..   that epdtot must be set to 0. to avoid double counting

      rectot=0.d0
      eectot=0.d0
      rpdtot=0.d0
      redtot=0.d0
      eedtot=0.d0
      do 110  j=1,nful
      iyed=icode(j)
      iyec=nrr(2,iyed)
      if (iyec.eq.0) go to 105
      rectot=rectot+y(iyec)*recful(j)
      rpdtot=rpdtot+y(iyec)*rpdful(j)
      eectot=eectot+y(iyec)*eec(j)
 105  redtot=redtot+y(iyed)*redful(j)
 110  eedtot=eedtot+y(iyed)*eed(j)
      rectot=rectot+y(i2+2)*rpen
      redtot=redtot+y(i2+1)*rnep
      eectot=eectot+y(i2+2)*eec(nful+2)
      eedtot=eedtot+y(i2+1)*eec(nful+1)
      eectot=eectot*6.02254d+23*8.18683d-7
      eedtot=eedtot*6.02254d+23*8.18683d-7
      spen=spen*y(i2+2)*6.02254d+23
      snep=snep*y(i2+1)*6.02254d+23

c..   calculate overall rates

      wrate=rectot+rpdtot-redtot
      snuw=eectot+eedtot

      return
      end
c..
c..
c..
      subroutine ecapnuc(etakep,temp)
      implicit real*8(a-h,o-z)
      save
c..   this subroutine calculates the rates for electron capture on
c..   protons (rpen, captures/sec/proton), positron capture on
c..   neutrons (rnep, captures/sec/neutron), and the associated
c..   neutrino loss rates (spen, ergs/sec/proton; snep, ergs/sec/neutron
c..   given the electron degeneracy parameter qeta (chemical potential
c..   excluding the electron's rest mass divided by kt)
c..   and the temperature tt in degk.
      common /ratpn/ rpen,rnep,spen,snep
c.
      data qn1,ft,twoln,cmk5,cmk6,bk,pi,qn2,c2me/-2.0716446d-06,
     &1083.9269d0,0.6931472d0,1.3635675d-49,2.2993864d-59,1.38062d-16,
     &3.1415927d0,2.0716446d-06,8.1872665d-07/
c.
c. tmean and qndeca are the mean lifetime and decay energy of the nuetro
c. xmp,xnp are masses of the p and n in grams.
      data xmp,xmn,qndeca,tmean/1.6726485d-24,1.6749543d-24,
     &1.2533036d-06,935.14d0/

      t9=temp/1.0e+9
c
c..c2me is the constant used to convert the neutrino energy
c..loss rate from mec2/s (as in the paper) to ergs/particle/sec.
c..c2me=8.1872665e-07
c
      pi2=pi*pi
      iflag=0
      qn=qn1
c..etaef uses the total fermi energy unlike keplers eta
       etaef = etakep + c2me/bk/temp
  502  iflag=iflag+1
       if(iflag.eq.1) etael=qn2/bk/temp
       if(iflag.eq.2) etael=c2me/bk/temp
c..  chanoes the electron etaf to that for positrons for p+ n -> e- + v-
       if(iflag.eq.2) etaef=-etaef
       t5=temp*temp*temp*temp*temp
       zetan = qn/bk/temp
       eta = etaef - etael

c..the significance of 680 is to protect the code from overflowing 
       if(eta.le.6.8e+02) exeta = dexp(eta)
       if(eta.gt.6.8e+02) exeta = 0.
       etael2 = etael*etael
       etael3=etael2*etael
       etael4=etael3*etael
       etael5=etael4*etael
       zetan2 = zetan*zetan
c..
       if(eta.le.6.8e+02) f0 = dlog(1.+exeta)
       if(eta.gt.6.8e+02) f0 = eta

c..if eta le. 0., the following fermi integrals apply:
       f1l = exeta
       f2l = 2.*f1l
       f3l = 6.*f1l
       f4l = 24.*f1l
       f5l = 120.*f1l
c..if eta gt. 0., the following fermi integrals apply:
       f1g=0.
       f2g=0.
       f3g=0.
       f4g=0.
       f5g=0.
       if(eta.le.0.) go to 515
       exmeta = dexp(-eta)
       eta2 = eta*eta
       eta3 = eta2*eta
       eta4 = eta3*eta
       f1g = eta2/2. + 2. - exmeta
       f2g = eta3/3. + 4.*eta + 2.*exmeta
       f3g = eta4/4. + pi2/2.*eta2 + 12. - 6.*exmeta
       f4g = eta4*eta/5. + 2.*pi2/3.*eta3 + 48.*eta
     &   + 24.*exmeta
       f5g = eta4*eta2/6. + 5./6.*pi2*eta4 + 7./6.*pi2*eta2
     &   + 240. -120.*exmeta
 515  continue
c..calculate factors which are multiplied by the fermi integrals
c..in order to compute ie.
      fac3 = 2.*zetan + 4.*etael
      fac2 = 6.*etael2 + 6.*etael*zetan + zetan2
      fac1 = 4.*etael3 + 6.*etael2*zetan + 2.*etael*zetan2
      fac0 = etael4 + 2.*zetan*etael3 + etael2*zetan2
c..compute the electron capture rates onto protons:no blocking(free-
      rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
      rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0
c..
c..
c...neutrino emission rate for electron capture:
      facv4 = 5.*etael + 3.*zetan
      facv3 = 10.*etael2 + 12.*etael*zetan + 3.*zetan2
      facv2 = 10.*etael3 + 18.*etael2*zetan +
     &   9.*etael*zetan2 + zetan2*zetan
      facv1 = 5.*etael4 + 12.*etael3*zetan +
     &   9.*etael2*zetan2 + 2.*etael*zetan2*zetan
      facv0 = etael5 + 3.*etael4*zetan
     &   + 3.*etael3*zetan2 + etael2*zetan2*zetan
      rjv1 = f5l + facv4*f4l + facv3*f3l +
     &   facv2*f2l + facv1*f1l + facv0*f0
      rjv2 = f5g + facv4*f4g + facv3*f3g +
     &   facv2*f2g + facv1*f1g + facv0*f0
c..
      if(iflag.eq.2) go to 503
      if(eta.gt.0.) go to 505
         rpen=twoln*cmk5*t5*rie1/ft
         spen=twoln*cmk6*t5*temp*rjv1/ft
         spenc=twoln*cmk6*t5*temp*rjv1/ft*c2me
         go to 504
  505 rpen=twoln*cmk5*t5*rie2/ft
      spen=twoln*cmk6*t5*temp*rjv2/ft
      spenc=twoln*cmk6*t5*temp*rjv2/ft*c2me
 504  continue
c..
      qn=qn2
c.. 
c..now compute positron capture rate onto neutron.
      go to 502
c..
 503  if(eta.gt.0.) go to 507
         rnep=twoln*cmk5*t5*rie1/ft
         snep=twoln*cmk6*t5*temp*rjv1/ft
         snepc=twoln*cmk6*t5*temp*rjv1/ft*c2me
c..       if(rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
         go to 506
  507  rnep=twoln*cmk5*t5*rie2/ft
       snep=twoln*cmk6*t5*temp*rjv2/ft
       snepc=twoln*cmk6*t5*temp*rjv2/ft*c2me
c..    if(rho.lt.1.0e+06) snep=snep+qndeca*xn(9)/xmn/tmean
506   continue

      return
      end
c..
c..
c..
      subroutine sneutx (t9,rho,rmu,q)
      implicit real*8 (a-h,o-z)
      save
c..   calculates the total neutrino energy loss rate, q,
c..   as a function of t9 (temp in 10**9 deg k) and
c..   rho (density in g/cc) in units of erg/g/s.
c..   rates taken from beaudet, petrosian, and salpeter (1967).
c..   rmu is the mean atomic weight per electron in amu.
c
      xl=t9/5.9302
      c=(rho/(rmu*1.0d+09))**.3333333
      c=c/xl
      xl2=xl*xl
      xl4=xl2*xl2
      xl8=xl4*xl4
      xl6=xl2*xl4
      xl5=xl*xl4
      xlm1=1./xl
      xlm2=xlm1*xlm1
      xlm3=xlm1*xlm2
      rm=rho/rmu
c
      xnum=(6.002d+19+2.084d+20*c+1.872d+21*c*c)*exp(-5.5924*c)
      xden=c**3+9.383d-01*xlm1-4.141d-01*xlm2+5.829d-02*xlm3
      fpa=xnum/xden
c
      xnum=(4.886d+10+7.580d+10*c+6.023d+10*c*c)*exp(-1.5654*c)
      xden=c**3+6.290d-03*xlm1+7.483d-03*xlm2+3.061d-04*xlm3
      fph=xnum/xden
c
      xnum=(2.320d-07+8.449d-08*c+1.787d-08*c*c)*exp(-.56457*c)
      xden=c**3+2.581d-02*xlm1+1.734d-02*xlm2+6.990d-04*xlm3
      fpl=xnum/xden
c
      gl=1.-13.04*xl2+133.5*xl4+1534.*xl6+918.6*xl8
c
      q=(rm**3)*fpl+rm*xl5*fph+gl*exp(-2.*xlm1)*fpa
      q=q/rho
c
      return
      end
c..
c..
c..
      subroutine es0(dd,tt,yye,eta)
      implicit real*8 (a-h,o-z)
      save
      real*8 n0,me,k

c..   fundamental constants

      parameter (n1=1,n2=2,n3=3,n4=4,n5=5,n6=6,n7=7,n8=8,n9=9,n10=10)

      parameter (   n0=6.02254d23,   sigt=6.65205d-25,  k=1.38054d-16,
     1               a=7.5648d-15,     me=9.10908d-28,  h=6.62559d-27,
     2          gee=6.670d-8,        c=2.997925d10,  pie=3.1415926536)

      parameter (  solmass=1.9892d+33,  solmassi=1./solmass  )

c the parameter list below has been moved up from the middle of the subr
      parameter(cpf=3.*h*h*h/(8.*pie*me*me*me*c*c*c),
     1 cbeta=k/(me*c*c))
      parameter(fa0=64./(9.*pie), cxne=2.488668349711d+30)
      parameter(c2mec2=2.*me*c*c)

c..   given the temperature (tt, degk), density (dd, g/cc), and
c..   electron abundance (yye, moles/g), this subroutine calculates
c..   the electron degeneracy parameter (eta, dimensionless-chiu's
c..   convention, i.e. doesn't include electron rest mass).
c..   it allows for degeneracy and relativity in its treatment
c..   of electrons and pairs.

      external funcne
c..
      real*8    twoth
      parameter (twoth = 2./3.)

c..common link to funcne
      common/compar/xncom,tcom,npflagcom
c..
c..
c..zero return arguments
       
      xne=0.d0
      pe=0.d0
      ee=0.d0
      xnebe=0.d0
      xnebt=0.d0
      pee=0.d0
      pet=0.d0
      eee=0.d0
      eet=0.d0
      xnp=0.d0
      pp=0.d0
      ep=0.d0
      xnpbt=0.d0
      xnpbep=0.d0
      ppep=0.d0
      ppt=0.d0
      epep=0.d0
      ept=0.d0

c..   restore input values

      t=tt
      d=dd
      ye=yye
      yebt=0.d0
      yebd=0.d0
      zlt=0.d0
      zld=0.d0

c..   calculate number of free matter electrons
c..   assumes that ionization is complete when pairs are present

      xnem=d*n0*ye
      xnembd=xnem*(1./d+zld+yebd/ye)
      xnembt=xnem*(zlt+yebt/ye)

c..   estimate the fermi degeneracy parameter, eta, from its
c..   asymptotic limits (c.f. cox and giuli, chapt. 24)

      beta=cbeta*t
      betai=1./beta

      if(beta.ge.1.)xnpair=2.345d+30*beta**2
      if(beta.lt.1.)xnpair=2.546d+30*beta*(1.+.75*beta)*
     1 exp(-betai)
      if(xnpair.lt.xnem) go to 10
      xeta=-0.5d0
      go to 100

 10   cdegl=log10(xnem**0.6/t)

      if(cdegl.gt.9.5) go to 20

      funcgx=(1.+fa0*beta)*sqrt(1.+fa0*beta*.5)-1.
      xeta=-log(cxne*beta*sqrt(beta*pie)*0.5*(1.+1.5*
     1 funcgx/fa0)/xnem)
      if(cdegl.lt.8.5) go to 100
      xetal=xeta

 20   pf2=(xnem*cpf)**twoth
      if(pf2.lt.1.d-6)xeta=pf2*(1.-pf2*0.25)*.5*betai
      if(pf2.ge.1.d-6)xeta=betai*(sqrt(pf2+1.)-1.)
      if(cdegl.gt.9.5) go to 100

      fxx=9.5-cdegl
      xeta=xetal*fxx+xeta*(1.-fxx)

c..   solve for the exact value of eta using the root-finding
c..   subroutine rootnr

 100  xncom=xnem
      tcom=t

c..   not carrying parameters from kepler
c..   set npflag=1 always
c..   that is, always include pairs in calculation

      npflagcom=1
      npflag=1
      jroot=49
      eta=rootnr(funcne,xeta,etaconv,jroot)

c..   calculate the negaton contribution
      call emeos(t,eta,n7,xne,pe,ee,xnebe,xnebt,pee,pet,eee,eet)
      sige=((pe+ee)/(k*t)-eta*xne)/(d*n0)

c..   bypass pairs if temperature is too low

      if((beta.le.0.02).or.(npflag.le.0)) go to 40

c..   calculate the positron contribution and add to negaton values

      etap=-eta-2./beta
      etapbeta=-1.d0
      etapbt=2./(beta*t)
      if(etap.le.-200.) etapbeta=0.d0
      if(etap.le.-200.) etapbt=0.d0
      if(etap.le.-200.) etap=-200.d0
      call emeos(t,etap,n7,xnp,pp,ep,xnpbep,xnpbt,ppep,ppt,epep,ept)

      xne=xne+xnp
      xnemc=xne-xnp
      pe=pe+pp
      ee=ee+ep+xnp*c2mec2
      xnpbe=xnpbep*etapbeta
      xnpbt=xnpbt+xnpbep*etapbt
      xnemcbe=xnebe-xnpbe
      xnebe=xnebe+xnpbe
      xnemcbt=xnebt-xnpbt
      xnebt=xnebt+xnpbt
      pee=pee+ppep*etapbeta
      pet=pet+ppt+ppep*etapbt
      eee=eee+epep*etapbeta+c2mec2*xnpbe
      eet=eet+ept+epep*etapbt+c2mec2*xnpbt
      sigp=((pp+ep)/(k*t)-etap*xnp)/(d*n0)
      go to 50

c..   no pair case

 40   xnemc=xne
      xnemcbe=xnebe
      xnemcbt=xnebt
      sigp=0.d0

c..   calculate the final density and temperature derivatives

 50   etabd=xnembd/xnemcbe
      etabt=(-xnemcbt+xnembt)/xnemcbe
      pet=pet+pee*etabt
      eet=eet+eee*etabt
      ped=pee*etabd
      eed=eee*etabd
      xnebt=xnebt+xnebe*etabt
      xnebd=xnebe*etabd

      return
      end
c..
c..
c..
      function funcne(eta,funcnebeta)
      implicit real*8 (a-h,o-z)
      save
c..   funcne returns the normalized value of the difference
c..   between the required value of xnem (xncom) and its value
c..   as a function of eta (xnemc).
c..   and funcnebeta, the derivative of funcne by eta.
c..   called by rootnr in es

      parameter (n1=1,n2=2,n3=3,n4=4,n5=5,n6=6,n7=7,n8=8,n9=9,n10=10)
      parameter(c=2.997925d+10,cbeta=1.38054d-16/(9.10908d-28*c*c))
      common/compar/xncom,tcom,npflagcom
      data xnp/0./, xnpbe/0./, xnebe/0./
      beta=cbeta*tcom
      xne=0.d0
      pd0=0.d0
      call emeos(tcom,eta,n4,xne,pd0,pd0,xnebe,pd0,pd0,pd0,pd0,pd0)
      xnemc=xne

      if((beta.le.0.02).or.(npflagcom.le.0)) go to 10
      etap=max(-eta-2./beta,-200.d0)
      call emeos(tcom,etap,n4,xnp,pd0,pd0,xnpbe,pd0,pd0,pd0,pd0,pd0)
      xnemc=xne-xnp
      if(etap.gt.-200.) xnebe=xnebe+xnpbe
 10   denom=1./(xncom+20.3*tcom*tcom*tcom)
      funcne=(xnemc-xncom)*denom
      funcnebeta=xnebe*denom
      return
      end
c..
c..
c..
      subroutine emeos(t,eta,ixne,xne,pe,ee,xnebe,xnebt,pee,pet,eee,eet)
      implicit real*8 (a-h,o-z)
      save

c..   for an arbitrary type of electron (negatron or positron), and
c..   given a value for the degeneracy parameter (=mu/kt) and temperatur
c..   this routine computes the electron number density, electron pressu
c..   and the electron kinetic energy density.
c..   the method is based on that given in  divine, ap. j. 142,1652(1965
c..   and is valid for arbitrary degeneracy and arbitrarily relativistic
c..   situations.
c..   n.b. eta is additive inverse of divine's and clayton's alpha.
c..   fe(1)=f(1/2),fe(2)=f(1),fe(3)=f(3/2), etc. are fermi integrals
c..   except fe(7)=f(-1/2).
c..   also fm(ii,7) and fx(ii,7) are tables for f(5) which
c..            is used in ecapnuc.

c..   required input
c..   t..  temperature (deg k)
c..   eta..degeneracy parameter (mu/kt) for electron type
c..           in question (follows chiu's convention)
c..   ixne..  if ixne=4, only xne is calculated.
c..           if ixne=7, all return variables are calculated.
c..           do not use other values.

c..   output quantities

c..   xne..number density (1/cc)
c..   pe.. electron pressure (erg/cc)
c..   ee.. electron energy density (erg/cc)
c..   xnebe.. dxne/deta (1/cc)
c..   xnebt.. dxne/dt (1/cc/degk)
c..   pee..dpe/deta (erg/cc)
c..   pet..dpe/dt (erg/cc/degk)
c..   eee..dee/deta (erg/cc)
c..   eet..dee/dt (erg/cc/degk)
c..
c..fermi integral declarations
      integer       nfpt,nfpt0
      parameter     (nfpt=67, nfpt0=41)
      integer       nfpt7,nfpt07
      parameter     (nfpt7=7*nfpt, nfpt07 = 7*nfpt0)
      real*8        fm(nfpt,7),fmbe(nfpt,7),etai(nfpt)
      real*8        fx(nfpt0,7),fxbe(nfpt0,7),etax(nfpt0)
      real*8        fmeq(nfpt7),fmbeeq(nfpt7),
     1              fxeq(nfpt07),fxbeeq(nfpt07)
      equivalence   (fm,fmeq)
      equivalence   (fmbe,fmbeeq)
      equivalence   (fx,fxeq)
      equivalence   (fxbe,fxbeeq)
c..
c..
      dimension gam(7),xfe(7),febe(7),fe(7),afac(7)
      real*8 n0,k,me
      parameter (n0=6.02254d+23,k=1.38054d-16,me=9.10908d-28,
     1   h=6.62559d-27,c=2.997925d+10,pie=3.1415926536,
     1   hbar=1.0544953d-27)

      parameter (pie2=pie*pie, c2=pie2/48., c4=c2*pie2*7./80.,
     1   c6=c4*pie2*1395./3528.)
      parameter (c10=2./3.,c12=6.*c2,c14=6.*c4,c16=14.*c6,
     1 c22=8.*c2,c32=30.*c2,c34=-10.*c4,c36=-10.*c6,
     2 c42=16.*c2,c50=5./7.,c52=70.*c2,c54=70.*c4,
     3 c56=14.*c6,c62=24.*c2,c64=64.*c4,c72=-2.*c2,
     4 c74=-10.*c4,c76=-42.*c6)

      parameter (cbeta=k/(me*c*c))

      parameter (cxne=2.488668349711d+30,cee=cxne*me*c*c,
     1 cpe=2.*cee/3.)
c..
c..
      data xfe/0.5966,0.7178,1.0000,1.5570,2.6527,4.8771,0./

c..
c..   tables of reduced fermi intergrals..  
c..   fm(i,j) is the fermi intergral of order j/2, evaluated
c..   at eta=etai(i), and divided by fbar=((eta+.16)**(j/2+1.))
c..   /(j/2+1.)+xfe(j) to remove most of the eta dependence.
c..
      data etai / -.16, -.14, -.12, -.10,  -.08, -.06,-.04,  -.02,  0.,
     1            .02,  .04,  .06,  .08,  .10,  .15,  .20,  .25,  .30,
     2            .35,  .40,  .45,  .50,  .60,  .70,  .80,  .90,  1.00,
     3          1.2,   1.4,  1.6,   1.8, 2.,  2.2, 2.4,  2.6,  2.8,
     4           3.,   3.2,  3.4,  3.6,  3.8,  4.,   4.2,  4.4,   4.6,
     5           4.8,  5.,   5.5,  6.,   6.5,  7.,   7.5,  8.,   8.5,
     6           9.,  9.5,  10.,  11.,  12.,  13.,  14.,  15.,  16.,
     7           17.,  18.,  19.,  20./
      data ((fm(ii,ij),ij=1,6),ii=1,8)/
     1  1.000000,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000,
     1  1.013117,  1.017019,  1.018019,  1.018598,  1.019019,  1.019332,
     1  1.023702,  1.033711,  1.036246,  1.037507,  1.038379,  1.039023,
     1  1.032646,  1.050053,  1.054641,  1.056720,  1.058084,  1.059080,
     1  1.040272,  1.066021,  1.073171,  1.076230,  1.078138,  1.079507,
     1  1.046775,  1.081595,  1.091809,  1.096028,  1.098541,  1.100310,
     1  1.052296,  1.096753,  1.110525,  1.116105,  1.119296,  1.121496,
     1  1.056942,  1.111478,  1.129293,  1.136452,  1.140403,  1.143068/
      data ((fm(ii,ij),ij=1,6),ii=9,27)/
     1  1.060804,  1.125753,  1.148086,  1.157058,  1.161862,  1.165032,
     1  1.063959,  1.139562,  1.166876,  1.177910,  1.183673,  1.187393,
     1  1.066474,  1.152892,  1.185638,  1.198998,  1.205836,  1.210155,
     1  1.068409,  1.165732,  1.204346,  1.220306,  1.228347,  1.233322,
     1  1.069819,  1.178072,  1.222971,  1.241821,  1.251204,  1.256898,
     1  1.070751,  1.189904,  1.241489,  1.263528,  1.274404,  1.280886,
     1  1.071282,  1.217227,  1.287143,  1.318522,  1.333873,  1.342678,
     1  1.069707,  1.241296,  1.331563,  1.374329,  1.395365,  1.407097,
     1  1.066516,  1.262129,  1.374367,  1.430636,  1.458755,  1.474156,
     1  1.062109,  1.279799,  1.415197,  1.487097,  1.523879,  1.543844,
     1  1.056811,  1.294431,  1.453728,  1.543342,  1.590531,  1.616121,
     1  1.050882,  1.306187,  1.489680,  1.598976,  1.658465,  1.690915,
     1  1.044534,  1.315257,  1.522814,  1.653592,  1.727391,  1.768118,
     1  1.037936,  1.321854,  1.552948,  1.706778,  1.796978,  1.847581,
     1  1.024496,  1.328517,  1.603747,  1.807238,  1.936605,  2.012468,
     1  1.011317,  1.327970,  1.641677,  1.897303,  2.073955,  2.183452,
     1  0.998845,  1.321904,  1.667116,  1.974436,  2.205241,  2.357597,
     1  0.987320,  1.311809,  1.681060,  2.036854,  2.326565,  2.531177,
     1  0.976851,  1.298932,  1.684912,  2.083663,  2.434260,  2.699802/
      data ((fm(ii,ij),ij=1,6),ii=28,46)/
     1  0.959138,  1.268614,  1.668821,  2.131264,  2.597315,  3.002927,
     1  0.945423,  1.236474,  1.631461,  2.125634,  2.681284,  3.230589,
     1  0.935134,  1.205512,  1.583147,  2.081449,  2.689582,  3.358925,
     1  0.927654,  1.177145,  1.531015,  2.013821,  2.637851,  3.384739,
     1  0.922422,  1.151882,  1.479436,  1.934962,  2.546243,  3.323774,
     1  0.918961,  1.129751,  1.430792,  1.853310,  2.433192,  3.201724,
     1  0.916878,  1.110550,  1.386189,  1.774017,  2.312554,  3.044725,
     1  0.915856,  1.093979,  1.345970,  1.699853,  2.193350,  2.873722,
     1  0.915647,  1.079714,  1.310064,  1.632041,  2.080732,  2.702955,
     1  0.916049,  1.067447,  1.278184,  1.570878,  1.977186,  2.540828,
     1  0.916909,  1.056893,  1.249958,  1.516152,  1.883535,  2.391537,
     1  0.918105,  1.047805,  1.224990,  1.467396,  1.799653,  2.256594,
     1  0.919542,  1.039966,  1.202896,  1.424043,  1.724925,  2.135948,
     1  0.921147,  1.033194,  1.183323,  1.385507,  1.658518,  2.028736,
     1  0.922863,  1.027331,  1.165956,  1.351227,  1.599544,  1.933733,
     1  0.924646,  1.022245,  1.150515,  1.320691,  1.547138,  1.849610,
     1  0.926464,  1.017824,  1.136755,  1.293438,  1.500503,  1.775072,
     1  0.928291,  1.013972,  1.124466,  1.269064,  1.458921,  1.708925,
     1  0.930108,  1.010610,  1.113464,  1.247213,  1.421762,  1.650096/
      data ((fm(ii,ij),ij=1,6),ii=47,65)/
     1  0.931901,  1.007670,  1.103591,  1.227577,  1.388472,  1.597648,
     1  0.936219,  1.001808,  1.082993,  1.186513,  1.319187,  1.489339,
     1  0.940240,  0.997555,  1.066964,  1.154422,  1.265376,  1.406086,
     1  0.943934,  0.994446,  1.054322,  1.128990,  1.222936,  1.340998,
     1  0.947305,  0.992161,  1.044231,  1.108573,  1.188987,  1.289306,
     1  0.950372,  0.990477,  1.036087,  1.091988,  1.161479,  1.247666,
     1  0.953160,  0.989239,  1.029449,  1.078374,  1.138932,  1.213693,
     1  0.955697,  0.988333,  1.023989,  1.067089,  1.120258,  1.185658,
     1  0.958009,  0.987678,  1.019463,  1.057654,  1.104645,  1.162284,
     1  0.960120,  0.987214,  1.015683,  1.049704,  1.091482,  1.142618,
     1  0.962052,  0.986897,  1.012507,  1.042958,  1.080298,  1.125934,
     1  0.965455,  0.986575,  1.007538,  1.032243,  1.062491,  1.099400,
     1  0.968346,  0.986528,  1.003918,  1.024252,  1.049147,  1.079518,
     1  0.970828,  0.986649,  1.001236,  1.018177,  1.038939,  1.064293,
     1  0.972976,  0.986871,  0.999224,  1.013483,  1.030995,  1.052418,
     1  0.974850,  0.987152,  0.997699,  1.009808,  1.024719,  1.043007,
     1  0.976499,  0.987466,  0.996535,  1.006895,  1.019698,  1.035448,
     1  0.977958,  0.987796,  0.995641,  1.004564,  1.015635,  1.029305,
     1  0.979258,  0.988131,  0.994954,  1.002682,  1.012316,  1.024260/
      data ((fm(ii,ij),ij=1,6),ii=66,67)/
     1  0.980422,  0.988464,  0.994426,  1.001152,  1.009582,  1.020078,
     1  0.981471,  0.988792,  0.994021,  0.999900,  1.007311,  1.016585/

      data ((fm(ii,ij),ij=7,7),ii=1,67)/
     + 1.0000000, 1.0199580, 1.0403098, 1.0610631, 1.0822256, 1.1038051,
     + 1.1258098, 1.1482479, 1.1711275, 1.1944574, 1.2182460, 1.2425021,
     + 1.2672348, 1.2924531, 1.3576847, 1.4261588, 1.4980308, 1.5734622,
     + 1.6526202, 1.7356781, 1.8228141, 1.9142114, 2.1105391, 2.3261794,
     + 2.5626401, 2.8213428, 3.1035136, 3.7412440, 4.4748891, 5.2861403,
     + 6.1300174, 6.9303090, 7.5886824, 8.0121513, 8.1485232, 8.0065453,
     + 7.6463371, 7.1499510, 6.5934455, 6.0326793, 5.5017276, 5.0175590,
     + 4.5859108, 4.2060742, 3.8741229, 3.5848527, 3.3328433, 2.8336372,
     + 2.4722872, 2.2047422, 2.0020915, 1.8453269, 1.7217470, 1.6226850,
     + 1.5420971, 1.4756789, 1.4203042, 1.3340125, 1.2706905, 1.2228846,
     + 1.1859375, 1.1568176, 1.1334818, 1.1145124, 1.0989003, 1.0859112,
     + 1.0750006/

c..   extended fermi integral table for eta=-4 to 0.

      data etax/-4.,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,
     1   -2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,
     2   -1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,
     3   -.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0./

      data ((fx(ii,ij),ij=1,6),ii=1,15)/
     1 0.9935882, 0.9954580, 0.9967836, 0.9977229, 0.9983882, 0.9988594,
     1 0.9929212, 0.9949845, 0.9964477, 0.9974848, 0.9982195, 0.9987399,
     1 0.9921857, 0.9944622, 0.9960771, 0.9972220, 0.9980333, 0.9986080,
     1 0.9913747, 0.9938861, 0.9956682, 0.9969320, 0.9978277, 0.9984623,
     1 0.9904808, 0.9932508, 0.9952171, 0.9966119, 0.9976007, 0.9983014,
     1 0.9894957, 0.9925503, 0.9947195, 0.9962587, 0.9973502, 0.9981238,
     1 0.9884106, 0.9917782, 0.9941707, 0.9958690, 0.9970737, 0.9979278,
     1 0.9872155, 0.9909273, 0.9935656, 0.9954392, 0.9967687, 0.9977114,
     1 0.9858999, 0.9899899, 0.9928987, 0.9949652, 0.9964321, 0.9974726,
     1 0.9844522, 0.9889576, 0.9921638, 0.9944426, 0.9960609, 0.9972091,
     1 0.9828598, 0.9878212, 0.9913541, 0.9938665, 0.9956515, 0.9969184,
     1 0.9811090, 0.9865705, 0.9904623, 0.9932317, 0.9952001, 0.9965978,
     1 0.9791851, 0.9851948, 0.9894806, 0.9925323, 0.9947025, 0.9962442,
     1 0.9770723, 0.9836823, 0.9884002, 0.9917620, 0.9941541, 0.9958543,
     1 0.9747533, 0.9820201, 0.9872118, 0.9909140, 0.9935499, 0.9954245/
      data ((fx(ii,ij),ij=1,6),ii=16,29)/
     1 0.9722100, 0.9801946, 0.9859050, 0.9899806, 0.9928846, 0.9949508,
     1 0.9694225, 0.9781908, 0.9844689, 0.9889539, 0.9921520, 0.9944289,
     1 0.9663700, 0.9759930, 0.9828916, 0.9878250, 0.9913457, 0.9938541,
     1 0.9630303, 0.9735839, 0.9811602, 0.9865842, 0.9904587, 0.9932213,
     1 0.9593799, 0.9709456, 0.9792609, 0.9852213, 0.9894833, 0.9925248,
     1 0.9553940, 0.9680586, 0.9771789, 0.9837252, 0.9884113, 0.9917585,
     1 0.9510467, 0.9649025, 0.9748985, 0.9820839, 0.9872338, 0.9909158,
     1 0.9463112, 0.9614558, 0.9724028, 0.9802845, 0.9859409, 0.9899896,
     1 0.9411597, 0.9576958, 0.9696741, 0.9783134, 0.9845225, 0.9889721,
     1 0.9355637, 0.9535990, 0.9666935, 0.9761560, 0.9829674, 0.9878550,
     1 0.9294942, 0.9491410, 0.9634414, 0.9737967, 0.9812637, 0.9866292,
     1 0.9229221, 0.9442967, 0.9598970, 0.9712192, 0.9793986, 0.9852851,
     1 0.9158184, 0.9390404, 0.9560391, 0.9684063, 0.9773587, 0.9838125,
     1 0.9081548, 0.9333462, 0.9518454, 0.9653399, 0.9751298, 0.9822001/
      data ((fx(ii,ij),ij=1,6),ii=30,41)/
     1 0.8999038, 0.9271882, 0.9472933, 0.9620012, 0.9726967, 0.9804364,
     1 0.8910396, 0.9205407, 0.9423598, 0.9583707, 0.9700437, 0.9785088,
     1 0.8815383, 0.9133786, 0.9370216, 0.9544284, 0.9671543, 0.9764043,
     1 0.8713786, 0.9056778, 0.9312555, 0.9501539, 0.9640114, 0.9741089,
     1 0.8605424, 0.8974157, 0.9250387, 0.9455263, 0.9605972, 0.9716083,
     1 0.8490153, 0.8885715, 0.9183489, 0.9405248, 0.9568935, 0.9688874,
     1 0.8367872, 0.8791267, 0.9111650, 0.9351287, 0.9528820, 0.9659305,
     1 0.8238531, 0.8690656, 0.9034670, 0.9293177, 0.9485439, 0.9627216,
     1 0.8102129, 0.8583758, 0.8952365, 0.9230718, 0.9438604, 0.9592442,
     1 0.7958728, 0.8470484, 0.8864572, 0.9163723, 0.9388129, 0.9554815,
     1 0.7808969, 0.8351   , 0.8771039, 0.9092181, 0.9333918, 0.9514126,
     1 0.7651953, 0.8225   , 0.8671884, 0.9015595, 0.9275623, 0.9470281/

      data ((fx(ii,ij),ij=7,7),ii=1,41)/
     1 0.9997143, 0.9996843, 0.9996511, 0.9996145, 0.9995741, 0.9995294,
     1 0.9994801, 0.9994256, 0.9993654, 0.9992989, 0.9992254, 0.9991444,
     1 0.9990549, 0.9989560, 0.9988469, 0.9987265, 0.9985936, 0.9984470,
     1 0.9982852, 0.9981068, 0.9979099, 0.9976929, 0.9974536, 0.9971899,
     1 0.9968993, 0.9965793, 0.9962269, 0.9958390, 0.9954121, 0.9949426,
     1 0.9944264, 0.9938592, 0.9932363, 0.9925526, 0.9918025, 0.9909804,
     1 0.9900797, 0.9890938, 0.9880156, 0.9868374, 0.9855511/


c
c..   calculate spline intepolation coefficients, gamma functions,
c..   and series constants, if this is the first pass.
      nzr0=0.d0
      if (it) 10,10,666
10    it=1
      itab=1

      gam(7)=sqrt(pie)
      gam(1)=0.5*gam(7)
      gam(2)=1.d0
      gam(3)=1.5*gam(1)
      gam(4)=2.d0
      gam(5)=2.5*gam(3)
      gam(6)=6.d0
      afac(7)=1./sqrt(2.d0)
      afac(1)=.5*afac(7)
      do 20 i1=2,6
20    afac(i1)=afac(7)*afac(i1-1)

      lf=0
      do 11 i1=1,7
      do 12 i2=1,nfpt
      fmbeeq(lf+i2)=fmeq(lf+i2)
12    continue
      call solve(etai,fmbeeq(lf+1),nfpt,nzr0,nzr0,0.d0,0.d0)
      lf=lf+nfpt
11    continue

      lf=0
      do 14 i1=1,7
      do 15 i2=1,nfpt0
      fxbeeq(lf+i2)=fxeq(lf+i2)
15    continue
      call solve(etax,fxbeeq(lf+1),nfpt0,nzr0,nzr0,0.d0,0.d0)
      lf=lf+nfpt0
14    continue
c..
c..bug fix
666   continue
c..   check range of eta
      if(eta+3.8) 31,30,30
30    if(eta+.14) 41,40,40
40    if (eta-18.) 50,50,60

c..    must use asymptotic series to calculate fermi integrals
c..   if eta lt -3.8
c..   see cox and giuli equ. 24.277.

31    eeta=exp(eta)
      do 17 l=1,7
      fe(l)=gam(l)*eeta*(1.-afac(l)*eeta)
17    continue
      go to 70

c..   calculate fermi integrals and derivatives using spline
c..   interpolated table for eta between -3.8 and -0.14

41    lf=0
      eeta=exp(eta)
      call search(eta,itab,etax,nfpt0)
      do 16 i3=1,ixne
      if(i3.eq.2) go to 80
      if(i3.eq.7)go to 16
      call seval2(eta,fe(i3),febe(i3),etax,fxeq(lf+1),
     1  fxbeeq(lf+1),nfpt0,itab)
      xk=dble(i3)*.5
      fbar=gam(i3)*eeta
      febe(i3)=(febe(i3)+fe(i3))*fbar
      fe(i3)=fe(i3)*fbar
80    lf=lf+nfpt0
16    continue
      go to 70

c..   calculate fermi integrals and derivatives using spline interpolate
c..   table for eta between -0.14 and 18.

50    lf=0
      call search(eta,itab,etai,nfpt)
      do 13 i3=1,ixne
      if(i3.eq.2) go to 90
      if(i3.eq.7) go to 13
      call seval2(eta,fe(i3),febe(i3),etai,fmeq(lf+1),fmbeeq(lf+1),
     1 nfpt,itab)
      xk=dble(i3)*.5
      fbar=((eta+.16)**(xk+1.))/(xk+1.)+xfe(i3)
      febe(i3)=febe(i3)*fbar+fe(i3)*(eta+0.16)**xk
      fe(i3)=fe(i3)*fbar
90    lf=lf+nfpt
13    continue
      go to 70

c..   eta gt 18. use asymptotic series(eq.24.38 and 24.39 of cox giuli )


60    etar=sqrt(eta)
      eta2=eta*eta
      eta2i=1./eta2
      fe1=c10*eta*etar
      fe3=0.6*eta*fe1
      fe(1)=fe1*(1.+eta2i*(c12+eta2i*(c14+eta2i*c16)))
      fe(2)=0.5*eta2+c22
      fe(3)=fe3*(1.+eta2i*(c32+eta2i*(c34+eta2i*c36)))
      fe(4)=eta*(eta2/3.+c42)
      fe(7)=2.*etar*(1.+eta2i*(c72+eta2i*(c74+eta2i*c76)))
      if(ixne.eq.4) go to 70

      fe(5)=c50*eta*fe3*(1.+eta2i*(c52+eta2i*(c54*eta2i*c56)))
      fe(6)=0.25*eta2*eta2+c62*eta2+c64

c..    must use nonrelativistic fermi functions to calculate state varia
c..       for arbitrary degeneracy and for relativistic case.
c..    must use approximation of  divine,n., ap. j. 142, 1652 (1965).
c..        see also chiu page 131.

70    beta=cbeta*t
      betar=sqrt(beta)
      beta32=beta*betar

      if(eta.le.18..and.eta.ge.-3.8) go to 100
      febe(1)=0.5*fe(7)
      febe(3)=1.5*fe(1)
      febe(4)=2.*fe(2)
      febe(5)=2.5*fe(3)
      febe(6)=3.*fe(4)

c..   cxne=sqrt(2.d0)*pie*(me*c/(hbar*pie))**3

100   u1=(fe(4)/fe(3))**2
      u1i=1./u1
      funu1r=sqrt(1.+0.5*u1*beta)
      funu1l=1.+u1*beta
      ffunu1=funu1l*funu1r
      xne = cxne * ( fe(1) + fe(3)*(ffunu1 - 1.)*u1i)*beta32
      if(ixne.eq.4) go to 110

      betai=1./beta
      beta52=beta32*beta
      u3=(fe(6)/fe(5))**2
      u3i=1./u3
      u3i2=u3i*u3i
      funu3r=sqrt(1.+0.5*u3*beta)
      funu3l=1.+u3*beta
      ffunu3=funu3l*funu3r
      gfunu3=funu3r**3
      ee  = cee  * ( fe(3) + fe(5)*(ffunu3 - 1.)*u3i)*beta52
      pe  = cpe  * ( fe(3) + fe(5)*(gfunu3 - 1.)*u3i)*beta52

c..   calculate derivatives: xnebe,xnebt,eee,eet,pee,pet

      betabt=cbeta
110   u1beta=2.*u1*(febe(4)/fe(4)-febe(3)/fe(3))
      f1bbeta=u1*(funu1r+funu1l/(4.*funu1r))
      f1beta=u1beta*beta*f1bbeta*u1i

      x=ffunu1-1.
      xnebe=cxne*beta32*(febe(1)+febe(3)*x*u1i+fe(3)*(f1beta*u1i
     1   -x*u1i*u1i*u1beta))
      if(ixne.eq.4) go to 200
      xnebt=betabt*(1.5*betai*xne+cxne*beta32*fe(3)*f1bbeta*u1i)

      u3beta=2.*u3*(febe(6)/fe(6)-febe(5)/fe(5))
      f3bbeta=u3*(funu3r+funu3l/(4.*funu3r))
      f3beta=u3beta*beta*f3bbeta*u3i
      g3beta=0.75*funu3r*beta*u3beta
      g3bbeta=0.75*funu3r*u3

      x=ffunu3-1.
      eee=cee*beta52*(febe(3)+febe(5)*x*u3i+fe(5)*(f3beta*u3i
     1   -x*u3i2*u3beta))
      eet=betabt*(2.5*betai*ee+cee*beta52*fe(5)*f3bbeta*u3i)

      x=gfunu3-1.
      pee=cpe*beta52*(febe(3)+febe(5)*x*u3i+fe(5)*(g3beta*u3i
     1   -x*u3i2*u3beta))
      pet=betabt*(2.5*betai*pe+cpe*beta52*fe(5)*g3bbeta*u3i)

200   return
      end
c..
c..
c..
      subroutine search(xval,i,x,n)
      implicit real*8 (a-h,o-z)
      save
c..   this subroutine calculates the table entry i which falls
c..   just below xval in array x(i) which is n points long.
c..   the x(i) array must increase monotonically.
c..   warning.. overflow protection removed.
      dimension x(1)
c..   if((xval.lt.x(1)).or.(xval.gt.x(n))) go to err
      if((i.le.0).or.(i.gt.n)) i=n/2+1
      if(xval-x(i)) 10,20,30
10    i=i-1
      if(xval-x(i)) 10,20,40
30    i=i+1
      if(xval-x(i))20,20,30
20    i=i-1
      if(i.eq.0)i=1
40    return
cerr  xx=0.d0
c..   y=1./xx
c..   return
      end
c..
c..
c..
      subroutine seval2(yval,x,dx,rv,hv,fyh,npts,ily)
      implicit real*8 (a-h,o-z)
      save

c..   this subroutine calculates the spline interpolated
c..   value,x, (and derivative, dx), at the point yval
c..   given the values of the function (rv--independent
c..   variable array, hv--dependent variable array) and
c..   its derivatives (fyh) at npts points. rv must
c..   increase monotonically and yval must be in the range
c..   defined by rv.  the derivative array (fyh) must be
c..   calculated beforehand using the subroutine solve.
c..   the first table entry,ily, less than yval must also
c..   be supplied.
c..
      integer   npts,ily
      real*8    yval,x,dx,rv(1),hv(1),fyh(1)
c..
      if(ily.eq.0) call search(yval,ily,rv,npts)
      ily1=ily+1
      dy=yval-rv(ily)
      hy=rv(ily1)-rv(ily)
      hyi=1./hy
      r=dy*hyi
      r2=r*r
      fh2=-2.*r2*r+3.*r2
      fh3=dy*(r2-2.*r+1.)
      fh4=dy*(r2-r)
      fh1d=6.*(r2-r)*hyi
      fh3d=3.*r2-4.*r+1.
      fh4d=3.*r2-2.*r
      x=hv(ily)*(1.-fh2)+hv(ily1)*fh2+fyh(ily)*fh3+fyh(ily1)*fh4
      dx=(hv(ily)-hv(ily1))*fh1d+fyh(ily)*fh3d+fyh(ily1)*fh4d
      return
      end
c..
c..
c..
      subroutine solve(x,y,n,if1,if2,d1,d2)
      implicit real*8 (a-h,o-z)
      save

c..   this subroutine calculates the derivatives needed for one
c..   dimensional spline interpolation.

c..   x..  array of independent variable values.
c..   y..  array of dependent variable values.
c..   n..  number of values in the x and y arrays.  note that x
c..        and y must be dimensioned at least this large in the
c..         calling program.
c..   the derivatives of y with respect to x are returned
c..   in array y.

      integer     if1,if2
      real*8      x(*),y(*),d1,d2
c..
      integer     l2,l3      
      parameter   (l2=301,l3=2*l2)
      real*8      xt(4),yt(4),al(l2),ad(l2),au(l2),b(l2),pq(l3)
c..
c..
c iflag initialized to avoid warning error
      iflag=0
      dd1=d1
      if((if1.eq.1).or.(n.lt.4)) go to 1
      dd1=funda(x,y)
    1 continue
      dd2=d2
      if((if2.eq.1).or.(n.lt.4)) go to 3
      do 2 i=1,4
      k=n+1-i
      xt(i)=x(k)
      yt(i)=y(k)
    2 continue
      dd2=funda(xt,yt)
    3 continue
      nm2=n-2
      do 4 i=1,nm2
      cl=(x(i+2)-x(i+1))/(x(i+2)-x(i))
      cm=1.-cl
      al(i)=cl
      ad(i)=2.d0
      au(i)=cm
      b(i)=3.*cl*(y(i+1)-y(i))/(x(i+1)-x(i))+
     13.*cm*(y(i+2)-y(i+1))/(x(i+2)-x(i+1))
    4 continue
      b(1)=b(1)-al(1)*dd1
      b(nm2)=b(nm2)-au(nm2)*dd2
      call trislv(al,ad,au,nm2,b,pq,y(2),iflag)
      if(iflag.eq.0) go to 6
c..   write(5,5)iflag
c5    format(///,' iflag=',i3,///)
      call exit(1)
    6 continue
      y(1)=dd1
      y(n)=dd2
      return
      end
c..
c..
c..
      function funda(x,y)
      implicit real*8 (a-h,o-z)
      save
c..   called by solve.
      real*8     x(*),y(*)
      h1=x(2)-x(1)
      h2=x(3)-x(2)
      h3=x(4)-x(3)
      h12=h1+h2
      h123=h12+h3
      h23=h2+h3
      t1=((y(2)-y(1))*h12*h123)/(h1*h2*h23)
      t2=((y(3)-y(1))*h123*h1)/(h12*h2*h3)
      t3=((y(4)-y(1))*h1*h12)/(h123*h3*h23)
      funda=t1-t2+t3
      return
      end
c..
c..
c..
      subroutine  trislv (al, ad, au, n, b, pq, x, iflag)
c..  called by solve.
      implicit real*8 (a-h,o-z)
      save
      real*8  al(n), ad(n), au(n), b(n), x(n), pq(1)
c..
c..  solution of the tridiagonal linear system
c..        a*x = b
c..  by gaussian elimination without pivoting.
c..
c..     al  is the lower diagonal of a.  ( al(1) is arbitrary.)
c..     ad  is the main diagonal of a.
c..     au  is the upper diagonal of a.  ( au(n) is arbitrary.)
c..      n  is the size if the system.
c..      b  is the right-hand side vector.
c..     pq  is a working vector of length at least  2*n .
c..  iflag  is an error flag:
c..          iflag = 0 for normal return.
c..          iflag.gt.0 for singular matrix indicated at step iflag.
c..
c..    input variables:   al, ad, au, n, b
c..    output variables:   x, iflag
c..    working variables:   pq
c..
c..  written by:         f. n. fritsch, lawrence livermore laboratory.
c..  date last changed:  13 march 1974.
c..  modified by t.a. weaver on 5/1/76 to give correct results
c..  for n=2
c..
c..
c..  .forward solution
c..
c..        notes on memory management:
c..           pq(i+1) = p(i),  i = o(1)n-1
c..           pq(n+i) = q(i),  i = 1(1)n
c..
      if(n-2) 90,50,51
 51   pq(1) = 0.d0
      pq(n) = 0.d0
      nm1 = n-1
      do 10  i = 1, nm1
         denom = al(i)*pq(i) + ad(i)
         if (denom.eq.0.)  go to 101
         pq(i+1) = -au(i)/denom
         npi = n+i
         pq(npi) = ( b(i) - al(i)*pq(npi-1) )/denom
   10 continue
      denom = al(n)*pq(n) + ad(n)
      if (denom.eq.0.)  go to 100
      pq(npi+1) = ( b(n) - al(n)*pq(npi) )/denom
c
c..  .back substitution
c
      x(n) = pq(2*n)
      do 20  k = 2, n
         i = n-k+1
         npi = n+i
         x(i) = pq(npi) + pq(i+1)*x(i+1)
   20 continue
c
c..normal return
c
      iflag = 0
      return
 50   denom=ad(1)*ad(2)-al(2)*au(1)
      i=222
      if(denom)21,101,21
 21   x(1)=(ad(2)*b(1)-au(1)*b(2))/denom
      x(2)=(ad(1)*b(2)-al(2)*b(1))/denom
      iflag=0
      return

cbad  wot 59,f1,n
cf1    format(1x,'n.lt.2 in trislv, n=',i10)
 90   return
c
c..error return
c
  100 continue
      i = n
  101 continue
      iflag = i
      return
      end
c..
c..
c..
      function rootnr(func,gess,xconv,j)
      implicit real*8 (a-h,o-z)
      save
c..   this routine finds the root of the function, func, to a
c..   relative accuracy of xconv, given an initial guess, gess,
c..   for the root.  if convergence is not reached in j cycles
c..   j is changed to minus j (initial) as a return error flag.
c..   otherwise the returned j is the number of iterations
c..   performed.  func must be a function of two arguments, the
c..   first being the independent variable, and the second being
c..   the derivative of func with respect to the indep variable.
c..   a newton-raphson scheme is used to find the root.

      external func
      double precision func
      dimension x(50),y(50),dy(50)
      data dydx/0./
      jroot=j
      x(1)=gess
      do 10 j=1,jroot
      y(j)=func(x(j),dydx)
      dy(j)=dydx
      dx=-y(j)/(dydx+1.d-199)
      x(j+1)=x(j)+dx
      if(abs(dx)-abs(x(1)*xconv)) 20,20,10
10    continue
      rootnr=x(jroot+1)
      j=-jroot
      return

20    rootnr=x(j+1)
      return
      end
c..
c..
c..
      subroutine spcrat
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'
c..
      real*8    fivsix,third
      parameter (fivsix = 5./6., third = 1./3.)
c..  
c..  
c..   spcrat computes special rates not computed by rate
c..   screening factors are computed in subroutine screen
c..   and applied here
c..  
c..   updated aug 20/1988 to add new reactions for low z elements li,be,
c..   updated june/88 cf88 revisions for 3-alpha; c12+c12; c12+o16; o16+
c..  
c..   compute necessary variables
c..  
      t139=t9**.333333
      t239=t139*t139
      t9m1=1./t9
      t9m23=1./t239
      t129=sqrt(t9)
      t569=t139*t129
      rt9=11.60485*t9m1
      t932=t9*t129
      t913=t139
      t923=t239
      t943=t9*t913
      t953=t9*t923
c
c..    triple alpha reaction rate
c..    rate is for formation of c12 compound nucleus
c..    (has been divided by 6)
c
c..   rate revised according to caughlan and fowler (nomoto ref.) 1988
c
c..   q = -0.092
      r2abe=(7.40d+05/t932)*exp(-1.0663*t9m1)
     1     +4.164d+09/t923*exp(-13.49/t913-(t9/0.098)**2)
     2     *(1.+0.031*t913+8.009*t923+1.732*t9+49.883
     3     *t943+27.426*t953)
c..   q = 7.367
      rbeac=(130./t932)*exp(-3.3364/t9)
     1     +2.510d+07/t923*exp(-23.57/t913-(t9/0.235)**2)
     2     *(1.+0.018*t913+5.249*t923+0.650*t9+19.176
     3     *t943+6.034*t953)
c..   q = 7.275
      if (t9.gt.0.08) then
         ra3=2.90d-16*(r2abe*rbeac)
     1   +0.1*1.35d-07/t932*exp(-24.811/t9)
      else
         ra3=2.90d-16*(r2abe*rbeac)
     1   *(0.01 + 0.2*(1. + 4.*exp(-(0.025/t9)**3.263))
     2   /(1. + 4.*exp(-(t9/0.025)**9.227)))
     3   +0.1*1.35d-07/t932*exp(-24.811/t9)
      end if
      ra3=ra3*(rho*rho)*sc3a/6.
      rev=2.00d+20*exp(-84.424*t9m1)
c..    ral=rev*(t9**3)*6.*ra3/(rho*rho)
c..  do not screen c12photo
      ral=rev*(t9**3)*6.*ra3/((rho*rho)*sc3a)
c
c..       c12 + c12  reaction
c..   rate revised according to caughlan and fowler 1988
c..   see cf88 references 47-49
c
      t9a=t9/(1.+0.0396*t9)
      t9a13=t9a**third
      t9a56=t9a**fivsix
      r24 = 4.27e+26*t9a56/t932*exp(-84.165/t9a13-2.12e-03*t9**3)
      r24 = 0.5*rho*r24*sc1212
c
c..   r24*pn(1)**2 is formation rate of mg24 cmpd. nucleus
c..    branching ratios from fcz 1975
c..      neutron branching from dayras switkowski and woosley
c..      1976
c
            if(t9.ge.1.5) go to 21
c..  
c..   repair typographical error in dayras et al neutron branching
c..   ratio. (t9**2 not t9**3)
c..  
      b24n = 0.859*exp(-((0.766/t9**3)*(1.+0.0789
     1    *t9+7.74*t9**2)))
      go to 22
   21 b24n = 0.055*(1.-exp(-(0.789*t9-0.976)))
   22 continue
      if(t9.gt.3.) go to 23
      b24p = (1.-b24n)/2.
      b24a = b24p
      go to 24
   23 b24p = (1.-b24n)/3.
      b24a = 2.*b24p
   24 continue
c
c..    reverse rates. ground state partition functions assumed.
c
      r23n=0.d0
      r23=0.d0
      r20=0.d0
      if (t9.lt.0.1) go to 25
      r23n=3.93*exp(2.599*rt9)*r24*b24n
      r23=3.93*exp(-2.239*rt9)*r24*b24p
      r20=2.42*exp(-4.6167*rt9)*r24*b24a
c..  
c..   16o+16o rate updated to cf88  q = 16.542
c..   y16*y16*r32 is rate of formation of 32s compound nucleus
c..  
   25 r32 = 7.10e+36/t923*exp(-135.93/t913-0.629*t923-0.445*t943
     1    +0.0103*t9**2)
      r32 = rho*r32*0.5*sc1616
c
c..    branching ratios highly uncertain
c..    guessed using fowler caughlan and zimmerman 1975
c..    deuteron channel is endoergic. apply error function cut-off.
c
      ezro=3.9*t239
      dlt=1.34*t569
      xt=2.*(2.406-ezro)/dlt
      b32d=0.05*thrs(xt)
      b32n=0.1d0
      b32a=0.25d0
      b32p=1.-b32d-b32a-b32n
c
c..    reverse rates. ground state partion functions assumed.
c
      r31=0.d0
      r30=0.d0
      r28=0.d0
      r31n=0.d0
      if (t9.lt.0.1) go to 26
      r31=5.92*exp(-7.676*rt9)*r32*b32p
      r30=0.984*exp(2.412*rt9)*r32*b32d
      r28=3.46*exp(-9.592*rt9)*r32*b32a
      r31n=5.92*exp(-1.448*rt9)*r32*b32n
c
c..       c12 + o16 reaction  q = 16.755
c..   rate revised according to caughlan and fowler 1988
c..   this rate only valid for t9 .gt. 0.5
c..   expression for rate from cf88
c
c..   y(nc12)*y(no16)*rc28 is the rate of formation
c..   of the si28 compound nucleus
c
   26 if (t9.ge.0.5) then
         t9a=t9/(1.+0.055*t9)
         t9a13=t9a**third
         t9a23=t9a13*t9a13
         t9a56=t9a**fivsix
         rc28 = 1.72e+31*t9a56/t932*exp(-106.594/t9a13)
     1        /(exp(-0.18*t9a**2)+1.06e-03*exp(2.562*t9a23))
         rc28 = rho*rc28*sc1216
      else
         rc28 = 0.
      end if
c
c..   branching ratios from pwnsz data
c
      bc27n=.1d0
      bc27p=.5d0
      bc24a=.4d0
c
c..   inverse reaction rates from detailed balance
c..   ground state partition functions assumed
c
      rc27p=0.d0
      rc27n=0.d0
      rc24a=0.d0
      if (t9.lt.0.1) return
      rc27p=1.58*rc28*bc27p*exp(-5.17*rt9)
      rc27n=1.58*rc28*bc27n*exp(0.422*rt9)
      rc24a=2.83*rc28*bc24a*exp(-6.77*rt9)
      return
      end
c..
c..
c..
      function thrs(x)
      implicit real*8 (a-h,o-z)
      save
c..willies threshold fudge function
c..error function rational approx. used (abramowitz p.299)7.1.25
      z=abs(x)
      z2=z*z
      t=1./(1.+0.47047*z)
      t2=t*t
      t3=t2*t
      tt=0.3480242*t-0.0958798*t2+0.7478556*t3
      er=1.-tt*exp(-z2)
c..   er(z) = error function
      if(z)10,11,10
   10 thrs=0.5*(1.-z/x*er)
      return
   11 thrs=0.5d0
      return
      end

      subroutine specl(c12agm)
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'

c..   updated 8-89 to include sevral reactions in cf88 not prev. in spec
c..   major update to expunge one line stmnt func defs 8-89
c..   major update to include light isotope rates pp-c11 12-88
c..   major update by r.d.h. during 3/88 using new cf88 rates.
c..   updated during 11-83 using fcziii by r.d.h.
c..   and 2-86 for o15(ag),ne19(pg):fowler, mg25(pg):parker
c
c..   subroutine generates special rates from fowler
c..   caughlan, and zimmerman 1975 ann. rev. astr. and ap.
c
      real*8    theta,fivsix,third,twoth
      parameter (theta = 0.1, fivsix = 5./6., third = 1./3.,
     1           twoth = 2./3.)
c..
c.. 10/5/95 limit frwd rates to use t9.le.10.
      t9sav=t9
      if (t9.ge.10.) t9=10.

      t9m1=1./t9
      t9m1rev=1./t9sav
      t912=sqrt(t9)
      t9m12=1./t912
      t932=t912*t9
      t9m32=1./t932
      t972=t912*t932*t932
      t913=t9**third
      t923=t913*t913
      t9m13=1./t913
      t9m23=1./t923
      t943=t9*t913
      t953=t9*t923
      t914=sqrt(t912)
      t934=t914*t914*t914
      t954=t914*t9
      t974=t934*t9 
      t915=t9**(1./5.)
      t9m15=1./t915
      t935=t915*t915*t915
      t945=t915*t935
      t965=t9*t915
      t985=t9*t935
      t976=t9**(7./6.)
      t927=t9**(2./7.)
      t947=t927*t927
      t918=t9**(1./8.)
      t938=t918*t918*t918
      t958=t918*t918*t938

      rli7pag=0.
      rb11pa=0.
      rb8epn=0.
      rc12np=0.
      rc12pa=0.
      rc12he3=0.
      rc11na=0.
      rpp=0.
      rpng=0.
      rpgn=0.
      rh2h3n=0.
      rhe3he3=0.
      rh3h3=0.
      rhe3dp=0.
      rbe9pd=0.

      if(nh2.eq.0) go to 19

c..  the fundamental reaction p(p,e+nu)2h  (cf88)
c..    forward rate only; divide by two for identical reactants
c..  q = 1.442  
c.. 
      rpp=rho*4.01e-15*t9m23*exp(-3.38*t9m13)*(1.+0.123*t913
     1 +1.09*t923+0.938*t9)
      rpp=rpp*0.5
c.. 
c..  2h(p,g)3he (cf88)
c..  q = 5.494
c.. 
      term=2.24e+03*t9m23*exp(-3.720*t9m13)*(1.+0.112*t913
     1     +3.38*t923+2.65*t9)
      sig(7,nh2)=rho*term
      sig(8,nh2)=1.63e+10*t932*exp(-63.750*t9m1rev)*term
c..  
c..  2h(n,g)3h schramm and wagoner 1977 ann. rev. nucl. sci. #27, 37
c.. 
      term=66.2*(1.+18.9*t9)
      sig(1,nh2)=rho*term
      sig(2,nh2)=1.63e+10*t932*exp(-72.62*t9m1rev)*term
c.. 
c..  1h(n,g)2h  neutron capture on the proton and inverse
c..  rate expression from schramm and wagoner 1977
c.. 
      term=4.4e+04*(1.-0.86*t912+0.429*t9)
      rpng=rho*term
      rpgn=4.71e+09*t932*exp(-25.82*t9m1rev)*term

   19 if (nh3.eq.0) go to 21

c.. 
c..  h3(p,n)he3 (cf88)
c..  q = -0.764 note: pos q in rev rate exp cancels here
c.. 
      term=7.07e+08*(1.-0.15*t912+0.098*t9)
      sig(3,nh3)=rho*term*exp(-8.863*t9m1)
      sig(4,nh3)=rho*term*0.998
c.. 
c..  h3(p,g)he4 (cf88)
c..  q = 19.814
c.. 
      term=2.20e+04*t9m23*exp(-3.869*t9m13)*(1.+0.108*t913
     1 +1.68*t923+1.26*t9+0.551*t943+1.06*t953)
      sig(7,nh3)=rho*term
      sig(8,nh3)=2.61e+10*t932*exp(-229.932*t9m1rev)*term
c.. 
c..  h2(h3,n)he4 i.e. the "dt" reaction (cf88)
c..  q = 17.589 note: no rev rate, it's term*rho*exp(-204.117/t9)=0.)
c.. 
      term=8.09e+10*t9m23*exp(-4.524*t9m13-(t9/0.120)**2)
     1 *(1.+0.092*t913+1.80*t923+1.16*t9+10.52*t943
     2 +17.24*t953)+8.73e+08*t9m23*exp(-0.523*t9m1)
      rh2h3n=rho*term
c.. 
c..  h3(h3,2n)he4; divide by two for identical reactants (cf88)
c..  q = 11.332 note: no rev rate, it's term*rho*exp(-131.504/t9) =0.)
c.. 
      term=1.67e+09*t9m23*exp(-4.872*t9m13)*(1.+0.086*t913
     1 -0.455*t923-0.272*t9+0.148*t943+0.225*t953)
      rh3h3=rho*term*0.5
c.. 
c..  h3(a,g)li7 (listed as he4(t,g)li7 in cf88)
c..  q = 2.468
c.. 
      term=8.67e+05*t9m23*exp(-8.08*t9m13)
     1 *(1.+0.052*t913-0.448*t923-0.165*t9+0.144*t943
     2 +0.134*t953)
      sig(13,nh3)=rho*term
      sig(14,nh3)=1.11e+10*t932*exp(-28.640*t9m1rev)*term
c.. 
c..  h3(a,n)li6 (listed as he4(t,n)li6 in cf88)
c..  q = -4.782 note: pos q in rev rate exp cancels here
c.. 
      t9ax=t9/(1.+49.18*t9)
      t9ax32=t9ax*sqrt(t9ax)
      term=1.80e+08*(1.-0.261*t9ax32*t9m32)
     1 +2.72e+09*t9m32*exp((55.494-57.884)*t9m1)
      sig(11,nh3)=rho*term*exp(-55.494*t9m1)
      sig(12,nh3)=rho*term*0.935

   21 if (nli6.eq.0) go to 29

c.. 
c..  li6(p,g)7be (cf88)
c..  q = 5.606
      t9ax=t9/(1.-0.0969*t9+0.0284*t953
     1 /(1.-0.0969*t9)**twoth)
      t9ax56=t9ax**fivsix
      t9axm13=1./(t9ax**0.3333333)
      term=6.69e+05*t9ax56*t9m32*exp(-8.413*t9axm13)
      sig(7,nli6)=term*rho
      sig(8,nli6)=term*1.19e+10*t932*exp(-65.054*t9m1rev)
c.. 
c..  li6(a,g)b10 (cf88)
c..  q = 4.460
      term=4.06e+06*t9m23*exp(-18.790*t9m13-(t9/1.326)**2)
     1 *(1.+0.022*t913+1.54*t923+0.239*t9+2.20*t943
     2 +0.869*t953)+1910.*t9m32*exp(-3.484*t9m1)
     3 +1.01e+04*t9m1*exp(-7.269*t9m1)
      sig(13,nli6)=rho*term
      sig(14,nli6)=term*1.58e+10*t932*exp(-51.753*t9m1rev)
c.. 
c..  li6(a,p)be9 (listed as be9(p,a)li6 in cf88)
c..  q = 2.126
c.. 
      term=2.11e+11*t9m23*exp(-10.359*t9m13-(t9/0.52)**2)
     1 *(1.+0.04*t913+1.09*t923+0.307*t9+3.21*t943+2.30
     2 *t953)+4.51e+08*t9m1*exp(-3.046*t9m1)+6.70e+08
     3 *exp(-5.16*t9m1)/t934
      sig(10,nli6)=rho*term
      sig(9,nli6)=rho*term*0.618*exp(-24.674*t9m1rev)

   29 if (nhe3.eq.0) go to 39

c.. 
c..  3he(n,g)4he schramm and wagoner 1977 ann. rev. nucl. sci. #27, 37
c
      term=6.62*(1.+905.*t9)
      sig(1,nhe3)=rho*term
      sig(2,nhe3)=2.61e+10*t932*exp(-238.81*t9m1rev)*term
c.. 
c..   3he(ag)7be (listed as he4(he3,g) in cf88)
c..   q=1.588 
c.. 
      t9ax=t9/(1.+0.0495*t9)
      t9a56x=t9ax**fivsix
      t9am13x=1./(t9ax**0.333333)
      term=5.61e+6*t9a56x*t9m32*exp(-12.826*t9am13x)
      sig(13,nhe3)=rho*term
      sig(14,nhe3)=term*1.11e+10*t932*exp(-18.423*t9m1rev)
c.. 
c..  he3(a,p)li6 and inverse (listed as li6(p,he3)he4 in cf88)
c..  q = 4.018
c.. 
      term=3.73e+10*t9m23*exp(-8.413*t9m13-(t9/5.50)**2)
     1 *(1.+0.05*t913-0.061*t923-0.021*t9+0.006*t943
     2 +0.005*t953)+1.33e+10*t9m32*exp(-17.763*t9m1)
     3 +1.29e+09*t9m1*exp(-21.82*t9m1)
      sig(10,nhe3)=rho*term
      sig(9,nhe3)=rho*term*1.07*exp(-46.631*t9m1rev)
c.. 
c..  he3+he3 makes 2p plus he4 (listed as he3(he3,2p)he4 in cf88)
c..  forward rate only; divide by two for identical reactants
c..  q = 12.860
c.. 
      term=6.04e+10*t9m23*exp(-12.276*t9m13)*(1.+0.034*t913
     1 -0.522*t923-0.124*t9+0.353*t943+0.213*t953)
      rhe3he3=rho*term*0.5
c.. 
c..  3he(d,p)4he; no inverse (cf88)
c..  q = 18.353 
c..  
      term=5.86e+10*t9m23*exp(-7.181*t9m13-(t9/0.315)**2)
     1 *(1.+0.058*t913+0.142*t923+0.0578*t9+2.25*t943
     2 +2.32*t953)+4.36e+08*t9m12*exp(-1.72*t9m1)
      rhe3dp=rho*term
c.. 
c..  3he(t,d)4he;  no inverse (cf88)
c..  q = 14.320
c.. 
      t9ax=t9/(1.+0.128*t9)
      t9ax56=t9ax**fivsix
      t9axm13=1./(t9ax**third)
      term=5.46e+09*t9ax56*t9m32*exp(-7.733*t9axm13)
      rhe3td=rho*term
c.. 
c..  3he(t,np)4he;  no inverse (cf88)
c..  q = 12.096
c.. 
      t9ax=t9/(1.+0.115*t9)
      t9ax56=t9ax**fivsix
      t9axm13=1./(t9ax**third)
      term=7.71e+09*t9ax56*t9m32*exp(-7.733*t9axm13)
      rhe3tnp=rho*term

   39 if(nb11.eq.0) go to 99

c..   
c..   11b(p,n)11c (cf88)
c..   q = -2.764 note: pos q in rev rate exp cancels here
c.. 
      sig(4,nb11)=rho*9.98e-01*1.69e+08*(1.-0.048*t912+0.010*t9)
      sig(3,nb11)=sig(4,nb11)/9.98e-01*exp(-32.080*t9m1)
c..   
c..   11b(p,g)12c  (cf88)
c..   q = 15.957 
c.. 
      b11pg = 4.62e+07*t9m23*exp(-12.095*t9m13-(t9/0.239)**2)
     1 *(1.+0.035*t913+3.0*t923+0.723*t9+9.91*t943+6.07*t953)
     2 +7.89e+03*t9m32*exp(-1.733*t9m1)
     3 +9.68e+04*t9m15*exp(-5.617*t9m1)
      sig(7,nb11)=rho*b11pg
      sig(8,nb11)=b11pg*7.01e+10*t932*exp(-185.173*t9m1rev)
c..   
c..   11b(p,a)2he4 (cf88)
c..   q = 8.628 
c..   special matrix element rb11pa used in spcomp
c.. 
      rb11pa = 2.20e+12*t9m23*exp(-12.095*t9m13-(t9/1.644)**2)
     1 *(1.+0.034*t913+0.14*t923+0.034*t9+0.19*t943+0.116*t953)
     2 +4.03e+06*t9m32*exp(-1.734*t9m1)
     3 +6.73e+09*t9m32*exp(-6.262*t9m1)
     4 +3.88e+09*t9m1*exp(-14.154*t9m1)
      rb11pa=rb11pa*rho
c.. 
c..   currently no reverse rate is in effect (though cf88 provides one)

c.. 
c..   11b(a,p)14c  (cf88)
c..   q = 0.784
c.. 
      sig(9,nb11)=rho*(5.37e+11*t9m23*exp(-28.234*t9m13-(t9/0.347)**2)
     1 *(1.+0.015*t913+5.575*t923+0.576*t9+15.888*t943+4.174*t953)
     2 +5.44e-03*t9m32*exp(-2.827*t9m1)
     3 +3.36e+02*t9m32*exp(-5.178*t9m1)
     4 +5.32e+06/t938*exp(-11.617*t9m1))
      sig(10,nb11)=sig(9,nb11)*1.10e+01*exp(-9.098*t9m1rev)
c.. 
c..   11b(a,n)14n  (cf88)
c..   q = 0.158    
c.. 
      sig(11,nb11)=(6.97e+12*t9m23*exp(-28.234*t9m13-(t9/0.140)**2) 
     1 *(1.+0.015*t913+8.115*t923+0.838*t9+39.804*t943+10.456*t953)
     2 +1.79*t9m32*exp(-2.827*t9m1)
     3 +1.71e+03*t9m32*exp(-5.178*t9m1)
     4 +4.49e+06*t935*exp(-8.596*t9m1))*rho
      sig(12,nb11)=sig(11,nb11)*3.67*exp(-1.835*t9m1rev)
c.. 

   99 if(nli7.eq.0) go to 101

c..   
c..  7li(p,g)8be  (cf88)
c..  q = 17.254
c.. 
      rli7pg = 1.56e+05*t9m23*exp(-8.472*t9m13-(t9/1.696)**2)
     1 *(1.+0.049*t913+2.498*t923+0.86*t9+3.518*t943+3.08*t953)
     2 +1.55e+06*t9m32*exp(-4.478*t9m1)
c.. 
c..  7li(p,n)7be (cf88)
c..  q = -1.644
c.. 
      term=5.15e+09*exp(-1.167*t913-19.081*t9m1)
     1 +7.84e+09*t9m32*exp(-22.832*t9m1)
      sig(3,nli7)=rho*term
      sig(4,nli7)=rho*0.998*(5.15e+09*exp(-1.167*t9sav**(1./3.))
     1 +7.84e+09*t9m32*exp((-22.832+19.081)*t9m1rev))
c.. 
c..   7li(p,a)4he  (cf88)
c..   q = 17.346
c.. 
      t9a1   = t9/(1.+0.759*t9)
      t9a13a = t9a1**third
      t9a56a = t9a1**fivsix
      rli7pa = 1.096e+09*t9m23*exp(-8.472*t9m13)
     1 -4.83e+08*t9a56a*t9m32*exp(-8.472/t9a13a)
     2 +1.06e+10*t9m32*exp(-30.442*t9m1)
c.. 
c..   7lipag = 7lipg + 7lipa  no reverse rates included as yet
c.. 
      rli7pag=rho*(rli7pa+rli7pg)  
c.. 
c..   li7(a,g)b11  (cf88)
c..   q = 8.664
c.. 
      sig(13,nli7)=rho*(3.55e+07*t9m23*exp(-19.161*t9m13-(t9/4.195)**2)
     1 *(1.+0.022*t913+0.775*t923+0.118*t9+0.884*t943+0.342*t953)
     2 +3.33e+02*t9m32*exp(-2.977*t9m1)
     3 +4.10e+04*t9m1*exp(-6.227*t9m1))
      sig(14,nli7)=sig(13,nli7)/rho*4.02e+10*t932*exp(-100.538*t9m1rev)
c..   
c..   li7(a,n)b10  (cf88)
c..   q = -2.790 note: pos q in rev rate exp cancels here
c.. 
      sig(12,nli7) = rho*1.32*3.84e+08 
      sig(11,nli7) = sig(12,nli7)/1.32*exp(-32.382*t9m1)
c..   
c..   be7(e-,nu)li7   weak decay (e.c.) from be7 into lithium 7
c..   7be(e-,nu+g)7li   (cf88)
c..   q = 0.862 
c.. 
      sig(6,nli7)=0.
      if (t9.le.3.0) then
      sig(6,nli7) = 1.34e-10*t9m12*(1.-0.537*t913+3.86*t923
     1 +0.0027*t9m1*exp(2.515e-03*t9m1))
      end if
c..             t9 less than or equal to 3.0
c..  q = 0.049  exclusive of nu energy
c..             rate must not exceed 1.51e-07/(rho*(1.+x)/2.)
c..             for t9 less than 0.001.
c..             we never get that cold!
c..            no rev factor given by cf88

  101 if (nbe7.eq.0) go to 1020

c..   
c..   7be(p,g)8b (cf88)
c..   q = 0.137
c.. 
      be7pg = 3.11e+05*t9m23*exp(-10.262*t9m13)
     1 +2.53e+03*t9m32*exp(-7.306*t9m1)
      sig(7,nbe7)=rho*be7pg      
      sig(8,nbe7)=be7pg*1.30e+10*t932*exp(-1.595*t9m1rev)
c..   
c..   10b(p,a)7be (cf88)
c..   q = 1.146
c.. 
      b10pa = 1.26e+11*t9m23*exp(-12.062*t9m13-(t9/4.402)**2)
     1 *(1.+0.035*t913-0.498*t923-0.121*t9+0.3*t943+0.184*t953)
     2 +2.59e+09*t9m1*exp(-12.260*t9m1)
      sig(10,nbe7)=rho*b10pa
      sig(9,nbe7)=rho*b10pa*7.54e-01*exp(-13.301*t9m1rev)
c.. 
c..  7be(a,g)11c (cf88)
c..  q = 7.544
c.. 
      term=8.45e+07*t9m23*exp(-23.212*t9m13-(t9/4.769)**2)
     1 *(1.+0.018*t913+0.488*t923+0.061*t9+0.296*t943
     2 +0.095*t953)+1.25e+04*t9m32*exp(-6.510*t9m1)
     3 +1.29e+05*exp(-10.039*t9m1)/t954  
      sig(13,nbe7)=rho*term
      sig(14,nbe7)=4.02e+10*t932*exp(-87.539*t9m1rev)*term

 1020 if (nbe9.eq.0) go to 102

c.. 
c..  a(an,g)9be  (cf88)
c.. 
      term=(2.59e-06/((t9*t9)*(1.+0.344*t9)))
     1 *exp(-1.062*t9m1)
      rev=5.84e+19*t9*t9*t9*exp(-18.260*t9m1rev)
      raan=term*rho*rho/2.
      rgaan=rev*term
c.. 
c..  be9(a,n)c12 (Wrean 94)
c..  q = 5.701
c... updated aug 95 - Phys Rev C (1994) vol 49, #2, 1205
c...
      term = 6.476e13*t9m23*exp(-23.8702*t9m13)*(1.0-0.3270*t913) 
     1 + 6.044e-3*t9m32*exp(-1.401*t9m1)+7.268*t9m32*exp(-2.063*t9m1)  
     2 + 3.256e4*t9m32*exp(-3.873*t9m1)+1.946e5*t9m32*exp(-4.966*t9m1)  
     3 + 1.838e9*t9m32*exp(-15.39*t9m1)  
      sig(11,nbe9)=rho*term
      sig(12,nbe9)=rho*term*10.3*exp(-66.160*t9m1rev)
c.. previously 12-5-92 with Kavanaugh & Wrean 1992
c      term = 4.909e13*t9m23*exp(-23.8702*t9m13)*(1.0-0.3266*t913) 
c     1 + 8.002e-4*t9m32*exp(-1.240*t9m1)+5.471*t9m32*exp(-2.037*t9m1)  
c     2 + 2.520e4*t9m32*exp(-3.849*t9m1)+1.749e5*t9m32*exp(-4.966*t9m1)  
c     3 + 1.432e9*t9m32*exp(-15.17*t9m1)  

c.. 
c..  be9(p,d)2 alpha; no inverse (cf88)
c..  q = 0.651
c.. 
      term=2.11e+11*t9m23*exp(-10.359*t9m13-(t9/0.52)**2)
     1 *(1.+0.04*t913+1.09*t923+0.307*t9+3.21*t943+2.30*t953)
     2 +5.79e+08*t9m1*exp(-3.046*t9m1)+8.50e+08
     3 *exp(-5.800*t9m1)/t934
      rbe9pd=rho*term
c.. 
c..  be9(p,g)b10 (cf88)
c..  q = 6.586
c.. 
      term=1.33e+07*t9m23*exp(-10.359*t9m13-(t9/0.846)**2)
     1 *(1.+0.04*t913+1.52*t923+0.428*t9+2.15*t943
     2 +1.54*t953)+9.64e+04*t9m32*exp(-3.445*t9m1)
     3 +2.72e+06*t9m32*exp(-10.62*t9m1)
      sig(7,nbe9)=rho*term
      sig(8,nbe9)=term*9.73e+09*t932*exp(-76.427*t9m1rev)
c.. 
c..  be9(p,n)b9 (cf88)
c..  q = -1.850 note: pos q in rev rate exp cancels here
c.. 
      term=5.58e+07*(1.+0.042*t912+0.985*t9)
     1 +1.02e+09*t9m32*exp((21.473-26.725)*t9m1)
      sig(3,nbe9)=rho*term*exp(-21.473*t9m1)
      sig(4,nbe9)=rho*term*0.998   

  102 if (nb8.eq.0) go to 103

c.. 
c..  b8(ap)c11 (cf88)
c..  q = 7.406
c.. 
      term=1.67e-09*exp(-1.079*t9m1)+9.55e-08*exp(-1.311*t9m1)
     1 +1.98e-01*exp(-2.704*t9m1)+1.34e+00*exp(-4.282*t9m1)
     2 +3.22e+04*exp(-6.650*t9m1)+2.33e+05*exp(-8.123*t9m1)
     3 +2.55e+06*exp(-11.99*t9m1)+9.90e+06*exp(-13.50*t9m1)
     4 +1.41e+06*exp(-16.51*t9m1)+1.99e+07*exp(-18.31*t9m1)
     5 +6.01e+07*exp(-20.63*t9m1)
      term=term*t9m32*rho
      sig(9,nb8)=term
      sig(10,nb8)=3.101*exp(-85.95*t9m1rev)*term
c.. 
c..   positron decay of boron 8 (?)
c.. 
      rb8epn=0.90
c.. 

  103 if(nb10.eq.0) go to 104

c..   
c..  10b(p,g)11c  (cf88)
c..  q = 8.690 
c.. 
      sig(7,nb10)=rho*(4.61e+05*t9m23*exp(-12.062*t9m13-(t9/4.402)**2)
     1 *(1.+0.035*t913+0.426*t923+0.103*t9+0.281*t943+0.173*t953)
     2 +1.93e+05*t9m32*exp(-12.041*t9m1)
     3 +1.14e+04*t9m32*exp(-16.164*t9m1))
      sig(8,nb10) = sig(7,nb10)/rho*3.03e+10*t932*exp(-100.840*t9m1rev)
c.. 
c..  10b(a,n)13n  (cf88) 
c..  q = 1.059
c..  error in matrix elem assignment - fixed 11-96 r.d.h, was 9,nb10, now 11,nb10

      sig(11,nb10) = rho*1.20e+13*t9m23*exp(-27.989*t9m13-(t9/9.589)**2)
      sig(12,nb10) = sig(11,nb10)*9.34*exp(-12.287*t9m1rev)
c.. 
  104 if(nc11.eq.0) go to 105

c..   
c..  n14(p,a)11c caughlan & fowler tnrr5 1988
c..  q = -2.923 note; pos q rev exp cancles here
c.. 
      t9a1=t9/(1.+4.78e-02*t9+7.56e-03*t953/(1.+4.78e-02*t9)**twoth)
      t9a13a = t9a1**third
      t9a56a = t9a1**fivsix
      sig(9,nc11) = 2.72e-01*2.63e+16*t9a56a*t9m32*exp(-31.883/t9a13a)
      sig(9,nc11) = rho*sig(9,nc11)
      sig(10,nc11) = sig(9,nc11)/2.72e-01*exp(-33.915*t9m1)
c.. 
c..  c11(p,g)n12 (cf88)
c..  q = 0.601
c.. 
      term=4.24e+04*t9m23*exp(-13.658*t9m13-(t9/1.627)**2)
     1 *(1.+0.031*t913+3.11*t923+0.665*t9+4.61*t943+2.50*t953)
     2 +8.84e+03*t9m32*exp(-7.021*t9m1)
      sig(7,nc11)=rho*term
      sig(8,nc11)=2.33e+10*t932*exp(-6.975*t9m1rev)*term
c.. 
c..  c11(na)2he4
c..  based upon hauser feshbach calculation by sew on aug 26, 1988
c..  include only reaction and not its inverse. use a constant
c..  rate over entire t range of 7.e+04.
c.. 
      rc11na=7.e+04*rho
c.. 

  105 if(nc12.eq.0) go to 106

c.. 
c..  12c(p,g)13n  fcz2
c..  q = 1.944
c.. 
       sig(7,nc12) = rho*(2.04e+07*t9m23*exp(-13.69*t9m13-(t9/1.5)**2)
     1 *(1.+0.03*t913+1.19*t923+0.254*t9+2.06*t943+1.12*t953)
     2 +1.08e+05*t9m32*exp(-4.925*t9m1)+2.15e+05*t9m32
     3 *exp(-18.179*t9m1))
       sig(8,nc12) = sig(7,nc12)/rho*8.84e+09*t932*exp(-22.553*t9m1rev)
c.. 
c..  12c(a,n)15o fcz3
c..  q = -8.502 note; pos q in rev rate exp cancels here
c.. 
      sig(12,nc12) = rho*1.41*2.48e+07*(1.+0.188*t912+0.015*t9)
      sig(11,nc12) = sig(12,nc12)/1.41*exp(-98.663*t9m1)
c.. 
c..  12c(a,g)16o cf88
c..  q = 7.162
c.. now has multiplier c12agm

      sig(13,nc12)=1.04e+08*t9m1**2/(1.+0.0489*t9m23)**2
     1 *exp(-32.120*t9m13-(t9/3.496)**2)
     2 +1.76e+08*t9m1**2/(1.+0.2654*t9m23)**2*exp(-32.120*t9m13)
     3 +1.25e+03*t9m32*exp(-27.499*t9m1)+1.43e-02*t9**5
     4 *exp(-15.541*t9m1)
      sig(13,nc12)=sig(13,nc12)*rho*c12agm
c..  see comment in cfhz85 
      sig(14,nc12) = sig(13,nc12)/rho*5.13e+10*t932*exp(-83.111*t9m1rev)
c..
c.. oap400 c12(a,g)016 rate
c..
c..    sig(13,nc12)=2.93e+08*t9m1**2/(1.+0.0489*t9m23)**2
c..   1 *exp(-32.120*t9m13-(t9/3.496)**2)
c..   2 +3.14e+08*t9m1**2/(1.+0.2654*t9m23)**2*exp(-32.120*t9m13)
c..   3 +1.25e+03*t9m32*exp(-27.499*t9m1)+1.43e-02*t9**5
c..   4 *exp(-15.541*t9m1)
c..    sig(13,nc12)=sig(13,nc12)*rho
c..    sig(14,nc12) = sig(13,nc12)/rho*5.13e+10*t932*exp(-83.111*t9m1)
c.. 
c..  15n(p,a)12c fcz3
c..  q = 4.966
c.. 
      sig(10,nc12) = rho*(1.08e+12*t9m23*exp(-15.251*t9m13-(t9
     1 /0.522)**2)* (1.+0.027*t913+2.62*t923+0.501*t9+5.36*t943
     2 +2.60*t953)+ 1.19e+08*t9m32*exp(-3.676*t9m1) + 5.41e+08*t9m12
     3 *exp(-8.926*t9m1)+theta*(4.72e+08*t9m32*exp(-7.721*t9m1) +
     4 2.20e+09*t9m32*exp(-11.418*t9m1)))
      sig(9,nc12) = sig(10,nc12)*7.06e-01*exp(-57.625*t9m1rev)

  106 continue
      if(nc13.eq.0) go to 201

c.. 
c..  13c(p,g)14n cf88
c..  q = 7.551
      sig(7,nc13) = rho*(8.01e+07*t9m23*exp(-13.717*t9m13-(t9/2.000)**2)
     1 *(1.+0.030*t913+0.958*t923+0.204*t9+1.39*t943+0.753*t953)
     2 +1.21e+06/t965*exp(-5.701*t9m1))
      sig(8,nc13) = sig(7,nc13)/rho*1.19e+10*t932*exp(-87.621*t9m1rev)
c.. 
c..  13c(a,n)16o fcz2
c..  q = 2.216
c.. 
      sig(11,nc13) = 6.77e+15*t9m23*exp(-32.329*t9m13-(t9/1.284)**2)
     1 *(1.+0.013*t913+2.04*t923+0.184*t9)
     2 +3.82e+05*t9m32*exp(-9.373*t9m1)
     3 +1.41e+06*t9m32*exp(-11.873*t9m1)
     4 +2.0e+09*t9m32*exp(-20.409*t9m1)
     5 +2.92e+09*t9m32*exp(-29.283*t9m1)
      sig(11,nc13) = sig(11,nc13)*rho
      sig(12,nc13) = sig(11,nc13)*5.79e+00*exp(-25.711*t9m1rev)
c.. 
c..  13c(p,n)13n fcz2
c..  q = -3.003 note: pos q in rev rate exp cancels here
c.. 
      sig(4,nc13) = rho*1.88e+08*0.998*(1.-0.167*t912+0.037*t9)
      sig(3,nc13) = sig(4,nc13)/0.998*exp(-34.846*t9m1)

  201 continue
      if(nc14.eq.0) go to 202

c.. 
c..  14c(p,g)15n cf88
c..  q = 10.207
c.. 
      sig(7,nc14) = 6.80e+06*t9m23*exp(-13.741*t9m13-(t9/5.721)**2)
     1 *(1.+0.03*t913+0.503*t923+0.107*t9+0.213*t943+0.115*t953)
     2 +5.36e+03*t9m32*exp(-3.811*t9m1)
     3 +9.82e+04*t9m13*exp(-4.739*t9m1)
      sig(7,nc14) = sig(7,nc14)*rho
      sig(8,nc14) = sig(7,nc14)/rho*9.00e+09*t932*exp(-118.452*t9m1rev)
c.. 
c..  14c(p,n)14n   cf88
c..  q = -0.626 note: pos q in rev rate exp cancels here
c.. 
      sig(4,nc14) = rho*(3.33e-01*7.19e+05
     1 *(1.+0.361*t912+0.502*t9)
     2 +3.33e-01*3.34e+08*t9m12*exp(-4.983*t9m1))
      sig(3,nc14) = rho*(7.19e+05*(1.+0.361*t912+0.502*t9)
     1 *exp(-7.263*t9m1)+3.34e+08*t9m12*exp(-12.246*t9m1))
c.. 
c..  14c(a,g)18o cf88
c..  q = 6.227
c.. 
      sig(13,nc14)=rho*(3.375e+08*t9m1**2*exp(-32.513*t9m13)
     1 +1.528e+09*t9m23*exp(-32.513*t9m13-(t9/2.662)**2)
     2 *(1.+0.0128*t913-0.869*t923-0.0779*t9+0.321*t943
     3 +0.0732*t953)+9.29e-08*t9m32*exp(-2.048*t9m1)
     4 +2.77e+03/t945*exp(-9.876*t9m1))
      sig(14,nc14)=sig(13,nc14)/rho*5.42e+10*t932*exp(-72.262*t9m1rev)

  202 continue
       if(n13.eq.0) go to 203

c.. 
c..  13n(p,g)14o cf88
c..  q = 4.628 
c.. 
      sig(7,n13) = rho*(4.04e+07*t9m23*exp(-15.202*t9m13-(t9/1.191)**2)
     1 *(1.+0.027*t913-0.803*t923-0.154*t9+5.00*t943+2.44*t953)
     2 +2.43e+05*t9m32*exp(-6.348*t9m1))
      sig(8,n13) = sig(7,n13)/rho*3.57e+10*t932*exp(-53.706*t9m1rev)
c.. 
c..  low q
c..  16o(p,a)13n fcz3
c..  q = -5.218 note; pos q in rev rate exp cancels here
c.. 
      t9a=t9/(1.+7.76e-02*t9+2.64e-02*t953/(1.+7.76e-02*t9)
     1**twoth)
      t9a13a=t9a**third
      t9a56a=t9a**fivsix
      sig(10,n13)=rho*(1.88e+18*t9a56a*t9m32*exp(-35.829
     1 /t9a13a-60.557*t9m1))
      sig(9,n13)=rho*(1.88e+18*t9a56a*t9m32*1.72e-01*exp(-35.829
     1/t9a13a))

  203 continue
      if(n14.eq.0) go to 204

c.. 
c..  14n(p,n)14o fcz2
c..  q = -5.925 note: pos q in rev rate exp cancels here
c.. 
      sig(4,n14) = rho*2.99*6.74e+07*(1.+0.658*t912+0.379*t9)
      sig(3,n14) = sig(4,n14)/2.99*exp(-68.762*t9m1)
c.. 
c..  14n(p,g)15o cf88
c..  q= 7.297
c.. 
      sig(7,n14)=rho*(4.90e+07*t9m23*exp(-15.228*t9m13-(t9/3.294)**2)
     1 *(1.+0.027*t913-0.778*t923-0.149*t9+0.261*t943+0.127*t953)
     2 +2.37e+03*t9m32*exp(-3.011*t9m1)+2.19e+04*exp(-12.530*t9m1))
      sig(8,n14)=sig(7,n14)/rho*2.70e+10*t932*exp(-84.678*t9m1rev)

c..  14n(a,n)17f cf88
c..  q = -4.735 note: pos q in rev rate exp cancels here

      term1=5.24e+09*(1.-1.15*t912+0.365*t9)
     1 *exp(-(t9/2.798)**2) 
      term2=3.28e+10*t9m32*exp(-1.5766e+01*t9m1)
      sig(11,n14)=rho*(term1*exp(-54.942*t9m1)
     1                +term2*exp(-54.942*t9m1))
      sig(12,n14)=rho*(term1+term2)*1.48e+00
c.. 
c..  14n(a,g)18f fcz2
c..  q = 4.415
c.. 
      sig(13,n14)=rho*(7.78e+09*t9m23*exp(-36.031*t9m13-(t9/0.881)**2)
     1 *(1.+0.012*t913+1.45*t923+0.117*t9+1.97*t943+0.406*t953)
     2 +2.36e-10*t9m32*exp(-2.798*t9m1)+2.03*t9m32*exp(-5.054*t9m1)
     3 +1.15e+04*t9m23*exp(-12.31*t9m1))
      sig(14,n14)=sig(13,n14)/rho*5.42e+10*t932*exp(-51.236*t9m1rev)
c.. 
c..  17o(p,a)14n fcz2
c..  q = 1.191
c.. 
      sig(10,n14)=rho*(1.53e+07*t9m23*exp(-16.712*t9m13-(t9/0.565)**2)
     1 *(1.+0.025*t913+5.39*t923+0.94*t9+13.5*t943+5.98*t953)
     2 +2.92e+06*t9*exp(-4.247*t9m1)
     3 +theta*(4.81e+10*t9*exp(-16.712*t9m13-(t9/0.04)**2)
     4 +5.05e-05*t9m32*exp(-0.723*t9m1))
     5 +theta*(13.1*t9m32*exp(-1.961*t9m1)))
      sig(9,n14)=sig(10,n14)*0.676*exp(-13.825*t9m1rev)

  204 continue
      if(n15.eq.0)  go to 205
c.. 
c..  15n(p,g)16o fcz2
c..  q = 12.128
c.. 
      sig(7,n15) = rho*(9.78e+08*t9m23*exp(-15.251*t9m13-(t9/0.45)**2)
     1 *(1.+0.027*t913+0.219*t923+0.042*t9+6.83*t943+3.32*t953)
     2 +1.11e+04*t9m32*exp(-3.328*t9m1)
     3 +1.49e+04*t9m32*exp(-4.665*t9m1)
     4 +3.8e+06*t9m32*exp(-11.048*t9m1))
      sig(8,n15) = sig(7,n15)/rho*3.62e+10*t932*exp(-140.734*t9m1rev)
c.. 
c..  15n(p,n)15o cf88
c..  q = -3.563 note: pos q in rev rate exp cancels here
c.. 
      sig(4,n15) = rho*0.998*3.51e+08*(1.+0.452*t912-0.191*t9)
      sig(3,n15) = sig(4,n15)/0.998*exp(-41.037*t9m1)
c.. 
c..  15n(a,n)18f cf88
c..  q = -6.418 note: pos q in rev rate exp cancels here

      term=3.14e+08*(1.-0.641*t912+0.108*t9)
      sig(11,n15)=rho*term*exp(-74.479*t9m1)
      sig(12,n15)=rho*term*2.0
c.. 
c..  15n(a,g)19f cf88
c..  q = 4.104
c.. 
      sig(13,n15)=rho*(2.54e+10*t9m23*exp(-36.211/t913-(t9/0.616)**2)
     1 *(1.+0.012*t913+1.69*t923+0.136*t9+1.91*t943+0.391*t953)
     2 +9.83e-03*t9m32*exp(-4.232*t9m1)+1.52e+03*t9*exp(-9.747*t9m1))
      sig(14,n15)=sig(13,n15)/rho*5.54e+10*t932*exp(-46.578*t9m1rev)
c.. 
c..  18o(p,a)15n cf88
c..  q = 3.980
c.. 
      sig(10,n15)=rho*(3.63e+11*t9m23*exp(-16.729*t9m13-(t9/1.361)**2)
     1 *(1.+0.025*t913+1.88*t923+0.327*t9+4.66*t943+2.06*t953)
     2 +9.90e-14*t9m32*exp(-0.231*t9m1)+2.66e+04*t9m32*exp(-1.670*t9m1)
     3 +2.41e+09*t9m32*exp(-7.638*t9m1)+1.46e+09*t9m1*exp(-8.310*t9m1))
      sig(9,n15)=sig(10,n15)*1.66e-01*exp(-46.191*t9m1rev)

  205 continue
      if (no14.eq.0) go to 4

c.. 
c..  14o(a,g)18ne cf88
c..  q = 5.112
c.. 
      sig(13,no14)=rho*(9.47e+08*t9m23*exp(-39.388*t9m13-(t9/0.717)**2)
     1 *(1.+0.011*t913+1.974*t923+0.146*t9+3.036*t943
     2 +0.572*t953)+1.16e-01*t9m32*exp(-11.733*t9m1)
     3 +3.39e+01*t9m32*exp(-22.609*t9m1)+9.10e-03*t9**5
     4 *exp(-12.159*t9m1))
      sig(14,no14)=sig(13,no14)/rho*5.42e+10*t932*exp(-59.328*t9m1rev)
c.. 
c..  14o(a,p)17f cf88
c..  q = 1.191
c.. 
      sig(9,no14)=rho*(1.68e+13*t9m23*exp(-39.388*t9m13-(t9/0.717)**2)
     1 *(1.+0.011*t913+13.117*t923+0.971*t9+85.295*t943
     2 +16.06*t953)+3.31e+04*t9m32*exp(-11.733*t9m1)
     3 +1.79e+07*t9m32*exp(-22.609*t9m1)+9.00e+03*t9**(11./3.)
     4 *exp(-12.517*t9m1))
      sig(10,no14)=sig(9,no14)*4.93e-01*exp(-13.820*t9m1rev)

    4 continue
      if(no15.eq.0) go to 5

c.. 
c..  15o(a,g)19ne cf88
c..  q = 3.529
c.. 
      sig(13,no15)=rho*(3.57e+11*t9m23*exp(-39.584*t9m13-(t9/3.000)**2)
     1 *(1.+0.011*t913-0.273*t923-0.020*t9)
     2 +5.10e+10*t9m23*exp(-39.584*t9m13-(t9/1.937)**2)
     3 *(1.+0.011*t913+1.59*t923+0.117*t9+1.81*t943+0.338*t953)
     4 +3.95e-01*t9m32*exp(-5.849*t9m1)
     5 +1.90e+01*t9**2.85*exp(-7.356*t9m1-(t9/8.000)**2))
      sig(14,no15)=sig(13,no15)/rho*5.54e+10*t932*exp(-40.957*t9m1rev)
c..  
c..  18f(p,a)15o  wiescher and kettner, ap. j., 263, 891 (1982)
c..  last three resonancesare very uncertain
c.. 
      sig(10,no15)=rho*(t9m32*(1.66e-10*exp(-0.302*t9m1)
     1 +1.56e+05*exp(-3.84*t9m1)+1.36e+06*exp(-5.22*t9m1)
     2 +8.1e-05*exp(-1.05*t9m1)+8.9e-04*exp(-1.51*t9m1)
     3 +3.0e+05*exp(-4.29*t9m1)))
      sig(9,no15)=sig(10,no15)*4.93e-01*exp(-3.3433e+01*t9m1rev)

    5 continue
      if(no16.eq.0) go to 10

c.. 
c..  16o(p,g)17f fcz2
c..  q = 0.0600
c.. 
      sig(7,no16)=rho*(1.50e+08/(t923*(1.+2.13*(1.-exp(-0.728
     1 *t923))))*exp(-16.692*t9m13))
      sig(8,no16)=sig(7,no16)/rho*3.03e+09*t932*exp(-6.968*t9m1rev)
c.. 
c..  16o(a,g)20ne cf88
c..  q = 4.734
c.. 
      sig(13,no16)=rho*(9.37e+09*t9m23*exp(-39.757*t9m13-(t9/1.586)**2)
     1 +6.21e+01*t9m32*exp(-10.297*t9m1)
     2 +5.38e+02*t9m32*exp(-12.226*t9m1)
     3 +1.30e+01*t9**2*exp(-20.093*t9m1))
      sig(14,no16)=sig(13,no16)/rho*5.65e+10*t932*exp(-54.937*t9m1rev)
c.. 
c..  19f(p,a)16o fcz3
c..  q = 8.114
c.. 
      gt9x=(1.+4.*exp(-2.090*t9m1)+7.*exp(-16.440*t9m1))
      sig(10,no16)=rho*(3.55e+11*t9m23*exp(-18.113*t9m13-(t9/0.845)**2)
     1 *(1.+0.023*t913+1.96*t923+0.316*t9+2.86*t943+1.17*t953)
     2 +3.67e+06*t9m32*exp(-3.752*t9m1) + 3.07e+08*exp(-6.019*t9m1))
     3 /gt9x
      sig(9,no16)=sig(10,no16)*6.54e-01*exp(-9.4159e+01*t9m1rev)
 
   10 continue
      if(no17.eq.0) go to 20

c.. 
c..  17o(p,g)18f fcz2
c..  q = 5.607
c.. 
      t9a=t9/(1.+2.69*t9)
      t9a13x=t9a**third
      t9a56x=t9a**fivsix
      sig(7,no17) = rho*(7.97e+07*t9a56x*t9m32*exp(-16.712/t9a13x)
     1 +1.51e+08*t9m23*exp(-16.712*t9m13)
     2 *(1.+0.025*t913-0.051*t923-8.82e-03*t9)
     3 +1.56e+05*t9m1*exp(-6.272*t9m1)
     4 +theta*(13.1*t9m32*exp(-1.961*t9m1)))
      sig(8,no17)=sig(7,no17)/rho*3.66e+10*t932*exp(-65.061*t9m1rev)
c.. 
c..  17o(a,g)21ne fcz2
c..  q = 7.351
c.. 
      t9a=t9/(1.+0.1646*t9)
      t9a13x=t9a**third
      t9a56x=t9a**fivsix
      ft9a=exp(-(0.786/t9a)**3.51)
      fpt9a=exp(-(t9a/1.084)**1.69)
      gt9x=1.+exp(-10.106*t9m1)/3.
      sig(13,no17)=1.73e+17*fpt9a/gt9x*t9a56x*t9m32*exp(-39.914/t9a13x)
     1 +3.50e+15*ft9a/gt9x*t9a56x*t9m32*exp(-39.914/t9a13x)
      sig(13,no17) = sig(13,no17)*rho
      sig(14,no17)=sig(13,no17)/rho*8.63e+10*t932*exp(-85.305*t9m1rev)
c.. 
c..  17o(a,n)20ne fcz2
c..  q = 0.590
c.. 
      t9a=t9/(1.+0.0268*t9+0.0232*t953/(1.+0.0268*t9)**twoth)
      t9a56x=t9a**fivsix
      t9a13x=t9a**third
      gt9x=1.+exp(-10.106*t9m1)/3.
      sig(11,no17)=rho*(1.03e+18/gt9x*t9a56x*t9m32*exp(-39.914/t9a13x))
      sig(12,no17)=sig(11,no17)*1.86e+01*exp(-6.852*t9m1rev)

   20 continue
      if(no18.eq.0) go to 50

c.. 
c..  18o(p,g)19f cf88
c..  q = 7.994
c..  
      sig(7,no18)=rho*(3.45e+08*t9m23*exp(-16.729*t9m13-(t9/0.139)**2)
     1 *(1.+0.025*t913+2.26*t923+0.394*t9+30.56*t943+13.55*t953)
     2 +1.25e-15*t9m32*exp(-0.231*t9m1) 
     3 +1.64e+02*t9m32*exp(-1.670*t9m1)
     4 +1.28e+04*t912*exp(-5.098*t9m1))
      sig(8,no18)=sig(7,no18)/rho*9.20e+09*t932*exp(-92.769*t9m1rev)
c.. 
c..  18o(a,g)22ne fcz3
c..  q = 9.669
c.. 
      sig(13,no18) = 1.82e+12*t9m23*exp(-40.057*t9m13-(t9/0.343)**2)
     1 *(1.+0.01*t913+0.988*t923+0.072*t9+3.17*t943+0.586*t953)
     2 +7.54*t9m32*exp(-6.228*t9m1)+34.8*t9m32*exp(-7.301*t9m1)
     3 +6.23e+03*t9*exp(-16.987*t9m1)
     4 +theta*1.0e-11*t9m32*exp(-1.994*t9m1)
      sig(13,no18) = sig(13,no18)*rho
      sig(14,no18)=sig(13,no18)/rho*5.85e+10*t932*exp(-112.208*t9m1rev)
c.. 
c..  18o(a,n)21ne cf88
c..  q = -0.693 note, small q value allows for pos exp in rev rate
c.. 
      gt91 = 1.+5.*exp(-23.002*t9m1)
      t9a = t9/(1.+0.0483*t9+0.00569*t953/(1.+0.0483*t9)**twoth)
      ft9a = exp(-(0.431/t9a)**3.89)
      sig(11,no18) = rho*(7.22e+17*ft9a/gt91*t9a**fivsix/t932
     1 *exp(-40.056/t9a**third)+150.31/gt91*exp(-8.045*t9m1))
      if(t9.lt.0.04) go to 50
      sig(12,no18) = sig(11,no18)*0.784*exp(8.101*t9m1rev)

   50 continue
      if (nf17.eq.0) go to 71

c.. 
c..  17f(p,g)18ne  wiescher and kettner, ap. j., 263, 891 (1982)
c.. 
      sig(7,nf17) = rho*(1.66e+07*t9m23*exp(-18.03*t9m13)*(2.194
     1 +0.050*t913-0.376*t923-0.061*t9+0.026*t943+0.011*t953)
     2 +839.0*t9m32*exp(-6.93*t9m1)+33.56*t9m32*exp(-7.75*t9m1))
      sig(8,nf17)=sig(7,nf17)/rho*1.087e+11*t932*exp(-45.501*t9m1rev)
c.. 
c..  20ne(p,a)17f cf88
c..  q = -4.134 note: pos q in rev rate exp cancels here

      t9a=t9/(1.+6.12e-02*t9+1.30e-02*t953
     1      /(1.+6.12e-02*t9)**twoth)
      t9a13=t9a**third
      t9a56=t9a**fivsix
      term=3.25e+19*(5.31+0.544*t9-0.0523*t9*t9)
     1     *t9a56*t9m32*exp(-43.176/t9a13)
      sig(10,nf17)=rho*term*exp(-47.969*t9m1)
      sig(9,nf17)=rho*term*5.37e-02

   71 continue
      if (nf18.eq.0) go to 72

c.. 
c..  18f(p,g)19ne  wiescher and kettner, ap. j., 263, 891 (1982)
c..  last resonance is very uncertain
c.. 
      sig(7,nf18) = rho*(1.658e+7*t9m23*exp(-18.06*t9m13)*(4.604
     1 +0.106*t913+0.053*t923+0.009*t9-0.036*t943-0.015*t953)
     2 +t9m32*(4.55e-14*exp(-0.302*t9m1)+327.*exp(-3.84*t9m1)
     3 +1.32e+04*exp(-5.22*t9m1)+93.*exp(-4.29*t9m1)))
      sig(8,nf18)=sig(7,nf18)/rho*2.73e+10*t932*exp(-7.4396e+01*t9m1rev)

   72 continue
      if(nf19.eq.0) go to 80

c.. 
c..  19f(p,g)20ne fcz2
c..  q = 12.848
c.. 
      gt9x=1.+4.0*exp(-2.09*t9m1)+7.0*exp(-16.44*t9m1)
      sig(7,nf19)=rho*(6.04e+07*t9m23*exp(-18.113*t9m13-(t9/0.416)**2)
     1 *(1.+0.023*t913+2.06*t923+0.332*t9+3.16*t943+1.30*t953)
     2 +6.32e+02*t9m32*exp(-3.752*t9m1)+7.56e+04/t927*exp(-5.722*t9m1))
      sig(7,nf19)=sig(7,nf19)/gt9x
      sig(8,nf19)=sig(7,nf19)/rho*3.7e+10*t932*exp(-149.093*t9m1rev)
c.. 
c..  19f(p,n)19ne cf88 
c..  q = -4.021 note: pos q in rev rate exp cancels here
c.. 
      sig(4,nf19) = rho*1.27e+08*0.998*(1.-0.147*t912+0.069*t9)
      sig(3,nf19) = sig(4,nf19)/0.998*exp(-46.659*t9m1)
c.. 
c..  19f(a,p)22ne cf88
c..  q = 1.675
c.. 
      sig(9,nf19) = rho*(4.50e+18*t9m23*exp(-43.467*t9m13
     1 -(t9/0.637)**2)+7.98e+04*t932*exp(-12.760*t9m1))
      sig(10,nf19) = sig(9,nf19)*6.36*exp(-19.435*t9m1rev)
c.. 
c..  na22(n,a)f19 cf88
c..  q = 1.949
c.. 
      sig(12,nf19)=rho*(1.21e+06*exp(1.+8.955e-01*t9-5.645e-02*t9**2
     1        +7.302e-04*t9**3))
      sig(11,nf19)=sig(12,nf19)*1.10*exp(-22.620*t9m1rev)

   80 continue
      if (ne19.eq.0) go to 95

c.. 
c..  19ne(p,g)20na cf88
c..  q = 2.199
c.. 
      sig(7,ne19)=rho*(1.71e+06*t9m23*exp(-19.431*t9m13)*(1.+0.021*t913
     1 +0.130*t923+1.95e-02*t9+3.86e-02*t943+1.47e-02*t953)
     2 +1.89e+05*t9m23*exp(-19.431*t9m13-(t9/1.142)**2)
     3 *(1.+0.021*t913+2.13*t923+0.320*t9+2.80*t943+1.07*t953)
     4 +8.45e+03/t954*exp(-7.64*t9m1))
      sig(8,ne19)=sig(7,ne19)/rho*7.39e+09*t932*exp(-25.519*t9m1rev)

   95 continue
      if(ne20.eq.0) go to 100

c.. 
c..  20ne(p,g)21na fcz2
c..  q = 2.431
c.. 
      sig(7,ne20) = 9.55e+06*exp(-19.447*t9m13)
     1 /(t9*t9*(1.+0.0127*t9m23)**2)
     1 +2.05e+08*t9m23*exp(-19.447*t9m13)*(1.+2.67*exp(-sqrt(t9/0.21)))
     2 +18.0*t9m32*exp(-4.242*t9m1)+10.2*t9m32*exp(-4.607*t9m1)
     3 +3.6e+04/t914*exp(-11.249*t9m1)
      sig(7,ne20) = sig(7,ne20)*rho
      sig(8,ne20) = sig(7,ne20)/rho*4.63e+09*t932*exp(-28.216*t9m1rev)
c.. 
c..  20ne(a,g)24mg fcz2
c..  q = 9.312
c.. 
      gt9x=1.+5.*exp(-18.96*t9m1)
      sig(13,ne20) = 4.11e+11*t9m23*exp(-46.766*t9m13-(t9/2.219)**2)
     1 *(1.+9.0e-3*t913+0.882*t923+0.055*t9+0.749*t943+0.119*t953)
     2 +5.27e+03*t9m32*exp(-15.869*t9m1)
     3 +6.51e+03*t912*exp(-16.223*t9m1)
     4 +theta*(42.1*t9m32*exp(-9.115*t9m1)+32.0*t9m23*exp(-9.383*t9m1))
      sig(13,ne20) = sig(13,ne20)*rho/gt9x
      sig(14,ne20)=sig(13,ne20)/rho*6.01e+10*t932*exp(-108.059*t9m1rev)
c.. 
c..  23na(p,a)20ne fcz3
c..  q = 2.379
c.. 
      sig(10,ne20)=rho*(8.56e+09*t9m23*exp(-20.766*t9m13-(t9/0.131)**2)
     1 *(1.+0.02*t913+8.21*t923+1.15*t9+44.36*t943+15.84*t953)
     2 +4.02*t9m32*exp(-1.99*t9m1)+1.18e+04/t954*exp(-3.148*t9m1)
     3 +8.59e+05*t943*exp(-4.375*t9m1)
     4 +theta*3.06e-12*t9m32*exp(-0.447*t9m1))
      sig(9,ne20)=sig(10,ne20)*1.25*exp(-27.606*t9m1rev)

  100 continue
      if(ne21.eq.0) go to 110

c.. 
c..  21ne(p,g)22na cf88
c..  q = 6.738
c..  
      sig(7,ne21)=rho*(4.37e+08*t9m23*exp(-19.462*t9m13)+5.85*t9m32
     1 *exp(-1.399*t9m1)+1.29e+04*t9m32*exp(-3.009*t9m1)+3.15e+05/t935
     2 *exp(-5.763*t9m1)+theta*(2.95e+08*t9m23*exp(-19.462*t9m13
     3 -(t9/0.058)**2)*(1.+0.021*t913+13.29*t923+1.99*t9+124.1
     4 *t943+47.29*t953)+7.80e-01*t9m32*exp(-1.085*t9m1)))
      sig(8,ne21)=sig(7,ne21)/rho*1.06e+10*t932*exp(-78.194*t9m1rev)
c.. 
c..  21ne(a,g)25mg fcz3
c..  q = 9.882
c.. 
      t9ax=t9/(1.+0.0537*t9)
      t9a13x=t9ax**third
      t9a56x=t9ax**fivsix
      tmult=1.52e-04*exp(-46.90*t9m13*(8.72e-03*t9-6.87e-04*t9**2
     1 +2.15e-05*t9**3)) / (1.+1.5*exp(-4.068*t9m1)
     2 +2.0*exp(-20.258*t9m1))
      sig(13,ne21) = 4.94e+19*t9a56x*t9m32*exp(-46.89/t9a13x)
     1 +2.66e+07*t9m32*exp(-22.049*t9m1)
      sig(13,ne21) = sig(13,ne21)*rho*tmult
      sig(14,ne21) = sig(13,ne21)/rho*4.06e+10*t932*
     1 exp(-114.676*t9m1rev)
c.. 
c..  21ne(a,n)24mg fcz2
c..  q = 2.551
c.. 
      gt9x=1.+1.5*exp(-4.068*t9m1)+2.0*exp(-20.258*t9m1)
      t9ax=t9/(1.+0.0537*t9)
      t9a13x=t9ax**third
      t9a56x=t9ax**fivsix
      sig(11,ne21) = 4.94e+19*t9a56x*t9m32*exp(-46.89/t9a13x)
     1 +2.66e+07*t9m32*exp(-22.049*t9m1)
      sig(11,ne21) = sig(11,ne21)*rho/gt9x
      sig(12,ne21)=sig(11,ne21)*12.9*exp(-29.606*t9m1rev)

  110 continue
      if(ne22.eq.0) go to 120

c.. 
c..  22ne(p,g)23na fcz3
c..  q = 8.794
c.. 
      sig(7,ne22)=rho*(1.15e+09*t9m23*exp(-19.475*t9m13)
     1 +9.77e-12*t9m32*exp(-0.348*t9m1)
     2 +8.96e+03*t9m32*exp(-4.84*t9m1)
     3 +6.52e+04*t9m32*exp(-5.319*t9m1)
     4 +7.97e+05*t9m12*exp(-7.418*t9m1)
     5 +theta*1.63e-01*t9m32*exp(-1.775*t9m1))
      sig(8,ne22)=sig(7,ne22)/rho*4.67e+09*t932*exp(-102.048*t9m1rev)
c.. 
c..  22ne(a,g)26mg fcz2
c..  q = 10.612
c.. 
      t9a=t9/(1.+0.0548*t9)
      t9a13x=t9a**third
      t9a56x=t9a**fivsix
      gt9x=1.+5.0*exp(-14.791*t9m1)
      ft9a=exp(-(0.197/t9a)**4.82)
      fpt9a=exp(-(t9a/0.249)**2.31)
      sig(13,ne22) = rho*(4.16e+19*fpt9a/gt9x*t9a56x*t9m32
     1 *exp(-47.004/t9a13x)+2.08e+16*ft9a/gt9x
     2 *t9a56x*t9m32*exp(-47.004/t9a13x))
      sig(14,ne22)=sig(13,ne22)/rho*6.15e+10*t932*exp(-123.151*t9m1rev)
c.. 
c..  22na(n,p)22ne cf88
c..  q = 3.624
c..  
      sig(4,ne22)=rho*(1.24e+08*exp(1.-3.037e-02*t9+8.380e-03*t9**2
     1        -7.101e-04*t9**3))
      sig(3,ne22)=sig(4,ne22)*7.01*exp(-42.059*t9m1rev)
c.. 
c..  22ne(a,n)25mg cf88
c..  q = -0.481
c.. 
      t9ax=t9/(1.+0.0548*t9)
      t9a13x=t9ax**third
      t9a56x=t9ax**fivsix
      ft9a=exp(-(0.197/t9ax)**4.82)
      gt9x=1.+5.0*exp(-14.791*t9m1)
      sig(11,ne22) = rho*(4.16e+19*ft9a/gt9x*t9a56x*t9m32
     1 *exp(-47.004/t9a13x)+1.44e-04/gt9x*exp(-5.577*t9m1))

c..   following is to prevent crashes on vax at low t

      sig(12,ne22)=7.83e-05  
      if(t9.lt.0.008) go to 120
      sig(12,ne22)=sig(11,ne22)*5.44e-01*exp(5.577*t9m1rev)

120   continue
      if(na21.eq.0) go to 124

c.. 
c..  21na(p,g)22mg cf88
c..  q = 5.497
c.. 
      sig(7,na21)=rho*(1.41e+05*t9m23*exp(-20.739*t9m13-(t9/0.366)**2)
     1 *(1.+0.020*t913+4.741*t923+0.667*t9+16.380*t943
     2 +5.858*t953)+6.72e+02/t934*exp(-2.436*t9m1))
      sig(8,na21)=sig(7,na21)/rho*7.44e+10*t932*exp(-63.790*t9m1rev)
c.. 
c..  24mg(p,a)21na cf88
c..  q = -6.880 note: pos q in rev rate exp cancels here
c..  see comment hfcz83

      t9a=t9/(1.+0.127*t9)
      t9a13=t9a**third
      t9a56=t9a**fivsix
      term=1.81e+21*(4.43+3.31*t9-0.229*t9*t9)
     1     *t9a56*t9m32*exp(-49.967/t9a13)
      sig(10,na21)=rho*term*exp(-79.843*t9m1)
      sig(9,na21) =rho*term*7.71e-02

  124 continue
      if(na22.eq.0) go to 128

c.. 
c..  22na(p,g)23mg cf88
c..  q = 7.578
c.. 
      sig(7,na22)=rho*(9.63e-05*t932*exp(-0.517*t9m1)+2.51e+04*t9
     1 *exp(-2.013/t9))
      sig(8,na22)=sig(7,na22)/rho*3.27e+10*t932*exp(-87.933*t9m1rev)

  128 continue
      if(na23.eq.0) go to 130

c.. 
c..  23na(p,g)24mg fcz2
c..  q = 11.691
c.. 
      gt9x = 1.+1.5*exp(-5.105*t9m1)
      sig(7,na23) = 2.93e+08*t9m23*exp(-20.766*t9m13-(t9/0.297)**2)
     1 *(1.+0.02*t913+1.61*t923+0.226*t9+4.94*t943+1.76*t953)
     2 +9.34e+01*t9m32*exp(-2.789*t9m1)+1.89e+04*t9m32*exp(-3.434*t9m1)
     3 +5.1e+04*t915*exp(-5.51*t9m1)
      sig(7,na23) = sig(7,na23)/gt9x*rho
      sig(8,na23)=sig(7,na23)/rho*7.49e+10*t932*exp(-135.665*t9m1rev)
c.. 
c..   23na(p,n)23mg cf88
c..   q = -4.841 note: pos q in rev rate exp cancels here
c.. 
      t9a=t9/(1.+0.141*t9)
      t9a32x=t9a**(3./2.)
      sig(4,na23)=rho*(9.29e+08*(1.-0.881*t9a32x*t9m32))
      sig(3,na23)=sig(4,na23)*0.998*exp(-56.173*t9m1)
c.. 
c..  26al(n,a)23na t cf88
c..  q = 2.968
c..  old rate multiplied by 0.3, fit param adjusted in bdat882.dat
c.. 

  130 continue
      if(mg24.eq.0) go to 150

c.. 
c..  24mg(p,g)25al fcz3
c..  q = 2.271
c.. 
      gt9x=(1.+5.*exp(-15.882*t9m1))
      sig(7,mg24)=rho*(5.60e+08*t9m23*exp(-22.019*t9m13)*(1.+0.019*t913
     1 -0.173*t923-0.023*t9)+1.48e+03*t9m32*exp(-2.484*t9m1)
     2 +4.00e+03*exp(-4.180*t9m1))/gt9x
      sig(8,mg24)=sig(7,mg24)/rho*3.13e+09*t932*exp(-26.358*t9m1rev)
c.. 
c..  24mg(a,g)28si fcz3
c..  q = 9.984
c.. 
      gt9x=(1.+5.0*exp(-15.882*t9m1))
      sig(13,mg24)=rho*(4.78e+01*t9m32*exp(-13.506*t9m1) 
     1 +2.38e+03*t9m32*exp(-15.218*t9m1)
     2 +2.47e+02*t932*exp(-15.147*t9m1)
     3 +theta*(1.72e-09*t9m32*exp(-5.028*t9m1)
     4 +1.25e-03*t9m32*exp(-7.929*t9m1)
     5 +2.43e+01*t9m1*exp(-11.523*t9m1)))
      sig(13,mg24)=sig(13,mg24)/gt9x
      sig(14,mg24)=sig(13,mg24)/rho*6.27e+10*t932*exp(-115.862*t9m1rev)
c.. 
c..  27al(p,a)24mg cf88
c..  q = 1.600
c.. 
      gt9x=1.+exp(-9.792*t9m1)/3.0+2.*exp(-11.773*t9m1)/3.
      sig(10,mg24)=1.1e+08*t9m23*exp(-23.261*t9m13-(t9/0.157)**2)
     1 *(1.+0.018*t913+12.85*t923+1.61*t9+89.87*t943+28.66*t953)
     2 +1.29e+02*t9m32*exp(-2.517*t9m1)+5.66e+03*t972*exp(-3.421*t9m1)
     3 +theta*(3.89e-08*t9m32*exp(-0.853*t9m1)
     4 +8.18e-09*t9m32*exp(-1.001*t9m1))
      sig(10,mg24)=sig(10,mg24)*rho/gt9x
      sig(9,mg24)=sig(10,mg24)*1.81*exp(-18.572*t9m1rev)

  150 continue
      if(mg25.eq.0) go to 160


c...
c... 25mg(p,g)26alt cf88
c... q = 6.306
c... altered sept. 28, 1990 to adjust for measurements of Iliadis
c...   et al. resonance terms divided by 90., 4.5, and 9.
c...   respectively. This fits their data to better than 20%.
c... reduction is by about 4.5 at 20 to 70 million K. No effect
c...   above 0.3 billion K.
c...
      sig(7,mg25)=rho*(3.57e+09*t9m23*exp(-22.031*t9m13-(t9/0.06)**2)
     1 *(1.+0.019*t913+7.669*t923+1.015*t9+167.4*t943+56.35*t953)
     2 +3.07e-13*t9m32*exp(-0.435*t9m1)/90.
     3 +1.94e-07*t9m32*exp(-0.673*t9m1)/4.5
     4 +3.15e-05*t9m1**(3.40)*exp(-1.342*t9m1-(t9/13.)**2)/9.
     5 +1.77e+04*t958*exp(-3.049*t9m1-(t9/13.)**2))
      sig(8,mg25)=sig(7,mg25)/rho*1.03e+10*t932*exp(-73.183*t9m1rev)
c.. 
c..  25mg(a,p)28al cf88
c..  q = -1.206 note: pos q in rev rate exp cancels here
      term=3.23e+08*t9m23*exp(-23.271*t9m13+6.46*t9
     1    -2.39*t9*t9+0.506*t9*t9*t9
     2    -6.04e-02*t9*t9*t9*t9
     3    +3.75e-03*t9*t9*t9*t9*t9
     4    -9.38e-05*t9*t9*t9*t9*t9*t9)
      sig(9,mg25)=rho*term*exp(-13.995*t9m1)
      sig(10,mg25)=rho*term*2.86
c.. 
c..  25mg(a,g)29si fcz3
c..  q = 11.127
c.. 
      t9a=t9/(1.+0.0630*t9)
      t9a13x=t9a**third
      t9a56x=t9a**fivsix
      gt9x=(1.+10.*exp(-13.180*t9m1)/3.)
      sig(13,mg25)=rho*(3.59e+20/gt9x*t9a56x*t9m32*exp(-53.41/t9a13x)
     1 *5.87e-04*exp(-53.42*t9m13*(0.0156*t9-1.79e-03*t9**2+9.08e-05
     2 *t9**3)))
      sig(14,mg25)=sig(13,mg25)/rho*1.90e+11*t932*exp(-129.139*t9m1rev)
c.. 
c..  25mg(a,n)28si fcz2
c..  q = 2.653
c.. 
      gt9x=1.+10.*exp(-13.18*t9m1)/3.0
      t9ax=t9/(1.+0.063*t9)
      t9a13x=t9ax**third
      t9a56x=t9ax**fivsix
      sig(11,mg25)= rho*(3.59e+20/gt9x*t9a56x*t9m32*exp(-53.41/t9a13x))
      sig(12,mg25)=sig(11,mg25)*20.0*exp(-30.792*t9m1rev)

  160 continue
      if(mg26.eq.0) go to 170

c.. 
c..  26mg(p,g)27al cf88
c..  q = 8.272
c.. 

      gt9x = (1.+5.*exp(-20.990*t9m1))
      sig(7,mg26) = 7.39e+08*t9m23*exp(-22.042*t9m13-(t9/0.299)**2)
     1 *(1.+0.019*t913+3.61*t923+0.478*t9+9.78*t943+3.29*t953)
     2 +1.32e-10*t9m32*exp(-0.603*t9m1)
     3 +2.90e-05*t9m32*exp(-1.056*t9m1)
     4 +6.45e-05*t9m32*exp(-1.230*t9m1)
     5 +5.64e-02*t9m32*exp(-1.694*t9m1)
     6 +2.86e+03*t9m32*exp(-3.265*t9m1)
     7 +7.99e+04*t9m32*exp(-3.781*t9m1)
     8 +4.23e+04*t912*exp(-3.661*t9m1)
      sig(7,mg26) = sig(7,mg26)*rho/gt9x
      sig(8,mg26) = sig(7,mg26)/rho*3.14e+09*t932*exp(-95.99*t9m1rev)
c.. 
c..  26mg(a,g)30si fcz3
c..  q = 10.644
c.. 
      gt9x=(1.+5.*exp(-20.990*t9m1))
      t9a=t9/(1.+0.0628*t9)
      t9a13x=t9a**third
      t9a56x=t9a**fivsix
      sig(13,mg26)=rho*(2.93e+20/gt9x*t9a56x*t9m32*exp(-53.505/t9a13x)
     1 *4.55e-02*exp(-53.51*t9m13*(0.0751*t9-0.0105*t9**2
     2 +5.57e-04*t9**3)))
      sig(14,mg26)=sig(13,mg26)/rho*6.38e+10*t932*exp(-123.52*t9m1rev)
c.. 
c..  26mg(a,n)si29 fcz2
c..  q = 0.035
c.. 
      t9a=t9/(1.+6.28e-02*t9)
      t9a13x=t9a**third
      t9a56x=t9a**fivsix
      gt9x=1.+5.*exp(-20.99*t9m1)
      sig(11,mg26)=rho*(2.93e+20/gt9x*t9a56x*t9m32*exp(-53.505/t9a13x))
      sig(12,mg26)=sig(11,mg26)*1.68*exp(-0.401*t9m1rev)
c.. 
c..  26al(n,p)26mg t cf88
c..   updated as specified by caughlan & fowler tnrr5 1988
c..  q = 4.786
c..   essentially identical to old fcz expression, leave it alone (s.e.w.)

  170 continue
      if(nal26.eq.0) go to 175

c.. 
c..  26al(p,g)27si cf88
c..  q = 7.464
c.. 
      sig(7,nal26)=rho*(6.78e+13*t9m23*exp((-23.261*t9m13)
     1 *(1.+2.004e-01*t9-1.538e-02*t9**2+5.723e-04*t9**3))
     2 +6.13e+02*t9m32*exp(-3.202*t9m1)+9.45e+03*t9m1*exp(-4.008*t9m1))
      sig(8,nal26)=sig(7,nal26)/rho*3.46e+10*t932*exp(-86.621*t9m1rev)

  175 continue
      if(nal27.eq.0) go to 180

c.. 
c..  27al(p,g)28si cf88
c..  q = 11.585
c..  
      gt9h=(1.+exp(-9.792*t9m1)/3.+2.*exp(-11.773*t9m1)/3.)
      sig(7,nal27)=rho*(1.67e+08*t9m23*exp(-23.261*t9m13-(t9/0.155)**2)
     1 *(1.+0.018*t913+5.81*t923+0.728*t9+27.31*t943+8.71*t953)
     2 +2.20e+00*t9m32*exp(-2.269*t9m1)+1.22e+01*t9m32*exp(-2.491*t9m1)
     3 +1.50e+04*t9*exp(-4.112*t9m1)+theta
     4 *(6.50e-10*t9m32*exp(-0.853*t9m1)+1.63e-10*t9m32
     5 *exp(-1.001*t9m1)))/gt9h
      sig(8,nal27)=sig(7,nal27)/rho*1.13e+11*t932*exp(-134.434*t9m1rev)
c.. 
c..  27al(a,n)30p cf88
c..  q = -2.636 note: pos q in rev rate exp cancels here
c.. 
      sig(12,nal27)=6.75*rho*(8.20e+04+5.21e+05*t974
     1 *exp(-2.966*t9m1)) 
      sig(11,nal27)=rho*(8.2e+04*exp(-30.588*t9m1)
     1 +5.21e+05*t974*exp(-33.554*t9m1))

  180 continue
      if (nsi27.eq.0) go to 186

c.. 
c..  27si(p,g)p28 cf88
c..  q = 2.065
c.. 
      sig(7,nsi27)=rho*(1.64e+09*t9m23*exp(-24.439*t9m13)
     1 +2.00e-08*t9m32*exp(-0.928*t9m1)+1.95e-02*t9m32*exp(-1.857*t9m1)
     2 +3.70e+02/t947*exp(-3.817*t9m1))
      sig(8,nsi27)=sig(7,nsi27)/rho*1.62e+10*t932*exp(-23.960*t9m1rev)

  186 continue
      if (nsi28.eq.0) go to 190

c.. 
c..  28si(p,g)29p cf88
c..  q = 2.747
c.. 
      sig(7,nsi28)=rho*(1.64e+08*t9m23*exp(-24.449*t9m13-(t9/2.91)**2)
     1 *(1.+0.017*t913-4.11*t923-0.491*t9+5.22*t943+1.58*t953)
     2 +3.52e+02*t9m32*exp(-4.152*t9m1)+6.3e+05*t9m32*exp(-18.505*t9m1)
     3 +1.69e+03*exp(-14.518*t9m1))
      sig(8,nsi28)=sig(7,nsi28)/rho*t932*9.46e+09*exp(-31.879*t9m1rev)

  190 continue
      if (nsi29.eq.0) go to 195

c.. 
c..  29si(p,g)30p cf88
c..  q = 5.601
c.. 
      sig(7,nsi29)=rho*(3.26e+09*t9m23*exp(-24.459*t9m13-(t9/0.256)**2)
     1 *(1.+0.017*t913+4.27*t923+0.509*t9+15.40*t943+4.67*t953)
     2 +2.98e+03*t9m32*exp(-3.667*t9m1)+3.94e+04*t9m32*exp(-4.665*t9m1)
     3 +2.08e+04*t912*exp(-8.657*t9m1))
      sig(8,nsi29)=sig(7,nsi29)/rho*1.26e+10*t932*exp(-65.002*t9m1rev)
  195 continue

      if (nsi30.eq.0) go to 196

c..  30si(p,g)31p cf88
c..  q = 7.297

      term=4.25e+08*t9m23*exp(-24.468*t9m13
     1    -(t9/0.670)**2)
     2    *(1.+0.017*t913+0.150*t923+0.018*t9
     3    +5.53*t943+1.68*t953)
     4    +1.86e+04*t9m32*exp(-5.601*t9m1)
     5    +3.15e+05*t9m32*exp(-6.961*t9m1)
     6    +2.75e+05*t9m12*exp(-10.062*t9m1)
      sig(7,nsi30)=rho*term
      sig(8,nsi30)=term*9.50e+09*t932*exp(-84.673*t9m1rev)

  196 continue 
c.. reset t9 to it's original value
      t9=t9sav
      return
      end

      subroutine qma0(w,u,v,y,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1),y(1)

      do 100 k=1,n

      w(k) = ( u(k) * v(k) ) + y(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qvdot(e,u,v,n)
      implicit real*8 (a-h,o-z)
      save

      dimension u(1),v(1)

      e=0.
      do 100 k=1,n

      e = e + u(k) * v(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qma1(w,u,v,c,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1)

      do 100 k=1,n

      w(k) = ( u(k) * v(k) ) + c

100   continue
      return
      end
c..
c..
c..
      subroutine qma2(w,u,c,y,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),y(1)

      do 100 k=1,n

      w(k) = ( u(k) * c ) + y(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qma3(w,u,c1,c2,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1)

      do 100 k=1,n

      w(k) = ( u(k) * c1 ) + c2

100   continue
      return
      end
c..
c..
c..
      subroutine qvmpy0(w,u,v,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1)
      nk=n
      do 100 k=1,nk

      w(k) =  u(k) * v(k)

100   continue
      nk=0
      return
      end
c..
c..
c..
      subroutine qvmpy1(w,u,c,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1)

      do 100 k=1,n

      w(k) =  u(k) * c

100   continue
      return
      end
c..
c..
c..
      subroutine qjsub0(m,i,j,n)
      implicit real*8 (a-h,o-z)
      save

      dimension m(1),i(1),j(1)

      do 100 k=1,n

      m(k) =  i(k) - j(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qvsub0(w,u,v,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1)

      do 100 k=1,n

      w(k) =  u(k) - v(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qvsub2(w,c,v,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),v(1)

      do 100 k=1,n

      w(k) =  c - v(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qmm0(w,u,v,y,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1),y(1)

      do 100 k=1,n

      w(k) = ( u(k) * v(k) ) * y(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qmm2(w,u,c,y,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),y(1)

      do 100 k=1,n

      w(k) = ( u(k) * c ) * y(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qam0(w,u,v,y,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1),y(1)

      do 100 k=1,n

      w(k) = ( u(k) + v(k) ) * y(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qaa0(w,u,v,y,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1),y(1)

      do 100 k=1,n

      w(k) = ( u(k) + v(k) ) + y(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qmr1(w,u,v,c,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),u(1),v(1)

      do 100 k=1,n

      w(k) =  u(k) - ( v(k) * c)

100   continue
      return
      end
c..
c..
c..
      subroutine qvset(c,w,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1)

      do 100 k=1,n

      w(k) = c

100   continue
      return
      end
c..
c..
c..
      subroutine qvseti(kc,ka,n)
      implicit real*8 (a-h,o-z)
      save

      dimension ka(1)

      do 100 k=1,n

      ka(k) = kc

100   continue
      return
      end
c..
c..
c..
      subroutine qvcopy(v,w,n)
      implicit real*8 (a-h,o-z)
      save

      dimension v(1),w(1)

      do 100 k=1,n

      w(k) = v(k)

100   continue
      return
      end
c..
c..
c..
      subroutine qvsums(e,v,n)
      implicit real*8 (a-h,o-z)
      save

      dimension v(1)

      do 100 k=1,n

      e = v(k) + e

100   continue
      return
      end
c..
c..
c..
      subroutine qvflot(w,i,n)
      implicit real*8 (a-h,o-z)
      save

      dimension w(1),i(1)

      do 100 k=1,n

      w(k) = dble( i(k) )

 100  continue
      return
      end
c..
c..
c..
      subroutine blockcopy(s,t,n)
      implicit real*8 (a-h,o-z)
      save

      dimension s(1),t(1)

      do 100 k=1,n

      t(k)=s(k)

 100  continue

      return
      end
c..
c..
c..
c..
c..
      subroutine compt(inu)
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'
c..
c..this routine forms the right hand side of the reaction equations and their
c..jacobian. actually the negative of the jacobian.
c..
c..declare
      integer          itot,iat,nt,i,j,k,in,ip,ia
      double precision a1,a2,a3,a4,a5,a6,b1,b2
c..
c..zero the a matrix and the b vector
      nt = 0
      do 01 i=1,nzo
       a(i) = 0.0d0
01    continue
      itot = i2 + 3
      do 02 i=1,itot
       b(i) = 0.0d0
02    continue
c..
c..set the n,p, alpha pointers       
      in = i2 + 1
      ip = i2 + 2
      ia = i2 + 3
c..
c..for every isotope
      do 11 j=1,i2
c..
c..(n,g) reactions
       k = nrr(1,j)
       if (k .ne. 0) then
        a1     = sig(1,j)*aan
        a2     = sig(2,j)
        a3     = sig(1,j)*y(j)
        b1     = -aan*sig(1,j)*y(j) + sig(2,j)*y(k)
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        b(j)   = b(j) + b1
        b(in)  = b(in) + b1
        b(k)   = b(k) - b1
       end if
c..
c..(p,n) and for beta-, beta+ decay
       k = nrr(2,j)
       if (k .ne. 0) then
        a1     = sig(3,j) * scfacp(j) * aap
        a2     = sig(4,j) * scfacp(j) * aan
        a3     = sig(3,j) * scfacp(j) * y(j)
        a4     = sig(4,j) * scfacp(j) * y(k)
        a5     = sig(5,j)
        a6     = sig(6,j)
        b1     = -aap*sig(3,j)*scfacp(j)*y(j)
     1           + aan*sig(4,j)*scfacp(j)*y(k)
        b2     = -y(j)*sig(5,j)+y(k)*sig(6,j)
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1 + a5
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2 - a6
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2 + a6
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1 - a5
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        b(j)   = b(j) + b1 + b2
        b(k)   = b(k) - b1 - b2
        b(ip)  = b(ip) + b1
        b(in)  = b(in) - b1
       end if
c..
c..(p,gam); do not screen (g,p)
       k = nrr(3,j)
       if (k .ne. 0) then 
        a1     = sig(7,j) * scfacp(j) * aap
c..        a2    = sig(8,j) * scfacp(j)
        a2     = sig(8,j)
        a3     = sig(7,j) * scfacp(j) * y(j)
c..        b1=-aap*sig(7,j)*scfacp(j)*y(j)+sig(8,j)*scfacp(j)*y(k)
        b1     = -aap * sig(7,j) * scfacp(j) * y(j) + sig(8,j) * y(k)
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        b(j)   = b(j) + b1
        b(k)   = b(k) - b1
        b(ip)  = b(ip) + b1
       end if
c..
c..(alp,p) reactions
       k = nrr(4,j)
       if (k .ne. 0) then 
        a1     = sig(9,j)*scfaca(j)*aaa
        a2     = sig(10,j)*scfaca(j)*aap
        a3     = sig(9,j)*scfaca(j)*y(j)
        a4     = sig(10,j)*scfaca(j)*y(k)
        b1     = -aaa*sig(9,j)*scfaca(j)*y(j)
     1           + aap*sig(10,j)*scfaca(j)*y(k)
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1 
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        b(j)   = b(j) + b1
        b(ia)  = b(ia) + b1
        b(k)   = b(k) - b1
        b(ip)  = b(ip) - b1
       end if
c..
c..(alp,n) reactions
       k = nrr(5,j)
       if (k .ne. 0) then
        a1     = sig(11,j)*scfaca(j)*aaa
        a2     = sig(11,j)*scfaca(j)*y(j)
        a3     = sig(12,j)*scfaca(j)*aan
        a4     = sig(12,j)*scfaca(j)*y(k)
        b1     = -aaa*sig(11,j)*scfaca(j)*y(j)
     1           + aan*sig(12,j)*scfaca(j)*y(k)
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1 
        iat    = eloc(nt)
        a(iat) = a(iat) - a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a4
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a4
        nt     = nt + 1 
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1 
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        b(j)   = b(j) + b1
        b(ia)  = b(ia) + b1
        b(k)   = b(k) - b1
        b(in)  = b(in) - b1
       end if
c..
c..(a,gam) reactions; do not screen (g,a)
       k = nrr(6,j)
       if (k .ne. 0) then
        a1     = sig(13,j)*scfaca(j)*aaa
        a2     = sig(13,j)*scfaca(j)*y(j)
c..        a3     = sig(14,j)*scfaca(j)
        a3     = sig(14,j)
c..        b1     = -aaa*sig(13,j)*scfaca(j)*y(j)+sig(14,j)*scfaca(j)*y(k)
        b1     = -aaa*sig(13,j)*scfaca(j)*y(j)+sig(14,j)*y(k)
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1 
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        b(j)   = b(j) + b1
        b(ia)  = b(ia) + b1
        b(k)   = b(k) - b1
       end if
11    continue
c..
c..
c..if we are doing neutrino interactions, add in those 
      if (inu .ne. 0) then
       do 110 j=1,i2
c..
c..(nu, e-) reaction
        k = nrrneut(1,j)
        if (k .ne. 0) then 
         a1     = signuc(1,j)
         b1     = signuc(1,j)*y(j)
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) + a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         b(j)   = b(j) - b1
         b(in)  = b(in) + b1
         b(k)   = b(k) + b1
        end if
c..
c..(nu, e-) reaction
        k = nrrneut(3,j)
        if (k .ne. 0) then 
         a1     = signuc(2,j)
         b1     = signuc(2,j)*y(j)
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) + a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         b(j)   = b(j) - b1
         b(ip)  = b(ip) + b1
         b(k)   = b(k) + b1
        end if
c..
c..(nu, e+) reaction
        k = nrrneut(4,j)
        if (k .ne. 0) then
         a1     = signuc(3,j)
         b1     = signuc(3,j)*y(j)
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) + a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         b(j)   = b(j) - b1
         b(in)  = b(in) + b1
         b(k)   = b(k) + b1
        end if
c..
c..(nu, e+) reaction
        k = nrrneut(6,j)
        if (k .ne. 0) then 
         a1     = signuc(4,j)
         b1     = signuc(4,j)*y(j)
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) + a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         nt     = nt + 1
         iat    = eloc(nt)
         a(iat) = a(iat) - a1
         b(j)   = b(j) - b1
         b(ip)  = b(ip) + b1
         b(k)   = b(k) + b1
        end if
110    continue
      end if
c..
c..now add in the special reactions; this was routine spcomp
c..
c..p(e-,nu)n and n(e+,nub)p reaction
      nt     = nt + 1
      iat    = eloc(nt)
      a(iat) = a(iat) + rpen
      nt     = nt + 1
      iat    = eloc(nt)
      a(iat) = a(iat) - rnep
      nt     = nt + 1
      iat    = eloc(nt)
      a(iat) = a(iat) - rpen
      nt     = nt + 1
      iat    = eloc(nt)
      a(iat) = a(iat)+rnep
      b(in)  = b(in) - aan*rnep + aap*rpen
      b(ip)  = b(ip) - aap*rpen + aan*rnep
c..
c..if nc12 is 0 do not include triple alpha or c burning
      if (nc12 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 9.0*aaa*ra3*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 3.0*ral
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat)+ral
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 3.0*aaa*ra3*aaa
       b(ia)   = b(ia) - 3.0*aaa*ra3*(aaa**2) + 3.0*ral*y(nc12)
       b(nc12) = b(nc12) + aaa*ra3*(aaa**2) - ral*y(nc12)
c..
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 4.0*y(nc12)*r24
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r20*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r23*aap
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r20*y(ne20)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r23*y(na23)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r20*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r24*y(nc12)*2.0*b24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r20*y(ne20)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r20*y(ne20)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r24*y(nc12)*2.0*b24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r20*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23*aap
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r24*y(nc12)*2.0*b24p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23*y(na23)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23*y(na23)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r24*y(nc12)*2.0*b24p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23*aap
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r23n*aan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r23n*y(mg23)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23n*aan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r24*y(nc12)*2.0*b24n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23n*y(mg23)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23n*y(mg23)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r24*y(nc12)*2.0*b24n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r23n*aan
       b(nc12) = b(nc12) - 2.0*y(nc12)*r24*y(nc12)
     1           + 2.0*aaa*r20*y(ne20)
     2           + 2.0*aap*r23*y(na23) 
     3           + 2.0*aan*r23n*y(mg23)
       b(ne20) = b(ne20) - aaa*r20*y(ne20) + y(nc12)*r24*y(nc12)*b24a
       b(ia)   = b(ia) - aaa*r20*y(ne20) + y(nc12)*r24*y(nc12)*b24a
       b(na23) = b(na23) - aap*r23*y(na23) + y(nc12)*r24*y(nc12)*b24p
       b(ip)   = b(ip) - aap*r23*y(na23) + y(nc12)*r24*y(nc12)*b24p
       b(mg23) = b(mg23) - aan*r23n*y(mg23) + y(nc12)*r24*y(nc12)*b24n
       b(in)   = b(in) - aan*r23n*y(mg23) + y(nc12)*r24*y(nc12)*b24n
      end if
c..
c..o16 burning
      if (no16 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 4.0*r32*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r28*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r31*aap
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r31n*aan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r28*y(nsi28)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r31*y(np31)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.0*r31n*y(ns31)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r28*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r28*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32a*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32a*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r28*y(nsi28)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r28*y(nsi28)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31*aap
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31*aap
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32p*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32p*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31*y(np31)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31*y(np31)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32d*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32d*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32d*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31n*aan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31n*aan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32n*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - r32*b32n*2.0*y(no16)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31n*y(ns31)
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + r31n*y(ns31)
       b(nsi28)=b(nsi28) + y(no16)*r32*b32a*y(no16) - aaa*r28*y(nsi28)
       b(ia)   = b(ia) + y(no16)*r32*b32a*y(no16) - aaa*r28*y(nsi28)
       b(np31) = b(np31) + y(no16)*r32*b32p*y(no16) - aap*r31*y(np31)
       b(ip)   = b(ip) + y(no16)*r32*b32p*y(no16) - aap*r31*y(np31)
       b(ns31) = b(ns31) + y(no16)*r32*b32n*y(no16) - aan*r31n*y(ns31)
       b(in)   = b(in) + y(no16)*r32*b32n*y(no16) - aan*r31n*y(ns31)
       b(np30) = b(np30) + y(no16)*r32*b32d*y(no16)
       b(in)   = b(in) + y(no16)*r32*b32d*y(no16)
       b(ip)   = b(ip) + y(no16)*r32*b32d*y(no16)
       b(no16) = b(no16) + 2.0*(aaa*r28*y(nsi28)+aap*r31*y(np31)
     1           + aan*r31n*y(ns31) - y(no16)*r32*y(no16))
      end if
c..
c..c12+o16 reaction
      if (nc12 .ne. 0 .and. no16 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(no16)*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nc12)*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aap*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aan*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aaa*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nal27)*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nsi27)*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(mg24)*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nc12)*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(no16)*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aap*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aan*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aaa*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nal27)*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nsi27)*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(mg24)*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aaa*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(no16)*bc24a*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nc12)*bc24a*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(mg24)*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aap*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(no16)*bc27p*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nc12)*bc27p*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nal27)*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aan*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(no16)*bc27n*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nc12)*bc27n*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nsi27)*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(mg24)*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aaa*rc24a
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(no16)*bc24a*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nc12)*bc24a*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nal27)*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aap*rc27p
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(no16)*bc27p*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nc12)*bc27p*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nsi27)*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aan*rc27n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nc12)*bc27n*rc28
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(no16)*bc27n*rc28
       b(nc12)  = b(nc12) - y(nc12)*y(no16)*rc28 + aap*y(nal27)*rc27p
     1            + aan*y(nsi27)*rc27n + aaa*y(mg24)*rc24a
       b(no16)  = b(no16) - y(nc12)*y(no16)*rc28+aap*y(nal27)*rc27p
     1            + aan*y(nsi27)*rc27n + aaa*y(mg24)*rc24a
       b(mg24)  = b(mg24) - aaa*y(mg24)*rc24a 
     1            + y(nc12)*y(no16)*bc24a*rc28
       b(nal27) = b(nal27) - aap*y(nal27)*rc27p
     1            + y(nc12)*y(no16)*bc27p*rc28
       b(nsi27) = b(nsi27) - aan*y(nsi27)*rc27n
     1            + y(nc12)*y(no16)*bc27n*rc28
       b(ia)    = b(ia) - aaa*y(mg24)*rc24a + y(nc12)*y(no16)*bc24a*rc28
       b(ip)    = b(ip) - aap*y(nal27)*rc27p 
     1            + y(nc12)*y(no16)*bc27p*rc28
       b(in)    = b(in) - aan*y(nsi27)*rc27n 
     1            + y(nc12)*y(no16)*bc27n*rc28
      end if
c..
c..7li(pa) reaction
      if (nli7 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aap*rli7pag
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nli7)*rli7pag
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aap*rli7pag
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nli7)*rli7pag
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*aap*rli7pag
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*y(nli7)*rli7pag
       b(nli7) = b(nli7) - aap*y(nli7)*rli7pag
       b(ip)   = b(ip) - aap*y(nli7)*rli7pag
       b(ia)   = b(ia) + 2.*aap*y(nli7)*rli7pag
      end if
c..
c..11b(pa) reaction
      if (nb11 .ne. 0) then
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + aap*rb11pa
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + y(nb11)*rb11pa
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + aap*rb11pa
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + y(nb11)*rb11pa
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - 3.*aap*rb11pa
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - 3.*y(nb11)*rb11pa
       b(nb11) = b(nb11) - aap*y(nb11)*rb11pa
       b(ip)   = b(ip) - aap*y(nb11)*rb11pa
       b(ia)   = b(ia) + 3.*aap*y(nb11)*rb11pa
      end if
c..
c..8b(e+ nu) reactin
      if (nb8 .ne. 0) then
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + rb8epn
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - 2.*rb8epn
       b(nb8) = b(nb8)-y(nb8)*rb8epn
       b(ia)  = b(ia)+2.*y(nb8)*rb8epn
      end if
c..
c..12c(nu,nup n p) reaction
      if (nb10 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + rc12np
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12np
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12np
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12np
       b(nc12) = b(nc12) - y(nc12)*rc12np
       b(nb10) = b(nb10) + y(nc12)*rc12np
       b(ip)   = b(ip) + y(nc12)*rc12np
       b(in)   = b(in) + y(nc12)*rc12np
      end if
c..
c..12c(nu,nup p a) reaction
      if (nli7 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + rc12pa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12pa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12pa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12pa
       b(nc12) = b(nc12) - y(nc12)*rc12pa
       b(nli7) = b(nli7) + y(nc12)*rc12pa
       b(ip)   = b(ip) + y(nc12)*rc12pa
       b(ia)   = b(ia) + y(nc12)*rc12pa
      end if
c..
c..12c(nu,nup n a) reaction
      if (nbe7 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + rc12na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12na
       b(nc12) = b(nc12) - y(nc12)*rc12na
       b(nbe7) = b(nbe7) + y(nc12)*rc12na
       b(in)   = b(in) + y(nc12)*rc12na
       b(ia)   = b(ia) + y(nc12)*rc12na
      end if
c..
c..c12(nu nup he3 n 2a) reaction
      if (nhe3 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + rc12he3n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12he3n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rc12he3n
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*rc12he3n
       b1      = y(nc12)*rc12he3n
       b(nc12) = b(nc12)-b1
       b(nhe3) = b(nhe3)+b1
       b(in)   = b(in)+b1
       b(ia)   = b(ia)+2.*b1
c..
c..12c(nu,nup he3)9be reaction
       if (nbe9 .ne. 0) then
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) + rc12he3
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - rc12he3
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - rc12he3
        b(nc12) = b(nc12) - y(nc12)*rc12he3
        b(nbe9) = b(nbe9) + y(nc12)*rc12he3
        b(nhe3) = b(nhe3) + y(nc12)*rc12he3
       end if
      end if
c..
c..9be(p,d)2 alpha reaction
      if (nbe9 .ne. 0 .and. nh2 .ne. 0) then
       a1      = aap*rbe9pd
       a2      = y(nbe9)*rbe9pd
       b1      = y(nbe9)*aap*rbe9pd
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a2
       nt      = nt + 1 
       iat     = eloc(nt)
       a(iat)  = a(iat) + a2
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - a2
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - a1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*a2
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*a1
       b(nbe9) = b(nbe9) - b1
       b(ip)   = b(ip) - b1
       b(nh2)  = b(nh2) + b1
       b(ia)   = b(ia) + 2.*b1
      end if
c..
c..11c(na) reaction
      if (nc11 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aan*rc11na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nc11)*rc11na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aan*rc11na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nc11)*rc11na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 3.*aan*rc11na
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 3.*y(nc11)*rc11na
       b(nc11) = b(nc11) - aan*y(nc11)*rc11na
       b(in)   = b(in) - aan*y(nc11)*rc11na
       b(ia)   = b(ia) + 3.*aan*y(nc11)*rc11na
      end if
c..
c..a(an,g)9be reaction
      if (nbe9 .ne. 0) then
       b1      = aaa*aaa*aan*raan - y(nbe9)*rgaan
       b(ia)   = b(ia) - 2.*b1
       b(in)   = b(in) - b1
       b(nbe9) = b(nbe9) + b1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 4.*aaa*aan*raan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*rgaan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 2.*aaa*aaa*raan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 2.*aaa*aan*raan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - rgaan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + aaa*aaa*raan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*aaa*aan*raan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + rgaan
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - aaa*aaa*raan
      end if
c..
c..proton-proton reaction
      if (nh2 .ne. 0) then
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + 4.*aap*rpp
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - 2.*aap*rpp
       b(ip)  = b(ip) - 2.*aap*aap*rpp
       b(nh2) = b(nh2) + aap*aap*rpp
c..
c..neutron capture on the proton
       a1     = rpng*aan
       a2     = rpgn
       a3     = rpng*aap
       b1     = -aan*rpng*aap + rpgn*y(nh2)
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + a1
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - a2
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + a3
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + a2
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - a1
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - a3
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + a3
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) + a1
       nt     = nt + 1
       iat    = eloc(nt)
       a(iat) = a(iat) - a2
       b(ip)  = b(ip) + b1
       b(in)  = b(in) + b1
       b(nh2) = b(nh2) - b1
c..
c..the "dt" reaction  2h(3h,n)4he
       if (nh3 .ne. 0) then
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + y(nh3)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + y(nh2)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + y(nh2)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + y(nh3)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - y(nh3)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - y(nh2)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - y(nh3)*rh2h3n
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) -y(nh2)*rh2h3n
        b(nh2) = b(nh2) - y(nh2)*y(nh3)*rh2h3n
        b(nh3) = b(nh3) - y(nh2)*y(nh3)*rh2h3n
        b(ia)  = b(ia) + y(nh2)*y(nh3)*rh2h3n
        b(in)  = b(in) + y(nh2)*y(nh3)*rh2h3n
c..
c..h3(p,g)he4   and inverse including neutrino irradiation
        a1     = sig(7,nh3)*aap
        a2     = sig(8,nh3)
        a3     = sig(7,nh3)*y(nh3)
        b1     = -aap*sig(7,nh3)*y(nh3) + sig(8,nh3)*aaa
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a2
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + a1
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - a2
        b(nh3) = b(nh3) + b1
        b(ia)  = b(ia) - b1
        b(ip)  = b(ip) + b1
c..
c..t(t,2n)he4; no inverse included
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) + 4.*y(nh3)*rh3h3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - 4.*y(nh3)*rh3h3
        nt     = nt + 1
        iat    = eloc(nt)
        a(iat) = a(iat) - 2.*y(nh3)*rh3h3
        b(nh3) = b(nh3) - 2.*y(nh3)*y(nh3)*rh3h3
        b(in)  = b(in) + 2.*y(nh3)*y(nh3)*rh3h3
        b(ia)  = b(ia) + y(nh3)*y(nh3)*rh3h3
c..
c..c12(nu nup h3 p 2a) reaction
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) + rc12h3p
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - rc12h3p
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - rc12h3p
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - 2.*rc12h3p
        b1      = y(nc12)*rc12h3p
        b(nc12) = b(nc12) - b1
        b(nh3)  = b(nh3) + b1
        b(ip)   = b(ip) + b1
        b(ia)   = b(ia) + 2.*b1
       end if
      end if
c..
c..he3 plus he3 reaction
      if (nhe3 .ne. 0) then
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + 4.*y(nhe3)*rhe3he3
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 4.*y(nhe3)*rhe3he3
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - 2.*y(nhe3)*rhe3he3
       b(nhe3) = b(nhe3)- 2.*y(nhe3)*y(nhe3)*rhe3he3
       b(ip)   = b(ip) + 2.*y(nhe3)*y(nhe3)*rhe3he3
       b(ia)   = b(ia) + y(nhe3)*y(nhe3)*rhe3he3
c..
c..he3(n,g) and inverse (including neutrino irradiation)
       a1      = sig(1,nhe3)*aan
       a2      = sig(2,nhe3)
       a3      = sig(1,nhe3)*y(nhe3)
       b1      = -aan*sig(1,nhe3)*y(nhe3) + sig(2,nhe3)*aaa
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - a2
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a3
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a2
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - a1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - a3
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a3
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + a1
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - a2
       b(nhe3) = b(nhe3) + b1
       b(in)   = b(in) + b1
       b(ia)   = b(ia) - b1
c..
c..he3(d,p)he4
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nhe3)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nh2)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nh2)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) + y(nhe3)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nhe3)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nh2)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nhe3)*rhe3dp
       nt      = nt + 1
       iat     = eloc(nt)
       a(iat)  = a(iat) - y(nh2)*rhe3dp
       b(nh2)  = b(nh2) - y(nh2)*y(nhe3)*rhe3dp
       b(nhe3) = b(nhe3) - y(nh2)*y(nhe3)*rhe3dp
       b(ia)   = b(ia) + y(nh2)*y(nhe3)*rhe3dp
       b(ip)   = b(ip) + y(nh2)*y(nhe3)*rhe3dp
c..
c..he3(t,d)4he
       if (nh3 .ne. 0 .and. nh2 .ne. 0) then
        a1     = y(nhe3)*rhe3td
        a2     = y(nh3)*rhe3td
        b1     = y(nh3)*y(nhe3)*rhe3td
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) + a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) + a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) + a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) + a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) - a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) - a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) - a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat) = a(iat) - a2
        b(nhe3)= b(nhe3) - b1
        b(nh3) = b(nh3) - b1
        b(nh2) = b(nh2) + b1
        b(ia)  = b(ia) + b1
c..
c..he3(t,np)4he
        a1      = y(nhe3)*rhe3tnp
        a2      = y(nh3)*rhe3tnp
        b1      = y(nh3)*y(nhe3)*rhe3tnp
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) + a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) + a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) + a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) + a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - a2
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - a1
        nt      = nt + 1
        iat     = eloc(nt)
        a(iat)  = a(iat) - a2
        b(nhe3) = b(nhe3) - b1
        b(nh3)  = b(nh3) - b1
        b(ip)   = b(ip) + b1
        b(in)   = b(in) + b1
        b(ia)   = b(ia) + b1
       end if
      end if
c..
c..bullet check the counting
      if (nt .ne. nterm) then
       write(6,*) 'counted ',nt,' terms'
       write(6,*) 'should be ',nterm,' terms'
       stop 'fatal counting error in routine compt'
      end if
      return
      end
c..
c..
c..
      subroutine dtnuc(nucleu)
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'
c..  
c..   modified 9-19-91 to not take dth .gt. 0.3*taud for hd freeze outs
c..   this subroutine computes the current maximum allowed nuclear
c..   time step for a given zone by examining the abundance
c..   changes that occured in the previous time step.
c..  
c..   delchi is the maximum fractional change allowed in an
c..   abundance whose "y" value exceeds chimin.
c..  
c
      itot   =  i2+3
      tau    = fdtn * dth
      nucleu = 0
      do 50  n=1,itot
      if (y(n).lt.(0.1*chimin)) go to 50
      fak=b(n)
      if (fak)  20,50,20
   20 taug=abs((y(n)+chimin)/fak)*delchi*dth
      if (taug) 50,50,30
   30 if (tau-taug) 50,50,40
   40 tau=taug
      nucleu=n
   50 continue
      dth=tau

      return
      end
c..
c..
c..
      subroutine bwrit(nucleu,tim,iprint)
      implicit real*8 (a-h,o-z)
      save
      include 'burnf.dek'
c..  
c..   this subroutine writes out a complete edit of flows and
c..   abundances for a given zone and time step.
c..  
c..  
      dimension  flow(13,nburn),xa(nburn),x(nburn)
      dimension floweq(13*nburn)
      equivalence (flow,floweq)
c..  iprint indicates which unit number bwrit prints to, 6 or 20
      integer iprint
c..

c.. stable nuclei H-Pu238

      character*2 isymb(0:98)

      data isymb/
     &     ' n',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     1     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc',
     2     'Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
     3     'As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc',
     4     'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe',
     5     'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb',
     6     'Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os',
     7     'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr',
     8     'Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf'/

c
c..1=n,g  2=g,n  3=p,n  4=n,p  5=p,g  6=g,p  7=a,p  8=p,a  9=a,n
c..10=n,a  11=a,g  12=g,a  13=beta-positron
c
      i213=13*nburn
      call qvset (0.d0,floweq,i213)
      in=i2+1
      ip=i2+2
      ia=i2+3
c
c..  apply screening factors to general flows here, same as in compt
      do 140  i=1,i2
      j=nrr(1,i)
      if (j) 40,40,30
   30 fact=sig(1,i)*(y(i)*aan-b(i)*b(in))
      flow(2,i)=fact
      flow(1,i)=fact-sig(2,i)*y(j)
   40 j=nrr(2,i)
      if (j) 60,60,50
   50 fact=sig(3,i)*scfacp(i)*(y(i)*aap-b(i)*b(ip))
      flow(4,i)=fact
      flow(3,i)=fact-sig(4,i)*scfacp(i)*(y(j)*aan-b(j)*b(in))
      flow(13,i)=y(i)*sig(5,i)-y(j)*sig(6,i)
   60 j=nrr(3,i)
      if (j) 80,80,70
   70 fact=sig(7,i)*scfacp(i)*(y(i)*aap-b(i)*b(ip))
      flow(6,i)=fact
      flow(5,i)=fact-sig(8,i)*y(j)
   80 j=nrr(4,i)
      if (j) 100,100,90
   90 fact=sig(9,i)*scfaca(i)*(y(i)*aaa-b(i)*b(ia))
      flow(8,i)=fact
      flow(7,i)=fact-sig(10,i)*scfaca(i)*(y(j)*aap-b(j)*b(ip))
  100 j=nrr(5,i)
      if (j) 120,120,110
  110 fact=sig(11,i)*scfaca(i)*(y(i)*aaa-b(i)*b(ia))
      flow(10,i)=fact
      flow(9,i)=fact-sig(12,i)*scfaca(i)*(y(j)*aan-b(j)*b(in))
  120 j=nrr(6,i)
      if (j) 140,140,130
  130 fact=sig(13,i)*scfaca(i)*(y(i)*aaa-b(i)*b(ia))
      flow(12,i)=fact
      flow(11,i)=fact-sig(14,i)*y(j)
  140 continue
c..  rezro this matrix element to keep from double counting on next pas
      sig(6,nb8) = 0.
c..  assign flows for h3(pg) and he3(ng) not treated in general
      if (nh3.ne.0) then
       flow(6,nh3)=ph3pg
       flow(5,nh3)=ph3pg-phe4gp
      end if
      if (nhe3.ne.0) then
       flow(2,nhe3)=phe3ng
       flow(1,nhe3)=phe3ng-phe4gn
      end if
c
c..   figure flows for anomolous reactions
c     bullet proof pointers
      pra3=aaa*ra3*(aaa**2)
      if (nc12.ne.0) then
       pral=ral*y(nc12)
       pr24=y(nc12)*y(nc12)*r24
       if (no16.ne.0) then
        prc28=y(nc12)*y(no16)*rc28
       end if
      end if
      if (no16.ne.0) pr32=y(no16)*y(no16)*r32
      if (ne20.ne.0) pr20=y(ne20)*r20*aaa
      if (na23.ne.0) pr23=y(na23)*r23*aap
      if (mg23.ne.0) pr23n=y(mg23)*r23n*aan
      if (mg24.ne.0) prc24a=y(mg24)*aaa*rc24a
      if (nal27.ne.0) prc27p=y(nal27)*aap*rc27p
      if (nsi27.ne.0) prc27n=y(nsi27)*aan*rc27n
      if (nsi28.ne.0) pr28=y(nsi28)*r28*aaa
      if (np31.ne.0) pr31=y(np31)*r31*aap
      if (ns31.ne.0) pr31n=y(ns31)*r31n*aan
      pecap=aap*rpen
      rnecap=aan*rnep
c
c..   figure flows for new anomolous reactions
c     bullet proof pointers
      prpp=rpp*aap*aap
      ppng=rpng*aap*aan
      if (nh2.ne.0) then
       ppgn=rpgn*y(nh2)
       if (nh3.ne.0) ph2h3n=rh2h3n*y(nh2)*y(nh3)
      end if
      if (nh3.ne.0) then
       ph3pg=sig(7,nh3)*y(nh3)*aap
       phe4gp=sig(8,nh3)*aaa
       ph3h3=rh3h3*y(nh3)*y(nh3)
      end if
      if (nhe3.ne.0) then
       phe3ng=sig(1,nhe3)*y(nhe3)*aan
       phe4gn=sig(2,nhe3)*aaa
       phe3he3=rhe3he3*y(nhe3)*y(nhe3)
       if (nh2.ne.0) phe3dp=rhe3dp*y(nhe3)*y(nh2)
       if (nh3.ne.0) then 
        phe3td=y(nhe3)*y(nh3)*rhe3td
        phe3tnp=y(nhe3)*y(nh3)*rhe3tnp
       end if
      end if
      if (nb11.ne.0) pb11pa=rb11pa*y(nb11)*aap
      if (nli7.ne.0) then
       pli7pag=rli7pag*y(nli7)*aap
       if (nbe7.ne.0) pbe7eng=sig(6,nli7)*y(nbe7)
      end if
      if (nbe9.ne.0) then 
       prgaan=y(nbe9)*rgaan
       pbe9pd=y(nbe9)*aap*rbe9pd
      end if
      if (nb8.ne.0) pb8epn=rb8epn*y(nb8)
      if (nc11.ne.0) pc11na=rc11na*y(nc11)*aan
      if (nc12.ne.0) then
       pc12pa=y(nc12)*rc12pa
       pc12np=y(nc12)*rc12np
       pc12he3=y(nc12)*rc12he3
      end if
      praan=raan*aaa*aaa*aan
c
      nuclz=0
      nucln=0
      rmu=1./ye
      if (nucleu.gt.0) then
       nuclz=nz(nucleu)
       nucln=nn(nucleu)
      else
       nucleu=0
      end if
      rtim=etim-btime
      write(iprint,161) ncycle,tim
  161 format (1x,'cycle number =',i6,3x,'time =',1pe12.4)
c..  
      write(iprint,179) t9,rho,eb
  179 format (2x,'t9 =',f9.6,16x,'rho = ',1pe12.6,6x,'eb =',e11.3)
      write(iprint,170) nucleu,nuclz,nucln,dth
  170 format (2x,'nucleus changes fastest: ',i4,' z=',i3,
     1    3h n=,i3,12x,5hdth =,1pe11.4,10x,2e10.2)
c..  
c..   sum mass fractions
c..  
      call qvflot (xa,na,i2)
      call qvdot (xnorm,xa,y,i2)
      xnorm=xnorm+aan+aap+4.0*aaa
      xnorm=xnorm-1.

c..  
c..   write out mass non-conservation and neutron excess
c..  
      write(iprint,280)  xnorm,etan
  280 format (2x,'mass non-conservation =',1pe14.5,16x,'eta =',e14.7)
c..  
c..   write neutrino losses and sum of weak flows from fuller rates.
c..  
      call sneutx (t9,rho,rmu,qq)
      write(iprint,290) qq,snuw,eectot,eedtot
  290 format(2x,'plasma nu loss =',1pe9.2,4x,'beta nu loss =',e9.2,
     1     ' ecap+epos =',e10.3,10x,'e- decay =',e10.3)
      write(iprint,295) wrate,rectot,redtot
  295 format (2x,'weak rate =',1pe10.3,8x,'ecap =',e10.3,8x,
     1     'beta =',e10.3,//)
c..  
c..   write out mass fractions
c..  
      call qvmpy0 (x,xa,y,i2)
      aa=aaa*4.
      write(iprint,300) aap,aan,aa,(isymb(nz(i)),na(i),x(i),i=1,5)
  300 format   (2x,'  p',  1pe10.3,'    n',e10.3,'    a',e10.3,5(1x,a2,
     1i2,     e10.3))
      write(iprint,310)  (isymb(nz(i)),na(i),x(i),i=6,i2)
  310 format (7(1x,a2,i3, 1pe10.3))
      write(iprint,310)  (isymb(nz(i)),na(i),g(i),i=6,i2)

c..   write out special flows due to heavy ion reactions
c..   do new special light nuclei flows first

      write(iprint,180) prpp,ppng,ppgn
  180 format(/,2x,'pp',1pe11.3,2x,'fpng',e11.3,2x,'fpgn',e11.3)
      write(iprint,181)
  181 format(/,2x,'h2h3n',3x,2x,'h3h32n',2x,3x,'h3pg',3x,2x,
     1         'he3dp',3x,2x,'he3ng',3x,2x,'he3td',3x,2x,'he3tnp',2x,
     2            2x,'he3he3',2x)
      rvno=0.
      write(iprint,182) ph2h3n,ph3h3,ph3pg,phe3dp,phe3ng,phe3td,
     1  phe3tnp,phe3he3
      write(iprint,182) rvno,rvno,phe4gp,rvno,phe4gn,rvno,rvno,rvno
  182 format(8(1pe10.2))
      write(iprint,183)
  183 format(/,2x,'li7pag',2x,2x,'be7eng',2x,2x,'b8epn',3x,2x,'be9pd',
     1         3x,2x,'b11pa',3x,2x,'c11na',3x,3x,'raan',3x,3x,'rgaan')
      write(iprint,184) pli7pag,pbe7eng,pb8epn,pbe9pd,pb11pa,pc11na,
     1 praan,prgaan
  184 format(8(1pe10.2))
      write(iprint,185) pc12pa,pc12np,pc12he3
  185 format(/,2x,'c12(nu,nup,pa)=',1pe10.2,2x,'c12(nu,nup,pn)=',
     1       e10.2,2x,'c12(nu,nup,he3)=',e10.3)

c.. continue with flows due to heavy ion reactions and their decay chann
      write(iprint,190)
  190 format (/,5x,'3a',5x,'c12photo',3x,'c12+c12',3x,'c12+o16',3x,
     1     'o16+o16',3x,'na23+p',4x,'ne20+a',4x,'mg23+n',4x,'p31+p',
     2     5x,'si28+a',4x,'s31+n')
      write(iprint,200) pra3,pral,pr24,prc28,pr32,pr23,pr20,pr23n,
     1     pr31,pr28,pr31n
  200 format (11(1pe10.2))
      write(iprint,210) sc3a,sc1212,sc1216,sc1616,b24p,b24a,b24n,
     1     b32p,b32a,b32n
  210 format (1pe10.2,10x,9e10.2)
      write(iprint,211)
  211 format (/,53x,'mg24+a',4x,'al27+p',4x,'si27+n',14x,
     1     'p(e-,nu)n',2x,'n(e+,nub)p')
      write(iprint,212) prc24a,prc27p,prc27n,pecap,rnecap
  212 format (50x,3(1pe10.2),11x,e10.2,e10.2)
      write(iprint,212) bc24a,bc27p,bc27n,spen,snep

c..   write out detailed flows in pairs (net flow and forward flow)

      write(iprint,220)
  220 format (/,' z/a ',8x,'ng',16x,'pn',16x,'pg',16x,'ap',
     1     16x,'an',16x,'ag',10x,'weak',/)
      write(iprint,230)  (isymb(nz(i)),na(i),(flow(j,i),j=1,13),i=1,i2)
  230 format (1x,a2,i3, 1pe9.2,11e9.2,e8.1)
      write(iprint,230)  (isymb(nz(i)),na(i),(sig(j,i),j=1,13),i=1,i2)
c
      return
      end
c..

c..
      subroutine builda(inu)
      implicit double precision (a-h,o-z)
      save
      include 'burnf.dek'
c..
c..this routine tags the nonzero matrix element locations.
c..
c..note: if routine compt is changed, then this routine
c..needs to be changed. the manner in which to change it should be clear.
c..
c..
c..declare
      integer          inu,in,ip,ia,j,k
c..
c..set the n, p, alpha pointers
      in = i2 + 1
      ip = i2 + 2
      ia = i2 + 3
c..
c..from routine compt:
      do 11 j=1,i2
c..
c..(n,gam) reactions
       k = nrr(1,j)
       if (k .ne. 0) then
        call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,k,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
c..
c..(p,n) and beta-, beta+ reactions
       k = nrr(2,j)
       if (k .ne. 0) then
        call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
c..
c..(p,gam) reactions
       k = nrr(3,j)
       if (k .ne. 0) then
        call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,k,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
c..
c..(alp,p) reactions
       k = nrr(4,j)
       if (k .ne. 0) then
        call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
c..
c..(alp,n) reaction
       k = nrr(5,j)
       if (k .ne. 0) then
        call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
c..
c..(a,gam) reactions
       k = nrr(6,j)
       if (k .ne. 0) then
        call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(j,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,k,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(k,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,k,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
11    continue
c..
c..branch around this if we are not doing neutrinos
      if (inu .ne. 0) then
       do 110 j=1,i2
c..
c..(nu, e-) recations
        k = nrrneut(1,j)
        if (k .ne. 0) then
         call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(in,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        end if
c..
c..(nu, e-) reactions
        k = nrrneut(3,j)
        if (k .ne. 0) then
         call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(ip,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        end if
c..
c..(nu, e+) reactions
        k = nrrneut(4,j)
        if (k .ne. 0) then
         call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(in,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        end if
c..
c..(nu, e+) reactions
        k = nrrneut(6,j)
        if (k .ne. 0) then
         call locat(j,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(k,j,irow,icol,nzo,naij,eloc,nterm,neloc)
         call locat(ip,j,irow,icol,nzo,naij,eloc,nterm,neloc)
        end if
110    continue
      end if
c..
c..
c..from routine spcomp:
c..
c..p(e-,nu)n and n(e+,nub)p reaction
      call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
      call locat(ip,in,irow,icol,nzo,naij,eloc,nterm,neloc)
      call locat(in,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
      call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..triple alpha
      if (nc12 .ne. 0) then
       call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..c12+c12
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,ne20,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,na23,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ne20,ne20,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ne20,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ne20,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ne20,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(na23,na23,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(na23,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(na23,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,na23,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,mg23,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg23,mg23,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg23,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg23,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,mg23,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..o16 burning
      if (no16 .ne. 0) then
       call locat(no16,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,nsi28,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,np31,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,ns31,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi28,nsi28,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nsi28,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi28,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi28,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(np31,np31,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,np31,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(np31,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(np31,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(np30,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ns31,ns31,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,ns31,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ns31,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ns31,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..c12+o16 reaction
      if (nc12 .ne. 0 .and. no16 .ne. 0) then
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,nal27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,nsi27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,mg24,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc12,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,nal27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,nsi27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,mg24,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(no16,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg24,mg24,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg24,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg24,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(mg24,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nal27,nal27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nal27,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nal27,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nal27,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi27,nsi27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi27,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi27,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nsi27,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,mg24,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nal27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nsi27,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,no16,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..7li(pa) reaction
      if (nli7 .ne. 0) then 
       call locat(nli7,nli7,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nli7,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nli7,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nli7,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..11b(pa) reaction
      if (nb11 .ne. 0) then
       call locat(nb11,nb11,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nb11,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nb11,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nb11,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c.. 8b(e+ nu) reaction
      if (nb8 .ne. 0) then
       call locat(nb8,nb8,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nb8,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..12c(nu,nup n p) reaction
      if (nb10 .ne. 0) then
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nb10,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..12c(nu,nup p a) reaction
      if (nli7 .ne. 0) then
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nli7,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..12c(nu,nup n a) reaction
      if (nbe7 .ne. 0) then
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nbe7,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..c12(nu nup he3 n 2a) reaction
      if (nhe3 .ne. 0) then
       call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nhe3,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..12c(nu,nup he3)9be reaction
       if (nbe9 .ne. 0) then
        call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nbe9,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nhe3,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
      end if
c..
c..9be(p,d)2 alpha reaction
      if (nbe9 .ne. 0 .and. nh2 .ne. 0) then
       call locat(nbe9,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nbe9,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..11c(na) reaction
      if (nc11 .ne. 0) then
       call locat(nc11,nc11,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nc11,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nc11,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nc11,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,in,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..a(an,g)9be reaction
      if (nbe9 .ne. 0) then
       call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nbe9,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nbe9,nbe9,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nbe9,in,irow,icol,nzo,naij,eloc,nterm,neloc)
      end if
c..
c..proton-proton reaction
      if (nh2. ne. 0) then
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..neutron capture on the proton reaction
       call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..the "dt" reaction  2h(3h,n)4he
       if (nh3 .ne. 0) then 
        call locat(nh2,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh2,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..h3(p,g)he4   and inverse including neutrino irradiation
        call locat(nh3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,ip,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..t(t,2n)he4; no inverse included
        call locat(nh3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..c12(nu nup h3 p 2a)
        call locat(nc12,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nc12,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
      end if
c..
c..he3 plus he3 reaction
      if (nhe3 .ne. 0) then 
       call locat(nhe3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..he3(n,g) and inverse (including neutrino irradiation)
       call locat(nhe3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nhe3,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nhe3,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,in,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(in,ia,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..he3(d,p)he4
       call locat(nh2,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nh2,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nhe3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(nhe3,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ia,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nh2,irow,icol,nzo,naij,eloc,nterm,neloc)
       call locat(ip,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..he3(t,d)4he
       if (nh3 .ne. 0 .and. nh2 .ne. 0) then 
        call locat(nhe3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nhe3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh2,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh2,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
c..
c..he3(t,np)4he
        call locat(nhe3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nhe3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(nh3,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ia,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(ip,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,nh3,irow,icol,nzo,naij,eloc,nterm,neloc)
        call locat(in,nhe3,irow,icol,nzo,naij,eloc,nterm,neloc)
       end if
      end if
      return
      end
c..
c..
c..
c..
c..
      subroutine locat(i,j,irow,icol,nzo,np,eloc,nterm,neloc)
      implicit double precision (a-h,o-z)
      save
c..
c..this routine check that each location is tagged only once.
c..
c..input is the matrix location i and j, and the arrays that store the
c..locations iloc andd jloc, both of physical dimension np.
c..
c..output is the current number of nonzeros and the filled iloc and
c..jloc arrays.
c..
c..declare
      integer          i,j,np,irow(np),icol(np),nzo,
     1                 neloc,eloc(neloc),nterm,
     2                 n,ifirst
      data             ifirst/0/
c..
c..initialize nzo on the first pass through
      if (ifirst .eq. 0) then
       ifirst = 1
       nzo    = 0
       nterm  = 0
      end if
c..
c..increment the number of terms
      nterm = nterm + 1
      if (nterm .gt. neloc) then
       write(6,*) 'nterm > neloc in routine locat',nterm,neloc
       call exit(1)
      end if
c..
c..see if we have this location already
      do 10 n=1,nzo
       if (irow(n) .eq. i  .and.  icol(n) .eq. j) then
        eloc(nterm) = n
        return
       end if
10    continue
c..
c..we dont have this location; add it in 
      nzo         = nzo + 1
      if (nzo .gt. np) stop 'nzo > np in routine locat'
      irow(nzo)   = i
      icol(nzo)   = j
      eloc(nterm) = nzo
      return
      end
c..
c..
c..
c..
c..
c..the rest of this file contains the harwell ma28 sparse matrix routines
c..
c..easy to use front end routines: 
c..routine ma28ad does the symbolic and numeric lu decomp 
c..routine ma28bd does the numeric lu decomp of ma28ad
c..routine ma28cd solves the system of equations directly
c..routine ma28id solves the system iteratively
c..routine ma28dd does some pointer work
c..
c..these are the hardball routines:
c..routine ma30ad does core symbolic and numeric lu decomp 
c..routine ma30bd does the numeric decomp the sparse pattern
c..routine ma30cd solves the linear system
c..routine ma30dd does garbage collection
c..
c..support hardball routines:
c..routine ma28int1 does some common block initialization
c..roytine ma28int2 does some common block initialization
c..routine ma28int3 does some common block initialization
c..routine mc20ad sort a matrix
c..routine mc23ad does the block triangularization pointers
c..routine mc22ad reorders the off diagonal blocks based on the pivot info
c..routine mc21a front end of mc21b
c..routine mc21b pernutes the rows to get a zero free diagonal
c..routine mc13d front end for mc13e
c..routine mc13e permutes a lower triangular block
c..routine mc24ad gets a growth rate of fillin
c..
c..
c..
c..
      subroutine ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag)
      implicit double precision (a-h,o-z)
      save
c..
c..this subroutine performs the lu factorization of a.
c..
c..input:
c..n     order of matrix ..  not altered by subroutine
c..nz    number of non-zeros in input matrix ..  not altered by subroutine
c..a     is a real array  length licn.  holds non-zeros of matrix on entry
c..      and non-zeros of factors on exit.  reordered by mc20a/ad and
c..      mc23a/ad and altered by ma30a/ad
c..licn  integer  length of arrays a and icn ..  not altered by subroutine
c..irn   integer array of length lirn.  holds row indices on input.
c..      used as workspace by ma30a/ad to hold column orientation of matrix
c..lirn  integer  length of array irn ..  not altered by the subroutine
c..icn   integer array of length licn.  holds column indices on entry
c..      and column indices of decomposed matrix on exit. reordered by
c..      mc20a/ad and mc23a/ad and altered by ma30a/ad.
c..u     real variable  set by user to control bias towards numeric or
c..      sparsity pivoting.  u=1.0 gives partial pivoting while u=0. does
c..      not check multipliers at all.  values of u greater than one are
c..      treated as one while negative values are treated as zero.  not
c..      altered by subroutine.
c..ikeep integer array of length 5*n  used as workspace by ma28a/ad
c..      it is not required to be set on entry and, on exit, it contains 
c..      information about the decomposition. it should be preserved between 
c..      this call and subsequent calls to ma28b/bd or ma28c/cd.
c..      ikeep(i,1),i=1,n  holds the total length of the part of row i
c..      in the diagonal block.
c..      row ikeep(i,2),i=1,n  of the input matrix is the ith row in
c..      pivot order.
c..      column ikeep(i,3),i=1,n  of the input matrix is the ith column
c..      in pivot order.
c..      ikeep(i,4),i=1,n  holds the length of the part of row i in
c..      the l part of the l/u decomposition.
c..      ikeep(i,5),i=1,n  holds the length of the part of row i in the
c..      off-diagonal blocks.  if there is only one diagonal block,
c..      ikeep(1,5) will be set to -1.
c..iw    integer array of length 8*n.  if the option nsrch.le.n is
c..      used, then the length of array iw can be reduced to 7*n.
c..w     real array  length n.  used by mc24a/ad both as workspace and to
c..      return growth estimate in w(1).  the use of this array by ma28a/ad
c..      is thus optional depending on common block logical variable grow.
c..iflag integer variable  used as error flag by routine.  a positive
c..      or zero value on exit indicates success.  possible negative
c..      values are -1 through -14.
c..
c..declare
      integer          n,nz,licn,lirn,iflag,irn(lirn),icn(licn),
     1                 ikeep(n,5),iw(n,8),i,j1,j2,jj,j,length,
     2                 newpos,move,newj1,jay,knum,ii,i1,iend
      double precision a(licn),u,w(n)
c..
c..common and private variables. common block ma28f/fd is used merely
c..to communicate with common block ma30f/fd  so that the user
c..need not declare this common block in his main program.
c..
c..the common block variables are:
c..lp,mp    default value 6 (line printer).  unit number
c..         for error messages and duplicate element warning resp.
c..nlp,mlp  unit number for messages from ma30a/ad and
c..         mc23a/ad resp.  set by ma28a/ad to value of lp.
c..lblock   logical  default value true.  if true mc23a/ad is used
c..         to first permute the matrix to block lower triangular form.
c..grow     logical  default value true.  if true then an estimate
c..         of the increase in size of matrix elements during l/u
c..         decomposition is given by mc24a/ad.
c..eps,rmin,resid  variables not referenced by ma28a/ad.
c..irncp,icncp  set to number of compresses on arrays irn and icn/a 
c..minirn,minicn  minimum length of arrays irn and icn/a; for success on 
c..               future runs.
c..irank  integer   estimated rank of matrix.
c..mirncp,micncp,mirank,mirn,micn communicate between ma30f/fd and ma28f/fd 
c..                               values of above named variables with 
c..                               somewhat similar names.
c..abort1,abort2  logical variables with default value true.  if false
c..               then decomposition will be performed even if the matrix is
c..               structurally or numerically singular respectively.
c..aborta,abortb  logical variables used to communicate values of
c..               abort1 and abort2 to ma30a/ad.
c..abort  logical  used to communicate value of abort1 to mc23a/ad.
c..abort3  logical variable not referenced by ma28a/ad.
c..idisp  integer array  length 2.  used to communicate information
c..       on decomposition between this call to ma28a/ad and subsequent
c..       calls to ma28b/bd and ma28c/cd.  on exit, idisp(1) and
c..       idisp(2) indicate position in arrays a and icn of the
c..       first and last elements in the l/u decomposition of the
c..       diagonal blocks, respectively.
c..numnz  integer  structural rank of matrix.
c..num    integer  number of diagonal blocks.
c..large  integer  size of largest diagonal block.
c..
c..
      logical          grow,lblock,abort,abort1,abort2,abort3,aborta,
     1                 abortb,lbig,lbig1
      integer          idisp(2),lp,mp,irncp,icncp,minirn,minicn,
     1                 irank,ndrop,maxit,noiter,nsrch,istart,
     2                 ndrop1,nsrch1,nlp,mirncp,micncp,mirank,
     3                 mirn,micn,mlp,numnz,num,large,lpiv(10),
     4                 lnpiv(10),mapiv,manpiv,iavpiv,ianpiv,kountl,
     5                 ifirst
      double precision tol,themax,big,dxmax,errmax,dres,cgce,
     1                 tol1,big1,upriv,rmin,eps,resid,zero
c..
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn,
     1                 irank,abort1,abort2
      common /ma28gd/  idisp
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce,
     1                 ndrop,maxit,noiter,nsrch,istart,lbig
      common /ma30id/  tol1,big1,ndrop1,nsrch1,lbig1
      common /ma30ed/  nlp,aborta,abortb,abort3
      common /ma30fd/  mirncp,micncp,mirank,mirn,micn
      common /mc23bd/  mlp,numnz,num,large,abort
      common /lpivot/  lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
      data zero        /0.0d0/, ifirst/0/
c..
c..format statements
99990 format(1x,'error return from ma28a/ad because')
99991 format(1x,'error return from ma30a/ad')
99992 format(1x,'error return from mc23a/ad')
99993 format(1x,'duplicate element in position',i8,' ',i8,
     1          'with value ',1pe22.14)
99994 format (1x,i6,'element with value',1pe22.14,'is out of range',/,
     1        1x,'with indices',i8,' ',i8)
99995 format(1x,'error return from ma28a/ad; indices out of range')
99996 format(1x,'lirn too small = ',i10)
99997 format(1x,'licn too small = ',i10)
99998 format(1x,'nz non positive = ',i10)
99999 format(1x,'n out of range = ',i10)
c..
c..

c..initialization and transfer of information between common blocks
      if (ifirst .eq. 0) then
       ifirst = 1
       call ma28int1
       call ma28int2
       call ma28int3
      end if
      iflag  = 0
      aborta = abort1
      abortb = abort2
      abort  = abort1
      mlp    = lp
      nlp    = lp
      tol1   = tol
      lbig1  = lbig
      nsrch1 = nsrch
c..
c..upriv private copy of u is used in case it is outside
      upriv = u
c..
c..simple data check on input variables and array dimensions.
      if (n.gt.0) goto 10
      iflag = -8
      if (lp.ne.0) write (lp,99999) n
      goto 210
10    if (nz.gt.0) goto 20
      iflag = -9
      if (lp.ne.0) write (lp,99998) nz
      goto 210
20    if (licn.ge.nz) goto 30
      iflag = -10
      if (lp.ne.0) write (lp,99997) licn
      goto 210
30    if (lirn.ge.nz) goto 40
      iflag = -11
      if (lp.ne.0) write (lp,99996) lirn
      goto 210
c..
c..data check to see if all indices lie between 1 and n.
40    do 50 i=1,nz
       if (irn(i).gt.0 .and. irn(i).le.n .and. icn(i).gt.0 .and.
     1     icn(i).le.n) goto 50
       if (iflag.eq.0 .and. lp.ne.0) write (lp,99995)
       iflag = -12
       if (lp.ne.0) write (lp,99994) i,a(i),irn(i),icn(i)
50    continue
      if (iflag.lt.0) goto 220
c..
c..sort matrix into row order.
      call mc20ad(n,nz,a,icn,iw,irn,0)
c..
c..part of ikeep is used here as a work-array.  ikeep(i,2) is the last row to 
c..have a non-zero in column i.  ikeep(i,3) is the off-set of column i from 
c..the start of the row.
      do 60 i=1,n
       ikeep(i,2) = 0
       ikeep(i,1) = 0
60    continue
c..
c..check for duplicate elements .. summing any such entries and printing a 
c..warning message on unit mp. move is equal to the number of duplicate 
c..elements found; largest element in the matrix is themax; j1 is position in 
c..arrays of first non-zero in row.
      move   = 0
      themax = zero
      j1     = iw(1,1)
      do 130 i=1,n
       iend = nz + 1
       if (i.ne.n) iend = iw(i+1,1)
       length = iend - j1
       if (length.eq.0) goto 130
       j2 = iend - 1
       newj1 = j1 - move
       do 120 jj=j1,j2
        j = icn(jj)
        themax = max(themax,abs(a(jj)))
        if (ikeep(j,2).eq.i) goto 110
c..
c..first time column has ocurred in current row.
        ikeep(j,2) = i
        ikeep(j,3) = jj - move - newj1
        if (move.eq.0) goto 120
c..
c..shift necessary because of previous duplicate element.
        newpos = jj - move
        a(newpos) = a(jj)
        icn(newpos) = icn(jj)
        goto 120
c..
c..duplicate element.
110     move = move + 1
        length = length - 1
        jay = ikeep(j,3) + newj1
        if (mp.ne.0) write (mp,99993) i,j,a(jj)
        a(jay) = a(jay) + a(jj)
        themax = max(themax,abs(a(jay)))
120    continue
       ikeep(i,1) = length
       j1 = iend
130    continue
c..
c..knum is actual number of non-zeros in matrix with any multiple entries 
c..counted only once
      knum = nz - move
      if (.not.lblock) goto 140
c..
c..perform block triangularisation.
      call mc23ad(n,icn,a,licn,ikeep,idisp,ikeep(1,2),
     1            ikeep(1,3),ikeep(1,5),iw(1,3),iw)
      if (idisp(1).gt.0) goto 170
      iflag = -7
      if (idisp(1).eq.-1) iflag = -1
      if (lp.ne.0) write (lp,99992)
      goto 210
c..
c..block triangularization not requested. move structure to end of data arrays 
c..in preparation for ma30a/ad; set lenoff(1) to -1 and set permutation arrays.
140   do 150 i=1,knum
       ii = knum - i + 1
       newpos = licn - i + 1
       icn(newpos) = icn(ii)
       a(newpos) = a(ii)
150   continue
      idisp(1) = 1
      idisp(2) = licn - knum + 1
      do 160 i=1,n
       ikeep(i,2) = i
       ikeep(i,3) = i
160   continue
      ikeep(1,5) = -1
170   if (lbig) big1 = themax
      if (nsrch.le.n) goto 180
c..
c..perform l/u decomosition on diagonal blocks.
      call ma30ad(n,icn,a,licn,ikeep,ikeep(1,4),idisp,
     1           ikeep(1,2),ikeep(1,3),irn,lirn,iw(1,2),iw(1,3),iw(1,4),
     2           iw(1,5),iw(1,6),iw(1,7),iw(1,8),iw,upriv,iflag)
      goto 190
c..
c..this call if used if nsrch has been set less than or equal n; in this case, 
c..two integer work arrays of length can be saved.
180    call ma30ad(n,icn,a,licn,ikeep,ikeep(1,4),idisp,
     1           ikeep(1,2),ikeep(1,3),irn,lirn,iw(1,2),iw(1,3),iw(1,4),
     2           iw(1,5),iw,iw,iw(1,6),iw,upriv,iflag)
c..
c..transfer common block information.
190   minirn = max0(mirn,nz)
      minicn = max0(micn,nz)
      irncp = mirncp
      icncp = micncp
      irank = mirank
      ndrop = ndrop1
      if (lbig) big = big1
      if (iflag.ge.0) goto 200
      if (lp.ne.0) write (lp,99991)
      goto 210
c..
c..reorder off-diagonal blocks according to pivot permutation.
200   i1 = idisp(1) - 1
      if (i1.ne.0) call mc22ad(n,icn,a,i1,ikeep(1,5),ikeep(1,2),
     1                         ikeep(1,3),iw,irn)
      i1 = idisp(1)
      iend = licn - i1 + 1
c..
c..optionally calculate element growth estimate.
      if (grow) call mc24ad(n,icn,a(i1),iend,ikeep,ikeep(1,4),w)
c..
c..increment growth estimate by original maximum element.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n.gt.1) w(2) = themax
c..
c..set flag if the only error is due to duplicate elements.
      if (iflag.ge.0 .and. move.ne.0) iflag = -14
      goto 220
210   if (lp.ne.0) write (lp,99990)
220   return
      end
c..
c..
c..
c..
c..
      subroutine ma28bd(n,nz,a,licn,ivect,jvect,icn,ikeep,iw,w,iflag)
      implicit double precision (a-h,o-z)
      save
c..
c..this subroutine factorizes a matrix with the same pattern as that
c..previously factorized by ma28a/ad.
c..
c..input is :
c..n      order of matrix  not altered by subroutine.
c..nz     number of non-zeros in input matrix  not altered by subroutine.
c..a      array  length licn.  holds non-zeros of matrix on entry and 
c..       non-zeros of factors on exit.  reordered by ma28d/dd and altered by 
c..       subroutine ma30b/bd.
c..licn   integer  length of arrays a and icn.  not altered by subroutine.
c..ivect,jvect  integer arrays of length nz.  hold row and column
c..       indices of non-zeros respectively.  not altered by subroutine.
c..icn    integer array of length licn.  same array as output from ma28a/ad.  
c..       unchanged by ma28b/bd.
c..ikeep  integer array of length 5*n.  same array as output from
c..       ma28a/ad.  unchanged by ma28b/bd.
c..iw     integer array  length 5*n.  used as workspace by ma28d/dd and
c..       ma30b/bd.
c..w      array  length n.  used as workspace by ma28d/dd,ma30b/bd and 
c..       (optionally) mc24a/ad.
c..iflag  error flag with positive or zero value indicating success.
c..
c..
c..declare
      integer          n,nz,licn,iw(n,5),iflag,ikeep(n,5),ivect(nz),
     1                 jvect(nz),icn(licn),i1,iend,idup
      double precision a(licn),w(n)
c..
c..private and common variables: unless otherwise stated common block 
c..variables are as in ma28a/ad. those variables referenced by ma28b/bd are 
c..mentioned below.
c..
c..lp,mp  used as in ma28a/ad as unit number for error and
c..       warning messages, respectively.
c..nlp    variable used to give value of lp to ma30e/ed.
c..eps    real/double precision  ma30b/bd will output a positive value
c..       for iflag if any modulus of the ratio of pivot element to the
c..       largest element in its row (u part only) is less than eps (unless
c..       eps is greater than 1.0 when no action takes place).
c..rmin   variable equal to the value of this minimum ratio in cases where 
c..       eps is less than or equal to 1.0.
c..meps,mrmin variables used by the subroutine to communicate between common 
c..        blocks ma28f/fd and ma30g/gd.
c..
c..declare
      logical          grow,lblock,aborta,abortb,abort1,abort2,abort3,
     1                 lbig,lbig1
      integer          idisp(2),mp,lp,irncp,icncp,minirn,minicn,irank,
     1                 ndrop,maxit,noiter,nsrch,istart,nlp,ndrop1,nsrch1
      double precision eps,meps,rmin,mrmin,resid,tol,themax,big,dxmax,
     1                 errmax,dres,cgce,tol1,big1
      common /ma28ed/  mp,lp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn,irank,
     1                 abort1,abort2
      common /ma28gd/  idisp
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce,ndrop,
     1                 maxit,noiter,nsrch,istart,lbig
      common /ma30ed/  nlp,aborta,abortb,abort3
      common /ma30gd/  meps,mrmin
      common /ma30id/  tol1,big1,ndrop1,nsrch1,lbig1
c..
c..formats
99994 format(1x,'error return from ma28b/bd because')
99995 format(1x,'error return from ma30b/bd')
99996 format(1x,'licn too small = ',i10)
99997 format(1x,'nz non positive = ',i10)
99998 format(1x,'n out of range = ',i10)
99999 format(1x,'error return from ma28b/bd with iflag=',i4,/,
     1       1x,i7,' entries dropped from structure by ma28a/ad')
c..
c..
c..check to see if elements were dropped in previous ma28a/ad call.
      if (ndrop.eq.0) goto 10
      iflag = -15
      write (6,99999) iflag,ndrop
      goto 70
10    iflag = 0
      meps  = eps
      nlp   = lp
c..
c..simple data check on variables.
      if (n.gt.0) goto 20
      iflag = -11
      if (lp.ne.0) write (lp,99998) n
      goto 60
20    if (nz.gt.0) goto 30
      iflag = -10
      if (lp.ne.0) write (lp,99997) nz
      goto 60
30    if (licn.ge.nz) goto 40
      iflag = -9
      if (lp.ne.0) write (lp,99996) licn
      goto 60
c..
c..
40     call ma28dd(n,a,licn,ivect,jvect,nz,icn,ikeep,ikeep(1,4),
     1             ikeep(1,5),ikeep(1,2),ikeep(1,3),iw(1,3),iw,
     2             w(1),iflag)
c..
c..themax is largest element in matrix
      themax = w(1)
      if (lbig) big1 = themax
c..
c..idup equals one if there were duplicate elements, zero otherwise.
      idup = 0
      if (iflag.eq.(n+1)) idup = 1
      if (iflag.lt.0) goto 60
c..
c..perform row-gauss elimination on the structure received from ma28d/dd
      call ma30bd(n,icn,a,licn,ikeep,ikeep(1,4),idisp,
     1            ikeep(1,2),ikeep(1,3),w,iw,iflag)
c..
c..transfer common block information.
      if (lbig) big1 = big
      rmin = mrmin
      if (iflag.ge.0) goto 50
      iflag = -2
      if (lp.ne.0) write (lp,99995)
      goto 60
c..
c..optionally calculate the growth parameter.
50    i1   = idisp(1)
      iend = licn - i1 + 1
      if (grow) call mc24ad(n,icn,a(i1),iend,ikeep,ikeep(1,4),w)
c..
c..increment estimate by largest element in input matrix.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n.gt.1) w(2) = themax
c..
c..set flag if the only error is due to duplicate elements.
      if (idup.eq.1 .and. iflag.ge.0) iflag = -14
      goto 70
60    if (lp.ne.0) write (lp,99994)
70    return
      end
c..
c..
c..
c..
c..
      subroutine ma28cd(n,a,licn,icn,ikeep,rhs,w,mtype)
      implicit double precision (a-h,o-z)
      save
c..
c..uses the factors from ma28a/ad or ma28b/bd to solve a system of equations
c..
c..input:
c..n     order of matrix  not altered by subroutine.
c..a     array  length licn.  the same array as most recent call to ma28a/ad 
c..      or ma28b/bd.
c..licn  length of arrays a and icn.  not altered by subroutine.
c..icn   integer array of length licn.  same array as output from ma28a/ad.  
c..      unchanged by ma28c/cd.
c..ikeep integer array of length 5*n.  same array as output from ma28a/ad.  
c..      unchanged by ma28c/cd.
c..rhs   array  length n.  on entry, it holds the right hand side.  
c..      on exit, the solution vector.
c..w     array  length n. used as workspace by ma30c/cd.
c..mtype integer  used to tell ma30c/cd to solve the direct equation
c..      (mtype=1) or its transpose (mtype.ne.1).
c..
c..resid  variable returns maximum residual of equations where pivot was zero.
c..mresid variable used by ma28c/cd to communicate with ma28f/fd and ma30h/hd.
c..idisp  integer array ; the same as that used by ma28a/ad. un changed.
c..
c..declare
      logical          abort1,abort2
      integer          n,licn,idisp(2),icn(licn),ikeep(n,5),
     1                 irncp,icncp,minirn,minicn,irank,mtype
      double precision a(licn),rhs(n),w(n),resid,mresid,eps,rmin
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn,
     1                 irank,abort1,abort2
      common /ma28gd/  idisp
      common /ma30hd/  mresid
c..
c..this call performs the solution of the set of equations.
      call ma30cd(n,icn,a,licn,ikeep,ikeep(1,4),ikeep(1,5),
     1            idisp,ikeep(1,2),ikeep(1,3),rhs,w,mtype)
c..
c..transfer common block information.
      resid = mresid
      return
      end
c..
c..
c..
c..
c..
      subroutine ma28id(n,nz,aorg,irnorg,icnorg,licn,a,icn,
     1                  ikeep,rhs,x,r,w,mtype,prec,iflag)
      implicit double precision (a-h,o-z)
      save
c..
c..this subroutine uses the factors from an earlier call to ma28a/ad
c..or ma28b/bd to solve the system of equations with iterative refinement.
c..
c..parameters are:
c..
c..n    order of the matrix. it is not altered by the subroutine.
c..nz   number of entries in the original matrix.  not altered by subroutine.
c..     for this entry the original matrix must have been saved in
c..     aorg,irnorg,icnorg where entry aorg(k) is in row irnorg(k) and
c..     column icnorg(k), k=1,.. nz.  information about the factors of a
c..     is communicated to this subroutine via the parameters licn, a, icn
c..     and ikeep where:
c..aorg   array of length nz.  not altered by ma28i/id.
c..irnorg array of length nz.  not altered by ma28i/id.
c..icnorg array of length nz.  not altered by ma28i/id.
c..licn   equal to the length of arrays a and icn. not altered
c..a    array of length licn. it must be unchanged since the last call
c..     to ma28a/ad or ma28b/bd. it is not altered by the subroutine.
c..icn, ikeep are the arrays (of lengths licn and 5*n, respectively) of
c..     the same names as in the previous all to ma28a/ad. they should be
c..     unchanged since this earlier call. not altered.
c..
c..other parameters are as follows:
c..rhs array of length n. the user must set rhs(i) to contain the
c..    value of the i th component of the right hand side. not altered.
c..
c..x   array of length n. if an initial guess of the solution is
c..    given (istart equal to 1), then the user must set x(i) to contain
c..    the value of the i th component of the estimated solution.  on
c..    exit, x(i) contains the i th component of the solution vector.
c..r   array of length n. it need not be set on entry.  on exit, r(i)
c..    contains the i th component of an estimate of the error if maxit
c..    is greater than 0.
c..w is an array of length n. it is used as workspace by ma28i/id.
c..mtype must be set to determine whether ma28i/id will solve a*x=rhs
c..     (mtype equal to 1) or at*x=rhs (mtype ne 1, zero say). not altered.
c..prec should be set by the user to the relative accuracy required. the
c..     iterative refinement will terminate if the magnitude of the
c..     largest component of the estimated error relative to the largest
c..     component in the solution is less than prec. not altered.
c..iflag is a diagnostic flag which will be set to zero on successful
c..      exit from ma28i/id, otherwise it will have a non-zero value. the
c..      non-zero value iflag can have on exit from ma28i/id are .. 
c..      -16    indicating that more than maxit iteartions are required.
c..      -17    indicating that more convergence was too slow.
c..
c..declare
      integer          n,nz,licn,mtype,iflag,icnorg(nz),irnorg(nz),
     1                 ikeep(n,5),icn(licn),i,iterat,nrow,ncol
      double precision a(licn),aorg(nz),rhs(n),r(n),x(n),w(n),prec,
     1                 d,dd,conver,zero
c..
c..common block communication
      logical          lblock,grow,lbig
      integer          lp,mp,ndrop,maxit,noiter,nsrch,istart
      double precision tol,themax,big,dxmax,errmax,dres,cgce
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce,
     1                 ndrop,maxit,noiter,nsrch,istart,lbig
      data             zero /0.0d0/
c..
c..
c..formats
99998 format(1x,'error return from ma28i with iflag = ', i3,/,
     1       1x,'convergence rate of',1pe9.2,'too slow',/,
     2       1x,'maximum acceptable rate set to ',1pe9.2)
99999 format(1x,'error return from ma28i/id with iflag = ',i3,/,
     1       1x,'more than',i5,'iterations required')
c..
c..
c..initialization of noiter,errmax and iflag.
      noiter = 0
      errmax = zero
      iflag   = 0
c..
c..jump if a starting vector has been supplied by the user.
      if (istart.eq.1) goto 20
c..
c..make a copy of the right-hand side vector.
      do 10 i=1,n
       x(i) = rhs(i)
10    continue
c..
c..find the first solution.
      call ma28cd(n,a,licn,icn,ikeep,x,w,mtype)
c..
c..stop the computations if   maxit=0.
20    if (maxit.eq.0) goto 160
c..
c..calculate the max-norm of the first solution.
      dd = 0.0
      do 30 i=1,n
       dd = max(dd,abs(x(i)))
30    continue
      dxmax = dd
c..
c..begin the iterative process.
      do 120 iterat=1,maxit
       d = dd
c..
c..calculate the residual vector.
       do 40 i=1,n
        r(i) = rhs(i)
40     continue
       if (mtype.eq.1) goto 60
       do 50 i=1,nz
        nrow = irnorg(i)
        ncol = icnorg(i)
        r(ncol) = r(ncol) - aorg(i)*x(nrow)
50     continue
       goto 80
c..
c..mtype=1.
60     do 70 i=1,nz
        nrow = irnorg(i)
        ncol = icnorg(i)
        r(nrow) = r(nrow) - aorg(i)*x(ncol)
70     continue
80     dres = 0.0
c..
c..find the max-norm of the residual vector.
       do 90 i=1,n
        dres = max(dres,abs(r(i)))
90     continue
c..
c..stop the calculations if the max-norm of the residual vector is zero.
       if (dres.eq.0.0) goto 150
c..
c..calculate the correction vector.
       noiter = noiter + 1
       call ma28cd(n,a,licn,icn,ikeep,r,w,mtype)
c..
c..find the max-norm of the correction vector.
       dd = 0.0
       do 100 i=1,n
        dd = max(dd,abs(r(i)))
100    continue
c..
c..check the convergence.
       if (dd.gt.d*cgce .and. iterat.ge.2) goto 130
       if (dxmax*10.0+dd.eq.dxmax*10.0) goto 140
c..
c..attempt to improve the solution.
       dxmax = 0.0
       do 110 i=1,n
        x(i) = x(i) + r(i)
        dxmax = max(dxmax,abs(x(i)))
110    continue
c..
c..check the stopping criterion; end of iteration loop
       if (dd.lt.prec*dxmax) goto 140
120   continue
c..
c..more than maxit iterations required.
      iflag = -16
      write (lp,99999) iflag,maxit
      goto 140
c..
c..convergence rate unacceptably slow.
130   iflag = -17
      conver = dd/d
      write (lp,99998) iflag,conver,cgce
c..
c..the iterative process is terminated.
140   errmax = dd
150   continue
160   return
      end
c..
c..
c..
c..
c..
      subroutine ma28dd(n,a,licn,ivect,jvect,nz,icn,lenr,lenrl,
     1                  lenoff,ip,iq,iw1,iw,w1,iflag)
      implicit double precision (a-h,o-z)
      save
c..
c..this subroutine need never be called by the user directly.
c..it sorts the user's matrix into the structure of the decomposed
c..form and checks for the presence of duplicate entries or
c..non-zeros lying outside the sparsity pattern of the decomposition
c..it also calculates the largest element in the input matrix.
c..
c..declare
      logical          lblock,grow,blockl
      integer          n,licn,nz,iw(n,2),idisp(2),icn(licn),ivect(nz),
     1                 jvect(nz),ip(n),iq(n),lenr(n),iw1(n,3),lenrl(n),
     2                 lenoff(n),iflag,lp,mp,i,ii,jj,inew,jnew,iblock,
     3                 iold,jold,j1,j2,idisp2,idummy,jdummy,midpt,jcomp
      double precision a(licn),zero,w1,aa
c..
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28gd/  idisp
      data             zero /0.0d0/
c..
c..formats
99997 format(1x,'element',i6,' ',i6,' was not in l/u pattern')
99998 format(1x,'non-zero',i7,' ',i6,' in zero off-diagonal block')
99999 format(1x,'element ',i6,' with value ',1pe22.14,/,
     1       1x,'has indices',i8,' ',i8,'out of range')
c..
c..iw1(i,3)  is set to the block in which row i lies and the inverse 
c..permutations to ip and iq are set in iw1(.,1) and iw1(.,2) resp.
c..pointers to beginning of the part of row i in diagonal and off-diagonal 
c..blocks are set in iw(i,2) and iw(i,1) resp.
      blockl  = lenoff(1).ge.0
      iblock  = 1
      iw(1,1) = 1
      iw(1,2) = idisp(1)
      do 10 i=1,n
       iw1(i,3) = iblock
       if (ip(i).lt.0) iblock = iblock + 1
       ii = iabs(ip(i)+0)
       iw1(ii,1) = i
       jj = iq(i)
       jj = iabs(jj)
       iw1(jj,2) = i
       if (i.eq.1) goto 10
       if (blockl) iw(i,1) = iw(i-1,1) + lenoff(i-1)
       iw(i,2) = iw(i-1,2) + lenr(i-1)
10    continue
c..
c..place each non-zero in turn into its correct location in the a/icn array.
      idisp2 = idisp(2)
      do 170 i=1,nz
       if (i.gt.idisp2) goto 20
       if (icn(i).lt.0) goto 170
20     iold = ivect(i)
       jold = jvect(i)
       aa = a(i)
c..
c..dummy loop for following a chain of interchanges. executed nz times.
       do 140 idummy=1,nz
        if (iold.le.n .and. iold.gt.0 .and. jold.le.n .and. 
     1      jold.gt.0) goto 30
        if (lp.ne.0) write (lp,99999) i, a(i), iold, jold
        iflag = -12
        goto 180
30      inew = iw1(iold,1)
        jnew = iw1(jold,2)
c..
c..are we in a valid block and is it diagonal or off-diagonal?
        if (iw1(inew,3)-iw1(jnew,3)) 40, 60, 50
40      iflag = -13
        if (lp.ne.0) write (lp,99998) iold, jold
        goto 180
50      j1 = iw(inew,1)
        j2 = j1 + lenoff(inew) - 1
        goto 110
c..
c..element is in diagonal block.
60      j1 = iw(inew,2)
        if (inew.gt.jnew) goto 70
        j2 = j1 + lenr(inew) - 1
        j1 = j1 + lenrl(inew)
        goto 110
70      j2 = j1 + lenrl(inew)
c..
c..binary search of ordered list  .. element in l part of row.
        do 100 jdummy=1,n
         midpt = (j1+j2)/2
         jcomp = iabs(icn(midpt)+0)
         if (jnew-jcomp) 80, 130, 90
80       j2 = midpt
         goto 100
90       j1 = midpt
100     continue
        iflag = -13
        if (lp.ne.0) write (lp,99997) iold, jold
        goto 180
c..
c..linear search ..  element in l part of row or off-diagonal blocks.
110     do 120 midpt=j1,j2
         if (iabs(icn(midpt)+0).eq.jnew) goto 130
120     continue
        iflag = -13
        if (lp.ne.0) write (lp,99997) iold, jold
        goto 180
c..
c..equivalent element of icn is in position midpt.
130     if (icn(midpt).lt.0) goto 160
        if (midpt.gt.nz .or. midpt.le.i) goto 150
        w1 = a(midpt)
        a(midpt) = aa
        aa = w1
        iold = ivect(midpt)
        jold = jvect(midpt)
        icn(midpt) = -icn(midpt)
140    continue
c..
150    a(midpt) = aa
       icn(midpt) = -icn(midpt)
       goto 170
160    a(midpt) = a(midpt) + aa
c..
c..set flag for duplicate elements; end of big loop
       iflag = n + 1
170   continue
c..
c..reset icn array  and zero elements in l/u but not in a. get max a element
180   w1 = zero
      do 200 i=1,idisp2
       if (icn(i).lt.0) goto 190
       a(i) = zero
       goto 200
190    icn(i) = -icn(i)
       w1 = max(w1,abs(a(i)))
200   continue
      return
      end
c..
c..
c..
c..
c..
      subroutine ma30ad(nn,icn,a,licn,lenr,lenrl,idisp,ip,iq,
     1                  irn,lirn,lenc,ifirst,lastr,nextr,lastc,
     2                  nextc,iptr,ipc,u,iflag)
      implicit double precision (a-h,o-z)
      save
c..
c..
c..if the user requires a more convenient data interface then the ma28
c..package should be used.  the ma28 subroutines call the ma30 routines after 
c..checking the user's input data and optionally using mc23a/ad to permute the 
c..matrix to block triangular form.
c..
c..this package of subroutines (ma30a/ad, ma30b/bd, ma30c/cd and ma30d/dd) 
c..performs operations pertinent to the solution of a general sparse n by n 
c..system of linear equations (i.e. solve ax=b). structually singular matrices 
c..are permitted including those with row or columns consisting entirely of 
c..zeros (i.e. including rectangular matrices).  it is assumed that the 
c..non-zeros of the matrix a do not differ widely in size. if necessary a 
c..prior call of the scaling subroutine mc19a/ad may be made.
c..
c..a discussion of the design of these subroutines is given by duff and reid 
c..(acm trans math software 5 pp 18-35,1979 (css 48)) while fuller details of 
c..the implementation are given in duff (harwell report aere-r 8730,1977).  
c..the additional pivoting option in ma30a/ad and the use of drop tolerances 
c..(see common block ma30i/id) were added to the package after joint work with 
c..duff,reid,schaumburg,wasniewski and zlatev, harwell report css 135, 1983.
c..
c..ma30a/ad performs the lu decomposition of the diagonal blocks of the
c..permutation paq of a sparse matrix a, where input permutations p1 and q1 
c..are used to define the diagonal blocks.  there may be non-zeros in the 
c..off-diagonal blocks but they are unaffected by ma30a/ad. p and p1 differ 
c..only within blocks as do q and q1. the permutations p1 and q1 may be found 
c..by calling mc23a/ad or the matrix may be treated as a single block by 
c..using p1=q1=i. the matrix non-zeros should be held compactly by rows, 
c..although it should be noted that the user can supply the matrix by columns
c..to get the lu decomposition of a transpose.
c..
c..this description of the following parameters should also be consulted for 
c..further information on most of the parameters of ma30b/bd and ma30c/cd:
c..
c..n    is an integer variable which must be set by the user to the order
c..     of the matrix.  it is not altered by ma30a/ad.
c..
c..icn  is an integer array of length licn. positions idisp(2) to
c..     licn must be set by the user to contain the column indices of
c..     the non-zeros in the diagonal blocks of p1*a*q1. those belonging
c..     to a single row must be contiguous but the ordering of column
c..     indices with each row is unimportant. the non-zeros of row i
c..     precede those of row i+1,i=1,.. ,n-1 and no wasted space is
c..     allowed between the rows.  on output the column indices of the
c..     lu decomposition of paq are held in positions idisp(1) to
c..     idisp(2), the rows are in pivotal order, and the column indices
c..     of the l part of each row are in pivotal order and precede those
c..     of u. again there is no wasted space either within a row or
c..     between the rows. icn(1) to icn(idisp(1)-1), are neither
c..     required nor altered. if mc23a/ad been called, these will hold
c..     information about the off-diagonal blocks.
c..
c..a    is a real/double precision array of length licn whose entries
c..     idisp(2) to licn must be set by the user to the  values of the
c..     non-zero entries of the matrix in the order indicated by  icn.
c..     on output a will hold the lu factors of the matrix where again
c..     the position in the matrix is determined by the corresponding
c..     values in icn. a(1) to a(idisp(1)-1) are neither required nor altered.
c..
c..licn is an integer variable which must be set by the user to the
c..     length of arrays icn and a. it must be big enough for a and icn
c..     to hold all the non-zeros of l and u and leave some "elbow
c..     room".  it is possible to calculate a minimum value for licn by
c..     a preliminary run of ma30a/ad. the adequacy of the elbow room
c..     can be judged by the size of the common block variable icncp. it
c..     is not altered by ma30a/ad.
c..
c..lenr is an integer array of length n.  on input, lenr(i) should
c..     equal the number of non-zeros in row i, i=1,.. ,n of the
c..     diagonal blocks of p1*a*q1. on output, lenr(i) will equal the
c..     total number of non-zeros in row i of l and row i of u.
c..
c..lenrl is an integer array of length n. on output from ma30a/ad,
c..      lenrl(i) will hold the number of non-zeros in row i of l.
c..
c..idisp is an integer array of length 2. the user should set idisp(1)
c..      to be the first available position in a/icn for the lu
c..      decomposition while idisp(2) is set to the position in a/icn of
c..      the first non-zero in the diagonal blocks of p1*a*q1. on output,
c..      idisp(1) will be unaltered while idisp(2) will be set to the
c..      position in a/icn of the last non-zero of the lu decomposition.
c..
c..ip    is an integer array of length n which holds a permutation of
c..      the integers 1 to n.  on input to ma30a/ad, the absolute value of
c..      ip(i) must be set to the row of a which is row i of p1*a*q1. a
c..      negative value for ip(i) indicates that row i is at the end of a
c..      diagonal block.  on output from ma30a/ad, ip(i) indicates the row
c..      of a which is the i th row in paq. ip(i) will still be negative
c..      for the last row of each block (except the last).
c..
c..iq    is an integer array of length n which again holds a
c..      permutation of the integers 1 to n.  on input to ma30a/ad, iq(j)
c..      must be set to the column of a which is column j of p1*a*q1. on
c..      output from ma30a/ad, the absolute value of iq(j) indicates the
c..      column of a which is the j th in paq.  for rows, i say, in which
c..      structural or numerical singularity is detected iq(i) is negated.
c..
c..irn  is an integer array of length lirn used as workspace by ma30a/ad.
c..
c..lirn is an integer variable. it should be greater than the
c..     largest number of non-zeros in a diagonal block of p1*a*q1 but
c..     need not be as large as licn. it is the length of array irn and
c..     should be large enough to hold the active part of any block,
c..     plus some "elbow room", the  a posteriori  adequacy of which can
c..     be estimated by examining the size of common block variable irncp.
c..
c..lenc,ifirst,lastr,nextr,lastc,nextc 
c..     are all integer arrays of length n which are used as workspace by 
c..     ma30a/ad.  if nsrch is set to a value less than or equal to n, then 
c..     arrays lastc and nextc are not referenced by ma30a/ad and so can be 
c..     dummied in the call to ma30a/ad.
c..
c..iptr,ipc are integer arrays of length n; used as workspace by ma30a/ad.
c..
c..u    is a real/double precision variable which should be set by the
c..     user to a value between 0. and 1.0. if less than zero it is
c..     reset to zero and if its value is 1.0 or greater it is reset to
c..     0.9999 (0.999999999 in d version).  it determines the balance
c..     between pivoting for sparsity and for stability, values near
c..     zero emphasizing sparsity and values near one emphasizing
c..     stability. we recommend u=0.1 as a posible first trial value.
c..     the stability can be judged by a later call to mc24a/ad or by
c..     setting lbig to .true.
c..
c..iflag is an integer variable. it will have a non-negative value if
c..      ma30a/ad is successful. negative values indicate error
c..      conditions while positive values indicate that the matrix has
c..      been successfully decomposed but is singular. for each non-zero
c..      value, an appropriate message is output on unit lp.  possible
c..      non-zero values for iflag are
c.. -1   the matrix is structually singular with rank given by irank in
c..      common block ma30f/fd.
c.. +1   if, however, the user wants the lu decomposition of a
c..      structurally singular matrix and sets the common block variable
c..      abort1 to .false., then, in the event of singularity and a
c..      successful decomposition, iflag is returned with the value +1
c..      and no message is output.
c.. -2   the matrix is numerically singular (it may also be structually
c..      singular) with estimated rank given by irank in common block ma30f/fd.
c.. +2   the  user can choose to continue the decomposition even when a
c..      zero pivot is encountered by setting common block variable
c..      abort2 to .false.  if a singularity is encountered, iflag will
c..      then return with a value of +2, and no message is output if the
c..      decomposition has been completed successfully.
c.. -3   lirn has not been large enough to continue with the
c..      decomposition.  if the stage was zero then common block variable
c..      minirn gives the length sufficient to start the decomposition on
c..      this block.  for a successful decomposition on this block the user
c..      should make lirn slightly (say about n/2) greater than this value.
c.. -4   licn not large enough to continue with the decomposition.
c.. -5   the decomposition has been completed but some of the lu factors
c..      have been discarded to create enough room in a/icn to continue
c..      the decomposition. the variable minicn in common block ma30f/fd
c..      then gives the size that licn should be to enable the
c..      factorization to be successful.  if the user sets common block
c..      variable abort3 to .true., then the subroutine will exit
c..      immediately instead of destroying any factors and continuing.
c.. -6   both licn and lirn are too small. termination has been caused by
c..      lack of space in irn (see error iflag= -3), but already some of
c..      the lu factors in a/icn have been lost (see error iflag= -5).
c..      minicn gives the minimum amount of space required in a/icn for
c..      decomposition up to this point.
c..
c..
c..declare
      logical          abort1,abort2,abort3,lbig
      integer          nn,licn,lirn,iptr(nn),pivot,pivend,dispc,
     1                 oldpiv,oldend,pivrow,rowi,ipc(nn),idisp(2),
     2                 colupd,icn(licn),lenr(nn),lenrl(nn),ip(nn),
     3                 iq(nn),lenc(nn),irn(lirn),ifirst(nn),lastr(nn),
     4                 nextr(nn),lastc(nn),nextc(nn),lpiv(10),lnpiv(10),
     5                 msrch,nsrch,ndrop,kk,mapiv,manpiv,iavpiv,ianpiv,
     6                 kountl,minirn,minicn,morei,irank,irncp,icncp,
     7                 iflag,ibeg,iactiv,nzrow,num,nnm1,i,ilast,nblock,
     8                 istart,irows,n,ising,lp,itop,ii,j1,jj,j,indrow,
     9                 j2,ipos,nzcol,nzmin,nz,isw,isw1,jcost,
     &                 isrch,ll,ijfir,idummy,kcost,ijpos,ipiv,jpiv,i1,
     1                 i2,jpos,k,lc,nc,lr,nr,lenpiv,ijp1,nz2,l,lenpp,
     2                 nzpc,iii,idrop,iend,iop,jnew,ifill,jdiff,jnpos,
     3                 jmore,jend,jroom,jbeg,jzero,idispc,kdrop,
     4                 ifir,jval,jzer,jcount,jdummy,jold
      double precision a(licn),u,au,umax,amax,zero,pivrat,pivr,
     1                 tol,big,anew,aanew,scale
c..
      common /ma30ed/  lp,abort1,abort2,abort3
      common /ma30fd/  irncp,icncp,irank,minirn,minicn
      common /ma30id/  tol,big,ndrop,nsrch,lbig
      common /lpivot/  lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
c..
      data             umax/0.999999999d0/
      data             zero /0.0d0/
c..
c..
c..formats
99992 format(1x,'to continue set lirn to at least',i8)
99993 format(1x,'at stage',i5,'in block',i5,'with first row',i5,
     1          'and last row',i5)
99994 format(1x,'error return from ma30a/ad lirn and licn too small')
99995 format(1x,'error return from ma30a/ad ; lirn not big enough')
99996 format(1x,'error return from ma30a/ad ; licn not big enough')
99997 format(1x,'lu decomposition destroyed to create more space')
99998 format(1x,'error return from ma30a/ad;',
     1          'matrix is numerically singular')
99999 format(1x,'error return from ma30a/ad;',
     1          'matrix is structurally singular')
c..
c..
c..
c..initialize
      msrch = nsrch
      ndrop = 0
      do 1272 kk=1,10
       lnpiv(kk) = 0
       lpiv(kk)  = 0
1272  continue
      mapiv  = 0
      manpiv = 0
      iavpiv = 0
      ianpiv = 0
      kountl = 0
      minirn = 0
      minicn = idisp(1) - 1
      morei  = 0
      irank  = nn
      irncp  = 0
      icncp  = 0
      iflag  = 0
      u      = min(u,umax)
      u      = max(u,zero)
c..
c..ibeg is the position of the next pivot row after elimination step using it.
c..iactiv is the position of the first entry in the active part of a/icn.
c..nzrow is current number of non-zeros in active and unprocessed part of row 
c..file icn.
      ibeg   = idisp(1)
      iactiv = idisp(2)
      nzrow  = licn - iactiv + 1
      minicn = nzrow + minicn
c..
c..count the number of diagonal blocks and set up pointers to the beginnings of 
c..the rows. num is the number of diagonal blocks.
      num = 1
      iptr(1) = iactiv
      if (nn.eq.1) goto 20
      nnm1 = nn - 1
      do 10 i=1,nnm1
       if (ip(i).lt.0) num = num + 1
       iptr(i+1) = iptr(i) + lenr(i)
10    continue
c..
c..ilast is the last row in the previous block.
20    ilast = 0
c..
c..lu decomposition of block nblock starts
c..each pass on this loop performs lu decomp on one of the diagonal blocks.
      do 1000 nblock=1,num
       istart = ilast + 1
       do 30 irows=istart,nn
        if (ip(irows).lt.0) goto 40
30     continue
       irows = nn
40     ilast = irows
c..
c..n is the number of rows in the current block.
c..istart is the index of the first row in the current block.
c..ilast is the index of the last row in the current block.
c..iactiv is the position of the first entry in the block.
c..itop is the position of the last entry in the block.
       n = ilast - istart + 1
       if (n.ne.1) goto 90
c..
c..code for dealing with 1x1 block.
       lenrl(ilast) = 0
       ising = istart
       if (lenr(ilast).ne.0) goto 50
c..
c..block is structurally singular.
       irank = irank - 1
       ising = -ising
       if (iflag.ne.2 .and. iflag.ne.-5) iflag = 1
       if (.not.abort1) goto 80
       idisp(2) = iactiv
       iflag = -1
       if (lp.ne.0) write (lp,99999)
       goto 1120
c..
50     scale = abs(a(iactiv))
       if (scale.eq.zero) goto 60
       if (lbig) big = max(big,scale)
       goto 70
60     ising = -ising
       irank = irank - 1
       iptr(ilast) = 0
       if (iflag.ne.-5) iflag = 2
       if (.not.abort2) goto 70
       idisp(2) = iactiv
       iflag = -2
       if (lp.ne.0) write (lp,99998)
       goto 1120
70     a(ibeg) = a(iactiv)
       icn(ibeg) = icn(iactiv)
       iactiv = iactiv + 1
       iptr(istart) = 0
       ibeg = ibeg + 1
       nzrow = nzrow - 1
80     lastr(istart) = istart
       ipc(istart) = -ising
       goto 1000
c..
c..non-trivial block.
90     itop = licn
       if (ilast.ne.nn) itop = iptr(ilast+1) - 1
c..
c..set up column oriented storage.
       do 100 i=istart,ilast
        lenrl(i) = 0
        lenc(i) = 0
100    continue
       if (itop-iactiv.lt.lirn) goto 110
       minirn = itop - iactiv + 1
       pivot = istart - 1
       goto 1100
c..
c..calculate column counts.
110    do 120 ii=iactiv,itop
        i = icn(ii)
        lenc(i) = lenc(i) + 1
120    continue
c..
c..set up column pointers so that ipc(j) points to position after end of 
c..column j in column file.
       ipc(ilast) = lirn + 1
       j1 = istart + 1
       do 130 jj=j1,ilast
        j = ilast - jj + j1 - 1
        ipc(j) = ipc(j+1) - lenc(j+1)
130    continue
       do 150 indrow=istart,ilast
        j1 = iptr(indrow)
        j2 = j1 + lenr(indrow) - 1
        if (j1.gt.j2) goto 150
        do 140 jj=j1,j2
         j = icn(jj)
         ipos = ipc(j) - 1
         irn(ipos) = indrow
         ipc(j) = ipos
140     continue
150    continue
c..
c..dispc is the lowest indexed active location in the column file.
       dispc = ipc(istart)
       nzcol = lirn - dispc + 1
       minirn = max0(nzcol,minirn)
       nzmin = 1
c..
c..initialize array ifirst.  ifirst(i) = +/- k indicates that row/col k has i 
c..non-zeros.  if ifirst(i) = 0,there is no row or column with i non zeros.
       do 160 i=1,n
        ifirst(i) = 0
160    continue
c..
c..compute ordering of row and column counts. first run through columns (from 
c..column n to column 1).
       do 180 jj=istart,ilast
        j = ilast - jj + istart
        nz = lenc(j)
        if (nz.ne.0) goto 170
        ipc(j) = 0
        goto 180
170     if (nsrch.le.nn) goto 180
        isw = ifirst(nz)
        ifirst(nz) = -j
        lastc(j) = 0
        nextc(j) = -isw
        isw1 = iabs(isw)
        if (isw.ne.0) lastc(isw1) = j
180    continue
c..
c..now run through rows (again from n to 1).
       do 210 ii=istart,ilast
        i = ilast - ii + istart
        nz = lenr(i)
        if (nz.ne.0) goto 190
        iptr(i) = 0
        lastr(i) = 0
        goto 210
190     isw = ifirst(nz)
        ifirst(nz) = i
        if (isw.gt.0) goto 200
        nextr(i) = 0
        lastr(i) = isw
        goto 210
200     nextr(i) = isw
        lastr(i) = lastr(isw)
        lastr(isw) = i
210    continue
c..
c..
c..start of main elimination loop
c..
c..first find the pivot using markowitz criterion with stability control.
c..jcost is the markowitz cost of the best pivot so far,.. this pivot is in 
c..row ipiv and column jpiv.
c..
       do 980 pivot=istart,ilast
        nz2 = nzmin
        jcost = n*n
c..
c..examine rows/columns in order of ascending count.
        do 340 l=1,2
         pivrat = zero
         isrch = 1
         ll = l
c..
c..a pass with l equal to 2 is only performed in the case of singularity.
         do 330 nz=nz2,n
          if (jcost.le.(nz-1)**2) goto 420
          ijfir = ifirst(nz)
          if (ijfir) 230, 220, 240
220       if (ll.eq.1) nzmin = nz + 1
          goto 330
230       ll = 2
          ijfir = -ijfir
          goto 290
240       ll = 2
c..
c.. scan rows with nz non-zeros.
          do 270 idummy=1,n
           if (jcost.le.(nz-1)**2) goto 420
           if (isrch.gt.msrch) goto 420
           if (ijfir.eq.0) goto 280
c..
c..row ijfir is now examined.
           i = ijfir
           ijfir = nextr(i)
c..
c..first calculate multiplier threshold level.
           amax = zero
           j1 = iptr(i) + lenrl(i)
           j2 = iptr(i) + lenr(i) - 1
           do 250 jj=j1,j2
            amax = max(amax,abs(a(jj)))
250        continue
           au = amax*u
           isrch = isrch + 1
c..
c..scan row for possible pivots
           do 260 jj=j1,j2
            if (abs(a(jj)).le.au .and. l.eq.1) goto 260
            j = icn(jj)
            kcost = (nz-1)*(lenc(j)-1)
            if (kcost.gt.jcost) goto 260
            pivr = zero
            if (amax.ne.zero) pivr = abs(a(jj))/amax
            if (kcost.eq.jcost .and. (pivr.le.pivrat .or. 
     1          nsrch.gt.nn+1)) goto 260
c..
c..best pivot so far is found.
            jcost = kcost
            ijpos = jj
            ipiv = i
            jpiv = j
            if (msrch.gt.nn+1 .and. jcost.le.(nz-1)**2) goto 420
            pivrat = pivr
260        continue
270       continue
c..
c..columns with nz non-zeros now examined.
280       ijfir = ifirst(nz)
          ijfir = -lastr(ijfir)
290       if (jcost.le.nz*(nz-1)) goto 420
          if (msrch.le.nn) goto 330
          do 320 idummy=1,n
           if (ijfir.eq.0) goto 330
           j = ijfir
           ijfir = nextc(ijfir)
           i1 = ipc(j)
           i2 = i1 + nz - 1
c..
c..scan column j
           do 310 ii=i1,i2
            i = irn(ii)
            kcost = (nz-1)*(lenr(i)-lenrl(i)-1)
            if (kcost.ge.jcost) goto 310
c..
c..pivot has best markowitz count so far ..  now check its suitability on 
c..numeric grounds by examining the other non-zeros in its row.
            j1 = iptr(i) + lenrl(i)
            j2 = iptr(i) + lenr(i) - 1
c..
c..we need a stability check on singleton columns because of possible problems 
c..with underdetermined systems.
            amax = zero
            do 300 jj=j1,j2
             amax = max(amax,abs(a(jj)))
             if (icn(jj).eq.j) jpos = jj
300         continue
            if (abs(a(jpos)).le.amax*u .and. l.eq.1) goto 310
            jcost = kcost
            ipiv = i
            jpiv = j
            ijpos = jpos
            if (amax.ne.zero) pivrat = abs(a(jpos))/amax
            if (jcost.le.nz*(nz-1)) goto 420
310        continue
320       continue
330      continue
c..
c..in the event of singularity; must make sure all rows and columns are tested.
c..matrix is numerically or structurally singular; it will be diagnosed later.
         msrch = n
         irank = irank - 1
340     continue
c..
c..assign rest of rows and columns to ordering array. matrix is singular.
        if (iflag.ne.2 .and. iflag.ne.-5) iflag = 1
        irank = irank - ilast + pivot + 1
        if (.not.abort1) goto 350
        idisp(2) = iactiv
        iflag = -1
        if (lp.ne.0) write (lp,99999)
        goto 1120
350     k = pivot - 1
        do 390 i=istart,ilast
         if (lastr(i).ne.0) goto 390
         k = k + 1
         lastr(i) = k
         if (lenrl(i).eq.0) goto 380
         minicn = max0(minicn,nzrow+ibeg-1+morei+lenrl(i))
         if (iactiv-ibeg.ge.lenrl(i)) goto 360
         call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
c..
c..check now to see if ma30d/dd has created enough available space.
         if (iactiv-ibeg.ge.lenrl(i)) goto 360
c..
c..create more space by destroying previously created lu factors.
         morei = morei + ibeg - idisp(1)
         ibeg = idisp(1)
         if (lp.ne.0) write (lp,99997)
         iflag = -5
         if (abort3) goto 1090
360      j1 = iptr(i)
         j2 = j1 + lenrl(i) - 1
         iptr(i) = 0
         do 370 jj=j1,j2
          a(ibeg) = a(jj)
          icn(ibeg) = icn(jj)
          icn(jj) = 0
          ibeg = ibeg + 1
370      continue
         nzrow = nzrow - lenrl(i)
380      if (k.eq.ilast) goto 400
390     continue
400     k = pivot - 1
        do 410 i=istart,ilast
         if (ipc(i).ne.0) goto 410
         k = k + 1
         ipc(i) = k
         if (k.eq.ilast) goto 990
410     continue
c..
c..the pivot has now been found in position (ipiv,jpiv) in location ijpos in 
c..row file. update column and row ordering arrays to correspond with removal
c..of the active part of the matrix.
420     ising = pivot
        if (a(ijpos).ne.zero) goto 430
c..
c..numerical singularity is recorded here.
        ising = -ising
        if (iflag.ne.-5) iflag = 2
        if (.not.abort2) goto 430
        idisp(2) = iactiv
        iflag = -2
        if (lp.ne.0) write (lp,99998)
        goto 1120
430     oldpiv = iptr(ipiv) + lenrl(ipiv)
        oldend = iptr(ipiv) + lenr(ipiv) - 1
c..
c..changes to column ordering.
        if (nsrch.le.nn) goto 460
        colupd = nn + 1
        lenpp = oldend-oldpiv+1
        if (lenpp.lt.4) lpiv(1) = lpiv(1) + 1
        if (lenpp.ge.4 .and. lenpp.le.6) lpiv(2) = lpiv(2) + 1
        if (lenpp.ge.7 .and. lenpp.le.10) lpiv(3) = lpiv(3) + 1
        if (lenpp.ge.11 .and. lenpp.le.15) lpiv(4) = lpiv(4) + 1
        if (lenpp.ge.16 .and. lenpp.le.20) lpiv(5) = lpiv(5) + 1
        if (lenpp.ge.21 .and. lenpp.le.30) lpiv(6) = lpiv(6) + 1
        if (lenpp.ge.31 .and. lenpp.le.50) lpiv(7) = lpiv(7) + 1
        if (lenpp.ge.51 .and. lenpp.le.70) lpiv(8) = lpiv(8) + 1
        if (lenpp.ge.71 .and. lenpp.le.100) lpiv(9) = lpiv(9) + 1
        if (lenpp.ge.101) lpiv(10) = lpiv(10) + 1
        mapiv = max0(mapiv,lenpp)
        iavpiv = iavpiv + lenpp
        do 450 jj=oldpiv,oldend
         j = icn(jj)
         lc = lastc(j)
         nc = nextc(j)
         nextc(j) = -colupd
         if (jj.ne.ijpos) colupd = j
         if (nc.ne.0) lastc(nc) = lc
         if (lc.eq.0) goto 440
         nextc(lc) = nc
         goto 450
440      nz = lenc(j)
         isw = ifirst(nz)
         if (isw.gt.0) lastr(isw) = -nc
         if (isw.lt.0) ifirst(nz) = -nc
450     continue
c..
c..changes to row ordering.
460     i1 = ipc(jpiv)
        i2 = i1 + lenc(jpiv) - 1
        do 480 ii=i1,i2
         i = irn(ii)
         lr = lastr(i)
         nr = nextr(i)
         if (nr.ne.0) lastr(nr) = lr
         if (lr.le.0) goto 470
         nextr(lr) = nr
         goto 480
470      nz = lenr(i) - lenrl(i)
         if (nr.ne.0) ifirst(nz) = nr
         if (nr.eq.0) ifirst(nz) = lr
480     continue
c..
c..move pivot to position lenrl+1 in pivot row and move pivot row to the 
c..beginning of the available storage. the l part and the pivot in the old 
c..copy of the pivot row is nullified while, in the strictly upper triangular 
c..part, the column indices, j say, are overwritten by the corresponding
c..entry of iq (iq(j)) and iq(j) is set to the negative of the displacement of 
c..the column index from the pivot entry.
        if (oldpiv.eq.ijpos) goto 490
        au = a(oldpiv)
        a(oldpiv) = a(ijpos)
        a(ijpos) = au
        icn(ijpos) = icn(oldpiv)
        icn(oldpiv) = jpiv
c..
c..check if there is space available in a/icn to hold new copy of pivot row.
490     minicn = max0(minicn,nzrow+ibeg-1+morei+lenr(ipiv))
        if (iactiv-ibeg.ge.lenr(ipiv)) goto 500
        call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
        oldpiv = iptr(ipiv) + lenrl(ipiv)
        oldend = iptr(ipiv) + lenr(ipiv) - 1
c..
c..check now to see if ma30d/dd has created enough available space.
        if (iactiv-ibeg.ge.lenr(ipiv)) goto 500
c..
c..create more space by destroying previously created lu factors.
        morei = morei + ibeg - idisp(1)
        ibeg = idisp(1)
        if (lp.ne.0) write (lp,99997)
        iflag = -5
        if (abort3) goto 1090
        if (iactiv-ibeg.ge.lenr(ipiv)) goto 500
c..
c..there is still not enough room in a/icn.
        iflag = -4
        goto 1090
c..
c..copy pivot row and set up iq array.
500     ijpos = 0
        j1 = iptr(ipiv)
        do 530 jj=j1,oldend
         a(ibeg) = a(jj)
         icn(ibeg) = icn(jj)
         if (ijpos.ne.0) goto 510
         if (icn(jj).eq.jpiv) ijpos = ibeg
         icn(jj) = 0
         goto 520
510      k = ibeg - ijpos
         j = icn(jj)
         icn(jj) = iq(j)
         iq(j) = -k
520      ibeg = ibeg + 1
530     continue
c..
        ijp1 = ijpos + 1
        pivend = ibeg - 1
        lenpiv = pivend - ijpos
        nzrow = nzrow - lenrl(ipiv) - 1
        iptr(ipiv) = oldpiv + 1
        if (lenpiv.eq.0) iptr(ipiv) = 0
c..
c..remove pivot row (including pivot) from column oriented file.
        do 560 jj=ijpos,pivend
         j = icn(jj)
         i1 = ipc(j)
         lenc(j) = lenc(j) - 1
c..
c..i2 is last position in new column.
         i2 = ipc(j) + lenc(j) - 1
         if (i2.lt.i1) goto 550
         do 540 ii=i1,i2
          if (irn(ii).ne.ipiv) goto 540
          irn(ii) = irn(i2+1)
          goto 550
540      continue
550      irn(i2+1) = 0
560     continue
        nzcol = nzcol - lenpiv - 1
c..
c..go down the pivot column and for each row with a non-zero add the 
c..appropriate multiple of the pivot row to it. we loop on the number of 
c..non-zeros in the pivot column since ma30d/dd may change its actual position.
        nzpc = lenc(jpiv)
        if (nzpc.eq.0) goto 900
        do 840 iii=1,nzpc
         ii = ipc(jpiv) + iii - 1
         i = irn(ii)
c..
c..search row i for non-zero to be eliminated, calculate multiplier, and place 
c..it in position lenrl+1 in its row. idrop is the number of non-zero entries 
c..dropped from row i because these fall beneath tolerance level.
         idrop = 0
         j1 = iptr(i) + lenrl(i)
         iend = iptr(i) + lenr(i) - 1
         do 570 jj=j1,iend
          if (icn(jj).ne.jpiv) goto 570
c..
c..if pivot is zero, rest of column is and so multiplier is zero.
          au = zero
          if (a(ijpos).ne.zero) au = -a(jj)/a(ijpos)
          if (lbig) big = max(big,abs(au))
          a(jj) = a(j1)
          a(j1) = au
          icn(jj) = icn(j1)
          icn(j1) = jpiv
          lenrl(i) = lenrl(i) + 1
          goto 580
570      continue
c..
c..jump if pivot row is a singleton.
580      if (lenpiv.eq.0) goto 840
c..
c..now perform necessary operations on rest of non-pivot row i.
         rowi = j1 + 1
         iop = 0
c..
c..jump if all the pivot row causes fill-in.
         if (rowi.gt.iend) goto 650
c..
c..perform operations on current non-zeros in row i. innermost loop.
         lenpp = iend-rowi+1
         if (lenpp.lt.4) lnpiv(1) = lnpiv(1) + 1
         if (lenpp.ge.4 .and. lenpp.le.6) lnpiv(2) = lnpiv(2) + 1
         if (lenpp.ge.7 .and. lenpp.le.10) lnpiv(3) = lnpiv(3) + 1
         if (lenpp.ge.11 .and. lenpp.le.15) lnpiv(4) = lnpiv(4) + 1
         if (lenpp.ge.16 .and. lenpp.le.20) lnpiv(5) = lnpiv(5) + 1
         if (lenpp.ge.21 .and. lenpp.le.30) lnpiv(6) = lnpiv(6) + 1
         if (lenpp.ge.31 .and. lenpp.le.50) lnpiv(7) = lnpiv(7) + 1
         if (lenpp.ge.51 .and. lenpp.le.70) lnpiv(8) = lnpiv(8) + 1
         if (lenpp.ge.71 .and. lenpp.le.100) lnpiv(9) = lnpiv(9) + 1
         if (lenpp.ge.101) lnpiv(10) = lnpiv(10) + 1
         manpiv = max0(manpiv,lenpp)
         ianpiv = ianpiv + lenpp
         kountl = kountl + 1
         do 590 jj=rowi,iend
          j = icn(jj)
          if (iq(j).gt.0) goto 590
          iop = iop + 1
          pivrow = ijpos - iq(j)
          a(jj) = a(jj) + au*a(pivrow)
          if (lbig) big = max(abs(a(jj)),big)
          icn(pivrow) = -icn(pivrow)
          if (abs(a(jj)).lt.tol) idrop = idrop + 1
590      continue
c..
c..jump if no non-zeros in non-pivot row have been removed because these are 
c..beneath the drop-tolerance  tol.
         if (idrop.eq.0) goto 650
c..
c..run through non-pivot row compressing row so that only non-zeros greater 
c..than tol are stored. all non-zeros less than tol are also removed from the 
c..column structure.
         jnew = rowi
         do 630 jj=rowi,iend
          if (abs(a(jj)).lt.tol) goto 600
          a(jnew) = a(jj)
          icn(jnew) = icn(jj)
          jnew = jnew + 1
          goto 630
c..
c..remove non-zero entry from column structure.
600       j = icn(jj)
          i1 = ipc(j)
          i2 = i1 + lenc(j) - 1
          do 610 ii=i1,i2
           if (irn(ii).eq.i) goto 620
610       continue
620       irn(ii) = irn(i2)
          irn(i2) = 0
          lenc(j) = lenc(j) - 1
          if (nsrch.le.nn) goto 630
c..
c..remove column from column chain and place in update chain.
          if (nextc(j).lt.0) goto 630
c..
c..jump if column already in update chain.
          lc = lastc(j)
          nc = nextc(j)
          nextc(j) = -colupd
          colupd = j
          if (nc.ne.0) lastc(nc) = lc
          if (lc.eq.0) goto 622
          nextc(lc) = nc
          goto 630
622       nz = lenc(j) + 1
          isw = ifirst(nz)
          if (isw.gt.0) lastr(isw) = -nc
          if (isw.lt.0) ifirst(nz) = -nc
630      continue
         do 640 jj=jnew,iend
          icn(jj) = 0
640      continue
c..
c..the value of idrop might be different from that calculated earlier because, 
c..we may have dropped some non-zeros which were not modified by the pivot row.
         idrop = iend + 1 - jnew
         iend = jnew - 1
         lenr(i) = lenr(i) - idrop
         nzrow = nzrow - idrop
         nzcol = nzcol - idrop
         ndrop = ndrop + idrop
650      ifill = lenpiv - iop
c..
c..jump is if there is no fill-in.
         if (ifill.eq.0) goto 750
c..
c..now for the fill-in.
         minicn = max0(minicn,morei+ibeg-1+nzrow+ifill+lenr(i))
c..
c..see if there is room for fill-in. get maximum space for row i in situ.
         do 660 jdiff=1,ifill
          jnpos = iend + jdiff
          if (jnpos.gt.licn) goto 670
          if (icn(jnpos).ne.0) goto 670
660      continue
c..
c..there is room for all the fill-in after the end of the row so it can be 
c..left in situ. next available space for fill-in.
         iend = iend + 1
         goto 750
c..
c..jmore spaces for fill-in are required in front of row.
670      jmore = ifill - jdiff + 1
         i1 = iptr(i)
c..
c..look in front of the row to see if there is space for rest of the fill-in.
         do 680 jdiff=1,jmore
          jnpos = i1 - jdiff
          if (jnpos.lt.iactiv) goto 690
          if (icn(jnpos).ne.0) goto 700
680      continue
690      jnpos = i1 - jmore
         goto 710
c..
c..whole row must be moved to the beginning of available storage.
700      jnpos = iactiv - lenr(i) - ifill
c..
c..jump if there is space immediately available for the shifted row.
710      if (jnpos.ge.ibeg) goto 730
         call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
         i1 = iptr(i)
         iend = i1 + lenr(i) - 1
         jnpos = iactiv - lenr(i) - ifill
         if (jnpos.ge.ibeg) goto 730
c..
c..no space available; try to create some by trashing previous lu decomposition.
         morei = morei + ibeg - idisp(1) - lenpiv - 1
         if (lp.ne.0) write (lp,99997)
         iflag = -5
         if (abort3) goto 1090
c..
c..keep record of current pivot row.
         ibeg = idisp(1)
         icn(ibeg) = jpiv
         a(ibeg) = a(ijpos)
         ijpos = ibeg
         do 720 jj=ijp1,pivend
          ibeg = ibeg + 1
          a(ibeg) = a(jj)
          icn(ibeg) = icn(jj)
  720    continue
         ijp1 = ijpos + 1
         pivend = ibeg
         ibeg = ibeg + 1
         if (jnpos.ge.ibeg) goto 730
c..
c..this still does not give enough room.
         iflag = -4
         goto 1090
730      iactiv = min0(iactiv,jnpos)
c..
c..move non-pivot row i.
         iptr(i) = jnpos
         do 740 jj=i1,iend
          a(jnpos) = a(jj)
          icn(jnpos) = icn(jj)
          jnpos = jnpos + 1
          icn(jj) = 0
740      continue
c..
c..first new available space.
         iend = jnpos
750      nzrow = nzrow + ifill
c..
c..innermost fill-in loop which also resets icn.
         idrop = 0
         do 830 jj=ijp1,pivend
          j = icn(jj)
          if (j.lt.0) goto 820
          anew = au*a(jj)
          aanew = abs(anew)
          if (aanew.ge.tol) goto 760
          idrop = idrop + 1
          ndrop = ndrop + 1
          nzrow = nzrow - 1
          minicn = minicn - 1
          ifill = ifill - 1
          goto 830
760       if (lbig) big = max(aanew,big)
          a(iend) = anew
          icn(iend) = j
          iend = iend + 1
c..
c..put new entry in column file.
          minirn = max0(minirn,nzcol+lenc(j)+1)
          jend = ipc(j) + lenc(j)
          jroom = nzpc - iii + 1 + lenc(j)
          if (jend.gt.lirn) goto 770
          if (irn(jend).eq.0) goto 810
770       if (jroom.lt.dispc) goto 780
c..
c..compress column file to obtain space for new copy of column.
          call ma30dd(a,irn,ipc(istart),n,dispc,lirn,.false.)
          if (jroom.lt.dispc) goto 780
          jroom = dispc - 1
          if (jroom.ge.lenc(j)+1) goto 780
c..
c..column file is not large enough.
          goto 1100
c..
c..copy column to beginning of file.
780       jbeg = ipc(j)
          jend = ipc(j) + lenc(j) - 1
          jzero = dispc - 1
          dispc = dispc - jroom
          idispc = dispc
          do 790 ii=jbeg,jend
           irn(idispc) = irn(ii)
           irn(ii) = 0
           idispc = idispc + 1
790       continue
          ipc(j) = dispc
          jend = idispc
          do 800 ii=jend,jzero
           irn(ii) = 0
800       continue
810       irn(jend) = i
          nzcol = nzcol + 1
          lenc(j) = lenc(j) + 1
c..
c..end of adjustment to column file.
          goto 830
c..
820       icn(jj) = -j
830      continue
         if (idrop.eq.0) goto 834
         do 832 kdrop=1,idrop
          icn(iend) = 0
          iend = iend + 1
832      continue
834      lenr(i) = lenr(i) + ifill
c..
c..end of scan of pivot column.
840     continue
c..
c..
c..remove pivot column from column oriented storage; update row ordering arrays.
        i1 = ipc(jpiv)
        i2 = ipc(jpiv) + lenc(jpiv) - 1
        nzcol = nzcol - lenc(jpiv)
        do 890 ii=i1,i2
         i = irn(ii)
         irn(ii) = 0
         nz = lenr(i) - lenrl(i)
         if (nz.ne.0) goto 850
         lastr(i) = 0
         goto 890
850      ifir = ifirst(nz)
         ifirst(nz) = i
         if (ifir) 860, 880, 870
860      lastr(i) = ifir
         nextr(i) = 0
         goto 890
870      lastr(i) = lastr(ifir)
         nextr(i) = ifir
         lastr(ifir) = i
         goto 890
880      lastr(i) = 0
         nextr(i) = 0
         nzmin = min0(nzmin,nz)
890     continue
c..
c..restore iq and nullify u part of old pivot row. record the column 
c..permutation in lastc(jpiv) and the row permutation in lastr(ipiv).
900     ipc(jpiv) = -ising
        lastr(ipiv) = pivot
        if (lenpiv.eq.0) goto 980
        nzrow = nzrow - lenpiv
        jval = ijp1
        jzer = iptr(ipiv)
        iptr(ipiv) = 0
        do 910 jcount=1,lenpiv
         j = icn(jval)
         iq(j) = icn(jzer)
         icn(jzer) = 0
         jval = jval + 1
         jzer = jzer + 1
910     continue
c..
c..adjust column ordering arrays.
        if (nsrch.gt.nn) goto 920
        do 916 jj=ijp1,pivend
         j = icn(jj)
         nz = lenc(j)
         if (nz.ne.0) goto 914
         ipc(j) = 0
         goto 916
914      nzmin = min0(nzmin,nz)
916     continue
        goto 980
920     jj = colupd
        do 970 jdummy=1,nn
         j = jj
         if (j.eq.nn+1) goto 980
         jj = -nextc(j)
         nz = lenc(j)
         if (nz.ne.0) goto 924
         ipc(j) = 0
         goto 970
924      ifir = ifirst(nz)
         lastc(j) = 0
         if (ifir) 930, 940, 950
930      ifirst(nz) = -j
         ifir = -ifir
         lastc(ifir) = j
         nextc(j) = ifir
         goto 970
940      ifirst(nz) = -j
         nextc(j) = 0
         goto 960
950      lc = -lastr(ifir)
         lastr(ifir) = -j
         nextc(j) = lc
         if (lc.ne.0) lastc(lc) = j
960      nzmin = min0(nzmin,nz)
970     continue
980    continue
c..
c..that was the end of main elimination loop
c..
c..
c..reset iactiv to point to the beginning of the next block.
990    if (ilast.ne.nn) iactiv = iptr(ilast+1)
1000  continue
c..
c..that was the end of deomposition of block   
c..
c..
c..record singularity (if any) in iq array.
      if (irank.eq.nn) goto 1020
      do 1010 i=1,nn
       if (ipc(i).lt.0) goto 1010
       ising = ipc(i)
       iq(ising) = -iq(ising)
       ipc(i) = -ising
1010  continue
c..
c..
c..run through lu decomposition changing column indices to that of new order 
c..and permuting lenr and lenrl arrays according to pivot permutations.
1020  istart = idisp(1)
      iend = ibeg - 1
      if (iend.lt.istart) goto 1040
      do 1030 jj=istart,iend
       jold = icn(jj)
       icn(jj) = -ipc(jold)
1030  continue
1040  do 1050 ii=1,nn
       i = lastr(ii)
       nextr(i) = lenr(ii)
       iptr(i) = lenrl(ii)
1050  continue
      do 1060 i=1,nn
       lenrl(i) = iptr(i)
       lenr(i) = nextr(i)
1060  continue
c..
c..update permutation arrays ip and iq.
      do 1070 ii=1,nn
       i = lastr(ii)
       j = -ipc(ii)
       nextr(i) = iabs(ip(ii)+0)
       iptr(j) = iabs(iq(ii)+0)
1070  continue
      do 1080 i=1,nn
       if (ip(i).lt.0) nextr(i) = -nextr(i)
       ip(i) = nextr(i)
       if (iq(i).lt.0) iptr(i) = -iptr(i)
       iq(i) = iptr(i)
1080  continue
      ip(nn) = iabs(ip(nn)+0)
      idisp(2) = iend
      goto 1120
c..
c..
c..error returns
1090  idisp(2) = iactiv
      if (lp.eq.0) goto 1120
      write (lp,99996)
      goto 1110
1100  if (iflag.eq.-5) iflag = -6
      if (iflag.ne.-6) iflag = -3
      idisp(2) = iactiv
      if (lp.eq.0) goto 1120
      if (iflag.eq.-3) write (lp,99995)
      if (iflag.eq.-6) write (lp,99994)
1110  pivot = pivot - istart + 1
      write (lp,99993) pivot,nblock,istart,ilast
      if (pivot.eq.0) write (lp,99992) minirn
1120  return
      end
c..
c..
c..
c..
c..
      subroutine ma30bd(n,icn,a,licn,lenr,lenrl,idisp,ip,iq,w,iw,iflag)
      implicit double precision (a-h,o-z)
      save
c..
c..ma30b/bd performs the lu decomposition of the diagonal blocks of a new 
c..matrix paq of the same sparsity pattern, using information from a previous 
c..call to ma30a/ad. the entries of the input matrix  must already be in their 
c..final positions in the lu decomposition structure.  this routine executes 
c..about five times faster than ma30a/ad.
c..
c..parameters (see also ma30ad):
c..n   is an integer variable set to the order of the matrix.
c..
c..icn is an integer array of length licn. it should be unchanged
c..    since the last call to ma30a/ad. it is not altered by ma30b/bd.
c..
c..a   is a real/double precision array of length licn the user must set
c..    entries idisp(1) to idisp(2) to contain the entries in the
c..    diagonal blocks of the matrix paq whose column numbers are held
c..    in icn, using corresponding positions. note that some zeros may
c..    need to be held explicitly. on output entries idisp(1) to
c..    idisp(2) of array a contain the lu decomposition of the diagonal
c..    blocks of paq. entries a(1) to a(idisp(1)-1) are neither
c..    required nor altered by ma30b/bd.
c..
c..licn is an integer variable which must be set by the user to the
c..     length of arrays a and icn. it is not altered by ma30b/bd.
c..
c..lenr,lenrl are integer arrays of length n. they should be
c..    unchanged since the last call to ma30a/ad. not altered by ma30b/bd.
c..
c..idisp is an integer array of length 2. it should be unchanged since
c..      the last call to ma30a/ad. it is not altered by ma30b/bd.
c..
c..ip,iq are integer arrays of length n. they should be unchanged
c..      since the last call to ma30a/ad. not altered by ma30b/bd.
c..
c..w   is a array of length n which is used as workspace by ma30b/bd.
c..
c..iw  is an integer array of length n which is used as workspace by ma30b/bd.
c..
c..iflag  is an integer variable. on output from ma30b/bd, iflag has
c..       the value zero if the factorization was successful, has the
c..       value i if pivot i was very small and has the value -i if an
c..       unexpected singularity was detected at stage i of the decomposition.
c..
c..declare
      logical          abort1,abort2,abort3,stab,lbig
      integer          n,licn,iflag,iw(n),idisp(2),pivpos,icn(licn),
     1                 lenr(n),lenrl(n),ip(n),iq(n),ndrop,nsrch,ising,
     2                 i,istart,ifin,ilend,j,ipivj,jfin,jay,jayjay,
     3                 jj,lp
      double precision a(licn),w(n),au,eps,rowmax,zero,one,rmin,tol,big
c..
      common /ma30ed/  lp,abort1,abort2,abort3
      common /ma30id/  tol,big,ndrop,nsrch,lbig
      common /ma30gd/  eps,rmin
      data             zero /0.0d0/,one /1.0d0/
c..
c..formats
99999 format(1x,'error return from ma30b/bd singularity in row',i8)
c..
c..initialize
      stab  = eps.le.one
      rmin  = eps
      ising = 0
      iflag = 0
      do 10 i=1,n
       w(i) = zero
10    continue
c..
c..set up pointers to the beginning of the rows.
      iw(1) = idisp(1)
      if (n.eq.1) goto 25
      do 20 i=2,n
       iw(i) = iw(i-1) + lenr(i-1)
20    continue
c..
c..start  of main loop
c..at step i, row i of a is transformed to row i of l/u by adding appropriate 
c..multiples of rows 1 to i-1. using row-gauss elimination.
c..istart is beginning of row i of a and row i of l.
c..ifin is end of row i of a and row i of u.
c..ilend is end of row i of l.
25    do 160 i=1,n
       istart = iw(i)
       ifin = istart + lenr(i) - 1
       ilend = istart + lenrl(i) - 1
       if (istart.gt.ilend) goto 90
c..
c..load row i of a into vector w.
       do 30 jj=istart,ifin
        j = icn(jj)
        w(j) = a(jj)
30     continue
c..
c..add multiples of appropriate rows of  i to i-1  to row i.
c..ipivj is position of pivot in row j. 
       do 70 jj=istart,ilend
        j = icn(jj)
        ipivj = iw(j) + lenrl(j)
        au = -w(j)/a(ipivj)
        if (lbig) big = max(abs(au),big)
        w(j) = au
c..
c..au * row j (u part) is added to row i.
        ipivj = ipivj + 1
        jfin = iw(j) + lenr(j) - 1
        if (ipivj.gt.jfin) goto 70
c..
c..innermost loop.
        if (lbig) goto 50
        do 40 jayjay=ipivj,jfin
         jay = icn(jayjay)
         w(jay) = w(jay) + au*a(jayjay)
40      continue
        goto 70
50      do 60 jayjay=ipivj,jfin
         jay    = icn(jayjay)
         w(jay) = w(jay) + au*a(jayjay)
         big    = max(abs(w(jay)),big)
60      continue
70     continue
c..
c..reload w back into a (now l/u)
       do 80 jj=istart,ifin
        j = icn(jj)
        a(jj) = w(j)
        w(j) = zero
80     continue
c..
c..now perform the stability checks.
90     pivpos = ilend + 1
       if (iq(i).gt.0) goto 140
c..
c..matrix had singularity at this point in ma30a/ad.
c..is it the first such pivot in current block ?
       if (ising.eq.0) ising = i
c..
c..does current matrix have a singularity in the same place ?
       if (pivpos.gt.ifin) goto 100
       if (a(pivpos).ne.zero) goto 170
c..
c..it does .. so set ising if it is not the end of the current block
c..check to see that appropriate part of l/u is zero or null.
100    if (istart.gt.ifin) goto 120
       do 110 jj=istart,ifin
        if (icn(jj).lt.ising) goto 110
        if (a(jj).ne.zero) goto 170
110    continue
120    if (pivpos.le.ifin) a(pivpos) = one
       if (ip(i).gt.0 .and. i.ne.n) goto 160
c..
c..end of current block ..  reset zero pivots and ising.
       do 130 j=ising,i
        if ((lenr(j)-lenrl(j)).eq.0) goto 130
        jj = iw(j) + lenrl(j)
        a(jj) = zero
130    continue
       ising = 0
       goto 160
c..
c..matrix had non-zero pivot in ma30a/ad at this stage.
140    if (pivpos.gt.ifin) goto 170
       if (a(pivpos).eq.zero) goto 170
       if (.not.stab) goto 160
       rowmax = zero
       do 150 jj=pivpos,ifin
        rowmax = max(rowmax,abs(a(jj)))
150    continue
       if (abs(a(pivpos))/rowmax.ge.rmin) goto 160
       iflag = i
       rmin = abs(a(pivpos))/rowmax
160   continue
      goto 180
c..
c..error return
170   if (lp.ne.0) write (lp,99999) i
      iflag = -i
180   return
      end
c..
c..
c..
c..
c..
      subroutine ma30cd(n,icn,a,licn,lenr,lenrl,lenoff,idisp,ip,
     1                  iq,x,w,mtype)
      implicit double precision (a-h,o-z)
      save
c..
c..ma30c/cd uses the factors produced by ma30a/ad or ma30b/bd to solve
c..ax=b or a transpose x=b when the matrix p1*a*q1 (paq) is block lower 
c..triangular (including the case of only one diagonal block).
c..
c..parameters: 
c..n  is an integer variable set to the order of the matrix. it is not
c..   altered by the subroutine.
c..
c..icn is an integer array of length licn. entries idisp(1) to
c..    idisp(2) should be unchanged since the last call to ma30a/ad. if
c..    the matrix has more than one diagonal block, then column indices
c..    corresponding to non-zeros in sub-diagonal blocks of paq must
c..    appear in positions 1 to idisp(1)-1. for the same row those
c..    entries must be contiguous, with those in row i preceding those
c..    in row i+1 (i=1,.. ,n-1) and no wasted space between rows.
c..    entries may be in any order within each row. not altered by ma30c/cd.
c..
c..a  is a real/double precision array of length licn.  entries
c..   idisp(1) to idisp(2) should be unchanged since the last call to
c..   ma30a/ad or ma30b/bd.  if the matrix has more than one diagonal
c..   block, then the values of the non-zeros in sub-diagonal blocks
c..   must be in positions 1 to idisp(1)-1 in the order given by icn.
c..   it is not altered by ma30c/cd.
c..
c..licn  is an integer variable set to the size of arrays icn and a.
c..      it is not altered by ma30c/cd.
c..
c..lenr,lenrl are integer arrays of length n which should be
c..     unchanged since the last call to ma30a/ad. not altered by ma30c/cd.
c..
c..lenoff  is an integer array of length n. if the matrix paq (or
c..        p1*a*q1) has more than one diagonal block, then lenoff(i),
c..        i=1,.. ,n should be set to the number of non-zeros in row i of
c..        the matrix paq which are in sub-diagonal blocks.  if there is
c..        only one diagonal block then lenoff(1) may be set to -1, in
c..        which case the other entries of lenoff are never accessed. it is
c..        not altered by ma30c/cd.
c..
c..idisp  is an integer array of length 2 which should be unchanged
c..       since the last call to ma30a/ad. it is not altered by ma30c/cd.
c..
cc..ip,iq are integer arrays of length n which should be unchanged
c..       since the last call to ma30a/ad. they are not altered by ma30c/cd.
c..
c..x   is a real/double precision array of length n. it must be set by
c..    the user to the values of the right hand side vector b for the
c..    equations being solved.  on exit from ma30c/cd it will be equal
c..    to the solution x required.
c..
c..w  is a real/double precision array of length n which is used as
c..   workspace by ma30c/cd.
c..
c mtype is an integer variable which must be set by the user. if
c..     mtype=1, then the solution to the system ax=b is returned; any
c..     other value for mtype will return the solution to the system a
c..     transpose x=b. it is not altered by ma30c/cd.
c..
c..declare
      logical          neg,nobloc
      integer          n,licn,idisp(2),icn(licn),lenr(n),lenrl(n),
     1                 lenoff(n),ip(n),iq(n),mtype,ii,i,lt,ifirst,
     2                 iblock,ltend,jj,j,iend,j1,ib,iii,j2,jpiv,
     3                 jpivp1,ilast,iblend,numblk,k,j3,iback,lj2,lj1
      double precision a(licn),x(n),w(n),wii,wi,resid,zero
      common /ma30hd/  resid
      data zero       /0.0d0/
c..
c..final value of resid is the max residual for inconsistent set of equations.
      resid = zero
c..
c..nobloc is .true. if subroutine block has been used previously and is .false. 
c..otherwise.  the value .false. means that lenoff will not be subsequently 
c..accessed.
      nobloc = lenoff(1).lt.0
      if (mtype.ne.1) goto 140
c..
c..now solve   a * x = b. neg is used to indicate when the last row in a block 
c..has been reached.  it is then set to true whereafter backsubstitution is
c..performed on the block.
      neg = .false.
c..
c..ip(n) is negated so that the last row of the last block can be recognised.  
c..it is reset to its positive value on exit.
      ip(n) = -ip(n)
c..
c..preorder vector ..  w(i) = x(ip(i))
      do 10 ii=1,n
       i = ip(ii)
       i = iabs(i)
       w(ii) = x(i)
10    continue
c..
c..lt is the position of first non-zero in current row of off-diagonal blocks.
c..ifirst holds the index of the first row in the current block.
c..iblock holds the position of the first non-zero in the current row
c..of the lu decomposition of the diagonal blocks.
      lt = 1
      ifirst = 1
      iblock = idisp(1)
c..
c..if i is not the last row of a block, then a pass through this loop adds the 
c..inner product of row i of the off-diagonal blocks and w to w and performs 
c..forward elimination using row i of the lu decomposition.   if i is the last 
c..row of a block then, after performing these aforementioned operations, 
c..backsubstitution is performed using the rows of the block.
      do 120 i=1,n
       wi = w(i)
       if (nobloc) goto 30
       if (lenoff(i).eq.0) goto 30
c..
c..operations using lower triangular blocks. 
c..ltend is the end of row i in the off-diagonal blocks.
       ltend = lt + lenoff(i) - 1
       do 20 jj=lt,ltend
        j = icn(jj)
        wi = wi - a(jj)*w(j)
20     continue
c..
c..lt is set the beginning of the next off-diagonal row.
c..set neg to .true. if we are on the last row of the block.
       lt = ltend + 1
30     if (ip(i).lt.0) neg = .true.
       if (lenrl(i).eq.0) goto 50
c..
c..forward elimination phase.
c..iend is the end of the l part of row i in the lu decomposition.
       iend = iblock + lenrl(i) - 1
       do 40 jj=iblock,iend
        j = icn(jj)
        wi = wi + a(jj)*w(j)
40     continue
c..
c..iblock is adjusted to point to the start of the next row.
50     iblock = iblock + lenr(i)
       w(i) = wi
       if (.not.neg) goto 120
c..
c..back substitution phase.
c..j1 is position in a/icn after end of block beginning in row ifirst
c..and ending in row i.
       j1 = iblock
c..
c..are there any singularities in this block?  if not, continue
       ib = i
       if (iq(i).gt.0) goto 70
       do 60 iii=ifirst,i
        ib = i - iii + ifirst
        if (iq(ib).gt.0) goto 70
        j1 = j1 - lenr(ib)
        resid = max(resid,abs(w(ib)))
        w(ib) = zero
60     continue
c..
c..entire block is singular.
       goto 110
c..
c..
c..each pass through this loop performs the back-substitution
c..operations for a single row, starting at the end of the block and
c..working through it in reverse order.
c..j2 is end of row ii. j1 is beginning of row ii. jpiv is the position of the 
c..pivot in row ii. jump out if row ii of u has no non-zeros.
70     do 100 iii=ifirst,ib
        ii = ib - iii + ifirst
        j2 = j1 - 1
        j1 = j1 - lenr(ii)
        jpiv = j1 + lenrl(ii)
        jpivp1 = jpiv + 1
        if (j2.lt.jpivp1) goto 90
        wii = w(ii)
        do 80 jj=jpivp1,j2
         j = icn(jj)
         wii = wii - a(jj)*w(j)
80      continue
        w(ii) = wii
90      w(ii) = w(ii)/a(jpiv)
100    continue
110    ifirst = i + 1
       neg = .false.
120   continue
c..
c..reorder solution vector ..  x(i) = w(iqinverse(i))
      do 130 ii=1,n
       i = iq(ii)
       i = iabs(i)
       x(i) = w(ii)
130   continue
      ip(n) = -ip(n)
      goto 320
c..
c..
c..now solve  atranspose * x = b. preorder vector ..  w(i)=x(iq(i))
140   do 150 ii=1,n
       i = iq(ii)
       i = iabs(i)
       w(ii) = x(i)
150   continue
c..
c..lj1 points to the beginning the current row in the off-diagonal blocks.
c..iblock is initialized to point to beginning of block after the last one
c..ilast is the last row in the current block.
c..iblend points to the position after the last non-zero in the current block.
      lj1 = idisp(1)
      iblock = idisp(2) + 1
      ilast = n
      iblend = iblock
c..
c..each pass through this loop operates with one diagonal block and
c..the off-diagonal part of the matrix corresponding to the rows
c..of this block.  the blocks are taken in reverse order and the
c..number of times the loop is entered is min(n,no. blocks+1).
      do 290 numblk=1,n
       if (ilast.eq.0) goto 300
       iblock = iblock - lenr(ilast)
c..
c..this loop finds the index of the first row in the current block. it is 
c..first and iblock is set to the position of the beginning of this first row.
       do 160 k=1,n
        ii = ilast - k
        if (ii.eq.0) goto 170
        if (ip(ii).lt.0) goto 170
        iblock = iblock - lenr(ii)
160    continue
170    ifirst = ii + 1
c..
c..j1 points to the position of the beginning of row i (lt part) or pivot
       j1 = iblock
c..
c..forward elimination. each pass through this loop performs the operations 
c..for one row of the block.  if the corresponding entry of w is zero then the
c..operations can be avoided.
       do 210 i=ifirst,ilast
        if (w(i).eq.zero) goto 200
c..
c.. jump if row i singular.
        if (iq(i).lt.0) goto 220
c..
c..j2 first points to the pivot in row i and then is made to point to the
c..first non-zero in the u transpose part of the row.
c..j3 points to the end of row i.
        j2 = j1 + lenrl(i)
        wi = w(i)/a(j2)
        if (lenr(i)-lenrl(i).eq.1) goto 190
        j2 = j2 + 1
        j3 = j1 + lenr(i) - 1
        do 180 jj=j2,j3
         j = icn(jj)
         w(j) = w(j) - a(jj)*wi
180     continue
190     w(i) = wi
200     j1 = j1 + lenr(i)
210    continue
       goto 240
c..
c..deals with rest of block which is singular.
220    do 230 ii=i,ilast
        resid = max(resid,abs(w(ii)))
        w(ii) = zero
230    continue
c..
c..back substitution. this loop does the back substitution on the rows of the 
c..block in the reverse order doing it simultaneously on the l transpose part
c..of the diagonal blocks and the off-diagonal blocks.
c..j1 points to the beginning of row i.
c..j2 points to the end of the l transpose part of row i.
240    j1 = iblend
       do 280 iback=ifirst,ilast
        i = ilast - iback + ifirst
        j1 = j1 - lenr(i)
        if (lenrl(i).eq.0) goto 260
        j2 = j1 + lenrl(i) - 1
        do 250 jj=j1,j2
         j = icn(jj)
         w(j) = w(j) + a(jj)*w(i)
250     continue
260     if (nobloc) goto 280
c..
c..operations using lower triangular blocks.
c..lj2 points to the end of row i of the off-diagonal blocks.
c..lj1 points to the beginning of row i of the off-diagonal blocks.
        if (lenoff(i).eq.0) goto 280
        lj2 = lj1 - 1
        lj1 = lj1 - lenoff(i)
        do 270 jj=lj1,lj2
         j = icn(jj)
         w(j) = w(j) - a(jj)*w(i)
270     continue
280    continue
       iblend = j1
       ilast = ifirst - 1
290   continue
c..
c..reorder solution vector ..  x(i)=w(ipinverse(i))
300   do 310 ii=1,n
       i = ip(ii)
       i = iabs(i)
       x(i) = w(ii)
310   continue
320   return
      end
c..
c..
c..
c..
c..
c..
      subroutine ma30dd(a,icn,iptr,n,iactiv,itop,reals)
      implicit double precision (a-h,o-z)
      save
c..
c..this subroutine performs garbage collection operations on the arrays a, 
c..icn and irn. iactiv is the first position in arrays a/icn from which the 
c..compress starts.  on exit, iactiv equals the position of the first entry
c..in the compressed part of a/icn
c..
      logical          reals
      integer          n,itop,iptr(n),icn(itop),
     1                 irncp,icncp,irank,minirn,minicn,j,k,kn,kl,
     2                 jpos,iactiv
      double precision a(itop)
      common /ma30fd/  irncp,icncp,irank,minirn,minicn
c..
c..
      if (reals) icncp = icncp + 1
      if (.not.reals) irncp = irncp + 1
c..
c..set the first non-zero entry in each row to the negative of the
c..row/col number and hold this row/col index in the row/col
c..pointer.  this is so that the beginning of each row/col can
c..be recognized in the subsequent scan.
      do 10 j=1,n
       k = iptr(j)
       if (k.lt.iactiv) goto 10
       iptr(j) = icn(k)
       icn(k) = -j
10    continue
      kn = itop + 1
      kl = itop - iactiv + 1
c..
c..go through arrays in reverse order compressing to the back so
c..that there are no zeros held in positions iactiv to itop in icn.
c..reset first entry of each row/col and pointer array iptr.
      do 30 k=1,kl
       jpos = itop - k + 1
       if (icn(jpos).eq.0) goto 30
       kn = kn - 1
       if (reals) a(kn) = a(jpos)
       if (icn(jpos).ge.0) goto 20
c..
c..first non-zero of row/col has been located
       j = -icn(jpos)
       icn(jpos) = iptr(j)
       iptr(j) = kn
20     icn(kn) = icn(jpos)
30    continue
      iactiv = kn
      return
      end
c..
c..
c..
c..
c..
      subroutine ma28int1
      implicit double precision (a-h,o-z)
      save
c..  
c..lp,mp are used by the subroutine as the unit numbers for its warning
c..      and diagnostic messages. default value for both is 6 (for line
c..      printer output). the user can either reset them to a different
c..      stream number or suppress the output by setting them to zero.
c..      while lp directs the output of error diagnostics from the
c..      principal subroutines and internally called subroutines, mp
c..      controls only the output of a message which warns the user that he
c..      has input two or more non-zeros a(i), .   ,a(k) with the same row
c..      and column indices.  the action taken in this case is to proceed
c..      using a numerical value of a(i)+.. +a(k). in the absence of other
c..      errors, iflag will equal -14 on exit.
c..lblock is a logical variable which controls an option of first
c..       preordering the matrix to block lower triangular form (using
c..       harwell subroutine mc23a). the preordering is performed if lblock
c..       is equal to its default value of .true. if lblock is set to
c..       .false. , the option is not invoked and the space allocated to
c..       ikeep can be reduced to 4*n+1.
c..grow is a logical variable. if it is left at its default value of
c..     .true. , then on return from ma28a/ad or ma28b/bd, w(1) will give
c..     an estimate (an upper bound) of the increase in size of elements
c..     encountered during the decomposition. if the matrix is well
c..     scaled, then a high value for w(1), relative to the largest entry
c..     in the input matrix, indicates that the lu decomposition may be
c..     inaccurate and the user should be wary of his results and perhaps
c..     increase u for subsequent runs.  we would like to emphasise that
c..     this value only relates to the accuracy of our lu decomposition
c..     and gives no indication as to the singularity of the matrix or the
c..     accuracy of the solution.  this upper bound can be a significant
c..     overestimate particularly if the matrix is badly scaled. if an
c..     accurate value for the growth is required, lbig (q.v.) should be
c..     set to .true.
c..eps,rmin are real variables. if, on entry to ma28b/bd, eps is less
c..     than one, then rmin will give the smallest ratio of the pivot to
c..     the largest element in the corresponding row of the upper
c..     triangular factor thus monitoring the stability of successive
c..     factorizations. if rmin becomes very large and w(1) from
c..     ma28b/bd is also very large, it may be advisable to perform a
c..     new decomposition using ma28a/ad.
c..resid is a real variable which on exit from ma28c/cd gives the value
c..      of the maximum residual over all the equations unsatisfied because
c..      of dependency (zero pivots).
c..irncp,icncp are integer variables which monitor the adequacy of "elbow
c..     room" in irn and a/icn respectively. if either is quite large (say
c..     greater than n/10), it will probably pay to increase the size of
c..     the corresponding array for subsequent runs. if either is very low
c..     or zero then one can perhaps save storage by reducing the size of
c..     the corresponding array.
c..minirn,minicn are integer variables which, in the event of a
c..     successful return (iflag ge 0 or iflag=-14) give the minimum size
c..     of irn and a/icn respectively which would enable a successful run
c..     on an identical matrix. on an exit with iflag equal to -5, minicn
c..     gives the minimum value of icn for success on subsequent runs on
c..     an identical matrix. in the event of failure with iflag= -6, -4,
c..     -3, -2, or -1, then minicn and minirn give the minimum value of
c..     licn and lirn respectively which would be required for a
c..     successful decomposition up to the point at which the failure occurred.
c..irank is an integer variable which gives an upper bound on the rank of
c..      the matrix.
c..abort1 is a logical variable with default value .true.  if abort1 is
c..       set to .false.  then ma28a/ad will decompose structurally singular
c..       matrices (including rectangular ones).
c..abort2 is a logical variable with default value .true.  if abort2 is
c..       set to .false. then ma28a/ad will decompose numerically singular
c..       matrices.
c..idisp is an integer array of length 2. on output from ma28a/ad, the
c..      indices of the diagonal blocks of the factors lie in positions
c..      idisp(1) to idisp(2) of a/icn. this array must be preserved
c..      between a call to ma28a/ad and subsequent calls to ma28b/bd,
c..      ma28c/cd or ma28i/id.
c..tol is a real variable.  if it is set to a positive value, then any
c..    non-zero whose modulus is less than tol will be dropped from the
c..    factorization.  the factorization will then require less storage
c..    but will be inaccurate.  after a run of ma28a/ad with tol positive
c..    it is not possible to use ma28b/bd and the user is recommended to
c..    use ma28i/id to obtain the solution.  the default value for tol is 0.0.
c..themax is a real variable.  on exit from ma28a/ad, it will hold the
c..       largest entry of the original matrix.
c..big is a real variable. if lbig has been set to .true., big will hold
c..    the largest entry encountered during the factorization by ma28a/ad
c..     or ma28b/bd.
c..dxmax is a real variable. on exit from ma28i/id, dxmax will be set to
c..      the largest component of the solution.
c..errmax is a real variable.  on exit from ma28i/id, if maxit is
c..       positive, errmax will be set to the largest component in the
c..       estimate of the error.
c..dres is a real variable.  on exit from ma28i/id, if maxit is positive,
c..     dres will be set to the largest component of the residual.
c..cgce is a real variable. it is used by ma28i/id to check the
c..     convergence rate.  if the ratio of successive corrections is
c..     not less than cgce then we terminate since the convergence
c..     rate is adjudged too slow.
c..ndrop is an integer variable. if tol has been set positive, on exit
c..     from ma28a/ad, ndrop will hold the number of entries dropped from
c..     the data structure.
c..maxit is an integer variable. it is the maximum number of iterations
c..     performed by ma28i/id. it has a default value of 16.
c..noiter is an integer variable. it is set by ma28i/id to the number of
c..     iterative refinement iterations actually used.
c..nsrch is an integer variable. if nsrch is set to a value less than n,
c..     then a different pivot option will be employed by ma28a/ad.  this
c..     may result in different fill-in and execution time for ma28a/ad.
c..     if nsrch is less than or equal to n, the workspace array iw can be
c..     reduced in length.  the default value for nsrch is 32768.
c..istart is an integer variable. if istart is set to a value other than
c..     zero, then the user must supply an estimate of the solution to
c..     ma28i/id.  the default value for istart is zero.
c..lbig is a logical variable. if lbig is set to .true., the value of the
c..    largest element encountered in the factorization by ma28a/ad or
c..    ma28b/bd is returned in big.  setting lbig to .true.  will
c..    increase the time for ma28a/ad marginally and that for ma28b/bd
c..    by about 20%.  the default value for lbig is .false.
c..
c..declare
      logical          lblock,grow,abort1,abort2,lbig
      integer          lp,mp,irncp,icncp,minirn,minicn,irank,
     1                 ndrop,maxit,noiter,nsrch,istart
      double precision eps,rmin,resid,tol,themax,big,dxmax,
     1                 errmax,dres,cgce
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn,
     1                 irank,abort1,abort2
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce,
     1                 ndrop,maxit,noiter,nsrch,istart,lbig
c..
      eps    = 1.0d-4
      tol    = 0.0d0
      cgce   = 0.5d0
      maxit  = 16
      lp     = 6
      mp     = 6
      nsrch  = 32768
      istart = 0
      lblock = .true.
      grow   = .false.
      lbig   = .false.
      abort1 = .true.
      abort2 = .true.
      return
      end
c..
c..
c..
c..
c..
      subroutine ma28int2
      implicit double precision (a-h,o-z)
      save
c..
c..common block ma30e/ed holds control parameters
c..    common /ma30ed/ lp, abort1, abort2, abort3
c..the integer lp is the unit number to which the error messages are
c..sent. lp has a default value of 6.  this default value can be
c..reset by the user, if desired.  a value of 0 suppresses all
c..messages.
c..the logical variables abort1,abort2,abort3 are used to control the
c..conditions under which the subroutine will terminate.
c..if abort1 is .true. then the subroutine will exit  immediately on
c..detecting structural singularity.
c..if abort2 is .true. then the subroutine will exit immediately on
c..detecting numerical singularity.
c..if abort3 is .true. then the subroutine will exit immediately when
c..the available space in a/icn is filled up by the previously decomposed, 
c..active, and undecomposed parts of the matrix.
c..
c..the default values for abort1,abort2,abort3 are set to .true.,.true.
c..and .false. respectively.
c..
c..
c..the variables in the common block ma30f/fd are used to provide the
c..user with information on the decomposition.
c..common /ma30fd/ irncp, icncp, irank, minirn, minicn
c..
c..irncp and icncp are integer variables used to monitor the adequacy
c..of the allocated space in arrays irn and a/icn respectively, by
c..taking account of the number of data management compresses
c..required on these arrays. if irncp or icncp is fairly large (say
c..greater than n/10), it may be advantageous to increase the size
c..of the corresponding array(s).  irncp and icncp are initialized
c..to zero on entry to ma30a/ad and are incremented each time the
c..compressing routine ma30d/dd is entered.
c..
c..icncp is the number of compresses on a/icn.
c..irncp is the number of compresses on irn.
c..
c..irank is an integer variable which gives an estimate (actually an
c..upper bound) of the rank of the matrix. on an exit with iflag
c..equal to 0, this will be equal to n.
c..
c minirn is an integer variable which, after a successful call to
c..ma30a/ad, indicates the minimum length to which irn can be
c..reduced while still permitting a successful decomposition of the
c..same matrix. if, however, the user were to decrease the length
c..of irn to that size, the number of compresses (irncp) may be
c..very high and quite costly. if lirn is not large enough to begin
c..the decomposition on a diagonal block, minirn will be equal to
c..the value required to continue the decomposition and iflag will
c..be set to -3 or -6. a value of lirn slightly greater than this
c..(say about n/2) will usually provide enough space to complete
c..the decomposition on that block. in the event of any other
c..failure minirn gives the minimum size of irn required for a
c..successful decomposition up to that point.
c..
c..minicn is an integer variable which after a successful call to
c..ma30a/ad, indicates the minimum size of licn required to enable
c..a successful decomposition. in the event of failure with iflag=
c..-5, minicn will, if abort3 is left set to .false., indicate the
c..minimum length that would be sufficient to prevent this error in
c..a subsequent run on an identical matrix. again the user may
c..prefer to use a value of icn slightly greater than minicn for
c..subsequent runs to avoid too many conpresses (icncp). in the
c..event of failure with iflag equal to any negative value except
c..-4, minicn will give the minimum length to which licn could be
c..reduced to enable a successful decomposition to the point at
c..which failure occurred.  notice that, on a successful entry
c..idisp(2) gives the amount of space in a/icn required for the
c..decomposition while minicn will usually be slightly greater
c..because of the need for "elbow room".  if the user is very
c..unsure how large to make licn, the variable minicn can be used
c..to provide that information. a preliminary run should be
c..performed with abort3 left set to .false. and licn about 3/2
c..times as big as the number of non-zeros in the original matrix.
c..unless the initial problem is very sparse (when the run will be
c..successful) or fills in extremely badly (giving an error return
c..with iflag equal to -4), an error return with iflag equal to -5
c..should result and minicn will give the amount of space required
c..for a successful decomposition.
c..
c..
c..common block ma30g/gd is used by the ma30b/bd entry only.
c..   common /ma30gd/ eps, rmin
c eps is a real/double precision variable. it is used to test for
c..small pivots. its default value is 1.0e-4 (1.0d-4 in d version).
c..if the user sets eps to any value greater than 1.0, then no
c..check is made on the size of the pivots. although the absence of
c..such a check would fail to warn the user of bad instability, its
c..absence will enable ma30b/bd to run slightly faster. an  a
c..posteriori  check on the stability of the factorization can be
c..obtained from mc24a/ad.
c..
c..rmin is a real/double precision variable which gives the user some
c..information about the stability of the decomposition.  at each
c..stage of the lu decomposition the magnitude of the pivot apiv
c..is compared with the largest off-diagonal entry currently in its
c..row (row of u), rowmax say. if the ratio min (apiv/rowmax)
c..where the minimum is taken over all the rows, is less than eps
c..then rmin is set to this minimum value and iflag is returned
c..with the value +i where i is the row in which this minimum
c..occurs.  if the user sets eps greater than one, then this test
c..is not performed. in this case, and when there are no small
c..pivots rmin will be set equal to eps.
c..
c..
c..common block ma30h/hd is used by ma30c/cd only.
c..   common /ma30hd/ resid
c..resid is a real/double precision variable. in the case of singular
c..or rectangular matrices its final value will be equal to the
c..maximum residual for the unsatisfied equations; otherwise its
c..value will be set to zero.
c..
c..
c..common  block ma30i/id controls the use of drop tolerances, the
c..modified pivot option and the the calculation of the largest
c..entry in the factorization process. this common block was added
c..to the ma30 package in february, 1983.
c..   common /ma30id/ tol, big, ndrop, nsrch, lbig
c..
c..tol is a real/double precision variable.  if it is set to a positive
c..value, then ma30a/ad will drop from the factors any non-zero
c..whose modulus is less than tol.  the factorization will then
c..require less storage but will be inaccurate.  after a run of
c..ma30a/ad where entries have been dropped, ma30b/bd  should not
c..be called.  the default value for tol is 0.0.
c..
c..big is a real/double precision variable.  if lbig has been set to
c.. true., big will be set to the largest entry encountered during
c..the factorization.
c..ndrop is an integer variable. if tol has been set positive, on exit
c..from ma30a/ad, ndrop will hold the number of entries dropped
c..from the data structure.
c..
c..nsrch is an integer variable. if nsrch is set to a value less than
c..or equal to n, then a different pivot option will be employed by
c..ma30a/ad.  this may result in different fill-in and execution
c..time for ma30a/ad. if nsrch is less than or equal to n, the
c..workspace arrays lastc and nextc are not referenced by ma30a/ad.
c..the default value for nsrch is 32768.
c..lbig is a logical variable. if lbig is set to .true., the value of
c..the largest entry encountered in the factorization by ma30a/ad
c..is returned in big.  setting lbig to .true.  will marginally
c..increase the factorization time for ma30a/ad and will increase
c..that for ma30b/bd by about 20%.  the default value for lbig is
c.. false.
c..
c..declare
      logical          abort1,abort2,abort3,lbig
      integer          lp,ndrop,nsrch  
      double precision eps,rmin,tol,big
c..
      common /ma30ed/ lp,abort1,abort2,abort3
      common /ma30gd/ eps,rmin
      common /ma30id/ tol,big,ndrop,nsrch,lbig
c..
      eps    = 1.0d-4
      tol    = 0.0d0
      big    = 0.0d0
      lp     = 6
      nsrch  = 32768
      lbig   = .false.
      abort1 = .true.
      abort2 = .true.
      abort3 = .false.
      return
      end
c..
c..
c..
c..
c..
      subroutine ma28int3
      logical         abort
      integer         lp,numnz,num,large
      common /mc23bd/ lp,numnz,num,large,abort
      lp       = 6
      abort    = .false.
      return
      end
c..
c..
c..
c..
c..
      subroutine mc20ad(nc,maxa,a,inum,jptr,jnum,jdisp)
      implicit double precision (a-h,o-z)
      save
c..
c..sorts a matrix into row order
c..
      integer          nc,maxa,inum(maxa),jnum(maxa),jptr(nc),jdisp,
     1                 null,j,k,kr,ice,jce,ja,jb,i,loc,icep,jcep
      double precision a(maxa),ace,acep
c..
c..go
      null=-jdisp
      do 60 j=1,nc
       jptr(j)=0
60    continue
c..
c..count the number of elements in each column.
      do 120 k=1,maxa
       j=jnum(k)+jdisp
       jptr(j)=jptr(j)+1
120   continue
c..
c..set the jptr array
      k=1
      do 150 j=1,nc
       kr=k+jptr(j)
       jptr(j)=k
       k=kr
150   continue
c..
c..reorder the elements into column order; an in-place sort of order maxa.
c.. jce is the current entry.
      do 230 i=1,maxa
       jce=jnum(i)+jdisp
       if(jce.eq.0) goto 230
       ace=a(i)
       ice=inum(i)
c..
c..clear the location vacated.
       jnum(i)=null
c..
c..chain from current entry to store items.
       do 200 j=1,maxa
        loc=jptr(jce)
        jptr(jce)=jptr(jce)+1
        acep=a(loc)
        icep=inum(loc)
        jcep=jnum(loc)
        a(loc)=ace
        inum(loc)=ice
        jnum(loc)=null
        if(jcep.eq.null) goto 230
        ace=acep
        ice=icep
        jce=jcep+jdisp
200    continue
230   continue
c..
c..reset jptr vector.
      ja = 1
      do 250 j=1,nc
       jb      = jptr(j)
       jptr(j) = ja
       ja      = jb
250   continue
      return
      end
c..
c..
c..
c..
c..
      subroutine mc23ad(n,icn,a,licn,lenr,idisp,ip,iq,lenoff,iw,iw1)
      implicit double precision (a-h,o-z)
      save
c..
c..performs the block triangularization
c..
c..declare
      logical          abort
      integer          n,licn,idisp(2),iw1(n,2),icn(licn),lenr(n),
     1                 ip(n),iq(n),lenoff(n),iw(n,5),lp,numnz,num,
     2                 large,i,ii,ibeg,iend,i1,i2,k,iblock,jnpos,
     3                 ilend,inew,irowe,irowb,leni,nz,j,jj,iold,jold,
     4                 jnew
      double precision a(licn)
      common /mc23bd/  lp,numnz,num,large,abort
c..
c..formats
180   format(1x,'matrix is structurally singular, rank = ',i6)
200   format(1x,'licn not big enough increase by ',i6)
220   format(1x,'error return from mc23ad because')
c..
c..set pointers iw(.,1) to beginning of the rows and set lenoff equal to lenr.
      iw1(1,1)=1
      lenoff(1)=lenr(1)
      if (n.eq.1) goto 20
      do 10 i=2,n
       lenoff(i)=lenr(i)
       iw1(i,1)=iw1(i-1,1)+lenr(i-1)
10    continue
c..
c..idisp(1) points to the first position in a/icn after the off-diagonal blocks 
c..and untreated rows.
20    idisp(1)=iw1(n,1)+lenr(n)
c..
c..find row permutation ip to make diagonal zero-free.
      call mc21a(n,icn,licn,iw1,lenr,ip,numnz,iw)
c..
c..possible error return for structurally singular matrices.
      if (numnz.ne.n.and.abort) goto 170
c..
c..iw1(.,2) and lenr are permutations of iw1(.,1) and lenr/lenoff suitable for 
c..entry to mc13d since matrix with these row pointer and length arrays has 
c..maximum number of non-zeros on the diagonal.
      do 30 ii=1,n
       i=ip(ii)
       iw1(ii,2)=iw1(i,1)
       lenr(ii)=lenoff(i)
30    continue
c..
c..find symmetric permutation iq to block lower triangular form.
      call mc13d(n,icn,licn,iw1(1,2),lenr,iq,iw(1,4),num,iw)
      if (num.ne.1) goto 60
c..
c..action taken if matrix is irreducible. whole matrix is just moved to the 
c..end of the storage.
      do 40 i=1,n
       lenr(i)=lenoff(i)
       ip(i)=i
       iq(i)=i
40    continue
      lenoff(1)=-1
c..
c..idisp(1) is the first position after the last element in the off-diagonal 
c..blocks and untreated rows.
      nz=idisp(1)-1
      idisp(1)=1
c..
c..idisp(2) is position in a/icn of the first element in the diagonal blocks.
      idisp(2)=licn-nz+1
      large=n
      if (nz.eq.licn) goto 230
      do 50 k=1,nz
       j=nz-k+1
       jj=licn-k+1
       a(jj)=a(j)
       icn(jj)=icn(j)
50    continue
      goto 230
c..
c..data structure reordered. form composite row permutation:ip(i) = ip(iq(i)).
60    do 70 ii=1,n
       i=iq(ii)
       iw(ii,1)=ip(i)
70    continue
      do 80 i=1,n
       ip(i)=iw(i,1)
80    continue
c..
c..run through blocks in reverse order separating diagonal blocks which are 
c..moved to the end of the storage.  elements in off-diagonal blocks are left 
c..in place unless a compress is necessary.
c..ibeg indicates the lowest value of j for which icn(j) has been
c..     set to zero when element in position j was moved to the
c..     diagonal block part of storage.
c..iend is position of first element of those treated rows which are in 
c..     diagonal blocks.
c..large is the dimension of the largest block encountered so far.
c..num is the number of diagonal blocks.
c..i1 is first row (in permuted form) of block iblock.
c..i2 is last row (in permuted form) of block iblock.
      ibeg=licn+1
      iend=licn+1
      large=0
      do 150 k=1,num
       iblock=num-k+1
       i1=iw(iblock,4)
       i2=n
       if (k.ne.1) i2=iw(iblock+1,4)-1
       large=max0(large,i2-i1+1)
c..
c..go through the rows of block iblock in the reverse order.
       do 140 ii=i1,i2
        inew=i2-ii+i1
c..
c..we now deal with row inew in permuted form (row iold in original matrix).
        iold=ip(inew)
c..
c..if there is space to move up diagonal block portion of row goto 110
        if (iend-idisp(1).ge.lenoff(iold)) goto 110
c..
c..in-line compress.; moves separated off-diagonal elements and untreated rows 
c..to front of storage.
        jnpos=ibeg
        ilend=idisp(1)-1
        if (ilend.lt.ibeg) goto 190
        do 90 j=ibeg,ilend
         if (icn(j).eq.0) goto 90
         icn(jnpos)=icn(j)
         a(jnpos)=a(j)
         jnpos=jnpos+1
90      continue
        idisp(1)=jnpos
        if (iend-jnpos.lt.lenoff(iold)) goto 190
        ibeg=licn+1
c..
c..reset pointers to the beginning of the rows.
        do 100 i=2,n
         iw1(i,1)=iw1(i-1,1)+lenoff(i-1)
100     continue
c..
c..row iold is now split into diag. and off-diag. parts.
110     irowb=iw1(iold,1)
        leni=0
        irowe=irowb+lenoff(iold)-1
c..
c..backward scan of whole of row iold (in original matrix).
        if (irowe.lt.irowb) goto 130
        do 120 jj=irowb,irowe
         j=irowe-jj+irowb
         jold=icn(j)
c..
c..iw(.,2) holds the inverse permutation to iq.; it was set to this in mc13d.
         jnew=iw(jold,2)
c..
c..if (jnew.lt.i1) then element is in off-diagonal block and so is left in situ.
         if (jnew.lt.i1) goto 120
c..
c..element is in diagonal block and is moved to the end of the storage.
         iend=iend-1
         a(iend)=a(j)
         icn(iend)=jnew
         ibeg=min0(ibeg,j)
         icn(j)=0
         leni=leni+1
120     continue
        lenoff(iold)=lenoff(iold)-leni
130     lenr(inew)=leni
140    continue
       ip(i2)=-ip(i2)
150   continue
c..
c..resets ip(n) to positive value.
c..idisp(2) is position of first element in diagonal blocks.
      ip(n)=-ip(n)
      idisp(2)=iend
c..
c..this compress used to move all off-diagonal elements to the front of storage.
      if (ibeg.gt.licn) goto 230
      jnpos=ibeg
      ilend=idisp(1)-1
      do 160 j=ibeg,ilend
       if (icn(j).eq.0) goto 160
       icn(jnpos)=icn(j)
       a(jnpos)=a(j)
       jnpos=jnpos+1
160   continue
c..
c..idisp(1) is first position after last element of off-diagonal blocks.
      idisp(1)=jnpos
      goto 230
c..
c..error return
170   if (lp.ne.0) write(lp,180) numnz
      idisp(1)=-1
      goto 210
190   if (lp.ne.0) write(lp,200) n
      idisp(1)=-2
210   if (lp.ne.0) write(lp,220)
230   return
      end
c..
c..
c..
c..
c..
      subroutine mc22ad(n,icn,a,nz,lenrow,ip,iq,iw,iw1)
      implicit double precision (a-h,o-z)
      save
c..
c..reorders the off diagonal blocks based on the pivot information
c..
c..declare
      integer          n,nz,iw(n,2),icn(nz),lenrow(n),ip(n),iq(n),
     1                 iw1(nz),i,jj,iold,j2,length,j,ipos,jval,
     2                 ichain,newpos,jnum
      double precision a(nz),aval
c..
c..go
      if (nz.le.0) goto 1000
      if (n.le.0) goto 1000
c..
c..set start of row i in iw(i,1) and lenrow(i) in iw(i,2)
      iw(1,1)=1
      iw(1,2)=lenrow(1)
      do 10 i=2,n
       iw(i,1)=iw(i-1,1)+lenrow(i-1)
       iw(i,2)=lenrow(i)
10    continue
c..
c..permute lenrow according to ip.  set off-sets for new position of row iold 
c..in iw(iold,1) and put old row indices in iw1 in positions corresponding to 
c..the new position of this row in a/icn.
      jj=1
      do 20 i=1,n
       iold=ip(i)
       iold=iabs(iold)
       length=iw(iold,2)
       lenrow(i)=length
       if (length.eq.0) goto 20
       iw(iold,1)=iw(iold,1)-jj
       j2=jj+length-1
       do 15 j=jj,j2
        iw1(j)=iold
15     continue
       jj=j2+1
20    continue
c..
c..set inverse permutation to iq in iw(.,2).
      do 30 i=1,n
       iold=iq(i)
       iold=iabs(iold)
       iw(iold,2)=i
30    continue
c..
c..permute a and icn in place, changing to new column numbers.
c..main loop; each pass through this loop places a closed chain of column 
c..indices in their new (and final) positions ..  this is recorded by
c..setting the iw1 entry to zero so that any which are subsequently
c..encountered during this major scan can be bypassed.
      do 200 i=1,nz
       iold=iw1(i)
       if (iold.eq.0) goto 200
       ipos=i
       jval=icn(i)
c..
c..if row iold is in same positions after permutation goto 150.
       if (iw(iold,1).eq.0) goto 150
       aval=a(i)
c..
c..chain loop; each pass through this loop places one (permuted) column index
c..in its final position  .. viz. ipos.
c..newpos is the original position in a/icn of the element to be placed
c..in position ipos.  it is also the position of the next element in the chain.
       do 100 ichain=1,nz
        newpos=ipos+iw(iold,1)
        if (newpos.eq.i) goto 130
        a(ipos)=a(newpos)
        jnum=icn(newpos)
        icn(ipos)=iw(jnum,2)
        ipos=newpos
        iold=iw1(ipos)
        iw1(ipos)=0
100    continue
130    a(ipos)=aval
150    icn(ipos)=iw(jval,2)
200   continue
1000  return
      end
c..
c..
c..
c..
c..
      subroutine mc21a(n,icn,licn,ip,lenr,iperm,numnz,iw)
      implicit double precision (a-h,o-z)
      save
      integer n,licn,ip(n),icn(licn),lenr(n),iperm(n),iw(n,4),numnz
      call mc21b(n,icn,licn,ip,lenr,iperm,numnz,iw(1,1),iw(1,2),iw(1,3),
     1           iw(1,4))
      return
      end
c..
c..
c..
c..
c..
      subroutine mc21b(n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      implicit double precision (a-h,o-z)
      save
c..
c..does a row permutation to make the diagonal zero free
c..
c..pr(i) is the previous row to i in the depth first search.
c..     it is used as a work array in the sorting algorithm.
c..     elements (iperm(i),i) i=1, ..  n  are non-zero at the end of the
c..     algorithm unless n assignments have not been made.  in which case
c..(iperm(i),i) will be zero for n-numnz entries.
c..cv(i)  is the most recent row extension at which column i was visited.
c..arp(i) is one less than the number of non-zeros in row i
c..       which have not been scanned when looking for a cheap assignment.
c..out(i) is one less than the number of non-zeros in row i
c..       which have not been scanned during one pass through the main loop.
c..
c..declare
      integer n,licn,ip(n),icn(licn),lenr(n),iperm(n),pr(n),cv(n),
     1        arp(n),out(n),i,jord,j,in1,in2,k,ii,ioutk,j1,kk,numnz
c..
c..initialization of arrays.
      do 10 i=1,n
       arp(i)=lenr(i)-1
       cv(i)=0
       iperm(i)=0
10    continue
      numnz=0
c..
c..main loop. each pass round this loop either results in a new assignment
c..or gives a row with no assignment.
      do 130 jord=1,n
       j=jord
       pr(j)=-1
       do 100 k=1,jord
c..
c..look for a cheap assignment
        in1=arp(j)
        if (in1.lt.0) goto 60
        in2=ip(j)+lenr(j)-1
        in1=in2-in1
        do 50 ii=in1,in2
         i=icn(ii)
         if (iperm(i).eq.0) goto 110
50      continue
c..
c..no cheap assignment in row.
c..begin looking for assignment chain starting with row j.
        arp(j)=-1
60      out(j)=lenr(j)-1
c..
c..c inner loop.  extends chain by one or backtracks.
        do 90 kk=1,jord
         in1=out(j)
         if (in1.lt.0) goto 80
         in2=ip(j)+lenr(j)-1
         in1=in2-in1
c..
c..forward scan.
         do 70 ii=in1,in2
          i=icn(ii)
          if (cv(i).eq.jord) goto 70
c..
c..column i has not yet been accessed during this pass.
          j1=j
          j=iperm(i)
          cv(i)=jord
          pr(j)=j1
          out(j1)=in2-ii-1
          goto 100
70       continue
c..
c..backtracking step.
80       j=pr(j)
         if (j.eq.-1) goto 130
90      continue
100    continue
c..
c..new assignment is made.
110    iperm(i)=j
       arp(j)=in2-ii-1
       numnz=numnz+1
       do 120 k=1,jord
        j=pr(j)
        if (j.eq.-1) goto 130
        ii=ip(j)+lenr(j)-out(j)-2
        i=icn(ii)
        iperm(i)=j
120    continue
130   continue
c..
c..if matrix is structurally singular, we now complete the permutation iperm.
      if (numnz.eq.n) return
      do 140 i=1,n
       arp(i)=0
140   continue
      k=0
      do 160 i=1,n
       if (iperm(i).ne.0) goto 150
       k=k+1
       out(k)=i
       goto 160
150    j=iperm(i)
       arp(j)=i
160   continue
      k=0
      do 170 i=1,n
       if (arp(i).ne.0) goto 170
       k=k+1
       ioutk=out(k)
       iperm(ioutk)=i
170   continue
      return
      end
c..
c..
c..
c..
c..
      subroutine mc13d(n,icn,licn,ip,lenr,ior,ib,num,iw)
      implicit double precision (a-h,o-z)
      save
      integer n,licn,ip(n),icn(licn),lenr(n),ior(n),ib(n),iw(n,3),num
      call mc13e(n,icn,licn,ip,lenr,ior,ib,num,iw(1,1),iw(1,2),iw(1,3))
      return
      end
c..
c..
c..
c..
c..
      subroutine mc13e(n,icn,licn,ip,lenr,arp,ib,num,lowl,numb,prev)
      implicit double precision (a-h,o-z)
      save
c..
c.. arp(i) is one less than the number of unsearched edges leaving
c..        node i.  at the end of the algorithm it is set to a
c..        permutation which puts the matrix in block lower
c..        triangular form.
c..ib(i)   is the position in the ordering of the start of the ith
c..        block.  ib(n+1-i) holds the node number of the ith node
c..        on the stack.
c..lowl(i) is the smallest stack position of any node to which a path
c..        from node i has been found.  it is set to n+1 when node i
c..        is removed from the stack.
c..numb(i) is the position of node i in the stack if it is on
c..        it, is the permuted order of node i for those nodes
c..        whose final position has been found and is otherwise zero.
c..prev(i) is the node at the end of the path when node i was
c..        placed on the stack.
c..
c..declare
      integer n,licn,stp,dummy,ip(n),icn(licn),lenr(n),arp(n),ib(n),
     1        lowl(n),numb(n),prev(n),icnt,num,nnm1,j,iv,ist,i1,i2,
     2        ii,iw,ist1,lcnt,i,isn,k
c..
c..
c..icnt is number of nodes whose positions in final ordering have been found.
c..num is the number of blocks that have been found.
      icnt = 0
      num  = 0
      nnm1 = n + n-1
c..
c..initialization of arrays.
      do 20 j=1,n
       numb(j)=0
       arp(j)=lenr(j)-1
20    continue
c..
c..look for a starting node
c..ist is the number of nodes on the stack ..  it is the stack pointer.
      do 120 isn=1,n
       if (numb(isn).ne.0) goto 120
       iv=isn
       ist=1
c..
c..put node iv at beginning of stack.
       lowl(iv)=1
       numb(iv)=1
       ib(n)=iv
c..
c..the body of this loop puts a new node on the stack or backtracks.
       do 110 dummy=1,nnm1
        i1=arp(iv)
c..
c..have all edges leaving node iv been searched.
        if (i1.lt.0) goto 60
        i2=ip(iv)+lenr(iv)-1
        i1=i2-i1
c..
c..look at edges leaving node iv until one enters a new node or all edges are 
c..exhausted.
        do 50 ii=i1,i2
         iw=icn(ii)
         if (numb(iw).eq.0) goto 100
         lowl(iv)=min0(lowl(iv),lowl(iw))
50      continue
c..
c..there are no more edges leaving node iv.
        arp(iv)=-1
c..
c..is node iv the root of a block.
60      if (lowl(iv).lt.numb(iv)) goto 90
c..
c..order nodes in a block.
        num=num+1
        ist1=n+1-ist
        lcnt=icnt+1
c..
c..peel block off the top of the stack starting at the top and working down to 
c..the root of the block.
        do 70 stp=ist1,n
         iw=ib(stp)
         lowl(iw)=n+1
         icnt=icnt+1
         numb(iw)=icnt
         if (iw.eq.iv) goto 80
70      continue
80      ist=n-stp
        ib(num)=lcnt
c..
c..are there any nodes left on the stack.
        if (ist.ne.0) goto 90
c..
c..have all the nodes been ordered.
        if (icnt.lt.n) goto 120
        goto 130
c..
c..backtrack to previous node on path.
90      iw=iv
        iv=prev(iv)
c..
c..update value of lowl(iv) if necessary.
        lowl(iv)=min0(lowl(iv),lowl(iw))
        goto 110
c..
c..put new node on the stack.
100     arp(iv)=i2-ii-1
        prev(iw)=iv
        iv=iw
        ist=ist+1
        lowl(iv)=ist
        numb(iv)=ist
        k=n+1-ist
        ib(k)=iv
110    continue
120   continue
c..
c..put permutation in the required form.
130   do 140 i=1,n
       ii=numb(i)
       arp(ii)=i
140   continue
      return
      end
c..
c..
c..
c..
c..
      subroutine mc24ad(n,icn,a,licn,lenr,lenrl,w)
      implicit double precision (a-h,o-z)
      save
c..
c..computes the gwoth rate of fill in
c..
      integer          n,licn,icn(licn),lenr(n),lenrl(n),i,j0,j2,j1,jj,j
      double precision a(licn),w(n),amaxl,wrowl,amaxu,zero
      data             zero/0.0d0/
c..
c..initialize
      amaxl=zero
      do 10 i=1,n
       w(i)=zero
10    continue
      j0=1
      do 100 i=1,n
       if (lenr(i).eq.0) goto 100
       j2=j0+lenr(i)-1
       if (lenrl(i).eq.0) goto 50
c..
c..calculation of 1-norm of l.
       j1=j0+lenrl(i)-1
       wrowl=zero
       do 30 jj=j0,j1
        wrowl=wrowl+abs(a(jj))
30     continue
c..
c..amaxl is the maximum norm of columns of l so far found.
       amaxl=max(amaxl,wrowl)
       j0=j1+1
c..
c..calculation of norms of columns of u (max-norms).
50     j0=j0+1
       if (j0.gt.j2) goto 90
       do 80 jj=j0,j2
        j=icn(jj)
        w(j)=max(abs(a(jj)),w(j))
80     continue
90     j0=j2+1
100   continue
c..
c..amaxu is set to maximum max-norm of columns of u.
      amaxu=zero
      do 200 i=1,n
       amaxu=max(amaxu,w(i))
200   continue
c..
c..grofac is max u max-norm times max l 1-norm.
      w(1) = amaxl*amaxu
      return
      end
c..
c..
c..
c..
c..
      subroutine mc20bd(nc,maxa,a,inum,jptr)
      implicit double precision (a-h,o-z)
      save
c..
c..never called
c..
c..
c..declare
      integer          nc,maxa,inum(maxa),jptr(nc),kmax,jj,j,klo,kor,
     1                 kdummy,ice,k,ik
      double precision a(maxa),ace
c..
c..go
      kmax=maxa
      do 35 jj=1,nc
       j=nc+1-jj
       klo=jptr(j)+1
       if(klo.gt.kmax)goto 30
       kor=kmax
c..
c..items kor, kor+1, ...  ,kmax are in order
       do 25 kdummy=klo,kmax
        ace=a(kor-1)
        ice=inum(kor-1)
        do 10 k=kor,kmax
         ik=inum(k)
         if (iabs(ice).le.iabs(ik)) goto 20
         inum(k-1)=ik
         a(k-1)=a(k)
10      continue
        k=kmax+1
20      inum(k-1)=ice
        a(k-1)=ace
        kor=kor-1
25     continue
30     kmax=klo-2
35    continue
      return
      end

