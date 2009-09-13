      subroutine helmeos(do_coulomb, npts, eosfail,
     &                   temp_row, den_row, abar_row, zbar_row,
     &                   etot_row, ptot_row, 
     &                   cv_row, cp_row, xne_row, xnp_row, etaele_row,
     &                   pele_row, ppos_row, 
     &                   dpd_row, dpt_row, dpa_row, dpz_row, 
     &                   ded_row, det_row, dea_row, dez_row, 
     &                   gam1_row, cs_row, stot_row,
     &                   dsd_row, dst_row)
      use bl_error_module
      implicit none
      include 'vector_eos.dek'

c  Frank Timmes Helmholtz based Equation of State
c  http://cococubed.asu.edu/


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

c..MLW:  Modified to pass explicit args and remove save variables
c        allowing for multiple threads to call eos simultaniously
c..input arguments
      logical do_coulomb, eosfail
      integer npts
      double precision temp_row(npts), den_row(npts), abar_row(npts),
     1                 zbar_row(npts), etot_row(npts), ptot_row(npts),
     2                 cv_row(npts), cp_row(npts), 
     4                 xne_row(npts), xnp_row(npts), etaele_row(npts),
     5                 pele_row(npts), ppos_row(npts), dpd_row(npts), 
     6                 dpt_row(npts), dpa_row(npts), dpz_row(npts), 
     7                 ded_row(npts), det_row(npts), dea_row(npts), 
     8                 dez_row(npts), 
     9                 stot_row(npts), dsd_row(npts), dst_row(npts)
      double precision gam1_row(npts), cs_row(npts)

c..declare some parameters
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


c..declare local variables
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
     9                 temp,den,abar,zbar,ytot1,ye

      double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
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



c..for the interpolations
      integer          iat,jat
      double precision tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
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




c..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     2                 ecoul,decouldd,decouldt,decoulda,decouldz,
     3                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     4                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     5                 tmelt,tfermi,rhocond,z2,x1,x2

      double precision  a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
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


c ======================================================================
c..Define Statement Functions

c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
c      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      psi0(z)   = z*z*z * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
c      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      dpsi0(z)  = z*z * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

c..psi1 and its derivatives
c      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      psi1(z)   = z* ( z*z * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
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



! Dont allow multi-streaming of this function on CRAY X1
!DIR$  NOSTREAM

c..popular format statements
01    format(1x,5(a,1pe11.3))

c..start of vectorization loop, normal executaion starts here
      eosfail = .false.
      do j =1, npts

         temp  = temp_row(j)
         den   = den_row(j)
         abar  = abar_row(j)
         zbar  = zbar_row(j)
         ytot1 = 1.0d0/abar
         ye    = ytot1 * zbar
         din = ye*den

         if (temp .le. 0.0D0) then
            call bl_error('EOS: temp less than 0')
         end if
         if (den  .le. 0.0D0) then
            call bl_error('EOS: den less than 0')
         end if

         if ( temp .lt. t(1) ) then
            print *, 'TEMP = ', temp
            call bl_warn("EOS: temp too cold")
            eosfail = .true.
            return
         end if
         if ( temp .gt. t(jmax) ) then
            print *, 'TEMP = ', temp
            call bl_warn("EOS: temp too hot")
            eosfail = .true.
            return
         end if
         if ( din .gt. d(imax) ) then
            call bl_warn("EOS: ye*den out of bounds")
            eosfail = .true.
            return
         end if
         if ( din .lt. d(1) ) then
            call bl_warn("EOS: ye*den out of bounds")
            eosfail = .true.
            return
         end if
      enddo

      do j = 1, npts

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
c       if (temp .gt. t(jmax)) then
c        write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
c        write(6,*) 'temp too hot, off grid'       
c        write(6,*) 'setting eosfail to true and returning'
cc       call flush(6)
c        eosfail = .true.
c        return
c       end if
c       if (temp .lt. t(1)) then
c        write(6,01) 'temp=',temp,' t(1)=',t(1)
c        write(6,*) 'temp too cold, off grid'
c        write(6,*) 'setting eosfail to true and returning'
cc       call flush(6)
c        eosfail = .true.
c        return
c       end if
c       if (din  .gt. d(imax)) then
c        write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
c        write(6,*) 'ye*den too big, off grid'
c        write(6,*) 'setting eosfail to true and returning'
cc       call flush(6)
c        eosfail = .true.
c        return
c       end if
c       if (din  .lt. d(1)) then
c        write(6,01) 'ye*den=',din,' d(1)=',d(1)
c        write(6,*) 'ye*den too small, off grid'
c        write(6,*) 'setting eosfail to true and returning'
cc       call flush(6)
c        eosfail = .true.
c        return
c       end if

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



!     TURN ON/OFF COULOMB
        IF ( do_coulomb ) THEN
cc..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
           if (plasg .ge. 1.0D0) then
              x        = plasg**(0.25d0) 
              y        = avo * ytot1 * kerg 
              ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
              pcoul    = third * den * ecoul
              scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1             + d1 * (log(plasg) - 1.0d0) - e1)
            
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
            
              y        = -avo*kerg/(abar*plasg)*
     1             (0.75d0*b1*x+1.25d0*c1/x+d1)
              dscouldd = y * plasgdd
              dscouldt = y * plasgdt
              dscoulda = y * plasgda - scoul/abar
              dscouldz = y * plasgdz
          
c     c..yakovlev & shalybkov 1989 equations 102, 103, 104
           else if (plasg .lt. 1.0D0) then
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
              
              s        = -avo*kerg/(abar*plasg)*
     1             (1.5d0*c2*x-a2*(b2-1.0d0)*y)
              dscouldd = s * plasgdd
              dscouldt = s * plasgdt
              dscoulda = s * plasgda - scoul/abar
              dscouldz = s * plasgdz
           end if
c
cc..bomb proof
           x   = prad + pion + pele + pcoul
           if (x .le. 0.0D0) then

              write(6,*) 
              write(6,*) 'coulomb corr. are causing a negative pressure'
              write(6,*) 'setting all coulomb corrections to zero'
              write(6,*) 

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
        END IF
!     TURN OFF COULOMB

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
c..MLW: comment out storage to *_row arrays we dont use 
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
c       dsa_row(j)    = dentrda        
c       dsz_row(j)    = dentrdz

c       prad_row(j)   = prad
c       erad_row(j)   = erad
c       srad_row(j)   = srad 

c       pion_row(j)   = pion
c       eion_row(j)   = eion
c       sion_row(j)   = sion 
c       xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
c       dpept_row(j)  = dpepdt
c       dpepd_row(j)  = dpepdd
c       dpepa_row(j)  = dpepda  
c       dpepz_row(j)  = dpepdz

c       eele_row(j)   = eele
c       epos_row(j)   = 0.0d0
c       deept_row(j)  = deepdt
c       deepd_row(j)  = deepdd
c       deepa_row(j)  = deepda   
c       deepz_row(j)  = deepdz

c       sele_row(j)   = sele 
c       spos_row(j)   = 0.0d0
c       dsept_row(j)  = dsepdt 
c       dsepd_row(j)  = dsepdd 
c       dsepa_row(j)  = dsepda        
c       dsepz_row(j)  = dsepdz

c       xnem_row(j)   = xnem
        xne_row(j)    = xnefer
c       dxnet_row(j)  = dxnedt
c       dxned_row(j)  = dxnedd
c       dxnea_row(j)  = dxneda
c       dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0

        etaele_row(j) = etaele
c       detat_row(j)  = detadt
c       detad_row(j)  = detadd
c       detaa_row(j)  = detada
c       detaz_row(j)  = detadz
c       etapos_row(j) = 0.0d0

c       pcou_row(j)   = pcoul
c       ecou_row(j)   = ecoul
c       scou_row(j)   = scoul 
c       plasg_row(j)  = plasg

c       dse_row(j)    = dse
c       dpe_row(j)    = dpe
c       dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
c       gam2_row(j)   = gam2
c       gam3_row(j)   = gam3
        cs_row(j)     = sound

c..end of vectorization loop
      enddo
      return
      end

c     -----------------------------------------------------------------
      subroutine helmeos_init
      use bl_error_module
      use parallel
      implicit none
      
      include 'vector_eos.dek'

      double precision dth, dt2, dti, dt2i
      double precision dd, dd2, ddi, dd2i
      double precision tsav, dsav
      integer i, j

      logical          first
      data             first/.false./ 

c..   popular format statements
 03   format(1x,4(a,1pe11.3))
 04   format(1x,4(a,i4))

      if ( first ) return

      first = .true.
      
c..   open the table
      open(unit=2,file='helm_table.dat',status='old',err=1010)
      GO TO 1011
 1010 CONTINUE
      call bl_error('EOSINIT: Failed to open helm_table.dat')
 1011 CONTINUE

c..   read the helmholtz free energy table
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
     &           fddt(i,j),fdtt(i,j),fddtt(i,j)
         end do
      end do

c..   read the pressure derivative with density table
      do j = 1, jmax
         do i = 1, imax
            read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
         end do
      end do

c..   read the electron chemical potential table
      do j = 1, jmax
         do i = 1, imax
            read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
         end do
      end do

c..   read the number density table
      do j = 1, jmax
         do i = 1, imax
            read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
         end do
      end do

c..   construct the temperature and density deltas and their inverses 
      do j = 1, jmax-1
         dth         = t(j+1) - t(j)
         dt2         = dth * dth
         dti         = 1.0d0/dth
         dt2i        = 1.0d0/dt2
         dt_sav(j)   = dth
         dt2_sav(j)  = dt2
         dti_sav(j)  = dti
         dt2i_sav(j) = dt2i
      end do
      do i = 1, imax-1
         dd          = d(i+1) - d(i)
         dd2         = dd * dd
         ddi         = 1.0d0/dd
         dd2i        = 1.0d0/dd2
         dd_sav(i)   = dd
         dd2_sav(i)  = dd2
         ddi_sav(i)  = ddi
         dd2i_sav(i) = dd2i
      end do

      close(unit=2)
      if ( parallel_ioprocessor() ) then
         write(6,*)
         write(6,*) 'finished reading eos table'
         write(6,04) 'imax=',imax,' jmax=',jmax
         write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
         write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
         write(6,*)
      end if

      end

