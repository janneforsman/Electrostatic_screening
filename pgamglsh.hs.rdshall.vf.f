      program platem
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
      dimension c(maxmon,0:maxel),chA(0:maxel),chB(0:maxel)
c     Generalized van der Waals theory of an HS solution
c     Specifically designed to calculate the interfacial tension
c     between fluid bulk phases
c     Note: this is an unshifted version!

c     THIS IS A RELATIVELY UNTESTED PROGRAM! IT APPEARS TO
c     WORK FOR HS SYSTEMS (FINAL FREE ENERGY EXPRESSION NOT
c     CHECKED), BUT ONLY IF LAGRANGE MULTIPLIERS ARE NOT USED
c     (OR SET EQUAL TO ONE)!

      do 1 i = 0,maxel
      fdmon(i) = 0.d0
      fem(i) = 0.d0
      fdsol(i) = 0.d0
      cdmonm(i) = 0.d0
      cdsolm(i) = 0.d0
      cdsols(i) = 0.d0
      cdmons(i) = 0.d0
      ae0(i) = 0.d0
      ae1(i) = 0.d0
      ae2(i) = 0.d0
      convpm(i) = 0.d0
      convps(i) = 0.d0
      convsm(i) = 0.d0
      convss(i) = 0.d0
      Ymix(i) = 0.d0
      chA(i) = 0.d0
      chB(i) = 0.d0
      do 2 j = 1,maxmon
 2    c(j,i) = 0.d0
 1    continue
      ifc = 38
      ins = 49
      iefil = 54
      issolv = 56
      ialp = 60
      ibdm = 65
      ieplag = 85
      pi = acos(-1.d0)
      bk = 1.38066D-23
      avno = 6.02214D23
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      rphi = fourpi*0.2d0 
      aphi = fourpi*0.5d0
      pis = pi/6.d0
      pit = pi/3.d0
      pif = pi/4.d0
      es22 = -32.d0*pi/9.d0    
c      dFmaxok = 1.E-5
      ddtol = 0.000002d0
c     CLOSE TO THE WALLS, THE DENSITY IS ASSUMED TO BE ZERO
      open (ifc,file='fcdfil',form='formatted')    
      open (ins,file='input.dftp',form='formatted')
      open (iefil,file='ekblamfil',form='formatted')
      open (issolv,file='sigmasolv',form='formatted')
      open (ialp,file='alpha',form='formatted')
      open(ieplag,file='eplagfil',form='formatted')
      open(ibdm,file='bdmfil',form='formatted')
      rewind ifc
      rewind ins
      rewind iefil
      rewind issolv
      rewind ialp
      rewind ieplag
      rewind ibdm
      read(ins,*) refbdm
      read(ins,*) refbds
c     bds written in sigma(solvent) units
      read(ins,*) nmon
      read(ins,*) trams
      read(ins,*) ioimaxm
      read(ins,*) ntrams
      read(ins,*) T
      read(ins,*) h
      read(ins,*) dz
      read(ins,*) closew
      read(ins,*) dmm
      read(ins,*) dms
      read(ins,*) kread
      read(ins,*) hwdist
      read(ins,*) awm
      read(ins,*) rwm
      read(ins,*) aws
      read(ins,*) rws
      read(ins,*) ufact
      read(ins,*) ess
      read(ins,*) ems
      ufact = 0.d0
      ess = 0.d0
      ems = 0.d0
      read(ins,*) refbdist
      read(issolv,*) q1
c      read(issolv,*) sclosew
      read(ibdm,*) bdmg
      read(ibdm,*) bdsg
      read(ibdm,*) bdml
      read(ibdm,*) bdsl
      read(ialp,*) alpha
      read(ialp,*) bdist
      if (q1.gt.1.00001d0) then
      write(*,*) 'Q1 IS NOT ALLOWED TO BE GT. 1!'
      goto 9999
      endif
      q2 = q1*q1
      q3 = q2*q1
      q4 = q2*q2
      q10 = q4*q4*q2
      p1 = (1.d0+q1)*0.5d0
      p2 = p1*p1
      p3 = p2*p1
      p4 = p2*p2
      p10 = p4*p4*p2
      rq3 = 1.d0/q3
      rp3 = 1.d0/p3
      write(*,*) 'q1,p1 = ',q1,p1
c     transformation to sigma(mon) units
      bdsl = bdsl*rq3
      bdsg = bdsg*rq3
      tdmm = 1.d0-dmm
      tdms = 1.d0-dms
      rrT = 1.d0/T
      es22 = 0.d0
      es22mm = 0.d0
      es22ss = 0.d0
      es22ms = 0.d0
      phzmm = 0.d0
      phzss = 0.d0
      phzms = 0.d0
      aphi = 0.d0
      rphi = 0.d0
      rdz = 1.d0/dz
      rnmon = real(nmon)
      twopidz = twopi*dz
      dmax = bdist+rnmon
      irdz = int(rdz+0.001d0)
      ism = int(rdz+0.01d0)
      isms = int(q1/dz+0.01d0)
      ismms = int(p1/dz+0.01d0)

      istp1 = int((rnmon-1.d0)*rdz+0.01d0)+1
c      istp1 = int((rnmon+2.d0)*rdz+0.01d0)+1
      ibd = int(bdist*rdz+0.01d0)
      imitt = istp1+ibd+2*ism
c      imitt = istp1+ibd+3*ism

      dmitt = real(imitt)*dz-0.5d0*dz
      distp1 = real(istp1)*dz-0.5d0*dz

      dgibbs = (real(imitt-ism)*dz+rnmon)*0.5d0
      idg = int(dgibbs*rdz+0.01d0)
      write(*,*) 'istp1,imitt = ',istp1,imitt
      write(*,*) 'ibd,ism = ',ibd,ism
      write(*,*) 'dgibbs,idg = ',dgibbs,idg
      write(*,*) 'dmitt,distp1 = ',dmitt,distp1
      write(*,*) 'alpha = ',alpha
      bfeml = bdml*2.d0/rnmon
      bfemg = bdmg*2.d0/rnmon

      pie = pi/8.d0
      dzpie = pie*dz
      rnmon = real(nmon)
c      Yfact = (rnmon-2.d0)*Y
      ve1 = pis*(1.d0+q1)**3
      ve2 = pis*(2.d0+6.d0*q1+4.5d0*q2+q3)
      ve3 = 1.57d0+4.75d0*q1+2.99d0*q2+0.52d0*q3
      venmon= ve3+(rnmon-3.d0)*(ve3-ve2)-
     *0.04915d0*(rnmon-3.d0)**1.09d0*q1**2.71d0
      veq = (venmon-ve2)/(ve2-ve1)
      sve1 = pis*(1.d0+1.d0)**3
      sve2 = pis*(2.d0+6.d0+4.5d0+1.d0)
      sve3 = 1.57d0+4.75d0+2.99d0+0.52d0
      svenmon= sve3+(rnmon-3.d0)*(sve3-sve2)-
     *0.04915d0*(rnmon-3.d0)**1.09d0
      Yfact = (svenmon-sve2)/(sve2-sve1)
      rrnmon = 1.d0/rnmon      
      rrcmon = 1.d0/(rnmon-2.d0)
c      checknm = abs(rnmon*0.5d0-real(int(rnmon*0.5d0+0.0001)))
c      if (checknm.gt.0.00001) then
c      write(*,*) 'NMON MUST BE AN EVEN NUMBER!!!'
c      goto 9999
c      endif
c      nhalfmon = nmon/2

      r2 = 0.75d0
      r1 = 0.5d0
      r0 = 0.5d0*q1
      s2 = twopi
      s1 = pi
      s0 = q2*pi
      b2 = pit
      b1 = pis
      b0 = pis*q3
      r2sq = r2*r2
      r1sq = r1*r1
      r0sq = r0*r0

      bdm = bdml
      bds = bdsl
      bdt = bds*q3+bdm
      bdpol = bdm/rnmon
      cdmlbulk = bdm
      cdslbulk = bds*q3

      aeta = pis*bdt
      xsib = 1.d0-aeta
      rxsib = 1.d0/xsib
      rxsib2 = rxsib*rxsib
      rxsib3 = rxsib*rxsib2
      raeta = 1.d0/aeta
      aeta2 = aeta*aeta
      aeta3 = aeta2*aeta
      useful = -dlog(xsib)*raeta-xsib*dlog(xsib)/aeta2-raeta

      xs = bds/(bds+bdpol)
      xp = 1.d0-xs
      xd = 0.5d0*bdm/(0.5d0*bdm+bds)
      xm = bdm/(bdm+bds)
      xsd = 1.d0-xd
      xsm = 1.d0-xm
      hr2 = xd*r2+xsd*r0
      hr1 = xm*r1+xsm*r0
      hs2 = xd*s2+xsd*s0
      hs1 = xm*s1+xsm*s0
      hb2 = xd*b2+xsd*b0
      hb1 = xm*b1+xsm*b0
      hc2 = xd*r2sq+xsd*r0sq
      hc1 = xm*r1sq+xsm*r0sq
      rhb2 = 1.d0/hb2
      rhb2sq = rhb2*rhb2
      rnhb2sq = rhb2sq/9.d0
      hs2sq = hs2*hs2
      hb2sq = hb2*hb2
      rhb1 = 1.d0/hb1
      rhb1sq = rhb1*rhb1
      rnhb1sq = rhb1sq/9.d0      
      hs1sq = hs1*hs1
      hb1sq = hb1*hb1
      write(*,*) 'xd,xm = ',xd,xm

      t1ae = xsib*dlog(xsib)*raeta
      aeT1 = dlog(xsib)*(xsib*raeta+1.d0)-rxsib-0.5d0*rxsib2+2.5d0
      aeTXZ = 0.5d0*rxsib2-0.5d0
      aeT2 = (rxsib-0.5d0*rxsib2-0.5d0)
      aeW = 0.5d0*rxsib2-2.d0*rxsib-dlog(xsib)+1.5d0

      psi2 = hc2*hs2sq*rnhb2sq
      ohm2 = hr2*hs2*rhb2

      tjoho1 = hc2*hs2sq*rhb2
      tjoho2 = hr2*hs2*rhb2
      psip2 = (r2sq*hs2sq+2.d0*(hc2*hs2*s2-b2*tjoho1))*rnhb2sq
      ohmp2 = (r2*hs2+hr2*s2-b2*tjoho2)*rhb2
      psip0 = (r0sq*hs2sq+2.d0*(hc2*hs2*s0-b0*tjoho1))*rnhb2sq
      ohmp0 = (r0*hs2+hr2*s0-b0*tjoho2)*rhb2

      WW2 = -psip2
      XX2 = -psip2-ohmp2-(ohm2+1.d0)*b2*rhb2
      YY2 = psip2+2.d0*ohmp2-(3.d0*psi2-ohm2-2.d0)*b2*rhb2
      ZZ2 = (psi2-1.d0)*b2*rhb2-ohmp2
      T12 = WW2+YY2+3.d0*ZZ2
      T22 = 0.5d0*XX2+2.5d0*YY2+5.5d0*ZZ2+3.d0*WW2
      ZET2 = -(1.d0+T12)*dlog(xsib)*raeta-T12*rxsib3+
     *(T22+2.d0*ZZ2-XX2)*aeta*rxsib3-T22*aeta2*rxsib3+
     *WW2*aeta3*rxsib3
      aex2 = t1ae+T12*aeT1+(T22+2.d0*ZZ2-XX2)*aeTXZ+
     *T22*aeT2+WW2*aeW

      WW0 = -psip0
      XX0 = -psip0-ohmp0-(ohm2+1.d0)*b0*rhb2
      YY0 = psip0+2.d0*ohmp0-(3.d0*psi2-ohm2-2.d0)*b0*rhb2
      ZZ0 = (psi2-1.d0)*b0*rhb2-ohmp0
      T10 = WW0+YY0+3.d0*ZZ0
      T20 = 0.5d0*XX0+2.5d0*YY0+5.5d0*ZZ0+3.d0*WW0
      ZET0 = -(1.d0+T10)*dlog(xsib)*raeta-T10*rxsib3+
     *(T20+2.d0*ZZ0-XX0)*aeta*rxsib3-T20*aeta2*rxsib3+
     *WW0*aeta3*rxsib3
      aex0 = t1ae+T10*aeT1+(T20+2.d0*ZZ0-XX0)*aeTXZ+
     *T20*aeT2+WW0*aeW

      psi1 = hc1*hs1sq*rnhb1sq
      ohm1 = hr1*hs1*rhb1
      psip1 = (r1sq*hs1sq+2.d0*hc1*hs1*s1-
     *2.d0*b1*hc1*hs1sq*rhb1)*rnhb1sq
      ohmp1 = (r1*hs1+hr1*s1-b1*hr1*hs1*rhb1)*rhb1
      WW1 = -psip1
      XX1 = -psip1-ohmp1-(ohm1+1.d0)*b1*rhb1
      YY1 = psip1+2.d0*ohmp1-(3.d0*psi1-ohm1-2.d0)*b1*rhb1
      ZZ1 = (psi1-1.d0)*b1*rhb1-ohmp1
      T11 = WW1+YY1+3.d0*ZZ1
      T21 = 0.5d0*XX1+2.5d0*YY1+5.5d0*ZZ1+3.d0*WW1
      ZET1 = -(1.d0+T11)*dlog(xsib)*raeta-T11*rxsib3+
     *(T21+2.d0*ZZ1-XX1)*aeta*rxsib3-T21*aeta2*rxsib3+
     *WW1*aeta3*rxsib3
      aex1 = t1ae+T11*aeT1+(T21+2.d0*ZZ1-XX1)*aeTXZ+
     *T21*aeT2+WW1*aeW

      phim = bdm/(bdm+bds*q3)
      phis = 1.d0-phim
      Ymb = phim*Yfact+phis*veq
      Z = xp*((Ymb+1.d0)*ZET2-Ymb*ZET1)+xs*ZET0
      write(*,*) 'Z = ',Z
      Pex = (bds+bdpol)*(Z-1.d0)
      write(*,*) 'Pex = ',Pex

      rrVt2 = 1.d0/(bdm+bds*q3)**2
      rNpdYmdNp = (bdm*bds*q3*Yfact-bdm*bds*q3*veq)*rrVt2
      rNpdYmdNs = (-bdpol*bdm*q3*Yfact+bdpol*bdm*q3*veq)*rrVt2

      rNt1 = bdm+bds
      rNt2 = 0.5d0*bdm+bds
      rrNt1 = 1.d0/rNt1
      rrNt2 = 1.d0/rNt2
      dxdNp1 = rnmon*xsm*rrNt1
      dxdNp2 = 0.5d0*rnmon*xsd*rrNt2
      dxdNs1 = -xm*rrNt1
      dxdNs2 = -xd*rrNt2

      dhcdNp1 = dxdNp1*(r1sq-r0sq) 
      dhcdNp2 = dxdNp2*(r2sq-r0sq) 
      dhcdNs1 = dxdNs1*(r1sq-r0sq) 
      dhcdNs2 = dxdNs2*(r2sq-r0sq) 
      dhbdNp1 = dxdNp1*(b1-b0) 
      dhbdNp2 = dxdNp2*(b2-b0) 
      dhbdNs1 = dxdNs1*(b1-b0) 
      dhbdNs2 = dxdNs2*(b2-b0) 
      dhsdNp1 = dxdNp1*(s1-s0) 
      dhsdNp2 = dxdNp2*(s2-s0) 
      dhsdNs1 = dxdNs1*(s1-s0) 
      dhsdNs2 = dxdNs2*(s2-s0) 
      dhrdNp1 = dxdNp1*(r1-r0) 
      dhrdNp2 = dxdNp2*(r2-r0) 
      dhrdNs1 = dxdNs1*(r1-r0) 
      dhrdNs2 = dxdNs2*(r2-r0) 

      tfact1 = dhsdNp2*hs2
      tfact2 = dhcdNp2*hs2+hc2*dhsdNp2
      tfact3 = (dhcdNp2*hs2sq+
     *2.d0*dhsdNp2*hs2*hc2-
     *dhbdNp2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdNp2*rhb2
      dpsipdNp0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdNp1 = 2.d0*((r1sq*dhsdNp1*hs1+
     *s1*(dhcdNp1*hs1+
     *hc1*dhsdNp1)-b1*(dhcdNp1*hs1sq+
     *2.d0*dhsdNp1*hs1*hc1-
     *dhbdNp1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdNp1*psip1*rhb1)
      dpsipdNp2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact1 = dhsdNs2*hs2
      tfact2 = dhcdNs2*hs2+hc2*dhsdNs2
      tfact3 = (dhcdNs2*hs2sq+
     *2.d0*dhsdNs2*hs2*hc2-
     *dhbdNs2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdNs2*rhb2
      dpsipdNs0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdNs1 = 2.d0*((r1sq*dhsdNs1*hs1+
     *s1*(dhcdNs1*hs1+
     *hc1*dhsdNs1)-b1*(dhcdNs1*hs1sq+
     *2.d0*dhsdNs1*hs1*hc1-
     *dhbdNs1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdNs1*psip1*rhb1)
      dpsipdNs2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact = rhb2*(dhrdNp2*hs2+dhsdNp2*hr2-
     *dhbdNp2*hr2*hs2*rhb2)
      dohmpdNp0 = (r0*dhsdNp2+s0*dhrdNp2-
     *b0*tfact)*rhb2-
     *dhbdNp2*ohmp0*rhb2
      dohmpdNp1 = (r1*dhsdNp1+s1*dhrdNp1-
     *b1*rhb1*(dhrdNp1*hs1+dhsdNp1*hr1-
     *dhbdNp1*hr1*hs1*rhb1))*rhb1-
     *dhbdNp1*ohmp1*rhb1
      dohmpdNp2 = (r2*dhsdNp2+s2*dhrdNp2-
     *b2*tfact)*rhb2-
     *dhbdNp2*ohmp2*rhb2

      tfact = rhb2*(dhrdNs2*hs2+dhsdNs2*hr2-
     *dhbdNs2*hr2*hs2*rhb2)
      dohmpdNs0 = (r0*dhsdNs2+s0*dhrdNs2-
     *b0*tfact)*rhb2-
     *dhbdNs2*ohmp0*rhb2
      dohmpdNs1 = (r1*dhsdNs1+s1*dhrdNs1-
     *b1*rhb1*(dhrdNs1*hs1+dhsdNs1*hr1-
     *dhbdNs1*hr1*hs1*rhb1))*rhb1-
     *dhbdNs1*ohmp1*rhb1
      dohmpdNs2 = (r2*dhsdNs2+s2*dhrdNs2-
     *b2*tfact)*rhb2-
     *dhbdNs2*ohmp2*rhb2

      tfact1 = hr1*hs1*rhb1
      tfact2 = hr2*hs2*rhb2
      dohmdNp1 = (dhrdNp1*hs1+dhsdNp1*hr1-
     *tfact1*dhbdNp1)*rhb1
      dohmdNp2 = (dhrdNp2*hs2+dhsdNp2*hr2-
     *tfact2*dhbdNp2)*rhb2
      dohmdNs1 = (dhrdNs1*hs1+dhsdNs1*hr1-
     *tfact1*dhbdNs1)*rhb1
      dohmdNs2 = (dhrdNs2*hs2+dhsdNs2*hr2-
     *tfact2*dhbdNs2)*rhb2

      tfact1 = 2.d0*hc1*hs1sq*rhb1
      tfact2 = 2.d0*hc2*hs2sq*rhb2
      dpsidNp1 = (dhcdNp1*hs1sq+2.d0*dhsdNp1*hc1*hs1-
     *tfact1*dhbdNp1)*rnhb1sq
      dpsidNp2 = (dhcdNp2*hs2sq+2.d0*dhsdNp2*hc2*hs2-
     *tfact2*dhbdNp2)*rnhb2sq
      dpsidNs1 = (dhcdNs1*hs1sq+2.d0*dhsdNs1*hc1*hs1-
     *tfact1*dhbdNs1)*rnhb1sq
      dpsidNs2 = (dhcdNs2*hs2sq+2.d0*dhsdNs2*hc2*hs2-
     *tfact2*dhbdNs2)*rnhb2sq

      dXXdNp0 = -dpsipdNp0-dohmpdNp0-
     *dohmdNp2*b0/hb2+(ohm2+1.d0)*dhbdNp2*b0*rhb2sq
      dYYdNp0 = dpsipdNp0+2.d0*dohmpdNp0-
     *(3.d0*dpsidNp2-dohmdNp2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdNp2*rhb2sq
      dZZdNp0 = dpsidNp2*b0*rhb2-
     *dhbdNp2*(psi2-1.d0)*b0*rhb2sq-dohmpdNp0
      dXXdNp1 = -dpsipdNp1-dohmpdNp1-
     *dohmdNp1*b1*rhb1+(ohm1+1.d0)*dhbdNp1*b1*rhb1sq
      dYYdNp1 = dpsipdNp1+2.d0*dohmpdNp1-
     *(3.d0*dpsidNp1-dohmdNp1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdNp1*rhb1sq
      dZZdNp1 = dpsidNp1*b1*rhb1-
     *dhbdNp1*(psi1-1.d0)*b1*rhb1sq-dohmpdNp1
      dXXdNp2 = -dpsipdNp2-dohmpdNp2-
     *dohmdNp2*b2*rhb2+(ohm2+1.d0)*dhbdNp2*b2*rhb2sq
      dYYdNp2 = dpsipdNp2+2.d0*dohmpdNp2-
     *(3.d0*dpsidNp2-dohmdNp2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdNp2*rhb2sq
      dZZdNp2 = dpsidNp2*b2*rhb2-
     *dhbdNp2*(psi2-1.d0)*b2*rhb2sq-dohmpdNp2
      dT1dNp0 = -dpsipdNp0+dYYdNp0+3.d0*dZZdNp0
      dT2dNp0 = 0.5d0*dXXdNp0+2.5d0*dYYdNp0+5.5d0*dZZdNp0-
     *3.d0*dpsipdNp0
      dT1dNp1 = -dpsipdNp1+dYYdNp1+3.d0*dZZdNp1
      dT2dNp1 = 0.5d0*dXXdNp1+2.5d0*dYYdNp1+5.5d0*dZZdNp1-
     *3.d0*dpsipdNp1
      dT1dNp2 = -dpsipdNp2+dYYdNp2+3.d0*dZZdNp2
      dT2dNp2 = 0.5d0*dXXdNp2+2.5d0*dYYdNp2+5.5d0*dZZdNp2-
     *3.d0*dpsipdNp2

      dXXdNs0 = -dpsipdNs0-dohmpdNs0-
     *dohmdNs2*b0*rhb2+(ohm2+1.d0)*dhbdNs2*b0*rhb2sq
      dYYdNs0 = dpsipdNs0+2.d0*dohmpdNs0-
     *(3.d0*dpsidNs2-dohmdNs2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdNs2*rhb2sq
      dZZdNs0 = dpsidNs2*b0*rhb2-
     *dhbdNs2*(psi2-1.d0)*b0*rhb2sq-dohmpdNs0
      dXXdNs1 = -dpsipdNs1-dohmpdNs1-
     *dohmdNs1*b1*rhb1+(ohm1+1.d0)*dhbdNs1*b1*rhb1sq
      dYYdNs1 = dpsipdNs1+2.d0*dohmpdNs1-
     *(3.d0*dpsidNs1-dohmdNs1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdNs1*rhb1sq
      dZZdNs1 = dpsidNs1*b1*rhb1-
     *dhbdNs1*(psi1-1.d0)*b1*rhb1sq-dohmpdNs1
      dXXdNs2 = -dpsipdNs2-dohmpdNs2-
     *dohmdNs2*b2*rhb2+(ohm2+1.d0)*dhbdNs2*b2*rhb2sq
      dYYdNs2 = dpsipdNs2+2.d0*dohmpdNs2-
     *(3.d0*dpsidNs2-dohmdNs2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdNs2*rhb2sq
      dZZdNs2 = dpsidNs2*b2*rhb2-
     *dhbdNs2*(psi2-1.d0)*b2*rhb2sq-dohmpdNs2
      dT1dNs0 = -dpsipdNs0+dYYdNs0+3.d0*dZZdNs0
      dT2dNs0 = 0.5d0*dXXdNs0+2.5d0*dYYdNs0+5.5d0*dZZdNs0-
     *3.d0*dpsipdNs0
      dT1dNs1 = -dpsipdNs1+dYYdNs1+3.d0*dZZdNs1
      dT2dNs1 = 0.5d0*dXXdNs1+2.5d0*dYYdNs1+5.5d0*dZZdNs1-
     *3.d0*dpsipdNs1
      dT1dNs2 = -dpsipdNs2+dYYdNs2+3.d0*dZZdNs2
      dT2dNs2 = 0.5d0*dXXdNs2+2.5d0*dYYdNs2+5.5d0*dZZdNs2-
     *3.d0*dpsipdNs2

      fdaex2dNp = pis*rnmon*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNp2*aeT1+(dT2dNp2+2.d0*dZZdNp2-dXXdNp2)*aeTXZ+
     *dT2dNp2*aeT2-dpsipdNp2*aeW

      fdaex1dNp = pis*rnmon*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNp1*aeT1+(dT2dNp1+2.d0*dZZdNp1-dXXdNp1)*aeTXZ+
     *dT2dNp1*aeT2-dpsipdNp1*aeW

      fdaex0dNp = pis*rnmon*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNp0*aeT1+(dT2dNp0+2.d0*dZZdNp0-dXXdNp0)*aeTXZ+
     *dT2dNp0*aeT2-dpsipdNp0*aeW

      exchempp = Ymb*(aex2-aex1)+aex2+
     *bdpol*(Ymb*(fdaex2dNp-fdaex1dNp)+fdaex2dNp)+
     *bds*fdaex0dNp+rNpdYmdNp*(aex2-aex1)

      fdaex2dNs = pis*q3*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNs2*aeT1+(dT2dNs2+2.d0*dZZdNs2-dXXdNs2)*aeTXZ+
     *dT2dNs2*aeT2-dpsipdNs2*aeW

      fdaex1dNs = pis*q3*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNs1*aeT1+(dT2dNs1+2.d0*dZZdNs1-dXXdNs1)*aeTXZ+
     *dT2dNs1*aeT2-dpsipdNs1*aeW

      fdaex0dNs = pis*q3*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNs0*aeT1+(dT2dNs0+2.d0*dZZdNs0-dXXdNs0)*aeTXZ+
     *dT2dNs0*aeT2-dpsipdNs0*aeW

      exchemps = aex0+
     *bdpol*(Ymb*(fdaex2dNs-fdaex1dNs)+fdaex2dNs)+
     *bds*fdaex0dNs+rNpdYmdNs*(aex2-aex1)

      csbl = bdpol*(Ymb*(fdaex2dNs-fdaex1dNs)+fdaex2dNs)+
     *bds*fdaex0dNs+rNpdYmdNs*(aex2-aex1)
      write(*,*) 'csbl = ',csbl
      write(*,*) 'exchempp = ',exchempp
      write(*,*) 'exchemps = ',exchemps
      Pbl = bdpol+bds+Pex
      chemppl = dlog(bdpol)+exchempp-(rnmon-1.d0)*dlog(fourpi)
      chempsl = dlog(bds)+exchemps

      scalem = chemppl/(2.d0*rnmon)
      emscale = 2.d0*scalem
      scales = chempsl
      bclamb = rrcmon*Ymb*(aex2-aex1)+
     *rrnmon*(bdpol*(Ymb*(fdaex2dNp-fdaex1dNp)+fdaex2dNp)+
     *bds*fdaex0dNp+rNpdYmdNp*(aex2-aex1))
      belamb = 0.5d0*aex2+
     *rrnmon*(bdpol*(Ymb*(fdaex2dNp-fdaex1dNp)+fdaex2dNp)+
     *bds*fdaex0dNp+rNpdYmdNp*(aex2-aex1))
      bslamb = exchemps
      write(*,*) 'Analytical bulk lambda values (liq. phase)'
      write(*,*) 'bclamb = ',bclamb
      write(*,*) 'belamb = ',belamb
      write(*,*) 'bslamb = ',bslamb
      elhblcmb = dexp(-0.5d0*bclamb+scalem)
      elblemb = dexp(-belamb+emscale)
      elblsmb = dexp(-bslamb+scales)
      bde = 2.d0*bdm*rrnmon
      exFreen = (bdm-bde)*rrcmon*Ymb*(aex2-aex1)+
     *0.5d0*bde*aex2+bds*aex0
      btrams = -(bdm-bde)*bclamb-bde*belamb-bdm*rrnmon+
     *bds*(dlog(bds)-1.d0-chempsl)+exFreen 
      write(*,*) 'exFreen,btrams (l) = ',exFreen,btrams 
      btramsl = btrams

      bdm = bdmg
      bds = bdsg
      bdt = bds*q3+bdm
      bdpol = bdm/rnmon
      cdmbulk = bdm
      cdsbulk = bds*q3

      aeta = pis*bdt
      xsib = 1.d0-aeta
      rxsib = 1.d0/xsib
      rxsib2 = rxsib*rxsib
      rxsib3 = rxsib*rxsib2
      raeta = 1.d0/aeta
      aeta2 = aeta*aeta
      aeta3 = aeta2*aeta
      useful = -dlog(xsib)*raeta-xsib*dlog(xsib)/aeta2-raeta

      xs = bds/(bds+bdpol)
      xp = 1.d0-xs
      xd = 0.5d0*bdm/(0.5d0*bdm+bds)
      xm = bdm/(bdm+bds)
      xsd = 1.d0-xd
      xsm = 1.d0-xm
      hr2 = xd*r2+xsd*r0
      hr1 = xm*r1+xsm*r0
      hs2 = xd*s2+xsd*s0
      hs1 = xm*s1+xsm*s0
      hb2 = xd*b2+xsd*b0
      hb1 = xm*b1+xsm*b0
      hc2 = xd*r2sq+xsd*r0sq
      hc1 = xm*r1sq+xsm*r0sq
      rhb2 = 1.d0/hb2
      rhb2sq = rhb2*rhb2
      rnhb2sq = rhb2sq/9.d0
      hs2sq = hs2*hs2
      hb2sq = hb2*hb2
      rhb1 = 1.d0/hb1
      rhb1sq = rhb1*rhb1
      rnhb1sq = rhb1sq/9.d0      
      hs1sq = hs1*hs1
      hb1sq = hb1*hb1
      write(*,*) 'xd,xm = ',xd,xm

      t1ae = xsib*dlog(xsib)*raeta
      aeT1 = dlog(xsib)*(xsib*raeta+1.d0)-rxsib-0.5d0*rxsib2+2.5d0
      aeTXZ = 0.5d0*rxsib2-0.5d0
      aeT2 = (rxsib-0.5d0*rxsib2-0.5d0)
      aeW = 0.5d0*rxsib2-2.d0*rxsib-dlog(xsib)+1.5d0

      psi2 = hc2*hs2sq*rnhb2sq
      ohm2 = hr2*hs2*rhb2

      tjoho1 = hc2*hs2sq*rhb2
      tjoho2 = hr2*hs2*rhb2
      psip2 = (r2sq*hs2sq+2.d0*(hc2*hs2*s2-b2*tjoho1))*rnhb2sq
      ohmp2 = (r2*hs2+hr2*s2-b2*tjoho2)*rhb2
      psip0 = (r0sq*hs2sq+2.d0*(hc2*hs2*s0-b0*tjoho1))*rnhb2sq
      ohmp0 = (r0*hs2+hr2*s0-b0*tjoho2)*rhb2

      WW2 = -psip2
      XX2 = -psip2-ohmp2-(ohm2+1.d0)*b2*rhb2
      YY2 = psip2+2.d0*ohmp2-(3.d0*psi2-ohm2-2.d0)*b2*rhb2
      ZZ2 = (psi2-1.d0)*b2*rhb2-ohmp2
      T12 = WW2+YY2+3.d0*ZZ2
      T22 = 0.5d0*XX2+2.5d0*YY2+5.5d0*ZZ2+3.d0*WW2
      ZET2 = -(1.d0+T12)*dlog(xsib)*raeta-T12*rxsib3+
     *(T22+2.d0*ZZ2-XX2)*aeta*rxsib3-T22*aeta2*rxsib3+
     *WW2*aeta3*rxsib3
      aex2 = t1ae+T12*aeT1+(T22+2.d0*ZZ2-XX2)*aeTXZ+
     *T22*aeT2+WW2*aeW

      WW0 = -psip0
      XX0 = -psip0-ohmp0-(ohm2+1.d0)*b0*rhb2
      YY0 = psip0+2.d0*ohmp0-(3.d0*psi2-ohm2-2.d0)*b0*rhb2
      ZZ0 = (psi2-1.d0)*b0*rhb2-ohmp0
      T10 = WW0+YY0+3.d0*ZZ0
      T20 = 0.5d0*XX0+2.5d0*YY0+5.5d0*ZZ0+3.d0*WW0
      ZET0 = -(1.d0+T10)*dlog(xsib)*raeta-T10*rxsib3+
     *(T20+2.d0*ZZ0-XX0)*aeta*rxsib3-T20*aeta2*rxsib3+
     *WW0*aeta3*rxsib3
      aex0 = t1ae+T10*aeT1+(T20+2.d0*ZZ0-XX0)*aeTXZ+
     *T20*aeT2+WW0*aeW

      psi1 = hc1*hs1sq*rnhb1sq
      ohm1 = hr1*hs1*rhb1
      psip1 = (r1sq*hs1sq+2.d0*hc1*hs1*s1-
     *2.d0*b1*hc1*hs1sq*rhb1)*rnhb1sq
      ohmp1 = (r1*hs1+hr1*s1-b1*hr1*hs1*rhb1)*rhb1
      WW1 = -psip1
      XX1 = -psip1-ohmp1-(ohm1+1.d0)*b1*rhb1
      YY1 = psip1+2.d0*ohmp1-(3.d0*psi1-ohm1-2.d0)*b1*rhb1
      ZZ1 = (psi1-1.d0)*b1*rhb1-ohmp1
      T11 = WW1+YY1+3.d0*ZZ1
      T21 = 0.5d0*XX1+2.5d0*YY1+5.5d0*ZZ1+3.d0*WW1
      ZET1 = -(1.d0+T11)*dlog(xsib)*raeta-T11*rxsib3+
     *(T21+2.d0*ZZ1-XX1)*aeta*rxsib3-T21*aeta2*rxsib3+
     *WW1*aeta3*rxsib3
      aex1 = t1ae+T11*aeT1+(T21+2.d0*ZZ1-XX1)*aeTXZ+
     *T21*aeT2+WW1*aeW

      phim = bdm/(bdm+bds*q3)
      phis = 1.d0-phim
      Ymb = phim*Yfact+phis*veq
      Z = xp*((Ymb+1.d0)*ZET2-Ymb*ZET1)+xs*ZET0
      write(*,*) 'Z = ',Z
      Pex = (bds+bdpol)*(Z-1.d0)
      write(*,*) 'Pex = ',Pex

      rrVt2 = 1.d0/(bdm+bds*q3)**2
      rNpdYmdNp = (bdm*bds*q3*Yfact-bdm*bds*q3*veq)*rrVt2
      rNpdYmdNs = (-bdpol*bdm*q3*Yfact+bdpol*bdm*q3*veq)*rrVt2

      rNt1 = bdm+bds
      rNt2 = 0.5d0*bdm+bds
      rrNt1 = 1.d0/rNt1
      rrNt2 = 1.d0/rNt2
      dxdNp1 = rnmon*xsm*rrNt1
      dxdNp2 = 0.5d0*rnmon*xsd*rrNt2
      dxdNs1 = -xm*rrNt1
      dxdNs2 = -xd*rrNt2

      dhcdNp1 = dxdNp1*(r1sq-r0sq) 
      dhcdNp2 = dxdNp2*(r2sq-r0sq) 
      dhcdNs1 = dxdNs1*(r1sq-r0sq) 
      dhcdNs2 = dxdNs2*(r2sq-r0sq) 
      dhbdNp1 = dxdNp1*(b1-b0) 
      dhbdNp2 = dxdNp2*(b2-b0) 
      dhbdNs1 = dxdNs1*(b1-b0) 
      dhbdNs2 = dxdNs2*(b2-b0) 
      dhsdNp1 = dxdNp1*(s1-s0) 
      dhsdNp2 = dxdNp2*(s2-s0) 
      dhsdNs1 = dxdNs1*(s1-s0) 
      dhsdNs2 = dxdNs2*(s2-s0) 
      dhrdNp1 = dxdNp1*(r1-r0) 
      dhrdNp2 = dxdNp2*(r2-r0) 
      dhrdNs1 = dxdNs1*(r1-r0) 
      dhrdNs2 = dxdNs2*(r2-r0) 

      tfact1 = dhsdNp2*hs2
      tfact2 = dhcdNp2*hs2+hc2*dhsdNp2
      tfact3 = (dhcdNp2*hs2sq+
     *2.d0*dhsdNp2*hs2*hc2-
     *dhbdNp2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdNp2*rhb2
      dpsipdNp0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdNp1 = 2.d0*((r1sq*dhsdNp1*hs1+
     *s1*(dhcdNp1*hs1+
     *hc1*dhsdNp1)-b1*(dhcdNp1*hs1sq+
     *2.d0*dhsdNp1*hs1*hc1-
     *dhbdNp1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdNp1*psip1*rhb1)
      dpsipdNp2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact1 = dhsdNs2*hs2
      tfact2 = dhcdNs2*hs2+hc2*dhsdNs2
      tfact3 = (dhcdNs2*hs2sq+
     *2.d0*dhsdNs2*hs2*hc2-
     *dhbdNs2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdNs2*rhb2
      dpsipdNs0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdNs1 = 2.d0*((r1sq*dhsdNs1*hs1+
     *s1*(dhcdNs1*hs1+
     *hc1*dhsdNs1)-b1*(dhcdNs1*hs1sq+
     *2.d0*dhsdNs1*hs1*hc1-
     *dhbdNs1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdNs1*psip1*rhb1)
      dpsipdNs2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact = rhb2*(dhrdNp2*hs2+dhsdNp2*hr2-
     *dhbdNp2*hr2*hs2*rhb2)
      dohmpdNp0 = (r0*dhsdNp2+s0*dhrdNp2-
     *b0*tfact)*rhb2-
     *dhbdNp2*ohmp0*rhb2
      dohmpdNp1 = (r1*dhsdNp1+s1*dhrdNp1-
     *b1*rhb1*(dhrdNp1*hs1+dhsdNp1*hr1-
     *dhbdNp1*hr1*hs1*rhb1))*rhb1-
     *dhbdNp1*ohmp1*rhb1
      dohmpdNp2 = (r2*dhsdNp2+s2*dhrdNp2-
     *b2*tfact)*rhb2-
     *dhbdNp2*ohmp2*rhb2

      tfact = rhb2*(dhrdNs2*hs2+dhsdNs2*hr2-
     *dhbdNs2*hr2*hs2*rhb2)
      dohmpdNs0 = (r0*dhsdNs2+s0*dhrdNs2-
     *b0*tfact)*rhb2-
     *dhbdNs2*ohmp0*rhb2
      dohmpdNs1 = (r1*dhsdNs1+s1*dhrdNs1-
     *b1*rhb1*(dhrdNs1*hs1+dhsdNs1*hr1-
     *dhbdNs1*hr1*hs1*rhb1))*rhb1-
     *dhbdNs1*ohmp1*rhb1
      dohmpdNs2 = (r2*dhsdNs2+s2*dhrdNs2-
     *b2*tfact)*rhb2-
     *dhbdNs2*ohmp2*rhb2

      tfact1 = hr1*hs1*rhb1
      tfact2 = hr2*hs2*rhb2
      dohmdNp1 = (dhrdNp1*hs1+dhsdNp1*hr1-
     *tfact1*dhbdNp1)*rhb1
      dohmdNp2 = (dhrdNp2*hs2+dhsdNp2*hr2-
     *tfact2*dhbdNp2)*rhb2
      dohmdNs1 = (dhrdNs1*hs1+dhsdNs1*hr1-
     *tfact1*dhbdNs1)*rhb1
      dohmdNs2 = (dhrdNs2*hs2+dhsdNs2*hr2-
     *tfact2*dhbdNs2)*rhb2

      tfact1 = 2.d0*hc1*hs1sq*rhb1
      tfact2 = 2.d0*hc2*hs2sq*rhb2
      dpsidNp1 = (dhcdNp1*hs1sq+2.d0*dhsdNp1*hc1*hs1-
     *tfact1*dhbdNp1)*rnhb1sq
      dpsidNp2 = (dhcdNp2*hs2sq+2.d0*dhsdNp2*hc2*hs2-
     *tfact2*dhbdNp2)*rnhb2sq
      dpsidNs1 = (dhcdNs1*hs1sq+2.d0*dhsdNs1*hc1*hs1-
     *tfact1*dhbdNs1)*rnhb1sq
      dpsidNs2 = (dhcdNs2*hs2sq+2.d0*dhsdNs2*hc2*hs2-
     *tfact2*dhbdNs2)*rnhb2sq

      dXXdNp0 = -dpsipdNp0-dohmpdNp0-
     *dohmdNp2*b0/hb2+(ohm2+1.d0)*dhbdNp2*b0*rhb2sq
      dYYdNp0 = dpsipdNp0+2.d0*dohmpdNp0-
     *(3.d0*dpsidNp2-dohmdNp2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdNp2*rhb2sq
      dZZdNp0 = dpsidNp2*b0*rhb2-
     *dhbdNp2*(psi2-1.d0)*b0*rhb2sq-dohmpdNp0
      dXXdNp1 = -dpsipdNp1-dohmpdNp1-
     *dohmdNp1*b1*rhb1+(ohm1+1.d0)*dhbdNp1*b1*rhb1sq
      dYYdNp1 = dpsipdNp1+2.d0*dohmpdNp1-
     *(3.d0*dpsidNp1-dohmdNp1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdNp1*rhb1sq
      dZZdNp1 = dpsidNp1*b1*rhb1-
     *dhbdNp1*(psi1-1.d0)*b1*rhb1sq-dohmpdNp1
      dXXdNp2 = -dpsipdNp2-dohmpdNp2-
     *dohmdNp2*b2*rhb2+(ohm2+1.d0)*dhbdNp2*b2*rhb2sq
      dYYdNp2 = dpsipdNp2+2.d0*dohmpdNp2-
     *(3.d0*dpsidNp2-dohmdNp2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdNp2*rhb2sq
      dZZdNp2 = dpsidNp2*b2*rhb2-
     *dhbdNp2*(psi2-1.d0)*b2*rhb2sq-dohmpdNp2
      dT1dNp0 = -dpsipdNp0+dYYdNp0+3.d0*dZZdNp0
      dT2dNp0 = 0.5d0*dXXdNp0+2.5d0*dYYdNp0+5.5d0*dZZdNp0-
     *3.d0*dpsipdNp0
      dT1dNp1 = -dpsipdNp1+dYYdNp1+3.d0*dZZdNp1
      dT2dNp1 = 0.5d0*dXXdNp1+2.5d0*dYYdNp1+5.5d0*dZZdNp1-
     *3.d0*dpsipdNp1
      dT1dNp2 = -dpsipdNp2+dYYdNp2+3.d0*dZZdNp2
      dT2dNp2 = 0.5d0*dXXdNp2+2.5d0*dYYdNp2+5.5d0*dZZdNp2-
     *3.d0*dpsipdNp2

      dXXdNs0 = -dpsipdNs0-dohmpdNs0-
     *dohmdNs2*b0*rhb2+(ohm2+1.d0)*dhbdNs2*b0*rhb2sq
      dYYdNs0 = dpsipdNs0+2.d0*dohmpdNs0-
     *(3.d0*dpsidNs2-dohmdNs2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdNs2*rhb2sq
      dZZdNs0 = dpsidNs2*b0*rhb2-
     *dhbdNs2*(psi2-1.d0)*b0*rhb2sq-dohmpdNs0
      dXXdNs1 = -dpsipdNs1-dohmpdNs1-
     *dohmdNs1*b1*rhb1+(ohm1+1.d0)*dhbdNs1*b1*rhb1sq
      dYYdNs1 = dpsipdNs1+2.d0*dohmpdNs1-
     *(3.d0*dpsidNs1-dohmdNs1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdNs1*rhb1sq
      dZZdNs1 = dpsidNs1*b1*rhb1-
     *dhbdNs1*(psi1-1.d0)*b1*rhb1sq-dohmpdNs1
      dXXdNs2 = -dpsipdNs2-dohmpdNs2-
     *dohmdNs2*b2*rhb2+(ohm2+1.d0)*dhbdNs2*b2*rhb2sq
      dYYdNs2 = dpsipdNs2+2.d0*dohmpdNs2-
     *(3.d0*dpsidNs2-dohmdNs2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdNs2*rhb2sq
      dZZdNs2 = dpsidNs2*b2*rhb2-
     *dhbdNs2*(psi2-1.d0)*b2*rhb2sq-dohmpdNs2
      dT1dNs0 = -dpsipdNs0+dYYdNs0+3.d0*dZZdNs0
      dT2dNs0 = 0.5d0*dXXdNs0+2.5d0*dYYdNs0+5.5d0*dZZdNs0-
     *3.d0*dpsipdNs0
      dT1dNs1 = -dpsipdNs1+dYYdNs1+3.d0*dZZdNs1
      dT2dNs1 = 0.5d0*dXXdNs1+2.5d0*dYYdNs1+5.5d0*dZZdNs1-
     *3.d0*dpsipdNs1
      dT1dNs2 = -dpsipdNs2+dYYdNs2+3.d0*dZZdNs2
      dT2dNs2 = 0.5d0*dXXdNs2+2.5d0*dYYdNs2+5.5d0*dZZdNs2-
     *3.d0*dpsipdNs2

      fdaex2dNp = pis*rnmon*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNp2*aeT1+(dT2dNp2+2.d0*dZZdNp2-dXXdNp2)*aeTXZ+
     *dT2dNp2*aeT2-dpsipdNp2*aeW

      fdaex1dNp = pis*rnmon*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNp1*aeT1+(dT2dNp1+2.d0*dZZdNp1-dXXdNp1)*aeTXZ+
     *dT2dNp1*aeT2-dpsipdNp1*aeW

      fdaex0dNp = pis*rnmon*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNp0*aeT1+(dT2dNp0+2.d0*dZZdNp0-dXXdNp0)*aeTXZ+
     *dT2dNp0*aeT2-dpsipdNp0*aeW

      exchempp = Ymb*(aex2-aex1)+aex2+
     *bdpol*(Ymb*(fdaex2dNp-fdaex1dNp)+fdaex2dNp)+
     *bds*fdaex0dNp+rNpdYmdNp*(aex2-aex1)

      fdaex2dNs = pis*q3*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNs2*aeT1+(dT2dNs2+2.d0*dZZdNs2-dXXdNs2)*aeTXZ+
     *dT2dNs2*aeT2-dpsipdNs2*aeW

      fdaex1dNs = pis*q3*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNs1*aeT1+(dT2dNs1+2.d0*dZZdNs1-dXXdNs1)*aeTXZ+
     *dT2dNs1*aeT2-dpsipdNs1*aeW

      fdaex0dNs = pis*q3*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dNs0*aeT1+(dT2dNs0+2.d0*dZZdNs0-dXXdNs0)*aeTXZ+
     *dT2dNs0*aeT2-dpsipdNs0*aeW

      exchemps = aex0+
     *bdpol*(Ymb*(fdaex2dNs-fdaex1dNs)+fdaex2dNs)+
     *bds*fdaex0dNs+rNpdYmdNs*(aex2-aex1)

      csb = bdpol*(Ymb*(fdaex2dNs-fdaex1dNs)+fdaex2dNs)+
     *bds*fdaex0dNs+rNpdYmdNs*(aex2-aex1)
      write(*,*) 'csb = ',csb
      write(*,*) 'exchempp = ',exchempp
      write(*,*) 'exchemps = ',exchemps
      Pb = bdpol+bds+Pex
      chempp = dlog(bdpol)+exchempp-(rnmon-1.d0)*dlog(fourpi)
      chemps = dlog(bds)+exchemps
c      Fexp = rNp*(Ymb*(aex2-aex1)+aex2)
c      Fexs = rNs*aex0
c      F = rNp*(dlog(bdpol)-1.d0)-rNp*(rnmon-1)*dlog(fourpi)+
c     *Fexp+rNs*(dlog(bds)-1.d0)+Fexs

      scalem = chempp/(2.d0*rnmon)
      emscale = 2.d0*scalem
      scales = chemps
      write(*,*) 'bdm,bdpol =',bdm,bdpol 
      write(*,*) 'monomer density in bulk = ',bdm
      write(*,*) 'solvent density in bulk (in sigmon units)= ',bds
      write(*,*) 'solvent density in bulk (in sigsolv units)= ',bds*q3
      write(*,*) 'no. of monomers/polymer = ',nmon
      write(*,*) 'max no. of iterations = ',ioimaxm
      write(*,*) 'polymer chemical pot. (betamu,g) = '
      write(*,'(1e25.14)') chempp
      write(*,*) 'polymer chemical pot. (betamu,l) = '
      write(*,'(1e25.14)') chemppl
      write(*,*) 'solvent chemical pot. (betamu,g) = '
      write(*,'(1e25.14)') chemps
      write(*,*) 'solvent chemical pot. (betamu,l) = '
      write(*,'(1e25.14)') chempsl
      write(*,*) 'total bulk pressure (g) = '
      write(*,'(1e25.14)') Pb
      write(*,*) 'total bulk pressure (l) = '
      write(*,'(1e25.14)') Pbl
      write(*,*) 'temperature = ',T
      write(*,*) 'bdist,dmax = ',bdist,dmax
      write(*,*) 'dz = ',dz
      write(*,*) 'closew = ',closew
      write(*,*) 'istp1 = ',istp1
      write(*,*) 'imitt = ',imitt
      write(*,*) 'ism,isms = ',ism,isms
      write(*,*) 'dmm,dms (density mixing param. mon.,solv.) = ',dmm,dms
      write(*,*) 'hwdist = ',hwdist

      bclamb = rrcmon*Ymb*(aex2-aex1)+
     *rrnmon*(bdpol*(Ymb*(fdaex2dNp-fdaex1dNp)+fdaex2dNp)+
     *bds*fdaex0dNp+rNpdYmdNp*(aex2-aex1))
      belamb = 0.5d0*aex2+
     *rrnmon*(bdpol*(Ymb*(fdaex2dNp-fdaex1dNp)+fdaex2dNp)+
     *bds*fdaex0dNp+rNpdYmdNp*(aex2-aex1))
      bslamb = exchemps
      bbclamb = bclamb
      bbelamb = belamb
      baex0 = aex0
      baex1 = aex1
      baex2 = aex2
      write(*,*) 'Analytical bulk lambda values'
      write(*,*) 'bclamb = ',bclamb
      write(*,*) 'belamb = ',belamb
      write(*,*) 'bslamb = ',bslamb
      ehblcmb = dexp(-0.5d0*bclamb+scalem)
      eblemb = dexp(-belamb+emscale)
      eblsmb = dexp(-bslamb+scales)
      write(*,*) 'elhblcmb,ehblcmb = ',elhblcmb,ehblcmb 
      write(*,*) 'elblemb,eblemb = ',elblemb,eblemb 
      write(*,*) 'elblsmb,eblsmb = ',elblsmb,eblsmb 
      bde = 2.d0*bdm*rrnmon
      exFreen = (bdm-bde)*rrcmon*Ymb*(aex2-aex1)+
     *0.5d0*bde*aex2+bds*aex0
      btrams = -(bdm-bde)*bclamb-bde*belamb-bdm*rrnmon+
     *bds*(dlog(bds)-1.d0-chemps)+exFreen 
      write(*,*) 'exFreen,btrams (g) = ',exFreen,btrams 
      btramsg = btrams


      bdld = bdml*real(idg-istp1-ism)
      bdgd = bdmg*real(imitt-ism-idg)
      bdsum = bdld+bdgd+0.5d0*(bdml+bdmg)
      sbdld = bdsl*real(idg-istp1-ism)
      sbdgd = bdsg*real(imitt-ism-idg)
      sbdsum = sbdld+sbdgd+0.5d0*(bdsl+bdsg)
      write(*,*) 'idg-istp1-ism,imitt-ism-idg = ',idg-istp1-ism,
     *imitt-ism-idg

cc     unshifted version
c      closew = closew-0.5d0*dz
      ifin = imitt+ism+1
      ist = imitt-ism
      kfin = istp1+ism
      kst = istp1-ism-1
      rrjkdiff = 1.d0/real(2*ism+1)



      iefin = imitt+ism+1
      iest = imitt-3*ism
      kefin = istp1+3*ism
      kest = istp1-ism-1
      rrejkdiff = 1.d0/real(4*ism+1)



      dsum = 0.d0
      sdsum = 0.d0
c      do 747 iz = minel,0
c      fdsol(iz) = bdsl
c 747  fdmon(iz) = bdml
      do 12 iz = 0,istp1-ism
      fdsol(iz) = bdsl
      fdmon(iz) = bdml
 12   fem(iz) = 2.d0*bdml*rrnmon
      do 24 iz = imitt+ism+1,maxel
      fdsol(iz) = bds
      fdmon(iz) = bdm
 24   fem(iz) = 2.d0*bdm*rrnmon
      if (kread.eq.0) then

c      z = distp1+1.d0-dz
      z = distp1+1.d0-0.5d0*dz

      do 412 iz = istp1+ism,imitt-ism
      z = z+dz
      if (abs((z-dgibbs)*alpha).lt.30.d0) then
      fem(iz) = bfemg+(bfeml-bfemg)/(1.d0+dexp(alpha*(z-dgibbs)))
      fdmon(iz) = bdmg+(bdml-bdmg)/(1.d0+dexp(alpha*(z-dgibbs)))
      fdsol(iz) = bdsg+(bdsl-bdsg)/(1.d0+dexp(alpha*(z-dgibbs)))
      else
      if (iz.lt.idg) then
      fem(iz) = bfeml
      fdmon(iz) = bdml
      fdsol(iz) = bdsl
      else
      fem(iz) = bfemg
      fdmon(iz) = bdmg
      fdsol(iz) = bdsg
      endif
      endif
      dsum = dsum+fdmon(iz)
 412  sdsum = sdsum+fdsol(iz)
      eplag = bdsum/dsum
      seplag = sbdsum/sdsum
      else
      read(ieplag,'(1e25.14)') eplag
      read(ieplag,'(1e25.14)') seplag
      do 13 iz = istp1-ism,imitt+ism
      read(ifc,'(1f11.4,3e25.14)') trams,fdmon(iz),fdsol(iz),fem(iz)
 13   fdsol(iz) = fdsol(iz)*rq3
      endif
      webelam = fem(istp1+ism)
      weblam = fdsol(istp1+ism)
      wehbclam = fdmon(istp1+ism)      
      do 776 iz = istp1-ism,istp1+ism-1
      rkk = real(iz-kst)
      rjj = real(kfin-iz)
      fem(iz) = (rkk*webelam+rjj*2.d0*bdml*rrnmon)*rrjkdiff
      fdsol(iz) = (rkk*weblam+rjj*bdsl)*rrjkdiff
 776  fdmon(iz) = (rkk*wehbclam+rjj*bdml)*rrjkdiff
      webelam = fem(imitt-ism)
      weblam = fdsol(imitt-ism)
      wehbclam = fdmon(imitt-ism)
      do 676 iz = imitt-ism+1,imitt+ism
      rkk = real(iz-ist)
      rjj = real(ifin-iz)
      fem(iz) = (rjj*webelam+rkk*2.d0*bdpol)*rrjkdiff
      fdsol(iz) = (rjj*weblam+rkk*bds)*rrjkdiff
 676  fdmon(iz) = (rjj*wehbclam+rkk*bdm)*rrjkdiff
      do 502 iz = 0,istp1-ism
      eblam(iz) = elblsmb
      ebelam(iz) = elblemb
 502  ehbclam(iz) = elhblcmb 
      do 504 iz = imitt+ism,maxel
      eblam(iz) = eblsmb
      ebelam(iz) = eblemb
 504  ehbclam(iz) = ehblcmb 
      do 518 iz = 0,istp1-ism
      cdmonm(iz) = bdml
      cdmons(iz) = bdml
      cdsolm(iz) = bdsl*q3
 518  cdsols(iz) = bdsl*q3
      do 519 iz = imitt+ism,maxel
      cdmonm(iz) = bdmg
      cdmons(iz) = bdmg
      cdsolm(iz) = bdsg*q3
 519  cdsols(iz) = bdsg*q3
      eplag = 1.d0
      seplag = 1.d0
      write(*,*) 'bdsum,dsum = ',bdsum,dsum
      write(*,*) 'eplag = ',eplag
      write(*,*) 'sbdsum,sdsum = ',sbdsum,sdsum
      write(*,*) 'seplag = ',seplag
c      call flush
      ddmax = 10000.
      niter = 0
 100  continue
      niter = niter+1    
      CALL CDMM
      CALL CDMS
      CALL CDSM
      CALL CDSS
c     SOME USEFUL VECTORS ARE OBTAINED IN AVEC
      CALL AVECa12
      CALL AVECa0
c     EVALUATION OF exp(-0.5*beta*lambda(iz)/2)=ehbclam(iz)
      CALL EBLMNEW
      CALL EBLSNEW

      if (ddmax.lt.ddtol) goto 200
      if (mod(niter,100).eq.0) then
      write(*,*) 'ddmax,niter = ',ddmax,niter
c      call flush
      endif
      if (niter.gt.ioimaxm) then
      write(*,*) 'NITER.GT.IOIMAXM !',niter
      goto 200
      endif
      ddmax = 0.
      nAB = 1
      do 445 iz = ism,imitt+(nmon-3)*ism
      sume = 0.5d0*ebelam(iz-ism)
      do 545 jz = iz-ism+1,iz+ism-1
 545  sume = ebelam(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*ebelam(iz+ism)+sume)*
     *twopidz
      chA(iz) = tuu*ehbclam(iz)
 445  c(nmon-1,iz) = tuu

      k = nmon-1
      do 745 mmm = 2,nmon-2 
      k = k-1
      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do 1045 iz = mmm*ism,imitt+(k-2)*ism
      sume = 0.5d0*chA(iz-ism)
      do 1145 jz = iz-ism+1,iz+ism-1
 1145 sume = chA(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
 1045 c(k,iz) = tuu
      else
      do 3045 iz = mmm*ism,imitt+(k-2)*ism
      sume = 0.5d0*chB(iz-ism)
      do 3145 jz = iz-ism+1,iz+ism-1
 3145 sume = chB(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 3045 c(k,iz) = tuu
      endif
 745  continue

      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do 1445 iz = istp1+ism,imitt-ism
      sume = 0.5d0*chA(iz-ism)
      do 1545 jz = iz-ism+1,iz+ism-1
 1545 sume = chA(jz)+sume
 1445 c(1,iz) = ebelam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      else
      do 3445 iz = istp1+ism,imitt-ism
      sume = 0.5d0*chB(iz-ism)
      do 3545 jz = iz-ism+1,iz+ism-1
 3545 sume = chB(jz)+sume
 3445 c(1,iz) = ebelam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      endif

      dsum = 0.d0
      sdsum = 0.d0
      do 9 i = istp1+ism,imitt-ism
      dumsum = 0.d0 
      do 10 k = 2,nmon-1
 10   dumsum = c(k,i)*c(nmon+1-k,i)+dumsum
      tfem = 2.d0*c(1,i)
      tfdm = dumsum+tfem
      ddiff = abs(tfdm-fdmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fem(i) = fem(i)*dmm+tdmm*tfem
      fdmon(i) = fdmon(i)*dmm+tdmm*tfdm
      tfds = eblam(i)
      ddiff = abs(tfds-fdsol(i))/tfds
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdsol(i) = fdsol(i)*dms+tdms*tfds
      dsum = fdmon(i)+dsum
 9    sdsum = fdsol(i)+sdsum
c      eplag = bdsum/dsum
c      seplag = sbdsum/sdsum
      webelam = fem(istp1+ism)
      weblam = fdsol(istp1+ism)
      wehbclam = fdmon(istp1+ism)      
      do 786 iz = istp1-ism,istp1+ism-1
      rkk = real(iz-kst)
      rjj = real(kfin-iz)
      fem(iz) = (rkk*webelam+rjj*2.d0*bdml*rrnmon)*rrjkdiff
      fdsol(iz) = (rkk*weblam+rjj*bdsl)*rrjkdiff
 786  fdmon(iz) = (rkk*wehbclam+rjj*bdml)*rrjkdiff
      webelam = fem(imitt-ism)
      weblam = fdsol(imitt-ism)
      wehbclam = fdmon(imitt-ism)
      do 76 iz = imitt-ism+1,imitt+ism
      rkk = real(iz-ist)
      rjj = real(ifin-iz)
      fem(iz) = (rjj*webelam+rkk*2.d0*bdpol)*rrjkdiff
c      fdsol(iz) = (rkk*weblam+rjj*bds)*rrjkdiff
      fdsol(iz) = (rjj*weblam+rkk*bds)*rrjkdiff
 76   fdmon(iz) = (rjj*wehbclam+rkk*bdm)*rrjkdiff

      goto 100
 200  continue
      rewind 25
      sumexFreen = 0.d0
      niz = 0
      sumfdm = 0.d0
      sumfds = 0.d0
      Freen = 0.d0
      plag = 0.d0
      splag = 0.d0
      iz = istp1+5*ism+1
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      fds = fdsol(iz)
      fdm = fdmon(iz)      
      fde = fem(iz)
      fdc = fdm-fde
      exFreen = (fdm-fem(iz))*rrcmon*Ymix(iz)*(ae2(iz)-ae1(iz))+
     *0.5d0*fde*ae2(iz)+fds*ae0(iz)
      sumexFreen = exFreen+sumexFreen
      weblam = fdc*bclamb+fde*belamb+fdm*(plag-rrnmon)+fds*(dlog(fds)-
     *1.d0-chemps+splag)+exFreen+Pb

      kefin = istp1+5*ism+1
      kest = istp1-ism
      rrejkdiff = 1.d0/real(6*ism+1)
c      btrams = -Pb
      btrams = 0.d0
      write(*,*) btrams,weblam
      do iz = istp1-ism+1,istp1+5*ism
      rkk = real(iz-kest)
      rjj = real(kefin-iz)
      fds = fdsol(iz)
      fdm = fdmon(iz)      
      sumfdm = fdm+sumfdm
      sumfds = fds+sumfds
      Fwc = (rkk*weblam+rjj*btrams)*rrejkdiff
      write(25,'(2e25.14)') real(iz)*dz,Fwc
      Freen = Fwc+Freen
      enddo

      do 500 iz = istp1+5*ism+1,imitt-5*ism-1
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      niz = niz+1
      fds = fdsol(iz)
      fdm = fdmon(iz)      
      fde = fem(iz)
      fdc = fdm-fde
      sumfdm = fdm+sumfdm
      sumfds = fds+sumfds
      exFreen = (fdm-fem(iz))*rrcmon*Ymix(iz)*(ae2(iz)-ae1(iz))+
     *0.5d0*fde*ae2(iz)+fds*ae0(iz)
      sumexFreen = exFreen+sumexFreen
      Fwc = fdc*bclamb+fde*belamb+fdm*(plag-rrnmon)+fds*(dlog(fds)-
     *1.d0-chemps+splag)+exFreen+Pb
      write(25,'(2e25.14)') real(iz)*dz,Fwc
 500  Freen = Fwc+Freen

      kefin = imitt+ism
      kest = imitt-5*ism-1
      weblam = Fwc
      do iz = imitt-5*ism,imitt+ism-1
      rkk = real(iz-kest)
      rjj = real(kefin-iz)
      fds = fdsol(iz)
      fdm = fdmon(iz)      
      sumfdm = fdm+sumfdm
      sumfds = fds+sumfds
      Fwc = (rjj*weblam+rkk*btrams)*rrejkdiff
      write(25,'(2e25.14)') real(iz)*dz,Fwc
      Freen = Fwc+Freen
      enddo

      write(*,*) 'exFreen = ',sumexFreen*dz
      write(*,*) 'sumfdm*plag = ',sumfds*plag*dz
      write(*,*) 'sumfds*splag = ',sumfds*splag*dz
      write(*,*) 'eplag = ',eplag
      write(*,*) 'seplag = ',seplag
      Freen = Freen*dz
c      write(*,*) '/beta*(grand pot.) = ',Freen
c      write(*,*) '(grand pot.)/epsmm = ',Freen*T      
      rewind ifc
      do 61 iz = istp1-ism,imitt+ism
      z = real(iz)*dz
 61   write(ifc,'(1f11.5,3e25.14)') z,fdmon(iz),fdsol(iz)*q3,fem(iz)
      rewind ieplag
      write(ieplag,'(1e25.14)') eplag
      write(ieplag,'(1e25.14)') seplag
      dch = real(imitt)*dz-real(istp1)*dz+dz
c      Fsture = Freen+Pb*dch
c      Fsture = Freen+Pb*(dch+2.d0)
      write(*,*)
      write(*,*) '/beta*gs:'
      write(*,'(1e25.14)') Freen
 9999 continue
      STOP
      END

      subroutine CDMM
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
c      z = distp1+1.d0-dz
      z = distp1+1.d0-0.5d0*dz
      do 36 iz = istp1+ism,imitt-ism
      z = z+dz
      zs = z-1.d0
      z1 = zs
      z2 = zs+dz
      f1 = fdmon(iz-ism)
      f2 = fdmon(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
 36   cdmonm(iz) = 0.75d0*sancint
      weblam = cdmonm(istp1+ism)
      do 676 iz = istp1-ism,istp1+ism-1
      rkk = real(iz-kst)
      rjj = real(kfin-iz)
 676  cdmonm(iz) = (rkk*weblam+rjj*cdmlbulk)*rrjkdiff
      weblam = cdmonm(imitt-ism)
      do 876 iz = imitt-ism+1,imitt+ism
      rkk = real(iz-ist)
      rjj = real(ifin-iz)
 876  cdmonm(iz) = (rjj*weblam+rkk*cdmbulk)*rrjkdiff
      return
      end



      subroutine CDMS
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
      z = distp1+1.d0-0.5d0*dz
      do 36 iz = istp1+ism,imitt-ism
      z = z+dz
      zs = z-p1
      z1 = zs
      z2 = zs+dz
      f1 = fdmon(iz-ismms)
      f2 = fdmon(iz-ismms+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-ismms+1,iz+ismms-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdmon(iz+ismms)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*p2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
 36   cdmons(iz) = 0.75d0*sancint*rp3

      weblam = cdmons(istp1+ism)
      do 676 iz = istp1-ism,istp1+ism-1
      rkk = real(iz-kst)
      rjj = real(kfin-iz)
 676  cdmons(iz) = (rkk*weblam+rjj*cdmlbulk)*rrjkdiff
      weblam = cdmons(imitt-ism)
      do 876 iz = imitt-ism+1,imitt+ism
      rkk = real(iz-ist)
      rjj = real(ifin-iz)
 876  cdmons(iz) = (rjj*weblam+rkk*cdmbulk)*rrjkdiff
      return
      end

      subroutine CDSS
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
c      z = distp1+1.d0-dz
      z = distp1+1.d0-0.5d0*dz
      do 36 iz = istp1+ism,imitt-ism
      z = z+dz
      zs = z-q1
      z1 = zs
      z2 = zs+dz
      f1 = fdsol(iz-isms)
      f2 = fdsol(iz-isms+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-isms+1,iz+isms-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdsol(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdsol(iz+isms)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*q2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
 36   cdsols(iz) = 0.75d0*sancint
      weblam = cdsols(istp1+ism)
      do 676 iz = istp1-ism,istp1+ism-1
      rkk = real(iz-kst)
      rjj = real(kfin-iz)
 676  cdsols(iz) = (rkk*weblam+rjj*cdslbulk)*rrjkdiff
      weblam = cdsols(imitt-ism)
      do 876 iz = imitt-ism+1,imitt+ism
      rkk = real(iz-ist)
      rjj = real(ifin-iz)
 876  cdsols(iz) = (rjj*weblam+rkk*cdsbulk)*rrjkdiff
      return
      end

      subroutine CDSM
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
c      z = distp1+1.d0-dz
      z = distp1+1.d0-0.5d0*dz
      do 36 iz = istp1+ism,imitt-ism
      z = z+dz
      zs = z-p1
      z1 = zs
      z2 = zs+dz
      f1 = fdsol(iz-ismms)
      f2 = fdsol(iz-ismms+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-ismms+1,iz+ismms-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdsol(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = fdsol(iz+ismms)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*p2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
 36   cdsolm(iz) = 0.75d0*sancint*rp3*q3
      weblam = cdsolm(istp1+ism)
      do 676 iz = istp1-ism,istp1+ism-1
      rkk = real(iz-kst)
      rjj = real(kfin-iz)
 676  cdsolm(iz) = (rkk*weblam+rjj*cdslbulk)*rrjkdiff
      weblam = cdsolm(imitt-ism)
      do 876 iz = imitt-ism+1,imitt+ism
      rkk = real(iz-ist)
      rjj = real(ifin-iz)
 876  cdsolm(iz) = (rjj*weblam+rkk*cdsbulk)*rrjkdiff
      return
      end


      subroutine AVECa12
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
      do 72 iz = istp1,imitt
      bdt = cdmonm(iz)+cdsolm(iz)
      bds = cdsolm(iz)*rq3
      bdm = cdmonm(iz)

      aeta = pis*bdt
      xsib = 1.d0-aeta
      rxsib = 1.d0/xsib
      rxsib2 = rxsib*rxsib
      rxsib3 = rxsib*rxsib2
      raeta = 1.d0/aeta
      aeta2 = aeta*aeta
      aeta3 = aeta2*aeta
      useful = -dlog(xsib)*raeta-xsib*dlog(xsib)/aeta2-raeta

      xd = 0.5d0*bdm/(0.5d0*bdm+bds)
      xm = bdm/(bdm+bds)
      xsd = 1.d0-xd
      xsm = 1.d0-xm
      hr2 = xd*r2+xsd*r0
      hr1 = xm*r1+xsm*r0
      hs2 = xd*s2+xsd*s0
      hs1 = xm*s1+xsm*s0
      hb2 = xd*b2+xsd*b0
      hb1 = xm*b1+xsm*b0
      hc2 = xd*r2sq+xsd*r0sq
      hc1 = xm*r1sq+xsm*r0sq
      rhb2 = 1.d0/hb2
      rhb2sq = rhb2*rhb2
      rnhb2sq = rhb2sq/9.d0
      hs2sq = hs2*hs2
      hb2sq = hb2*hb2
      rhb1 = 1.d0/hb1
      rhb1sq = rhb1*rhb1
      rnhb1sq = rhb1sq/9.d0      
      hs1sq = hs1*hs1
      hb1sq = hb1*hb1

      t1ae = xsib*dlog(xsib)*raeta
      aeT1 = dlog(xsib)*(xsib*raeta+1.d0)-rxsib-0.5d0*rxsib2+2.5d0
      aeTXZ = 0.5d0*rxsib2-0.5d0
      aeT2 = (rxsib-0.5d0*rxsib2-0.5d0)
      aeW = 0.5d0*rxsib2-2.d0*rxsib-dlog(xsib)+1.5d0

      psi2 = hc2*hs2sq*rnhb2sq
      ohm2 = hr2*hs2*rhb2

      tjoho1 = hc2*hs2sq*rhb2
      tjoho2 = hr2*hs2*rhb2
      psip2 = (r2sq*hs2sq+2.d0*(hc2*hs2*s2-b2*tjoho1))*rnhb2sq
      ohmp2 = (r2*hs2+hr2*s2-b2*tjoho2)*rhb2
      psip0 = (r0sq*hs2sq+2.d0*(hc2*hs2*s0-b0*tjoho1))*rnhb2sq
      ohmp0 = (r0*hs2+hr2*s0-b0*tjoho2)*rhb2
      WW2 = -psip2
      XX2 = -psip2-ohmp2-(ohm2+1.d0)*b2*rhb2
      YY2 = psip2+2.d0*ohmp2-(3.d0*psi2-ohm2-2.d0)*b2*rhb2
      ZZ2 = (psi2-1.d0)*b2*rhb2-ohmp2
      T12 = WW2+YY2+3.d0*ZZ2
      T22 = 0.5d0*XX2+2.5d0*YY2+5.5d0*ZZ2+3.d0*WW2
      aex2 = t1ae+T12*aeT1+(T22+2.d0*ZZ2-XX2)*aeTXZ+
     *T22*aeT2+WW2*aeW

      WW0 = -psip0
      XX0 = -psip0-ohmp0-(ohm2+1.d0)*b0*rhb2
      YY0 = psip0+2.d0*ohmp0-(3.d0*psi2-ohm2-2.d0)*b0*rhb2
      ZZ0 = (psi2-1.d0)*b0*rhb2-ohmp0
      T10 = WW0+YY0+3.d0*ZZ0
      T20 = 0.5d0*XX0+2.5d0*YY0+5.5d0*ZZ0+3.d0*WW0
      aex0 = t1ae+T10*aeT1+(T20+2.d0*ZZ0-XX0)*aeTXZ+
     *T20*aeT2+WW0*aeW

      psi1 = hc1*hs1sq*rnhb1sq
      ohm1 = hr1*hs1*rhb1
      psip1 = (r1sq*hs1sq+2.d0*hc1*hs1*s1-
     *2.d0*b1*hc1*hs1sq*rhb1)*rnhb1sq
      ohmp1 = (r1*hs1+hr1*s1-b1*hr1*hs1*rhb1)*rhb1
      WW1 = -psip1
      XX1 = -psip1-ohmp1-(ohm1+1.d0)*b1*rhb1
      YY1 = psip1+2.d0*ohmp1-(3.d0*psi1-ohm1-2.d0)*b1*rhb1
      ZZ1 = (psi1-1.d0)*b1*rhb1-ohmp1
      T11 = WW1+YY1+3.d0*ZZ1
      T21 = 0.5d0*XX1+2.5d0*YY1+5.5d0*ZZ1+3.d0*WW1
      aex1 = t1ae+T11*aeT1+(T21+2.d0*ZZ1-XX1)*aeTXZ+
     *T21*aeT2+WW1*aeW

      phim = bdm/(bdm+bds*q3)
      phis = 1.d0-phim
      Ymb = phim*Yfact+phis*veq
      rrVt2 = 1.d0/(bdm+bds*q3)**2
      dYmdbdm = (bds*q3*Yfact-bds*q3*veq)*rrVt2
      dYmdbds = (-bdm*q3*Yfact+bdm*q3*veq)*rrVt2

      rNt1 = bdm+bds
      rNt2 = 0.5d0*bdm+bds
      rrNt1 = 1.d0/rNt1
      rrNt2 = 1.d0/rNt2
      dxdbdm1 = xsm*rrNt1
      dxdbdm2 = 0.5d0*xsd*rrNt2
      dxdbds1 = -xm*rrNt1
      dxdbds2 = -xd*rrNt2

      dhcdbdm1 = dxdbdm1*(r1sq-r0sq) 
      dhcdbdm2 = dxdbdm2*(r2sq-r0sq) 
      dhcdbds1 = dxdbds1*(r1sq-r0sq) 
      dhcdbds2 = dxdbds2*(r2sq-r0sq) 
      dhbdbdm1 = dxdbdm1*(b1-b0) 
      dhbdbdm2 = dxdbdm2*(b2-b0) 
      dhbdbds1 = dxdbds1*(b1-b0) 
      dhbdbds2 = dxdbds2*(b2-b0) 
      dhsdbdm1 = dxdbdm1*(s1-s0) 
      dhsdbdm2 = dxdbdm2*(s2-s0) 
      dhsdbds1 = dxdbds1*(s1-s0) 
      dhsdbds2 = dxdbds2*(s2-s0) 
      dhrdbdm1 = dxdbdm1*(r1-r0) 
      dhrdbdm2 = dxdbdm2*(r2-r0) 
      dhrdbds1 = dxdbds1*(r1-r0) 
      dhrdbds2 = dxdbds2*(r2-r0) 

      tfact1 = dhsdbdm2*hs2
      tfact2 = dhcdbdm2*hs2+hc2*dhsdbdm2
      tfact3 = (dhcdbdm2*hs2sq+
     *2.d0*dhsdbdm2*hs2*hc2-
     *dhbdbdm2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdbdm2*rhb2
      dpsipdbdm0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdbdm1 = 2.d0*((r1sq*dhsdbdm1*hs1+
     *s1*(dhcdbdm1*hs1+
     *hc1*dhsdbdm1)-b1*(dhcdbdm1*hs1sq+
     *2.d0*dhsdbdm1*hs1*hc1-
     *dhbdbdm1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdbdm1*psip1*rhb1)
      dpsipdbdm2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact1 = dhsdbds2*hs2
      tfact2 = dhcdbds2*hs2+hc2*dhsdbds2
      tfact3 = (dhcdbds2*hs2sq+
     *2.d0*dhsdbds2*hs2*hc2-
     *dhbdbds2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdbds2*rhb2
      dpsipdbds0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdbds1 = 2.d0*((r1sq*dhsdbds1*hs1+
     *s1*(dhcdbds1*hs1+
     *hc1*dhsdbds1)-b1*(dhcdbds1*hs1sq+
     *2.d0*dhsdbds1*hs1*hc1-
     *dhbdbds1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdbds1*psip1*rhb1)
      dpsipdbds2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact = rhb2*(dhrdbdm2*hs2+dhsdbdm2*hr2-
     *dhbdbdm2*hr2*hs2*rhb2)
      dohmpdbdm0 = (r0*dhsdbdm2+s0*dhrdbdm2-
     *b0*tfact)*rhb2-
     *dhbdbdm2*ohmp0*rhb2
      dohmpdbdm1 = (r1*dhsdbdm1+s1*dhrdbdm1-
     *b1*rhb1*(dhrdbdm1*hs1+dhsdbdm1*hr1-
     *dhbdbdm1*hr1*hs1*rhb1))*rhb1-
     *dhbdbdm1*ohmp1*rhb1
      dohmpdbdm2 = (r2*dhsdbdm2+s2*dhrdbdm2-
     *b2*tfact)*rhb2-
     *dhbdbdm2*ohmp2*rhb2

      tfact = rhb2*(dhrdbds2*hs2+dhsdbds2*hr2-
     *dhbdbds2*hr2*hs2*rhb2)
      dohmpdbds0 = (r0*dhsdbds2+s0*dhrdbds2-
     *b0*tfact)*rhb2-
     *dhbdbds2*ohmp0*rhb2
      dohmpdbds1 = (r1*dhsdbds1+s1*dhrdbds1-
     *b1*rhb1*(dhrdbds1*hs1+dhsdbds1*hr1-
     *dhbdbds1*hr1*hs1*rhb1))*rhb1-
     *dhbdbds1*ohmp1*rhb1
      dohmpdbds2 = (r2*dhsdbds2+s2*dhrdbds2-
     *b2*tfact)*rhb2-
     *dhbdbds2*ohmp2*rhb2

      tfact1 = hr1*hs1*rhb1
      tfact2 = hr2*hs2*rhb2
      dohmdbdm1 = (dhrdbdm1*hs1+dhsdbdm1*hr1-
     *tfact1*dhbdbdm1)*rhb1
      dohmdbdm2 = (dhrdbdm2*hs2+dhsdbdm2*hr2-
     *tfact2*dhbdbdm2)*rhb2
      dohmdbds1 = (dhrdbds1*hs1+dhsdbds1*hr1-
     *tfact1*dhbdbds1)*rhb1
      dohmdbds2 = (dhrdbds2*hs2+dhsdbds2*hr2-
     *tfact2*dhbdbds2)*rhb2

      tfact1 = 2.d0*hc1*hs1sq*rhb1
      tfact2 = 2.d0*hc2*hs2sq*rhb2
      dpsidbdm1 = (dhcdbdm1*hs1sq+2.d0*dhsdbdm1*hc1*hs1-
     *tfact1*dhbdbdm1)*rnhb1sq
      dpsidbdm2 = (dhcdbdm2*hs2sq+2.d0*dhsdbdm2*hc2*hs2-
     *tfact2*dhbdbdm2)*rnhb2sq
      dpsidbds1 = (dhcdbds1*hs1sq+2.d0*dhsdbds1*hc1*hs1-
     *tfact1*dhbdbds1)*rnhb1sq
      dpsidbds2 = (dhcdbds2*hs2sq+2.d0*dhsdbds2*hc2*hs2-
     *tfact2*dhbdbds2)*rnhb2sq

      dXXdbdm0 = -dpsipdbdm0-dohmpdbdm0-
     *dohmdbdm2*b0/hb2+(ohm2+1.d0)*dhbdbdm2*b0*rhb2sq
      dYYdbdm0 = dpsipdbdm0+2.d0*dohmpdbdm0-
     *(3.d0*dpsidbdm2-dohmdbdm2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdbdm2*rhb2sq
      dZZdbdm0 = dpsidbdm2*b0*rhb2-
     *dhbdbdm2*(psi2-1.d0)*b0*rhb2sq-dohmpdbdm0
      dXXdbdm1 = -dpsipdbdm1-dohmpdbdm1-
     *dohmdbdm1*b1*rhb1+(ohm1+1.d0)*dhbdbdm1*b1*rhb1sq
      dYYdbdm1 = dpsipdbdm1+2.d0*dohmpdbdm1-
     *(3.d0*dpsidbdm1-dohmdbdm1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdbdm1*rhb1sq
      dZZdbdm1 = dpsidbdm1*b1*rhb1-
     *dhbdbdm1*(psi1-1.d0)*b1*rhb1sq-dohmpdbdm1
      dXXdbdm2 = -dpsipdbdm2-dohmpdbdm2-
     *dohmdbdm2*b2*rhb2+(ohm2+1.d0)*dhbdbdm2*b2*rhb2sq
      dYYdbdm2 = dpsipdbdm2+2.d0*dohmpdbdm2-
     *(3.d0*dpsidbdm2-dohmdbdm2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdbdm2*rhb2sq
      dZZdbdm2 = dpsidbdm2*b2*rhb2-
     *dhbdbdm2*(psi2-1.d0)*b2*rhb2sq-dohmpdbdm2
      dT1dbdm0 = -dpsipdbdm0+dYYdbdm0+3.d0*dZZdbdm0
      dT2dbdm0 = 0.5d0*dXXdbdm0+2.5d0*dYYdbdm0+5.5d0*dZZdbdm0-
     *3.d0*dpsipdbdm0
      dT1dbdm1 = -dpsipdbdm1+dYYdbdm1+3.d0*dZZdbdm1
      dT2dbdm1 = 0.5d0*dXXdbdm1+2.5d0*dYYdbdm1+5.5d0*dZZdbdm1-
     *3.d0*dpsipdbdm1
      dT1dbdm2 = -dpsipdbdm2+dYYdbdm2+3.d0*dZZdbdm2
      dT2dbdm2 = 0.5d0*dXXdbdm2+2.5d0*dYYdbdm2+5.5d0*dZZdbdm2-
     *3.d0*dpsipdbdm2

      dXXdbds0 = -dpsipdbds0-dohmpdbds0-
     *dohmdbds2*b0*rhb2+(ohm2+1.d0)*dhbdbds2*b0*rhb2sq
      dYYdbds0 = dpsipdbds0+2.d0*dohmpdbds0-
     *(3.d0*dpsidbds2-dohmdbds2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdbds2*rhb2sq
      dZZdbds0 = dpsidbds2*b0*rhb2-
     *dhbdbds2*(psi2-1.d0)*b0*rhb2sq-dohmpdbds0
      dXXdbds1 = -dpsipdbds1-dohmpdbds1-
     *dohmdbds1*b1*rhb1+(ohm1+1.d0)*dhbdbds1*b1*rhb1sq
      dYYdbds1 = dpsipdbds1+2.d0*dohmpdbds1-
     *(3.d0*dpsidbds1-dohmdbds1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdbds1*rhb1sq
      dZZdbds1 = dpsidbds1*b1*rhb1-
     *dhbdbds1*(psi1-1.d0)*b1*rhb1sq-dohmpdbds1
      dXXdbds2 = -dpsipdbds2-dohmpdbds2-
     *dohmdbds2*b2*rhb2+(ohm2+1.d0)*dhbdbds2*b2*rhb2sq
      dYYdbds2 = dpsipdbds2+2.d0*dohmpdbds2-
     *(3.d0*dpsidbds2-dohmdbds2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdbds2*rhb2sq
      dZZdbds2 = dpsidbds2*b2*rhb2-
     *dhbdbds2*(psi2-1.d0)*b2*rhb2sq-dohmpdbds2
      dT1dbds0 = -dpsipdbds0+dYYdbds0+3.d0*dZZdbds0
      dT2dbds0 = 0.5d0*dXXdbds0+2.5d0*dYYdbds0+5.5d0*dZZdbds0-
     *3.d0*dpsipdbds0
      dT1dbds1 = -dpsipdbds1+dYYdbds1+3.d0*dZZdbds1
      dT2dbds1 = 0.5d0*dXXdbds1+2.5d0*dYYdbds1+5.5d0*dZZdbds1-
     *3.d0*dpsipdbds1
      dT1dbds2 = -dpsipdbds2+dYYdbds2+3.d0*dZZdbds2
      dT2dbds2 = 0.5d0*dXXdbds2+2.5d0*dYYdbds2+5.5d0*dZZdbds2-
     *3.d0*dpsipdbds2

      fdaex2dbdm = pis*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbdm2*aeT1+(dT2dbdm2+2.d0*dZZdbdm2-dXXdbdm2)*aeTXZ+
     *dT2dbdm2*aeT2-dpsipdbdm2*aeW

      fdaex1dbdm = pis*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbdm1*aeT1+(dT2dbdm1+2.d0*dZZdbdm1-dXXdbdm1)*aeTXZ+
     *dT2dbdm1*aeT2-dpsipdbdm1*aeW

      fdaex0dbdm = pis*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbdm0*aeT1+(dT2dbdm0+2.d0*dZZdbdm0-dXXdbdm0)*aeTXZ+
     *dT2dbdm0*aeT2-dpsipdbdm0*aeW

      convpm(iz) = (fdmon(iz)-fem(iz))*rrcmon*(Ymb*(fdaex2dbdm-
     *fdaex1dbdm)+dYmdbdm*(aex2-aex1))+0.5d0*fem(iz)*fdaex2dbdm

      fdaex2dbds = pis*q3*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbds2*aeT1+(dT2dbds2+2.d0*dZZdbds2-dXXdbds2)*aeTXZ+
     *dT2dbds2*aeT2-dpsipdbds2*aeW

      fdaex1dbds = pis*q3*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbds1*aeT1+(dT2dbds1+2.d0*dZZdbds1-dXXdbds1)*aeTXZ+
     *dT2dbds1*aeT2-dpsipdbds1*aeW

      fdaex0dbds = pis*q3*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbds0*aeT1+(dT2dbds0+2.d0*dZZdbds0-dXXdbds0)*aeTXZ+
     *dT2dbds0*aeT2-dpsipdbds0*aeW

      convsm(iz) = (fdmon(iz)-fem(iz))*rrcmon*(Ymb*(fdaex2dbds-
     *fdaex1dbds)+dYmdbds*(aex2-aex1))+0.5d0*fem(iz)*fdaex2dbds

      ae1(iz) = aex1
      ae2(iz) = aex2
 72   Ymix(iz) = Ymb
      return
      end


      subroutine AVECa0
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
      do 72 iz = istp1,imitt
      bdt = cdmons(iz)+cdsols(iz)
      bds = cdsols(iz)*rq3
      bdm = cdmons(iz)

      aeta = pis*bdt
      xsib = 1.d0-aeta
      rxsib = 1.d0/xsib
      rxsib2 = rxsib*rxsib
      rxsib3 = rxsib*rxsib2
      raeta = 1.d0/aeta
      aeta2 = aeta*aeta
      aeta3 = aeta2*aeta
      useful = -dlog(xsib)*raeta-xsib*dlog(xsib)/aeta2-raeta

      xd = 0.5d0*bdm/(0.5d0*bdm+bds)
      xm = bdm/(bdm+bds)
      xsd = 1.d0-xd
      xsm = 1.d0-xm
      hr2 = xd*r2+xsd*r0
      hr1 = xm*r1+xsm*r0
      hs2 = xd*s2+xsd*s0
      hs1 = xm*s1+xsm*s0
      hb2 = xd*b2+xsd*b0
      hb1 = xm*b1+xsm*b0
      hc2 = xd*r2sq+xsd*r0sq
      hc1 = xm*r1sq+xsm*r0sq
      rhb2 = 1.d0/hb2
      rhb2sq = rhb2*rhb2
      rnhb2sq = rhb2sq/9.d0
      hs2sq = hs2*hs2
      hb2sq = hb2*hb2
      rhb1 = 1.d0/hb1
      rhb1sq = rhb1*rhb1
      rnhb1sq = rhb1sq/9.d0      
      hs1sq = hs1*hs1
      hb1sq = hb1*hb1

      t1ae = xsib*dlog(xsib)*raeta
      aeT1 = dlog(xsib)*(xsib*raeta+1.d0)-rxsib-0.5d0*rxsib2+2.5d0
      aeTXZ = 0.5d0*rxsib2-0.5d0
      aeT2 = (rxsib-0.5d0*rxsib2-0.5d0)
      aeW = 0.5d0*rxsib2-2.d0*rxsib-dlog(xsib)+1.5d0

      psi2 = hc2*hs2sq*rnhb2sq
      ohm2 = hr2*hs2*rhb2

      tjoho1 = hc2*hs2sq*rhb2
      tjoho2 = hr2*hs2*rhb2
      psip2 = (r2sq*hs2sq+2.d0*(hc2*hs2*s2-b2*tjoho1))*rnhb2sq
      ohmp2 = (r2*hs2+hr2*s2-b2*tjoho2)*rhb2
      psip0 = (r0sq*hs2sq+2.d0*(hc2*hs2*s0-b0*tjoho1))*rnhb2sq
      ohmp0 = (r0*hs2+hr2*s0-b0*tjoho2)*rhb2
      WW2 = -psip2
      XX2 = -psip2-ohmp2-(ohm2+1.d0)*b2*rhb2
      YY2 = psip2+2.d0*ohmp2-(3.d0*psi2-ohm2-2.d0)*b2*rhb2
      ZZ2 = (psi2-1.d0)*b2*rhb2-ohmp2
      T12 = WW2+YY2+3.d0*ZZ2
      T22 = 0.5d0*XX2+2.5d0*YY2+5.5d0*ZZ2+3.d0*WW2
      aex2 = t1ae+T12*aeT1+(T22+2.d0*ZZ2-XX2)*aeTXZ+
     *T22*aeT2+WW2*aeW

      WW0 = -psip0
      XX0 = -psip0-ohmp0-(ohm2+1.d0)*b0*rhb2
      YY0 = psip0+2.d0*ohmp0-(3.d0*psi2-ohm2-2.d0)*b0*rhb2
      ZZ0 = (psi2-1.d0)*b0*rhb2-ohmp0
      T10 = WW0+YY0+3.d0*ZZ0
      T20 = 0.5d0*XX0+2.5d0*YY0+5.5d0*ZZ0+3.d0*WW0
      aex0 = t1ae+T10*aeT1+(T20+2.d0*ZZ0-XX0)*aeTXZ+
     *T20*aeT2+WW0*aeW

      psi1 = hc1*hs1sq*rnhb1sq
      ohm1 = hr1*hs1*rhb1
      psip1 = (r1sq*hs1sq+2.d0*hc1*hs1*s1-
     *2.d0*b1*hc1*hs1sq*rhb1)*rnhb1sq
      ohmp1 = (r1*hs1+hr1*s1-b1*hr1*hs1*rhb1)*rhb1
      WW1 = -psip1
      XX1 = -psip1-ohmp1-(ohm1+1.d0)*b1*rhb1
      YY1 = psip1+2.d0*ohmp1-(3.d0*psi1-ohm1-2.d0)*b1*rhb1
      ZZ1 = (psi1-1.d0)*b1*rhb1-ohmp1
      T11 = WW1+YY1+3.d0*ZZ1
      T21 = 0.5d0*XX1+2.5d0*YY1+5.5d0*ZZ1+3.d0*WW1
      aex1 = t1ae+T11*aeT1+(T21+2.d0*ZZ1-XX1)*aeTXZ+
     *T21*aeT2+WW1*aeW

      phim = bdm/(bdm+bds*q3)
      phis = 1.d0-phim
      Ymb = phim*Yfact+phis*veq
      rrVt2 = 1.d0/(bdm+bds*q3)**2
      dYmdbdm = (bds*q3*Yfact-bds*q3*veq)*rrVt2
      dYmdbds = (-bdm*q3*Yfact+bdm*q3*veq)*rrVt2

      rNt1 = bdm+bds
      rNt2 = 0.5d0*bdm+bds
      rrNt1 = 1.d0/rNt1
      rrNt2 = 1.d0/rNt2
      dxdbdm1 = xsm*rrNt1
      dxdbdm2 = 0.5d0*xsd*rrNt2
      dxdbds1 = -xm*rrNt1
      dxdbds2 = -xd*rrNt2

      dhcdbdm1 = dxdbdm1*(r1sq-r0sq) 
      dhcdbdm2 = dxdbdm2*(r2sq-r0sq) 
      dhcdbds1 = dxdbds1*(r1sq-r0sq) 
      dhcdbds2 = dxdbds2*(r2sq-r0sq) 
      dhbdbdm1 = dxdbdm1*(b1-b0) 
      dhbdbdm2 = dxdbdm2*(b2-b0) 
      dhbdbds1 = dxdbds1*(b1-b0) 
      dhbdbds2 = dxdbds2*(b2-b0) 
      dhsdbdm1 = dxdbdm1*(s1-s0) 
      dhsdbdm2 = dxdbdm2*(s2-s0) 
      dhsdbds1 = dxdbds1*(s1-s0) 
      dhsdbds2 = dxdbds2*(s2-s0) 
      dhrdbdm1 = dxdbdm1*(r1-r0) 
      dhrdbdm2 = dxdbdm2*(r2-r0) 
      dhrdbds1 = dxdbds1*(r1-r0) 
      dhrdbds2 = dxdbds2*(r2-r0) 

      tfact1 = dhsdbdm2*hs2
      tfact2 = dhcdbdm2*hs2+hc2*dhsdbdm2
      tfact3 = (dhcdbdm2*hs2sq+
     *2.d0*dhsdbdm2*hs2*hc2-
     *dhbdbdm2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdbdm2*rhb2
      dpsipdbdm0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdbdm1 = 2.d0*((r1sq*dhsdbdm1*hs1+
     *s1*(dhcdbdm1*hs1+
     *hc1*dhsdbdm1)-b1*(dhcdbdm1*hs1sq+
     *2.d0*dhsdbdm1*hs1*hc1-
     *dhbdbdm1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdbdm1*psip1*rhb1)
      dpsipdbdm2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact1 = dhsdbds2*hs2
      tfact2 = dhcdbds2*hs2+hc2*dhsdbds2
      tfact3 = (dhcdbds2*hs2sq+
     *2.d0*dhsdbds2*hs2*hc2-
     *dhbdbds2*hc2*hs2sq*rhb2)*rhb2
      tfact4 = dhbdbds2*rhb2
      dpsipdbds0 = 2.d0*((r0sq*tfact1+s0*tfact2-b0*tfact3)*rnhb2sq-
     *psip0*tfact4)
      dpsipdbds1 = 2.d0*((r1sq*dhsdbds1*hs1+
     *s1*(dhcdbds1*hs1+
     *hc1*dhsdbds1)-b1*(dhcdbds1*hs1sq+
     *2.d0*dhsdbds1*hs1*hc1-
     *dhbdbds1*hc1*hs1sq*rhb1)*rhb1)*rnhb1sq-
     *dhbdbds1*psip1*rhb1)
      dpsipdbds2 = 2.d0*((r2sq*tfact1+s2*tfact2-b2*tfact3)*rnhb2sq-
     *psip2*tfact4)

      tfact = rhb2*(dhrdbdm2*hs2+dhsdbdm2*hr2-
     *dhbdbdm2*hr2*hs2*rhb2)
      dohmpdbdm0 = (r0*dhsdbdm2+s0*dhrdbdm2-
     *b0*tfact)*rhb2-
     *dhbdbdm2*ohmp0*rhb2
      dohmpdbdm1 = (r1*dhsdbdm1+s1*dhrdbdm1-
     *b1*rhb1*(dhrdbdm1*hs1+dhsdbdm1*hr1-
     *dhbdbdm1*hr1*hs1*rhb1))*rhb1-
     *dhbdbdm1*ohmp1*rhb1
      dohmpdbdm2 = (r2*dhsdbdm2+s2*dhrdbdm2-
     *b2*tfact)*rhb2-
     *dhbdbdm2*ohmp2*rhb2

      tfact = rhb2*(dhrdbds2*hs2+dhsdbds2*hr2-
     *dhbdbds2*hr2*hs2*rhb2)
      dohmpdbds0 = (r0*dhsdbds2+s0*dhrdbds2-
     *b0*tfact)*rhb2-
     *dhbdbds2*ohmp0*rhb2
      dohmpdbds1 = (r1*dhsdbds1+s1*dhrdbds1-
     *b1*rhb1*(dhrdbds1*hs1+dhsdbds1*hr1-
     *dhbdbds1*hr1*hs1*rhb1))*rhb1-
     *dhbdbds1*ohmp1*rhb1
      dohmpdbds2 = (r2*dhsdbds2+s2*dhrdbds2-
     *b2*tfact)*rhb2-
     *dhbdbds2*ohmp2*rhb2

      tfact1 = hr1*hs1*rhb1
      tfact2 = hr2*hs2*rhb2
      dohmdbdm1 = (dhrdbdm1*hs1+dhsdbdm1*hr1-
     *tfact1*dhbdbdm1)*rhb1
      dohmdbdm2 = (dhrdbdm2*hs2+dhsdbdm2*hr2-
     *tfact2*dhbdbdm2)*rhb2
      dohmdbds1 = (dhrdbds1*hs1+dhsdbds1*hr1-
     *tfact1*dhbdbds1)*rhb1
      dohmdbds2 = (dhrdbds2*hs2+dhsdbds2*hr2-
     *tfact2*dhbdbds2)*rhb2

      tfact1 = 2.d0*hc1*hs1sq*rhb1
      tfact2 = 2.d0*hc2*hs2sq*rhb2
      dpsidbdm1 = (dhcdbdm1*hs1sq+2.d0*dhsdbdm1*hc1*hs1-
     *tfact1*dhbdbdm1)*rnhb1sq
      dpsidbdm2 = (dhcdbdm2*hs2sq+2.d0*dhsdbdm2*hc2*hs2-
     *tfact2*dhbdbdm2)*rnhb2sq
      dpsidbds1 = (dhcdbds1*hs1sq+2.d0*dhsdbds1*hc1*hs1-
     *tfact1*dhbdbds1)*rnhb1sq
      dpsidbds2 = (dhcdbds2*hs2sq+2.d0*dhsdbds2*hc2*hs2-
     *tfact2*dhbdbds2)*rnhb2sq

      dXXdbdm0 = -dpsipdbdm0-dohmpdbdm0-
     *dohmdbdm2*b0/hb2+(ohm2+1.d0)*dhbdbdm2*b0*rhb2sq
      dYYdbdm0 = dpsipdbdm0+2.d0*dohmpdbdm0-
     *(3.d0*dpsidbdm2-dohmdbdm2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdbdm2*rhb2sq
      dZZdbdm0 = dpsidbdm2*b0*rhb2-
     *dhbdbdm2*(psi2-1.d0)*b0*rhb2sq-dohmpdbdm0
      dXXdbdm1 = -dpsipdbdm1-dohmpdbdm1-
     *dohmdbdm1*b1*rhb1+(ohm1+1.d0)*dhbdbdm1*b1*rhb1sq
      dYYdbdm1 = dpsipdbdm1+2.d0*dohmpdbdm1-
     *(3.d0*dpsidbdm1-dohmdbdm1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdbdm1*rhb1sq
      dZZdbdm1 = dpsidbdm1*b1*rhb1-
     *dhbdbdm1*(psi1-1.d0)*b1*rhb1sq-dohmpdbdm1
      dXXdbdm2 = -dpsipdbdm2-dohmpdbdm2-
     *dohmdbdm2*b2*rhb2+(ohm2+1.d0)*dhbdbdm2*b2*rhb2sq
      dYYdbdm2 = dpsipdbdm2+2.d0*dohmpdbdm2-
     *(3.d0*dpsidbdm2-dohmdbdm2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdbdm2*rhb2sq
      dZZdbdm2 = dpsidbdm2*b2*rhb2-
     *dhbdbdm2*(psi2-1.d0)*b2*rhb2sq-dohmpdbdm2
      dT1dbdm0 = -dpsipdbdm0+dYYdbdm0+3.d0*dZZdbdm0
      dT2dbdm0 = 0.5d0*dXXdbdm0+2.5d0*dYYdbdm0+5.5d0*dZZdbdm0-
     *3.d0*dpsipdbdm0
      dT1dbdm1 = -dpsipdbdm1+dYYdbdm1+3.d0*dZZdbdm1
      dT2dbdm1 = 0.5d0*dXXdbdm1+2.5d0*dYYdbdm1+5.5d0*dZZdbdm1-
     *3.d0*dpsipdbdm1
      dT1dbdm2 = -dpsipdbdm2+dYYdbdm2+3.d0*dZZdbdm2
      dT2dbdm2 = 0.5d0*dXXdbdm2+2.5d0*dYYdbdm2+5.5d0*dZZdbdm2-
     *3.d0*dpsipdbdm2

      dXXdbds0 = -dpsipdbds0-dohmpdbds0-
     *dohmdbds2*b0*rhb2+(ohm2+1.d0)*dhbdbds2*b0*rhb2sq
      dYYdbds0 = dpsipdbds0+2.d0*dohmpdbds0-
     *(3.d0*dpsidbds2-dohmdbds2)*b0*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b0*dhbdbds2*rhb2sq
      dZZdbds0 = dpsidbds2*b0*rhb2-
     *dhbdbds2*(psi2-1.d0)*b0*rhb2sq-dohmpdbds0
      dXXdbds1 = -dpsipdbds1-dohmpdbds1-
     *dohmdbds1*b1*rhb1+(ohm1+1.d0)*dhbdbds1*b1*rhb1sq
      dYYdbds1 = dpsipdbds1+2.d0*dohmpdbds1-
     *(3.d0*dpsidbds1-dohmdbds1)*b1*rhb1+
     *(3.d0*psi1-ohm1-2.d0)*b1*dhbdbds1*rhb1sq
      dZZdbds1 = dpsidbds1*b1*rhb1-
     *dhbdbds1*(psi1-1.d0)*b1*rhb1sq-dohmpdbds1
      dXXdbds2 = -dpsipdbds2-dohmpdbds2-
     *dohmdbds2*b2*rhb2+(ohm2+1.d0)*dhbdbds2*b2*rhb2sq
      dYYdbds2 = dpsipdbds2+2.d0*dohmpdbds2-
     *(3.d0*dpsidbds2-dohmdbds2)*b2*rhb2+
     *(3.d0*psi2-ohm2-2.d0)*b2*dhbdbds2*rhb2sq
      dZZdbds2 = dpsidbds2*b2*rhb2-
     *dhbdbds2*(psi2-1.d0)*b2*rhb2sq-dohmpdbds2
      dT1dbds0 = -dpsipdbds0+dYYdbds0+3.d0*dZZdbds0
      dT2dbds0 = 0.5d0*dXXdbds0+2.5d0*dYYdbds0+5.5d0*dZZdbds0-
     *3.d0*dpsipdbds0
      dT1dbds1 = -dpsipdbds1+dYYdbds1+3.d0*dZZdbds1
      dT2dbds1 = 0.5d0*dXXdbds1+2.5d0*dYYdbds1+5.5d0*dZZdbds1-
     *3.d0*dpsipdbds1
      dT1dbds2 = -dpsipdbds2+dYYdbds2+3.d0*dZZdbds2
      dT2dbds2 = 0.5d0*dXXdbds2+2.5d0*dYYdbds2+5.5d0*dZZdbds2-
     *3.d0*dpsipdbds2

      fdaex2dbdm = pis*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbdm2*aeT1+(dT2dbdm2+2.d0*dZZdbdm2-dXXdbdm2)*aeTXZ+
     *dT2dbdm2*aeT2-dpsipdbdm2*aeW

      fdaex1dbdm = pis*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbdm1*aeT1+(dT2dbdm1+2.d0*dZZdbdm1-dXXdbdm1)*aeTXZ+
     *dT2dbdm1*aeT2-dpsipdbdm1*aeW

      fdaex0dbdm = pis*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbdm0*aeT1+(dT2dbdm0+2.d0*dZZdbdm0-dXXdbdm0)*aeTXZ+
     *dT2dbdm0*aeT2-dpsipdbdm0*aeW

      convps(iz) =  fdsol(iz)*fdaex0dbdm

      fdaex2dbds = pis*q3*(useful+
     *T12*(useful-rxsib-rxsib2-rxsib3)+
     *(T22+2.d0*ZZ2-XX2)*rxsib3+
     *T22*(rxsib2-rxsib3)+
     *WW2*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbds2*aeT1+(dT2dbds2+2.d0*dZZdbds2-dXXdbds2)*aeTXZ+
     *dT2dbds2*aeT2-dpsipdbds2*aeW

      fdaex1dbds = pis*q3*(useful+
     *T11*(useful-rxsib-rxsib2-rxsib3)+
     *(T21+2.d0*ZZ1-XX1)*rxsib3+
     *T21*(rxsib2-rxsib3)+
     *WW1*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbds1*aeT1+(dT2dbds1+2.d0*dZZdbds1-dXXdbds1)*aeTXZ+
     *dT2dbds1*aeT2-dpsipdbds1*aeW

      fdaex0dbds = pis*q3*(useful+
     *T10*(useful-rxsib-rxsib2-rxsib3)+
     *(T20+2.d0*ZZ0-XX0)*rxsib3+
     *T20*(rxsib2-rxsib3)+
     *WW0*(rxsib3-2.d0*rxsib2+rxsib))+
     *dT1dbds0*aeT1+(dT2dbds0+2.d0*dZZdbds0-dXXdbds0)*aeTXZ+
     *dT2dbds0*aeT2-dpsipdbds0*aeW

      convss(iz) = fdsol(iz)*fdaex0dbds

 72   ae0(iz) = aex0
      return
      end



      subroutine EBLMNEW
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
c      z = distp1+1.d0-dz
      z = distp1+1.d0-0.5d0*dz
      do 36 iz = istp1+ism,imitt-ism
      z = z+dz
      zs = z-1.d0
      z1 = zs
      z2 = zs+dz
      f1 = convpm(iz-ism)
      f2 = convpm(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convpm(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convpm(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
c 36   trams(iz) = sancint*0.75d0
      trams = sancint*0.75d0

c      do 96 iz = istp1s+ismms,imitt-ism
c      z = z+dz
      zs = z-p1
      z1 = zs
      z2 = zs+dz
      f1 = convps(iz-ismms)
      f2 = convps(iz-ismms+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 97 jz = iz-ismms+1,iz+ismms-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convps(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 97   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convps(iz+ismms)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*p2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
c      ttrams = trams(iz)+sancint*0.75d0*rp3
      ttrams = trams+sancint*0.75d0*rp3
      emtrams = ttrams+0.5d0*ae2(iz)
      cmtrams = ttrams+rrcmon*Ymix(iz)*(ae2(iz)-ae1(iz))
      ebelam(iz) = dexp(-emtrams+emscale)
 36   ehbclam(iz) = dexp(-0.5d0*cmtrams+scalem)

      webelam = ebelam(istp1+3*ism)
      wehbclam = ehbclam(istp1+3*ism)
      do 676 iz = istp1-ism,istp1+3*ism-1
      rkk = real(iz-kest)
      rjj = real(kefin-iz)
      ebelam(iz) = (rkk*webelam+rjj*elblemb)*rrejkdiff
 676  ehbclam(iz) = (rkk*wehbclam+rjj*elhblcmb)*rrejkdiff

      webelam = ebelam(imitt-3*ism)
      wehbclam = ehbclam(imitt-3*ism)
      do 76 iz = imitt-3*ism+1,imitt+ism
      rkk = real(iz-iest)
      rjj = real(iefin-iz)
      ebelam(iz) = (rjj*webelam+rkk*eblemb)*rrejkdiff
 76   ehbclam(iz) = (rjj*wehbclam+rkk*ehblcmb)*rrejkdiff

      return 
      end

      subroutine EBLSNEW
      implicit double precision (a-h,o-z)
      include 'dftpol.rdshall.inc'
c      z = distp1+1.d0-dz
      z = distp1+1.d0-0.5d0*dz
      do 36 iz = istp1+ism,imitt-ism
c      do 46 iz = istp1s+ismms,imitt-ism
      z = z+dz
      zs = z-p1
      z1 = zs
      z2 = zs+dz
      f1 = convsm(iz-ismms)
      f2 = convsm(iz-ismms+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 47 jz = iz-ismms+1,iz+ismms-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convsm(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 47   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convsm(iz+ismms)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*p2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
c 46   trams(iz) = sancint*0.75d0*rp3
      trams = sancint*0.75d0*rp3

c      do 36 iz = istp1s+isms,imitt-ism
c      z = z+dz
      zs = z-q1
      z1 = zs
      z2 = zs+dz
      f1 = convss(iz-isms)
      f2 = convss(iz-isms+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do 37 jz = iz-isms+1,iz+isms-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convss(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
 37   sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convss(iz+isms)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*q2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
c      ttrams = trams(iz)+sancint*0.75d0*rq3
      ttrams = trams+sancint*0.75d0*rq3
      strams = ttrams+ae0(iz)
 36   eblam(iz) = dexp(-strams+scales)

      weblam = eblam(istp1+3*ism)
      do 676 iz = istp1-ism,istp1+3*ism-1
      rkk = real(iz-kest)
      rjj = real(kefin-iz)
 676  eblam(iz) = (rkk*weblam+rjj*elblsmb)*rrejkdiff
      weblam = eblam(imitt-3*ism)
      do 76 iz = imitt-3*ism+1,imitt+ism
      rkk = real(iz-iest)
      rjj = real(iefin-iz)
 76   eblam(iz) = (rjj*weblam+rkk*eblsmb)*rrejkdiff
      return 
      end
      
