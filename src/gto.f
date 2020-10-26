c ==============================================================================
c =                                                                            =
c = Author: Sylwester Wysocki <sw143@wp.pl>                                    =
c = Created on: 2008                                                           =
c =                                                                            =
c = This is free and unencumbered software released into the public domain.    =
c =                                                                            =
c = Anyone is free to copy, modify, publish, use, compile, sell, or            =
c = distribute this software, either in source code form or as a compiled      =
c = binary, for any purpose, commercial or non-commercial, and by any          =
c = means.                                                                     =
c =                                                                            =
c = In jurisdictions that recognize copyright laws, the author or authors      =
c = of this software dedicate any and all copyright interest in the            =
c = software to the public domain. We make this dedication for the benefit     =
c = of the public at large and to the detriment of our heirs and               =
c = successors. We intend this dedication to be an overt act of                =
c = relinquishment in perpetuity of all present and future rights to this      =
c = software under copyright law.                                              =
c =                                                                            =
c = THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,            =
c = EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF         =
c = MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.     =
c = IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR          =
c = OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,      =
c = ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR      =
c = OTHER DEALINGS IN THE SOFTWARE.                                            =
c =                                                                            =
c = For more information, please refer to <https://unlicense.org>              =
c =                                                                            =
c ==============================================================================

c123456
c=============================================
c== Oblicza macierz calek nakrywania <p|q>  ==
c==                                         ==
c==     S = macierz calek nakrywania (WY)   ==
c=============================================

      subroutine ComputeOverlap(S)
        IMPLICIT NONE
        real*8 Base(1000,2)                           ! {wyk,C}
        integer BaseIndex(100)                        ! index w bazie dla danej
                                                      ! liczby atomowej

        integer Norb                                  ! ile orbitali AO
        real*8 AO(1000,5)                             ! lista orbitali AO
                                                      ! (index,m,x,y,z)
        COMMON /BASE/ Base, BaseIndex
        COMMON /AO/ AO,Norb

        integer  p,q
        integer  i,j

        real*8  suma
        real*8  S(Norb,Norb)                         ! macierz cal. nakryw.

        real*8  wp,wq                                ! wykladniki p i q
        real*8  cp,cq                                ! wsp. kontrakcji
        real*8  X                                    ! zm. pomocnicza
        real*8  Rpq,Rpq2                             ! odl. miedzy centrami
        real*8  dx,dy,dz                             ! roznice wsp. cetrow

        integer  lp,lq
        integer  np,nq
        integer  Ip,Iq
        integer  Ilep, Ileq


        do p=1,Norb
          do q=p,Norb
            Ip = AO(p,1)                             ! znajdz indeksy
            Iq = AO(q,1)                             ! orbitali w bazie

            lp = Base(Ip,1)                          ! pobierz typ
            lq = Base(Iq,1)                          ! orbitalu (l)

            Ilep = Base(ip,2)                        ! ile prymitywnych
            Ileq = Base(iq,2)                        ! gaussow?

            ip = ip + 1                              ! przeskocz do
            iq = iq + 1                              ! danych prymitywow

c           ===========================
c           = odleglosc miedzy centr. =
c           ===========================

            dx = AO(p,3)-AO(q,3)
            dy = AO(p,4)-AO(q,4)
            dz = AO(p,5)-AO(q,5)
            Rpq2 = (dx*dx + dy*dy + dz*dz)

c           =========================
c           = licz calke nakrywania =
c           =========================

            suma = 0.0
            do i=ip, ip+Ilep -1
              wp = Base(i,1)
              cp = Base(i,2)

              do j=iq, iq+Ileq -1
                wq = Base(j,1)
                cq = Base(j,2)

                X = (2.*wp*wq)/(wp+wq)

                suma = suma + cp*cq*(wp*wq)**(-0.75)
     &                      * X**(1.5)
     &                      * exp(-0.5*X*Rpq2)
              enddo
            enddo

            S(p,q) = suma
            S(q,p) = suma
          enddo
        enddo
      end subroutine

c=============================================
c== Oblicza macierz calek jednoelekt.<p|h|q>==
c==                                         ==
c==     h = macierz calek nakrywania (WY)   ==
c=============================================

      subroutine ComputeOneElectron (h)
        IMPLICIT NONE
        integer Nat                                  ! ile atomow
        real*8 GEOM(1000,4)                          ! geometria (3 wsp. + Z)
        real*8 Base(1000,2)                          ! {wyk,C}
        integer BaseIndex(100)                       ! index w bazie dla danej
                                                     ! liczby atomowej

        integer Norb                                 ! ile orbitali AO
        real*8 AO(1000,5)                            ! lista orbitali AO
                                                     ! (index,m,x,y,z)

        COMMON /BASE/ Base, BaseIndex
        COMMON /GEOM/ GEOM,Nat
        COMMON /AO/ AO,Norb

        real*8 F0
        integer p,q
        integer i,j

        real*8 suma
        real*8 h(Norb,Norb)                         ! calki jednoelektr.

        real*8 wp,wq                                ! wykladniki p i q
        real*8 cp,cq                                ! wsp. kontrakcji
        real*8 X                                    ! zm. pomocnicza
        real*8 Rpq,Rpq2, R2                         ! odl. miedzy centrami
        real*8 dx,dy,dz                             ! roznice wsp. cetrow
        real*8 xa,ya,za                             ! odl. od atomu
        real*8 rox,roy,roz

        integer lp,lq
        integer np,nq
        integer Ip,Iq
        integer Ilep, Ileq
        integer a                                   ! numer atomu
        real*8 Tpq,Vpq, S
        real*8 t                                    ! pomocnicza do F0
        real*8 c


        do p=1,Norb
          do q=p,Norb
            Ip = AO(p,1)                            ! znajdz indeksy
            Iq = AO(q,1)                            ! orbitali w bazie

            lp = Base(Ip,1)                         ! pobierz typ
            lq = Base(Iq,1)                         ! orbitalu (l)

            Ilep = Base(ip,2)                       ! ile prymitywnych
            Ileq = Base(iq,2)                       ! gaussow?

            ip = ip + 1                             ! przeskocz do
            iq = iq + 1                             ! danych prymitywow

c           ===========================
c           = odleglosc miedzy centr. =
c           ===========================

            dx = AO(p,3)-AO(q,3)
            dy = AO(p,4)-AO(q,4)
            dz = AO(p,5)-AO(q,5)
            Rpq2 = (dx*dx + dy*dy + dz*dz)

c           ==============
c           = licz calki =
c           ==============

            Tpq = 0.0
            Vpq = 0.0

            do i=ip, ip+Ilep -1
              wp = Base(i,1)
              cp = Base(i,2)

              do j=iq, iq+Ileq -1
                wq = Base(j,1)
                cq = Base(j,2)

c               =========================
c               = licz calke nakrywania =
c               =========================

                X = (2.*wp*wq)/(wp+wq)
                S = cp*cq* (2.*X/(wp+wq))**(0.75)*exp(-0.5*X*Rpq2)

c               =========================
c               = licz calke en. kinet. =
c               =========================

                Tpq = Tpq + X*(1.5-0.5*X*Rpq2)*S

c               =========================
c               = licz calke prz. jadr. =
c               =========================

                rox = (wp*AO(p,3)+wq*AO(q,3))/(wp+wq)
                roy = (wp*AO(p,4)+wq*AO(q,4))/(wp+wq)
                roz = (wp*AO(p,5)+wq*AO(q,5))/(wp+wq)

                do a=1,Nat
                  xa = GEOM(a,1) - rox
                  ya = GEOM(a,2) - roy
                  za = GEOM(a,3) - roz
                  R2 = (xa*xa + ya*ya + za*za)

                  t = (wp+wq)*R2
                  c = - 1.12837916709551257*GEOM(a,4)*sqrt(wp+wq)

                  Vpq= Vpq + c*F0(t)*S
                enddo
              enddo
            enddo
            h(p,q) = Vpq + Tpq
            h(q,p) = h(p,q)
          enddo
        enddo
      end subroutine

c=============================================
c== Oblicza macierz calek dwuelektr..<p|G|q>==
c==                                         ==
c==     G = macierz calek nakrywania (WY)   ==
c=============================================

      subroutine ComputeTwoElectron (G)
        IMPLICIT NONE
        integer Nat                                   ! ile atomow
        real*8 GEOM(1000,4)                           ! geometria (3 wsp. + Z)
        real*8 Base(1000,2)                           ! {wyk,C}
        integer BaseIndex(100)                        ! index w bazie dla danej
                                                      ! liczby atomowej

        integer Norb                                  ! ile orbitali AO
        real*8 AO(1000,5)                             ! lista orbitali AO
                                                      ! (index,m,x,y,z)
        COMMON /BASE/ Base, BaseIndex
        COMMON /GEOM/ GEOM,Nat
        COMMON /AO/ AO,Norb

        real*8 F0
        integer r,s,p,q
        integer i,j,k,l

        real*8 suma
        real*8 G(Norb,Norb,Norb,Norb)               ! calki dwuelektr.

        real*8 wr,ws,wp,wq                          ! wykladniki p i q
        real*8 cr,cs,cp,cq                          ! wsp. kontrakcji
        real*8 X                                    ! zm. pomocnicza
        real*8 Rpq2,Rrs2                            ! odl. miedzy centrami
        real*8 Ro2                                  ! odl. miedyz nowymi c.
        real*8 dx,dy,dz                             ! roznice wsp. cetrow
        real*8 xrs,yrs,zrs                          ! wsp. nowego centrum
        real*8 xpq,ypq,zpq                          ! wsp. nowego centrum
        real*8 Rpq,Rrs

        integer lr,ls,lp,lq
        integer nr,ns,np,nq
        integer Ir,Is,Ip,Iq
        integer Iler,Iles,Ilep, Ileq

        real*8 Srs, Spq
        real*8 t                                    ! pomocnicza do F0
        real*8 c
        real*8 PI

        parameter (PI=3.14159265358979324)


        do r=1,Norb
          do s=1,Norb
            do p=1,Norb
              do q=1,Norb
                Ir = AO(r,1)                        ! znajdz indeksy
                Is = AO(s,1)                        ! orbitali w bazie
                Ip = AO(p,1)                        ! znajdz indeksy
                Iq = AO(q,1)                        ! orbitali w bazie

                lr = Base(Ir,1)                     ! pobierz typ
                ls = Base(Is,1)                     ! orbitalu (l)
                lp = Base(Ip,1)                     ! pobierz typ
                lq = Base(Iq,1)                     ! orbitalu (l)

                Iler = Base(ir,2)                   ! ile prymitywnych
                Iles = Base(is,2)                   ! gaussow?
                Ilep = Base(ip,2)                   ! ile prymitywnych
                Ileq = Base(iq,2)                   ! gaussow?

                ir = ir + 1                         ! przeskocz do
                is = is + 1                         ! danych prymitywow
                ip = ip + 1                         ! przeskocz do
                iq = iq + 1                         ! danych prymitywow

c               ===========================
c               = licz odl. miedzy centr. =
c               ===========================

                dx = AO(p,3)-AO(q,3)
                dy = AO(p,4)-AO(q,4)
                dz = AO(p,5)-AO(q,5)
                Rpq2 = (dx*dx + dy*dy +dz*dz)

                dx = AO(r,3)-AO(s,3)
                dy = AO(r,4)-AO(s,4)
                dz = AO(r,5)-AO(s,5)
                Rrs2 = (dx*dx + dy*dy +dz*dz)

c               ==============
c               = licz calki =
c               ==============

                suma = 0.0

                do i=ir, ir+Iler -1
                  wr = Base(i,1)
                  cr = Base(i,2)

                  do j=is, is+Iles -1
                    ws = Base(j,1)
                    cs = Base(j,2)

                    do k=ip, ip+Ilep -1
                      wp = Base(k,1)
                      cp = Base(k,2)

                      do l=iq, iq+Ileq -1
                        wq = Base(l,1)
                        cq = Base(l,2)

c                       =========================
c                       = licz calki nakrywania =
c                       =========================

                        X = (2.*wp*wq)/(wp+wq)

                        Spq = cp*cq*(wp*wq)**(-0.75)*X**(1.5)
     &                      * exp(-0.5*X*Rpq2)


                        X = (2.*wr*ws)/(wr+ws)

                        Srs = cr*cs*(wr*ws)**(-0.75)*X**(1.5)
     &                      * exp(-0.5*X*Rrs2)

c                       =========================
c                       = licz calke dwuelektr. =
c                       =========================

                        xrs = (wr*AO(r,3)+ws*AO(s,3))/(wr+ws)
                        yrs = (wr*AO(r,4)+ws*AO(s,4))/(wr+ws)
                        zrs = (wr*AO(r,5)+ws*AO(s,5))/(wr+ws)

                        xpq = (wp*AO(p,3)+wq*AO(q,3))/(wp+wq)
                        ypq = (wp*AO(p,4)+wq*AO(q,4))/(wp+wq)
                        zpq = (wp*AO(p,5)+wq*AO(q,5))/(wp+wq)

                        Ro2 = (xrs-xpq)**2 + (yrs-ypq)**2 + (zrs-zpq)**2


                        X = (wr+ws)*(wp+wq)/(wr+ws+wp+wq)
                        c = 1.12837916709551257 * sqrt(X)          ! 2*sqr(X/pi)
                        t = X*Ro2

                        suma = suma + c*F0(t)*Srs*Spq
                      enddo
                    enddo
                  enddo
                enddo

                G(p,r,q,s) = suma
              enddo
            enddo
          enddo
        enddo
      end subroutine
