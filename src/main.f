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

c**************************************************
c*                                                *
c*            Hartree-Fock-Rothan                 *
c*          -----------------------               *
c*                                                *
c**************************************************

      program hfr
        IMPLICIT NONE

c       =====================
c       = Zmienne globalne  =
c       =====================

        character*2 Mendel(100)              ! symbole pierwiastkow

        integer Nat                          ! ile atomow
        real*8 GEOM(1000,4)                  ! geometria (3 wsp. + Z)

        real*8 Base(1000,2)                  ! {wyk,C}
        integer BaseIndex(100)               ! index w bazie dla danej
                                             ! liczby atomowej

        integer Norb                         ! ile orbitali AO
        real*8 AO(1000,5)                    ! lista orbitali AO
                                             ! (index,m,x,y,z)

        integer  occupied                    ! ile zajetych orbitali

        COMMON /occupied/ occupied
        COMMON /MENDEL/ Mendel
        COMMON /GEOM/ GEOM,Nat
        COMMON /BASE/ Base, BaseIndex
        COMMON /AO/ AO,Norb

c       =====================
c       == Zmienne lokalne  =
c       =====================

        character*10 BaseSymbol                ! symbol bazy

        real*8 S(100*100)                      ! macierz nakrywania
        real*8 C(100*100)                      ! wspolczyniki LCAO
        real*8 h(100*100)                      ! calki jednoelektr.
        real*8 G(50*50*50*50)                  ! calki dwuelektr.

        integer u_in, u_out
        parameter (u_in=31, u_out=32)

        character*20 nazwa                     ! nazwa projektu
        integer nazwaLEN                       ! dlugosc nazwy

        character temp
        integer i,j

        integer SymToZ
        real*8 Etotal
        real*8 SCF

c       ========================
c       == sprawdz parametry  ==
c       ========================

        call GetArg (1,nazwa)

        if (nazwa(1:1).eq.' ') then
          write (*,*) 'Blad: Nie podano nazwy zbioru wejsciowego'
          write (*,*) 'skladnia: hfr <nazwa zbioru>'
          stop
        endif

c       ================================
c       ====== zlicz dlugosc nazwy =====
c       ================================

        i=0                                           !
 200 0  i=i+1                                         !
        if (nazwa(i:i).ne.' ') goto 200               ! zlicz dlugosc nazwy
        nazwaLEN = i-1;                               !

        call LoadMendel()                             ! wczytaj symbole

c       ============================
c       ===== czytaj geometrie =====
c       ============================

        open (u_in, file=nazwa(1:nazwaLEN)//'.inp')
        call ReadGeometry (u_in)                      ! wczytaj geometrie

c       ==================================
c       ===== wczytaj baze funkcyjna  ====
c       ==================================

        read (u_in, *) BaseSymbol
        call ReadBase (BaseSymbol)
        close (u_in)

c       ==================================
c       ===== Tworz liste orbitali AO   ==
c       ==================================

        call MakeListAO ()

        call ComputeOverlap (S)
        call ComputeOneElectron (h)
        call ComputeTwoElectron (G)

        open (u_out, file=nazwa(1:nazwaLEN)//'.out')

        Etotal = SCF (C,h,G,S)
        write (u_out,*) 'Etotal=', Etotal
        write (u_out,*)

        write (u_out,*) '--- Funkcje falowe ---'      ! nieposortowane!
        call WriteMatrix (u_out,C,Norb,Norb)
        write (u_out,*)

        close (u_out)

      end

c==============================================
c==  SCF                                     ==
c==                                          ==
c==  C - wektory wlasne (WY)                 ==
c==  E - energie orbitalne (WY)              ==
c==  h - macierz calek 1-elektr. (WE)        ==
c==  G - macierz calek 2-elektr. (WE)        ==
c==  S - macierz calek nakrywania (WE)       ==
c==                                          ==
c==  zwraca:  Etotal                         ==
c==============================================

      real*8 function SCF (C,h,G,S)
        IMPLICIT NONE
        integer   Norb                                ! ile orbitali AO
        real*8   AO(1000,5)                           ! lista orbitali AO
                                                      ! (index,m,x,y,z)
        COMMON /AO/ AO,Norb

        real*8  C(Norb,Norb)
        real*8  h(Norb,Norb)
        real*8  G(Norb,Norb,Norb,Norb)
        real*8  S(Norb,Norb)

        real*8  P(Norb,Norb)
        real*8  F(Norb,Norb)                          ! macierz Focka

        real*8  FxS(Norb,Norb)
        real*8  SFS(Norb,Norb)
        real*8  Cij_(Norb,Norb)
        real*8  Cij(Norb,Norb)
        real*8  E(1000), Eold,E0, deltaE
        real*8  THRE
        real*8  Etotal                                ! E0 = NuclearRepul

        real*8  NuclearRepulsion
        integer  index(Norb)                          ! kolejnosc orbitali
        integer  iter, i,j


c     ============================
c     == na poczatku wypelnij 0  =
c     ============================
        do i=1,Norb
         index(i) = i
         do j=1,Norb
           C(i,j)=0.0
         enddo
        enddo

        THRE = 1.0D-7                                ! prog dokladnosci
        call PowerMatrix (S, -0.5D0, Norb)           ! S = S^(-1/2)

        iter = 1                                     ! zliczaj iteracje
 100    call SetUpP (P,C, index)
        call SetUpFock (F,h,G,P)

c   =================================
c   == rozwiazuje rownanie wiekowe ==
c   ==     F x C = S x C x E       ==
c   =================================

        call mxm (FxS, F,S, Norb,Norb,Norb,Norb)     ! FS =  F x S^-1/2
        call mxm (SFS, S,FxS, Norb,Norb,Norb,Norb)   ! SFS=S^-1/2 x H x S^-1/2
        call DiagMatrixSort (SFS, E, Cij_,Norb,index)
        call mxm (C, S,Cij_, Norb,Norb,Norb,Norb)    !

c   =================================
c   === licz energie calkowita ======
c   =================================

        call SetUpP(P,C, index)
        call SetUpFock (F,h,G,P)

        E0 = 0.0
        do i=1,Norb
         do j=1,Norb
           E0 = E0 + P(j,i)*(h(i,j)+F(i,j))
         enddo
        enddo
        E0 = 0.5*E0

c   ====================================
c   === sprawdz czy juz koniec?  =======
c   ====================================
        Etotal = E0 + NuclearRepulsion()
        deltaE = Eold-Etotal

        if (iter.eq.1) then
          Write (*,*)'iter=',iter,'   Etotal=',Etotal
        else
          Write (*,*)'iter=',iter,'   Etotal=',Etotal,' delta=',deltaE
        endif

        Eold = Etotal
        iter = iter + 1
        if (deltaE.gt.THRE) goto 100

        SCF = E0 + NuclearRepulsion ()                     ! zwroc Etotal
      end function


c====================================================
c== Oblicza macierz rzedow wiazan                  ==
c==                                                ==
c==      P - macierz rzedow wiazan  (WY)           ==
c==      C - wspolczynniki LCAO (WE)               ==
c==  index - kolejnosc zapelniania orbitali MO (WE)==
c====================================================

      subroutine SetUpP (P,C, index)
        IMPLICIT NONE
          integer  Norb                        ! ile orbitali AO
          real*8   AO(1000,5)                  ! lista orbitali AO
                                               ! (index,m,x,y,z)
          integer  occupied                    ! ile zajetych orbitali

          COMMON /AO/ AO,Norb
          COMMON /occupied/ occupied

          integer index(Norb)                  ! kolejnosc orbitali MO
          real*8  P(Norb,Norb)                 ! macierz rz. wiazan
          real*8  C(Norb,Norb)                 ! wspolczynniki LCAO

          integer r,s,i,j
          real*8  suma

          do r=1,Norb
           do s=1,Norb
             suma = 0.0
             do i=1,occupied
               j = index(i)
               suma = suma + C(r,j)*C(s,j)
             enddo
             P(s,r) = 2.* suma
           enddo
          enddo
      end subroutine

c====================================================
c== Oblicza macierz Focka                          ==
c==                                                ==
c==   Fock - macierz Focka (WY)                    ==
c==     AO - lista orbitali (WE)                   ==
c==      P - macierz rzedow wiazan (WE)            ==
c==   Norb - ile orbitali (WE)                     ==
c====================================================

      subroutine SetUpFock (F,h,G,Pr)
        IMPLICIT NONE
        integer  Norb                                 ! ile orbitali AO
        real*8   AO(1000,5)                           ! lista orbitali AO
                                                      ! (index,m,x,y,z)
        COMMON /AO/ AO,Norb

         real*8  F(Norb,Norb)
         real*8  h(Norb,Norb)
         real*8  G(Norb,Norb,Norb,Norb)
         real*8  Pr(Norb,Norb)
         real*8  suma

        integer  r,s,p,q

        do r=1,Norb
         do s=r, Norb
          suma = h(r,s)

          do p=1,Norb
           do q=1,Norb
            suma = suma + Pr(q,p)*(G(r,p,s,q)-0.5*G(r,p,q,s))
           enddo
          enddo

          F(r,s) = suma
          F(s,r) = suma

         enddo
        enddo

      end subroutine

c====================================================
c== Generuje liste orbitali na podstawie           ==
c==  geometrii molekuly                            ==
c==                                                ==
c==    AO = tablica orbitali (WY)                  ==
c==       podana w postaci (index, m, x,y,z)       ==
c==   index = index orbitalu w bazie               ==
c==       m = magnetyczna liczba kwantowa          ==
c==   x,y,z = wsp. centrum                         ==
c==                                                ==
c==  Norb = liczba orbitali (WY)                   ==
c==   Nat = liczba atomow (WE)                     ==
c==  GEOM = Geometria molekuly (3N + L.at) (WE)    ==
c====================================================

      subroutine MakeListAO ()
          IMPLICIT NONE

          integer Norb, Nat
           real*8 GEOM(1000,4)                   ! 3N wsp. + La
           real*8 AO(1000,5)                     ! lista orbitali (WY)
           real*8 Base(1000,2)                   ! baza funkcyjna
          integer BaseIndex(100)                 ! rozdzial funkcji
                                                 ! do pierwiastkow

          integer La                             ! liczba atomowa
           real*8 x,y,z                          ! wsp. atomow
          integer OrbInBase                      ! ile orb. ma pierw.
          integer Index                          ! index orb. w bazie
          integer i,j,k
          integer l,m                            ! liczby kwantowe

           COMMON /GEOM/ GEOM,Nat
           COMMON /BASE/ Base, BaseIndex
           COMMON /AO/ AO,Norb

          Norb = 0                               ! zliczaj orbitale AO
          do i=1, Nat                            ! jedz po atomach
            x = GEOM(i,1)                        !
            y = GEOM(i,2)                        !  wspolrzedne atomu
            z = GEOM(i,3)                        !  i jego liczba
            La = GEOM(i,4)                       !  atomowa

            index = BaseIndex (La)               ! znajdz pierw. w bazie

            OrbInBase = Base(index,2)            ! ile orbitali do wczyt.?
            index = index + 1                    ! przesun sie do orbitali

            do j=1,OrbInBase
              l = Base(index,1)                  ! typ orbitalu

              do m=-l,l
                Norb = Norb + 1                  ! nastepny orbital AO

                AO(Norb,1) = index               ! polozenie w bazie
                AO(Norb,2) = m                   ! liczba magnetyczna
                AO(Norb,3) = x                   !
                AO(Norb,4) = y                   ! wspolrzedne
                AO(Norb,5) = z                   ! centrum
              enddo

              index = index + Base(index,2) + 1  ! przesun sie do nast. orb.
            enddo
          enddo
      end subroutine

c================================================
c==  Oblicza energie odpychania jadrowego      ==
c==                                            ==
c================================================

       real*8 function NuclearRepulsion ()
            IMPLICIT NONE
            integer  Nat                          ! ile atomow
             real*8  GEOM(1000,4)                 ! geometria (3 wsp. + Z)

             COMMON /GEOM/ GEOM,Nat

             real*8  Ejj
            integer  i,j
            integer  Za,Zb
             real*8  dx,dy,dz
             real*8  Rab

            real*8  k
          parameter (k=0.52917693824475)

          Ejj = 0.0
          do i=1,Nat
            do j=i+1,Nat

               dx = GEOM(i,1) - GEOM(j,1)
               dy = GEOM(i,2) - GEOM(j,2)
               dz = GEOM(i,3) - GEOM(j,3)
               Rab = sqrt (dx*dx + dy*dy + dz*dz)

               Za = GEOM(i,4)
               Zb = GEOM(j,4)

                Ejj = Ejj + Za*Zb/Rab
              enddo
          enddo
c          Ejj = Ejj !* k !* 27.2097
          NuclearRepulsion = Ejj !* k
       end function
