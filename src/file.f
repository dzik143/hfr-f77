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

c================================================
c==   Wczytuje geometrie molekuly z pliku      ==
c==                                            ==
c==     u_in = uchwyt do pliku (WE)            ==
c==      Nat = ile atomow (WY)                 ==
c==   GEOM = wspolrzedne atomow i l.at. (WY)   ==
c================================================

      subroutine ReadGeometry (u_in)
        IMPLICIT NONE

        integer u_in                          ! uchwyt do pliku
        integer Nat
        real*8 GEOM(1000,4)                   ! 3 wsp. i liczba at.

        character*10 tag
        real*8 x,y,z
        character*2  sym                      ! symbol pierwiastka
        integer i
        integer SymToZ
        integer occupied                      ! ile zajetych orbitali

        COMMON /GEOM/ GEOM,Nat
        COMMON /OCCUPIED/ occupied

c       ======================
c       = sprawdz tag 'GEOM' =
c       ======================

        read (u_in,'(A4)') tag

        if (tag.ne.'GEOM') then
          write (*,*) 'BLAD! Niepoprawna geometria!'
          stop
        endif

c       ====================
c       = czytaj geometrie =
c       ====================

        Nat = 0                                    ! zliczaj atomy
        occupied = 0                               ! licz orbitale
  100   read (u_in,'(A8)') tag
        if (tag.eq.'END_GEOM') goto 200

        backspace (u_in)
        read (u_in,*) sym, x,y,z

        Nat = Nat + 1

        GEOM(Nat,1) = x/0.5292
        GEOM(Nat,2) = y/0.5292
        GEOM(Nat,3) = z/0.5292
        GEOM(Nat,4) = SymToZ (sym)
        occupied = occupied + GEOM(Nat,4)          ! dodaj liczbe atom.
        goto 100                                   ! nastepny atom
  200   continue

c       ==========================
c       = sprwadz liczbe elektr. =
c       ==========================

        if (MOD(occupied,2).eq.1) then
          write (*,*) 'Blad! Nieparzysta liczba elektronow'
          stop
        endif
        occupied = occupied / 2
      end subroutine

c================================================
c==   Wczytuje baze funkcyjna z pliku          ==
c==                                            ==
c==     u_in = uchwyt do pliku z baza  (WE)    ==
c==     Base = baza funkcyjna w postaci: (WY)  ==
c==            (wyk,C)                         ==
c==              C = wsp. w rozwinieciu        ==
c==            wyk = wykladnik                 ==
c== BaseIndex = zakres funkcji w bazie dla (WY)==
c==             zadanej liczby atomowej        ==
c================================================

      subroutine ReadBase (BaseSymbol)
        IMPLICIT NONE

        real*8 Base(1000,2)
        integer BaseIndex(100)
        character*10 BaseSymbol

        integer u_in
        parameter (u_in=35)

        character bufor                            ! bufor do czytania
        integer La                                 ! liczba atomowa
        character*2 Sym                            ! symbol
        integer charge

        character*2 typ                            ! typ orbitalu
        integer l                                  ! poboczna licz. kw.
        integer IleFun                             ! ile funkcji
        real*8 cos

        integer IleOrb                             ! ile orbitali

        real*8 wyk(20), c(20),c2(20)               ! prym. GTO
        integer index                              ! index w bazie
        integer SymToZ                             ! symbol to l.at.
        integer i,k

        character*10 symbol                        ! symbol bazy
        character*10 plik                          ! nazwa pliku

        COMMON /BASE/ Base, BaseIndex

c       ============================
c       = szukaj bazy w bibliotece =
c       ============================

        open (u_in, file='basis.lib')
  500   read (u_in,*,end=510) symbol, plik
        if (symbol.eq.BaseSymbol) goto 520
        goto 500

  510   close (u_in)
        write (*,*) 'Blad! Brak bazy w biliotece!'
        stop

  520   close (u_in)
        open (u_in, file=plik)                      ! otworz plik z baza

        index = 1                                   ! index w bazie

c       ========================================
c       ====== kolejny pierwiastek  ============
c       ========================================

  100   read (u_in, '(A)',end=200) bufor            ! czytaj 1 znak
        if (bufor.ne.'!') then                      ! czy to komentarz?

          backspace (u_in)                          ! wroc na pocz. linijki
          read (u_in,*) sym, charge

          La = SymToZ (sym)                         ! pobierz liczbe at.
          BaseIndex (La) = index                    ! zachowaj polozenie
                                                    ! pierw. w bazie
          IleOrb = 0                                ! zliczaj orbitale
          index = index + 1                         ! zostaw miejsce
                                                    ! na liczbe orbitali
c         ====================
c         = kolejny orbital  =
c         ====================

  300     read (u_in,*) bufor
          if (bufor.eq.'*') then                    ! czy nastepny pierw?
            k = BaseIndex (La)                      ! poczatek pierwiastka
            Base(k,1) = La                          ! zapisz liczb. atm.
            Base(k,2) = IleOrb                      ! oraz ilosc orbiali
            goto 100                                ! czytaj nas. pierw.
          endif

          backspace (u_in)
          read (u_in,*) typ, IleFun, cos            ! kolejny orbital

c         ========================
c         = orbitale S i P razem =
c         ========================

          if (typ.eq.'SP') then
            IleOrb = IleOrb + 2
            do i=1,IleFun
              read (u_in,*) wyk(i), c(i),c2(i)
            enddo

            Base(index,1) = 0
            Base(index,2) = IleFun
            index = index + 1

            do i=1,IleFun
              Base(index,1) = wyk(i)
              Base(index,2) = c(i)
              index = index + 1
            enddo

            Base(index,1) = 1
            Base(index,2) = IleFun
            index = index + 1

            do i=1,IleFun
              Base(index,1) = wyk(i)
              Base(index,2) = c2(i)
              index = index + 1
            enddo

          else

c         =======================
c         = pozostale typy orb. =
c         =======================

            IleOrb = IleOrb + 1
            do i=1,IleFun
              read (u_in,*) wyk(i), c(i)
            enddo

            if (typ.eq.'S') l=0
            if (typ.eq.'P') l=1
            if (typ.eq.'D') l=2
            if (typ.eq.'F') l=3

            Base(index,1) = l
            Base(index,2) = IleFun
            index = index + 1

            do i=1,IleFun
              Base(index,1) = wyk(i)
              Base(index,2) = c(i)
              index = index + 1
            enddo
          endif
          goto 300                                 ! nastepny orbital
        endif

        goto 100
  200   continue

        close (u_in)
      end subroutine
