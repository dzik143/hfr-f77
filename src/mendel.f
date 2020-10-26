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

c=====================================================
c==      = Wczytuje liste symboli =                 ==
c==                                                 ==
c==  Mendel(100) - lista symboli wg. liczby Z (WY)  ==
c=====================================================

      subroutine LoadMendel ()
        IMPLICIT NONE
        character*2 Mendel(100)              ! symbole pierwiastkow
        integer u_in                         ! uchwyt do pliku
        parameter (u_in=35)

        integer La                           ! liczba atomowa
        character*2 Sym                      ! symbol

        COMMON /MENDEL/ Mendel

c       ==========================
c       = wczytaj uklad okresowy =
c       ==========================

        open (u_in, file='sym.dat')
  100   read (u_in,*,end=110) La, sym
        Mendel(La) = sym
        goto 100
  110   close (u_in)

      end subroutine

c=====================================================
c==  Zwraca liczbe atomowa na podstawie symbolu     ==
c==                                                 ==
c==  Sym - 2 znakowy symbol (WE)                    ==
c=====================================================

      integer function SymToZ (sym)
        IMPLICIT NONE
        character*2 Sym                          ! szukany symbol
        character*2 Mendel(100)                  ! symbole pierwiastkow
        integer i

        COMMON /MENDEL/ Mendel

        do i=1,10
          if (sym.eq.Mendel(i)) then
            SymToZ = i
            return
          endif
        enddo
        Write (*,*) 'Blad! Brak pierwiastka w tablicy!'
        stop
      end function
