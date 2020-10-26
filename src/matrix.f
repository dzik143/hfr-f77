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

c12345678

c================================================
c= odwraca macierz A met. Gaussa-Jordana        =
c=       A = macierz wejsciowa  (WE)            =
c=   A_inv = macierz odwrotna do A (WY)         =
c=     wym = wymiar macierzy                    =
c================================================

      subroutine InvertMatrix (Ai, A, wym)
        IMPLICIT  NONE

        integer wym                        ! wymiar macierzy
        real*8 A                           ! macierz wejsciowa
        real*8 Ai                          ! macierz odwrotna

        DIMENSION A(wym,wym), Ai(wym,wym)  !

        integer i,j, k,l                   ! wskazuja na element (aij)
        real*8  skala


c       ==============================
c       == przygotuj macierz Ai jako =
c       == macierz jednostkowa       =
c       ==============================

        do i=1,wym
          do j=1,wym
            if (i.eq.j) then
              Ai(i,j)=1
            else
              Ai(i,j)=0
            endif
          enddo
        enddo


        do j=1, wym
c         ====================================
c         = gdy 0 na przek.to zamien wiersze =
c         ====================================

          if (A(j,j).eq.0) then

            if (j.eq.wym) then
              write (32,*) 'detS=0'
              stop
            endif

            k=j                              !
  100       if (A(k,j).ne.0) goto 200        ! szukaj wiersza bez 0
            k=k+1
            if (k.lt.wym) goto 100           ! czy koniec macierzy?
  200       continue                         !

            do l=1, wym                      ! zamien wiersze
              call xchg (A(j,l),A(k,l))      !    j <=> l
              call xchg (Ai(j,l),Ai(k,l))    !
            enddo
          endif

c         ====================================
c         = gdy dalej 0 to zamien kolumny    =
c         ====================================
          if (A(j,j).eq.0) then

            k=j                              !
  300       if (A(j,k).ne.0) goto 400        ! szukaj kolumny bez 0
            k=k+1
            if (k.lt.wym) goto 300           ! czy koniec macierzy?
  400       continue                         !

            do l=1, wym                      ! zamien kolumny
              call xchg (A(l,j),A(l,k))      !    j <=> l
              call xchg (Ai(l,j),Ai(l,k))    !
            enddo
          endif

c         ====================================
c         == przeskaluj wiersz zeby dostac 1 =
c         == na przekatnej                   =
c         ====================================

          skala = 1/A(j,j)

          do l=1, wym
            A(j,l) = A(j,l)*skala
            Ai(j,l) = Ai(j,l)*skala
          enddo

          do i=1, wym
c           =================================
c           == likwiduj wszystkie elementy  =
c           == w kolumnie procz elementu j  =
c           =================================

            if (i.ne.j) then
              skala = A(i,j)/A(j,j)
              do l=1,wym
                A(i,l)=A(i,l)-A(j,l)* skala
                Ai(i,l)=Ai(i,l)-Ai(j,l)* skala
              enddo
            endif
          enddo
        enddo
      end subroutine

c========================================================
c=      procedura mnozy 2 macierze                      =
c=            C = AxB                                   =
c=        w1,k1 = wymiar macierzy A                     =
c=        w2,k2 = wymiar macierzy B                     =
c========================================================

      subroutine mxm (C, A,B, w1,k1, w2,k2)
        IMPLICIT NONE

        integer wym, i,j,k
        integer w1,k1, w2,k2

        real*8 suma
        real*8 A,B,C

        dimension A(w1,k1)
        dimension B(w2,k2)
        dimension C(w1,k2)

c       ===================
c       = sprawdz wymiary =
c       ===================

        if (k1.ne.w2) then
          write (*,*)'Niepoprawny rozmiar!'
          stop
        endif

c       ========
c       = mnoz =
c       ========

        do j=1, k2      ! kolumny (=>B)
          do i=1, w1    ! wiersz  (-A)
            suma=0

            do k=1, k1
              suma = suma + A(i,k)*B(k,j)
            enddo

            C(i,j)=suma
          enddo
        enddo

        k1 = k2
      end subroutine

c==========================================
c=  transponuje macierz A                 =
c==========================================

      subroutine TranMatrix (At, A, wym)
        integer wym

        real*8 At(wym,wym)
        real*8 A(wym,wym)

        integer i,j

        do i=1,wym
          do j=1,wym
            At(i,j) = A(j,i)
          enddo
        enddo
      end subroutine

c====================================================
c=  Podnosi macierz do wybranej potegi              =
c=        T = macierz wejsciowa                     =
c=      wym = wymiar macierzy wejsciowej            =
c====================================================

      subroutine PowerMatrix (T, wykladnik, wym)

         IMPLICIT NONE
         integer  wym                              ! wymiar macierzy
         integer  i,j                              !

          real*8  wykladnik
          real*8  T,S
       DIMENSION  T(wym,wym), S(wym,wym)           ! bufory robocze
          real*8  buf(wym,wym), S_t(wym,wym)

          real*8  u,v,lambda                       !
          real*8  sinx, cosx, sinxcos, sin2, cos2  ! zmienne pomocnicze
          real*8  Tpq,Tqq,Tpp,Tpj,Tip,Tmax, Sip
          real*8  sgn                              ! funkcja do znaku

          real*8  thre                             ! prog czulosci

          real*8  evalue(wym), evec(wym,wym)

          call DiagMatrix (T, evalue, evec, wym,0)

          do i=1, wym
           do j=1, wym
            if (i.eq.j) then
              if (T(i,j).le.0.0) then
                write (*,*) 'BLAD!!! war. wlasne S<0!'
              else
                T(i,j) = T(i,j) ** wykladnik
              endif
            else
              T(i,j)=0.0
            endif
          enddo
         enddo

         call TranMatrix (S_t, evec, wym)
         call MxM (buf, evec, T, wym,wym,wym,wym)
         call MxM (T, buf, S_t, wym,wym,wym,wym)

         end subroutine


c==========================================
c=  procedura zamienia zmienne A,B        =
c==========================================

      subroutine xchg (A,B)
        real*8  A,B, Temp

        Temp=A
        A=B
        B=Temp
      end subroutine

c ================================================
c ==        laduje z pliku macierz              ==
c ==                                            ==
c ==   u_in = uchwyt do otwartego pliku (WE)    ==
c ==      A = macierz docelowa  (WY)            ==
c ==   wier = liczba wierszy  (WY)              ==
c ==    kol = liczba kolumn   (WY)              ==
c=================================================
       subroutine LoadMatrix (u_in, A, wier, kol)
              IMPLICIT NONE
              integer  wier, kol               ! wymiar macierzy
              integer  u_in                    ! uchwyt do pliku
               real*8  A
             dimension A(wier,kol)             ! macierz docelowa
              integer  i,j                     ! liczniki

              do i=1, wier
                read (u_in,*) (A(i,j), j=1,kol)
              enddo

       end subroutine

c==================================================
c=        Wypisuje macierz                        =
c=   u_out = ID urzadzenia wyjsciowego            =
c=  Matrix = macierz do wypisania                 =
c=  wier,kol = wymiary macierzy                   =
c==================================================
      subroutine WriteMatrix (u_out,A,wier,kol)
          integer  wier,kol
          real*8   A
       dimension   A(wier,kol)

          integer  u_out
          integer  i,j

          do i=1,wier
             write (u_out,'(12f12.6)') (A(i,j),j=1,kol)
          enddo
      end Subroutine

c=================================================
c    procedura kopiuje dane pomiedzy buforami    =
c        D = bufor docelowy (destiny)            =
c        S = bufor zrodlowy (source)             =
c      LEN = ilosc komorek do przeniesienia      =
c=================================================
       subroutine copy (D,S, LEN)
         IMPLICIT NONE
          real*8  D, S
         integer  i, LEN
       DIMENSION  D(*), S(*)

         do i=1, LEN
            D(i)=S(i)
         enddo

      end subroutine


c=========================================
c==    funkcja zwraca znak argumentu    ==
c=========================================
      real*8 function sgn(x)
         real*8  x

         if (x.lt.0.) then
           sgn =-1.
         else
           sgn = 1.
         endif
      end function

c =============================================
c ==         sortuje tablice TAB()           ==
c ==                                         ==
c ==     TAB = tablica do sortowania         ==
c ==     len = dlugosc tablicy               ==
c ==   index = nowa kolejnosc elementow      ==
c ==  direct = kierunek sortowania (1=rosn.) ==
c =============================================
       subroutine SortIndex (TAB, len, index, direct)
         IMPLICIT NONE
          real*8  TAB(1000), a0, tmp           ! a0 = "srodek podzialu"
         integer  len, i,j                     ! left, right
         integer  l,r
         integer  left(1000), right(1000)      ! kolejka zadan
         integer  sp                           ! wierzcholek kolejki
         integer  index(1000),tmp_i            ! do zapisu kolejnosci
         integer  direct                       ! kierunek sortowania

         sp = 1                                ! dodaj pierwsze zadanie
         left(1) = 1                           ! do kolejki:
         right(1) = len                        ! sortuj przedzial (a,b)


         do i=1, len                           ! poczatkowa kolejnosc
           index(i) = i
         enddo

c ====================================
c ====  pobierz zadanie ze stosu   ===
c ====================================
  100    l = left(sp)                          ! pobierz przedzial
         r = right(sp)                         ! do sortowania
         sp = sp-1                             ! usun zadanie z kolejki
         a0 = TAB((l+r)/2)                     ! element odniesienia a0
c         write (*,*)'a0=', a0

         i = l-1
         j = r+1
c ==================================
c ====    kierunek rosnacy    ======
c ==================================
  200    if (direct.eq.-1) goto 500

  300    i=i+1
         if (TAB(i).lt.a0) goto 300            ! szukaj a>=a0 (po lewej)
  400    j=j-1
         if (TAB(j).gt.a0) goto 400            ! szukaj a<=a0 (po prawej)
         goto 700

c ==================================
c ====    kierunek malejacy    =====
c ==================================

  500    i=i+1
         if (TAB(i).gt.a0) goto 500            ! szukaj a<=a0 (po lewej)
  600    j=j-1
         if (TAB(j).lt.a0) goto 600            ! szukaj a>=a0 (po prawej)

c ==================================
c ==   zamien miejscami i<->j    ===
c ==================================
  700      if (i.lt.j) then
c           write (*,*)'swap=', i,j
           tmp = TAB(i)                        ! zamien
           TAB(i) = TAB(j)                     ! miejscami
           TAB(j) = tmp                        ! (i)<->(j)

           tmp_i = index(i)                    ! zapisz nowa
           index(i) = index(j)                 ! kolejnosc
           index(j) = tmp_i                    !

           goto 200                            ! szukaj nast. pary (i,j)
         endif


c =================================
c ==    generuj nowe zadania     ==
c =================================
         if ((r-l).gt.1) then                  ! czy juz 1 element?

            if (j.gt.l) then                   !
c             write (*,*) 'add task:', l,j
              sp = sp+1                        ! lewy podprzedzial:
              left(sp) = l                     !   (l, j)
              right(sp) = j                    !
            endif

            if (j.lt.r) then
c             write (*,*) 'add task:', j+1, r
              sp = sp+1                        ! prawy podprzedzial:
              left(sp) = j+1                   !   (j+1, p)
              right(sp) = r                    !
            endif

         endif

         if (sp.gt.0) goto 100                 ! nastepne zadanie

       end subroutine


c====================================================
c=  Diagonalizuje macierz metoda Jacobiego          =
c=        T = macierz wejsciowa                     =
c=      wym = wymiar macierzy wejsciowej            =
c=   evalue = bufor na wartosci wlasne              =
c=        S = bufor na wektory wlasne               =
c====================================================
      subroutine DiagMatrix (T, evalue, evec, wym,sort)

         IMPLICIT NONE
         integer  wym                              ! wymiar macierzy
         integer  p,q                              ! (p,q) = max
c         integer  pq,qp,pp,qq                     ! indeksy do tablic
c         integer  pj,qj,iq,ip                     !
         integer  i,j                              !
         integer  index(1000)                      ! do sortowania wektorow
         integer  sort

          real*8  T,S
       DIMENSION  T(wym,wym), S(wym,wym)           ! bufory robocze
          real*8  u,v,lambda                       !
          real*8  sinx, cosx, sinxcos, sin2, cos2  ! zmienne pomocnicze
          real*8  Tpq,Tqq,Tpp,Tpj,Tip,Tmax,Sip,Siq
          real*8  sgn                              ! funkcja do znaku
          real*8  evalue(100), evec                !
       DIMENSION  evec(wym,wym)                    ! wektory wlasne
          real*8  thre                             ! prog czulosci

          real*8  teta, tau,te

          integer iter

c    ======================
c    ===    S0 = (1)  =====
c    ======================

        do i=1, wym
           do j=1, wym
             if (i.eq.j) then
               S(i,j) = 1.                         ! jedynki na
             else                                  ! na przekatnej
               S(i,j) = 0.
             endif
           enddo
        enddo


       thre=0.000000001
c    ==============================
c    ======= szukaj (p,q) max =====
c    ==============================
 100    tmax=0
        do i=1, wym
          do j=i+1, wym                            ! omin przekatna
                if (abs(T(i,j)).gt.tmax) then
                    tmax=abs(T(i,j))               ! nowe maksimum
                    p=i
                    q=j
               endif
          enddo
        enddo

c   =============================
c   ==   zmienne pomocnicze    ==
c   =============================
          teta = (T(q,q)-T(p,p))/(2.*T(p,q))
            te = sgn(teta)/(abs(teta)+sqrt(teta*teta+1.))
          cosx = 1. / sqrt(te*te + 1.)
          sinx = te * cosx
           tau = sinx/(1.+cosx)

       sinxcos = sinx*cosx                           !
          sin2 = sinx*sinx                           !
          cos2 = cosx*cosx                           !

          Tpp = cos2*T(p,p)+sin2*T(q,q)-2.*sinxcos*T(p,q)
          Tqq = sin2*T(p,p)+cos2*T(q,q)+2.*sinxcos*T(p,q)
          Tpq = (cos2-sin2)*T(p,q)+sinxcos*(T(p,p)-T(q,q))

c      =================================
c      == jedz po wierszu (P,j),(Q,j) ==
c      =================================

         do j=1, wym
            Tpj = T(p,j)
            T(p,j) = Tpj*cosx - T(q,j)*sinx
            T(q,j) = Tpj*sinx + T(q,j)*cosx
         enddo

c      ==================================
c      === jedz po kolumn.(i,P),(i,Q) ===
c      ==================================
         do i=1, wym
            Tip = T(i,p)
            T(i,p) = cosx*Tip - sinx*T(i,q)
            T(i,q) = sinx*Tip + cosx*T(i,q)

            Sip = S(i,p)
            Siq = S(i,q)
            S(i,p) = Sip*cosx - Siq*sinx            ! wektory
            S(i,q) = Sip*sinx + Siq*cosx            ! wlasne
         enddo

         T(p,q) = Tpq
         T(q,p) = Tpq
         T(p,p) = Tpp
         T(q,q) = Tqq
                                                     ! sprawdz
         if (tmax.gt.thre) goto 100                  ! czy koniec

c      ==================================
c      ==  wypelnij tablice na wyjscie  =
c      ==================================
         do i=1, wym
           evalue(i) = T(i,i)
         enddo

c      ==================================
c      ==         sortuj dane          ==
c      ==================================
         if (sort.eq.1) then
           call SortIndex (evalue, wym, index,-1)       ! sortuj wartosci wl.
c
           do i=1, wym                                 ! wypelnij macierz
             j = index(i)                              ! wektorow w takiej
             call copy (evec(1,i), S(1,j), wym)        ! kolejnosci jak wart.
           enddo                                       ! wlasne

         else

           do i=1,wym
            do j=1,wym
              evec(i,j) = S(i,j)
            enddo
           enddo
         endif

      end subroutine

c====================================================
c=  Diagonalizuje macierz metoda Jacobiego          =
c=        T = macierz wejsciowa                     =
c=      wym = wymiar macierzy wejsciowej            =
c=   evalue = bufor na wartosci wlasne              =
c=        S = bufor na wektory wlasne               =
c=    index = kolejnosc sortowania (WY)             =
c====================================================

      subroutine DiagMatrixSort (T, evalue, S, wym, index)

        IMPLICIT NONE
        integer wym                             ! wymiar macierzy
        integer p,q                             ! (p,q) = max
        integer i,j                             !
        integer index(wym)                      ! do sortowania wektorow
        integer sort

        real*8 T(wym,wym),S(wym,wym)            ! bufory robocze
        real*8 u,v,lambda                       !
        real*8 sinx, cosx, sinxcos, sin2, cos2  ! zmienne pomocnicze
        real*8 Tpq,Tqq,Tpp,Tpj,Tip,Tmax,Sip,Siq
        real*8 sgn                              ! funkcja do znaku
        real*8 evalue(wym), evec(wym,wym)       !
        real*8 evalue2(wym)
        real*8 thre                             ! prog czulosci

        integer iter
        real*8  teta, tau,te

c    ======================
c    ===    S0 = (1)  =====
c    ======================

        do i=1, wym
          do j=1, wym
            if (i.eq.j) then
              S(i,j) = 1.                         ! jedynki na
             else                                 ! na przekatnej
               S(i,j) = 0.
             endif
           enddo
        enddo

        thre= 1.0D-8

c    ==============================
c    ======= szukaj (p,q) max =====
c    ==============================
 100    tmax=0
        do i=1, wym
          do j=i+1, wym                        ! omin przekatna
            if (abs(T(i,j)).gt.tmax) then
              tmax=abs(T(i,j))                 ! nowe maksimum
              p=i
              q=j
            endif
          enddo
        enddo

c       =============================
c       ==   zmienne pomocnicze    ==
c       =============================

        teta = (T(q,q)-T(p,p))/(2.*T(p,q))
        te   = sgn(teta)/(abs(teta)+sqrt(teta*teta+1.))
        cosx = 1. / sqrt(te*te + 1.)
        sinx = te * cosx
        tau  = sinx/(1.+cosx)

        sinxcos = sinx*cosx                           !
        sin2    = sinx*sinx                           !
        cos2    = cosx*cosx                           !


        Tpp = cos2*T(p,p)+sin2*T(q,q)-2.*sinxcos*T(p,q)
        Tqq = sin2*T(p,p)+cos2*T(q,q)+2.*sinxcos*T(p,q)
        Tpq = (cos2-sin2)*T(p,q)+sinxcos*(T(p,p)-T(q,q))

c       =================================
c       == jedz po wierszu (P,j),(Q,j) ==
c       =================================

        do j=1, wym
          Tpj = T(p,j)
          T(p,j) = Tpj*cosx - T(q,j)*sinx          !
          T(q,j) = Tpj*sinx + T(q,j)*cosx          !
        enddo

c       ==================================
c       === jedz po kolumn.(i,P),(i,Q) ===
c       ==================================

        do i=1, wym
          Tip = T(i,p)                             !
          T(i,p) = cosx*Tip - sinx*T(i,q)
          T(i,q) = sinx*Tip + cosx*T(i,q)

          Sip = S(i,p)                             !
          Siq = S(i,q)
          S(i,p) = Sip*cosx - Siq*sinx             ! wektory
          S(i,q) = Sip*sinx + Siq*cosx             ! wlasne
        enddo

        T(p,q) = Tpq
        T(q,p) = Tpq
        T(p,p) = Tpp
        T(q,q) = Tqq
                                                   ! sprawdz
        if (tmax.gt.thre) goto 100                 ! czy koniec

c       ==================================
c       ==  wypelnij tablice na wyjscie  =
c       ==================================
        do i=1, wym
          evalue(i) = T(i,i)
          evalue2(i) = evalue(i)
        enddo

c       ==================================
c       ==         sortuj dane          ==
c       ==================================

        call SortIndex (evalue2, wym, index,1)        ! sortuj wartosci wl.

      end subroutine

c=========================================
c= funkcja zwraca wartosc wyznacznika    =
c= macierzy A o wymairze wym             =
c=========================================
      real*8 function det (A,wym)
        integer wym
        real*8 A(wym,wym)
        integer n, i,j,l
        real*8 suma, skala
        real*8 znak           !do zamiany wierszy

        znak=1

        do j=1, wym-1

c       =========================
c       = likwiduj j-ta kolumne =
c       =========================

          do i=(j+1), wym

c           ===============================
c           = sprawdz czy 0 na przekatnej =
c           ===============================

            if (A(j,j).eq.0) then
              if (znak.eq.1) then
               znak=-1
              else
                znak=1
              endif

              do l=1, wym
               call xchg (A(j,l),A(j+1,l))
              enddo
            endif

c           ===================
c           = odejmij wiersze =
c           ===================

            skala = A(i,j)/A(j,j)

            do l=1,wym
              A(i,l)=A(i,l)-A(j,l)* skala
            enddo
          enddo
        enddo

        det =1
        do i=1, wym
          det = det*A(i,i)
        enddo
        det = det*znak

      end function
