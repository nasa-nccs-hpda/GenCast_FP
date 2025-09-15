! ETA2XPRS (ECMWF method)
! ----------------------
!
! Python extension via f2py to vertically interpolate from model levels (hybrid sigma pressure corrdinates)
! to constant pressure levels, with an option to extrapolate below the surface. This extrapolation is
! based on the so-called "ECMWF method". A similar interface is available in NCL.
!
! Note: Rountine eta2xprs is the entry point for the python interface. The core of the calculation is 
! performed in the legacy routine VINTH2PECMWF.
!
! Arlindo da Silva, April 2025.
!
      
      subroutine eta2xprs ( im, jm, km, kp, A_prs, A_eta, &
                           ak, bk, ps, ts, phis, plevs,        &
                           xflag, vflag, method, undef )

      implicit NONE

      integer,      intent(in)  :: im, jm, km, kp   ! dimensions
      real(kind=8), intent(in)  :: A_eta(im,jm,km)  ! Input array

                                                    ! p = ak + bk * psfc
      real(kind=8), INTENT(IN)  :: ak(km)           !    ak [Pa]
      real(kind=8), INTENT(IN)  :: bk(km)           !    bk, non-dimentional
                                                    ! NOTE: if input A_eta is mid-layer then
                                                    ! ak/bk must be midlayer as well.
      
      real(kind=8), intent(in)  :: ps(im,jm)        ! surface pressure [Pa]
      real(kind=8), intent(in)  :: ts(im,jm)        ! temperature next to the ground [K]
      real(kind=8), intent(in)  :: phis(im,jm)      ! surface geopotential [m2 s-2]
      real(kind=8), intent(in)  :: plevs(kp)        ! pressure levels to interpolate to [hPa]

      integer, intent(in)  :: xflag ! Extrapolation flag: 
                                    !        0 - do not extrapolate
                                    !        1 - extrapolate using ECMWF method
      
      integer, intent(in)  :: vflag ! Variable flag: 
                                    !        -1  means geopotential height (Z)
                                    !        +1  means temperature (T)
                                    !         0  any other variable
      
      integer, intent(in)  :: method ! Method for vertical interpolation
                                     ! 1 - LINEAR, 2 - LOG, 3 - LOG LOG

      real(kind=8), intent(in)  :: undef           ! value to use when not extrapolating 
      
      real(kind=8), intent(out) :: A_prs(im,jm,kp) ! output in pressure coordinates

      ! This is a simple wrapper so that f2py interface list sane variable names;
      ! the order of arguments have also been adjusted.

      !                   ----
      

#if 0
      print *
      print *, 'Hello, World! I am inside the Fortran.'

      print *, 'im, jm, km, kp: ', im, jm, km, kp
      print *, 'xflag, vflag, method: ', xflag, vflag, method
      print *, 'A_eta:', shape(A_eta)

      !print *, 'A_eta:', shape(A_eta), minval(A_eta), maxval(A_eta)
      !print *, '   ps:', shape(ps),    minval(ps), maxval(ps)
      !print *, '   ts:', shape(ts),    minval(ts), maxval(ts)
      !print *, ' phis:', shape(phis),  minval(phis), maxval(phis)
      
      print *, 'plevs:', shape(plevs), minval(plevs), maxval(plevs)
      print *, '   ak:', shape(ak),    minval(ak), maxval(ak)
      print *, '   bk:', shape(bk),    minval(bk), maxval(bk)

      print *, '   ak:', ak(1:3), '...', ak(70:km)
      print *, '   bk:', bk(1:3), '...', bk(70:km)

#endif

      ! Legacy interpolator
      ! -------------------
      call VINTH2PECMWF ( A_eta, A_prs, ak, bk, plevs,   &
                          method, ps, undef, xflag,         &
                          im, jm, km, kp, vflag, ts, phis) 

      end subroutine eta2xprs

!...................................................................................................
!
!    --------------------------- Legacy code Follows ---------------------------------
!    Origin: https://github.com/yyr/ncl/blob/master/ni/src/lib/nfpfort/vinth2p_ecmwf.f
!    Converted to f90 with to_f90, then hand adjusted by A. da Silva:
!    - redefined hybrid coordinates ak/bk, eliminating p0 (like in GEOS)
!    - eliminated input parameter nlevip1 = nlevi+1
!    - eliminated worked space plevi which is now allocated internally
!    - fixed intents, eliminated DIMENSIONS, move declaration to top.


SUBROUTINE VINTH2PECMWF ( dati, dato, hbcofa, hbcofb, plevo, &
                          intyp, psfc, spvl, kxtrp,          &
                          imax, nlat, nlevi, nlevo, varflg, tbot, phis)

  IMPLICIT NONE

! Input parameters
! ----------------
real(kind=8), INTENT(IN) :: dati(imax,nlat,nlevi) ! Input array on eta coords
  
                                                  ! p = ak + bk * psfc
real(kind=8), INTENT(IN) :: hbcofa(nlevi)         !    ak [Pa]
real(kind=8), INTENT(IN) :: hbcofb(nlevi)         !    bk, non-dimentional
                                                  ! NOTE: if input is mid-layer then
                                                  ! HBCOF must be midlayer as well.

real(kind=8), INTENT(IN) :: plevo(nlevo)          ! pressure levels [hPa]
integer,      INTENT(IN) :: intyp                 ! Interpolation: 1=linear,2=log,3=loglog
real(kind=8), INTENT(IN) :: psfc(imax,nlat)       ! surface pressure [Pa]
real(kind=8), INTENT(IN) :: spvl                  ! value for UNDEFs
integer,      INTENT(IN) :: kxtrp                 ! 0: don't extrapolate, 1: extrapolate
integer,      INTENT(IN) :: imax                  ! number of longitudes
integer,      INTENT(IN) :: nlat                  ! number of latitudes
integer,      INTENT(IN) :: nlevi                 ! number of hybrid (eta) levels
integer,      INTENT(IN) :: nlevo                 ! number of output pressure levels
integer,      INTENT(IN) :: varflg                ! Input variable: -1: H, 1: T, 0: other
real(kind=8), INTENT(IN) :: tbot(imax,nlat)       ! Temperature closest to ground [K]
real(kind=8), INTENT(IN) :: phis(imax,nlat)       ! Surface geopotential [m2 s-2]

! Output parameters
! -----------------
real(kind=8), INTENT(OUT):: dato(imax,nlat,nlevo)

          
!         THIS ROUTINE INTERPLOATES CCM2/3 HYBRID COORDINATE DATA
!         TO PRESSURE COORDINATES USING PRESSURE SURFACES AS THE
!         COORDINATE SURFACE WHERE THE INTERPOLATION IS DONE.  THE
!         TYPE OF INTERPOLATION IS CURRENTLY A VARIANT OF TRANSFORMED
!         PRESSURE COORDINATES WITH THE  INTERPOLATION TYPE
!         SPECIFIED BY INTYP.  ALL HYBRID COORDINATE VALUES ARE
!         TRANSFORMED TO PRESSURE VALUES. WHERE THE
!         FORMULA FOR THE PRESSURE OF A HYBRID SURFACE IS;
!              P(K) = HBCOFA(K) + HBCOFB(K)*PSFC
!         WHERE,
!              HBCOFA - IS THE "A" OR PRESSURE HYBRID COEF [Pa]
!              K      - THE LEVEL INDEX (RUNNING FROM TOP TO BOTTOM)
!              HBCOFB - IS THE "B" OR SIGMA COEFICIENT
!              P(K)   - IS THE PRESSURE OF A HYBRID SURFACE IN MB.
!              PSFC   - IS THE SURFACE PRESSURE IN PASCALS
!                       (will later be converted to MB = .01*PASCALS)

!         FOR HYBRID DATA AT LEVEL INTERFACES SINCE THERE IS ONE
!         MORE VERTICAL LEVEL FOR INTERFACES THAN FOR LEVEL MIDPOINTS
!         IT IS ASSUNMED THAT THE FIRST INTERFACE LEVEL WITH A DATA
!         VALUE IS THE SECOND LEVEL FROM THE TOP.

!         ON INPUT-
!            DATI    - 3 DIMENSIONAL ARRAY (I,J,KI) CONTAINING DATA
!                      ON HYBRID SURFACES  WHERE I IS LONGTIUDE, J
!                      IS LATITUDE AND K IS THE VERTICAL HYBRID
!                      COORDINATE.  THE VERTICAL DATA RUN TOP TO BOTTOM.
!                      SIGMA DATA WITH THE DATA ORDERED TOP TO BOTTOM.
!            HBCOFA  - ARRAY CONTAINING "A" OR PRESSURE
!                      COEFICIENTS FOR COMPUTING PRESSURE AT A LEVEL.
!                      ARRAY IS 2XNLEVIP1.  THE 1ST INDEX TAKES ON
!                      THE VALUE OF EITHER
!                      NOTE THAT INTERNALLY COEFICIENTS ARE SCALED TO YIELD A
!                      PRESSURE IN MB.  THEY ARE ORDERED FROM TOP
!                      OF THE MODEL TO THE BOTTOM.
!            HBCOFB  - SAME AS HCOFA BUT FOR THE "B" OR SIGMA COEFICIENT
!            PLEVO   - LIST OF OUTPUT PRESSURE SURFACES IN MB
!                      LOW TO HIGH PRESSURE
!            INTYP   - A FLAG INDICATING INTERPOLATION FOR EACH
!                      FIELD (1 - LINEAR,2 - LOG ,3 - LOG LOG)
!                      WHERE EACH INTERPOLATION IS DONE IN TRANSFORMED
!                      PRESSURE COORDINATES.
!            PSFC    - MODEL SFC PRESSURE IN PASCALS (WILL BE CONVERTED
!                      TO MB)
!            VCOLI   - ARRAY TO STORE A LONGITUDINAL VERTICAL SLICE OF
!                      INPUT DATA (IMAX BY NLEVI).
!            VCOLO   - SAME BUT FOR OUTPUT DATA (IMAX BY NLEVO)
!            IMAX    - LONGITUDINAL DIMENSION OF THE DATA.
!            NLAT    - LATITUDINAL DIMENSION OF THE DATA.
!            NLEVI   - NO. OF LEVELS FOR THE HYBRID DATA
!            NLEVIP1 - NLEVI + 1
!            NLEVO   - NUMBER OF OUTPUT LEVELS FOR PRESSURE DATA
!            KXTRP   - FLAG WHICH INDICATES WHETHER OR NOT
!                      EXTRAPOLATION WILL BE USED WHEN THE OUTPUT
!                      PRESSURE SURFACE IS BELOW THE LOWEST LEVEL
!                      OF THE MODEL.
!                         0 - DON'T EXTRAPOLATE USE SPECIAL VALUE SPVL
!                         1 - EXTRAPOLATE DATA using ECMWF formulation
!                             below PSFC
!            SPVL    - SPECIAL VALUE TO USE WHEN DATA IS NOT
!                      EXTRAPOLATED
!            varflg  - flag which indicates the name of the variable
!                      -1 means geopotential (Z)
!                      +1 means geopotential (T)
!                       0 any other variable
!            tbot    - temperature at level closest to ground
!            phis    - surface geopotential

!         ON OUTPUT-
!            DATO  - 3 DIMENSIONAL ARRAY TO HOLD DATA INTERPOLATED
!                    TO PRESSURE SURFACES.


!               ------

! Local space
! -----------
! PLEVI -  1 DIMENSIONAL ARRAY TO HOLD PRESSURE VALUES
!          OF HYBRID SURFACES FOR A VERTICAL COLUMN
!          SLICE
!integer :: nlevip1 = nlevi+1 
real(kind=8) :: plevi(nlevi+1)


REAL(kind=8) :: tstar,hgt,alnp,t0,tplat,tprime0,alpha  
REAL(kind=8) :: a2ln,a1,alph,psfcmb
INTEGER      :: i,j,k,kp,kpi,iprint

! NCLEND

! for ecmwf extrapolation
REAL(kind=8) :: rd,ginv
PARAMETER (rd=287.04D0)
PARAMETER (ginv=1.d0/9.80616D0)
PARAMETER (alpha=0.0065D0*rd*ginv)

!         STATEMENT FCN. FOR DOUBLE LOG. INTERP ON PRESSURE SURFACES
!         PRESUMES PRESSURE IS IN MB

a2ln(a1) = LOG(LOG(a1+2.72D0))

!         STATEMENT FCN. FOR DOUBLE LOG. INTERP ON SIGMA SURFACES.
!         SETS UP ROUGH UPPER BPOUND SIMILAR TO STATEMENT FCN FOR
!         PRESSURE. I.E.    FIXED VALUE LN(LN(P) = LN(LN(FIXED VAL)
!         AT .001 SIGMA OR ABOUT 1 MB

!     A2LN(A1)=LOG(LOG(A1+1.001))

!    

DO  j = 1,nlat
   DO  i = 1,imax
      
! =======================================DJS special case===
    IF (psfc(i,j) == spvl) THEN
      DO k = 1,nlevo
        dato(i,j,k) = spvl
      END DO
      CYCLE
    END IF
! =========================================================
    
    
!         GET PRESSURE VALUES FOR HYBRID SURFACES FOR THIS POINT
!         AND FOR THE TYPE OF MODEL SURFACE THE DATA IS ON.
!         INTERFACE DATA STARTS AT THE SECOND INTERFACE LEVEL SO
!         IF THE DATA IS ON THOSE LEVELS START THE
    
    DO k = 1,nlevi
      kpi = k
!ams  plevi(k) = (hbcofa(kpi)*p0) + hbcofb(kpi)* (psfc(i,j)*.01D0)
      plevi(k) = (hbcofa(kpi) + hbcofb(kpi)* psfc(i,j) ) * 0.01 ! hPa
    END DO
    
!         CALL P2HBD TO PERFORM VERTICAL INTERP. THEN TRANSFER DATA TO
!         THE OUTPUT ARRAY
    
    DO  k = 1,nlevo
      
!         CHECK FOR BRACKETING LEVEL KP WILL BE THE INPUT LEVEL THAT
!         IS THE UPPER PORTION OF 2 INPUT BRACKETING LEVELS.
      
!         IF BRANCH FOR MODEL TOP
      
      IF (plevo(k) <= plevi(1)) THEN
        kp = 1
        GO TO 30
        
!         IF BRANCH FOR LEVEL BELOW LOWEST HYBRID LEVEL
        
      ELSE IF (plevo(k) > plevi(nlevi)) THEN
        IF (kxtrp == 0) THEN
          dato(i,j,k) = spvl
          GO TO 40
        ELSE IF (varflg > 0) THEN
! Variable is "T" and ECMWF extrapolation is desired
          psfcmb = psfc(i,j)*0.01D0
          tstar = dati(i,j,nlevi)* (1.d0+alpha* (psfcmb/plevi(nlevi)-1))
          hgt = phis(i,j)*ginv
          IF (hgt < 2000.d0) THEN
            alnp = alpha*LOG(plevo(k)/psfcmb)
          ELSE
            t0 = tstar + 0.0065D0*hgt
            tplat = MIN(t0,298.d0)
            IF (hgt <= 2500.d0) THEN
              tprime0 = 0.002D0* ((2500.d0-hgt)*t0+ (hgt-  &
                  2000.d0)*tplat)
            ELSE
              tprime0 = tplat
            END IF
            IF (tprime0 < tstar) THEN
              alnp = 0.d0
            ELSE
              alnp = rd* (tprime0-tstar)/phis(i,j)* LOG(plevo(k)/psfcmb)
            END IF
          END IF
          dato(i,j,k) = tstar* (1.d0+alnp+.5D0*alnp**2+ 1.d0/6.d0*alnp**3)
          GO TO 40
          
        ELSE IF (varflg < 0) THEN
! Variable is "Z" and ECMWF extrapolation is desired
          psfcmb = psfc(i,j)*0.01D0
          hgt = phis(i,j)*ginv
          tstar = tbot(i,j)* (1.d0+ alpha* (psfcmb/plevi(nlevi)-1.d0))
          t0 = tstar + 0.0065D0*hgt
          
          IF (tstar <= 290.5D0 .AND. t0 > 290.5D0) THEN
            alph = rd/phis(i,j)* (290.5D0-tstar)
          ELSE IF (tstar > 290.5D0 .AND.  &
                t0 > 290.5D0) THEN
            alph = 0
            tstar = 0.5D0* (290.5D0+tstar)
          ELSE
            alph = alpha
          END IF
          
          IF (tstar < 255.d0) THEN
            tstar = 0.5D0* (tstar+255.d0)
          END IF
          alnp = alph*LOG(plevo(k)/psfcmb)
          dato(i,j,k) = hgt - rd*tstar*ginv* LOG(plevo(k)/psfcmb)*  &
              (1.d0+.5D0*alnp+ 1.d0/6.d0*alnp**2)
          GO TO 40
        ELSE
! Use lowest sigma layer
          dato(i,j,k) = dati(i,j,nlevi)
          GO TO 40
        END IF
        
!         IF BRANCH FOR TO CHECK IF OUTPUT LEVEL IN BETWEEN
!         2 LOWEST HYBRID LEVELS
        
      ELSE IF (plevo(k) >= plevi(nlevi-1)) THEN
        kp = nlevi - 1
        GO TO 30
        
!         IF BRANCH FOR MODEL INTERIOR
!         LOOP THROUGH INPUT LEVELS TILL YOU ARE BRACKETING
!         OUTPUT LEVEL
        
      ELSE
        kp = 0
        20                 CONTINUE
        kp = kp + 1
        IF (plevo(k) <= plevi(kp+1)) GO TO 30
        IF (kp > nlevi) THEN
          WRITE (6,FMT=25) kp,nlevi
          25                     FORMAT (' KP.GT.KLEVI IN P2HBD.  KP,KLEVI= ',  &
              2I5)
!            CALL DPRNT(' PLEVI',PLEVI,NLEVIP1,1,1,1)
!            CALL DPRNT(' PLEVO',PLEVO,NLEVO,1,1,1)
!            CALL ABORT(' KP.GT.NLEVI IN P2HBD')
        END IF
        GO TO 20
      END IF
      30             CONTINUE
      
!         LEVEL BRACKETED PICK TYPE OF INTERP.
      
      
!         LINEAR INTERP.
      
      IF (intyp == 1) THEN
        dato(i,j,k) = dati(i,j,kp) + (dati(i,j,kp+1)-dati(i,j,kp))*  &
            (plevo(k)-plevi(kp))/ (plevi(kp+1)-plevi(kp))
        
!         LOG INTERPOLATION.
        
      ELSE IF (intyp == 2) THEN
        iprint = 1
!      IF (I.EQ.1.AND.IPRINT.EQ.1) THEN
!         PRINT 101,I,J,K,KP,ILEV
!  101    FORMAT('  IN S2HBD I,J,K,KP,ILEV ',5I3)
!         PRINT 102,DATI(I,J,KP),DATI(I,J,KP+1),PLEVO(K),
!     *             PLEVI(KP),PLEVI(KP+1)
!  102    FORMAT(' DATI(KP),DATI(KP+1),PLEVO(K),',
!     *          'PLEVI(KP),PLEVI(KP+1) ',
!     *          /,1X,1P5E12.5)
!      ENDIF
        dato(i,j,k) = dati(i,j,kp) + (dati(i,j,kp+1)-dati(i,j,kp))*  &
            LOG(plevo(k)/plevi(kp))/ LOG(plevi(kp+1)/plevi(kp))
        
!        FOR LOG LOG INTERP. NOTE A2LN IS A STATEMENT FCN.
        
      ELSE IF (intyp == 3) THEN
        dato(i,j,k) = dati(i,j,kp) + (dati(i,j,kp+1)-dati(i,j,kp))*  &
            (a2ln(plevo(k))-a2ln(plevi(kp)))/ (a2ln(plevi(kp+1))-a2ln(plevi(kp)))
      END IF
      40             CONTINUE
      
    END DO
  END DO
END DO
RETURN

END SUBROUTINE vinth2pecmwf

#if 0

SUBROUTINE dprnt(ifld,a,im,jm,is,js)

CHARACTER(LEN=*), intent(in) :: ifld
REAL(kind=8), INTENT(IN OUT) :: a
INTEGER, INTENT(IN)          :: im
INTEGER, INTENT(IN)          :: jm
INTEGER, INTENT(IN OUT)      :: is
INTEGER, INTENT(IN)          :: js

DIMENSION a(im,jm)

PRINT 10,ifld
10 FORMAT (1X,/,' FIELD ',a)
DO  j = 1,jm,js
  PRINT 13,j
  13     FORMAT (' J=',i4)
  PRINT 15, (a(i,j),i=1,im,is)
  15     FORMAT (1X,1P,10D12.5)
END DO
RETURN
END SUBROUTINE dprnt

#endif

