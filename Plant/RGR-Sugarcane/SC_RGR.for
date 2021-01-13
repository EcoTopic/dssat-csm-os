c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     SASRI CANEGRO
c     &------------
c     Matthew Jones, September - December 2006
c     Gainesville, Florida, USA
c     &----------------------------------------------------
c     The SASRI Canegro was modularised according to the 
c     DSSAT structure and incorporated into the DSSAT
c     CSM.  This serves both as an update to the legacy
c     DSSAT sugarcane model (also a CANEGRO, from years
c     ago), as well as a restructured CANEGRO, a platform
c     upon which further development will hopefully be
c     easier and more effective.
c     ::::::::::::::::::::::::::::::::::::::::
c     This work was funded by the International Consortium 
c     for Sugarcane Modelling (ICSM):
c     https://sasri.sasa.org.za/agronomy/icsm/
c     ::::::::::::::::::::::::::::::::::::::::
c     Acknowledgements and thanks to:
c     Dr. Abraham Singels (SASRI), Dr. Maurits van den Berg 
c     (SASRI), Dr. Jim Jones (UF), Cheryl Porter (UF), 
c     Jim Shine (SGC / UF), Gerrit Hoogenboom (UGA) 
c     ... and Geoff Inman-Bamber (CSIRO?), for creating the 
c     Sugarcane infrastructure in DSSAT
c     -----------------------------------------------------
c     Changes for v4.7
c     Gainesville, Florida, January 2018.
c     Matthew Jones
c
c     1. Temperature responses (cardinal temperatures for 
c        many processes).
c     2. Shoot population (linear tillering rate per 
c        primary shoot with light interception-based
c        slowing).
c     3. Respiration (growth and maintenance)
c     4. CO2 impacts (lookup function now has better support
c        for low CO2, no direct benefit to high CO2.)
c     5. [optional: aquacrop water uptake] 
c     &----------------------------------------------------



c     SUBROUTINE: SC_CNGRO()
c     &---------------------
c     This sub interfaces with the rest of the DSSAT CSM.
c     Inputs and outputs are passed solely through this
c     interface.  This forms the PLANT GROWTH MODULE for
c     sugarcane.  Soil and atmospheric processes, along 
c     with external calculations like irrigation are
c     handled OUTSIDE this module.
c     [interface copied from MZ_CERES.for]
!  Added IRRAMT July 2015
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SC_RGR (
     &    CONTROL, ISWITCH,                                   !Input
     &    CO2, DAYL, EOP, EP, EO, ES, HARVFRAC, NH4, NO3, SNOW,   !Input
     &    SOILPROP, SRAD, SW, TMAX, TMIN, TRWUP, TRWU, EOS,   !Input
     &    RWUEP1, TWILEN, YREND, YRPLT, WEATHER, IRRAMT,      !Input
     &    CANHT, HARVRES, KCAN, KTRANS, MDATE, NSTRES,        !Output
     &    PORMIN, RLV, RWUMX,SENESCE, STGDOY, UNH4,           !Output
     &    UNO3, XLAI, XHLAI, EORATIO)                 !Output
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     Define DSSAT composite variables:
c     [Taken from MZ_CERES.for]
      USE ModuleDefs

c     :::::::::::::::::::::::::::::::::::::::::::::::::::::

      IMPLICIT NONE
      SAVE
      
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     DSSAT variable declarations (i.e. subroutine 
c     argument declarations)
c     Explanation of variables in Appendix 1 at the 
c     end of this subroutine.
c     :::::::::::::::::::::::
c     & Input variables -
c     :::::::::::::::::::
      REAL     CO2
      REAL     DAYL
      REAL     EOP
      REAL     EP, EP1, RWUEP1, RWUEP2
      REAL     EO, EOS, ES
      REAL     HARVFRAC(2)
      REAL     NH4(NL)
      REAL     NO3(NL)
      REAL     SNOW
      REAL     SRAD
      REAL     SW(NL)
      REAL     TMAX
      REAL     TMIN
      REAL     TRWUP
      REAL     TRWU
      REAL     TWILEN
      INTEGER  YREND
      INTEGER  YRPLT
c     Daily irrigation amount (mm)
      REAL IRRAMT
c     & Output variables: -
c     :::::::::::::::::::::
      REAL     CANHT
      REAL     KCAN
      REAL     KTRANS
      INTEGER  MDATE
      REAL     NSTRES
      REAL     PORMIN
      REAL     RLV(NL)
! c     Maximum rate of root water uptake per unit root length, RWUMX (today) and RWUMX_REF (limited by water stress)
      REAL     RWUMX ! , RWUMX_REF
      INTEGER  STGDOY(20)
      REAL     UNH4(NL)
      REAL     UNO3(NL)
      REAL     XLAI
      REAL     XHLAI 
      REAL     EORATIO
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     DSSAT composite variables:
c     [Taken from MZ_CERES.for]
c     ::::::::::::::::::::::::::
      TYPE (ControlType) CONTROL
      TYPE (SoilType)    SOILPROP
      TYPE (SwitchType)  ISWITCH
      Type (ResidueType) HARVRES 
      Type (ResidueType) SENESCE
      Type (WeatherType) WEATHER


c     & Local variables -
c     :::::::::::::::::::
c     local copy of CONTROL%DYNAMIC:
      INTEGER DYNAMIC
c     Output file unit:
c      INTEGER NOUTDG
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     CANEGRO variable declarations (variables internal
c     to this subroutine and sub-subroutines)
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     Heat units for today: base10, base16, cultivar-
c     defined heat units for leaf dynamics & emergence
c      & population
      REAL HU10, HU16, HUBaseLeaf, HUBaseEm, HUBasePop
c     Cumulative versions of above:
      REAL CHU10, CHU16, CHUBaseLeaf, CHUBaseEm, CHUBasePop

c     Thermal time measures for canesim canopy:
      REAL  HU_CANESIM, CHU_CANESIM

c     Average temperature for today:
c      REAL AVGTEMP

c     Base temperatures for leaf and emergence
      REAL TBaseLeaf, TBaseEm, TBasePop

c     Yesterday's heat cum. thermal time, base16
      REAL T0

c     Canopy/population variables
      REAL    pplast
      INTEGER newbat
      REAL    rowspc, stdaye,stddiv, RATNUM

c     Local variables required for Photosynthesis
c      REAL SCV
      REAL CWSI
      REAL RECTIM 
      REAL SWDFMN
      REAL CRITSW
      REAL SLPF

c     Thermal time required to for photosynthesis to recover
c     from water stress (cultivar param):
      REAL HuRecover

      REAL ALAT
      INTEGER NEGGRO
      REAL DWDT
      REAL CODERD
      REAL RUNOFF
      INTEGER ALLOWLODG

c     Lodging variables (note: look in canop3 for lodge coeff.)
      REAL LODGE_LI

      INTEGER IWS
      REAL STDAYC, SENESF, EXPLA, PER, RESET, BOOTLI
c     BEST is a cultivar input
      REAL STRPOP, BEST
      INTEGER MDL

c      Average RM FRACtion and corres. mean temperature so far this sim run
       REAL RM_AVG, TMEAN_AVG

c     Cellulosic DM (t/ha)
      REAL  CELLSE_DM

c     Thermal time requirements for phenology
c     Germination:
!      REAL TT_GERM
c     Emergence
      REAL TT_EMERG
      LOGICAL EMERGED, FILE_EXISTS

c     General local vars:
      INTEGER I

c     Choice of Canopy routine
      INTEGER CANOPY_ROUTINE !, CANESIM, CANEGRO

c     

c     Error variable for cultivar coefficient read
!      LOGICAL CF_ERR
c     Error variable for species coefficient read
      LOGICAL SPC_ERROR
!      REAL HI

c     temp output file (remove)
      CHARACTER TFILENAME*25
      
c     Row-spacing warning text:
      CHARACTER*78 WARNINGS(2)      

c     REMOVE??
!      INTEGER INIDOY 
!      REAL INSWDF1, INSWDF2
!      INTEGER TDOY
!      REAL    TempLAI

!     CHP ADDED 2/13/2007
      REAL MAXLAI

c     Should the ratoon info be carried over?
c      LOGICAL CARRY - replaced by CaneCrop%CARRY_OVER

c     Partitioning:
c     To calc fresh mass
      REAL STKDMC, STKWAT

c     Temp:
!      REAL ES, ES_mm, ASWC_mm
      REAL cCOEFF
      LOGICAL CERROR

!     Unit number for output
      INTEGER SCLUN   !CHP


c     Yesterday's value of delta sucrose mass:
c     t/ha/d
      REAL DSUCyest

c     Soil water stress for tillering (not used in this sub):
      REAL SWDF30
      
c     Sunlit LAI fraction
      REAL F_SL      

c     Variables for calculating waterlogging stress
c     (saturation factor)      

      REAL  SATFAC
      REAL SUMEX, SUMRL, SWEXF
      INTEGER L
      REAL TSS(NL)
      

c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ~~~~~~~~ SUBROUTINE CODE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     Initialisation for DYNAMIC = all
c     ::::::::::::::::::::::::::::::::
c     Init DYNAMIC:
      DYNAMIC = CONTROL%DYNAMIC

c     Copy Soil properties from DSSAT object to CANEGRO
c     object:
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ! DO I=1,NL
      !     Soil%DUL(I)   = SOILPROP%DUL(I)
      !     Soil%LL(I)    = SOILPROP%LL(I)
      !     Soil%SW(I)    = SW(I)
      !     Soil%DLAYR(I) = SOILPROP%DLAYR(I)
      !     Soil%NLAYR    = SOILPROP%NLAYR
      !     Soil%SAT      = SOILPROP%SAT
      ! ENDDO
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c    
c              DYNAMIC = RUNINIT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      IF (DYNAMIC.EQ.RUNINIT) THEN
            WRITE(*,*) 'SC_RGR CALLED IN RUNINIT MODE...'
     

c     Call the output routine, to initialise output
c     ::::::::::


c     create work.out
      IF (INDEX('YDA',ISWITCH % IDETL) > 0) THEN
        TFILENAME = 'Work.OUT'
        CALL GETLUN('WORK.OUT', SCLUN)
        INQUIRE(FILE=TFILENAME, EXIST=FILE_EXISTS)

c       Open the file
        IF (FILE_EXISTS) THEN
c         In append mode if the file already exists
          OPEN(UNIT=SCLUN,FILE=TFILENAME,STATUS='OLD',POSITION='APPEND')
        ELSE
c         A new file if not existing
          OPEN (UNIT=SCLUN, FILE=TFILENAME, STATUS='NEW')
          WRITE(SCLUN,'("*CaneGro Supplemental Output File")')
        ENDIF

c       Output a header (treatment / run info, etc)
c       ::::::::::::::::::::::::
c        Use the predefined one:
        CALL HEADER(SEASINIT, SCLUN, CONTROL%RUN)
      ENDIF

c     I don't know what this ought to be... 1 seems fine.
      KCAN = 1.




c     Call the cultivar routine at least once:
      CALL GET_CULTIVAR_COEFF(cCOEFF, 'dummy', CONTROL, cERROR)


c     END of RUNINIT
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = SEASINIT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF (DYNAMIC.EQ.SEASINIT) THEN
            WRITE(*,*) 'SC_RGR CALLED IN SEASINIT MODE...'

c     EORATIO: this seems to be the closest thing to a crop coefficient
c     that DSSAT offers for the PENMON model at least.
c     ::::::::::::::::::::::::::::::::::::::::::::::::
      EORATIO = 1.15
      CALL GET_SPECIES_COEFF(EORATIO,'EORATIO', CONTROL, SPC_ERROR) 
c     ::::::::::::::::::::::::::::::::::::::::::::::::



c     Until N is explcitly modeled, N stress will be 1 (no stress)
      NSTRES = 1.

c     Canopy height
      CANHT = 0.

c     Leaf area index:
      XLAI   = 0.
      XHLAI  = 0.



c     Coefficients for SWDF1, SWDF2
c     From the species file
c     :::::::::::::::::
          RWUEP1         = 1.
          RWUEP2         = 2.

          CALL GET_SPECIES_COEFF(RWUEP1,'RWUEP1', CONTROL, SPC_ERROR)
          CALL GET_SPECIES_COEFF(RWUEP2,'RWUEP2', CONTROL, SPC_ERROR)


c         Critical soil water content?
c         ::::::::::::::::::::::::::::
c         Read from species file:
          CRITSW = 0.2
          CALL GET_SPECIES_COEFF(CRITSW,'CRITSW', CONTROL, SPC_ERROR)

c     Initialise population variables:
c     :::::::::::::::::::::::::::::::::::::::
c         READ from fileX!!!      
c         Set the standard rowspacing
          rowspc = 1.4
          CALL FIND_INP(ROWSPC, 'ROWSPC', Control)
          
c         MJ, March 2010: check rowspacing
c         ::::::::::::::::::::::::::::::::
c         Issue a warning if row-spacing is less than 3 cm.  The reason
c         being that someone might inadvertantly enter 1.5, thinking the input
c         is metres.  Rowspacing is unlikely to exceed 3 m, hence this threshold.
c         Note: rowspc is converted as a special case from cm to m by the FIND_INP
c         routine.
          IF (ROWSPC .LT. .03) THEN
            WARNINGS(1) = "Warning: please check that " 
     &      // "row-spacing is specified in CENTIMETRES"
            WARNINGS(2) = "[Row-spacing appears to be very low]"
            CALL WARNING(2, "SCCAN"//ModelVerTxt//": SCCNGRO", WARNINGS)
          ENDIF          



c     Is this a plant/ratoon crop?
          CALL FIND_INP(RATNUM, 'RATOON', Control)


c	Water logging stress 
      PORMIN = 0.05 
      CALL GET_SPECIES_COEFF(PORMIN,'PORM', CONTROL, SPC_ERROR) 
      ! init waterlogging stress
      !WRITE(*, '(A, F10.5)') 'PORMIN is ', PORMIN
      DO I=1, NL
        TSS(I) = 0.0
      ENDDO
      SATFAC = 0.0

c     END of SEASINIT
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     &----------------------------------------------------
c
c              DYNAMIC = RATE
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.RATE) THEN
            WRITE(*,*) 'SC_RGR CALLED IN RATE MODE...'


c     Calculate water stresses:
c     :::::::::::::::::::::::::::
c     Only call this sub if the crop has emerged:
c     :::::::::::::::::::::::::::::::::::::::::::
!      IF (EMERGED) THEN
!        IF (ISWITCH%ISWWAT .EQ. 'N') THEN
c         Water balance is NOT used
!          WaterBal%SWDF1 = 1.0
!          WaterBal%SWDF2 = 1.0
!          SATFAC = 0.0
!        ELSE
!          WaterBal%EOP = EOP
!          WaterBal%EOS = EOS
!          WaterBal%TRWUP = TRWUP

c         If (1) the water balance is simulated
c            (2) Potential transpiration > 0.
c            (3) LAI > 1.E-4, calculate stresses:
c         ::::::::::::::::::::::::::::::::::::::::::::::::::::::
!          IF ((EOP .LT. 0.00001) .OR. 
!     &      .NOT.(Growth%LAI .GT. 0.0001)) THEN
c           If there is NO potential transp., there can be no water stress...
!            WaterBal%SWDF1 = 1.0
!            WaterBal%SWDF2 = 1.0
!            WaterBal%ANAERF = 1.0
!            SATFAC = 0.0
!          ELSE
c           Init:
!            WaterBal%SWDF1 = 1.0
!            WaterBal%SWDF2 = 1.0
!            EP1 = EOP * 0.1

c           SWDF1 
!            WaterBal%SWDF1 = MIN(1., (1./RWUEP1* TRWUP/EP1))
c           SWDF2
!            WaterBal%SWDF2 = MIN(1., (1./RWUEP2 * TRWUP/EP1))
            
c           Find an equivalent value from DSSAT!
            ! WaterBal%ANAERF = 1.
            ! MJ, Jan 2018: Calculate anaerobic stress factor
            PORMIN = 0.20
!            CALL SC_ANAERF(
!     &        Soil%NLAYR, RLV, SW, Soil%DUL, Soil%SAT, 
!     &        Soil%DLAYR, CONTROL, WaterBal%ANAERF)
!c           SCV    = 0.2
!            WRITE(*, '(F10.5)') WaterBal%ANAERF

          !-------------------------------------------------------------
          !      Compute Water Saturation Factors       
          ! ------------------------------------------------------------
          ! taken from the CERES-Maize model, MZ_GROSUB.for
          ! Added by MJ, Jan 2018
          ! This is meant as a sort of fudge.  I will adjust PORMIN 
          ! until I get a similar fit to v4.5 with the AquaCrop 
          ! waterlogging stress routine.
          

          SATFAC = 0.0    
          SUMEX = 0.0
          SUMRL = 0.0
      
          DO L = 1, SoilProp%NLAYR

          !------------------------------------------------------------
          !PORMIN = Minimum pore space required for supplying oxygen to 
          !         roots for optimum growth and function    
          !TSS(L) = Number of days soil layer L has been saturated 
          !         above PORMIN
          !------------------------------------------------------------
              IF ((SOILPROP%SAT(L)-SW(L)) .GE. PORMIN) THEN
                  TSS(L) = 0.
              ELSE
                  TSS(L) = TSS(L) + 1.
              ENDIF
          !------------------------------------------------------------
          ! Delay of 2 days after soil layer is saturated before root
          ! water uptake is affected
          !------------------------------------------------------------
              IF (TSS(L) .GT. 2.) THEN
                  SWEXF = (SOILPROP%SAT(L)-SW(L))/PORMIN
                  SWEXF = MAX(SWEXF,0.0)
              ELSE
                  SWEXF = 1.0
              ENDIF

              SWEXF = MIN(SWEXF,1.0)
              SUMEX  = SUMEX + SOILPROP%DLAYR(L)*RLV(L)*(1.0 - SWEXF)
              SUMRL  = SUMRL + SOILPROP%DLAYR(L)*RLV(L)
          ENDDO

          IF (SUMRL .GT. 0.0) THEN
              SATFAC = SUMEX/SUMRL
          ELSE
              SATFAC = 0.0
          ENDIF
          SATFAC = AMAX1(SATFAC,0.0)
          SATFAC = AMIN1(SATFAC,1.0)
          ! WaterBal%ANAERF is passed to SC_PHOTOS and reduces
          ! photosynthesis rates.
          ! The same calculation is done in the ROOTWU routine
          ! (SATFAC reduces water uptake), and SWDF1 (and SWDF2)
          ! are calculated a function of total root water uptake
          ! and potential transp. - so this is maybe
          ! double-counting?  (Comment in MZ_GROSUB confirms this)
          ! WaterBal%ANAERF = 1.0 - SATFAC
          ! In CERES Maize, SATFAC is used in the calculation of
          ! leaf area growth, but not in photosynthesis.  Leaf
          ! area in growth in Canegro is limited by SWDF2, which
          ! includes some effect of SATFAC in ROOTWU... but there
          ! might be an argument for an additional effect of
          ! SATFAC on leaf elongation rate?
                      
!          ENDIF
!        ENDIF
!      ENDIF



c     END of RATE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = INTEGRATE
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.INTEGR) THEN
            WRITE(*,*) 'SC_RGR CALLED IN INTEGRATE MODE...'


c     END of INTEGRATE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
c              DYNAMIC = OUTPUT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.OUTPUT) THEN

c     Write output of daily values:
c     :::::::::::::::::::::::::::::



!         chp 4/7/2009
          IF (CONTROL%YRDOY == YREND) THEN
            STGDOY(16) = YREND
            MDATE = YREND
          ENDIF


c     END of OUTPUT
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = FINAL 
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.SEASEND) THEN
            WRITE(*,*) 'SC_RGR CALLED IN FINAL MODE...'



c     Temporary!!
      CLOSE(SCLUN)



c     END of FINAL
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c     END of DYNAMIC conditional statement
c     ::::::::::::::::::::::::::::::::::::
      ENDIF
c     ::::::::::::::::::::::::::::::::::::


      END
c     -----------------------------------------------------
c     END of SUBROUTINE SC_CNGRO
c     -----------------------------------------------------
c     APPENDIX 1: DSSAT interface variables:
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     - Input variables -
c     :::::::::::::::::::
c      CONTROL:  �init�, �rate�, etc 
c      ISWITCH:  Indicates which processes (Soil, N, etc) are to be modelled 
c      CO2:      Atmospheric CO2, ppm 
c      DAYL:     Day length on day of simulation (from sunrise to sunset) (hr) 
c      EOP:      Potential plant transpiration, mm/day 
c      EP:       Actual plant transpiration mm/day
c      HARVFRAC: Two-element array containing fractions of (1) yield harvested 
c                and (2) by-product harvested (fraction) 
c      NH4:      Ammonium in soil layer L, mg elemental N/kg soil 
c      NO3:      Nitrate in soil layer L (mg elemental N/kg soil) 
c      SNOW:     Snow accumulation today, mm 
c      SOILPROP: Soil properties (DUL, LL, etc) 
c      SRAD:     Total solar radiation today, MJ/m2   
c      Sw:       Volumetric soil water content of soil layer, cm3 water/cm3 soil    
c      TMAX:     Max temp., degrees C 
c      TMIN:     Min temp., degrees C 
c      TRWUP:    Total potential root water uptake, cm/day 
c      TWILEN:   TWILEN is Twilight to Twilight daylength 
c      YREND:    Year and day of year for harvest 
c      YRPLT:    Year and day of year of planting 
c     :::::::::::::::::::
c     - Output variables -
c     :::::::::::::::::::
c      CANHT:    Canopy height (m)
c      HARVRES:  Composite variable containing harvest residue amounts for total 
c                dry matter, lignin, and N amounts (kg/ha)
c      KCAN:     Canopy light extinction coefficient for daily PAR, for equidistant 
c                plant spacing, modified when in-row and between  row spacing are not equal
c      KTRANS:   Energy extinction coefficient for partitioning EO (potential plant 
c                transpiration) to EP (actual plant transpiration) [previously KEP]
c      MDATE:    Year and day of year of maturity
c      NSTRES:   Nitrogen stress factor (1=no stress, 0=max stress)
c      PORMIN:   Minimum pore space required for supplying oxygen to roots for optimal 
c                growth and function (cm3/cm3)
c      RLV:      Root length density for soil layer L (cm[root] / cm3[soil])
c      RWUMX:    Maximum root water uptake parameter from species file
c      SENESCE:  Composite variable containing data about daily senesced plant matter. 
c      STGDOY:   Year and day of year that a growth stage occurred on
c      UNH4:     Plant uptake of ammonium from layer (kg N/ha/day)
c      UNO3:     Plant uptake of nitrate from a layer (kg N/ha/day)
c      XLAI:     Leaf area index, m2/m2
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
