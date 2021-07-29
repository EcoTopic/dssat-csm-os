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

c     JvdM: RGR genetic input parameters
c     -----------------
      REAL InitialGLAI  ! Initial green leaf area index (m2/m2)
      REAL TBase_LAI  ! Cardinal temperatures for canopy development  (°C)
      REAL TOpt_LAI  ! Cardinal temperatures for canopy development  (°C)
      REAL TFin_LAI  ! Cardinal temperatures for canopy development  (°C)
      REAL TBase_Photos  ! Cardinal temperatures for photosynthesis  (°C)
      REAL TOpt1_Photos  ! Cardinal temperatures for photosynthesis  (°C)
      REAL TOpt2_Photos  ! Cardinal temperatures for photosynthesis  (°C)
      REAL TFin_Photos  ! Cardinal temperatures for photosynthesis  (°C)
      REAL TBase_LFAPP  ! Cardinal temperatures for leaf appearance (°C)
      REAL TOpt_LFAPP  ! Cardinal temperatures for leaf appearance (°C)
      REAL TFin_LFAPP  ! Cardinal temperatures for leaf appearance (°C)
      REAL LeafPI1  ! Leaf phyllocron intervals (°Cd):
      REAL LeafPI2  ! Leaf phyllocron intervals (°Cd):
      REAL KeMin  ! Minimum and maximum PAR extinction coefficients
      REAL KeMax  ! Minimum and maximum PAR extinction coefficients
      REAL KeMaxLf  ! Number of leaves/stalk at KeMax
      REAL SLAMin  ! Minimum and maximum specific leaf area (cm2/g)
      REAL SLAMax  ! Minimum and maximum specific leaf area (cm2/g)
      REAL Suc_LfNum_Delay  ! delay in number of leaves after onset of stalk growth when sucrose is permitted to start accumulating
      REAL RGRglaiMin  ! Minimum and Maximum relative growth rate of green leaf area index (RGRglai), and the slope
      REAL RGRglaiMax  ! Minimum and Maximum relative growth rate of green leaf area index (RGRglai), and the slope
      REAL RGRglaiSlope  ! coefficient describing the transition from max to min as the crop transitions to stalk growth
      REAL RUEo  ! maximum theoretical radiation use efficiency, under optimal conditions (ADM PAR basis, g/MJ)
      REAL MAX_ROOTPF  ! Root partitioning: average partitioning of mass to roots, maximum partitioning fraction
      REAL AvRootDMFrac  ! Root partitioning: average partitioning of mass to roots, maximum partitioning fraction
      REAL APFMX  ! Root partitioning: average partitioning of mass to roots, maximum partitioning fraction
      REAL PCB  ! Root partitioning: average partitioning of mass to roots, maximum partitioning fraction
      REAL RUE_FT_c ! MARCH 2021
      REAL FI_OSG  ! Fractional PAR interception at which 50% of shoots have started stalk growth
      REAL FIinOSG ! JULY 2021 FI_OSG independent of ROWSPC
      REAL OSG_log_c1  ! Parameter controlling the rate at which the crop transitions to stalk growth in response to PAR interception
      REAL STKPFmax  ! Maximum fraction of aerial biomass that can be allocated to stalks each day
      REAL SERo  !  max stalk elongation rate per day, under optimal temperature and water conditions (cm/d)
      REAL SSH  ! specific stalk height' = cm length/g (over 1 m2 for a mature crop i.e. 8-15 stalks)
      REAL SXRo ! May 2021
      REAL SSV ! May 2021
      REAL TBase_SXR ! May 2021
      REAL TOpt_SXR ! May 2021
      REAL TFin_SXR ! May 2021
      REAL lai_sen_light  ! APSIM-Sugar maximum GLAI for light interception reasons / LAI at which senescence starts (m2/m2)
      REAL sen_light_slope  ! APSIM-Sugar maximum GLAI for light interception reasons / LAI at which senescence starts (m2/m2)
      REAL ResExFrac ! May 2021
      REAL FT_SXR ! May 2021
      REAL SLA_prop ! May 2021
      REAL SLA_prop_raw ! JULY 2021


c     JvdM: compiled from genetic parameters
      REAL AvgSLA ! AVG (min+max/2) specific leaf area (cm2/g)


c     JvdM: RGR model local state and daily rate variables
      REAL TAVE  ! daily mean temperature
      REAL ADM  ! Above-ground dry mass t/ha (sum of components)
      REAL CANma  ! Conversion from leaf to dry mass
      REAL CTT_LFEM  ! Cumulative thermal time for leaf appearance °Cd
      REAL dADM ! MAY 2021
      REAL dADMSource  ! Daily change in Above-ground dry mass
      REAL dSourceLimit !
      REAL dGLA  ! New daily green leaf area m2/m2
      REAL dGLA_pot  ! Potential (source-unlimited) change in green leaf area m2/m2
      REAL dGLAI  ! Daily change in green leaf area index, new growth - senescence
      REAL dLeafDM  ! Daily change in total leaf dry mass t/ha/d
      REAL dlt_slai_light  ! Daily leaf area senescence m2/m2/d
      REAL dRootDM  ! Daily change in Root dry mass
      REAL dSDM  ! Daily change in stalk dry mass t/ha/d
      REAL dSenDM  ! Daily leaf area senescence m2/m2/d
      REAL dRootDM_Sett ! MAY 2021
      REAL dSFibDM  ! Daily change in stalk fibre dry mass t/ha/d
      REAL dSourceShortfall  ! daily sink - source, t/ha/d
      REAL LeafSourceShortfall  ! daily leaf source shortfall
      REAL StalkSourceShortfall  ! daily stalk source shortfall
      REAL AdditionalSource  ! additional source that could not be allocated (t/ha)
      REAL dAdditionalSource  ! daily additional source that could not be allocated (t/ha/d)
      REAL SourceShortfall
      REAL dSSuc  ! Daily change in stalk sucrose t/ha/d
      REAL dTTLf  ! Daily thermal time accumulation for driving leaf appearance °Cd
      REAL FIinter  ! Inter-row PAR fractional interception
      REAL FT_LAI  ! Temperature factor for LAI growth
      REAL FT_photos  ! Daily temperature factor for photosynthesis
      REAL GLAI  ! Green leaf area index m2/m2
      REAL GLeafDM  ! Green leaf canopy dry mass t/ha
      REAL CarbReserveDM  ! Sett dry mass (energy source) t/ha
      REAL KePAR  ! PAR canopy extinction coefficient
      REAL LAIsen  ! Daily GLAI senesced m2/m2/d
      REAL LeafDM  ! Total (living + dead) leaf dry mass t/ha
      REAL LeafPI  ! Current leaf phyllocron interval °Cd
      REAL CTT_LFEM_prev ! CTT when the most recent leaf appeared
      REAL LeafSource  ! Leaf carbon supply (source strength) t/ha/d
      REAL LeafSink  ! Leaf carbon demand (sink strength) t/ha/d
      REAL LFNUM_OSG  ! The number of leaves per shoot at which onset of stalk growth started l/s
      REAL LFNUM_SUCStart  ! Leaf number per stalk at which sucrose accumulation can start. l/s
      REAL NumLF  ! Number of leaves per shoot (reference) l/s
      REAL NumLFprev ! added February 2021
      REAL dNumLF  ! change in number of leaves l/s
      REAL RootFrac  ! Daily fraction of biomass allocated to roots, g/g
      REAL RootDM  ! Root dry mass
      REAL RGR_LAI_max  ! Today's maximum relative LAI growth rate
      REAL RGR_LAI_FTSW  ! Today's maximum relative LAI growth rate, considering temperature and soil moisture
      REAL RGR_LAI_act  ! Today's actual relative LAI growth rate, considering source limitations
      REAL SDM  ! Stalk dry mass t/ha
      REAL SenDM  ! Dry mass of senesced leaves t/ha
      REAL SXRpot  ! Stalk elongation rate cm/d
      REAL RSRlai  ! Daily fraction of GLAI senesced
      REAL PIAvgSSRv2  ! leaf PI average source:sink ratio (2nd approach)
      REAL SLA  ! specific leaf area (cm2/g) (daily / instantaneous value)
      REAL SLA_avg  ! Average SLA (cm2/g)
      REAL Source  ! Daily biomass increase (source strength) t/ha/d
      REAL StalkSource  ! stalk source strength
      REAL OSGfrac  ! fraction of shoots that have started stalk elongation
      REAL SPF_act  ! actual stalk partitioning fracton g/g
      REAL STKFRAC_avg  ! average stalk part frac, i.e. SDM/ADM
      REAL StalkSink  ! Stalk fibre growth sink strength t/ha/d
      REAL StalkVolume  ! stalk volume (t water / ha, equiv to fresh cane yield)
      REAL SFibDM  ! Stalk fibre dry mass t/ha
      REAL SWDF2  ! Soil water deficit (stress) factor affecting expansive growth
      REAL SUCDM  ! stalk sucrose (t/ha)
      REAL TotalDM  ! Total crop dry mass t/ha
      REAL vADM  ! Verify ADM t/ha
      REAL vSource  ! Verify Source t/ha
      REAL ADMSource  ! sum of ADM source
      REAL dExcessLeafSource ! MAY 2021
      REAL dReserveUsed ! MAY 2021
      REAL TotShortFall ! MAY 2021
      REAL LeafAllocRes ! MAY 2021
      REAL StalkAllocRes ! MAY 2021
      REAL dTotalSourceShortfall ! MAY 2021
      REAL ExtraSource  ! MAY 2021
      REAL PINumDays ! MAY 2021
      REAL PIsumSource ! MAY 2021
      REAL PIsumSink ! MAY 2021


c     to be used in seasinit      
      REAL LeafTT
      REAL LeafTTmax
      REAL StalkSucPart
      REAL SXR_EffTemp
      REAL Suc_LfNum_Max
      REAL, DIMENSION(5) :: LeafTTx, LeafTTy ! array
      REAL, DIMENSION(6) :: FT_photosx, FT_photosy ! array
      REAL, DIMENSION(5) :: FT_LAIx,FT_LAIy ! array
      REAL, DIMENSION(2) :: KePARx,KePARy ! array
      REAL, DIMENSION(3) :: RelativeSLAx,RelativeSLAy ! array
      REAL, DIMENSION(2) :: StalkSucPartx, StalkSucParty ! array
      REAL, DIMENSION(5) :: SXR_EffTempx, SXR_EffTempy ! array
      REAL TABEX
      REAL SWDF1
      
      REAL PIAvgSSR

      REAL dMaxReserveSource
      REAL dSourceFromSett
      REAL RootSource
      REAL SettDM

c     need for SPF calculation in rate calculations.
      REAL SPF
c     Need to fix digits to avoid overrun: REAL(8)
      REAL(8) inx

c     Emergence date (forced / override sim)
      INTEGER F_EDATE ! forced emergence date
      REAL F_EDATEIN ! same as above; reads as REAL


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

c     JvdM: JUN_2021 v13_6 additional functions
c     LeafMassFromArea
      REAL LMFAra
c     LeafAreaFromMass
      REAL LAFMa
c     SLAFromAreaAndMass
      REAL SLAFAM

      REAL LeafSink_SLAmin, LeafSink_SLAmax
     
c     JvdM: JULY_2021 v13_8 additional rowspacing var.
      REAL FII_c1, FII_c2
      REAL MaxResCont
      REAL NumPI_SSR


      

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
        !WRITE(*,*) 'SC_RGR CALLED IN RUNINIT MODE...'


c     need for OSGfrac calculation in rate calculations
      inx = 0.0
     

c     Call the output routine, to initialise output
c     ::::::::::
      CALL SC_RGOUTPUT(CONTROL, WEATHER,
     &  SW, SoilProp,
     &  YRPLT, CELLSE_DM,
     & TAVE, 
     & ADM, CANma, CTT_LFEM, dADMSource, 
     & dGLA, dGLA_pot, dGLAI, dLeafDM, 
     & dlt_slai_light, dRootDM, 
     & dSDM, dSenDM, dSFibDM, dSourceShortfall, 
     & LeafSourceShortfall, StalkSourceShortfall, 
     & AdditionalSource, dAdditionalSource, 
     & dSSuc, dTTLf, FIinter, FT_LAI, 
     & FT_photos, GLAI, GLeafDM, 
     & CarbReserveDM, KePAR, LAIsen, 
     & LeafDM, LeafPI, LeafSource, 
     & LeafSink, LFNUM_OSG, LFNUM_SUCStart, 
     & NumLF, dNumLF, RootFrac, 
     & RootDM, RGR_LAI_max, RGR_LAI_FTSW, 
     & RGR_LAI_act, SDM, SenDM, 
     & SXRpot, RSRlai, PIAvgSSRv2, 
     & SLA, SLA_avg, Source, StalkSource, 
     & OSGfrac, SPF_act, STKFRAC_avg, 
     & StalkSink, StalkVolume, 
     & SFibDM, SWDF2, SUCDM, TotalDM, 
     & vADM, vSource, ADMSource, 
     & SourceShortfall, SettDM)  

     


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
            !WRITE(*,*) 'SC_RGR CALLED IN SEASINIT MODE...'


      CALL GET_CULTIVAR_COEFF(InitialGLAI, 'InitialGLAI', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TBase_LAI, 'TBase_LAI', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TOpt_LAI, 'TOpt_LAI', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TFin_LAI, 'TFin_LAI', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TBase_Photos, 'TBase_Photos', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TOpt1_Photos, 'TOpt1_Photos', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TOpt2_Photos, 'TOpt2_Photos', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TFin_Photos, 'TFin_Photos', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TBase_LFAPP, 'TBase_LFAPP', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TOpt_LFAPP, 'TOpt_LFAPP', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TFin_LFAPP, 'TFin_LFAPP', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(LeafPI1, 'LeafPI1', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(LeafPI2, 'LeafPI2', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(KeMin, 'KeMin', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(KeMax, 'KeMax', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(KeMaxLf, 'KeMaxLf', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(SLAMin, 'SLAMin', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(SLAMax, 'SLAMax', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(Suc_LfNum_Delay, 'Suc_LfNum_Delay'
     & , CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(RGRglaiMin, 'RGRglaiMin', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(RGRglaiMax, 'RGRglaiMax', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(RUEo, 'RUEo', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(MAX_ROOTPF, 'MAX_ROOTPF', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(AvRootDMFrac, 'AvRootDMFrac', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(APFMX, 'APFMX', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(PCB, 'PCB', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(RUE_FT_c, 'RUE_FT_c', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(FIinOSG, 'FIinOSG', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(OSG_log_c1, 'OSG_log_c1', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(STKPFmax, 'STKPFmax', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(SXRo, 'SXRo', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(SSV, 'SSV', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TBase_SXR, 'TBase_SXR', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TOpt_SXR, 'TOpt_SXR', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(TFin_SXR, 'TFin_SXR', CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(lai_sen_light, 'lai_sen_light', 
     & CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(sen_light_slope, 'sen_light_slope', 
     & CONTROL, cERROR)
      CALL GET_CULTIVAR_COEFF(ResExFrac, 'ResExFrac', CONTROL, cERROR)


c     set lai_sen_light to test model behaviour
!      lai_sen_light = 2.3

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
          ! if RATNUM 0: PLANTCROP, OTHERWISE RATOON

c     Read reported date of emergence
            CALL FIND_INP(F_EDATEIN, 'EDATE', Control)
            F_EDATE = F_EDATEIN / 1 ! convert to integer


c	Water logging stress 
      PORMIN = 0.05 
      CALL GET_SPECIES_COEFF(PORMIN,'PORM', CONTROL, SPC_ERROR) 
      ! init waterlogging stress
      !WRITE(*, '(A, F10.5)') 'PORMIN is ', PORMIN
      DO I=1, NL
        TSS(I) = 0.0
      ENDDO
      SATFAC = 0.0

      CALL SC_RGOUTPUT(CONTROL, WEATHER,
     &  SW, SoilProp,
     &  YRPLT, CELLSE_DM,
     & TAVE, 
     & ADM, CANma, CTT_LFEM, dADMSource, 
     & dGLA, dGLA_pot, dGLAI, dLeafDM, 
     & dlt_slai_light, dRootDM, 
     & dSDM, dSenDM, dSFibDM, dSourceShortfall, 
     & LeafSourceShortfall, StalkSourceShortfall, 
     & AdditionalSource, dAdditionalSource, 
     & dSSuc, dTTLf, FIinter, FT_LAI, 
     & FT_photos, GLAI, GLeafDM, 
     & CarbReserveDM, KePAR, LAIsen, 
     & LeafDM, LeafPI, LeafSource, 
     & LeafSink, LFNUM_OSG, LFNUM_SUCStart, 
     & NumLF, dNumLF, RootFrac, 
     & RootDM, RGR_LAI_max, RGR_LAI_FTSW, 
     & RGR_LAI_act, SDM, SenDM, 
     & SXRpot, RSRlai, PIAvgSSRv2, 
     & SLA, SLA_avg, Source, StalkSource, 
     & OSGfrac, SPF_act, STKFRAC_avg, 
     & StalkSink, StalkVolume, 
     & SFibDM, SWDF2, SUCDM, TotalDM, 
     & vADM, vSource, ADMSource, 
     & SourceShortfall, SettDM)   
     



c      set initial GLAI
c      set initialGLAI to a higher value for ratoon crops
c      1.62 from Bezuidenhout
      IF (RATNUM .GT. 0.0) THEN
            InitialGLAI = InitialGLAI * 1.3
      ELSE
            InitialGLAI = InitialGLAI 
      ENDIF

c      Initial mass of energy in setts (which becomes roots) (t/ha)
      SettDM = 0.0

c      initialise the "new" reserve (development of a stool)
      CarbReserveDM = 0.0
      
c       ---- "lookup" functions --------
c      These functions interpolate linearly between points.
      FT_photosx = (/ 0.0, TBase_Photos, TOpt1_Photos, 
     & TOpt2_Photos, TFin_Photos, 0.0 /)
      FT_photosy = (/ 0.0,          0.0,          1.0, 
     &          1.0,         0.0, 0.0 /)

      !WRITE(*,'(A,F5.2)') "FT_photosx: ", FT_photosx(2)
      !WRITE(*,'(A,F5.2)') "MAX OF FT_photosy: ", MAXVAL(FT_photosy)

c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c     format: TABEX(yArray, xArray, valueToBeEstimated, dimensions of Arrays)
c     FT_photos = TABEX(FT_photosy, FT_photosx, TAVE, 6)

      !WRITE(*,'(A,F8.3)') "TFin_LAI: ", TFin_LAI

c     Max relative growth rate of GLAI
      FT_LAIx= (/ 0.0, TBase_LAI, TOpt_LAI, TFin_LAI, 0.0 /)
      FT_LAIy= (/ 0.0,       0.0,      1.0,      0.0, 0.0 /)

c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c     format: TABEX(xArray, yArray, valueToBeEstimated, dimensions of Arrays)
c      FT_LAI = TABEX(FT_LAIy, FT_LAIx, TAVE, 5)
      
c      WRITE(*,'(A,F8.3)') "FT_LAI: ", FT_LAI
      
c      radiation extinction coefficient:
c      linear interpolation between KeMin when there are no leaves
c      and KeMax when there are KeMaxLf leaves
      KePARx=(/ 0.0, KeMaxLf /)
      KePARy=(/ KeMin, KeMax /)

c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c     format: TABEX(yArray, xArray, valueToBeEstimated, dimensions of Arrays)
c      KePAR = TABEX(KePARy, KePARx, NumLF, 2)
      
c      WRITE(*,'(A,F8.3)') "KePAR: ", KePAR

c      OSG fraction - NOT IN USE
c       FI_OSGx = (/ FI_OSG - 0.05, min(1.0, FI_OSG + 0.05) /)
c      FI_OSGy = (/ 0.0, 1.0 /)

c     call TABEX function later (daily loop, rates of change) with arrays defined above:
C      OSGfrac = TABEX(FI_OSGy, FI_OSGx, FIinter, 2)

c     define the leaf appearance dTT function
      LeafTTmax = TOpt_LFAPP - TBase_LFAPP
      
c     prepare arrays for TABEX function
      LeafTTx = (/ 0.0, TBase_LFAPP, TOpt_LFAPP, TFin_LFAPP,  60.0 /)
      LeafTTy = (/ 0.0,         0.0,  LeafTTmax,        0.0,  0.0 /)

c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c     format: TABEX(yArray, xArray, valueToBeEstimated, dimensions of Arrays)
c      LeafTT = TABEX(LeafTTy, LeafTTx, TAVE, 5)

c      WRITE(*,'(A,F10.4)') "LeafTT: ", LeafTT

c      Stalk elongation dTT function
       SXR_EffTempx = (/ 0.0, TBase_SXR, TOpt_SXR, TFin_SXR, 60.0 /)
       SXR_EffTempy = (/ 0.0,       0.0,      1.0,      0.0,  0.0 /)
c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c      SXR_EffTemp = TABEX(SXR_EffTempy, SXR_EffTempx, TAVE, 5)

c      define the function relating relative specific leaf area to source sink ratio
      AvgSLA =  0.5 * SLAMax + 0.5 * SLAMin

      !RelativeSLAx = (/    0.5,    1.5 /)
      !RelativeSLAy = (/ SLAMax, SLAMin /)
c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c      RelativeSLA = TABEX(RelativeSLAy, RelativeSLAx, PIAvgSSR, 3)

c      define function describing transition from stalk fibre only to stalk fibre + sucrose
c      this starts at mpars$Suc_LfNum_Delay (leaf 3 or 4 or whatever after OSG)
      Suc_LfNum_Max = 7.0  ! move to parameters sometime

      StalkSucPartx = (/ 0.0, Suc_LfNum_Max /)
      StalkSucParty = (/ 0.0, 1.0 /)
c     call TABEX function later (daily loop, rates of change) with arrays defined above:
c      StalkSucPart = TABEX(StalkSucPartx, StalkSucParty, NumLF-LFNUM_SUCStart , 2)

c     convert input FInOSG to FI_OSG using relationships based on data from 
c     Singels and Smit 2009 (wagon wheel)
c     print(paste0('Row-spacing is ', RowSpc))
      
      FII_c1 =  0.07 * ROWSPC + 0.33
      FII_c2 = 0.28 * ROWSPC + 0.65
      FI_OSG = MIN(0.98, exp((FIinOSG - FII_c2)/FII_c1))
      WRITE(*,'(A,F8.3)') "FI_OSG: ", FI_OSG 
      

      
c       state variables (set to 0 or other initial value)
      GLAI      = InitialGLAI ! initial green leaf area index
      LeafPI    = LeafPI1  ! Current Leaf phyllocron interval, °Cd
      ADM = 0.0
      ADMSource = 0.0
      CANma = 0.0
      CTT_LFEM = 0.0
      dADM = 0.0
      dADMSource = 0.0
      dGLA = 0.0
      dGLA_pot = 0.0
      dGLAI = 0.0
      dLeafDM = 0.0
      dlt_slai_light = 0.0
      dRootDM = 0.0
      dSDM = 0.0
      dSenDM = 0.0
      !dRootDM_Sett = 0.0
      dSFibDM = 0.0
      dSSuc = 0.0
      dTTLf = 0.0
      FIinter = 0.0
      FT_LAI = 0.0
      FT_photos = 0.0
      GLeafDM = 0.0
      KePAR = 0.0
      LAIsen = 0.0
      LeafDM = 0.0
      LeafSink = 0.0
      LFNUM_OSG = 0.0
      LFNUM_SUCStart = 0.0
      NumLF = 1.0 ! updated Apil 2021
      NumLFprev = 1.0 ! added Feb 2021
      RootDM = 0.0
      RGR_LAI_max = 0.0
      RGR_LAI_act = 0.0
      RGR_LAI_FTSW = 0.0
      SDM = 0.0
      SenDM = 0.0
      SXRpot = 0.0
      !SER = 0.0
      RSRlai = 0.0
      !SLSR = 0.0
      Source = 0.0
      OSGfrac = 0.0
      !SPF = 0.0
      StalkSink = 0.0
      StalkVolume = 0.0
      SFibDM = 0.0
      SWDF2 = 0.0
      SUCDM = 0.0
      TotalDM = 0.0
      vADM = 0.0
      vSource = 0.0
      RootFrac = 0.0
      !dMaxReserveSource = 0.0
      !dSourceFromSett = 0.0
      !SettDM = 0.0
      !RootSource = 0.0
      !FT_SXR = 0.0
      !dExcessLeafSource = 0.0
      !SLA_prop = 0.0
      !dReserveUsed = 0.0
      !TotShortFall = 0.0
      !LeafAllocRes = 0.0
      !StalkAllocRes = 0.0
      !dTotalSourceShortfall = 0.0
      !ExtraSource = 0.0
      

  

c      leaf phyllocron interval starts at minimum value, not 0
       LeafPI = LeafPI1
c       CTT when the most recent leaf appeared
       CTT_LFEM_prev = 0.0

       SLA = AvgSLA

c      variables for keeping track of average source sink ratio over 
c      the duration of a leaf phyllocron interval
      PINumDays = 0.0
       !PIAvgSSR <- 2.0
      PIAvgSSRv2 = 2.0
      PIsumSource = 0.0
      PIsumSink = 1.0

c       place holders for soil water stress deficit (stress) factors
       SWDF1 = 1.0
       SWDF2 = 1.0

c       number of leaves per shoot at start of stalk growth (i.e. when SPF > 0.001)
       LFNUM_OSG = 0.0
c       number of leaves at which sucrose accumulation can start
c       (this is recalculated at runtime)
       LFNUM_SUCStart = 20.0

c       unallocated source (sink-limited growth)
       AdditionalSource = 0.0
c       accummulated source shortfall (sink > source)
       SourceShortfall = 0.0

c      maximum energy content of reserve
       MaxResCont = 1.0

c       increment iteration counter
       !Inti = 0 ! PROBALY ONLY NEEDED IN R SCRIPT

      !WRITE(*,*) "End of SEASINIT"

c     END of SEASINIT
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     &----------------------------------------------------
c
c              DYNAMIC = RATE
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.RATE) THEN
            !WRITE(*,*) 'SC_RGR CALLED IN RATE MODE...'

            !WRITE(*,'(A,I8)') "YRDOY: ", Control%YRDOY

      IF (Control%YRDOY >= F_EDATE) THEN
            !WRITE(*,'(A,I8)') "F_EDATE: ", F_EDATE
            !WRITE(*,'(A,I8)') "YRDOY: ", Control%YRDOY
c     set rates of change to 0 at start            
      dADM = 0.0
      dGLA = 0.0
      dGLA_pot = 0.0
      dGLAI = 0.0
      dLeafDM = 0.0
!      dLeafPI = 0.0
      dlt_slai_light = 0.0
      dRootDM = 0.0
      dSDM = 0.0
      dSenDM = 0.0
      dSFibDM = 0.0
      dSSuc = 0.0
      dTTLf = 0.0
      dNumLF = 0.0
      dSourceLimit = 0.0
      dADMSource = 0.0
      
c     Average daily temperature
      TAVE = (0.5 * TMAX) + (0.5 * TMIN)
c     SRAD already defined
      
c     -----------  CALCULATE RATES OF CHANGE ------------------
c     ------------- Phenology-related -------------------
      
c     convert Specific leaf area input to t/ha/m2
c     (SLA might be modified at runtime)
c      CANma = 1.0 / SLA * 100.0


c     radiation extinction coefficient:
        IF (NumLF .GT. MAXVAL(KePARx)) THEN
            KePAR = MAXVAL(KePARy)
        ELSEIF (NumLF .LT. MINVAL(KePARx)) THEN
            KePAR = MINVAL(KePARy)
        ELSE
          KePAR = TABEX(KePARy, KePARx, NumLF, 2)
        ENDIF

      !WRITE(*,'(A,F8.3)') "KePAR: ", KePAR

c      leaf phyllocron interval (transition from short to long PI linked to transition to stalk growth):
      LeafPI = LeafPI1 + ((LeafPI2 - LeafPI1) * OSGfrac)

c     fractional interception of radiation
      inx = -1.0 * KePAR * GLAI
      FIinter = 1.0 - EXP(inx)
      inx = 0.0

      !WRITE(*,'(A,F8.3)') "FIinter: ", FIinter
      !WRITE(*,'(A,F8.3)') "OSG_log_c1: ", OSG_log_c1
      !WRITE(*,'(A,F8.3)') "(FIinter - FI_OSG): ", (FIinter - FI_OSG)

c      stalk partitioning fraction
      inx = -1 * OSG_log_c1 * (FIinter - FI_OSG)
      OSGfrac = MAX( 1 / (1 + EXP(inx)), OSGfrac)
      inx = 0.0
     
       !WRITE(*,'(A,F8.3)') "OSGfrac: ", OSGfrac

! c     stalk partitioning fraction
! c     can OSG_log_c1 be linked to LAI_RGR, i.e. a faster canopy --> faster transition?
!       inx = (-1.0 * OSG_log_c1 * (FIinter - FI_OSG))
!       SPF = MAX((STKPFmax / (1.0 + EXP(inx))), SPF)
!       inx = 0.0

!       WRITE(*,'(A,F8.3)') "SPF: ", SPF

c     how many leaves are there? (first primary cohort)
c     new leaf appears if LeafPI Cd have elapsed since current leaf appeared
      !NumLF = (CTT_LFEM / LeafPI) + 1.0
c     FLOOR returns an integer
      dNumLF = floor((CTT_LFEM - CTT_LFEM_prev) / LeafPI) * 1.0 !convert int -> real

      !WRITE(*,'(A,F8.3)') "dNumLF: ", dNumLF

c     if SPF > 0.001, record number of leaves:

      IF (LFNUM_OSG .LT. 1.0) THEN
            IF (OSGfrac .GT. 0.001) THEN
                  LFNUM_OSG = NumLF
                  LFNUM_SUCStart = LFNUM_OSG + Suc_LfNum_Delay
            ENDIF
      ENDIF

      !WRITE(*,'(A,F8.3)') "LFNUM_SUCStart: ", LFNUM_SUCStart
      !WRITE(*,'(A,F8.3)') "GLAI: ", GLAI

c     daily thermal time for leaf appearance
c     dTTLf = ifelse(GLAI>0.0, LeafTTf(TAVE), 0.0)
      IF (GLAI .GT. 0.0) THEN
            dTTLf =  TABEX(LeafTTy, LeafTTx, TAVE, 5)
      ELSE
            dTTLf =  0.0
      ENDIF
      
      !WRITE(*,'(A,F8.3)') "dTTLf: ", dTTLf

c     leaf phyllocron interval transitions from short to long as crop
c     transitions to stalk growth phase
!      dLeafPI = (LeafPI2 - LeafPI1) * (SPF / STKPFmax)

c     ----------- Photosynthesis (source strength) --------------
c     temperature factor for photosynthesis
      FT_photos = TABEX(FT_photosy, FT_photosx, TAVE, 6)

c     maximum source strength that can come from underground reserves:
c     described as maximum carb reserves used; any unused goes back into
c     the reserve each day
c     the rhizome has structure, absorbing 40% of the dry mass.  The remaining
c     60% is available for use, but the rate is limited to ResExFrac (a 
c     fraction) of the remaining total per day.
      dMaxReserveSource = CarbReserveDM * MaxResCont * ResExFrac
      dSourceFromSett = SettDM * ResExFrac
      
c     daily source strength    
c     source strength is increases by the average fraction allocated to roots.
c     RUEo is defined in terms of above-ground biomass per unit intercepted 
c     photosynthetically-active radiation (PAR).
      Source = (FT_photos * FIinter * RUEo * 0.01 * (1.0/(1.0 - AvRootDMFrac)) 
     & * SRAD * 0.5) + dSourceFromSett

      !WRITE(*,'(A,F8.3)') "Source: ", Source
c      divide up source strength into source pools for leaves and stalks:
      inx = -1.0 * PCB * TotalDM
      RootSource  = Source * min(MAX_ROOTPF, 1.0 - APFMX * ( 1 - EXP(inx)))
      inx = 0.0
      !WRITE(*,'(A,F8.3)') "RootSource: ", RootSource
      dADMSource  = Source - RootSource
      !WRITE(*,'(A,F8.3)') "dADMSource: ", dADMSource
      LeafSource  = (1.0 - OSGfrac * STKPFmax) * dADMSource
      !WRITE(*,'(A,F8.3)') "LeafSource: ", LeafSource
      StalkSource = OSGfrac * STKPFmax * dADMSource
      !WRITE(*,'(A,F8.3)') "StalkSource: ", StalkSource

c     ----------- Sink strengths --------------
c     daily allocation to root DM
      dRootDM = RootSource

c     calculate relative growth rate of green leaf area index
      RGR_LAI_max = RGRglaiMin + 
     & ((RGRglaiMax - RGRglaiMin) * (1.0 - OSGfrac))
      !WRITE(*,'(A,F8.3)') "RGR_LAI_max: ", RGR_LAI_max
c     potential (source-unlimited) growth in GLAI today
      FT_LAI = TABEX(FT_LAIy, FT_LAIx, TAVE, 5)
      !WRITE(*,'(A,F8.3)') "FT_LAI: ", FT_LAI

c     temperature- and water-stress limited RGR_LAI:
      RGR_LAI_FTSW = RGR_LAI_max * FT_LAI * SWDF2

c     source-unlimited change in LAI
      dGLA_pot = RGR_LAI_FTSW * GLAI

c     associated potential sink strengths for min and max SLA
      !WRITE(*,'(A,F8.3)') "dGLA_pot: ", dGLA_pot
      LeafSink_SLAmin = LMFAra(dGLA_pot, SLAMin)
      LeafSink_SLAmax = LMFAra(dGLA_pot, SLAMax)

      !WRITE(*,'(A,F8.3)') "LeafSink_SLAmin: ", LeafSink_SLAmin
      !WRITE(*,'(A,F8.3)') "LeafSink_SLAmax: ", LeafSink_SLAmax

c     stalk elongation rate (cm/d) for the average stalk
c     using SXRpot-specific cardinal temperatures for stalk growth
      FT_SXR = TABEX(SXR_EffTempy, SXR_EffTempx, TAVE, 5)
      SXRpot = SXRo * FT_SXR * SWDF2 * OSGfrac

      !WRITE(*,'(A,F8.3)') "TBase_SXR: ", TBase_SXR

      !WRITE(*,'(A,F8.3)') "FT_SXR: ", FT_SXR
      !WRITE(*,'(A,F8.3)') "SXRpot: ", SXRpot

c     stalk structural growth sink strength
c     maybe delay by a few PIs before sucrose can start accumulating?
c     SSV = 'specific stalk volume' 
      StalkSink = SXRpot * 1.0 / SSV * 0.01

c     ------------- resolution of source and sink strengths ----------
C     demand for leaf dry mass (source and sink-limited)
      dLeafDM = 0.0     

c     spare source after meeting leaf growth demand at minimum SLA (dense leaves)
      dExcessLeafSource = 0.0

c     amount of source the leaves would need to achieve full expansion
c     today at maxSLA (i.e. thinnest leaves)
      LeafSourceShortfall = 0.0

      IF (LeafSource .LE. LeafSink_SLAmax) THEN
c       Source-limited scenario --> leaf area growth is reduced  
        SLA = SLAMax
        dLeafDM = LeafSource
c        dGLA = dLeafDM * SLA * 0.01
c       how much reserve would be required to meet shortfall?
        LeafSourceShortfall = LeafSink_SLAmax - LeafSource
      
      ELSEIF (LeafSource .GT. LeafSink_SLAmin) THEN
c       Excess-source scenario  
        SLA  = SLAMin
        dLeafDM = LeafSink_SLAmin
c        dGLA =  dLeafDM * SLA *0.01
        dExcessLeafSource = LeafSource - dLeafDM
      
      ELSE 
c       Leaf source strength is equivalent to sink strength within the SLA range  
        dLeafDM = LeafSource  
c       select an SLA proportional to where source is relative to sink
c       if source is close to sink with SLA_min, then SLA is close to SLA_min
c       if source is close to sink with SLA_max, then SLA is close to SLA_max
        SLA_prop_raw = 1.0 - ((LeafSource - LeafSink_SLAmax)/(LeafSink_SLAmin - LeafSink_SLAmax))
        SLA_prop = MIN(LeafSink_SLAmax, MAX(LeafSink_SLAmin, SLA_prop_raw))
        SLA = SLAMin + ((SLAMax - SLAMin) * SLA_prop)
c        dGLA = LeafSource * SLA * 0.01
      ENDIF

c     add anything left over from leaves to stalks:
      StalkSource = StalkSource + dExcessLeafSource

c     calculate the actual RGR LAI:
      RGR_LAI_act = dGLA / GLAI

c     convert Specific leaf area input to t/ha/m2
c     (SLA might be modified at runtime)
      CANma = 1.0 / SLA * 100.0

c     change in stalk dry fibre mass, limited by:
c     1. Source partitioned to source
c     2. Stalk sink strength
      dSFibDM = min(StalkSource, StalkSink)

c     how much additional source is required for stalk fibre to achieve
c     full expansion today?
      StalkSourceShortfall = max(0.0, StalkSink - StalkSource)

c     The maximum amount of source available in reserve is dMaxReserveSource
      dReserveUsed = 0.0
      TotShortFall = LeafSourceShortfall + StalkSourceShortfall
      LeafAllocRes = 0.0
      StalkAllocRes = 0.0

c     if the shortfalls for leaf and stalk are both less than the available
c     reserve, we can add this and allow full expansion of leaves and stalks
    
      IF (TotShortFall .GT. 0.00001) THEN
        IF (TotShortFall .LE. dMaxReserveSource ) THEN
c         if the shortfall is less than the available reserve, meet
c         the shortfall for leaves and stalks
          dReserveUsed = TotShortFall
          LeafAllocRes = LeafSourceShortfall
          StalkAllocRes = StalkSourceShortfall
        ELSE 
c         if the shortfall greater than the available reserve, partitiong
c         according to OSGfrac and STKPFMAX:
c         v4
          LeafAllocRes = (1 - OSGfrac * STKPFmax) * dMaxReserveSource
          StalkAllocRes = OSGfrac * STKPFmax * dMaxReserveSource
          dReserveUsed = LeafAllocRes + StalkAllocRes
        ENDIF
      ENDIF

c      issue a warning if necessary
c      if (dReserveUsed > dMaxReserveSource) warning('Allocation of reserves exceeds available reserve.')

c     update remaining reserve:
      dMaxReserveSource = dMaxReserveSource - dReserveUsed

c     Total shortfall today:
      dTotalSourceShortfall = LeafSourceShortfall + StalkSourceShortfall 
     & - dReserveUsed

c     update leaf delta mass:
      dLeafDM = dLeafDM + LeafAllocRes

c     calculate actual source-limited change in leaf area today
      !dGLA = LAFMa(dLeafDM, SLA)  !dLeafDM * SLA * 0.01
      dGLA = dLeafDM * SLA * 0.01

c     update stalk dry mass delta:
      dSFibDM = dSFibDM + StalkAllocRes

c     allocation to sucrose is whatever is left over after roots, leaves and stalks
c     C demands are met.  Sucrose can only accumulate after a certain number of
c     internodes have appeared, i.e. at leaf number LFNUM_SUCStart
c     I suspect this necessary, otherwise sucrose content will be very high
c     shortly after OSG, where the canopy is photosynthesising strongly and
c     sink strengths are not that great during the transition to stalk growth.
      dSSuc = 0.0

      IF (NumLF .GE. LFNUM_SUCStart) THEN
c       allocate excess to sucrose if enough internodes have developed
c       note: stalkallocres is added to stalksource, but in principle
c       there should not be anything leftover if reserves had to be used
c       to support structural growth today.
        ExtraSource = StalkSource + StalkAllocRes - dSFibDM  !+ dExcessLeafSource
c       the amount allocated to sucrose is limited by internode maturity
c     ! dSSuc = ExtraSource * 
c     !&   TABEX(StalkSucParty, StalkSucPartx, (NumLF-LFNUM_SUCStart) , 2)
        IF ((NumLF-LFNUM_SUCStart) .GT. MAXVAL(StalkSucPartx)) THEN
            dSSuc = ExtraSource * MAXVAL(StalkSucParty)
        ELSEIF ((NumLF-LFNUM_SUCStart) .LT. MINVAL(StalkSucPartx)) THEN
            dSSuc = ExtraSource * MINVAL(StalkSucParty)
        ELSE
            dSSuc = ExtraSource * TABEX(StalkSucParty, StalkSucPartx, 
     &       (NumLF-LFNUM_SUCStart), 2)
        ENDIF

      ENDIF

c     was there a source shortfall?
      dSourceLimit = (LeafSink + StalkSink) - (dLeafDM + dSFibDM)

c     actual GLAI growth, based on new leaf dry mass and
c     yesterday's senescence
      dGLAI = dGLA - LAIsen

c     senescence: based on light environment only:
c     this basic concept comes from APSIM-Sugar
      IF (GLAI .GT. lai_sen_light) THEN
            RSRlai = sen_light_slope * (GLAI - lai_sen_light)
      ELSE
            RSRlai = 0.0
      ENDIF

c     daily senescence rate (using APSIM_Sugar variable name)
      dlt_slai_light = GLAI * RSRlai
      LAIsen = dlt_slai_light

c     mass senesced
c     ##todo: senesce based on SLA from when these leaves appeared.
      dSenDM  = dlt_slai_light * CANma


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

      ENDIF ! END OF IF (Control%YRDOY >= F_EDATE)

      !WRITE(*,*) "End of RATE MODE..."

c     END of RATE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = INTEGRATE
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.INTEGR) THEN
            !WRITE(*,*) 'SC_RGR CALLED IN INTEGRATE MODE...'

            IF (Control%YRDOY >= F_EDATE) THEN
      
c     update number of leaves
      NumLF = NumLF + dNumLF

c     update cumulative CTT since appearance of this leaf
      IF (dNumLF .GT. 0.0) THEN
c       reset
        CTT_LFEM_prev = CTT_LFEM
      ENDIF
    
c      check if the whole number of leaves changed?  In which case, calculate average
c      source/sink ratio over the previous phyllocron interval (PI)
      NumPI_SSR = 3.0

c       update cumulative source strength over this PI
      PIsumSource = PIsumSource + dADMSource
      PIsumSink = PIsumSink + LeafSink + StalkSink
c      update the number of days in this phyllocron interval
      PINumDays = PINumDays + 1.0     


c     integrate GLAI
      GLAI = GLAI + dGLAI

c     stalk volume (expressed in t water / ha, equiv to fresh cane yield)
      StalkVolume = StalkVolume + (dSFibDM * SSV)

c     thermal time for leaf appearance:
      CTT_LFEM = CTT_LFEM + dTTLf

c     integrate ADM
      ADMSource = ADMSource + dADMSource
      dADM = dLeafDM + dSFibDM + dSSuc
      ADM =  ADM + dADM
    
c     calculate the change today:
c      dAdditionalSource = dADMSource - (dLeafDM + dSFibDM + dSSuc)
      dAdditionalSource = MAX(0.0, dADMSource - (LeafSink_SLAmin + StalkSink + dSSuc))
      AdditionalSource = AdditionalSource + dAdditionalSource

c     SettDM
      SettDM = SettDM - dSourceFromSett
    
c     to what extent was source under-supplying the sink?
c     SourceShortfall <- SourceShortfall + dSourceLimit
      SourceShortfall = SourceShortfall + dTotalSourceShortfall
c     energy in sett / stool / rhizome
c     SettDM <- SettDM - dRootDM_Sett
      CarbReserveDM = CarbReserveDM - dReserveUsed + dAdditionalSource

c     integrate stalk dry mass
      SDM = SDM + dSFibDM + dSSuc

      !WRITE(*,'(A,F8.3)') "SDM: ", SDM

c     actual stalk partitioning fraction
      IF ((dLeafDM + dSFibDM + dSSuc) .NE. 0.0) THEN
        SPF_act = (dSFibDM + dSSuc) / (dLeafDM + dSFibDM + dSSuc)
      ENDIF

      !WRITE(*,'(A,F8.3)') "SPF_act: ", SPF_act


c     average stalk fraction of ADM:
      IF (ADM .NE. 0.0 ) THEN
        STKFRAC_avg = SDM / ADM
      ENDIF

      !WRITE(*,'(A,F8.3)') "STKFRAC_avg: ", STKFRAC_avg
    
c     stalk fibre dry mass
      SFibDM = SFibDM + dSFibDM
c     stalk sucrose (t/ha)
      SUCDM = SUCDM + dSSuc
      LeafDM = LeafDM + dLeafDM
c     integrate green leaf dry mass
      GLeafDM = GLeafDM + dLeafDM - dSenDM
c     calculate average SLA; clamp to known max
      !WRITE(*,'(A,F8.3)') "GLeafDM: ", GLeafDM
      !WRITE(*,'(A,F8.3)') "GLAI: ", GLAI
      ! check for GLeadFM > 0 otherwise SLAFAM function returns error
      IF(GLeafDM .NE. 0.0) THEN
        SLA_avg = MIN(SLAMax, SLAFAM(GLeafDM,  GLAI))
      ENDIF
      !WRITE(*,'(A,F8.3)') "SLA_avg: ", SLA_avg
    
c     integrate total dry mass
      TotalDM = TotalDM + Source
c     integrate root dry mass
      RootDM = RootDM + dRootDM

c     integrate senesced DM
      SenDM = SenDM + dSenDM
c     verify ADM
      vADM = GLeafDM + SenDM + dSFibDM + dSSuc
    
c     verify source
      vSource = dRootDM + dLeafDM + dSFibDM + dSSuc      



      ENDIF ! end ofIF (Control%YRDOY >= F_EDATE)

      !WRITE(*,*) "End of INTEGRATE MODE..."


c     END of INTEGRATE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
c              DYNAMIC = OUTPUT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.OUTPUT) THEN
        
            !WRITE(*,*) "SC_RGR CALLED IN OUTPUT MODE..."

            IF (Control%YRDOY >= F_EDATE) THEN

c     Write output of daily values:
c     :::::::::::::::::::::::::::::

      CALL SC_RGOUTPUT(CONTROL, WEATHER,
     &  SW, SoilProp,
     &  YRPLT, CELLSE_DM,
     & TAVE, 
     & ADM, CANma, CTT_LFEM, dADMSource, 
     & dGLA, dGLA_pot, dGLAI, dLeafDM, 
     & dlt_slai_light, dRootDM, 
     & dSDM, dSenDM, dSFibDM, dSourceShortfall, 
     & LeafSourceShortfall, StalkSourceShortfall, 
     & AdditionalSource, dAdditionalSource, 
     & dSSuc, dTTLf, FIinter, FT_LAI, 
     & FT_photos, GLAI, GLeafDM, 
     & CarbReserveDM, KePAR, LAIsen, 
     & LeafDM, LeafPI, LeafSource, 
     & LeafSink, LFNUM_OSG, LFNUM_SUCStart, 
     & NumLF, dNumLF, RootFrac, 
     & RootDM, RGR_LAI_max, RGR_LAI_FTSW, 
     & RGR_LAI_act, SDM, SenDM, 
     & SXRpot, RSRlai, PIAvgSSRv2, 
     & SLA, SLA_avg, Source, StalkSource, 
     & OSGfrac, SPF_act, STKFRAC_avg, 
     & StalkSink, StalkVolume, 
     & SFibDM, SWDF2, SUCDM, TotalDM, 
     & vADM, vSource, ADMSource, 
     & SourceShortfall, SettDM)  
     
            ENDIF ! end of IF (Control%YRDOY >= F_EDATE)

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
            !WRITE(*,*) 'SC_RGR CALLED IN FINAL MODE...'
            
            !WRITE(*,'(A,F8.3)') "ROWSPC: ", ROWSPC
            !WRITE(*,'(A,F8.3)') "FI_OSG: ", FI_OSG




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
c     END of SUBROUTINE SC_RGR
c     -----------------------------------------------------

c     FUNCTION LeafMassFromArea(Area, SLA)
      FUNCTION LMFAra(Area_f1, SLA_f1)
        IMPLICIT NONE
        REAL Area_f1, SLA_f1, LMFAra
        LMFAra = Area_f1 * 1.0/SLA_f1 * 100.0
c      RETURN 
      END FUNCTION LMFAra


c     FUNCTION LeafAreaFromMass(Mass, SLA)
      FUNCTION LAFMa(Mass_f2, SLA_f2)
          IMPLICIT NONE
          REAL Mass_f2, SLA_f2, LAFMa
          LAFMa = Mass_f2 * SLA_f2 * 0.01
c      RETURN 
      END FUNCTION LAFMa


c     FUNCTION SLAFromAreaAndMass(Mass, Area)
      FUNCTION SLAFAM(Mass_f3, Area_f3)
        IMPLICIT NONE
        REAL Mass_f3, Area_f3, SLAFAM
        SLAFAM = Area_f3 / (Mass_f3 * 0.01)
c      RETURN 
      END FUNCTION SLAFAM



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
