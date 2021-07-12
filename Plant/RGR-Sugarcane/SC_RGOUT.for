c     -----------------------------------------------------
c     -----------------------------------------------------

c     SUBROUTINE: SC_RGOUTPUT(), was SC_OPGROW()
c     ----------------------
c     This subroutine handles the output of growth
c     aspects (yields, roots, etc)

c     Matthew Jones, Feb 2008:
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     1. I have made a number of changes to this routine.
c     Output is now controlled by variable (runtime)
c     format strings.  Set column width with the 
c     VAR_WIDTH variable.  The code will automagically
c     adjust formatting (do not let VAR_WIDTH go below
c     4, otherwise this will clash with the floating
c     point decimals).  This might not work so well in 
c     a different compiler, because there is more 
c     whitespace in the format strings than strictly
c     necessary (or, um, correct).  It works with the
c     CVF compiler; if this does not work with Intel,
c     modify the code to output the format string, then
c     copy and paste this from the output window into the
c     code.  Then clean this up and use it as a hard-
c     coded format string.
c
c     2. I have also included a routine for right-aligning
c     text strings so that it still works with GBuild
c     
c     3. I have changed many of the ouput variables from
c     DSSAT-standard kg/ha measures to sugarcane-
c     conventional t/ha values.  I have updated the 
c     list of data codes in DATA.CDE with new variables
c     for these measures.
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SC_RGOUTPUT (CONTROL, WEATHER,
     & SW, SoilProp,
     & YRPLT, CELLSE_DM,
     &  TAVE,   ! daily mean temperature
     &  ADM,   ! Above-ground dry mass t/ha (sum of components)
     &  CANma,   ! Conversion from leaf to dry mass
     &  CTT_LFEM,   ! Cumulative thermal time for leaf appearance °Cd
     &  dADMSource,   ! Daily change in Above-ground dry mass
     &  dGLA,   ! New daily green leaf area m2/m2
     &  dGLA_pot,   ! Potential (source-unlimited) change in green leaf area m2/m2
     &  dGLAI,   ! Daily change in green leaf area index, new growth - senescence
     &  dLeafDM,   ! Daily change in total leaf dry mass t/ha/d
     &  dlt_slai_light,   ! Daily leaf area senescence m2/m2/d
     &  dRootDM,   ! Daily change in Root dry mass
     &  dSDM,   ! Daily change in stalk dry mass t/ha/d
     &  dSenDM,   ! Daily leaf area senescence m2/m2/d
     &  dSFibDM,   ! Daily change in stalk fibre dry mass t/ha/d
     &  dTotalSourceShortfall,   ! DAILY TotalSourceShortfall
     &  CarbReserveDM,   ! Sett dry mass (energy source) t/ha
     &  SourceShortfall,   ! SourceShortfall
     &  LeafSourceShortfall,   ! daily leaf source shortfall
     &  StalkSourceShortfall,   ! daily stalk source shortfall
     &  AdditionalSource,   ! additional source that could not be allocated (t/ha)
     &  dAdditionalSource,   ! daily additional source that could not be allocated (t/ha/d)
     &  dSSuc,   ! Daily change in stalk sucrose t/ha/d
     &  dTTLf,   ! Daily thermal time accumulation for driving leaf appearance °Cd
     &  FIinter,   ! Inter-row PAR fractional interception
     &  FT_LAI,   ! Temperature factor for LAI growth
     &  FT_photos,   ! Daily temperature factor for photosynthesis
     &  GLAI,   ! Green leaf area index m2/m2
     &  GLeafDM,   ! Green leaf canopy dry mass t/ha
     &  SettDM,   ! SettDM
     &  KePAR,   ! PAR canopy extinction coefficient
     &  LAIsen,   ! Daily GLAI senesced m2/m2/d
     &  LeafDM,   ! Total (living + dead) leaf dry mass t/ha
     &  LeafPI,   ! Current leaf phyllocron interval °Cd
     &  LeafSource,   ! Leaf carbon supply (source strength) t/ha/d
     &  LeafSink,   ! Leaf carbon demand (sink strength) t/ha/d
     &  LFNUM_OSG,   ! The number of leaves per shoot at which onset of stalk growth started l/s
     &  LFNUM_SUCStart,   ! Leaf number per stalk at which sucrose accumulation can start. l/s
     &  NumLF,   ! Number of leaves per shoot (reference) l/s
     &  dNumLF,   ! change in number of leaves l/s
     &  PIAvgSSRv2,   ! leaf PI average source:sink ratio (2nd approach)
     &  RootDM,   ! Root dry mass
     &  RootFrac,   ! Daily fraction of biomass allocated to roots, g/g
     &  RGR_LAI_max,   ! Today's maximum relative LAI growth rate
     &  RGR_LAI_FTSW,   ! Today's maximum relative LAI growth rate, considering temperature and soil moisture
     &  RGR_LAI_act,   ! Today's actual relative LAI growth rate, considering source limitations
     &  SDM,   ! Stalk dry mass t/ha
     &  SenDM,   ! Dry mass of senesced leaves t/ha
     &  SXR,   ! Stalk elongation rate cm/d
     &  SFibDM,   ! Stalk fibre dry mass t/ha
     &  slai_light_fac,   ! Daily fraction of GLAI senesced
     &  SLA,   ! specific leaf area (cm2/g) (daily / instantaneous value)
     &  SLA_avg,   ! Average SLA (cm2/g)
     &  Source,   ! Daily biomass increase (source strength) t/ha/d
     &  StalkSource,   ! stalk source strength
     &  OSGfrac,   ! fraction of shoots that have started stalk elongation
     &  SPF_act,   ! actual stalk partitioning fracton g/g
     &  STKFRAC_avg,   ! average stalk part frac, i.e. SDM/ADM
     &  StalkSink,   ! Stalk fibre growth sink strength t/ha/d
     &  StalkLength,   ! Stalk length (height to top visible dewlap) cm
     &  SWDF2,   ! Soil water deficit (stress) factor affecting expansive growth
     &  SUCDM,   ! stalk sucrose (t/ha)
     &  TotalDM,   ! Total crop dry mass t/ha
     &  vADM,   ! Verify ADM t/ha
     &  vSource,   ! Verify Source t/ha
     &  ADMSource)   ! sum of ADM source





     

      
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::


c     Define DSSAT composite variables:
c     [Taken from MZ_CERES.for]
      USE ModuleDefs
      USE ModuleData

c     Define CANEGRO composite variables:
      USE CNG_ModuleDefs

c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     No implicit typing - yuck!
      IMPLICIT NONE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     Maintain the value of internal variables between 
c     calls to this subroutine (supports the modular 
c     control structure)
      SAVE
c     DSSAT composite variables:
c     [Taken from MZ_CERES.for]
c     ::::::::::::::::::::::::::
      TYPE (ControlType) CONTROL
      TYPE (SoilType)    SOILPROP !chp not actually used yet
!      TYPE (SwitchType)  ISWITCH  !chp not actually used yet
!      Type (ResidueType) HARVRES  !chp not actually used yet
!      Type (ResidueType) SENESCE  !chp not actually used yet
      Type (WeatherType) WEATHER

c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     DSSAT inputs:
c     :::::::::::::
      REAL TMIN, TMAX

c     The file unit number for the growth output file.
      INTEGER NOUTDG
c     local run number
      INTEGER RUN
c     local ?
      CHARACTER*20 FILEIO

c     Soil water content
      REAL SW(NL)

c     CANEGRO variables:
c     ::::::::::::::::::
c     The CaneCrop 'object', an input parameter
c     :::::::::::::::::::::::::::::::::::::::::
      TYPE (CaneCropType) CaneCrop
c     The growth 'object'
      TYPE (GrothType)    Growth
c     The Partitioning object
      TYPE (PartType) Part
c     The output object:
      TYPE (OutType) Out
c     Water balance object:
      TYPE (WaterType) WaterBal

c     General inputs:
c     :::::::::::::::
c     General variables:
c     :::::::::::::::::::
c     local copy of CONTROL%DYNAMIC:
      INTEGER DYNAMIC
c     Filename
      CHARACTER*20 OFILE
c     Does the file exist, is this the first run?:
      LOGICAL FILE_EXISTS !, FIRST
c     Error status:
      INTEGER ERRNUM

c     Date/planting variables (local):
c     ::::::::::::::::::::::::::::::::
      INTEGER DAP, DAS, FROP, YRDOY, YRPLT, YEAR, DOY
      INTEGER TIMDIF  !, YRSIM
!      INTEGER MDATE, 
c     Intermediate output variables (to accom.
c     different units, etc)
c     ::::::::::::::::::::::
c     Stalk population per m2
      REAL STKPOPm2 

c     Cellulosic DM (t/ha)
      REAL, INTENT(IN) ::  CELLSE_DM

c     SC_RGR outputs
      REAL, INTENT(IN) :: TAVE !daily mean temperature
      REAL, INTENT(IN) :: ADM  ! Above-ground dry mass t/ha (sum of components)
      REAL, INTENT(IN) :: CANma  ! Conversion from leaf to dry mass
      REAL, INTENT(IN) :: CTT_LFEM  ! Cumulative thermal time for leaf appearance °Cd
      REAL, INTENT(IN) :: dADMSource  ! Daily change in Above-ground dry mass
      REAL, INTENT(IN) :: dGLA  ! New daily green leaf area m2/m2
      REAL, INTENT(IN) :: dGLA_pot  ! Potential (source-unlimited) change in green leaf area m2/m2
      REAL, INTENT(IN) :: dGLAI  ! Daily change in green leaf area index, new growth - senescence
      REAL, INTENT(IN) :: dLeafDM  ! Daily change in total leaf dry mass t/ha/d
      REAL, INTENT(IN) :: dlt_slai_light  ! Daily leaf area senescence m2/m2/d
      REAL, INTENT(IN) :: dRootDM  ! Daily change in Root dry mass
      REAL, INTENT(IN) :: dSDM  ! Daily change in stalk dry mass t/ha/d
      REAL, INTENT(IN) :: dSenDM  ! Daily leaf area senescence m2/m2/d
      REAL, INTENT(IN) :: dSFibDM  ! Daily change in stalk fibre dry mass t/ha/d
      REAL, INTENT(IN) :: dTotalSourceShortfall  ! DAILY TotalSourceShortfall
      REAL, INTENT(IN) :: CarbReserveDM  ! Sett dry mass (energy source) t/ha
      REAL, INTENT(IN) :: SourceShortfall  ! SourceShortfall
      REAL, INTENT(IN) :: LeafSourceShortfall  ! daily leaf source shortfall
      REAL, INTENT(IN) :: StalkSourceShortfall  ! daily stalk source shortfall
      REAL, INTENT(IN) :: AdditionalSource  ! additional source that could not be allocated (t/ha)
      REAL, INTENT(IN) :: dAdditionalSource  ! daily additional source that could not be allocated (t/ha/d)
      REAL, INTENT(IN) :: dSSuc  ! Daily change in stalk sucrose t/ha/d
      REAL, INTENT(IN) :: dTTLf  ! Daily thermal time accumulation for driving leaf appearance °Cd
      REAL, INTENT(IN) :: FIinter  ! Inter-row PAR fractional interception
      REAL, INTENT(IN) :: FT_LAI  ! Temperature factor for LAI growth
      REAL, INTENT(IN) :: FT_photos  ! Daily temperature factor for photosynthesis
      REAL, INTENT(IN) :: GLAI  ! Green leaf area index m2/m2
      REAL, INTENT(IN) :: GLeafDM  ! Green leaf canopy dry mass t/ha
      REAL, INTENT(IN) :: SettDM  ! SettDM
      REAL, INTENT(IN) :: KePAR  ! PAR canopy extinction coefficient
      REAL, INTENT(IN) :: LAIsen  ! Daily GLAI senesced m2/m2/d
      REAL, INTENT(IN) :: LeafDM  ! Total (living + dead) leaf dry mass t/ha
      REAL, INTENT(IN) :: LeafPI  ! Current leaf phyllocron interval °Cd
      REAL, INTENT(IN) :: LeafSource  ! Leaf carbon supply (source strength) t/ha/d
      REAL, INTENT(IN) :: LeafSink  ! Leaf carbon demand (sink strength) t/ha/d
      REAL, INTENT(IN) :: LFNUM_OSG  ! The number of leaves per shoot at which onset of stalk growth started l/s
      REAL, INTENT(IN) :: LFNUM_SUCStart  ! Leaf number per stalk at which sucrose accumulation can start. l/s
      REAL, INTENT(IN) :: NumLF  ! Number of leaves per shoot (reference) l/s
      REAL, INTENT(IN) :: dNumLF  ! change in number of leaves l/s
      REAL, INTENT(IN) :: PIAvgSSRv2  ! leaf PI average source:sink ratio (2nd approach)
      REAL, INTENT(IN) :: RootDM  ! Root dry mass
      REAL, INTENT(IN) :: RootFrac  ! Daily fraction of biomass allocated to roots, g/g
      REAL, INTENT(IN) :: RGR_LAI_max  ! Today's maximum relative LAI growth rate
      REAL, INTENT(IN) :: RGR_LAI_FTSW  ! Today's maximum relative LAI growth rate, considering temperature and soil moisture
      REAL, INTENT(IN) :: RGR_LAI_act  ! Today's actual relative LAI growth rate, considering source limitations
      REAL, INTENT(IN) :: SDM  ! Stalk dry mass t/ha
      REAL, INTENT(IN) :: SenDM  ! Dry mass of senesced leaves t/ha
      REAL, INTENT(IN) :: SXR  ! Stalk elongation rate cm/d
      REAL, INTENT(IN) :: SFibDM  ! Stalk fibre dry mass t/ha
      REAL, INTENT(IN) :: slai_light_fac  ! Daily fraction of GLAI senesced
      REAL, INTENT(IN) :: SLA  ! specific leaf area (cm2/g) (daily / instantaneous value)
      REAL, INTENT(IN) :: SLA_avg  ! Average SLA (cm2/g)
      REAL, INTENT(IN) :: Source  ! Daily biomass increase (source strength) t/ha/d
      REAL, INTENT(IN) :: StalkSource  ! stalk source strength
      REAL, INTENT(IN) :: OSGfrac  ! fraction of shoots that have started stalk elongation
      REAL, INTENT(IN) :: SPF_act  ! actual stalk partitioning fracton g/g
      REAL, INTENT(IN) :: STKFRAC_avg  ! average stalk part frac, i.e. SDM/ADM
      REAL, INTENT(IN) :: StalkSink  ! Stalk fibre growth sink strength t/ha/d
      REAL, INTENT(IN) :: StalkLength  ! Stalk length (height to top visible dewlap) cm
      REAL, INTENT(IN) :: SWDF2  ! Soil water deficit (stress) factor affecting expansive growth
      REAL, INTENT(IN) :: SUCDM  ! stalk sucrose (t/ha)
      REAL, INTENT(IN) :: TotalDM  ! Total crop dry mass t/ha
      REAL, INTENT(IN) :: vADM  ! Verify ADM t/ha
      REAL, INTENT(IN) :: vSource  ! Verify Source t/ha
      REAL, INTENT(IN) :: ADMSource  ! sum of ADM source







c     CANEGRO 3.5 variables:
c     ::::::::::::::::::::::
c     GROHEAD - this is an array variable that stores the 
c     text headings in the plantgro.out file.
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     MJ, Feb 2008: Grohead was previously an array of four long
c     strings.  However, this was a pain to maintain, so it is
c     now a 2-dimensional array with each column header stored
c     separately.
      CHARACTER*15 GROHEAD(4, 69)
      CHARACTER*15 GRO_OUT(69)


c     Loop counter variables
      INTEGER I, J, NUM_OVARS

!     CHP added 08Jan10
!     Print only if output switch (IDETG) is 'Y'
      Type (SwitchType) ISWITCH
      CHARACTER*1 IDETG

c     OUTPUT variables (as they are in PlantGro.out)
c     ::::::::::::::::::::::::::::::::::::::::::::::
!      INTEGER GSTD
c     Leaf area index, DM, leafDM, stalkDM, sucroseDM, rootDM, crop DM  
c     The DM values are in kg/ha
      REAL LAID, LGDMD,  SMDMD, SUCMD,  RWAD, RDMD, BADMD
c     Stalk fresh mass, t/ha
      REAL SMFMD
c     Number of leaf tips, number of green leaves:
      REAL T_SD, L_SD
c     Harvest index
      REAL SUID
c     Trash
      REAL LDDMD  
c     Total (living + dead) LAI
      REAL LAITD  
c     Water stress
      REAL WSPD
      REAL WSGD
c     Specific leaf area
      REAL SLAD  
c     Canopy height
      REAL CHTD
c     Root depth
      REAL  RDPD
c     Potential soil evaporation
      REAL EOSA
c     Total root length density (weighted by layer thickness)
      REAL RTLD
c     Fractional light interception
      REAL LIPD
c     Water stress in the top 30 cm of soil, affecting tiller popn
      REAL SWDF30
c     Photosynthesis
      REAL GROSSP, BRDMD, PARCE
c     Thermal time / heat units for emergence, stalk population and
c     Leaf extension respectively (deg C[ > Tbase].day)
      REAL TTEBC, TTSPC, TTLEC
c     Sucrose content (%) of dry and fresh mass respectively
      REAL SUDMD, SUFMD
!----------------------------------------------------------------------
c             DATA (for headings) - peviously from DSSAT 3.5, modified
c             to be like DSSAT 4, and then changed to individual
c             variable approach
!----------------------------------------------------------------------
c     Length of variable name string (excluding leading and 
c     trailing whitespace)
      INTEGER VLEN

c     How many spaces need to be skipped?
      INTEGER SKIP

c     How wide should the columns be spaced?
      INTEGER VAR_WIDTH

c     String equivalents of VLEN and SKIP
      CHARACTER*10 SKIP_STR, VLEN_STR, WIDTH_STR, NUM_OVARS_STR

c     Runtime format statement:
      CHARACTER*1024 FMT_STR, T_FMT_STR

c     General format statments for outputting heading comments
c     and daily values
      CHARACTER*100 G_FMT_STR, D_FMT_STR

c     String lengths:
      INTEGER TLEN, FLEN, ENDI

c     SUGARCANE-specific output variables:
c     (See C:\DSSAT45\DATA.CDE)
c     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     @CDE   LABEL           DESCRIPTION............................................  SYNONYMS
c     EOSA   Pot. soil evap. Average potential soil evaporation per day  (mm/day)     .
c     SMFMD  Stalk fresh massStalk (millable) fresh mass (t/ha)                                 .
c     WSTD   Tiller stress   Water stress affecting tiller population                 .
c     SHTD   Stalk height m  Stalk height (m)                                         .
c     RTLD   Total RLV (cm2) Total root length per cm2 (cm/cm2)                       .
c     LI%D   % intercpt      Canopy light interception (%)                            .
c     PGRD   Gross photos.   Gross photosynthesis rate (t/ha/day)                     .
c     BRDMD  Bioms increase  Biomass accumulation rate (kg/ha/day)                     .
c     LDGFD  Lodging         Extent of lodging                                        .
c     ! New codes, MJ Feb 2008 (mainly t/ha measures as conv. for sugarcane)
c     LGDMD  Green tops DM   Green leaf canopy + meristem dry mass(t/ha)['green tops'].
c     SMDMD  Stalk DM t/ha   Stalk (millable) dry mass, t/ha)                         .
c     SUCMD  Sucrose DM t/ha Sucrose dry mass (t/ha)                                  .
c     LDDMD  Trash DM  t/ha  Trash (residue) dry mass (t/ha)                          .
c     BADMD  Aerial DM t/ha  Aerial dry biomass (t/ha)                                .
c     LAIGD  Green LAI       Green leaf area index (m2/m2)                            .
c     RDMD   Root DM t/ha    Root dry mass (t/ha)                                     .
c     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


c      Headings for output variables:
c      ::::::::::::::::::::::::::::::
c      Date:
           DATA GROHEAD(1, 1) /'!DATE'/
           DATA GROHEAD(2, 1) /'! '/
           DATA GROHEAD(3, 1) /'!YEAR'/
           DATA GROHEAD(4, 1) /'@YEAR'/

c      Day of year:
           DATA GROHEAD(1, 2) /'Day '/
           DATA GROHEAD(2, 2) /'of'/
           DATA GROHEAD(3, 2) /'year'/
           DATA GROHEAD(4, 2) /'DOY'/

c      Days after 'simulation start' (I think)
           DATA GROHEAD(1, 3) /'Days'/
           DATA GROHEAD(2, 3) /'after'/
           DATA GROHEAD(3, 3) /'sim.start'/
           DATA GROHEAD(4, 3) /'DAS'/

c      Days after 'planting' (/harvest?)
           DATA GROHEAD(1, 4) /'Days'/
           DATA GROHEAD(2, 4) /'after'/
           DATA GROHEAD(3, 4) /'plant'/
           DATA GROHEAD(4, 4) /'DAP'/

c        daily mean temperature
         DATA GROHEAD(1, 5) /'DAILY'/
         DATA GROHEAD(2, 5) /'MEAN'/
         DATA GROHEAD(3, 5) /'TEMPERATURE'/
         DATA GROHEAD(4, 5) /'TAVE'/
 
c        Above-ground dry mass t/ha (sum of components)
         DATA GROHEAD(1, 6) /'ABOVE-GROUND'/
         DATA GROHEAD(2, 6) /'DRY'/
         DATA GROHEAD(3, 6) /'MASS'/
         DATA GROHEAD(4, 6) /'BADMD'/
 
c        Conversion from leaf to dry mass
         DATA GROHEAD(1, 7) /'CONVERSION'/
         DATA GROHEAD(2, 7) /'FROM'/
         DATA GROHEAD(3, 7) /'LEAF'/
         DATA GROHEAD(4, 7) /'CANma'/
 
c        Cumulative thermal time for leaf appearance °Cd
         DATA GROHEAD(1, 8) /'CUMULATIVE'/
         DATA GROHEAD(2, 8) /'THERMAL'/
         DATA GROHEAD(3, 8) /'TIME'/
         DATA GROHEAD(4, 8) /'CTT_LFEM'/
 
c        Daily change in Above-ground dry mass
         DATA GROHEAD(1, 9) /'DAILY'/
         DATA GROHEAD(2, 9) /'CHANGE'/
         DATA GROHEAD(3, 9) /'IN'/
         DATA GROHEAD(4, 9) /'dADMSource'/
 
c        New daily green leaf area m2/m2
         DATA GROHEAD(1, 10) /'NEW'/
         DATA GROHEAD(2, 10) /'DAILY'/
         DATA GROHEAD(3, 10) /'GREEN'/
         DATA GROHEAD(4, 10) /'dGLA'/
 
c        Potential (source-unlimited) change in green leaf area m2/m2
         DATA GROHEAD(1, 11) /'POTENTIAL'/
         DATA GROHEAD(2, 11) /'(SOURCE-UNLIMITED)'/
         DATA GROHEAD(3, 11) /'CHANGE'/
         DATA GROHEAD(4, 11) /'dGLA_pot'/
 
c        Daily change in green leaf area index, new growth - senescence
         DATA GROHEAD(1, 12) /'DAILY'/
         DATA GROHEAD(2, 12) /'CHANGE'/
         DATA GROHEAD(3, 12) /'IN'/
         DATA GROHEAD(4, 12) /'dGLAI'/
 
c        Daily change in total leaf dry mass t/ha/d
         DATA GROHEAD(1, 13) /'DAILY'/
         DATA GROHEAD(2, 13) /'CHANGE'/
         DATA GROHEAD(3, 13) /'IN'/
         DATA GROHEAD(4, 13) /'dLeafDM'/
 
c        Daily leaf area senescence m2/m2/d
         DATA GROHEAD(1, 14) /'DAILY'/
         DATA GROHEAD(2, 14) /'LEAF'/
         DATA GROHEAD(3, 14) /'AREA'/
         DATA GROHEAD(4, 14) /'dlt_slai_light'/
 
c        Daily change in Root dry mass
         DATA GROHEAD(1, 15) /'DAILY'/
         DATA GROHEAD(2, 15) /'CHANGE'/
         DATA GROHEAD(3, 15) /'IN'/
         DATA GROHEAD(4, 15) /'dRootDM'/
 
c        Daily change in stalk dry mass t/ha/d
         DATA GROHEAD(1, 16) /'DAILY'/
         DATA GROHEAD(2, 16) /'CHANGE'/
         DATA GROHEAD(3, 16) /'IN'/
         DATA GROHEAD(4, 16) /'dSDM'/
 
c        Daily leaf area senescence m2/m2/d
         DATA GROHEAD(1, 17) /'DAILY'/
         DATA GROHEAD(2, 17) /'LEAF'/
         DATA GROHEAD(3, 17) /'AREA'/
         DATA GROHEAD(4, 17) /'dSenDM'/
 
c        Daily change in stalk fibre dry mass t/ha/d
         DATA GROHEAD(1, 18) /'DAILY'/
         DATA GROHEAD(2, 18) /'CHANGE'/
         DATA GROHEAD(3, 18) /'IN'/
         DATA GROHEAD(4, 18) /'dSFibDM'/
 
c        DAILY TotalSourceShortfall
         DATA GROHEAD(1, 19) /'DAILY'/
         DATA GROHEAD(2, 19) /'TOTALSOURCESHORTFALL'/
         DATA GROHEAD(3, 19) /'NA'/
         DATA GROHEAD(4, 19) /'dTotalSourceShortfall'/
 
c        Sett dry mass (energy source) t/ha
         DATA GROHEAD(1, 20) /'SETT'/
         DATA GROHEAD(2, 20) /'DRY'/
         DATA GROHEAD(3, 20) /'MASS'/
         DATA GROHEAD(4, 20) /'CarbReserveDM'/
 
c        SourceShortfall
         DATA GROHEAD(1, 21) /'SOURCESHORTFALL'/
         DATA GROHEAD(2, 21) /'NA'/
         DATA GROHEAD(3, 21) /'NA'/
         DATA GROHEAD(4, 21) /'SourceShortfall'/
 
c        daily leaf source shortfall
         DATA GROHEAD(1, 22) /'DAILY'/
         DATA GROHEAD(2, 22) /'LEAF'/
         DATA GROHEAD(3, 22) /'SOURCE'/
         DATA GROHEAD(4, 22) /'LeafSourceShortfall'/
 
c        daily stalk source shortfall
         DATA GROHEAD(1, 23) /'DAILY'/
         DATA GROHEAD(2, 23) /'STALK'/
         DATA GROHEAD(3, 23) /'SOURCE'/
         DATA GROHEAD(4, 23) /'StalkSourceShortfall'/
 
c        additional source that could not be allocated (t/ha)
         DATA GROHEAD(1, 24) /'ADDITIONAL'/
         DATA GROHEAD(2, 24) /'SOURCE'/
         DATA GROHEAD(3, 24) /'THAT'/
         DATA GROHEAD(4, 24) /'AdditionalSource'/
 
c        daily additional source that could not be allocated (t/ha/d)
         DATA GROHEAD(1, 25) /'DAILY'/
         DATA GROHEAD(2, 25) /'ADDITIONAL'/
         DATA GROHEAD(3, 25) /'SOURCE'/
         DATA GROHEAD(4, 25) /'dAdditionalSource'/
 
c        Daily change in stalk sucrose t/ha/d
         DATA GROHEAD(1, 26) /'DAILY'/
         DATA GROHEAD(2, 26) /'CHANGE'/
         DATA GROHEAD(3, 26) /'IN'/
         DATA GROHEAD(4, 26) /'dSSuc'/
 
c        Daily thermal time accumulation for driving leaf appearance °Cd
         DATA GROHEAD(1, 27) /'DAILY'/
         DATA GROHEAD(2, 27) /'THERMAL'/
         DATA GROHEAD(3, 27) /'TIME'/
         DATA GROHEAD(4, 27) /'dTTLf'/
 
c        Inter-row PAR fractional interception
         DATA GROHEAD(1, 28) /'INTER-ROW'/
         DATA GROHEAD(2, 28) /'PAR'/
         DATA GROHEAD(3, 28) /'FRACTIONAL'/
         DATA GROHEAD(4, 28) /'LI%D'/
 
c        Temperature factor for LAI growth
         DATA GROHEAD(1, 29) /'TEMPERATURE'/
         DATA GROHEAD(2, 29) /'FACTOR'/
         DATA GROHEAD(3, 29) /'FOR'/
         DATA GROHEAD(4, 29) /'FT_LAI'/
 
c        Daily temperature factor for photosynthesis
         DATA GROHEAD(1, 30) /'DAILY'/
         DATA GROHEAD(2, 30) /'TEMPERATURE'/
         DATA GROHEAD(3, 30) /'FACTOR'/
         DATA GROHEAD(4, 30) /'FT_photos'/
 
c        Green leaf area index m2/m2
         DATA GROHEAD(1, 31) /'GREEN'/
         DATA GROHEAD(2, 31) /'LEAF'/
         DATA GROHEAD(3, 31) /'AREA'/
         DATA GROHEAD(4, 31) /'LAIGD'/
 
c        Green leaf canopy dry mass t/ha
         DATA GROHEAD(1, 32) /'GREEN'/
         DATA GROHEAD(2, 32) /'LEAF'/
         DATA GROHEAD(3, 32) /'CANOPY'/
         DATA GROHEAD(4, 32) /'LGDMD'/
 
c        SettDM
         DATA GROHEAD(1, 33) /'SETTDM'/
         DATA GROHEAD(2, 33) /'NA'/
         DATA GROHEAD(3, 33) /'NA'/
         DATA GROHEAD(4, 33) /'SettDM'/
 
c        PAR canopy extinction coefficient
         DATA GROHEAD(1, 34) /'PAR'/
         DATA GROHEAD(2, 34) /'CANOPY'/
         DATA GROHEAD(3, 34) /'EXTINCTION'/
         DATA GROHEAD(4, 34) /'KePAR'/
 
c        Daily GLAI senesced m2/m2/d
         DATA GROHEAD(1, 35) /'DAILY'/
         DATA GROHEAD(2, 35) /'GLAI'/
         DATA GROHEAD(3, 35) /'SENESCED'/
         DATA GROHEAD(4, 35) /'LAIsen'/
 
c        Total (living + dead) leaf dry mass t/ha
         DATA GROHEAD(1, 36) /'TOTAL'/
         DATA GROHEAD(2, 36) /'(LIVING'/
         DATA GROHEAD(3, 36) /'+'/
         DATA GROHEAD(4, 36) /'LTDMD'/
 
c        Current leaf phyllocron interval °Cd
         DATA GROHEAD(1, 37) /'CURRENT'/
         DATA GROHEAD(2, 37) /'LEAF'/
         DATA GROHEAD(3, 37) /'PHYLLOCRON'/
         DATA GROHEAD(4, 37) /'LeafPI'/
 
c        Leaf carbon supply (source strength) t/ha/d
         DATA GROHEAD(1, 38) /'LEAF'/
         DATA GROHEAD(2, 38) /'CARBON'/
         DATA GROHEAD(3, 38) /'SUPPLY'/
         DATA GROHEAD(4, 38) /'LeafSource'/
 
c        Leaf carbon demand (sink strength) t/ha/d
         DATA GROHEAD(1, 39) /'LEAF'/
         DATA GROHEAD(2, 39) /'CARBON'/
         DATA GROHEAD(3, 39) /'DEMAND'/
         DATA GROHEAD(4, 39) /'LeafSink'/
 
c        The number of leaves per shoot at which onset of stalk growth started l/s
         DATA GROHEAD(1, 40) /'THE'/
         DATA GROHEAD(2, 40) /'NUMBER'/
         DATA GROHEAD(3, 40) /'OF'/
         DATA GROHEAD(4, 40) /'LFNUM_OSG'/
 
c        Leaf number per stalk at which sucrose accumulation can start. l/s
         DATA GROHEAD(1, 41) /'LEAF'/
         DATA GROHEAD(2, 41) /'NUMBER'/
         DATA GROHEAD(3, 41) /'PER'/
         DATA GROHEAD(4, 41) /'LFNUM_SUCStart'/
 
c        Number of leaves per shoot (reference) l/s
         DATA GROHEAD(1, 42) /'NUMBER'/
         DATA GROHEAD(2, 42) /'OF'/
         DATA GROHEAD(3, 42) /'LEAVES'/
         DATA GROHEAD(4, 42) /'NumLF'/
 
c        change in number of leaves l/s
         DATA GROHEAD(1, 43) /'CHANGE'/
         DATA GROHEAD(2, 43) /'IN'/
         DATA GROHEAD(3, 43) /'NUMBER'/
         DATA GROHEAD(4, 43) /'dNumLF'/
 
c        leaf PI average source:sink ratio (2nd approach)
         DATA GROHEAD(1, 44) /'LEAF'/
         DATA GROHEAD(2, 44) /'PI'/
         DATA GROHEAD(3, 44) /'AVERAGE'/
         DATA GROHEAD(4, 44) /'PIAvgSSRv2'/
 
c        Root dry mass
         DATA GROHEAD(1, 45) /'ROOT'/
         DATA GROHEAD(2, 45) /'DRY'/
         DATA GROHEAD(3, 45) /'MASS'/
         DATA GROHEAD(4, 45) /'RootDM'/
 
c        Daily fraction of biomass allocated to roots, g/g
         DATA GROHEAD(1, 46) /'DAILY'/
         DATA GROHEAD(2, 46) /'FRACTION'/
         DATA GROHEAD(3, 46) /'OF'/
         DATA GROHEAD(4, 46) /'RootFrac'/
 
c        Today_s maximum relative LAI growth rate
         DATA GROHEAD(1, 47) /'TODAY_S'/
         DATA GROHEAD(2, 47) /'MAXIMUM'/
         DATA GROHEAD(3, 47) /'RELATIVE'/
         DATA GROHEAD(4, 47) /'RGR_LAI_max'/
 
c        Today_s maximum relative LAI growth rate, considering temperature and soil moisture
         DATA GROHEAD(1, 48) /'TODAY_S'/
         DATA GROHEAD(2, 48) /'MAXIMUM'/
         DATA GROHEAD(3, 48) /'RELATIVE'/
         DATA GROHEAD(4, 48) /'RGR_LAI_FTSW'/
 
c        Today_s actual relative LAI growth rate, considering source limitations
         DATA GROHEAD(1, 49) /'TODAY_S'/
         DATA GROHEAD(2, 49) /'ACTUAL'/
         DATA GROHEAD(3, 49) /'RELATIVE'/
         DATA GROHEAD(4, 49) /'RGR_LAI_act'/
 
c        Stalk dry mass t/ha
         DATA GROHEAD(1, 50) /'STALK'/
         DATA GROHEAD(2, 50) /'DRY'/
         DATA GROHEAD(3, 50) /'MASS'/
         DATA GROHEAD(4, 50) /'SMDMD'/
 
c        Dry mass of senesced leaves t/ha
         DATA GROHEAD(1, 51) /'DRY'/
         DATA GROHEAD(2, 51) /'MASS'/
         DATA GROHEAD(3, 51) /'OF'/
         DATA GROHEAD(4, 51) /'LDDMD'/
 
c        Stalk elongation rate cm/d
         DATA GROHEAD(1, 52) /'STALK'/
         DATA GROHEAD(2, 52) /'ELONGATION'/
         DATA GROHEAD(3, 52) /'RATE'/
         DATA GROHEAD(4, 52) /'SXR'/
 
c        Stalk fibre dry mass t/ha
         DATA GROHEAD(1, 53) /'STALK'/
         DATA GROHEAD(2, 53) /'FIBRE'/
         DATA GROHEAD(3, 53) /'DRY'/
         DATA GROHEAD(4, 53) /'SFibDM'/
 
c        Daily fraction of GLAI senesced
         DATA GROHEAD(1, 54) /'DAILY'/
         DATA GROHEAD(2, 54) /'FRACTION'/
         DATA GROHEAD(3, 54) /'OF'/
         DATA GROHEAD(4, 54) /'slai_light_fac'/
 
c        specific leaf area (cm2/g) (daily / instantaneous value)
         DATA GROHEAD(1, 55) /'SPECIFIC'/
         DATA GROHEAD(2, 55) /'LEAF'/
         DATA GROHEAD(3, 55) /'AREA'/
         DATA GROHEAD(4, 55) /'SLAD'/
 
c        Average SLA (cm2/g)
         DATA GROHEAD(1, 56) /'AVERAGE'/
         DATA GROHEAD(2, 56) /'SLA'/
         DATA GROHEAD(3, 56) /'(CM2/G)'/
         DATA GROHEAD(4, 56) /'SLA_avg'/
 
c        Daily biomass increase (source strength) t/ha/d
         DATA GROHEAD(1, 57) /'DAILY'/
         DATA GROHEAD(2, 57) /'BIOMASS'/
         DATA GROHEAD(3, 57) /'INCREASE'/
         DATA GROHEAD(4, 57) /'Source'/
 
c        stalk source strength
         DATA GROHEAD(1, 58) /'STALK'/
         DATA GROHEAD(2, 58) /'SOURCE'/
         DATA GROHEAD(3, 58) /'STRENGTH'/
         DATA GROHEAD(4, 58) /'StalkSource'/
 
c        fraction of shoots that have started stalk elongation
         DATA GROHEAD(1, 59) /'FRACTION'/
         DATA GROHEAD(2, 59) /'OF'/
         DATA GROHEAD(3, 59) /'SHOOTS'/
         DATA GROHEAD(4, 59) /'OSGfrac'/
 
c        actual stalk partitioning fracton g/g
         DATA GROHEAD(1, 60) /'ACTUAL'/
         DATA GROHEAD(2, 60) /'STALK'/
         DATA GROHEAD(3, 60) /'PARTITIONING'/
         DATA GROHEAD(4, 60) /'SPF_act'/
 
c        average stalk part frac, i.e. SDM/ADM
         DATA GROHEAD(1, 61) /'AVERAGE'/
         DATA GROHEAD(2, 61) /'STALK'/
         DATA GROHEAD(3, 61) /'PART'/
         DATA GROHEAD(4, 61) /'STKFRAC_avg'/
 
c        Stalk fibre growth sink strength t/ha/d
         DATA GROHEAD(1, 62) /'STALK'/
         DATA GROHEAD(2, 62) /'FIBRE'/
         DATA GROHEAD(3, 62) /'GROWTH'/
         DATA GROHEAD(4, 62) /'StalkSink'/
 
c        Stalk length (height to top visible dewlap) cm
         DATA GROHEAD(1, 63) /'STALK'/
         DATA GROHEAD(2, 63) /'LENGTH'/
         DATA GROHEAD(3, 63) /'(HEIGHT'/
         DATA GROHEAD(4, 63) /'StalkLength'/
 
c        Soil water deficit (stress) factor affecting expansive growth
         DATA GROHEAD(1, 64) /'SOIL'/
         DATA GROHEAD(2, 64) /'WATER'/
         DATA GROHEAD(3, 64) /'DEFICIT'/
         DATA GROHEAD(4, 64) /'WSGD'/
 
c        stalk sucrose (t/ha)
         DATA GROHEAD(1, 65) /'STALK'/
         DATA GROHEAD(2, 65) /'SUCROSE'/
         DATA GROHEAD(3, 65) /'(T/HA)'/
         DATA GROHEAD(4, 65) /'SUCMD'/
 
c        Total crop dry mass t/ha
         DATA GROHEAD(1, 66) /'TOTAL'/
         DATA GROHEAD(2, 66) /'CROP'/
         DATA GROHEAD(3, 66) /'DRY'/
         DATA GROHEAD(4, 66) /'TotalDM'/
 
c        Verify ADM t/ha
         DATA GROHEAD(1, 67) /'VERIFY'/
         DATA GROHEAD(2, 67) /'ADM'/
         DATA GROHEAD(3, 67) /'T/HA'/
         DATA GROHEAD(4, 67) /'vADM'/
 
c        Verify Source t/ha
         DATA GROHEAD(1, 68) /'VERIFY'/
         DATA GROHEAD(2, 68) /'SOURCE'/
         DATA GROHEAD(3, 68) /'T/HA'/
         DATA GROHEAD(4, 68) /'vSource'/
 
c        sum of ADM source
         DATA GROHEAD(1, 69) /'SUM'/
         DATA GROHEAD(2, 69) /'OF'/
         DATA GROHEAD(3, 69) /'ADM'/
         DATA GROHEAD(4, 69) /'ADMSource'/
 



c          ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c          Number of output variables (important for output!)
           DATA NUM_OVARS /69/

c          Width of output columns:
           DATA VAR_WIDTH /21/
c         :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ~~~~~~~~ SUBROUTINE CODE  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     Initialisation for DYNAMIC = all
c     ::::::::::::::::::::::::::::::::
c     Init DYNAMIC:
      DYNAMIC = CONTROL%DYNAMIC
c     Init RUN
      RUN     = CONTROL%RUN
c     Fileio
      FILEIO  = CONTROL%FILEIO

c     Alert that sugarcane model is running
c      WRITE(*,*) 'Called the sugarcane output module'


c     -----------------------------------------------------
c
c              DYNAMIC = RUNINIT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      IF (DYNAMIC.EQ.RUNINIT) THEN


c     END of RUNINIT
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = SEASINIT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF (DYNAMIC.EQ.SEASINIT) THEN

          CALL GET(ISWITCH)
          IDETG = ISWITCH % IDETG
          IF (IDETG .NE. 'Y') RETURN

c     Open growth aspects output file:
c     ::::::::::::::::::::::::::::::::
c         Set file name:
          OFILE = 'PlantGro.OUT'

c         Get file unit number:
          CALL GETLUN('OUTG', NOUTDG)

c         Check that the file exists:
          FILE_EXISTS = .FALSE.
          INQUIRE(FILE=OFILE, EXIST=FILE_EXISTS)

c         Open the file
          IF (FILE_EXISTS) THEN
c             In append mode if the file already exists
              OPEN (UNIT=NOUTDG, FILE=OFILE, STATUS='OLD',
     &        IOSTAT=ERRNUM, POSITION='APPEND')
          ELSE
c             A new file if not existing
              OPEN (UNIT=NOUTDG, FILE=OFILE, STATUS='NEW',
     &        IOSTAT = ERRNUM)
              WRITE(NOUTDG,'("*GROWTH ASPECTS OUTPUT FILE")')
          ENDIF

c     Output a header (treatment / run info, etc)
c     ::::::::::::::::::::::::
c         Use the predefined one:
          CALL HEADER(SEASINIT, NOUTDG, RUN)

c     Now write the variable names:
c     :::::::::::::::::::::::::::::
c         (taken from DSSAT CANEGRO 3.5)
c         ::::::::::::::::::::::::::::::::::::::::::::
c         MJ, Feb 2008: changed from simply outputting the 
c         4 array elements (long strings) of GROHEAD to the
c         'finer-grained' one element per column and row
c         GROHEAD:
c         old:
c          WRITE (NOUTDG,2190) GROHEAD(1)
c          WRITE (NOUTDG,2190) GROHEAD(2)
c          WRITE (NOUTDG,2190) GROHEAD(3)
c          WRITE (NOUTDG,2190) GROHEAD(4)
c 2190     FORMAT (A)

c         new:
c              DO J=1,NUM_OVARS    
c                  GROHEAD(5, J) = '!--------------'  
c              ENDDO

c         Now, there is a problem.  GBuild expects column labels
c         to be right-aligned. Fortran, so far as I can find out,
c         only supports left-aligned text, BUT always right-
c         aligns Fx.x floating point numbers.  I do not know
c         why.  I have to then create a runtime format string
c         to correctly right-align the text in the headings:
          FMT_STR = '(A5, 1X, A3, '

          DO J=3, NUM_OVARS
c             Get length of variable name
              VLEN = LEN_TRIM(GROHEAD(4, J))
c             How many spaces need to be skipped?
              SKIP = VAR_WIDTH - VLEN
c             Create strings of these:
              WRITE(SKIP_STR, '(I3)') SKIP
              WRITE(VLEN_STR, '(I3)') VLEN
c              WRITE(*, '(A)') SKIP_STR, VLEN_STR

c             Add to format statement:
              IF (J .EQ. NUM_OVARS) THEN
c                 If this is format info for the last variable:
                  T_FMT_STR = TRIM(SKIP_STR) // 'X,'
     &                      // 'A' // TRIM(VLEN_STR)
     &                      // ')'
              ELSE
c                 For any variable
                  T_FMT_STR = TRIM(SKIP_STR) // 'X,'
     &                      // 'A' // TRIM(VLEN_STR)
     &                      // ','
              ENDIF

              TLEN = LEN_TRIM(T_FMT_STR) + 1
              FLEN = LEN_TRIM(FMT_STR) + 1
              ENDI = FLEN + TLEN - 1

c              WRITE(*, '(A)') T_FMT_STR
c              WRITE(*, '(3(I3, 2H, ))') TLEN, FLEN, ENDI
              WRITE(FMT_STR(FLEN:ENDI), '(A)') T_FMT_STR(1:TLEN)
              
          ENDDO

c         Proudly output the format string
c              WRITE(*, '(A)') FMT_STR

c         Format string for general headings:
c           added flexibility by num_ovars
          WRITE(WIDTH_STR, '(I3)') VAR_WIDTH-1
          WRITE(NUM_OVARS_STR, '(I3)') NUM_OVARS - 2
          G_FMT_STR = '(A5,1X,A3,1X,'
     &                 // TRIM(NUM_OVARS_STR) // '(1H|, A' 
     &                 // TRIM(WIDTH_STR) // '))'


c         Loop through each row of GROHEAD, outputting
c         each column value per row
          DO I=1, 4
              DO J=1,NUM_OVARS    
                  GRO_OUT(J) = GROHEAD(I, J)  
              ENDDO
               IF (I .LT. 4) THEN
c                 Write any old column headings
c                  WRITE (NOUTDG,'(A5,1X,A3,1X,4(A5,1X), 44(1X,A12))') 
c     &                       GRO_OUT(1:NUM_OVARS)
c                  WRITE (NOUTDG,'(A5,1X,A3,1X,50(1H|, A11))') 
c     &                       GRO_OUT(1:NUM_OVARS)
                  WRITE (NOUTDG,FMT=G_FMT_STR) 
     &                       GRO_OUT(1:NUM_OVARS)
               ELSE
c                   Write right-aligned column headings
                    WRITE (NOUTDG,FMT=FMT_STR) 
     &                       GRO_OUT(1:NUM_OVARS)
               ENDIF 
          ENDDO


c         Output format string for daily values
!          D_FMT_STR = '(1X, I4, 1X, I3, 4(1X, I' 
!     &                // TRIM(WIDTH_STR) // 
!     &                '), 95(1X, F' // 
!     &                TRIM(WIDTH_STR) // '.3))'

c          D_FMT_STR = '(1X, I4, 1X, I3, 2(9X, I3), 2(7X, F5.2))'
           D_FMT_STR = '(1X, I4, 1X, I3, 2(18X, I3), 95(1X, F' // 
     &                TRIM(WIDTH_STR) // '.3))'

c          WRITE(*, '(A)') D_FMT_STR

c         ::::::::::::::::::::::::::::::::::::::::::::

 
c     Initialise some date variables:
c     :::::::::::::::::::::::::::::::
      FROP   = CONTROL%FROP
!      YRSIM  = CONTROL%YRSIM


c     END of SEASINIT
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = RATE
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.RATE) THEN
c     END of RATE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = INTEGRATE
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.INTEGR) THEN

c     END of INTEGRATE
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     -----------------------------------------------------
c
c              DYNAMIC = OUTPUT
c
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSEIF(DYNAMIC.EQ.OUTPUT) THEN
      
          IF (IDETG .NE. 'Y') RETURN

c     Print daily output:
c     :::::::::::::::::::
      YRDOY  = CONTROL%YRDOY

      IF (YRDOY .GE. YRPLT) THEN
          DAP = MAX(0, TIMDIF(YRPLT, YRDOY))
!         DAS = MAX(0, TIMDIF(YRSIM, YRDOY))
          DAS = CONTROL % DAS

          !-------------------------------------------------------------
          !  Write output based on user specified frequency 
          !-------------------------------------------------------------
          IF ((MOD(DAS,FROP) .EQ. 0)      !Daily output every FROP days,
     &      .OR. (YRDOY .EQ. YRPLT)) THEN !on planting date

            CALL YR_DOY(YRDOY, YEAR, DOY)

            TMIN     = WEATHER%TMIN
            TMAX     = WEATHER%TMAX
            T_SD     = CaneCrop%LT
            LAID     = Growth%LAI
            LGDMD    = Part%TOPDM
            SMDMD    = Part%STKDM
            SMFMD    = Part%STKWM
            SUCMD    = Part%SUCMAS  
            RWAD     = Out%ROOTDM * 1000.
            RDMD     = Out%ROOTDM
            BADMD    = Part%AERLDM
            STKPOPm2 = CaneCrop%TOTPOP / 10000.
            L_SD     = CaneCrop%TMEANLF - CaneCrop%TMEANDEDLF
c           Harvest index = product / biomass
            IF (Part%STKDM .GT. 0.) THEN
               SUID     = Part%SUCMAS / Part%STKDM
            ELSE 
               SUID     = 0.
            ENDIF
            LDDMD    = Out%TRASDM
            LAITD    = Growth%TLAI
            WSPD     = WaterBal%SWDF1
            WSGD     = WaterBal%SWDF2
            
c           Calculate sucrose %, on a dry and fresh mass basis
c           STKWM is a function of STKDM, so it is only
c           necessary to check one.
            IF (Part%STKDM .GT. 0.00001) THEN
              SUDMD    = Part%SUCMAS / Part%STKDM * 100.
              SUFMD    = Part%SUCMAS / Part%STKWM * 100.
            ELSE
              SUDMD = 0.
              SUFMD = 0.  
            ENDIF

c           Specific leaf area (cm2/gram)
c           [calculated from TOPS, i.e. including meristem]
c           [assuming LI is an indicator of field coverage]
c                     ((area (m2) * LI * number of leaves/area) / total leaf DM)
            SLAD     = 0.
            IF (LGDMD .GT. 0.) SLAD = ((10000. * Growth%LI * 
     -                                Growth%LAI) / LGDMD)
            CHTD     = CaneCrop%SHGT ! / 100.
            
            
            RDPD     = WaterBal%RTDEP / 100.
            EOSA     = WaterBal%EOS
            
            LIPD     = Growth%LI * 100.
c           Total root length density
            RTLD = DOT_PRODUCT(SoilProp%DLAYR(1:SoilProp%NLAYR),
     -                         WaterBal%RLV(1:SoilProp%NLAYR))

            SWDF30 = WaterBal%SWDF30
            GROSSP = Out%GROSSP
            BRDMD  = Out%DWDT * 1000.
            PARCE  = Part%PARCE

c         Thermal time / heat units
            TTEBC= Out%CHU_EM
            TTSPC= Out%CHUPOP
            TTLEC= Out%CHUPI

 
c         MJ, Feb 2008:
c         Previously:
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c            WRITE(NOUTDG,400)YEAR, DOY,  DAS, DAP, T_SD, 0,
c     -                       LAID, LGDMD, SMFMD, SMDMD, SUCMD, RWAD, 
c     -                       BADMD, STKPOPm2, L_SD, SUID,
c     -                       NINT(LDDMD), LAITD, WSPD, WSGD, NSTD, LN_D,
c     -                       SH_D, SLAD, CHTD, CWID, NWAD, RDPD,
c     -                       WaterBal%RLV(1:10), EOSA, EOPA, TWUP, RTLD,
c     -                       LIPD, SWDF30, GROSSP, BRDMD, PARCE, Out%PAR,
c     -                       Part%FLODGE


c 400        FORMAT (1X, I4,1X, I3,2(1X, I5),1X, F5.0,1X, I5,1X, F5.3,1X,
c     -              5(F6.0,1X), F6.0,1X, F5.1, 1X, F5.2, 1X, F5.3, 
c     -              1X, I5,1X, F5.2,1X, F5.3, 9(1X, F5.3), 10(1X, F5.3),
c     -              1X, F5.2, 1X, F5.2, 1X, F5.2, 1X, F5.2, 1X, F5.2, 
c     -              1X, F5.2, 1X, F6.4, 1X, F6.4, 1X, F6.4, 1X, F6.2,
c     -              1X, F6.2)
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c         NOW:
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!       chp 9/26/2008 changed stress output to (1.0 - stress)

          WRITE(NOUTDG,FMT=D_FMT_STR) 
     &                     YEAR, DOY,  DAS, DAP, 
     & TAVE, ADM, CANma, CTT_LFEM, dADMSource, 
     & dGLA, dGLA_pot, dGLAI, dLeafDM, 
     & dlt_slai_light, dRootDM, 
     & dSDM, dSenDM, dSFibDM, dTotalSourceShortfall, 
     & CarbReserveDM, SourceShortfall, 
     & LeafSourceShortfall, StalkSourceShortfall, 
     & AdditionalSource, dAdditionalSource, 
     & dSSuc, dTTLf, FIinter, FT_LAI, 
     & FT_photos, GLAI, GLeafDM, 
     & SettDM, KePAR, LAIsen, LeafDM, 
     & LeafPI, LeafSource, LeafSink, 
     & LFNUM_OSG, LFNUM_SUCStart, 
     & NumLF, dNumLF, PIAvgSSRv2, 
     & RootDM, RootFrac, RGR_LAI_max, 
     & RGR_LAI_FTSW, RGR_LAI_act, 
     & SDM, SenDM, SXR, SFibDM, 
     & slai_light_fac, SLA, SLA_avg, 
     & Source, StalkSource, OSGfrac, 
     & SPF_act, STKFRAC_avg, StalkSink, 
     & StalkLength, SWDF2, SUCDM, 
     & TotalDM, vADM, vSource, 
     & ADMSource
!     &                     TMIN, TMAX, 
!     -                     NINT(T_SD), 
!     -                     CaneCrop%GROPHASE,
!     -                     LAID, LGDMD, SMFMD, SMDMD, SUCMD, RDMD, 
!     -                     BADMD, STKPOPm2, L_SD,   SUID,
!     -                     LDDMD, LAITD, 1.0-WSPD, 1.0-WSGD, 1.0-SWDF30,
!     -                     CHTD, RDPD,
!     -                     WaterBal%RLV(1:10), EOSA, RTLD,
!     -                     LIPD, SLAD, GROSSP, BRDMD, 
!     -                     Part%FLODGE, TTEBC, TTSPC, TTLEC, SUDMD,
!     &              SUFMD, Canecrop%CWSI, PARCE, Part%RESP_G,
!     &              Part%RESP_M, CELLSE_DM

c         Format statement
c         Original one
c         Dates, etc
c  500     FORMAT(I5,1X,I3,4(1X, I5),1X,
c         Up to 97 (total 100 vars) 12-character floating point nums
c     &           95(F12.3, 1X))

c         Latest one:
c  500     FORMAT(1X, I4, 1X, I3, 4(1X, I11),
c         Up to 95 (total 100 vars) 11-character floating point nums
c     &           95(1X, F11.3))


c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          ENDIF
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
          IF (IDETG .NE. 'Y') RETURN

c     Close the output file
          CLOSE(UNIT=NOUTDG)
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
c     END of SUBROUTINE SC_RGOUTPUT()
c     -----------------------------------------------------