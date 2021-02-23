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
     & YRPLT, CELLSE_DM, TAVE,
     & ADM,   ! Above-ground dry mass t/ha
     &  CANma,   ! Conversion from leaf to dry mass
     &  CTT_LFEM,   ! Cumulative thermal time for leaf appearance °Cd
     &  dADM,   ! Daily change in Above-ground dry mass
     &  dGLA,   ! New daily green leaf area m2/m2
     &  dGLA_pot,   ! Potential (source-unlimited) change in green leaf area m2/m2
     &  dGLAI,   ! Daily change in green leaf area index, new growth - senescence
     &  dLeafDM,   ! Daily change in total leaf dry mass t/ha/d
     &  dLeafPI,   ! Daily change in leaf phyllocron interval l/s/d
     &  dlt_slai_light,   ! Daily leaf area senescence m2/m2/d
     &  dRootDM,   ! Daily change in Root dry mass
     &  dSDM,   ! Daily change in stalk dry mass t/ha/d
     &  dSenDM,   ! Daily leaf area senescence m2/m2/d
     &  dSFibDM,   ! Daily change in stalk fibre dry mass t/ha/d
     &  dSSuc,   ! Daily change in stalk sucrose t/ha/d
     &  dTTLf,   ! Daily thermal time accumulation for driving leaf appearance °Cd
     &  FIinter,   ! Inter-row PAR fractional interception
     &  FT_LAI,   ! Temperature factor for LAI growth
     &  FT_photos,   ! Daily temperature factor for photosynthesis
     &  GLAI,   ! Green leaf area index m2/m2
     &  GLeafDM,   ! Green leaf canopy dry mass t/ha
     &  KePAR,   ! PAR canopy extinction coefficient
     &  LAIsen,   ! Daily GLAI senesced m2/m2/d
     &  LeafDM,   ! Total (living + dead) leaf dry mass t/ha
     &  LeafPI,   ! Current leaf phyllocron interval °Cd
     &  LeafSink,   ! Leaf carbon demand (sink strength) t/ha/d
     &  LFNUM_OSG,   ! The number of leaves per shoot at which onset of stalk growth started l/s
     &  LFNUM_SUCStart,   ! Leaf number per stalk at which sucrose accumulation can start. l/s
     &  NumLF,   ! Number of leaves per shoot (reference) l/s
     &  RootFrac,   ! Daily fraction of biomass allocated to roots, g/g
     &  RootDM,   ! Root dry mass
     &  RGR_LAI_max,   ! Today's maximum relative LAI growth rate, considering temperature
     &  SDM,   ! Stalk dry mass t/ha
     &  SenDM,   ! Dry mass of senesced leaves t/ha
     &  SER,   ! Stalk elongation rate cm/d
     &  slai_light_fac,   ! Daily fraction of GLAI senesced
     &  SLSR,   ! Source:sink ratio
     &  SLA,   ! specific leaf area (cm2/g)
     &  Source,   ! Daily biomass increase (source strength) t/ha/d
     &  SPF,   ! Stalk partitioning fraction t/t
     &  StalkSink,   ! Stalk fibre growth sink strength t/ha/d
     &  StalkLength,   ! Stalk length (height to top visible dewlap) cm
     &  SFibDM,   ! Stalk fibre dry mass t/ha
     &  SWDF2,   ! Soil water deficit (stress) factor affecting expansive growth
     &  SUCDM,   ! stalk sucrose (t/ha)
     &  TotalDM,   ! Total crop dry mass t/ha
     &  vADM,   ! Verify ADM t/ha
     &  vSource)   ! Verify Source t/ha

     

      
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
      REAL, INTENT(IN) :: TAVE  ! daily mean temperature
      REAL, INTENT(IN) :: ADM  ! Above-ground dry mass t/ha
      REAL, INTENT(IN) :: CANma  ! Conversion from leaf to dry mass
      REAL, INTENT(IN) :: CTT_LFEM  ! Cumulative thermal time for leaf appearance °Cd
      REAL, INTENT(IN) :: dADM  ! Daily change in Above-ground dry mass
      REAL, INTENT(IN) :: dGLA  ! New daily green leaf area m2/m2
      REAL, INTENT(IN) :: dGLA_pot  ! Potential (source-unlimited) change in green leaf area m2/m2
      REAL, INTENT(IN) :: dGLAI  ! Daily change in green leaf area index, new growth - senescence
      REAL, INTENT(IN) :: dLeafDM  ! Daily change in total leaf dry mass t/ha/d
      REAL, INTENT(IN) :: dLeafPI  ! Daily change in leaf phyllocron interval l/s/d
      REAL, INTENT(IN) :: dlt_slai_light  ! Daily leaf area senescence m2/m2/d
      REAL, INTENT(IN) :: dRootDM  ! Daily change in Root dry mass
      REAL, INTENT(IN) :: dSDM  ! Daily change in stalk dry mass t/ha/d
      REAL, INTENT(IN) :: dSenDM  ! Daily leaf area senescence m2/m2/d
      REAL, INTENT(IN) :: dSFibDM  ! Daily change in stalk fibre dry mass t/ha/d
      REAL, INTENT(IN) :: dSSuc  ! Daily change in stalk sucrose t/ha/d
      REAL, INTENT(IN) :: dTTLf  ! Daily thermal time accumulation for driving leaf appearance °Cd
      REAL, INTENT(IN) :: FIinter  ! Inter-row PAR fractional interception
      REAL, INTENT(IN) :: FT_LAI  ! Temperature factor for LAI growth
      REAL, INTENT(IN) :: FT_photos  ! Daily temperature factor for photosynthesis
      REAL, INTENT(IN) :: GLAI  ! Green leaf area index m2/m2
      REAL, INTENT(IN) :: GLeafDM  ! Green leaf canopy dry mass t/ha
      REAL, INTENT(IN) :: KePAR  ! PAR canopy extinction coefficient
      REAL, INTENT(IN) :: LAIsen  ! Daily GLAI senesced m2/m2/d
      REAL, INTENT(IN) :: LeafDM  ! Total (living + dead) leaf dry mass t/ha
      REAL, INTENT(IN) :: LeafPI  ! Current leaf phyllocron interval °Cd
      REAL, INTENT(IN) :: LeafSink  ! Leaf carbon demand (sink strength) t/ha/d
      REAL, INTENT(IN) :: LFNUM_OSG  ! The number of leaves per shoot at which onset of stalk growth started l/s
      REAL, INTENT(IN) :: LFNUM_SUCStart  ! Leaf number per stalk at which sucrose accumulation can start. l/s
      REAL, INTENT(IN) :: NumLF  ! Number of leaves per shoot (reference) l/s
      REAL, INTENT(IN) :: RootFrac  ! Daily fraction of biomass allocated to roots, g/g
      REAL, INTENT(IN) :: RootDM  ! Root dry mass
      REAL, INTENT(IN) :: RGR_LAI_max  ! Today's maximum relative LAI growth rate, considering temperature
      REAL, INTENT(IN) :: SDM  ! Stalk dry mass t/ha
      REAL, INTENT(IN) :: SenDM  ! Dry mass of senesced leaves t/ha
      REAL, INTENT(IN) :: SER  ! Stalk elongation rate cm/d
      REAL, INTENT(IN) :: slai_light_fac  ! Daily fraction of GLAI senesced
      REAL, INTENT(IN) :: SLSR  ! Source:sink ratio
      REAL, INTENT(IN) :: SLA  ! specific leaf area (cm2/g)
      REAL, INTENT(IN) :: Source  ! Daily biomass increase (source strength) t/ha/d
      REAL, INTENT(IN) :: SPF  ! Stalk partitioning fraction t/t
      REAL, INTENT(IN) :: StalkSink  ! Stalk fibre growth sink strength t/ha/d
      REAL, INTENT(IN) :: StalkLength  ! Stalk length (height to top visible dewlap) cm
      REAL, INTENT(IN) :: SFibDM  ! Stalk fibre dry mass t/ha
      REAL, INTENT(IN) :: SWDF2  ! Soil water deficit (stress) factor affecting expansive growth
      REAL, INTENT(IN) :: SUCDM  ! stalk sucrose (t/ha)
      REAL, INTENT(IN) :: TotalDM  ! Total crop dry mass t/ha
      REAL, INTENT(IN) :: vADM  ! Verify ADM t/ha
      REAL, INTENT(IN) :: vSource  ! Verify Source t/ha




c     CANEGRO 3.5 variables:
c     ::::::::::::::::::::::
c     GROHEAD - this is an array variable that stores the 
c     text headings in the plantgro.out file.
c     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     MJ, Feb 2008: Grohead was previously an array of four long
c     strings.  However, this was a pain to maintain, so it is
c     now a 2-dimensional array with each column header stored
c     separately.
      CHARACTER*15 GROHEAD(4, 53)
      CHARACTER*15 GRO_OUT(53)


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
 
c        Above-ground dry mass t/ha
         DATA GROHEAD(1, 6) /'ABOVE-GROUN'/
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
         DATA GROHEAD(4, 9) /'dADM'/
 
c        New daily green leaf area m2/m2
         DATA GROHEAD(1, 10) /'NEW'/
         DATA GROHEAD(2, 10) /'DAILY'/
         DATA GROHEAD(3, 10) /'GREEN'/
         DATA GROHEAD(4, 10) /'dGLA'/
 
c        Potential (source-unlimited) change in green leaf area m2/m2
         DATA GROHEAD(1, 11) /'POTENTIAL'/
         DATA GROHEAD(2, 11) /'(SOURCE-UNL'/
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
 
c        Daily change in leaf phyllocron interval l/s/d
         DATA GROHEAD(1, 14) /'DAILY'/
         DATA GROHEAD(2, 14) /'CHANGE'/
         DATA GROHEAD(3, 14) /'IN'/
         DATA GROHEAD(4, 14) /'dLeafPI'/
 
c        Daily leaf area senescence m2/m2/d
         DATA GROHEAD(1, 15) /'DAILY'/
         DATA GROHEAD(2, 15) /'LEAF'/
         DATA GROHEAD(3, 15) /'AREA'/
         DATA GROHEAD(4, 15) /'dlt_slai_li'/
 
c        Daily change in Root dry mass
         DATA GROHEAD(1, 16) /'DAILY'/
         DATA GROHEAD(2, 16) /'CHANGE'/
         DATA GROHEAD(3, 16) /'IN'/
         DATA GROHEAD(4, 16) /'dRootDM'/
 
c        Daily change in stalk dry mass t/ha/d
         DATA GROHEAD(1, 17) /'DAILY'/
         DATA GROHEAD(2, 17) /'CHANGE'/
         DATA GROHEAD(3, 17) /'IN'/
         DATA GROHEAD(4, 17) /'dSDM'/
 
c        Daily leaf area senescence m2/m2/d
         DATA GROHEAD(1, 18) /'DAILY'/
         DATA GROHEAD(2, 18) /'LEAF'/
         DATA GROHEAD(3, 18) /'AREA'/
         DATA GROHEAD(4, 18) /'dSenDM'/
 
c        Daily change in stalk fibre dry mass t/ha/d
         DATA GROHEAD(1, 19) /'DAILY'/
         DATA GROHEAD(2, 19) /'CHANGE'/
         DATA GROHEAD(3, 19) /'IN'/
         DATA GROHEAD(4, 19) /'dSFibDM'/
 
c        Daily change in stalk sucrose t/ha/d
         DATA GROHEAD(1, 20) /'DAILY'/
         DATA GROHEAD(2, 20) /'CHANGE'/
         DATA GROHEAD(3, 20) /'IN'/
         DATA GROHEAD(4, 20) /'dSSuc'/
 
c        Daily thermal time accumulation for driving leaf appearance °Cd
         DATA GROHEAD(1, 21) /'DAILY'/
         DATA GROHEAD(2, 21) /'THERMAL'/
         DATA GROHEAD(3, 21) /'TIME'/
         DATA GROHEAD(4, 21) /'dTTLf'/
 
c        Inter-row PAR fractional interception
         DATA GROHEAD(1, 22) /'INTER-ROW'/
         DATA GROHEAD(2, 22) /'PAR'/
         DATA GROHEAD(3, 22) /'FRACTIONAL'/
         DATA GROHEAD(4, 22) /'LI%D'/
 
c        Temperature factor for LAI growth
         DATA GROHEAD(1, 23) /'TEMPERATURE'/
         DATA GROHEAD(2, 23) /'FACTOR'/
         DATA GROHEAD(3, 23) /'FOR'/
         DATA GROHEAD(4, 23) /'FT_LAI'/
 
c        Daily temperature factor for photosynthesis
         DATA GROHEAD(1, 24) /'DAILY'/
         DATA GROHEAD(2, 24) /'TEMPERATURE'/
         DATA GROHEAD(3, 24) /'FACTOR'/
         DATA GROHEAD(4, 24) /'FT_photos'/
 
c        Green leaf area index m2/m2
         DATA GROHEAD(1, 25) /'GREEN'/
         DATA GROHEAD(2, 25) /'LEAF'/
         DATA GROHEAD(3, 25) /'AREA'/
         DATA GROHEAD(4, 25) /'LAIGD'/
 
c        Green leaf canopy dry mass t/ha
         DATA GROHEAD(1, 26) /'GREEN'/
         DATA GROHEAD(2, 26) /'LEAF'/
         DATA GROHEAD(3, 26) /'CANOPY'/
         DATA GROHEAD(4, 26) /'LGDMD'/
 
c        PAR canopy extinction coefficient
         DATA GROHEAD(1, 27) /'PAR'/
         DATA GROHEAD(2, 27) /'CANOPY'/
         DATA GROHEAD(3, 27) /'EXTINCTION'/
         DATA GROHEAD(4, 27) /'KePAR'/
 
c        Daily GLAI senesced m2/m2/d
         DATA GROHEAD(1, 28) /'DAILY'/
         DATA GROHEAD(2, 28) /'GLAI'/
         DATA GROHEAD(3, 28) /'SENESCED'/
         DATA GROHEAD(4, 28) /'LAIsen'/
 
c        Total (living + dead) leaf dry mass t/ha
         DATA GROHEAD(1, 29) /'TOTAL'/
         DATA GROHEAD(2, 29) /'(LIVING'/
         DATA GROHEAD(3, 29) /'+'/
         DATA GROHEAD(4, 29) /'LTDMD'/
 
c        Current leaf phyllocron interval °Cd
         DATA GROHEAD(1, 30) /'CURRENT'/
         DATA GROHEAD(2, 30) /'LEAF'/
         DATA GROHEAD(3, 30) /'PHYLLOCRON'/
         DATA GROHEAD(4, 30) /'LeafPI'/
 
c        Leaf carbon demand (sink strength) t/ha/d
         DATA GROHEAD(1, 31) /'LEAF'/
         DATA GROHEAD(2, 31) /'CARBON'/
         DATA GROHEAD(3, 31) /'DEMAND'/
         DATA GROHEAD(4, 31) /'LeafSink'/
 
c        The number of leaves per shoot at which onset of stalk growth started l/s
         DATA GROHEAD(1, 32) /'THE'/
         DATA GROHEAD(2, 32) /'NUMBER'/
         DATA GROHEAD(3, 32) /'OF'/
         DATA GROHEAD(4, 32) /'LFNUM_OSG'/
 
c        Leaf number per stalk at which sucrose accumulation can start. l/s
         DATA GROHEAD(1, 33) /'LEAF'/
         DATA GROHEAD(2, 33) /'NUMBER'/
         DATA GROHEAD(3, 33) /'PER'/
         DATA GROHEAD(4, 33) /'LFNUM_SUCSt'/
 
c        Number of leaves per shoot (reference) l/s
         DATA GROHEAD(1, 34) /'NUMBER'/
         DATA GROHEAD(2, 34) /'OF'/
         DATA GROHEAD(3, 34) /'LEAVES'/
         DATA GROHEAD(4, 34) /'NumLF'/
 
c        Daily fraction of biomass allocated to roots, g/g
         DATA GROHEAD(1, 35) /'DAILY'/
         DATA GROHEAD(2, 35) /'FRACTION'/
         DATA GROHEAD(3, 35) /'OF'/
         DATA GROHEAD(4, 35) /'RootFrac'/
 
c        Root dry mass
         DATA GROHEAD(1, 36) /'ROOT'/
         DATA GROHEAD(2, 36) /'DRY'/
         DATA GROHEAD(3, 36) /'MASS'/
         DATA GROHEAD(4, 36) /'RootDM'/
 
c        Today_s maximum relative LAI growth rate, considering temperature
         DATA GROHEAD(1, 37) /'TODAY_S'/
         DATA GROHEAD(2, 37) /'MAXIMUM'/
         DATA GROHEAD(3, 37) /'RELATIVE'/
         DATA GROHEAD(4, 37) /'RGR_LAI_max'/
 
c        Stalk dry mass t/ha
         DATA GROHEAD(1, 38) /'STALK'/
         DATA GROHEAD(2, 38) /'DRY'/
         DATA GROHEAD(3, 38) /'MASS'/
         DATA GROHEAD(4, 38) /'SMDMD'/
 
c        Dry mass of senesced leaves t/ha
         DATA GROHEAD(1, 39) /'DRY'/
         DATA GROHEAD(2, 39) /'MASS'/
         DATA GROHEAD(3, 39) /'OF'/
         DATA GROHEAD(4, 39) /'LDDMD'/
 
c        Stalk elongation rate cm/d
         DATA GROHEAD(1, 40) /'STALK'/
         DATA GROHEAD(2, 40) /'ELONGATION'/
         DATA GROHEAD(3, 40) /'RATE'/
         DATA GROHEAD(4, 40) /'SER'/
 
c        Daily fraction of GLAI senesced
         DATA GROHEAD(1, 41) /'DAILY'/
         DATA GROHEAD(2, 41) /'FRACTION'/
         DATA GROHEAD(3, 41) /'OF'/
         DATA GROHEAD(4, 41) /'slai_light_'/
 
c        Source:sink ratio
         DATA GROHEAD(1, 42) /'SOURCE:SINK'/
         DATA GROHEAD(2, 42) /'RATIO'/
         DATA GROHEAD(3, 42) /'NA'/
         DATA GROHEAD(4, 42) /'SLSR'/
 
c        specific leaf area (cm2/g)
         DATA GROHEAD(1, 43) /'SPECIFIC'/
         DATA GROHEAD(2, 43) /'LEAF'/
         DATA GROHEAD(3, 43) /'AREA'/
         DATA GROHEAD(4, 43) /'SLAD'/
 
c        Daily biomass increase (source strength) t/ha/d
         DATA GROHEAD(1, 44) /'DAILY'/
         DATA GROHEAD(2, 44) /'BIOMASS'/
         DATA GROHEAD(3, 44) /'INCREASE'/
         DATA GROHEAD(4, 44) /'Source'/
 
c        Stalk partitioning fraction t/t
         DATA GROHEAD(1, 45) /'STALK'/
         DATA GROHEAD(2, 45) /'PARTITIONIN'/
         DATA GROHEAD(3, 45) /'FRACTION'/
         DATA GROHEAD(4, 45) /'SPF'/
 
c        Stalk fibre growth sink strength t/ha/d
         DATA GROHEAD(1, 46) /'STALK'/
         DATA GROHEAD(2, 46) /'FIBRE'/
         DATA GROHEAD(3, 46) /'GROWTH'/
         DATA GROHEAD(4, 46) /'StalkSink'/
 
c        Stalk length (height to top visible dewlap) cm
         DATA GROHEAD(1, 47) /'STALK'/
         DATA GROHEAD(2, 47) /'LENGTH'/
         DATA GROHEAD(3, 47) /'(HEIGHT'/
         DATA GROHEAD(4, 47) /'StalkLength'/
 
c        Stalk fibre dry mass t/ha
         DATA GROHEAD(1, 48) /'STALK'/
         DATA GROHEAD(2, 48) /'FIBRE'/
         DATA GROHEAD(3, 48) /'DRY'/
         DATA GROHEAD(4, 48) /'SFibDM'/
 
c        Soil water deficit (stress) factor affecting expansive growth
         DATA GROHEAD(1, 49) /'SOIL'/
         DATA GROHEAD(2, 49) /'WATER'/
         DATA GROHEAD(3, 49) /'DEFICIT'/
         DATA GROHEAD(4, 49) /'WSGD'/
 
c        stalk sucrose (t/ha)
         DATA GROHEAD(1, 50) /'STALK'/
         DATA GROHEAD(2, 50) /'SUCROSE'/
         DATA GROHEAD(3, 50) /'(T/HA)'/
         DATA GROHEAD(4, 50) /'SUCMD'/
 
c        Total crop dry mass t/ha
         DATA GROHEAD(1, 51) /'TOTAL'/
         DATA GROHEAD(2, 51) /'CROP'/
         DATA GROHEAD(3, 51) /'DRY'/
         DATA GROHEAD(4, 51) /'TotalDM'/
 
c        Verify ADM t/ha
         DATA GROHEAD(1, 52) /'VERIFY'/
         DATA GROHEAD(2, 52) /'ADM'/
         DATA GROHEAD(3, 52) /'T/HA'/
         DATA GROHEAD(4, 52) /'vADM'/
 
c        Verify Source t/ha
         DATA GROHEAD(1, 53) /'VERIFY'/
         DATA GROHEAD(2, 53) /'SOURCE'/
         DATA GROHEAD(3, 53) /'T/HA'/
         DATA GROHEAD(4, 53) /'vSource'/


c          ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c          Number of output variables (important for output!)
           DATA NUM_OVARS /53/

c          Width of output columns:
           DATA VAR_WIDTH /12/
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
           D_FMT_STR = '(1X, I4, 1X, I3, 2(9X, I3), 95(1X, F' // 
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
     & TAVE, ADM,
     & CANma, CTT_LFEM,
     & dADM, dGLA, dGLA_pot, dGLAI,
     & dLeafDM, dLeafPI, dlt_slai_light, 
     & dRootDM, dSDM, dSenDM, dSFibDM, 
     & dSSuc, dTTLf, FIinter, FT_LAI, 
     & FT_photos, GLAI, GLeafDM,
     & KePAR, LAIsen, LeafDM, LeafPI,
     & LeafSink, LFNUM_OSG, LFNUM_SUCStart,
     & NumLF, RootFrac, RootDM,
     & RGR_LAI_max, SDM, SenDM, 
     & SER, slai_light_fac, SLSR,
     & SLA, Source, SPF, StalkSink, 
     & StalkLength, SFibDM, SWDF2, 
     & SUCDM, TotalDM, vADM, vSource
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