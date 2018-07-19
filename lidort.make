#------------------------------------------- LIDORT master modules

LIDORT_V23E_MASTER.o: $(MPATH_E)LIDORT_V23E_MASTER.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_RESULTS.VARS        \
         $(IPATH_E)LIDORT_L_CONTROL.VARS      \
         $(IPATH_E)LIDORT_L_RESULTS.VARS
	$(FF) $(MPATH_E)LIDORT_V23E_MASTER.f

LIDORT_V23E_INPUT.o: $(MPATH_E)LIDORT_V23E_INPUT.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS
	$(FF) $(MPATH_E)LIDORT_V23E_INPUT.f

LIDORT_V23E_FOURIER.o: $(MPATH_E)LIDORT_V23E_FOURIER.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_E)LIDORT_L_CONTROL.VARS      \
         $(IPATH_E)LIDORT_L_GEOPHYS.VARS      \
         $(IPATH_E)LIDORT_L_SOLUTION.VARS
	$(FF) $(MPATH_E)LIDORT_V23E_FOURIER.f

###################
#  Standard code  #
###################

#------------------------------------------- LIDORT Setups modules

LIDORT_MISCSETUPS.o: $(SPATH_S)LIDORT_MISCSETUPS.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS
	$(FF) $(SPATH_S)LIDORT_MISCSETUPS.f

LIDORT_CHECKINPUT.o: $(SPATH_S)LIDORT_CHECKINPUT.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS
	$(FF) $(SPATH_S)LIDORT_CHECKINPUT.f

LIDORT_DERIVEINPUT.o: $(SPATH_S)LIDORT_DERIVEINPUT.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS
	$(FF) $(SPATH_S)LIDORT_DERIVEINPUT.f

#------------------------------------------- LIDORT solutions modules

LIDORT_RTSOLUTION.o: $(SPATH_S)LIDORT_RTSOLUTION.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_LEGENDRE.VARS       \
         $(IPATH_S)LIDORT_SOLUTION.VARS
	$(FF) $(SPATH_S)LIDORT_RTSOLUTION.f

LIDORT_MULTIPLIERS.o: $(SPATH_S)LIDORT_MULTIPLIERS.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_LEGENDRE.VARS       \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS
	$(FF) $(SPATH_S)LIDORT_MULTIPLIERS.f

#------------------------------------------- LIDORT Intensity modules

LIDORT_INTENSITY.o: $(SPATH_S)LIDORT_INTENSITY.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(FF) $(SPATH_S)LIDORT_INTENSITY.f

LIDORT_CONVERGE.o: $(SPATH_S)LIDORT_CONVERGE.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(FF) $(SPATH_S)LIDORT_CONVERGE.f

LIDORT_BVPROBLEM.o: $(SPATH_S)LIDORT_BVPROBLEM.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS
	$(FF) $(SPATH_S)LIDORT_BVPROBLEM.f

LIDORT_SSCORRECTION.o: $(SPATH_S)LIDORT_SSCORRECTION.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SSCORRECTION.VARS   \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(FF) $(SPATH_S)LIDORT_SSCORRECTION.f
#------------------------------------------- LIDORT Auxiliary modules

LIDORT_AUX.o: $(SPATH_S)LIDORT_AUX.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS
	$(FF) $(SPATH_S)LIDORT_AUX.f

LIDORT_READINPUT.o: $(SPATH_S)LIDORT_READINPUT.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS
	$(FF) $(SPATH_S)LIDORT_READINPUT.f

LIDORT_WRITEMODULES.o: $(SPATH_S)LIDORT_WRITEMODULES.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_RESULTS.VARS
	$(FF) $(SPATH_S)LIDORT_WRITEMODULES.f

###################
#  Extended code  #
###################

#------------------------------------------- LIDORT Setups modules

LIDORT_L_MISCSETUPS.o: $(SPATH_E)LIDORT_L_MISCSETUPS.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_E)LIDORT_L_CONTROL.VARS      \
         $(IPATH_E)LIDORT_L_GEOPHYS.VARS      \
         $(IPATH_E)LIDORT_L_SETUPS.VARS       \
         $(IPATH_E)LIDORT_L_MULTIPLIERS.VARS
	$(FF) $(SPATH_E)LIDORT_L_MISCSETUPS.f

LIDORT_L_CHECKINPUT.o: $(SPATH_E)LIDORT_L_CHECKINPUT.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_E)LIDORT_L_CONTROL.VARS
	$(FF) $(SPATH_E)LIDORT_L_CHECKINPUT.f

#------------------------------------------- LIDORT linearized solution modules

LIDORT_L_RTSOLUTION.o: $(SPATH_E)LIDORT_L_RTSOLUTION.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_LEGENDRE.VARS       \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_E)LIDORT_L_SETUPS.VARS       \
         $(IPATH_E)LIDORT_L_SOLUTION.VARS
	$(FF) $(SPATH_E)LIDORT_L_RTSOLUTION.f

LIDORT_L_MULTIPLIERS.o: $(SPATH_E)LIDORT_L_MULTIPLIERS.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_E)LIDORT_L_SETUPS.VARS       \
         $(IPATH_E)LIDORT_L_MULTIPLIERS.VARS
	$(FF) $(SPATH_E)LIDORT_L_MULTIPLIERS.f

#------------------------------------------- LIDORT Weighting function modules

LIDORT_L_BVSETUPS.o: $(SPATH_E)LIDORT_L_BVSETUPS.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_E)LIDORT_L_SETUPS.VARS       \
         $(IPATH_E)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_E)LIDORT_L_MULTIPLIERS.VARS
	$(FF) $(SPATH_E)LIDORT_L_BVSETUPS.f

LIDORT_L_WFCALC.o: $(SPATH_E)LIDORT_L_WFCALC.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_E)LIDORT_L_SETUPS.VARS       \
         $(IPATH_E)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_E)LIDORT_L_RESULTS.VARS
	$(FF) $(SPATH_E)LIDORT_L_WFCALC.f

LIDORT_L_SSCORRECTION.o: $(SPATH_E)LIDORT_L_SSCORRECTION.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_S)LIDORT_SETUPS.VARS         \
         $(IPATH_S)LIDORT_SOLUTION.VARS       \
         $(IPATH_S)LIDORT_MULTIPLIERS.VARS    \
         $(IPATH_S)LIDORT_SSCORRECTION.VARS   \
         $(IPATH_E)LIDORT_L_CONTROL.VARS      \
         $(IPATH_E)LIDORT_L_GEOPHYS.VARS      \
         $(IPATH_E)LIDORT_L_SETUPS.VARS       \
         $(IPATH_E)LIDORT_L_SOLUTION.VARS     \
         $(IPATH_E)LIDORT_L_MULTIPLIERS.VARS  \
         $(IPATH_E)LIDORT_L_RESULTS.VARS
	$(FF) $(SPATH_E)LIDORT_L_SSCORRECTION.f

LIDORT_L_FOURIERADD.o: $(SPATH_E)LIDORT_L_FOURIERADD.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_E)LIDORT_L_CONTROL.VARS      \
         $(IPATH_E)LIDORT_L_GEOPHYS.VARS      \
         $(IPATH_E)LIDORT_L_RESULTS.VARS
	$(FF) $(SPATH_E)LIDORT_L_FOURIERADD.f

#------------------------------------------- LIDORT Auxiliary modules

LIDORT_L_INPUTREAD.o: $(SPATH_E)LIDORT_L_INPUTREAD.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L_CONTROL.VARS
	$(FF) $(SPATH_E)LIDORT_L_INPUTREAD.f

LIDORT_L_WRITEMODULES.o: $(SPATH_E)LIDORT_L_WRITEMODULES.f    \
         $(IPATH_S)LIDORT.PARS                \
         $(IPATH_E)LIDORT_L.PARS              \
         $(IPATH_S)LIDORT_CONTROL.VARS        \
         $(IPATH_S)LIDORT_GEOPHYS.VARS        \
         $(IPATH_S)LIDORT_MODEL.VARS          \
         $(IPATH_E)LIDORT_L_CONTROL.VARS      \
         $(IPATH_E)LIDORT_L_GEOPHYS.VARS      \
         $(IPATH_E)LIDORT_L_RESULTS.VARS
	$(FF) $(SPATH_E)LIDORT_L_WRITEMODULES.f
