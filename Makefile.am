bin_PROGRAMS = afid
afid_SOURCES = \
DeallocateVariables.F90       InitVariables.F90    ResetLogs.F90                \
AuxiliaryRoutines.F90         DebugRoutines.F90        InitPressureSolver.F90     SetTempBCs.F90               \
CalcDissipationNu.F90         decomp_2d.F90            interp.F90                 SlabDumpRoutines.F90         \
CalcMaxCFL.F90                decomp_2d_fft.F90        HdfRoutines.F90                LocateLargeDivergence.F90  SolveImpEqnUpdate_Temp.F90   \
CalcLocalDivergence.F90  ExplicitTermsTemp.F90    main.F90                   SolveImpEqnUpdate_X.F90     \
CalcPlateNu.F90          ExplicitTermsVX.F90      ImplicitAndUpdateTemp.F90      SolveImpEqnUpdate_YZ.F90      WriteFlowField.F90 \
CheckDivergence.F90           ExplicitTermsVY.F90      ImplicitAndUpdateVX.F90        MpiAuxRoutines.F90         SolvePressureCorrection.F90  WriteGridInfo.F90 \
CorrectPressure.F90           ExplicitTermsVZ.F90      ImplicitAndUpdateVY.F90        StatReadReduceWrite.F90 \
CorrectVelocity.F90           ImplicitAndUpdateVZ.F90        param.F90                  StatRoutines.F90 \
CreateGrid.F90                QuitRoutine.F90            TimeMarcher.F90 \
CreateInitialConditions.F90   GlobalQuantities.F90     InitTimeMarchScheme.F90  ReadInputFile.F90         ReadFlowField.F90 HdfReadContinua.F90

afid_LDADD = $(FFTW3_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LIBS)

EXTRA_DIST = bou.in transpose_z_to_x.F90 transpose_x_to_z.F90 transpose_x_to_y.F90 transpose_y_to_x.F90 transpose_y_to_z.F90 transpose_z_to_y.F90 factor.F90 halo.F90 fft_common.F90 alloc.F90 halo_common.F90

# Change preprocessor flag. This is specific to some Fortran compilers.
# This hack shouldn't be used if C/C++ code is included!
dist-hook:
	chmod u+w $(distdir)/configure
	sed -i 's/^DEFS=-DHAVE_CONFIG_H$$/DEFS=$${FC_DEFINE}HAVE_CONFIG_H/' $(distdir)/configure	

# Clean module files
mostlyclean-local:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

SolveImpEqnUpdate_X.o: SolveImpEqnUpdate_X.F90  param.o
CalcDissipationNu.o: CalcDissipationNu.F90  param.o
ExplicitTermsVZ.o: ExplicitTermsVZ.F90  param.o
TimeMarcher.o: TimeMarcher.F90  decomp_2d.o param.o
InitPressureSolver.o: InitPressureSolver.F90  decomp_2d_fft.o param.o
LocateLargeDivergence.o: LocateLargeDivergence.F90  param.o
CalcMaxCFL.o: CalcMaxCFL.F90  decomp_2d.o param.o
fft_common.o: fft_common.F90  
InitTimeMarchScheme.o: InitTimeMarchScheme.F90  param.o
AuxiliaryRoutines.o: AuxiliaryRoutines.F90  decomp_2d.o
ImplicitAndUpdateVX.o: ImplicitAndUpdateVX.F90  param.o
QuitRoutine.o: QuitRoutine.F90  decomp_2d_fft.o param.o
SetTempBCs.o: SetTempBCs.F90  param.o
ResetLogs.o: ResetLogs.F90  param.o
CorrectPressure.o: CorrectPressure.F90  param.o
SolvePressureCorrection.o: SolvePressureCorrection.F90  decomp_2d.o decomp_2d_fft.o param.o
param.o: param.F90
InitVariables.o: InitVariables.F90  AuxiliaryRoutines.o decomp_2d.o param.o
SolveImpEqnUpdate_YZ.o: SolveImpEqnUpdate_YZ.F90  param.o
ExplicitTermsVX.o: ExplicitTermsVX.F90  param.o
ExplicitTermsTemp.o: ExplicitTermsTemp.F90  param.o
DebugRoutines.o: DebugRoutines.F90  param.o
WriteGridInfo.o: WriteGridInfo.F90  param.o
CalcPlateNu.o: CalcPlateNu.F90  param.o
ImplicitAndUpdateVY.o: ImplicitAndUpdateVY.F90  param.o
CreateGrid.o: CreateGrid.F90  AuxiliaryRoutines.o param.o
CalcLocalDivergence.o: CalcLocalDivergence.F90  param.o
ReadFlowField.o: ReadFlowField.F90  decomp_2d.o param.o
CheckDivergence.o: CheckDivergence.F90  param.o
ReadInputFile.o: ReadInputFile.F90  param.o
GlobalQuantities.o: GlobalQuantities.F90  param.o
CreateInitialConditions.o: CreateInitialConditions.F90  param.o
decomp_2d_fft.o: decomp_2d_fft.F90  decomp_2d.o
ExplicitTermsVY.o: ExplicitTermsVY.F90  param.o
HdfRoutines.o: HdfRoutines.F90  param.o
SlabDumpRoutines.o: SlabDumpRoutines.F90  param.o
DeallocateVariables.o: DeallocateVariables.F90  AuxiliaryRoutines.o param.o
main.o: main.F90  decomp_2d.o decomp_2d_fft.o param.o
SolveImpEqnUpdate_Temp.o: SolveImpEqnUpdate_Temp.F90  param.o
WriteFlowField.o: WriteFlowField.F90  param.o
CorrectVelocity.o: CorrectVelocity.F90  param.o
MpiAuxRoutines.o: MpiAuxRoutines.F90  param.o
ImplicitAndUpdateVZ.o: ImplicitAndUpdateVZ.F90  param.o
StatRoutines.o: StatRoutines.F90  param.o
StatReadReduceWrite.o: StatReadReduceWrite.F90  param.o
HdfReadContinua.o: HdfReadContinua.F90  param.o
transpose_y_to_x.o: transpose_y_to_x.F90  
ImplicitAndUpdateTemp.o: ImplicitAndUpdateTemp.F90  param.o
decomp_2d.o: decomp_2d.F90 transpose_z_to_x.F90 transpose_x_to_z.F90 transpose_x_to_y.F90 transpose_y_to_x.F90 transpose_y_to_z.F90 transpose_z_to_y.F90
interp.o: interp.F90  param.o

