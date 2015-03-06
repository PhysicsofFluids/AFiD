bin_PROGRAMS = boutnp
boutnp_SOURCES = \
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

boutnp_LDADD = $(FFTW3_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) 

include dependencies
