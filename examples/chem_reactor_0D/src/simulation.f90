!> chem_reactor_0D: 0D chemistry solver (homogeneous adiabatic isobaric).
!> Solves mass fraction and temperature temporal evolution at constant pressure, adiabatic.
!> Uses CVODE (stiff ODE solver) with adaptive internal stepping.
!> Output: FM T0DIsoChor-style — "Number of outputs" = buffer capacity, step-strided SaveSolution + Reduce,
!> optional early stop after ignition + plateau (FlameMaster), optional "Impose max time" to integrate to Max time anyway.
!> Negative "Max time" = |value| and implies Impose max time (FlameMaster-compatible).
!> Chemical source terms: dY_i/dt = ydot_i from fcmech_get_ydot (mass-based production rate).
!> Uses NGA2 param parser and libraries.
!
!> Naming: hr_ib_ = homogeneous reactor isobaric (constant P). hr_ic_ = isochoric (constant V).

!> Module to pass pressure and constants to CVODE RHS callback (C binding cannot pass extra arguments).
module hr_ib_cvode_data
   use precision, only: WP
   real(WP), save :: hr_ib_P = 0.0_WP
   real(WP), parameter :: Y_floor = 1.0e-50_WP, Y_floor_OUT = 1.0e-20_WP
end module hr_ib_cvode_data

!> Module to pass density for isochoric CVODE RHS (constant volume).
module hr_ic_cvode_data
   use precision, only: WP
   real(WP), save :: hr_ic_rho = 0.0_WP
end module hr_ic_cvode_data

!> Module containing the RHS function for CVODE: dy/dt = f(t,y) (isobaric).
module hr_ib_rhs_mod
   use precision, only: WP
   use fcmech
   use hr_ib_cvode_data
   use fsundials_core_mod
   use fnvector_serial_mod
   use, intrinsic :: ISO_C_BINDING
   implicit none

   contains

   !> CVODE RHS callback: computes dy/dt = f(t,y).
   !> State y = [Y_1..Y_nS, T]; output f = [dY/dt, dT/dt].
   integer(c_int) function hr_ib_rhs_wrapper(t, sunvec_y, sunvec_f, user_data) result(ierr) bind(C, name='hr_ib_rhs_wrapper')
      real(c_double), value :: t
      type(N_Vector)        :: sunvec_y
      type(N_Vector)        :: sunvec_f
      type(c_ptr), value    :: user_data
      real(c_double), pointer :: yval(:), fval(:)
      real(WP), dimension(nS) :: h, cp, ydot
      real(WP) :: Cp_mix, W_mix
      ierr = 0_c_int
      ! Get pointers to CVODE state and output arrays
      yval => FN_VGetArrayPointer(sunvec_y)
      fval => FN_VGetArrayPointer(sunvec_f)
      ! Mixture molar mass: 1/W_mix = sum(Y_i/W_i)
      W_mix = 1.0_WP / sum(real(yval(1:nS), WP) / W_sp(1:nS))
      ! Enthalpy and Cp from NASA polynomials (h, cp in J/mol, J/(mol·K))
      call fcmech_get_thermodata(h, cp, real(yval(nS + 1), WP))
      ! Cp_mix in J/(kg·K): sum(Y_i * cp_i/W_i) for mass-based mixture Cp
      Cp_mix = sum(real(yval(1:nS), WP) * cp(1:nS) / W_sp(1:nS))
      ! Mass-based production rates: dY_i/dt = ydot_i (1/s)
      call fcmech_get_ydot(hr_ib_P, real(yval(nS + 1), WP), real(yval(1:nS), WP), ydot)
      fval(1:nS) = real(ydot(1:nS), c_double)
      ! Adiabatic: dT/dt = -sum(h_i/W_i * dY_i/dt) / Cp_mix, with Cp_mix in J/(kg·K)
      fval(nS+1) = real(-sum(h(1:nS) * ydot(1:nS) / W_sp(1:nS)) / Cp_mix, c_double)
   end function hr_ib_rhs_wrapper
end module hr_ib_rhs_mod

!> Module containing the RHS function for CVODE: dy/dt = f(t,y) (isochoric).
module hr_ic_rhs_mod
   use precision, only: WP
   use fcmech
   use hr_ic_cvode_data
   use fsundials_core_mod
   use fnvector_serial_mod
   use, intrinsic :: ISO_C_BINDING
   implicit none

   contains

   !> CVODE RHS callback: isochoric adiabatic. State y = [Y_1..Y_nS, T]; P = rho*R*T/W_mix.
   integer(c_int) function hr_ic_rhs_wrapper(t, sunvec_y, sunvec_f, user_data) result(ierr) bind(C, name='hr_ic_rhs_wrapper')
      real(c_double), value :: t
      type(N_Vector)        :: sunvec_y
      type(N_Vector)        :: sunvec_f
      type(c_ptr), value    :: user_data
      real(c_double), pointer :: yval(:), fval(:)
      real(WP), dimension(nS) :: h, cp, ydot
      real(WP) :: Cp_mix, Cv_mix, W_mix, Ploc
      ierr = 0_c_int
      yval => FN_VGetArrayPointer(sunvec_y)
      fval => FN_VGetArrayPointer(sunvec_f)
      W_mix = 1.0_WP / sum(real(yval(1:nS), WP) / W_sp(1:nS))
      call fcmech_get_thermodata(h, cp, real(yval(nS + 1), WP))
      Cp_mix = sum(real(yval(1:nS), WP) * cp(1:nS) / W_sp(1:nS))
      ! P = rho * R * T / W_mix (ideal gas, constant volume)
      Ploc = hr_ic_rho * Rcst * real(yval(nS + 1), WP) / W_mix
      call fcmech_get_ydot(Ploc, real(yval(nS + 1), WP), real(yval(1:nS), WP), ydot)
      fval(1:nS) = real(ydot(1:nS), c_double)
      ! Adiabatic isochoric: dT/dt = -sum((h_i - R*T)/W_i * ydot_i) / Cv_mix
      Cv_mix = Cp_mix - Rcst / W_mix
      if (Cv_mix .lt. 1.0e-30_WP) Cv_mix = 1.0_WP
      fval(nS+1) = real(-sum((h(1:nS) - Rcst * real(yval(nS + 1), WP)) * ydot(1:nS) / W_sp(1:nS)) / Cv_mix, c_double)
   end function hr_ic_rhs_wrapper
end module hr_ic_rhs_mod

!> Main program: 0D adiabatic chemistry reactor (isobaric or isochoric).
program chem_reactor_0D
   use precision, only: WP
   use string, only: str_medium
   use param, only: param_init, param_final, param_read, param_exists
   use parallel, only: parallel_init, parallel_final, amRoot
   use messager, only: messager_init, messager_final, die, warn
   use fcmech
   use hr_ib_cvode_data
   use hr_ib_rhs_mod
   use hr_ic_cvode_data
   use hr_ic_rhs_mod
   use fcvode_mod
   use fnvector_serial_mod
   use fsunmatrix_dense_mod
   use fsunlinsol_dense_mod
   use fsundials_core_mod
   use, intrinsic :: ISO_C_BINDING
   implicit none

   ! State vector size: Y(1:nS) = mass fractions, T = temperature; nT = nS+1
   integer, parameter :: nT = nS + 1

   ! CVODE tolerances and step limits
   real(C_DOUBLE), parameter :: cvode_rtol = 1.0e-12_WP
   real(C_DOUBLE), parameter :: cvode_atol = 1.0e-15_WP
   real(C_DOUBLE), parameter :: cvode_init_step = 1.0e-12_C_DOUBLE
   real(C_DOUBLE), parameter :: cvode_min_step = 1.0e-18_C_DOUBLE
   ! FM T0DIsoChor early-exit / plateau (same thresholds as FlameMaster T0DIsoChor.C)
   real(WP), parameter :: fm_rel_T_plateau = 1.0e-3_WP
   real(WP), parameter :: fm_dTdt_plateau_max = 8229.0_WP
   integer, parameter :: fm_nothing_happens_exit = 3
   real(WP), parameter :: fm_t_after_ign_mult = 3.0_WP
   real(WP), parameter :: fm_dT_ignition_K = 500.0_WP

   ! =======================================
   ! Variable declarations =================
   ! =======================================
   real(WP), allocatable :: Y(:), Y0(:), h(:), cp(:)
   real(WP), allocatable :: atom_masses_arr(:)
   integer, allocatable :: comp(:,:)
   real(WP) :: T, T0, P, rho, W_mix, time, time_end
   real(WP) :: Cp_mix, Ysum, phi, F_A_st
   real(WP) :: n_O2_st, n_C, n_H, n_O, W_fuel, W_O2
   real(WP) :: Y_O2_air, Y_N2_air, W_O2_air, W_N2_air, W_air
   integer :: i, j, iu, iO2, iN2, ifuel, iC, iH, iO, a
   integer :: reactor_type  ! 1=isobar, 2=isochor
   ! FM-style output buffer (FlameMaster T0DIsoChor): sized in Step 1; used in later steps
   integer :: Noutputs              ! user "Number of outputs" (max trajectory slots)
   integer :: N_slot                ! allocated length: Noutputs or 2*Noutputs+1 if .not. equidistant
   integer :: act_length            ! filled slots (FM fActLength); Step 2+ will update
   integer :: delta_step_save       ! FM fDeltaStepSave; doubles after each Reduce
   integer :: alloc_stat
   logical :: equidistant           ! FM fEquidistant; if .false., expand buffer (2*N+1)
   logical :: impose_max_time       ! FM fImposeTEnd: if .true., do not stop early after ignition
   logical :: ignited               ! FM ignition detected
   logical :: leave_run              ! FM leave: exit integration loop
   logical :: exited_fm_early        ! stopped on ignition+plateau (not Max time)
   character(len=str_medium) :: output_file, fuel_name, tag, reactor_str, eq_str
   character(len=str_medium) :: impose_str
   character(len=str_medium), dimension(nS) :: species_names
   character(len=2), dimension(:), allocatable :: atom_names_arr
   real(WP), allocatable :: time_hist(:)   ! stored physical time per slot
   real(WP), allocatable :: T_hist(:)      ! temperature history
   real(WP), allocatable :: Y_hist(:, :)   ! mass fractions (1:nS, 1:N_slot)

   ! CVODE variables
   type(C_PTR) :: cvode_mem, sunctx
   type(N_Vector), pointer :: yvec
   type(SUNMatrix), pointer :: sunmat_A
   type(SUNLinearSolver), pointer :: sunlinsol_LS
   real(C_DOUBLE), dimension(nT), target :: ydata
   real(C_DOUBLE) :: tstart
   real(C_DOUBLE) :: tret(1)
   integer(C_INT) :: ierr
   integer(C_LONG) :: flag
   integer(C_LONG), target :: nst_arr(1)
   integer :: irow
   integer :: nothing_happens       ! FM consecutive plateau segments
   real(WP) :: t_ignition           ! FM tIgnition [s]
   real(WP) :: dtemp_seg, delt_time_seg
   real(WP) :: T_prev_cv, time_prev_cv   ! previous CVODE state for FM plateau (consecutive steps)

   ! =======================================
   ! NGA2 initialization ====================
   ! =======================================
   call parallel_init
   call messager_init
   call param_init

   ! =======================================
   ! Read parameters from input file =======
   ! =======================================
   allocate (Y(nS), Y0(nS), h(nS), cp(nS))
   call fcmech_get_speciesnames(species_names)
   if (nA .gt. 0) then
      allocate (atom_names_arr(nA), atom_masses_arr(nA), comp(nA, nS))
      call fcmech_get_atomnames(atom_names_arr)
      call fcmech_get_atommasses(atom_masses_arr)
      call fcmech_get_composition(comp)
   end if

   ! =======================================
   ! Initial composition ====================
   ! Either: Fuel + Equivalence ratio (phi), or explicit Initial Y per species
   ! =======================================
   Y0 = 0.0_WP
   if (param_exists('Fuel') .and. param_exists('Equivalence ratio')) then
      call param_read('Fuel', fuel_name)
      call param_read('Equivalence ratio', phi)
      ! Find fuel, O2, N2 indices
      ifuel = -1
      iO2 = -1
      iN2 = -1
      do i = 1, nS
         if (trim(adjustl(species_names(i))) .eq. trim(adjustl(fuel_name))) ifuel = i
         if (trim(adjustl(species_names(i))) .eq. 'O2') iO2 = i
         if (trim(adjustl(species_names(i))) .eq. 'N2') iN2 = i
      end do
      ! Air mass fractions from mechanism molar masses (matches Cantera)
      W_O2_air = W_sp(iO2)
      W_N2_air = W_sp(iN2)
      W_air = 0.21_WP * W_O2_air + 0.79_WP * W_N2_air
      Y_O2_air = 0.21_WP * W_O2_air / W_air
      Y_N2_air = 0.79_WP * W_N2_air / W_air
      if (ifuel .le. 0 .or. iO2 .le. 0 .or. iN2 .le. 0) &
         call die('[chem_reactor_0D] Fuel, O2, or N2 not found in mechanism')
      ! (F/A)_st from species composition: n_O2_st = (2*n_C + n_H/2 - n_O)/2 moles O2 per mole fuel
      if (nA .le. 0) &
         call die('[chem_reactor_0D] Mechanism has no atom data; cannot compute F/A_st from composition')
      iC = -1
      iH = -1
      iO = -1
      do a = 1, nA
         if (trim(adjustl(atom_names_arr(a))) .eq. 'C') iC = a
         if (trim(adjustl(atom_names_arr(a))) .eq. 'H') iH = a
         if (trim(adjustl(atom_names_arr(a))) .eq. 'O') iO = a
      end do
      n_C = 0.0_WP
      n_H = 0.0_WP
      n_O = 0.0_WP
      if (iC .gt. 0) n_C = real(comp(iC, ifuel), WP)
      if (iH .gt. 0) n_H = real(comp(iH, ifuel), WP)
      if (iO .gt. 0) n_O = real(comp(iO, ifuel), WP)
      n_O2_st = (2.0_WP * n_C + n_H / 2.0_WP - n_O) / 2.0_WP
      if (n_O2_st .le. 0.0_WP) &
         call die('[chem_reactor_0D] Fuel has no combustible content or invalid composition')
      W_fuel = W_sp(ifuel)
      W_O2 = W_sp(iO2)
      F_A_st = W_fuel * Y_O2_air / (n_O2_st * W_O2)
      ! Y_fuel = phi*(F/A)_st / (1 + phi*(F/A)_st), Y_O2 = Y_O2_air/(1+phi*(F/A)_st), Y_N2 = Y_N2_air/(1+phi*(F/A)_st)
      Y0(ifuel) = phi * F_A_st / (1.0_WP + phi * F_A_st)
      Y0(iO2) = Y_O2_air / (1.0_WP + phi * F_A_st)
      Y0(iN2) = Y_N2_air / (1.0_WP + phi * F_A_st)
   else
      ! Per-species Initial Y: only specify non-zero species
      do i = 1, nS
         tag = 'Initial Y '//trim(adjustl(species_names(i)))
         if (param_exists(tag)) then
            call param_read(tag, Y0(i))
         end if
      end do
      Ysum = sum(Y0)
      Y0 = Y0 / Ysum
   end if

   Y = Y0
   call param_read('Temperature', T0)
   call param_read('Pressure', P)
   reactor_str = 'isobar'
   if (param_exists('Reactor')) call param_read('Reactor', reactor_str)
   reactor_str = trim(adjustl(reactor_str))
   do i = 1, len_trim(reactor_str)
      j = ichar(reactor_str(i:i))
      if (j .ge. 65 .and. j .le. 90) reactor_str(i:i) = char(j + 32)
   end do
   reactor_type = 1
   if (index(reactor_str, 'isochor') .gt. 0) reactor_type = 2

   ! Integration horizon: "Max time" (alias: "End time"). Negative => |Max time| + Impose max time (FM).
   if (param_exists('Max time')) then
      call param_read('Max time', time_end)
   else if (param_exists('End time')) then
      call param_read('End time', time_end)
   else
      call die('[chem_reactor_0D] Must specify "Max time" (or alias "End time")')
   end if
   impose_max_time = .false.
   if (time_end .lt. 0.0_WP) then
      time_end = abs(time_end)
      impose_max_time = .true.
      if (amRoot) write (*, '(A,ES18.10)') '[chem_reactor_0D] Negative Max time: imposing integration until TEnd = ', time_end
   end if
   if (param_exists('Impose max time')) then
      impose_str = 'no'
      call param_read('Impose max time', impose_str)
      impose_str = trim(adjustl(impose_str))
      if (len_trim(impose_str) .gt. 0) then
         j = ichar(impose_str(1:1))
         if (j .ge. 65 .and. j .le. 90) j = j + 32
         if (j .eq. ichar('y') .or. j .eq. ichar('t') .or. trim(impose_str) .eq. '1') impose_max_time = .true.
         if (j .eq. ichar('f') .or. j .eq. ichar('n') .or. trim(impose_str) .eq. '0') impose_max_time = .false.
      end if
   end if

   if (.not. param_exists('Number of outputs')) &
      call die('[chem_reactor_0D] Must specify "Number of outputs" (FlameMaster-style control)')
   call param_read('Number of outputs', Noutputs)
   if (Noutputs .lt. 1) call die('[chem_reactor_0D] Number of outputs must be >= 1')

   equidistant = .true.
   eq_str = 'yes'
   if (param_exists('Equidistant')) then
      call param_read('Equidistant', eq_str)
      eq_str = trim(adjustl(eq_str))
      if (len_trim(eq_str) .gt. 0) then
         j = ichar(eq_str(1:1))
         if (j .ge. 65 .and. j .le. 90) j = j + 32
         if (j .eq. ichar('f') .or. j .eq. ichar('n') .or. trim(eq_str) .eq. '0') equidistant = .false.
      end if
   end if
   if (equidistant) then
      N_slot = Noutputs
   else
      N_slot = 2 * Noutputs + 1
   end if
   if (N_slot .lt. 2) call die('[chem_reactor_0D] Number of outputs must be >= 2 (FM buffer / Reduce)')
   allocate (time_hist(N_slot), T_hist(N_slot), Y_hist(nS, N_slot), stat=alloc_stat)
   if (alloc_stat .ne. 0) call die('[chem_reactor_0D] Allocation of FM output buffers failed')
   act_length = 0
   delta_step_save = 1
   nothing_happens = 0
   ignited = .false.
   t_ignition = 0.0_WP

   call param_read('Output file', output_file, short='o', default='results_hr.out')
   T = T0
   ! Y < Y_floor => 0; renormalize
   Y(1:nS) = merge(0.0_WP, Y(1:nS), Y(1:nS) < Y_floor)
   Ysum = sum(Y(1:nS))
   Y(1:nS) = Y(1:nS) / Ysum
   if (reactor_type .eq. 1) then
      hr_ib_P = P
   else
      W_mix = 1.0_WP / sum(Y(1:nS) / W_sp(1:nS))
      hr_ic_rho = P * W_mix / (Rcst * T)
   end if

   ! =======================================
   ! Open output file ======================
   ! Header: time, species-Y1..Y_nS, T, P, rho
   ! =======================================
   if (amRoot) then
      iu = 50
      open (unit=iu, file=trim(output_file), status='replace', action='write')
      ! Header: A18 format to match ES18.10 column width (18 chars per column)
      write (iu, '(A18)', advance='no') 'time'
      do i = 1, nS
         write (tag, '(A,I0)') trim(adjustl(species_names(i)))//'-Y', i
         write (iu, '(1X,A18)', advance='no') trim(adjustl(tag))
      end do
      write (iu, '(1X,A18,1X,A18,1X,A18)') 'T', 'P', 'rho'
   end if

   ! =======================================
   ! CVODE integration ======================
   ! State vector: ydata(1:nS)=Y (mass fractions), ydata(nT)=T (temperature)
   ! =======================================
   ydata(1:nS) = real(Y(1:nS), C_DOUBLE)
   ydata(nT) = real(T, C_DOUBLE)

   ierr = FSUNContext_Create(SUN_COMM_NULL, sunctx)
   if (ierr .ne. 0) call die('[chem_reactor_0D] FSUNContext_Create failed')

   yvec => FN_VMake_Serial(int(nT, c_int64_t), ydata, sunctx)
   if (.not. associated(yvec)) call die('[chem_reactor_0D] FN_VMake_Serial failed')

   cvode_mem = FCVodeCreate(CV_BDF, sunctx)
   if (.not. c_associated(cvode_mem)) call die('[chem_reactor_0D] FCVodeCreate failed')

   tstart = 0.0_C_DOUBLE
   if (reactor_type .eq. 1) then
      ierr = FCVodeInit(cvode_mem, c_funloc(hr_ib_rhs_wrapper), tstart, yvec)
   else
      ierr = FCVodeInit(cvode_mem, c_funloc(hr_ic_rhs_wrapper), tstart, yvec)
   end if
   if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeInit failed')

   ierr = FCVodeSStolerances(cvode_mem, cvode_rtol, cvode_atol)
   if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeSStolerances failed')

   sunmat_A => FSUNDenseMatrix(int(nT, c_int64_t), int(nT, c_int64_t), sunctx)
   if (.not. associated(sunmat_A)) call die('[chem_reactor_0D] FSUNDenseMatrix failed')
   sunlinsol_LS => FSUNLinSol_Dense(yvec, sunmat_A, sunctx)
   if (.not. associated(sunlinsol_LS)) call die('[chem_reactor_0D] FSUNLinSol_Dense failed')
   ierr = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
   if (ierr .ne. 0) call die('[chem_reactor_0D] FCVodeSetLinearSolver failed')

   ierr = FCVodeSetMaxNumSteps(cvode_mem, 500000_C_LONG)
   if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeSetMaxNumSteps failed')
   ierr = FCVodeSetInitStep(cvode_mem, cvode_init_step)
   if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeSetInitStep failed')
   ierr = FCVodeSetMinStep(cvode_mem, cvode_min_step)
   if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeSetMinStep failed')

   ierr = FCVodeSetStopTime(cvode_mem, real(time_end, C_DOUBLE))
   if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeSetStopTime failed')

   ! FM-style output: store initial state (t=0), then CV_ONE_STEP with snapshot / Reduce
   ! Do NOT modify ydata before CVODE step: Cantera integrates the ODE as-is without
   ! clipping/renormalization. Pre-step modification caused divergence from Cantera.
   !
   ! Plateau / early-exit: FM uses consecutive *stored* rows; with a small buffer, Reduce
   ! makes stores sparse and that breaks plateau detection. Using consecutive *CVODE*
   ! steps for the same FM inequalities matches “consecutive solution points” along the
   ! actual trajectory (see discussion in source comments above).
   time = 0.0_WP
   call fm_save_snapshot(0.0_WP, Y, T, 0, .true.)
   T_prev_cv = T_hist(1)
   time_prev_cv = time_hist(1)
   exited_fm_early = .false.

   do
      flag = FCVode(cvode_mem, real(time_end, C_DOUBLE), yvec, tret(1), CV_ONE_STEP)
      if (flag .lt. 0) then
         if (amRoot) then
            write (*, '(A,I0,A)') '[chem_reactor_0D] CVODE failed with flag: ', flag, ' (see below)'
            call print_cvode_flag(flag)
         end if
         call die('[chem_reactor_0D] CVODE integration failed')
      end if

      time = real(tret(1), WP)
      Y(1:nS) = real(ydata(1:nS), WP)
      T = real(ydata(nT), WP)

      ierr = FCVodeGetNumSteps(cvode_mem, nst_arr)
      if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeGetNumSteps failed')

      ! Plateau: consecutive CVODE steps (same FM thresholds as T0DIsoChor.C)
      dtemp_seg = T - T_prev_cv
      delt_time_seg = time - time_prev_cv
      if (delt_time_seg .gt. 0.0_WP) then
         if (abs(dtemp_seg / T) .lt. fm_rel_T_plateau .and. &
             abs(dtemp_seg / delt_time_seg) .lt. fm_dTdt_plateau_max) then
            nothing_happens = nothing_happens + 1
         else
            nothing_happens = 0
         end if
      end if
      T_prev_cv = T
      time_prev_cv = time

      ! Ignition: temperature vs initial stored T (no CH criterion)
      if (.not. ignited) then
         if (T .gt. T_hist(1) + fm_dT_ignition_K) then
            ignited = .true.
            t_ignition = time
            if (amRoot) write (*, '(A,ES14.6)') '[chem_reactor_0D] ignition (dT > 500 K) at t [s] = ', time
         end if
      end if

      ! FM OneStep saves before leave check
      call fm_save_snapshot(time, Y, T, int(nst_arr(1)), .false.)

      leave_run = .false.
      if (time .ge. time_end) leave_run = .true.
      if (ignited .and. nothing_happens .gt. fm_nothing_happens_exit .and. &
          time .gt. fm_t_after_ign_mult * t_ignition .and. .not. impose_max_time) then
         leave_run = .true.
         exited_fm_early = .true.
      end if

      if (leave_run) exit

      if (int(flag) .eq. CV_TSTOP_RETURN) exit

      if (act_length .eq. N_slot .and. time .lt. time_end) call fm_reduce_buffers
   end do

   ! FM post-loop append if integrator time passed last stored point. Skip when we already
   ! stopped on ignition+plateau: otherwise a huge delta_step_save can leave last store at
   ! ~1 ms while time is 0.2 s and this append recreates the spurious tail the user saw.
   if (act_length .gt. 0 .and. .not. exited_fm_early) then
      if (time .gt. time_hist(act_length)) then
         do while (act_length .ge. N_slot)
            call fm_reduce_buffers
         end do
         ierr = FCVodeGetNumSteps(cvode_mem, nst_arr)
         if (ierr .ne. CV_SUCCESS) call die('[chem_reactor_0D] FCVodeGetNumSteps failed (post-loop)')
         call fm_save_snapshot(time, Y, T, int(nst_arr(1)), .true.)
      end if
   end if

   if (amRoot .and. act_length .gt. 0) then
      write (*, '(A,ES14.6)') '[chem_reactor_0D] final time in output [s] = ', time_hist(act_length)
   end if

   ! Write trajectory from FM buffers (same columns as before)
   if (amRoot) then
      do irow = 1, act_length
         time = time_hist(irow)
         Y(1:nS) = Y_hist(1:nS, irow)
         T = T_hist(irow)
         W_mix = 1.0_WP / sum(Y(1:nS) / W_sp(1:nS))
         if (reactor_type .eq. 1) then
            rho = P * W_mix / (Rcst * T)
         else
            rho = hr_ic_rho
            P = rho * Rcst * T / W_mix
         end if
         write (iu, '(ES18.10)', advance='no') time
         do i = 1, nS
            write (iu, '(1X,ES18.10)', advance='no') merge(0.0_WP, Y(i), Y(i) < Y_floor_OUT)
         end do
         write (iu, '(1X,ES18.10,1X,ES18.10,1X,ES18.10)') T, P, rho
      end do
   end if

   ! Free CVODE and SUNDIALS resources
   call FCVodeFree(cvode_mem)
   ierr = FSUNLinSolFree(sunlinsol_LS)
   call FSUNMatDestroy(sunmat_A)
   call FN_VDestroy(yvec)
   ierr = FSUNContext_Free(sunctx)

   if (amRoot) close (iu)

   ! =======================================
   ! NGA2 termination ======================
   ! =======================================
   deallocate (Y, Y0, h, cp)
   if (allocated(time_hist)) deallocate (time_hist, T_hist, Y_hist)
   if (allocated(atom_names_arr)) deallocate (atom_names_arr)
   if (allocated(atom_masses_arr)) deallocate (atom_masses_arr)
   if (allocated(comp)) deallocate (comp)
   call param_final
   call messager_final
   call parallel_final

contains

   !> FlameMaster T0DIsoChor::SaveSolution — append one trajectory row if mod(n_step, delta_step_save)==0 or force_save.
   !> Host: time_hist, T_hist, Y_hist, act_length, N_slot, delta_step_save.
   subroutine fm_save_snapshot(t_phys, Y, temp_k, n_step, force_save)
      implicit none
      real(WP), intent(in) :: t_phys, temp_k
      real(WP), intent(in) :: Y(nS)
      integer, intent(in) :: n_step
      logical, intent(in), optional :: force_save
      integer :: idx
      logical :: do_save
      if (delta_step_save .lt. 1) call die('[chem_reactor_0D] fm_save_snapshot: delta_step_save < 1')
      do_save = (mod(n_step, delta_step_save) .eq. 0)
      if (present(force_save)) then
         if (force_save) do_save = .true.
      end if
      if (.not. do_save) return
      if (act_length .ge. N_slot) &
         call die('[chem_reactor_0D] FM output buffer full: call fm_reduce_buffers before fm_save_snapshot')
      idx = act_length + 1
      time_hist(idx) = t_phys
      T_hist(idx) = temp_k
      Y_hist(1:nS, idx) = Y(1:nS)
      act_length = act_length + 1
   end subroutine fm_save_snapshot

   !> FlameMaster T0DIsoChor::Reduce — in-place decimation: keep indices 1,3,5,... (every other slot).
   !> Sets act_length = (N_slot+1)/2 and doubles delta_step_save.
   subroutine fm_reduce_buffers
      implicit none
      integer :: to_len, i, j_src
      if (act_length .ne. N_slot) &
         call die('[chem_reactor_0D] fm_reduce_buffers: require act_length == N_slot (buffer full)')
      to_len = (N_slot + 1)/2
      do i = 1, to_len
         j_src = 2*i - 1
         time_hist(i) = time_hist(j_src)
         T_hist(i) = T_hist(j_src)
         Y_hist(1:nS, i) = Y_hist(1:nS, j_src)
      end do
      act_length = to_len
      delta_step_save = delta_step_save*2
   end subroutine fm_reduce_buffers

end program chem_reactor_0D

!> Print human-readable CVODE return flag.
subroutine print_cvode_flag(flag)
   use, intrinsic :: ISO_C_BINDING
   integer(C_LONG), intent(in) :: flag
   select case (int(flag))
   case (-9)
      write (*, '(A)') '  CV_FIRST_RHSFUNC_ERR: RHS failed at first call'
   case (-8)
      write (*, '(A)') '  CV_RHSFUNC_FAIL: RHS failed unrecoverably'
   case (-4)
      write (*, '(A)') '  CV_CONV_FAILURE: Newton convergence failure or min step size reached'
   case (-3)
      write (*, '(A)') '  CV_ERR_FAILURE: Error test failure or min step size reached'
   case (-1)
      write (*, '(A)') '  CV_TOO_MUCH_WORK: Max steps exceeded'
   case default
      write (*, '(A,I0,A)') '  Unknown CVODE flag: ', flag, ' (consult SUNDIALS docs)'
   end select
end subroutine print_cvode_flag
