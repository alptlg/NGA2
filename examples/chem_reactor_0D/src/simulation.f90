!> chem_reactor_0D: 0D chemistry solver (homogeneous adiabatic isobaric).
!> Solves mass fraction and temperature temporal evolution at constant pressure, adiabatic.
!> Uses CVODE (stiff ODE solver) with adaptive internal stepping; output at user-specified "Time step" intervals.
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

   ! =======================================
   ! Variable declarations =================
   ! =======================================
   real(WP), allocatable :: Y(:), Y0(:), h(:), cp(:)
   real(WP), allocatable :: atom_masses_arr(:)
   integer, allocatable :: comp(:,:)
   real(WP) :: T, T0, P, rho, W_mix, dt, time, time_end
   real(WP) :: Cp_mix, Ysum, phi, F_A_st
   real(WP) :: n_O2_st, n_C, n_H, n_O, W_fuel, W_O2
   real(WP) :: Y_O2_air, Y_N2_air, W_O2_air, W_N2_air, W_air
   integer :: i, j, iu, iO2, iN2, ifuel, iC, iH, iO, a
   integer :: reactor_type  ! 1=isobar, 2=isochor
   character(len=str_medium) :: output_file, fuel_name, tag, reactor_str
   character(len=str_medium), dimension(nS) :: species_names
   character(len=2), dimension(:), allocatable :: atom_names_arr

   ! CVODE variables
   type(C_PTR) :: cvode_mem, sunctx
   type(N_Vector), pointer :: yvec
   type(SUNMatrix), pointer :: sunmat_A
   type(SUNLinearSolver), pointer :: sunlinsol_LS
   real(C_DOUBLE), dimension(nT), target :: ydata
   real(C_DOUBLE) :: tstart, tstop
   real(C_DOUBLE) :: tret(1)
   integer(C_INT) :: ierr
   integer(C_LONG) :: flag

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
   call param_read('End time', time_end)
   call param_read('Time step', dt)
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

   ! Write initial state to output file
   if (amRoot) then
      W_mix = 1.0_WP / sum(Y(1:nS) / W_sp(1:nS))
      if (reactor_type .eq. 1) then
         rho = P * W_mix / (Rcst * T)
      else
         rho = hr_ic_rho
         P = rho * Rcst * T / W_mix
      end if
      write (iu, '(ES18.10)', advance='no') 0.0_WP
      do i = 1, nS
         write (iu, '(1X,ES18.10)', advance='no') merge(0.0_WP, Y(i), Y(i) < Y_floor_OUT)
      end do
      write (iu, '(1X,ES18.10,1X,ES18.10,1X,ES18.10)') T, P, rho
   end if

   ! Integrate: advance to each output time; CVODE uses adaptive internal stepping
   ! Do NOT modify ydata before CVODE step: Cantera integrates the ODE as-is without
   ! clipping/renormalization. Pre-step modification caused divergence from Cantera.
   time = 0.0_WP
   do while (time .lt. time_end)
      time = min(time + dt, time_end)
      tstop = real(time, C_DOUBLE)
      flag = FCVode(cvode_mem, tstop, yvec, tret(1), CV_NORMAL)
      if (flag .lt. 0) then
         if (amRoot) then
            write (*, '(A,I0,A)') '[chem_reactor_0D] CVODE failed with flag: ', flag, ' (see below)'
            call print_cvode_flag(flag)
         end if
         call die('[chem_reactor_0D] CVODE integration failed')
      end if

      ! Extract solution (ydata is pointer target of yvec, updated by CVODE)
      ! Use raw CVODE output; no clipping/renormalization (matches Cantera behavior)
      Y(1:nS) = real(ydata(1:nS), WP)
      T = real(ydata(nT), WP)

      ! Write this time step to output file
      if (amRoot) then
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
      end if
   end do

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
   if (allocated(atom_names_arr)) deallocate (atom_names_arr)
   if (allocated(atom_masses_arr)) deallocate (atom_masses_arr)
   if (allocated(comp)) deallocate (comp)
   call param_final
   call messager_final
   call parallel_final
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
