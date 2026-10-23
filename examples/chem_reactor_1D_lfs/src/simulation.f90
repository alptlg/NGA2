!> chem_reactor_1D_lfs: steady freely-propagating premixed flame -> laminar flame speed S_L = u(inlet).
!>
!> Port of Cantera 3.0 FreeFlame (StFlow free flow, Inlet1D/Outlet1D, MultiNewton, MultiJac, Refiner, Sim1D and the
!> python auto-solve), with the source terms, thermo and transport of module fcmech (yaml2nga_lib_cf.py):
!>   unknowns per point: u, T, Y_1..Y_nS; mixture-averaged diffusion with correction flux, no Soret, no radiation;
!>   convection first-order upwind (Cantera) or second order ('Convection scheme : second-order', see deferred_source),
!>   central diffusion with transport properties at mid-points;
!>   energy: rho cp u dT/dz = d/dz(lambda dT/dz) - sum(wdot_k h_k) - sum(j_k dh_k/dz / W_k)  (dh_k/dz upwind; central
!>   with second order);
!>   flame anchored by a fixed-temperature point (continuity replaced by T = T_fixed there);
!>   damped Newton with finite-difference block-tridiagonal Jacobian (frozen transport), time-stepping fallback,
!>   grid refinement (ratio / slope / curve / prune), automatic start from a coarse grid (as FreeFlame.solve(auto=True))
!>   or start from a given profile (as set_initial_guess(data=...) + solve(refine_grid=True)).
!> Standalone program: NGA2 libraries (precision, param, messager, parallel) + LAPACK, not NGA2's core/main.f90.
!> Input (NGA2 param file): see input_single. Output: 'Result file' (key value) and 'Solution file' (z u T Y_k).

module flame1d
   use precision, only: WP
   use fcmech
   implicit none
   integer, parameter :: iU = 1, iT = 2, iY = 3, nv = nS + 2
   integer, parameter :: ncomp_ct = nv + 3      ! Cantera also carries V, lambda and E (always zero): enters the norm size
   real(WP), parameter :: Tiny = 1.0e-20_WP, Undef = -1.0e300_WP

   ! ---- problem
   real(WP) :: Pres, Tin, Yin(nS), Teq, Yeq(nS), Thigh = 4000.0_WP
   integer :: np = 0
   real(WP), allocatable :: zg(:), dz(:), xs(:,:), xprev(:,:), Tfix(:)
   logical :: energy = .true., width_check = .false., too_narrow = .false.
   real(WP) :: zfixed = Undef, tfixed = -1.0_WP, rdt = 0.0_WP
   ! ---- properties of the current residual evaluation
   real(WP), allocatable :: rho(:), wtm(:), cpm(:), hk(:,:), wdm(:,:), dmix(:,:), tcon(:), flux(:,:)
   real(WP), allocatable :: dcs(:,:)       ! deferred second-order convection source (zero for upwind)
   integer, allocatable :: mask(:,:)
   integer :: kexL = 1, kexR = 1
   ! ---- tolerances (Cantera defaults) and Jacobian
   real(WP) :: rtol_ss = 1.0e-4_WP, atol_ss = 1.0e-9_WP, rtol_ts = 1.0e-4_WP, atol_ts = 1.0e-11_WP
   real(WP), parameter :: jac_rtol = 1.0e-5_WP
   real(WP) :: jac_atol
   real(WP), allocatable :: JA(:,:,:), JB(:,:,:), JC(:,:,:), diag0(:,:), LUf(:,:,:), Gf(:,:,:)
   integer, allocatable :: ipiv(:,:)
   integer :: jac_age = 1000, max_age = 20, n_jac = 0, n_newton = 0, n_tsteps = 0
   logical :: jac_ok = .false., factored = .false.
   ! ---- refinement (set from input)
   real(WP) :: ratio = 3.0_WP, gslope = 0.08_WP, gcurve = 0.08_WP, gprune = 0.05_WP
   real(WP), parameter :: min_range = 0.01_WP, gridmin = 1.0e-10_WP
   integer :: npmax = 400
   ! ---- time stepping (Cantera defaults)
   real(WP), parameter :: tstep0 = 1.0e-5_WP, tmin = 1.0e-16_WP, tmax = 1.0e8_WP, tfactor = 0.5_WP
   integer, parameter :: nsteps_ts = 10, nsteps_max = 5000
   integer :: verbose = 0
   ! ---- convection: 0 = first-order upwind (Cantera); 1 = second order by deferred correction: Newton solves the
   !      upwind system plus the source dcs = rho u (upwind - second-order derivative) of the previous solution, repeated
   !      until the solution stops changing (deferred_source, sim_solve); dh_k/dz central
   integer :: conv_scheme = 0
   integer, parameter :: max_dc = 200      ! deferred-correction passes per grid
   real(WP) :: dc_alpha = 1.0_WP           ! fraction of the source applied (ramped 0.25 -> 1 when switching from upwind)

contains

   ! ============================================================================================ storage
   subroutine resize(n)
      integer, intent(in) :: n
      if (allocated(zg)) deallocate(zg, dz, xprev, Tfix, rho, wtm, cpm, hk, wdm, dmix, tcon, flux, mask, dcs)
      if (allocated(JA)) deallocate(JA, JB, JC, diag0, LUf, Gf, ipiv)
      np = n
      allocate(zg(n), dz(n-1), xprev(nv,n), Tfix(n), rho(n), wtm(n), cpm(n), hk(nS,n), wdm(nS,n), dmix(nS,n-1), &
               tcon(n-1), flux(nS,n-1), mask(nv,n), dcs(nv,n))
      allocate(JA(nv,nv,n), JB(nv,nv,n), JC(nv,nv,n), diag0(nv,n), LUf(nv,nv,n), Gf(nv,nv,n), ipiv(nv,n))
      mask = 0; dcs = 0.0_WP; jac_ok = .false.; factored = .false.
   end subroutine resize

   subroutine set_grid(znew, xnew)
      real(WP), intent(in) :: znew(:), xnew(:,:)
      if (allocated(xs)) deallocate(xs)
      call resize(size(znew))
      zg = znew
      dz = zg(2:np) - zg(1:np-1)
      allocate(xs(nv,np)); xs = xnew
   end subroutine set_grid

   !> StFlow::_finalize: fixed-temperature profile and fixed point bookkeeping
   subroutine finalize()
      integer :: j
      if (energy) Tfix = xs(iT,:)
      if (tfixed > 0.0_WP) then
         if (any(zg == zfixed)) return
         do j = 1, np - 1
            if ((xs(iT,j) - tfixed)*(xs(iT,j+1) - tfixed) <= 0.0_WP) then
               tfixed = xs(iT,j+1); zfixed = zg(j+1); return
            end if
         end do
      end if
   end subroutine finalize

   ! ============================================================================================ residual
   !> Properties: node thermo and source terms; midpoint transport and excess species only if full
   subroutine properties(xx, full)
      real(WP), intent(in) :: xx(nv,np)
      logical, intent(in) :: full
      real(WP) :: h(nS), cp(nS), yd(nS), mu(nS), lam(nS), invD(nS,nS), Ym(nS), Xm(nS), Tm, Wm, sm
      real(WP) :: Xj(nS), Xj1(nS), corr
      integer :: j, k
      do j = 1, np
         wtm(j) = 1.0_WP/sum(xx(iY:,j)/W_sp)
         rho(j) = Pres*wtm(j)/(Rcst*xx(iT,j))
         call fcmech_get_thermodata(h, cp, xx(iT,j))
         hk(:,j) = h
         cpm(j) = sum(xx(iY:,j)*cp/W_sp)
         call fcmech_get_ydot(Pres, xx(iT,j), xx(iY:,j), yd)
         wdm(:,j) = rho(j)*yd                                   ! W_k * wdot_k  [kg/(m^3 s)]
      end do
      if (full) then
         do j = 1, np - 1                                       ! StFlow::updateTransport (setGasAtMidpoint)
            Tm = 0.5_WP*(xx(iT,j) + xx(iT,j+1))
            Ym = 0.5_WP*(xx(iY:,j) + xx(iY:,j+1))
            Wm = 1.0_WP/sum(Ym/W_sp)
            Xm = max(Ym/W_sp*Wm, Tiny)                           ! MixTransport::update_C
            mu = 0.0_WP                                          ! (not used by fcmech_get_conductivity)
            call fcmech_get_conductivity(lam, Tm, mu)
            tcon(j) = 0.5_WP*(sum(Xm*lam) + 1.0_WP/sum(Xm/lam))
            call fcmech_get_invDij(invD, Tm, Pres)
            do k = 1, nS                                        ! GasTransport::getMixDiffCoeffs
               sm = sum(Xm*invD(:,k)) - Xm(k)*invD(k,k)
               dmix(k,j) = (Wm - Xm(k)*W_sp(k))/(Wm*sm)
            end do
         end do
         kexL = maxloc(xx(iY:,1), 1); kexR = maxloc(xx(iY:,np), 1)
      end if
      do j = 1, np - 1                                          ! StFlow::updateDiffFluxes (mixture-averaged)
         Xj = xx(iY:,j)*wtm(j)/W_sp
         Xj1 = xx(iY:,j+1)*wtm(j+1)/W_sp
         flux(:,j) = W_sp*(rho(j)*dmix(:,j)/wtm(j))*(Xj - Xj1)/dz(j)
         corr = -sum(flux(:,j))
         flux(:,j) = flux(:,j) + corr*xx(iY:,j)
      end do
   end subroutine properties

   !> StFlow::evalResidual + Inlet1D/Outlet1D (free flame). full: update transport, excess species and mask.
   subroutine residual(xx, rr, rdt_l, full)
      real(WP), intent(in) :: xx(nv,np), rdt_l
      real(WP), intent(out) :: rr(nv,np)
      logical, intent(in) :: full
      real(WP) :: dYdz, dTdz, divq, sm1, sm2, dhk, convec, diffus, mdot, cm, c0, cq, hl, hr
      integer :: j, k, jl, jr
      call properties(xx, full)
      if (full) mask = 0
      ! ---- left boundary (j = 1) + inlet
      rr(iU,1) = -(rho(2)*xx(iU,2) - rho(1)*xx(iU,1))/dz(1)
      rr(iT,1) = xx(iT,1) - merge(Tin, Tfix(1), energy)
      mdot = rho(1)*xx(iU,1)
      do k = 1, nS
         rr(iY+k-1,1) = -(flux(k,1) + rho(1)*xx(iU,1)*xx(iY+k-1,1))
         if (k /= kexL) rr(iY+k-1,1) = rr(iY+k-1,1) + mdot*Yin(k)
      end do
      rr(iY+kexL-1,1) = 1.0_WP - sum(xx(iY:,1))
      ! ---- interior
      do j = 2, np - 1
         ! continuity (free flame: anchored at the fixed-temperature point)
         if (zg(j) > zfixed) then
            rr(iU,j) = -(rho(j)*xx(iU,j) - rho(j-1)*xx(iU,j-1))/dz(j-1)
         else if (zg(j) == zfixed) then
            if (energy) then
               rr(iU,j) = xx(iT,j) - tfixed
            else
               rr(iU,j) = rho(j)*xx(iU,j) - rho(1)*0.3_WP        ! as Cantera StFlow: fixed mass flux, energy off
            end if
         else
            rr(iU,j) = -(rho(j+1)*xx(iU,j+1) - rho(j)*xx(iU,j))/dz(j)
         end if
         ! convection: upwind (Cantera); second order adds the deferred source dcs (see deferred_source)
         if (xx(iU,j) > 0.0_WP) then
            jl = j - 1; jr = j
         else
            jl = j; jr = j + 1
         end if
         ! species
         do k = 1, nS
            dYdz = (xx(iY+k-1,jr) - xx(iY+k-1,jl))/dz(jl)
            convec = rho(j)*xx(iU,j)*dYdz
            diffus = 2.0_WP*(flux(k,j) - flux(k,j-1))/(zg(j+1) - zg(j-1))
            if (conv_scheme == 1) then                       ! (separate expressions keep upwind bit-identical to Cantera's)
               rr(iY+k-1,j) = (wdm(k,j) - convec - diffus + dcs(iY+k-1,j))/rho(j) - rdt_l*(xx(iY+k-1,j) - xprev(iY+k-1,j))
            else
               rr(iY+k-1,j) = (wdm(k,j) - convec - diffus)/rho(j) - rdt_l*(xx(iY+k-1,j) - xprev(iY+k-1,j))
            end if
         end do
         if (full) mask(iY:,j) = 1
         ! energy
         if (energy) then
            dTdz = (xx(iT,jr) - xx(iT,jl))/dz(jl)
            if (conv_scheme == 1) then
               hl = dz(j-1); hr = dz(j)
               cm = -hr/(hl*(hl + hr)); c0 = (hr - hl)/(hl*hr); cq = hl/(hr*(hl + hr))
            end if
            sm1 = 0.0_WP; sm2 = 0.0_WP
            do k = 1, nS
               if (conv_scheme == 1) then
                  dhk = cm*hk(k,j-1) + c0*hk(k,j) + cq*hk(k,j+1)
               else
                  dhk = (hk(k,jr) - hk(k,jl))/dz(jl)
               end if
               sm1 = sm1 + wdm(k,j)/W_sp(k)*hk(k,j)
               sm2 = sm2 + 0.5_WP*(flux(k,j-1) + flux(k,j))*dhk/W_sp(k)
            end do
            divq = -2.0_WP*(tcon(j)*(xx(iT,j+1) - xx(iT,j))/dz(j) - tcon(j-1)*(xx(iT,j) - xx(iT,j-1))/dz(j-1)) &
                   /(zg(j+1) - zg(j-1))
            if (conv_scheme == 1) then
               rr(iT,j) = (-cpm(j)*rho(j)*xx(iU,j)*dTdz - divq - sm1 - sm2 + dcs(iT,j))/(rho(j)*cpm(j)) &
                          - rdt_l*(xx(iT,j) - xprev(iT,j))
            else
               rr(iT,j) = (-cpm(j)*rho(j)*xx(iU,j)*dTdz - divq - sm1 - sm2)/(rho(j)*cpm(j)) &
                          - rdt_l*(xx(iT,j) - xprev(iT,j))
            end if
            if (full) mask(iT,j) = 1
         else
            rr(iT,j) = xx(iT,j) - Tfix(j)
         end if
      end do
      ! ---- right boundary (j = np) + outlet: zero gradients
      j = np
      rr(iU,j) = rho(j)*xx(iU,j) - rho(j-1)*xx(iU,j-1)
      if (energy) then
         rr(iT,j) = xx(iT,j) - xx(iT,j-1)
      else
         rr(iT,j) = xx(iT,j) - Tfix(j)
      end if
      do k = 1, nS
         rr(iY+k-1,j) = xx(iY+k-1,j) - xx(iY+k-1,j-1)
      end do
      rr(iY+kexR-1,j) = 1.0_WP - sum(xx(iY:,j))
   end subroutine residual

   !> Weight of the central derivative for cell Peclet number pe = |u| h / D: w = 1/(1 + (pe/2)^2).
   !> pe << 1 (resolved flame): w = 1 - pe^2/4, central (added diffusion ~ h^3); large pe (coarse, convection-dominated
   !> cells): -> upwind. w <= 1/pe (equality at pe = 2) keeps the downstream coefficient non-negative (no wiggles).
   pure real(WP) function fitted_weight(pe)
      real(WP), intent(in) :: pe
      fitted_weight = 1.0_WP/(1.0_WP + 0.25_WP*pe*pe)
   end function fitted_weight

   !> Second order by deferred correction: dcs = rho u (upwind - second-order derivative) at the current solution xs, so
   !> that upwind residual + dcs = second-order residual. Second-order derivative: 3-point central (non-uniform grid)
   !> blended toward upwind per variable by its cell Peclet number (fitted_weight) and a slope limiter (upwind at local
   !> extrema, upwind in under-resolved tails: no undershoot), then shifted by Y_k x (sum of the changes) so that the
   !> source adds no total mass (sum_k dcs_k = 0).
   subroutine deferred_source()
      real(WP) :: dupw(nS), d2(nS), cm, c0, cq, hl, hr, w, dTu, dTc
      integer :: j, k, jl, jr
      dcs = 0.0_WP
      if (conv_scheme /= 1) return
      call properties(xs, .true.)
      do j = 2, np - 1
         if (xs(iU,j) > 0.0_WP) then
            jl = j - 1; jr = j
         else
            jl = j; jr = j + 1
         end if
         hl = dz(j-1); hr = dz(j)
         cm = -hr/(hl*(hl + hr)); c0 = (hr - hl)/(hl*hr); cq = hl/(hr*(hl + hr))
         do k = 1, nS
            dupw(k) = (xs(iY+k-1,jr) - xs(iY+k-1,jl))/dz(jl)
            w = fitted_weight(abs(xs(iU,j))*dz(jl)/dmix(k,jl))*limiter(xs(iY+k-1,j-1:j+1))
            d2(k) = dupw(k) + w*(cm*xs(iY+k-1,j-1) + c0*xs(iY+k-1,j) + cq*xs(iY+k-1,j+1) - dupw(k))
         end do
         d2 = d2 - xs(iY:,j)*(sum(d2) - sum(dupw))
         dcs(iY:,j) = rho(j)*xs(iU,j)*(dupw - d2)
         if (energy) then
            dTu = (xs(iT,jr) - xs(iT,jl))/dz(jl)
            dTc = cm*xs(iT,j-1) + c0*xs(iT,j) + cq*xs(iT,j+1)
            w = fitted_weight(abs(xs(iU,j))*dz(jl)*rho(j)*cpm(j)/tcon(jl))*limiter(xs(iT,j-1:j+1))
            dcs(iT,j) = rho(j)*cpm(j)*xs(iU,j)*w*(dTu - dTc)
         end if
      end do
      dcs = dc_alpha*dcs

   contains

      !> min(1, van Leer) of the ratio of downstream to upstream slope (0 at extrema, 1 where the profile is smooth),
      !> times 1/(1 + q^2), q = largest change to a neighbour / |value| (-> upwind in tails resolved only relative to
      !> the variable's full range, e.g. tiny upstream radical tails: no undershoot, no negative mass fractions)
      real(WP) function limiter(f)
         real(WP), intent(in) :: f(3)
         real(WP) :: gu, gd, r, q
         q = max(abs(f(2) - f(1)), abs(f(3) - f(2)))/max(abs(f(2)), 1.0e-300_WP)
         gu = (f(2) - f(1))/hl; gd = (f(3) - f(2))/hr
         if (jl == j) then
            r = gu; gu = gd; gd = r                          ! u < 0: upstream is the right side
         end if
         if (gu*gd <= 0.0_WP) then
            limiter = 0.0_WP
         else
            r = gd/gu
            limiter = min(1.0_WP, 2.0_WP*r/(1.0_WP + r))/(1.0_WP + q*q)
         end if
      end function limiter

   end subroutine deferred_source

   ! ============================================================================================ Jacobian
   !> MultiJac::eval at the steady state (rdt = 0), transport frozen at xx; columns perturbed three points at a time
   subroutine jacobian(xx)
      real(WP), intent(in) :: xx(nv,np)
      real(WP) :: r0(nv,np), r1(nv,np), x1(nv,np), dxj(np), xv, d
      integer :: n, c, j
      call residual(xx, r0, 0.0_WP, .true.)
      JA = 0.0_WP; JB = 0.0_WP; JC = 0.0_WP
      do n = 1, nv
         do c = 1, 3
            x1 = xx
            do j = c, np, 3
               xv = xx(n,j)
               d = merge(xv*jac_rtol + jac_atol, xv*jac_rtol - jac_atol, xv >= 0.0_WP)
               x1(n,j) = xv + d
               dxj(j) = x1(n,j) - xv
            end do
            call residual(x1, r1, 0.0_WP, .false.)
            do j = c, np, 3
               JA(:,n,j) = (r1(:,j) - r0(:,j))/dxj(j)
               if (j > 1) JB(:,n,j-1) = (r1(:,j-1) - r0(:,j-1))/dxj(j)
               if (j < np) JC(:,n,j+1) = (r1(:,j+1) - r0(:,j+1))/dxj(j)
            end do
         end do
      end do
      do j = 1, np
         do n = 1, nv
            diag0(n,j) = JA(n,n,j)
         end do
      end do
      n_jac = n_jac + 1; jac_age = 0; factored = .false.
   end subroutine jacobian

   !> MultiJac::updateTransient
   subroutine update_transient()
      integer :: j, n
      do j = 1, np
         do n = 1, nv
            JA(n,n,j) = diag0(n,j) - mask(n,j)*rdt
         end do
      end do
      factored = .false.
   end subroutine update_transient

   !> Block-tridiagonal LU (block Thomas): LUf_j = LU(A_j - C_j G_{j-1}), G_j = LUf_j^-1 B_j
   subroutine factor(info)
      integer, intent(out) :: info
      integer :: j
      external :: dgetrf, dgetrs, dgemm
      LUf(:,:,1) = JA(:,:,1)
      call dgetrf(nv, nv, LUf(:,:,1), nv, ipiv(:,1), info)
      if (info /= 0) return
      do j = 1, np - 1
         Gf(:,:,j) = JB(:,:,j)
         call dgetrs('N', nv, nv, LUf(:,:,j), nv, ipiv(:,j), Gf(:,:,j), nv, info)
         LUf(:,:,j+1) = JA(:,:,j+1)
         call dgemm('N', 'N', nv, nv, nv, -1.0_WP, JC(:,:,j+1), nv, Gf(:,:,j), nv, 1.0_WP, LUf(:,:,j+1), nv)
         call dgetrf(nv, nv, LUf(:,:,j+1), nv, ipiv(:,j+1), info)
         if (info /= 0) return
      end do
      factored = .true.
   end subroutine factor

   subroutine block_solve(b, info)
      real(WP), intent(inout) :: b(nv,np)
      integer, intent(out) :: info
      integer :: j
      external :: dgetrs, dgemv
      info = 0
      if (.not. factored) then
         call factor(info)
         if (info /= 0) return
      end if
      call dgetrs('N', nv, 1, LUf(:,:,1), nv, ipiv(:,1), b(:,1), nv, info)
      do j = 2, np
         call dgemv('N', nv, nv, -1.0_WP, JC(:,:,j), nv, b(:,j-1), 1, 1.0_WP, b(:,j), 1)
         call dgetrs('N', nv, 1, LUf(:,:,j), nv, ipiv(:,j), b(:,j), nv, info)
      end do
      do j = np - 1, 1, -1
         call dgemv('N', nv, nv, -1.0_WP, Gf(:,:,j), nv, b(:,j+1), 1, 1.0_WP, b(:,j), 1)
      end do
   end subroutine block_solve

   ! ============================================================================================ Newton
   !> weighted norm of a step (MultiNewton::wnorm)
   real(WP) function wnorm(xx, st)
      real(WP), intent(in) :: xx(nv,np), st(nv,np)
      real(WP) :: ewt, rt, at, ssum
      integer :: n
      rt = merge(rtol_ss, rtol_ts, rdt == 0.0_WP); at = merge(atol_ss, atol_ts, rdt == 0.0_WP)
      ssum = 0.0_WP
      do n = 1, nv
         ewt = rt*sum(abs(xx(n,:)))/np + at
         ssum = ssum + sum((st(n,:)/ewt)**2)
      end do
      wnorm = sqrt(ssum/real(np*ncomp_ct + 2, WP))
   end function wnorm

   !> Newton step: st = -J^-1 F(xx)
   subroutine newton_step(xx, st, info)
      real(WP), intent(in) :: xx(nv,np)
      real(WP), intent(out) :: st(nv,np)
      integer, intent(out) :: info
      call residual(xx, st, rdt, .true.)
      st = -st
      call block_solve(st, info)
   end subroutine newton_step

   !> damping factor that keeps the step within the bounds (bound_step)
   real(WP) function bound_step(xx, st)
      real(WP), intent(in) :: xx(nv,np), st(nv,np)
      real(WP) :: lo(nv), hi(nv), v, vn
      integer :: n, j
      lo = -1.0e-7_WP; hi = 1.0e5_WP
      lo(iU) = -1.0e20_WP; hi(iU) = 1.0e20_WP
      lo(iT) = 200.0_WP; hi(iT) = Thigh
      bound_step = 1.0_WP
      do j = 1, np
         do n = 1, nv
            v = xx(n,j); vn = v + st(n,j)
            if (vn > hi(n)) then
               bound_step = max(0.0_WP, min(bound_step, (hi(n) - v)/(vn - v)))
            else if (vn < lo(n)) then
               bound_step = min(bound_step, (v - lo(n))/(v - vn))
            end if
         end do
      end do
   end function bound_step

   !> Verbose 2: the variable that blocks a Newton step at its bound
   subroutine report_bound(xx, st)
      real(WP), intent(in) :: xx(nv,np), st(nv,np)
      integer :: n, j
      do j = 1, np
         do n = iT, nv
            if (xx(n,j) + st(n,j) < merge(200.0_WP, -1.0e-7_WP, n == iT) .and. &
                xx(n,j) <= merge(200.0_WP, -1.0e-7_WP, n == iT) + 1.0e-12_WP) &
               write(*,'(a,i5,a,i4,a,es10.3,a,es10.3,a,es10.3,a,es10.3)') '  at bound: point ', j, ' var ', n, &
                  ' z=', zg(j), ' x=', xx(n,j), ' step=', st(n,j), ' dcs=', dcs(n,j)
         end do
      end do
   end subroutine report_bound

   !> MultiNewton::dampStep: 1 converged, 0 keep iterating, <0 failure
   integer function damp_step(x0, st0, x1, st1, s1)
      real(WP), intent(in) :: x0(nv,np), st0(nv,np)
      real(WP), intent(out) :: x1(nv,np), st1(nv,np), s1
      real(WP) :: s0, fb, damp
      integer :: m, info
      s0 = wnorm(x0, st0)
      fb = bound_step(x0, st0)
      if (fb < 1.0e-10_WP) then
         if (verbose > 1) call report_bound(x0, st0)
         damp_step = -3; return
      end if
      damp = 1.0_WP
      do m = 1, 7
         x1 = x0 + fb*damp*st0
         call newton_step(x1, st1, info)
         if (info /= 0) then
            damp_step = -4; return
         end if
         s1 = wnorm(x1, st1)
         if (s1 < 1.0_WP .or. s1 < s0) exit
         damp = damp/sqrt(2.0_WP)
      end do
      if (m <= 7) then
         damp_step = merge(0, 1, s1 > 1.0_WP)
      else
         damp_step = -2
         if (verbose > 1) call report_step(x0, st0)
      end if
   end function damp_step

   !> Verbose 2: where the (failed) Newton step is largest, in units of the error weight
   subroutine report_step(xx, st)
      real(WP), intent(in) :: xx(nv,np), st(nv,np)
      real(WP) :: e(nv,np), rt, at
      integer :: n, loc(2)
      rt = merge(rtol_ss, rtol_ts, rdt == 0.0_WP); at = merge(atol_ss, atol_ts, rdt == 0.0_WP)
      do n = 1, nv
         e(n,:) = abs(st(n,:))/(rt*sum(abs(xx(n,:)))/np + at)
      end do
      loc = maxloc(e)
      write(*,'(a,i5,a,i4,a,es10.3,a,es10.3,a,es10.3,a,es10.3)') '  max step: point ', loc(2), ' var ', loc(1), &
         ' z=', zg(loc(2)), ' x=', xx(loc(1),loc(2)), ' step=', st(loc(1),loc(2)), ' /ewt=', maxval(e)
   end subroutine report_step

   !> OneDim::solve + MultiNewton::solve. Returns >= 0 on success (x1 = solution), < 0 on failure.
   integer function newton_solve(x0, x1)
      real(WP), intent(in) :: x0(nv,np)
      real(WP), intent(out) :: x1(nv,np)
      real(WP) :: xm(nv,np), st(nv,np), st1(nv,np), s1
      integer :: m, info, j0, nre
      logical :: force
      if (.not. jac_ok) then
         call jacobian(x0); call update_transient(); jac_ok = .true.
      end if
      xm = x0; force = .false.; nre = 0; j0 = n_jac; m = -1
      do
         if (jac_age > max_age) force = .true.
         if (force) then
            call jacobian(xm); call update_transient(); force = .false.
         end if
         call newton_step(xm, st, info)
         if (info /= 0) then
            m = -4; exit
         end if
         jac_age = jac_age + 1
         n_newton = n_newton + 1
         m = damp_step(xm, st, x1, st1, s1)
         if (m == 0) then
            xm = x1
         else if (m == 1) then
            exit
         else
            if (jac_age > 1) then
               force = .true.
               if (nre > 3) exit
               nre = nre + 1
            else
               exit
            end if
         end if
      end do
      if (m < 0) x1 = xm
      if (m > 0 .and. n_jac == j0) m = 100
      newton_solve = m
   end function newton_solve

   ! ============================================================================================ time stepping
   subroutine set_rdt(new_rdt)
      real(WP), intent(in) :: new_rdt
      if (abs(new_rdt - rdt) > Tiny) then
         rdt = new_rdt
         if (jac_ok) call update_transient()
      end if
   end subroutine set_rdt

   !> OneDim::timeStep. Returns the last time step; ok = .false. if the integration failed.
   real(WP) function time_step(nsteps, dt0, ok)
      integer, intent(in) :: nsteps
      real(WP), intent(in) :: dt0
      logical, intent(out) :: ok
      real(WP) :: dt, x1(nv,np), ysum
      integer :: n, nfail, m, j
      dt = dt0; n = 0; nfail = 0; ok = .true.
      do while (n < nsteps)
         call set_rdt(1.0_WP/dt)
         xprev = xs
         m = newton_solve(xs, x1)
         if (m >= 0) then
            nfail = 0; n = n + 1; n_tsteps = n_tsteps + 1
            xs = x1
            if (m == 100) dt = dt*1.5_WP
            dt = min(dt, tmax)
            if (n_tsteps >= nsteps_max) then
               ok = .false.; exit
            end if
         else
            nfail = nfail + 1
            if (nfail > 2) then                                  ! resetBadValues
               do j = 1, np
                  xs(iY:,j) = max(xs(iY:,j), 0.0_WP)
                  ysum = sum(xs(iY:,j)); xs(iY:,j) = xs(iY:,j)/ysum
               end do
               nfail = 0
            else
               dt = dt*tfactor
               if (dt < tmin) then
                  ok = .false.; exit
               end if
            end if
         end if
      end do
      time_step = dt
   end function time_step

   ! ============================================================================================ grid
   !> Refiner::analyze + Sim1D::refine. Returns the number of new points (< 0: too many points).
   integer function refine()
      integer :: keep(np), n, j, i, nn
      logical :: loc(np)
      real(WP) :: v(np), s(np-1), vmin, vmax, smin, smax, aa, ss, dmax, r, thresh
      real(WP), allocatable :: zn(:), xn(:,:)
      if (np >= npmax) then
         refine = -1; return
      end if
      thresh = sqrt(epsilon(1.0_WP))
      keep = 0; keep(1) = 1; keep(np) = 1; loc = .false.
      do i = 1, nv
         v = xs(i,:)
         s = (v(2:np) - v(1:np-1))/dz
         vmin = minval(v); vmax = maxval(v); smin = minval(s); smax = maxval(s)
         aa = max(abs(vmax), abs(vmin)); ss = max(abs(smax), abs(smin))
         if ((vmax - vmin) > min_range*aa) then
            dmax = gslope*(vmax - vmin) + thresh
            do j = 1, np - 1
               r = abs(v(j+1) - v(j))/dmax
               if (r > 1.0_WP .and. dz(j) >= 2.0_WP*gridmin) loc(j) = .true.
               if (r >= gprune) then
                  keep(j) = 1; keep(j+1) = 1
               else if (keep(j) == 0) then
                  keep(j) = -1
               end if
            end do
         end if
         if ((smax - smin) > min_range*ss) then
            dmax = gcurve*(smax - smin)
            do j = 1, np - 2
               r = abs(s(j+1) - s(j))/(dmax + thresh/dz(j))
               if (r > 1.0_WP .and. dz(j) >= 2.0_WP*gridmin .and. dz(j+1) >= 2.0_WP*gridmin) then
                  loc(j) = .true.; loc(j+1) = .true.
               end if
               if (r >= gprune) then
                  keep(j+1) = 1
               else if (keep(j+1) == 0) then
                  keep(j+1) = -1
               end if
            end do
         end if
      end do
      do j = 2, np - 1                                           ! grid-ratio criteria (0-based j-1 in Cantera)
         if (dz(j) > ratio*dz(j-1)) then
            loc(j) = .true.
            keep(max(j-1,1):min(j+2,np)) = 1
         end if
         if (dz(j) < dz(j-1)/ratio) then
            loc(j-1) = .true.
            keep(max(j-2,1):min(j+1,np)) = 1
         end if
         if (j > 2) then
            if (zg(j+1) - zg(j-1) > ratio*dz(j-2)) keep(j) = 1
         end if
         if (j < np - 1) then
            if (zg(j+1) - zg(j-1) > ratio*dz(j+1)) keep(j) = 1
         end if
         if (zg(j) == zfixed) keep(j) = 1
      end do
      do j = 3, np - 1
         if (keep(j) == -1 .and. keep(j-1) == -1) keep(j) = 1
      end do
      nn = count(loc)                              ! as Cantera Sim1D::refine: flagged intervals + inserted points
      allocate(zn(2*np), xn(nv,2*np))
      n = 0
      do j = 1, np
         if (keep(j) /= -1) then
            n = n + 1; zn(n) = zg(j); xn(:,n) = xs(:,j)
            if (loc(j) .and. j < np) then
               n = n + 1; nn = nn + 1
               zn(n) = 0.5_WP*(zg(j) + zg(j+1)); xn(:,n) = 0.5_WP*(xs(:,j) + xs(:,j+1))
            end if
         end if
      end do
      call set_grid(zn(1:n), xn(:,1:n))
      call finalize()
      refine = nn
   end function refine

   !> Sim1D::setFixedTemperature (free flame): fixed point at an existing point or a new interpolated point
   subroutine set_fixed_temperature(t)
      real(WP), intent(in) :: t
      real(WP) :: t1, t2, th, zf, z1, z2, f
      real(WP), allocatable :: zn(:), xn(:,:)
      integer :: m, mf
      logical :: found
      mf = 0; found = .false.
      do m = 1, np - 1
         t1 = xs(iT,m); t2 = xs(iT,m+1)
         th = min(1.0_WP, 0.1_WP*(t2 - t1))
         z1 = zg(m); z2 = zg(m+1)
         if (abs(t - t1) <= th) then
            zf = z1; found = .true.
         else if (abs(t2 - t) <= th) then
            zf = z2; found = .true.
         else if (t1 < t .and. t < t2) then
            mf = m; zf = (z1 - z2)/(t1 - t2)*(t - t2) + z2; found = .true.
         end if
         if (found) exit
      end do
      if (.not. found) return
      zfixed = zf; tfixed = t
      if (mf > 0) then
         allocate(zn(np+1), xn(nv,np+1))
         zn(1:mf) = zg(1:mf); xn(:,1:mf) = xs(:,1:mf)
         f = (zf - zg(mf+1))/(zg(mf) - zg(mf+1))
         zn(mf+1) = zf; xn(:,mf+1) = f*(xs(:,mf) - xs(:,mf+1)) + xs(:,mf+1)
         zn(mf+2:) = zg(mf+1:); xn(:,mf+2:) = xs(:,mf+1:)
         call set_grid(zn, xn)
      end if
      call finalize()
   end subroutine set_fixed_temperature

   !> FreeFlame auto 'check_width': flame touching the domain edges
   logical function domain_too_narrow()
      real(WP) :: mref, ml, mr
      mref = (xs(iT,np) - xs(iT,1))/(zg(np) - zg(1))
      ml = (xs(iT,2) - xs(iT,1))/(zg(2) - zg(1))/mref
      mr = (xs(iT,np-2) - xs(iT,np))/(zg(np-2) - zg(np))/mref
      domain_too_narrow = ml > 0.02_WP .or. mr > 0.02_WP
   end function domain_too_narrow

   ! ============================================================================================ solve
   !> Sim1D::solve. status 0: converged; 1: failed; 2: domain too narrow (auto mode)
   integer function sim_solve(do_refine)
      logical, intent(in) :: do_refine
      real(WP) :: dt, ch, chp, relax, frac
      real(WP), allocatable :: x1(:,:), xo(:,:), zo(:), so(:,:), dprev(:,:)
      integer :: newpts, st, it, j, n, ilo
      sim_solve = 1
      n_tsteps = 0; dt = tstep0
      call finalize()
      call deferred_source()
      newpts = 1
      do while (newpts > 0)
         st = solve_grid()
         if (st /= 0) then
            sim_solve = st; return
         end if
         if (conv_scheme == 1) then                          ! deferred correction: update the source until converged
            relax = 1.0_WP; chp = huge(1.0_WP)
            do it = 1, max_dc
               if (allocated(xo)) deallocate(xo, dprev)
               allocate(xo(nv,np), dprev(nv,np)); xo = xs; dprev = dcs
               if (dc_alpha < 1.0_WP) dc_alpha = min(1.0_WP, dc_alpha + 0.25_WP)
               call deferred_source()
               dcs = dprev + relax*(dcs - dprev)              ! under-relaxed when the passes stop contracting
               st = solve_grid()
               if (st /= 0) then
                  sim_solve = st; return
               end if
               ch = wnorm(xo, xs - xo)
               if (verbose > 0) write(*,'(a,i4,a,es10.3,a,f10.4)') '  deferred correction ', it, ': change ', ch, &
                  ' u0[cm/s]=', 100.0_WP*xs(iU,1)
               if (ch <= 1.0_WP .and. dc_alpha >= 1.0_WP) exit
               if (dc_alpha >= 1.0_WP .and. ch > 0.7_WP*chp) relax = max(0.125_WP, 0.5_WP*relax)
               chp = ch
            end do
            if (it > max_dc) return
         end if
         if (do_refine) then
            if (conv_scheme == 1) then                       ! keep the converged source (per unit upwind spacing:
               if (allocated(zo)) deallocate(zo, so)        ! it is a truncation error, proportional to h)
               allocate(zo(np), so(nv,np)); zo = zg; so = 0.0_WP
               do j = 2, np - 1
                  so(:,j) = dcs(:,j)/dz(merge(j - 1, j, xs(iU,j) > 0.0_WP))
               end do
            end if
            newpts = refine()
            if (verbose > 0) write(*,'(a,i0,a,i0)') '  refine: new points ', newpts, ' -> np = ', np
            if (newpts < 0) return
            if (newpts > 0) dt = tstep0
            if (conv_scheme == 1 .and. newpts > 0) then      ! source on the new grid: linear interpolation x new h
               do j = 2, np - 1
                  ilo = min(max(count(zo <= zg(j)), 1), size(zo) - 1)
                  frac = min(max((zg(j) - zo(ilo))/(zo(ilo+1) - zo(ilo)), 0.0_WP), 1.0_WP)
                  do n = 1, nv
                     dcs(n,j) = (so(n,ilo) + frac*(so(n,ilo+1) - so(n,ilo)))*dz(merge(j - 1, j, xs(iU,j) > 0.0_WP))
                  end do
               end do
               dcs(:,1) = 0.0_WP; dcs(:,np) = 0.0_WP
            end if
         else
            newpts = 0
         end if
      end do
      sim_solve = 0

   contains

      !> Steady solution on the current grid: damped Newton, time steps when it fails. 0 ok, 1 failed, 2 too narrow
      integer function solve_grid()
         integer :: m
         logical :: ok, tsok
         solve_grid = 1
         ok = .false.
         do while (.not. ok)
            call set_rdt(0.0_WP); max_age = 20
            if (allocated(x1)) deallocate(x1)
            allocate(x1(nv,np))                              ! np changes with refinement
            m = newton_solve(xs, x1)
            if (verbose > 0) write(*,'(a,i5,a,l2,a,i3,a,f10.4,a,i0)') '  newton: np=', np, ' energy=', energy, &
               ' m=', m, ' u0[cm/s]=', 100.0_WP*merge(x1(iU,1), xs(iU,1), m >= 0), ' jac=', n_jac
            if (m >= 0) then
               xs = x1; ok = .true.
               if (width_check) then
                  if (domain_too_narrow()) then
                     solve_grid = 2; return
                  end if
               end if
            else
               dt = time_step(nsteps_ts, dt, tsok)
               if (verbose > 0) write(*,'(a,es10.3,a,l2,a,i0)') '  timesteps: dt=', dt, ' ok=', tsok, ' total=', n_tsteps
               if (.not. tsok) return
               dt = min(dt, tmax)
            end if
         end do
         solve_grid = 0
      end function solve_grid

   end function sim_solve

   subroutine set_energy(on)
      logical, intent(in) :: on
      if (on .neqv. energy) then
         energy = on; jac_ok = .false.
      end if
   end subroutine set_energy

   !> FreeFlame.set_initial_guess (no data): linear profiles between locs 0.3 and 0.5 of the domain
   subroutine initial_guess_auto()
      real(WP) :: locs(4), frac, u0, u1, rho_in, rho_eq, Tmid, t
      real(WP) :: x1(nv,np)
      integer :: j, i
      locs = [0.0_WP, 0.3_WP, 0.5_WP, 1.0_WP]
      rho_in = Pres/(Rcst*Tin*sum(Yin/W_sp)); rho_eq = Pres/(Rcst*Teq*sum(Yeq/W_sp))
      u0 = 1.0_WP; u1 = rho_in*1.0_WP/rho_eq
      do j = 1, np
         frac = (zg(j) - zg(1))/(zg(np) - zg(1))
         t = min(max((frac - locs(2))/(locs(3) - locs(2)), 0.0_WP), 1.0_WP)
         x1(iU,j) = u0 + t*(u1 - u0)
         x1(iT,j) = Tin + t*(Teq - Tin)
         x1(iY:,j) = Yin + t*(Yeq - Yin)
      end do
      xs = x1
      zfixed = Undef; tfixed = -1.0_WP
      Tmid = 0.75_WP*Tin + 0.25_WP*Teq
      i = 1
      do j = 1, np
         if (xs(iT,j) < Tmid) i = j
      end do
      if (Tmid - xs(iT,i) < 0.2_WP*(Tmid - Tin)) then
         call set_fixed_temperature(xs(iT,i))
      else if (xs(iT,i+1) - Tmid < 0.2_WP*(Teq - Tmid)) then
         call set_fixed_temperature(xs(iT,i+1))
      else
         call set_fixed_temperature(Tmid)
      end if
      jac_ok = .false.
   end subroutine initial_guess_auto

   !> Adiabatic constant-pressure equilibrium of (Tin, Pres, Yin) -> (Teq, Yeq), replacing Cantera equilibrate('HP').
   !> Element-potential method (as Cantera ChemEquil / NASA CEA): n_k = N exp(-g_k/RT - ln(P/Pref) + sum_j a_jk lam_j),
   !> unknowns lam_j (present elements), ln N, ln T; equations: element balance, sum x_k = 1, h = h(Tin).
   !> Start: complete-combustion products (T from h(T) = h0), element potentials from a basis of major species, a
   !> fixed-T pass, then the full system. Newton with analytic Jacobian, step limits and backtracking.
   !> info: 0 converged, 1 no convergence, 2 no estimate (missing CO2/H2O/CO/H2/O2 or too rich). iters: Newton iterations.
   subroutine equilibrate_hp(info, iters)
      integer, intent(out) :: info, iters
      integer, parameter :: maxit = 200
      character(len=2), dimension(nA) :: anames
      integer, dimension(nA,nS) :: comp
      integer :: el(nA), ne, kb, k, j, it, pass, nb, nvar, ipv(nA+2), linfo, ls, basis(nA)
      logical :: on(nS)
      real(WP) :: a(nA,nS), b(nA), ncc(nS), hort(nS), cpor(nS), sor(nS), nk(nS), xk(nS)
      real(WP) :: lam(nA), lnN, lnTe, h0, Tg, fac, fnorm, fnew, lnP, d, v(nA), bm(nA,nA)
      real(WP) :: Jm(nA+2,nA+2), F(nA+2), dx(nA+2), x0(nA+2)
      external :: dgesv
      info = 2; iters = 0; Teq = -1.0_WP; Yeq = 0.0_WP
      call fcmech_get_atomnames(anames)
      call fcmech_get_composition(comp)
      b = matmul(real(comp, WP), Yin/W_sp)                   ! element moles per kg
      ne = 0
      do j = 1, nA
         if (b(j) > 1.0e-30_WP*sum(b)) then
            ne = ne + 1; el(ne) = j
         end if
      end do
      do k = 1, nS                                           ! species of absent elements stay zero
         on(k) = all(comp(:,k) == 0 .or. b > 1.0e-30_WP*sum(b))
         a(1:ne,k) = real(comp(el(1:ne),k), WP)
      end do
      b(1:ne) = b(el(1:ne))
      lnP = log(Pres/Pref)
      call fcmech_nasa7(Tin, hort, cpor, sor)
      h0 = sum(Yin/W_sp*hort)*Tin                            ! h(Tin)/R, K mol/kg

      ! complete combustion: C -> CO2 (CO), H -> H2O (H2), excess O -> O2, other elements -> most stable pure form
      ncc = 0.0_WP
      call complete_combustion(ls)
      if (ls /= 0) return
      Tg = Tin                                               ! frozen products: h(Tg) = h0
      do it = 1, 50
         call fcmech_nasa7(Tg, hort, cpor, sor)
         d = (sum(ncc*hort)*Tg - h0)/sum(ncc*cpor)
         Tg = min(max(Tg - d, Tin), 6000.0_WP)
         if (abs(d) < 1.0e-3_WP) exit
      end do
      Tg = min(Tg, 3500.0_WP)

      ! element potentials from the basis species (largest estimates first, linearly independent)
      call fcmech_nasa7(Tg, hort, cpor, sor)
      xk = ncc
      do k = 1, nS                                           ! small dissociation seeds so a basis exists at phi = 1
         if (on(k) .and. ncc(k) == 0.0_WP .and. (is(k, 'O', 2, 0, 0) .or. is(k, 'C', 1, 'O', 1) .or. is(k, 'H', 2, 0, 0))) &
            xk(k) = 1.0e-3_WP*sum(ncc)
      end do
      nb = 0
      do while (nb < ne)
         kb = 0
         do k = 1, nS
            if (.not. on(k) .or. xk(k) < 0.0_WP) cycle
            if (kb == 0) then
               kb = k
            else if (xk(k) > xk(kb) .or. (xk(k) == xk(kb) .and. hort(k) - sor(k) < hort(kb) - sor(kb))) then
               kb = k
            end if
         end do
         if (kb == 0) return
         v = a(:,kb)                                         ! Gram-Schmidt against the chosen basis
         do j = 1, nb
            v(1:ne) = v(1:ne) - dot_product(v(1:ne), bm(1:ne,j))*bm(1:ne,j)
         end do
         if (sqrt(sum(v(1:ne)**2)) > 1.0e-8_WP) then
            nb = nb + 1; basis(nb) = kb; bm(1:ne,nb) = v(1:ne)/sqrt(sum(v(1:ne)**2))
         end if
         xk(kb) = -1.0_WP                                    ! tried
      end do
      xk = ncc
      do k = 1, nS
         if (on(k) .and. ncc(k) == 0.0_WP .and. (is(k, 'O', 2, 0, 0) .or. is(k, 'C', 1, 'O', 1) .or. is(k, 'H', 2, 0, 0))) &
            xk(k) = 1.0e-3_WP*sum(ncc)
      end do
      do j = 1, ne
         bm(j,1:ne) = a(1:ne,basis(j))
         lam(j) = hort(basis(j)) - sor(basis(j)) + lnP + log(max(xk(basis(j))/sum(xk), 1.0e-8_WP))
      end do
      call dgesv(ne, 1, bm, nA, ipv, lam, nA, linfo)
      if (linfo /= 0) return
      lnN = log(sum(ncc)); lnTe = log(Tg)

      ! Newton: pass 1 at fixed T (ne+1 unknowns), pass 2 adiabatic (ne+2)
      info = 1
      do pass = 1, 2
         nvar = ne + pass
         call eval(F, .true.)
         fnorm = sqrt(sum(F(1:nvar)**2))
         do it = 1, maxit
            iters = iters + 1
            dx(1:nvar) = -F(1:nvar)
            call dgesv(nvar, 1, Jm, nA+2, ipv, dx, nA+2, linfo)
            if (linfo /= 0) return
            fac = 1.0_WP                                     ! step limits: |d lam| <= 2, |d lnN| <= 0.5, |d lnT| <= 0.1
            if (maxval(abs(dx(1:ne))) > 2.0_WP) fac = 2.0_WP/maxval(abs(dx(1:ne)))
            if (abs(dx(ne+1)) > 0.5_WP) fac = min(fac, 0.5_WP/abs(dx(ne+1)))
            if (pass == 2) then
               if (abs(dx(ne+2)) > 0.1_WP) fac = min(fac, 0.1_WP/abs(dx(ne+2)))
            end if
            x0(1:ne) = lam(1:ne); x0(ne+1) = lnN; x0(ne+2) = lnTe
            do ls = 1, 10                                    ! backtrack until the residual decreases
               lam(1:ne) = x0(1:ne) + fac*dx(1:ne); lnN = x0(ne+1) + fac*dx(ne+1)
               if (pass == 2) lnTe = min(max(x0(ne+2) + fac*dx(ne+2), log(200.0_WP)), log(6000.0_WP))
               call eval(F, .false.)
               fnew = sqrt(sum(F(1:nvar)**2))
               if (fnew < fnorm .or. fnew < 1.0e-12_WP) exit
               fac = 0.5_WP*fac
            end do
            fnorm = fnew
            call eval(F, .true.)
            if (maxval(abs(F(1:nvar))) < 1.0e-12_WP .or. maxval(abs(fac*dx(1:nvar))) < 1.0e-14_WP) exit
         end do
         if (maxval(abs(F(1:nvar))) > 1.0e-9_WP) return
      end do
      Teq = exp(lnTe)
      Yeq = nk*W_sp/sum(nk*W_sp)
      info = 0

   contains

      !> residuals F (element balance relative to b_j, sum x_k - 1, (h - h0)/(R T N)) and, if jac, the Jacobian Jm
      subroutine eval(F, jac)
         real(WP), intent(out) :: F(nA+2)
         logical, intent(in) :: jac
         real(WP) :: TT, hk_(nS)
         integer :: i, m
         TT = exp(lnTe)
         call fcmech_nasa7(TT, hort, cpor, sor)
         nk = 0.0_WP; xk = 0.0_WP
         do i = 1, nS
            if (on(i)) xk(i) = exp(min(-(hort(i) - sor(i)) - lnP + dot_product(a(1:ne,i), lam(1:ne)), 300.0_WP))
         end do
         nk = exp(lnN)*xk
         do m = 1, ne
            F(m) = sum(a(m,:)*nk)/b(m) - 1.0_WP
         end do
         F(ne+1) = sum(xk) - 1.0_WP
         F(ne+2) = sum(xk*hort) - h0/TT/exp(lnN)
         if (.not. jac) return
         hk_ = hort
         Jm = 0.0_WP
         do m = 1, ne
            do i = 1, ne
               Jm(m,i) = sum(a(m,:)*a(i,:)*nk)/b(m)
            end do
            Jm(m,ne+1) = sum(a(m,:)*nk)/b(m)
            Jm(m,ne+2) = sum(a(m,:)*nk*hk_)/b(m)
            Jm(ne+1,m) = sum(a(m,:)*xk)
            Jm(ne+2,m) = sum(a(m,:)*nk*hk_)/exp(lnN)
         end do
         Jm(ne+1,ne+1) = 0.0_WP
         Jm(ne+1,ne+2) = sum(xk*hk_)
         Jm(ne+2,ne+1) = h0/TT/exp(lnN)
         Jm(ne+2,ne+2) = (sum(nk*(hk_**2 + cpor)) - sum(nk*hk_))/exp(lnN) + h0/TT/exp(lnN)
      end subroutine eval

      !> true if species k is made of exactly n1 atoms of e1 (and n2 of e2; e2 = 0 for a pure species)
      logical function is(k, e1, n1, e2, n2)
         integer, intent(in) :: k, n1, n2
         character(len=*), intent(in) :: e1
         class(*), intent(in) :: e2
         integer :: q, want
         is = .true.
         do q = 1, nA
            want = 0
            if (trim(anames(q)) == e1) want = n1
            select type (e2)
            type is (character(*))
               if (trim(anames(q)) == e2) want = n2
            end select
            if (comp(q,k) /= want) is = .false.
         end do
      end function is

      !> complete-combustion products ncc (mol/kg) from the element moles; ls /= 0 if not possible
      subroutine complete_combustion(ls)
         integer, intent(out) :: ls
         integer :: q, i, kco2, kco, kh2o, kh2, ko2, kp
         real(WP) :: nC, nH, nO, dO, c
         ls = 1
         nC = 0.0_WP; nH = 0.0_WP; nO = 0.0_WP
         kco2 = find('C', 1, 'O', 2); kco = find('C', 1, 'O', 1); kh2o = find('H', 2, 'O', 1)
         kh2 = find('H', 2, 0, 0); ko2 = find('O', 2, 0, 0)
         do i = 1, ne
            q = el(i)
            select case (trim(anames(q)))
            case ('C'); nC = b(i)
            case ('H'); nH = b(i)
            case ('O'); nO = b(i)
            case default                                     ! N -> N2, AR -> AR, ...: lowest g/RT per atom at 2000 K
               call fcmech_nasa7(2000.0_WP, hort, cpor, sor)
               kp = 0
               do k = 1, nS
                  if (.not. on(k) .or. comp(q,k) == 0 .or. sum(comp(:,k)) /= comp(q,k)) cycle
                  if (kp == 0) then
                     kp = k
                  else if ((hort(k) - sor(k))/comp(q,k) < (hort(kp) - sor(kp))/comp(q,kp)) then
                     kp = k
                  end if
               end do
               if (kp == 0) return
               ncc(kp) = b(i)/comp(q,kp)
            end select
         end do
         if ((nC > 0.0_WP .and. kco2 == 0) .or. (nH > 0.0_WP .and. kh2o == 0) .or. nO <= 0.0_WP) return
         if (kco2 > 0) ncc(kco2) = ncc(kco2) + nC
         if (kh2o > 0) ncc(kh2o) = ncc(kh2o) + 0.5_WP*nH
         dO = nO - 2.0_WP*nC - 0.5_WP*nH                     ! excess (> 0) or missing (< 0) O atoms
         if (dO >= 0.0_WP) then
            if (ko2 == 0) return
            ncc(ko2) = ncc(ko2) + 0.5_WP*dO
         else
            c = min(-dO, nC)                                 ! CO2 -> CO first, then H2O -> H2
            if (c > 0.0_WP) then
               if (kco == 0) return
               ncc(kco2) = ncc(kco2) - c; ncc(kco) = ncc(kco) + c; dO = dO + c
            end if
            c = min(-dO, 0.5_WP*nH)
            if (c > 0.0_WP) then
               if (kh2 == 0) return
               ncc(kh2o) = ncc(kh2o) - c; ncc(kh2) = ncc(kh2) + c; dO = dO + c
            end if
            if (dO < -1.0e-12_WP*nO) return                  ! richer than CO + H2 (solid carbon region)
         end if
         ls = 0
      end subroutine complete_combustion

      integer function find(e1, n1, e2, n2)
         character(len=*), intent(in) :: e1
         integer, intent(in) :: n1, n2
         class(*), intent(in) :: e2
         integer :: i
         find = 0
         do i = nS, 1, -1
            if (on(i) .and. is(i, e1, n1, e2, n2)) find = i
         end do
      end function find

   end subroutine equilibrate_hp

   !> One pass of Sim1D.solve(auto=True). 0 ok, 1 failed, 2 too narrow
   integer function auto_once()
      integer :: Ns(4), k, st, j
      real(WP) :: z0, z1
      real(WP), allocatable :: zn(:), xn(:,:)
      logical :: solved
      Ns = [np, 12, 24, 48]
      z0 = zg(1); z1 = zg(np)
      auto_once = 1
      do k = 1, 4
         if (Ns(k) > npmax) return
         if (Ns(k) /= np) then
            allocate(zn(Ns(k)), xn(nv,Ns(k)))
            zn = [(z0 + (z1 - z0)*real(j-1, WP)/real(Ns(k)-1, WP), j = 1, Ns(k))]
            xn = 0.0_WP
            call set_grid(zn, xn)
            deallocate(zn, xn)
         end if
         call initial_guess_auto()
         call set_energy(.true.)
         st = sim_solve(.false.)
         if (st == 2) then
            auto_once = 2; return
         end if
         solved = st == 0
         if (.not. solved) then
            call set_energy(.false.)
            st = sim_solve(.false.)
            if (st == 2) then
               auto_once = 2; return
            end if
            if (st == 0) then
               call set_energy(.true.)
               st = sim_solve(.false.)
               if (st == 2) then
                  auto_once = 2; return
               end if
            end if
            solved = st == 0
         end if
         if (solved) then
            st = sim_solve(.true.)
            if (st == 2) then
               auto_once = 2; return
            end if
            if (st == 0) then
               auto_once = 0; return
            end if
         end if
      end do
   end function auto_once

   !> FreeFlame.solve(auto=True): retry with a doubled domain if the flame touches the edges
   !> Second-order convection: the coarse cold-start stages use upwind (robust), then the flame is re-solved and refined.
   integer function solve_auto(width)
      real(WP), intent(in) :: width
      real(WP) :: z8(8), x8(nv,8)
      integer :: it, st, scheme
      scheme = conv_scheme; conv_scheme = 0
      z8 = [0.0_WP, 0.2_WP, 0.3_WP, 0.4_WP, 0.5_WP, 0.6_WP, 0.8_WP, 1.0_WP]*width
      x8 = 0.0_WP
      call set_grid(z8, x8)
      width_check = .true.
      solve_auto = 1
      do it = 1, 12
         st = auto_once()
         if (st == 2) then
            zg = 2.0_WP*zg; dz = zg(2:np) - zg(1:np-1)
            st = refine()
            if (st < 0) exit
            cycle
         end if
         solve_auto = st
         exit
      end do
      width_check = .false.
      conv_scheme = scheme
      if (solve_auto == 0 .and. scheme == 1) then
         jac_ok = .false.; dc_alpha = 0.25_WP                ! switch from the upwind flame: ramp the source
         solve_auto = sim_solve(.true.)
         dc_alpha = 1.0_WP
      end if
   end function solve_auto

   !> set_initial_guess(data=profile) + solve(refine_grid=True)
   integer function solve_from_profile(zp, xp)
      real(WP), intent(in) :: zp(:), xp(:,:)
      real(WP) :: Tmid
      integer :: i, j, scheme
      call set_grid(zp, xp)
      zfixed = Undef; tfixed = -1.0_WP
      Tmid = 0.75_WP*xs(iT,1) + 0.25_WP*xs(iT,np)
      i = 1
      do j = 1, np
         if (xs(iT,j) < Tmid) i = j
      end do
      call set_fixed_temperature(xs(iT,i))
      call set_energy(.true.); jac_ok = .false.
      width_check = .false.
      if (conv_scheme == 1) dc_alpha = 0.0_WP              ! second order: upwind first, then ramp the source in
      solve_from_profile = sim_solve(.true.)
      dc_alpha = 1.0_WP
      if (solve_from_profile /= 0 .and. conv_scheme == 1) then   ! second order failed: upwind from the profile first
         call set_grid(zp, xp)
         zfixed = Undef; tfixed = -1.0_WP
         call set_fixed_temperature(xs(iT,i))
         jac_ok = .false.
         scheme = conv_scheme; conv_scheme = 0
         solve_from_profile = sim_solve(.true.)
         conv_scheme = scheme; jac_ok = .false.; dc_alpha = 0.25_WP
         if (solve_from_profile == 0) solve_from_profile = sim_solve(.true.)
         dc_alpha = 1.0_WP
      end if
   end function solve_from_profile

end module flame1d


program chem_reactor_1D_lfs
   use precision, only: WP
   use string, only: str_medium, str_long
   use param, only: param_init, param_final, param_read, param_exists, param_getsize
   use parallel, only: parallel_init, parallel_final, amRoot
   use messager, only: messager_init, messager_final, die
   use fcmech
   use flame1d
   implicit none
   character(len=str_medium) :: conv_name
   character(len=str_long) :: profile_file, result_file, solution_file
   character(len=str_medium), dimension(nS) :: names
   real(WP) :: width, cpu0, cpu1, eq_cpu
   real(WP), allocatable :: zp(:), xp(:,:)
   integer :: st, iunit, j, npr, k, eq_info, eq_iters

   call parallel_init
   call messager_init
   call param_init
   call fcmech_get_speciesnames(names)
   jac_atol = sqrt(epsilon(1.0_WP))

   ! Unburnt gas: same keywords and mixture rule as the 0D reactor (classical phi, as Cantera set_equivalence_ratio)
   call param_read('Pressure', Pres)
   call param_read('Temperature', Tin)
   call read_composition(Yin)
   call cpu_time(cpu0)
   call equilibrate_hp(eq_info, eq_iters)                   ! burnt state for cold starts (as Cantera equilibrate('HP'))
   call cpu_time(cpu1)
   eq_cpu = cpu1 - cpu0
   call param_read('Refine ratio', ratio, default=3.0_WP)
   call param_read('Refine slope', gslope, default=0.08_WP)
   call param_read('Refine curve', gcurve, default=0.08_WP)
   call param_read('Refine prune', gprune, default=0.05_WP)
   call param_read('Max grid points', npmax, default=400)
   call param_read('T upper bound', Thigh, default=4000.0_WP)
   call param_read('Width', width, default=0.04_WP)
   call param_read('Result file', result_file, default='lfs.res')
   call param_read('Solution file', solution_file, default='flame.dat')
   call param_read('Verbose', verbose, default=0)
   call param_read('Convection scheme', conv_name, default='upwind')
   select case (trim(conv_name))
   case ('upwind'); conv_scheme = 0
   case ('second-order'); conv_scheme = 1
   case default; call die('[chem_reactor_1D_lfs] Convection scheme: upwind or second-order')
   end select

   call cpu_time(cpu0)
   if (param_exists('Initial profile')) then
      call param_read('Initial profile', profile_file)       ! rows: z u T Y_1..Y_nS (mechanism order)
      open(newunit=iunit, file=trim(profile_file), status='old', action='read')
      read(iunit,*) npr
      allocate(zp(npr), xp(nv,npr))
      do j = 1, npr
         read(iunit,*) zp(j), xp(:,j)
      end do
      close(iunit)
      st = solve_from_profile(zp, xp)
   else
      if (eq_info /= 0) call die('[chem_reactor_1D_lfs] equilibrium (adiabatic flame state) not found: no cold start')
      st = solve_auto(width)
   end if
   call cpu_time(cpu1)

   if (amRoot) then
      open(newunit=iunit, file=trim(result_file), status='replace', action='write')
      write(iunit,'(a,1x,i0)') 'status', st
      if (st == 0) then
         write(iunit,'(a,1x,es25.16e3)') 'sl_m_s', xs(iU,1)
      else
         write(iunit,'(a,1x,a)') 'sl_m_s', 'nan'
      end if
      write(iunit,'(a,1x,i0)') 'npts', np
      do k = 1, nS
         write(iunit,'(a,1x,es25.16e3)') 'Y_in_'//trim(names(k)), Yin(k)
      end do
      write(iunit,'(a,1x,es25.16e3)') 'cpu_s', cpu1 - cpu0
      write(iunit,'(a,1x,i0)') 'eq_status', eq_info
      write(iunit,'(a,1x,i0)') 'eq_iterations', eq_iters
      write(iunit,'(a,1x,es25.16e3)') 'eq_cpu_s', eq_cpu
      write(iunit,'(a,1x,es25.16e3)') 'T_eq', Teq
      do k = 1, nS
         write(iunit,'(a,1x,es25.16e3)') 'Y_eq_'//trim(names(k)), Yeq(k)
      end do
      write(iunit,'(a,1x,i0)') 'n_jacobians', n_jac
      write(iunit,'(a,1x,i0)') 'n_newton', n_newton
      write(iunit,'(a,1x,es25.16e3)') 'T_max', maxval(xs(iT,:))
      write(iunit,'(a,1x,es25.16e3)') 'z_fixed', zfixed
      close(iunit)
      if (st == 0) then
         open(newunit=iunit, file=trim(solution_file), status='replace', action='write')
         write(iunit,'(a)', advance='no') 'z u T'
         do k = 1, nS
            write(iunit,'(1x,a)', advance='no') trim(names(k))
         end do
         write(iunit,*)
         do j = 1, np
            write(iunit,'(*(es25.16e3,1x))') zg(j), xs(:,j)
         end do
         close(iunit)
      end if
      write(*,'(a,i0,a,f10.4,a,i0,a,f8.1,a)') '[chem_reactor_1D_lfs] status ', st, ', S_L = ', &
         merge(100.0_WP*xs(iU,1), -1.0_WP, st == 0), ' cm/s, ', np, ' points, cpu ', cpu1 - cpu0, ' s'
   end if
   call param_final
   call messager_final
   call parallel_final

contains

   subroutine read_composition(Yout)
      real(WP), dimension(nS), intent(out) :: Yout
      real(WP), dimension(nS) :: Xin, xf, xo
      real(WP) :: phi, o_f, o_o
      integer :: k
      logical :: anyX, anyY
      Xin = 0.0_WP; anyX = .false.; anyY = .false.
      if (param_exists('Equivalence ratio')) then
         call param_read('Equivalence ratio', phi)
         if (.not. param_exists('Fuel')) call die('[chem_reactor_1D_lfs] Equivalence ratio needs Fuel : <species ...>')
         call read_stream('Fuel', 'Fuel mole fractions', '', xf)
         call read_stream('Oxidizer', 'Oxidizer mole fractions', 'O2 N2', xo)
         o_f = o2_need(xf); o_o = o2_need(xo)
         if (o_f .le. 0.0_WP .or. o_o .ge. 0.0_WP) call die('[chem_reactor_1D_lfs] fuel/oxidizer cannot define phi')
         Xin = phi*(-o_o)*xf + o_f*xo
      else
         do k = 1, nS
            if (param_exists('Initial X '//trim(names(k)))) then
               call param_read('Initial X '//trim(names(k)), Xin(k)); anyX = .true.
            end if
            if (param_exists('Initial Y '//trim(names(k)))) then
               call param_read('Initial Y '//trim(names(k)), Xin(k)); anyY = .true.
            end if
         end do
         if (anyX .and. anyY) call die('[chem_reactor_1D_lfs] use either Initial X or Initial Y lines, not both')
         if (sum(Xin) .le. 0.0_WP) call die('[chem_reactor_1D_lfs] no initial composition given')
         if (anyY) then
            Yout = Xin/sum(Xin); return
         end if
      end if
      Xin = Xin/sum(Xin)
      Yout = Xin*W_sp/sum(Xin*W_sp)
   end subroutine read_composition

   !> Stream composition (normalized mole fractions) from 'tag : A B ..' and 'tagx : xA xB ..'
   subroutine read_stream(tag, tagx, def_names, x)
      character(len=*), intent(in) :: tag, tagx, def_names
      real(WP), dimension(nS), intent(out) :: x
      character(len=str_medium), dimension(nS) :: toks
      character(len=str_long) :: line
      real(WP), allocatable :: v(:)
      integer :: j, n, p, q
      logical :: last
      call param_read(tag, line, default=def_names)
      n = 0; p = 0
      do j = 1, len_trim(line)                                 ! split on blanks (NGA2 param turns commas into blanks)
         if (line(j:j) .eq. ' ') cycle
         if (p .eq. 0) p = j
         last = j .eq. len_trim(line)
         if (.not. last) last = line(j+1:j+1) .eq. ' '
         if (last) then
            n = n + 1
            if (n .gt. nS) call die('[chem_reactor_1D_lfs] too many species in '//trim(tag))
            toks(n) = line(p:j); p = 0
         end if
      end do
      allocate(v(n))
      if (param_exists(tagx)) then
         if (param_getsize(tagx) .ne. n) call die('[chem_reactor_1D_lfs] '//trim(tagx)//': one value per species')
         call param_read(tagx, v)
      else if (n .eq. 1) then
         v = 1.0_WP
      else if (trim(tag) .eq. 'Oxidizer' .and. n .eq. 2 .and. trim(toks(1)) .eq. 'O2' .and. trim(toks(2)) .eq. 'N2') then
         v = [1.0_WP, 3.76_WP]                                 ! default oxidizer: air, O2:N2 = 1:3.76
      else
         call die('[chem_reactor_1D_lfs] missing '//trim(tagx))
      end if
      x = 0.0_WP
      do j = 1, n
         q = species_index(toks(j))
         if (q .le. 0) call die('[chem_reactor_1D_lfs] species not in mechanism: '//trim(toks(j)))
         x(q) = x(q) + v(j)
      end do
      x = x/sum(x)
   end subroutine read_stream

   !> O2 needed for complete oxidation minus O2 contained, per mole of stream (Cantera convention)
   real(WP) function o2_need(x)
      real(WP), dimension(nS), intent(in) :: x
      character(len=2), dimension(nA) :: anames
      integer, dimension(nA,nS) :: comp
      integer :: a
      real(WP) :: w
      call fcmech_get_atomnames(anames)
      call fcmech_get_composition(comp)
      o2_need = 0.0_WP
      do a = 1, nA
         select case (trim(anames(a)))
         case ('C', 'S'); w = 1.0_WP
         case ('H'); w = 0.25_WP
         case ('O'); w = -0.5_WP
         case default; w = 0.0_WP
         end select
         o2_need = o2_need + w*sum(x*real(comp(a,:), WP))
      end do
   end function o2_need

   integer function species_index(name)
      character(len=*), intent(in) :: name
      integer :: k
      species_index = -1
      do k = 1, nS
         if (trim(names(k)) .eq. trim(adjustl(name))) species_index = k
      end do
   end function species_index

end program chem_reactor_1D_lfs
