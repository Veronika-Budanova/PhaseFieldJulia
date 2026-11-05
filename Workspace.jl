#=
Файл, содержащий функции выделения рабочих массивов для функций, вызываемых на каждом шаге по времени
(Solver.jl, PiNS!, Mu1Hat!, Mu2Hat!). Массивы выделяются один раз до цикла по времени и перезаписываются.
Включен в Main.jl через include.

Формы массивов:
    c  - центры ячеек
    fx - грани по x (2D), fy - грани по y
    n  - узлы (2D)
Массивы tmp_* - временные, их содержимое используется только внутри одной функции.
=#

function Workspace2D()
    sz_c  = (Ny + 2 * ficty - 1, Nx + 2 * fictx - 1)
    sz_fx = (Ny + 2 * ficty - 1, Nx + 2 * fictx)
    sz_fy = (Ny + 2 * ficty,     Nx + 2 * fictx - 1)
    sz_n  = (Ny + 2 * ficty,     Nx + 2 * fictx)

    c()  = zeros(sz_c)
    fx() = zeros(sz_fx)
    fy() = zeros(sz_fy)
    n()  = zeros(sz_n)

    return (
        # --- ConservMass! ---
        flux1_x = fx(), flux1_y = fy(), flux2_x = fx(), flux2_y = fy(),
        d_rho1_u = c(), d_rho2_u = c(), theta_mu_diff = c(),
        flux_M_x = fx(), flux_M_y = fy(), d_M_Mu = c(),
        # --- ConservMomentum! ---
        rho = c(), sStarx_rho = fx(), sStary_rho = fy(), Jx = fx(), Jy = fy(),
        fluxJ_xx = c(), fluxJ_yy = c(), fluxJ_yx = n(), fluxJ_xy = n(),
        d_J_u = fx(), d_J_v = fy(), rho_Mu_sum_x = fx(), rho_Mu_sum_y = fy(),
        d_E_K_x = fx(), d_E_K_y = fy(), E_rq = c(), K_rq = c(),
        Pi_xx = c(), Pi_yy = c(), Pi_el_xy = c(), Pi_xy = n(),
        d_Pi_x = fx(), d_Pi_y = fy(), force_x = fx(), force_y = fy(),
        rho_new = c(), denom_x = fx(), denom_y = fy(),
        # --- ConservAlmansiTensor! ---
        E_xx = c(), E_yy = c(), E_zz = c(), E_xy = c(),
        du_dx = c(), dv_dy = c(), dv_dx = c(), du_dy = c(),
        conv_xx = c(), conv_yy = c(), conv_zz = c(), conv_xy = c(),
        # --- PiNS! ---
        pins_du_dx = c(), pins_dv_dy = c(),
        # --- Mu1Hat!, Mu2Hat! (LapWall!) ---
        lap1 = c(), lap2 = c(), lap_c = c(), lap_fx = fx(), lap_fy = fy(),
        # --- временные ---
        tmp_c = c(), tmp_fx = fx(), tmp_fx2 = fx(), tmp_fy = fy(), tmp_fy2 = fy(),
        tmp_n = n(), tmp_n2 = n(),
    )
end

function Workspace1D()
    c() = zeros(Ny + 2 * ficty - 1)
    f() = zeros(Ny + 2 * ficty)

    return (
        # --- ConservMass! ---
        flux1_y = f(), flux2_y = f(),
        d_rho1_u = c(), d_rho2_u = c(), theta_mu_diff = c(),
        flux_M_y = f(), d_M_Mu = c(),
        # --- ConservMomentum! ---
        rho = c(), sStary_rho = f(), Jy = f(),
        fluxJ_yy = c(), d_J_v = f(), rho_Mu_sum_y = f(),
        d_E_K_y = f(), E_rq = c(), K_rq = c(),
        Pi_yy = c(), d_Pi_y = f(), force_y = f(),
        rho_new = c(), denom_y = f(),
        # --- ConservAlmansiTensor! ---
        E_xx = c(), E_yy = c(), E_zz = c(), E_xy = c(), dv_dy = c(),
        conv_xx = c(), conv_yy = c(), conv_zz = c(), conv_xy = c(),
        # --- PiNS! ---
        pins_dv_dy = c(),
        # --- Mu1Hat!, Mu2Hat! (LapWall!) ---
        lap1 = c(), lap2 = c(), lap_fy = f(),
        # --- временные ---
        tmp_c = c(), tmp_fy = f(), tmp_fy2 = f(),
    )
end

# --- выбор по размерности поля (для вызовов вне цикла по времени, например при записи фрейма) ---
Workspace(rho::AbstractMatrix) = Workspace2D()
Workspace(rho::AbstractVector) = Workspace1D()
