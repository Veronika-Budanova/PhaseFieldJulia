#=
Файл со вспомогательными функциями, не связанными непосредственно с численным решением.
=#

# --- загрузка конфигурационного файла ---
function DownloadConfig()
    global config_file = ARGS[1] 
    println("Загружаю конфигурацию: $config_file")
    include(config_file)
end

# --- расчет и сохранение результатов одного фрейма ---
function MakeFrame(data_io, xc, x, rho1, rho2, rho, u, E, phi,
                   t, step, frame_count, energy_hist, time_hist, frames_data)
    mu_el  = MuEl!(zeros(length(rho1)), rho1, rho2)
    lam_el = LamEl!(zeros(length(rho1)), rho1, rho2)
    p      = Pressure(E, mu_el, lam_el, rho1, rho2)

    save_state_frame!(xc, x, rho1, rho2, rho, u, p, E, phi, fict,
                      t, step, frame_count, energy_hist, time_hist, frames_data)
    serialize(data_io, frames_data[end])
    flush(data_io)

    return frame_count + 1, mu_el, lam_el
end

# --- создание временных массивов для хранения данных ---
function CreateBuffers()
    N_centr = N + 2 * fict - 1
    N_face  = N + 2 * fict
    return (
        flux1            = zeros(N_face),
        flux2            = zeros(N_face),
        sStar_rho1       = zeros(N_face),
        sStar_rho2       = zeros(N_face),
        d_rho1_u         = zeros(N_centr),
        d_rho2_u         = zeros(N_centr),
        sStar_M          = zeros(N_face),
        theta_mu_diff    = zeros(N_centr),
        dStar_theta_mu   = zeros(N_face),
        flux_M           = zeros(N_face),
        d_M_Mu           = zeros(N_centr),
        u_x              = zeros(N_face),
        Jx               = zeros(N_face),
        s_J              = zeros(N_centr),
        s_u_x            = zeros(N_centr),
        fluxJ            = zeros(N_centr),
        d_J_u            = zeros(N_face),
        sStar_rho        = zeros(N_face),
        dStar_mu1        = zeros(N_face),
        dStar_mu2        = zeros(N_face),
        rho_Mu_sum       = zeros(N_face),
        d_E_K_x          = zeros(N_face),
        E_rq             = zeros(N_centr),
        K_rq             = zeros(N_centr),
        dStar_E_rq       = zeros(N_face),
        sStar_K_rq       = zeros(N_face),
        Pi_diag          = zeros(N_centr),
        d_Pi_x           = zeros(N_face),
        dStar_Phi        = zeros(N_face),
        force            = zeros(N_face),
        rho_new_sum      = zeros(N_centr),
        denom            = zeros(N_face),
        u_x_new          = zeros(N_face),
        rho              = zeros(N_centr),
        E_jj             = zeros(N_centr),
        dstar_E_jj       = zeros(N_face),
        u_x_dstar_E      = zeros(N_face),
        conv             = zeros(N_centr),
        du_dx            = zeros(N_centr),
        mu_el            = zeros(N_centr),
        lam_el           = zeros(N_centr),
        H_arr            = zeros(N_centr),
        mu1              = zeros(N_centr),
        mu2              = zeros(N_centr),
        mu1h             = zeros(N_centr),
        mu2h             = zeros(N_centr),
        M_arr            = zeros(N_centr),
        lap_rho1         = zeros(N_centr),
        lap_rho2         = zeros(N_centr),
        dStar_rho1       = zeros(N_face),
        dStar_rho2       = zeros(N_face),
        K_current        = [zeros(3, 3) for _ in 1:N_centr],
        Pi_el_current    = [zeros(3, 3) for _ in 1:N_centr],
        Pi_ns_current    = [zeros(3, 3) for _ in 1:N_centr],
    )
end

# --- проверка на NaN и расхождение решений
function CheckStability(rho1, rho2, u, step, t, E, mu_el, lam_el, phi)
    u_x = [u[i][1] for i in eachindex(u)]

    if any(isnan.(rho1)) || any(isnan.(rho2)) || any(isnan.(u_x))
        @warn "Обнаружен NaN, step = $step, t = $t, rho1 = $(any(isnan.(rho1))), rho2 = $(any(isnan.(rho2))), u_x = $(any(isnan.(u_x)))"
        return false
    end

    if maximum(abs.(rho1)) > 1e6 || maximum(abs.(rho2)) > 1e6
        @warn "Взрыв rho, step = $step, t = $t, max|rho1| = $(maximum(abs.(rho1))), max|rho2| = $(maximum(abs.(rho2)))"
        return false
    end

    if maximum(abs.(u_x)) > 1e4
        @warn "Взрыв u, step = $step, t = $t, max|u| = $(maximum(abs.(u_x)))"
        return false
    end
    return true
end

# --- проверка сохранения масс компонентов и полной массы ---
function CheckMassConservation(rho1, rho2, init_mass1, init_mass2, step, t, rtol)
    mass1 = sum(rho1[fict + 1:N + fict - 1]) * h
    mass2 = sum(rho2[fict + 1:N + fict - 1]) * h
    init_total = init_mass1 + init_mass2
    total = mass1 + mass2

    ok = true
    if abs(mass1 - init_mass1) > rtol * abs(init_mass1)
        @warn "Нарушено сохранение массы вещества 1, step = $step, t = $t, mass1 = $mass1, init1 = $init_mass1"
        ok = false
    end
    if abs(mass2 - init_mass2) > rtol * abs(init_mass2)
        @warn "Нарушено сохранение массы вещества 2, step = $step, t = $t, mass2 = $mass2, init2 = $init_mass2"
        ok = false
    end
    if abs(total - init_total) > rtol * abs(init_total)
        @warn "Нарушено сохранение полной массы, step = $step, t = $t, total_mass = $total, init_mass = $init_total"
        ok = false
    end
    return ok
end

# --- проверка невозрастания полной энергии ---
function CheckEnergyNonIncrease(energy_hist, step, t, rtol)
    n = length(energy_hist)
    if n < 3
        return true
    end
    prev = energy_hist[n - 1]
    curr = energy_hist[n]
    if curr > prev + rtol * abs(prev)
        @warn "Нарушено невозрастание энергии, step = $step, t = $t, E = $curr, E_prev = $prev"
        return false
    end
    return true
end