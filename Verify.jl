#=
Файл с функциями, проверяющими корректность и адекватность численного решения.
Включен в Main.jl через include.

Методы для матриц - двумерная задача, для векторов - одномерная задача вдоль y.
=#

####################################################################################################
# --- двумерная задача ---

function CheckStability(rho1::AbstractMatrix, rho2, u, step, t, E, mu_el, lam_el, phi)
    u_x = u[1]
    u_y = u[2]

    if any(isnan.(rho1)) || any(isnan.(rho2)) || any(isnan.(u_x)) || any(isnan.(u_y))
        @warn "Обнаружен NaN, step = $step, t = $t, rho1 = $(any(isnan.(rho1))), rho2 = $(any(isnan.(rho2))), u_x = $(any(isnan.(u_x))), u_y = $(any(isnan.(u_y)))"
        return false
    end

    if maximum(abs.(rho1)) > 1e6 || maximum(abs.(rho2)) > 1e6
        @warn "Взрыв rho, step = $step, t = $t, max|rho1| = $(maximum(abs.(rho1))), max|rho2| = $(maximum(abs.(rho2)))"
        return false
    end

    if maximum(abs.(u_x)) > 1e4 || maximum(abs.(u_y)) > 1e4
        @warn "Взрыв u, step = $step, t = $t, max|u_x| = $(maximum(abs.(u_x))), max|u_y| = $(maximum(abs.(u_y)))"
        return false
    end
    return true
end

function CheckMassConservation(rho1::AbstractMatrix, rho2, init_mass1, init_mass2, step, t, rtol)
    mass1 = sum(rho1[ficty + 1:Ny + ficty - 1, fictx + 1:Nx + fictx - 1]) * hx * hy
    mass2 = sum(rho2[ficty + 1:Ny + ficty - 1, fictx + 1:Nx + fictx - 1]) * hx * hy
    return CheckMasses(mass1, mass2, init_mass1, init_mass2, step, t, rtol)
end

####################################################################################################
# --- одномерная задача (вдоль y) ---

function CheckStability(rho1::AbstractVector, rho2, u, step, t, E, mu_el, lam_el, phi)
    u_y = u

    if any(isnan.(rho1)) || any(isnan.(rho2)) || any(isnan.(u_y))
        @warn "Обнаружен NaN, step = $step, t = $t, rho1 = $(any(isnan.(rho1))), rho2 = $(any(isnan.(rho2))), u_y = $(any(isnan.(u_y)))"
        return false
    end

    if maximum(abs.(rho1)) > 1e6 || maximum(abs.(rho2)) > 1e6
        @warn "Взрыв rho, step = $step, t = $t, max|rho1| = $(maximum(abs.(rho1))), max|rho2| = $(maximum(abs.(rho2)))"
        return false
    end

    if maximum(abs.(u_y)) > 1e4
        @warn "Взрыв u, step = $step, t = $t, max|u_y| = $(maximum(abs.(u_y)))"
        return false
    end
    return true
end

function CheckMassConservation(rho1::AbstractVector, rho2, init_mass1, init_mass2, step, t, rtol)
    mass1 = sum(rho1[ficty + 1:Ny + ficty - 1]) * hy
    mass2 = sum(rho2[ficty + 1:Ny + ficty - 1]) * hy
    return CheckMasses(mass1, mass2, init_mass1, init_mass2, step, t, rtol)
end

####################################################################################################
# --- общие для обеих размерностей ---

function CheckMasses(mass1, mass2, init_mass1, init_mass2, step, t, rtol)
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
