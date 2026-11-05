#=
Файл, содержащий функции тензоров напряжений, давления и вспомогательных к ним величин.
Включен в Main.jl через include.

Методы для матриц - двумерная задача, для векторов - одномерная задача вдоль y.
ws - рабочие массивы (Workspace.jl), выделенные один раз до цикла по времени.
Двумерные методы параллельны (Polyester, @batch), одномерные - последовательные.
Тензоры в ячейках - Tensor3 (неизменяемые), поэтому тензор ячейки собирается целиком функциями *Cell.
=#

####################################################################################################
# --- тензоры в одной ячейке (общие для обеих размерностей) ---

# --- Pi_ns: диагональный, компоненты xx и yy (в 1D первая из них считается по своей формуле) ---
@inline PiNSCell(p11, p22) = Tensor3(p11, 0.0, 0.0, 0.0, p22, 0.0, 0.0, 0.0, 0.0)

# --- Pi_el: Pi_lk = K_lk - 2 * sum_r E_lr * K_rk ---
@inline function PiElCell(Ei, Ki)
    f(l, k) = Ki[l, k] - 2 * (Ei[l, 1] * Ki[1, k] + Ei[l, 2] * Ki[2, k] + Ei[l, 3] * Ki[3, k])
    return Tensor3(f(1, 1), f(2, 1), f(3, 1),
                   f(1, 2), f(2, 2), f(3, 2),
                   f(1, 3), f(2, 3), f(3, 3))
end

# --- K: K_lk = 2 * mu_el * E_lk + lam_el * tr(E) * delta_lk ---
@inline function KCell(Ei, mui, lami)
    trE = Ei[1,1] + Ei[2,2] + Ei[3,3]
    f(l, k) = 2*mui * Ei[l,k] + (l == k ? lami * trE : 0.0)
    return Tensor3(f(1, 1), f(2, 1), f(3, 1),
                   f(1, 2), f(2, 2), f(3, 2),
                   f(1, 3), f(2, 3), f(3, 3))
end

####################################################################################################
# --- двумерная задача ---

function PiNS!(Pi, U::AbstractVector{<:AbstractMatrix}, ws)
    du_dx = dx!(ws.pins_du_dx, U[1])
    dv_dy = dy!(ws.pins_dv_dy, U[2])
    @batch for idx in eachindex(Pi)
        @inbounds begin
            div_u = du_dx[idx] + dv_dy[idx]
            Pi[idx] = PiNSCell(2 * eta * du_dx[idx] + (zeta - (2 / 3) * eta) * div_u,
                               2 * eta * dv_dy[idx] + (zeta - (2 / 3) * eta) * div_u)
        end
    end
    return Pi
end

function PiEl!(Pi, E::AbstractMatrix{<:AbstractMatrix}, K)
    rows = ficty + 1:Ny + ficty - 1
    @batch for j in fictx + 1:Nx + fictx - 1
        @inbounds for i in rows
            Pi[i, j] = PiElCell(E[i, j], K[i, j])
        end
    end
    return Pi
end

function K!(K_ar, E::AbstractMatrix{<:AbstractMatrix}, mu_el, lam_el)
    rows = ficty + 1:Ny + ficty - 1
    @batch for j in fictx + 1:Nx + fictx - 1
        @inbounds for i in rows
            K_ar[i, j] = KCell(E[i, j], mu_el[i, j], lam_el[i, j])
        end
    end
    return K_ar
end

####################################################################################################
# --- одномерная задача (вдоль y): u = 0, d/dx = 0, поэтому div u = dv/dy ---

function PiNS!(Pi, v::AbstractVector{<:Real}, ws)
    dv_dy = dy!(ws.pins_dv_dy, v)
    @inbounds for idx in eachindex(Pi)
        div_u = dv_dy[idx]
        Pi[idx] = PiNSCell((zeta - (2 / 3) * eta) * div_u,
                           2 * eta * dv_dy[idx] + (zeta - (2 / 3) * eta) * div_u)
    end
    return Pi
end

function PiEl!(Pi, E::AbstractVector{<:AbstractMatrix}, K)
    @inbounds for i in ficty + 1:Ny + ficty - 1
        Pi[i] = PiElCell(E[i], K[i])
    end
    return Pi
end

function K!(K_ar, E::AbstractVector{<:AbstractMatrix}, mu_el, lam_el)
    @inbounds for i in ficty + 1:Ny + ficty - 1
        K_ar[i] = KCell(E[i], mu_el[i], lam_el[i])
    end
    return K_ar
end

####################################################################################################
# --- общие для обеих размерностей ---

function Pressure(E, mu_el, lam_el, rho1, rho2)
	mu1 = Mu1!(zeros(size(rho1)), E, rho1, rho2)
	mu2 = Mu2!(zeros(size(rho1)), E, rho1, rho2)
	psi0 = Psi0(E, mu_el, lam_el, rho1, rho2)
	return rho1 .* mu1 .+ rho2 .* mu2 .- psi0
end
