#=
Файл, содержащий функции термодинамических переменных, потенциалов и вспомогательных к ним величин.
Включен в Main.jl через include.

Поточечные функции общие для обеих размерностей.
Функции с производными по пространству имеют методы для матриц (2D) и векторов (1D, вдоль y).
ws - рабочие массивы (Workspace.jl); если ws не передан (запись фреймов), он выделяется заново.
Поточечные циклы с @batch2d параллельны при dim == 2 и последовательны при dim == 1 (Parallel.jl).
=#

####################################################################################################
# --- общие для обеих размерностей ---

PsiA(rho1, rho2) = Apsi .* ((rho1 .- rho1A).^2 .+ (rho2 .- rho2A).^2)
PsiB(rho1, rho2) = Bpsi .* ((rho1 .- rho1B).^2 .+ (rho2 .- rho2B).^2)
function Psi0m(rho1, rho2)
    A = PsiA(rho1, rho2)
    B = PsiB(rho1, rho2)
    return (A .* B) ./ (A .+ B)
end
H(rho1, rho2) = 0.5 .+ 0.5 .* tanh.(30 .* ((rho1 ./ (rho1 .+ rho2)) .- 0.5))
I1(matrix) = tr(matrix)
# --- квадрат тензора считается для Matrix{Float64}, как до перехода на Tensor3: так энергия и давление ---
# --- во фреймах совпадают с прежними побитово (у SMatrix свое умножение с другим порядком операций) ---
# --- используется только при записи фреймов, поэтому выделение памяти здесь не влияет на скорость шага ---
I2(matrix) = 0.5 * ((I1(matrix))^2 - I1(Matrix(matrix)^2))
PsiEl(E, mu_el, lam_el, rho1, rho2) = [(mu_el[i] + 0.5 * lam_el[i]) * (I1(E[i]))^2 - 2 * mu_el[i] * I2(E[i]) for i in CartesianIndices(rho1)]
Psi0(E, mu_el, lam_el, rho1, rho2) = PsiEl(E, mu_el, lam_el, rho1, rho2) .+ Psi0m(rho1, rho2)

function M!(out, M0_val, rho1, rho2)
    @batch2d for i in eachindex(rho1)
        @inbounds out[i] = M0_val
    end
    return out
end

function MuEl!(out, rho1, rho2)
    @batch2d for i in eachindex(rho1)
        @inbounds begin
            c = rho1[i] / (rho1[i] + rho2[i])
            out[i] = muel * (0.5 + 0.5 * tanh(30 * (c - 0.5)))
        end
    end
    return out
end
function LamEl!(out, rho1, rho2)
    @batch2d for i in eachindex(rho1)
        @inbounds begin
            c = rho1[i] / (rho1[i] + rho2[i])
            out[i] = lamel * (0.5 + 0.5 * tanh(30 * (c - 0.5)))
        end
    end
    return out
end

function Mu1!(out, E, rho1, rho2)
    @batch2d for i in eachindex(rho1)
        @inbounds begin
            r1 = rho1[i]
            r2 = rho2[i]

            A = Apsi * ((r1 - rho1A)^2 + (r2 - rho2A)^2)
            B = Bpsi * ((r1 - rho1B)^2 + (r2 - rho2B)^2)

            dA = 2 * Apsi * (r1 - rho1A)
            dB = 2 * Bpsi * (r1 - rho1B)

            denom = (A + B)^2
            mu1_m = denom > 1e-30 ? (dA * B^2 + dB * A^2) / denom : 0.0

            c_i = r1 / (r1 + r2)
            H_prime = 15.0 / cosh(30.0 * (c_i - 0.5))^2
            dH_drho1 = H_prime * r2 / (r1 + r2)^2

            Ei = E[i]
            E_sq = Ei[1,1]^2 + Ei[1,2]^2 + Ei[1,3]^2 +
                   Ei[2,1]^2 + Ei[2,2]^2 + Ei[2,3]^2 +
                   Ei[3,1]^2 + Ei[3,2]^2 + Ei[3,3]^2
            trE = Ei[1,1] + Ei[2,2] + Ei[3,3]
            mu1_el = (muel * E_sq + 0.5 * lamel * trE^2) * dH_drho1

            out[i] = mu1_m + mu1_el
        end
    end
    return out
end
function Mu2!(out, E, rho1, rho2)
    @batch2d for i in eachindex(rho1)
        @inbounds begin
            r1 = rho1[i]
            r2 = rho2[i]

            A = Apsi * ((r1 - rho1A)^2 + (r2 - rho2A)^2)
            B = Bpsi * ((r1 - rho1B)^2 + (r2 - rho2B)^2)

            dA = 2 * Apsi * (r2 - rho2A)
            dB = 2 * Bpsi * (r2 - rho2B)

            denom = (A + B)^2
            mu2_m = denom > 1e-30 ? (dA * B^2 + dB * A^2) / denom : 0.0

            c_i = r1 / (r1 + r2)
            H_prime = 15.0 / cosh(30.0 * (c_i - 0.5))^2
            dH_drho2 = -H_prime * r1 / (r1 + r2)^2

            Ei = E[i]
            E_sq = Ei[1,1]^2 + Ei[1,2]^2 + Ei[1,3]^2 +
                   Ei[2,1]^2 + Ei[2,2]^2 + Ei[2,3]^2 +
                   Ei[3,1]^2 + Ei[3,2]^2 + Ei[3,3]^2
            trE = Ei[1,1] + Ei[2,2] + Ei[3,3]
            mu2_el = (muel * E_sq + 0.5 * lamel * trE^2) * dH_drho2

            out[i] = mu2_m + mu2_el
        end
    end
    return out
end

# --- лапласиан посчитан только во внутренних ячейках, фиктивные ячейки mu_hat заполняются граничными условиями ---
function Mu1Hat!(out, E, rho1, rho2, ws = Workspace(rho1))
    Mu1!(out, E, rho1, rho2)
    lap_rho1 = LapWall!(ws.lap1, rho1, ws)
    lap_rho2 = LapWall!(ws.lap2, rho2, ws)
    @batch2d for i in eachindex(rho1)
        @inbounds out[i] -= lam11 * lap_rho1[i] + lam12 * lap_rho2[i]
    end
    BoundCondScalar(out)
    return out
end
function Mu2Hat!(out, E, rho1, rho2, ws = Workspace(rho1))
    Mu2!(out, E, rho1, rho2)
    lap_rho1 = LapWall!(ws.lap1, rho1, ws)
    lap_rho2 = LapWall!(ws.lap2, rho2, ws)
    @batch2d for i in eachindex(rho1)
        @inbounds out[i] -= lam21 * lap_rho1[i] + lam22 * lap_rho2[i]
    end
    BoundCondScalar(out)
    return out
end

####################################################################################################
# --- двумерная задача ---

# --- дискретный лапласиан в центрах ячеек ---
# Условие стенки (нулевая нормальная производная) обеспечивается зеркальными фиктивными ячейками
# из BoundCondScalar, поэтому в приграничных ячейках лапласиан не обнуляется: так mu_hat остается
# вариационной производной дискретной градиентной энергии E_lam из ComputeEnergy.
function LapWall!(lap, rho::AbstractMatrix, ws)
    dx!(lap, dStarx!(ws.lap_fx, rho))
    dy!(ws.lap_c, dStary!(ws.lap_fy, rho))
    Map2D!(+, lap, lap, ws.lap_c)
    return lap
end

function Phi(yc, rho1::AbstractMatrix)
    out = zeros(size(rho1))
    @inbounds for j in axes(rho1, 2), i in axes(rho1, 1)
        out[i, j] = -g * yc[i]
    end
    return out
end

# --- плотность градиентной энергии в центрах ячеек: sx(dStarx a * dStarx b) + sy(dStary a * dStary b) ---
function GradProduct(a::AbstractMatrix, b)
    return sx(dStarx(a) .* dStarx(b)) .+ sy(dStary(a) .* dStary(b))
end

# --- полная энергия на двумерной сетке (интеграл по внутренним ячейкам, элемент площади hx * hy) ---
function ComputeEnergy(rho1::AbstractMatrix, rho2, U, E, Phi)
    I  = ficty + 1:Ny + ficty - 1
    J  = fictx + 1:Nx + fictx - 1
    dS = hx * hy

    mu_el  = MuEl!(zeros(size(rho1)), rho1, rho2)
    lam_el = LamEl!(zeros(size(rho1)), rho1, rho2)

    E_psi = sum(Psi0(E[I, J], mu_el[I, J], lam_el[I, J], rho1[I, J], rho2[I, J])) * dS

    u     = U[1]
    v     = U[2]
    s_u_2 = sx(u .* u) .+ sy(v .* v)
    rho   = rho1 .+ rho2
    E_kinetic = 0.5 * sum(rho[I, J] .* s_u_2[I, J]) * dS

    E_lam = 0.5 * sum(
        lam11 .* GradProduct(rho1, rho1)[I, J] .+
        lam12 .* GradProduct(rho2, rho1)[I, J] .+
        lam21 .* GradProduct(rho1, rho2)[I, J] .+
        lam22 .* GradProduct(rho2, rho2)[I, J]
    ) * dS

    E_force = sum(rho[I, J] .* Phi[I, J]) * dS
    E_total = E_psi + E_kinetic + E_lam - E_force

    return E_total, E_psi, E_kinetic, E_lam, E_force
end

####################################################################################################
# --- одномерная задача (вдоль y) ---

# --- дискретный лапласиан в центрах ячеек (комментарий к двумерному методу выше справедлив и здесь) ---
function LapWall!(lap, rho::AbstractVector, ws)
    dy!(lap, dStary!(ws.lap_fy, rho))
    return lap
end

function Phi(yc, rho1::AbstractVector)
    out = zeros(size(rho1))
    @inbounds for i in eachindex(rho1)
        out[i] = -g * yc[i]
    end
    return out
end

# --- плотность градиентной энергии в центрах ячеек: sy(dStary a * dStary b) ---
function GradProduct(a::AbstractVector, b)
    return sy(dStary(a) .* dStary(b))
end

# --- полная энергия на одномерной сетке (интеграл по внутренним ячейкам, элемент длины hy) ---
function ComputeEnergy(rho1::AbstractVector, rho2, U, E, Phi)
    I  = ficty + 1:Ny + ficty - 1
    dS = hy

    mu_el  = MuEl!(zeros(size(rho1)), rho1, rho2)
    lam_el = LamEl!(zeros(size(rho1)), rho1, rho2)

    E_psi = sum(Psi0(E[I], mu_el[I], lam_el[I], rho1[I], rho2[I])) * dS

    v     = U
    s_u_2 = sy(v .* v)
    rho   = rho1 .+ rho2
    E_kinetic = 0.5 * sum(rho[I] .* s_u_2[I]) * dS

    E_lam = 0.5 * sum(
        lam11 .* GradProduct(rho1, rho1)[I] .+
        lam12 .* GradProduct(rho2, rho1)[I] .+
        lam21 .* GradProduct(rho1, rho2)[I] .+
        lam22 .* GradProduct(rho2, rho2)[I]
    ) * dS

    E_force = sum(rho[I] .* Phi[I]) * dS
    E_total = E_psi + E_kinetic + E_lam - E_force

    return E_total, E_psi, E_kinetic, E_lam, E_force
end
