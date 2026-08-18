# J
J(rho, u) = sStar!(zeros(N + 2 * fict), rho) .* u

# Psi
PsiA(rho1, rho2) = Apsi .* ((rho1 .- rho1A).^2 .+ (rho2 .- rho2A).^2)
PsiB(rho1, rho2) = Bpsi .* ((rho1 .- rho1B).^2 .+ (rho2 .- rho2B).^2)

function Psi0m(rho1, rho2)
    A = PsiA(rho1, rho2)
    B = PsiB(rho1, rho2)
    return (A .* B) ./ (A .+ B)
end

H(rho1, rho2) = 0.5 .+ 0.5 .* tanh.(30 .* ((rho1 ./ (rho1 .+ rho2)) .- 0.5))

I1(matrix) = tr(matrix)
I2(matrix) = 0.5 * ((I1(matrix))^2 - I1(matrix^2))

PsiEl(E, mu_el, lam_el, rho1, rho2) = [
    (mu_el[i] + 0.5 * lam_el[i]) * (I1(E[i]))^2 - 2 * mu_el[i] * I2(E[i])
    for i in eachindex(rho1)
]

Psi0(E, mu_el, lam_el, rho1, rho2) = PsiEl(E, mu_el, lam_el, rho1, rho2) .+ Psi0m(rho1, rho2)

# M
function M!(out, M0_val, rho1, rho2)
    @inbounds for i in eachindex(rho1)
        out[i] = M0_val
    end
    return out
end

# Mu1, Mu2, Mu1Hat, Mu2Hat
function MuEl!(out, rho1, rho2)
    @inbounds for i in eachindex(rho1)
        c = rho1[i] / (rho1[i] + rho2[i])
        out[i] = muel * (0.5 + 0.5 * tanh(30 * (c - 0.5)))
    end
    return out
end
 
function LamEl!(out, rho1, rho2)
    @inbounds for i in eachindex(rho1)
        c = rho1[i] / (rho1[i] + rho2[i])
        out[i] = lamel * (0.5 + 0.5 * tanh(30 * (c - 0.5)))
    end
    return out
end

function Mu1!(out, E, rho1, rho2)
    @inbounds for i in eachindex(rho1)
        r1 = rho1[i]; r2 = rho2[i]
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
    return out
end
 
function Mu2!(out, E, rho1, rho2)
    @inbounds for i in eachindex(rho1)
        r1 = rho1[i]; r2 = rho2[i]
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
    return out
end
 
function LapWall!(out, rho, dStar_buf)
    dStar!(dStar_buf, rho)
    d!(out, dStar_buf)
    if bound == "wall"
        out[fict + 1] = 0.0
        out[N + fict - 1] = 0.0
    end
    return out
end
 
function Mu1Hat!(out, buf, E, rho1, rho2)
    Mu1!(out, E, rho1, rho2)
    LapWall!(buf.lap_rho1, rho1, buf.dStar_rho1)
    LapWall!(buf.lap_rho2, rho2, buf.dStar_rho2)
    @inbounds for i in eachindex(rho1)
        out[i] -= lam11 * buf.lap_rho1[i] + lam12 * buf.lap_rho2[i]
    end
    return out
end
 
function Mu2Hat!(out, buf, E, rho1, rho2)
    Mu2!(out, E, rho1, rho2)
    @inbounds for i in eachindex(rho1)
        out[i] -= lam21 * buf.lap_rho1[i] + lam22 * buf.lap_rho2[i]
    end
    return out
end

function Phi!(out, xc, rho1)
    @inbounds for i in eachindex(rho1)
        out[i] = -g * xc[i]
    end
    return out
end 