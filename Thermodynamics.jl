##########
# J
##########

J(rho, u) = sStar(rho) .* u

##########
# Psi
##########

PsiA(rho1, rho2) = Apsi .* ((rho1 .- rho1A).^2 .+ (rho2 .- rho2A).^2)
PsiB(rho1, rho2) = Bpsi .* ((rho1 .- rho1B).^2 .+ (rho2 .- rho2B).^2)

function Psi0m(rho1, rho2)
    A = PsiA(rho1, rho2)
    B = PsiB(rho1, rho2)
    return (A .* B) ./ (A .+ B)
end

H(rho1, rho2) = 0.5 .+ 0.5 .* tanh.(30 .* ((rho1 ./ (rho1 .+ rho2)) .- 0.5))

MuEl(muel, rho1, rho2) = muel * H(rho1, rho2)
LamEl(lamel, rho1, rho2) = lamel * H(rho1, rho2)

I1(matrix) = tr(matrix)
I2(matrix) = 0.5 * ((I1(matrix))^2 - I1(matrix^2))

PsiEl(E, muel, lamel, rho1, rho2) = [
    (muel[i] + 0.5 * lamel[i]) * (I1(E[i]))^2 - 2 * muel[i] * I2(E[i])
    for i in eachindex(rho1)
]

Psi0(E, muel, lamel, rho1, rho2) = PsiEl(E, muel, lamel, rho1, rho2) .+ Psi0m(rho1, rho2)

##########
# Mu1, Mu2, Mu1Hat, Mu2Hat
##########

function Mu1(E, muel, lamel, rho1, rho2)
    mu1 = zeros(length(rho1))
    for i in eachindex(rho1)

        A = Apsi * ((rho1[i]-rho1A)^2 + (rho2[i]-rho2A)^2)
        B = Bpsi * ((rho1[i]-rho1B)^2 + (rho2[i]-rho2B)^2)
        dA = 2*Apsi*(rho1[i]-rho1A)
        dB = 2*Bpsi*(rho1[i]-rho1B)
        denom = (A + B)^2
        mu1_m = denom > 1e-30 ? (dA*B^2 + dB*A^2) / denom : 0.0


        c_i = rho1[i] / (rho1[i] + rho2[i])
        H_prime = 15.0 / cosh(30.0*(c_i - 0.5))^2       
        dH_drho1 = H_prime * rho2[i] / (rho1[i]+rho2[i])^2 
        E_sq = sum(E[i][j,k]^2 for j in 1:3, k in 1:3)  
        trE = E[i][1,1] + E[i][2,2] + E[i][3,3]
        mu1_el = (muelB * E_sq + 0.5 * lamelB * trE^2) * dH_drho1

        mu1[i] = mu1_m + mu1_el
    end
    return mu1
end

function Mu2(E, muel, lamel, rho1, rho2)
    mu2 = zeros(length(rho1))
    for i in eachindex(rho1)

        A = Apsi * ((rho1[i]-rho1A)^2 + (rho2[i]-rho2A)^2)
        B = Bpsi * ((rho1[i]-rho1B)^2 + (rho2[i]-rho2B)^2)
        dA = 2*Apsi*(rho2[i]-rho2A)
        dB = 2*Bpsi*(rho2[i]-rho2B)
        denom = (A + B)^2
        mu2_m = denom > 1e-30 ? (dA*B^2 + dB*A^2) / denom : 0.0


        c_i = rho1[i] / (rho1[i] + rho2[i])
        H_prime = 15.0 / cosh(30.0*(c_i - 0.5))^2
        dH_drho2 = -H_prime * rho1[i] / (rho1[i]+rho2[i])^2  
        E_sq = sum(E[i][j,k]^2 for j in 1:3, k in 1:3)
        trE = E[i][1,1] + E[i][2,2] + E[i][3,3]
        mu2_el = (muelB * E_sq + 0.5 * lamelB * trE^2) * dH_drho2

        mu2[i] = mu2_m + mu2_el
    end
    return mu2
end

function LapWall(rho)
    lap = d(dStar(rho))
    if bound == "wall"
        lap[fict + 1] = 0.0
        lap[N + fict - 1] = 0.0
    elseif bound == "periodic"
    end
    return lap
end

Mu1Hat(E, muel, lamel, rho1, rho2) = Mu1(E, muel, lamel, rho1, rho2) .- (lam11 .* LapWall(rho1) .+ lam12 .* LapWall(rho2))

Mu2Hat(E, muel, lamel, rho1, rho2) = Mu2(E, muel, lamel, rho1, rho2) .- (lam21 .* LapWall(rho1) .+ lam22 .* LapWall(rho2))


##########
# M
##########

M(M0, rho1, rho2) = [M0 for i in eachindex(rho1)]
