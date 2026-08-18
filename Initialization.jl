# INITIAL DENSITIES
function PistonFraction(x_c, xc_piston)
    if bound == "periodic"
        dist = mod(x_c - xc_piston + 0.5 * L, L) - 0.5 * L
    elseif bound == "wall"
        dist = x_c - xc_piston
    end
    if rho_init == "sharp"
        return abs(dist) <= 0.5 * piston_width ? 1.0 : 0.0
    elseif rho_init == "smooth"
        return 0.5 * (1.0 + tanh(b * (0.5 * piston_width + dist))) *
               0.5 * (1.0 + tanh(b * (0.5 * piston_width - dist)))
    end
end
 
function BuildDensity(x_points, xc_piston)
    n = length(x_points)
    rho1 = zeros(n)
    rho2 = zeros(n)
    for i in 1:n
        f = PistonFraction(x_points[i], xc_piston)
        rho1[i] = rho1A + (rho1B - rho1A) * f
        rho2[i] = rho2A + (rho2B - rho2A) * f
    end
    return rho1, rho2
end
 
function InitDensitiesSpinodal(xc)
    rho1 = zeros(N + 2 * fict - 1)
    rho2 = zeros(N + 2 * fict - 1)
 
    rho1_mean = 0.5 * (rho1A + rho1B)
    rho2_mean = 0.5 * (rho2A + rho2B)
 
    amp = 0.02
 
    Random.seed!(123)
 
    for i in fict + 1:N - 1 + fict
        rho1[i] = rho1_mean * (1.0 + amp * (2.0 * rand() - 1.0))
        rho2[i] = 1000 - rho1[i]
    end
 
    rho = rho1 .+ rho2
    return rho1, rho2, rho
 
end
 
function InitDensitiesPistonImitation(xc)
    rho1 = zeros(N + 2 * fict - 1)
    rho2 = zeros(N + 2 * fict - 1)
 
    r1, r2 = BuildDensity(xc[fict + 1:N + fict - 1], x0_piston)
    rho1[fict + 1:N + fict - 1] .= r1
    rho2[fict + 1:N + fict - 1] .= r2
 
    rho = rho1 .+ rho2
    return rho1, rho2, rho
end
 
function InitDensityChoise(xc)
    if rho_init == "sharp" || rho_init == "smooth"
        return InitDensitiesPistonImitation(xc)
    elseif rho_init == "spinodal"
        return InitDensitiesSpinodal(xc)
    end
end

# INITIAL VELOCITY
function InitVelocity()
    u = [zeros(3) for i in 1:(N + 2 * fict)]
    for i in fict + 1:N + fict
        u[i][1] = u_init
    end
    return u
end

# INITIAL ALMANSI TENSOR
function InitAlmansiTens()
    E = [zeros(3,3) for i in 1:(N + 2 * fict - 1)]
    return E
end