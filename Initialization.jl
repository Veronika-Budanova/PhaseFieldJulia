##########
# INITIAL DENSITIES
##########

function InitDensitiesSpinodal(xc)
    rho1 = zeros(N + 2 * fict - 1)
    rho2 = zeros(N + 2 * fict - 1)

    rho1_mean = 0.5 * (rho1A + rho1B)
    rho2_mean = 0.5 * (rho2A + rho2B)

    amp = 0.02

    Random.seed!(123)

    for i in N + fict:N - 1 + fict
        rho1[i] = rho1_mean * (1.0 + amp * (2.0 * rand() - 1.0))
        rho2[i] = 1000 - rho1[i]
    end

    rho = rho1 .+ rho2

    return rho1, rho2, rho

end

function InitDensitiesPistonImitationSmooth(xc)
    rho1 = zeros(N + 2 * fict - 1)
    rho2 = zeros(N + 2 * fict - 1)

    for i in fict + 1: N - 1 + fict
        x_c = xc[i]
        if bound == "periodic"
            dist = mod(x_c - x0_piston + 0.5 * L, L) - 0.5 * L
        elseif bound == "wall"
            dist = x_c - x0_piston
        end
        f = 0.5 * (1.0 + tanh(b * (0.5 * piston_width + dist))) *
            0.5 * (1.0 + tanh(b * (0.5 * piston_width - dist)))

        rho1[i] = rho1A + (rho1B - rho1A) * f
        rho2[i] = rho2A + (rho2B - rho2A) * f
    end

    flush(stdout)

    rho = rho1 .+ rho2
    return rho1, rho2, rho
end

function InitDensitiesPistonImitationSharp(xc)
    rho1 = zeros(N + 2 * fict - 1)
    rho2 = zeros(N + 2 * fict - 1)

    for i in fict + 1: N - 1 + fict
        x_c = xc[i]
        if bound == "periodic"
            dist = mod(x_c - x0_piston + 0.5 * L, L) - 0.5 * L
        elseif bound == "wall"
            dist = x_c - x0_piston
        end
        if abs(dist) <= 0.5 * piston_width
            f = 1.0 
        else
            f = 0.0 
        end

        rho1[i] = rho1A + (rho1B - rho1A) * f
        rho2[i] = rho2A + (rho2B - rho2A) * f
    end

    flush(stdout)

    rho = rho1 .+ rho2
    return rho1, rho2, rho
end

##########
# INITIAL VELOCITY
##########

function InitVelocity()
    u = [zeros(3) for i in 1:(N + 2 * fict)]
    for i in fict + 1:N + fict
        u[i][1] = u_init
    end
    return u

end

function InitVelocityPistonImitation(x, rho1, rho2)
    u = [zeros(3) for i in 1:(N + 2 * fict)]
    rho_total = rho1 .+ rho2

    rho_min = minimum(rho_total)
    rho_max = maximum(rho_total)
    ch_centers = (rho_total .- rho_min) ./ (rho_max - rho_min + 1e-12)
    ch_nodes = sStar(ch_centers)

    for i in eachindex(u)
        u[i][1] = v0_piston * ch_nodes[i]
    end

    return u

end   



##########
# INITIAL ALMANSI TENSOR
##########

function InitAlmansiTens()

    E = [zeros(3,3) for i in 1:(N + 2 * fict - 1)]

    return E

end