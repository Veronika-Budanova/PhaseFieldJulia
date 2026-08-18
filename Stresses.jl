# PiNS, PiEl
function PiNS!(Pi, u, buf)
    u_x = buf.u_x
    @inbounds for j in 1:N + 2 * fict
        u_x[j] = u[j][1]
    end
    d!(buf.du_dx, u_x)
    @inbounds for i in 1:N + 2 * fict - 1
        div_u = buf.du_dx[i]
        Pii = Pi[i]
        Pii[1, 1] = 2 * eta * div_u + (zeta - (2 / 3) * eta) * div_u
        Pii[1, 2] = 0.0
        Pii[1, 3] = 0.0
        Pii[2, 1] = 0.0
        Pii[2, 2] = 0.0
        Pii[2, 3] = 0.0
        Pii[3, 1] = 0.0
        Pii[3, 2] = 0.0
        Pii[3, 3] = 0.0
    end
    return Pi
end

function PiEl!(Pi, E, K)
    @inbounds for i in fict + 1:N + fict - 1
        Ei = E[i]
        Ki = K[i]
        Pii = Pi[i]
        for j in 1:3
            for k in 1:3
                sum_r = 2 * (Ei[j, 1] * Ki[1, k] + Ei[j, 2] * Ki[2, k] + Ei[j, 3] * Ki[3, k])
                Pii[j, k] = Ki[j, k] - sum_r
            end
        end
    end
    return Pi
end

# K
function K!(K_ar, E, mu_el, lam_el)
    @inbounds for i in fict+1:N+fict-1
        mui = mu_el[i]
        lami = lam_el[i]
        Ei = E[i]
        Ki = K_ar[i]
        trE = Ei[1,1] + Ei[2,2] + Ei[3,3]
        for j in 1:3
            for k in 1:3
                Ki[j,k] = 2*mui * Ei[j,k] + (j == k ? lami * trE : 0.0)
            end
        end
    end
    return K_ar
end

# p
function Pressure(E, mu_el, lam_el, rho1, rho2)
	mu1 = Mu1!(zeros(length(rho1)), E, rho1, rho2)
	mu2 = Mu2!(zeros(length(rho1)), E, rho1, rho2)
	psi0 = Psi0(E, mu_el, lam_el, rho1, rho2)
	return rho1 .* mu1 .+ rho2 .* mu2 .- psi0
end
