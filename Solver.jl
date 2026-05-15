heaviside(u) = 0.5 * (sign(u) + 1.0)

function ConservMassAndMomentumStep!(rho1_new, rho2_new, u_new, buf,
                                     rho1, rho2, u, E, Phi, dt,
                                     mu_el, lam_el, mu1h, mu2h, M_arr,
                                     K_current, Pi_el_current, Pi_ns_current)

    Ncentr = N + 2 * fict - 1
    Nface  = N + 2 * fict

    u_x = buf.u_x
    @inbounds for j in 1:Nface
        u_x[j] = u[j][1]
    end

    sStar!(buf.sStar_rho1, rho1)
    sStar!(buf.sStar_rho2, rho2)

    @inbounds for i in 1:Nface
        buf.flux1[i] = buf.sStar_rho1[i] * u_x[i]
        buf.flux2[i] = buf.sStar_rho2[i] * u_x[i]
    end
    BoundCondFlux(buf.flux1)
    BoundCondFlux(buf.flux2)

    #=
    flux1 = zeros(N + 2*fict)
    flux2 = zeros(N + 2*fict)
    for i in fict+1:N+fict
        if u_x[i] >= 0
            flux1[i] = rho1[i-1] * u_x[i]
            flux2[i] = rho2[i-1] * u_x[i]
        else
            flux1[i] = rho1[i] * u_x[i]
            flux2[i] = rho2[i] * u_x[i]
        end
    end
    flux1 = BoundCondFlux(flux1)
    flux2 = BoundCondFlux(flux2)
    =#

    d!(buf.d_rho1_u, buf.flux1)
    d!(buf.d_rho2_u, buf.flux2)

    BoundCondScalar(mu1h)
    BoundCondScalar(mu2h)

    sStar!(buf.sStar_M, M_arr)

    inv_theta = 1 / theta
    @inbounds for i in 1:Ncentr
        buf.theta_mu_diff[i] = inv_theta * (mu1h[i] - mu2h[i])
    end
    dStar!(buf.dStar_theta_mu, buf.theta_mu_diff)

    @inbounds for i in 1:Nface
        buf.flux_M[i] = buf.sStar_M[i] * buf.dStar_theta_mu[i]
    end
    BoundCondFlux(buf.flux_M)

    d!(buf.d_M_Mu, buf.flux_M)

    @inbounds for i in 1:Ncentr
        rho1_new[i] = rho1[i] + dt * (-buf.d_rho1_u[i] + buf.d_M_Mu[i])
        rho2_new[i] = rho2[i] + dt * (-buf.d_rho2_u[i] - buf.d_M_Mu[i])
    end
    BoundCondScalar(rho1_new)
    BoundCondScalar(rho2_new)

    rho = buf.rho
    @inbounds for i in 1:Ncentr
        rho[i] = rho1[i] + rho2[i]
    end

    Jx = buf.Jx
    sStar!(Jx, rho)
    @inbounds for i in 1:Nface
        Jx[i] *= u_x[i]
    end

    s!(buf.s_J, Jx)
    s!(buf.s_u_x, u_x)
    @inbounds for i in 1:Ncentr
        buf.fluxJ[i] = buf.s_J[i] * buf.s_u_x[i]
    end
    BoundCondScalar(buf.fluxJ)
    dStar!(buf.d_J_u, buf.fluxJ)

    sStar!(buf.sStar_rho, rho)
    dStar!(buf.dStar_mu1, mu1h)
    dStar!(buf.dStar_mu2, mu2h)

    @inbounds for i in 1:Nface
        buf.rho_Mu_sum[i] = buf.sStar_rho1[i] * buf.dStar_mu1[i] +
                            buf.sStar_rho2[i] * buf.dStar_mu2[i]
    end

    fill!(buf.d_E_K_x, 0.0)
    @inbounds for r in 1:3, q in 1:3
        for i in 1:Ncentr
            buf.E_rq[i] = E[i][r, q]
            buf.K_rq[i] = K_current[i][r, q]
        end
        dStar!(buf.dStar_E_rq, buf.E_rq)
        sStar!(buf.sStar_K_rq, buf.K_rq)
        for i in 1:Nface
            buf.d_E_K_x[i] += buf.dStar_E_rq[i] * buf.sStar_K_rq[i]
        end
    end

    @inbounds for i in 1:Ncentr
        buf.Pi_diag[i] = Pi_ns_current[i][1, 1] + Pi_el_current[i][1, 1]
    end
    dStar!(buf.d_Pi_x, buf.Pi_diag)

    dStar!(buf.dStar_Phi, Phi)
    @inbounds for i in 1:Nface
        buf.force[i] = buf.sStar_rho[i] * buf.dStar_Phi[i]
    end

    @inbounds for i in 1:Ncentr
        buf.rho_new_sum[i] = rho1_new[i] + rho2_new[i]
    end
    sStar!(buf.denom, buf.rho_new_sum)

    @inbounds for i in 1:Nface
        buf.u_x_new[i] = (buf.sStar_rho[i] * u_x[i] +
                         dt * (-buf.d_J_u[i] - buf.rho_Mu_sum[i] +
                               buf.d_E_K_x[i] + buf.d_Pi_x[i] + buf.force[i])) / buf.denom[i]
    end

    @inbounds for i in 1:Nface
        u_new[i][1] = buf.u_x_new[i]
    end

    return nothing
end

function ConservAlmansiTensor!(E_new, buf, E, u, dt)
    Ncentr = N + 2 * fict - 1
    Nface  = N + 2 * fict

    u_x = buf.u_x
    @inbounds for j in 1:Nface
        u_x[j] = u[j][1]
    end

    d!(buf.du_dx, u_x)

    @inbounds for j in 1:3
        for l in 1:Ncentr
            buf.E_jj[l] = E[l][j, j]
        end

        dStar!(buf.dstar_E_jj, buf.E_jj)

        for i in 1:Nface
            buf.u_x_dstar_E[i] = u_x[i] * buf.dstar_E_jj[i]
        end
        s!(buf.conv, buf.u_x_dstar_E)

        if j == 1
            for i in fict + 1:N + fict - 1
                E_new[i][j, j] = E[i][j, j] + dt * (-buf.conv[i] + buf.du_dx[i] - 2 * E[i][j, j] * buf.du_dx[i])
            end
        else
            for i in fict + 1:N + fict - 1
                E_new[i][j, j] = E[i][j, j] + dt * (-buf.conv[i])
            end
        end
    end

    return nothing
end


function CheckStability(rho1, rho2, u, step, t, E, mu_el, lam_el, Phi)
    u_x = [u[i][1] for i in eachindex(u)]

    if any(isnan.(rho1)) || any(isnan.(rho2)) || any(isnan.(u_x))
        println("NaN, step = $step, t = $t")
        return false
    end

    if maximum(abs.(rho1)) > 1e6 || maximum(abs.(rho2)) > 1e6
        println("Blow-up rho, step = $step, t = $t, max|rho1| = $(maximum(abs.(rho1))), max|rho2| = $(maximum(abs.(rho2)))")
        return false
    end

	if maximum(abs.(u_x)) > 1e4
        println("Blow-up u, step = $step, t = $t, max|u| = $(maximum(abs.(u_x)))")
        return false
    end


    return true

end