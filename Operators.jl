##########
# s, d, s*, d*
##########

function s(array)
    s_array = zeros(N + 2 * fict - 1)
    @inbounds for i in (fict + 1):(N + fict - 1)
        s_array[i] = 0.5 * (array[i + 1] + array[i])
    end
    return s_array
end

function d(array)
    d_array = zeros(N + 2 * fict - 1)
    inv_h = 1 / h
    @inbounds for i in (fict + 1):(N + fict - 1)
        d_array[i] = inv_h * (array[i + 1] - array[i])
    end
    return d_array
end

function sStar(array)
    sStar_array = zeros(N + 2 * fict)
    @inbounds for i in (fict + 1):(N + fict)
        sStar_array[i] = 0.5 * (array[i] + array[i - 1])
    end
    return sStar_array
end

function dStar(array)
    dStar_array = zeros(N + 2 * fict)
    inv_h = 1 / h
    @inbounds for i in (fict + 1):(N + fict)
        dStar_array[i] = inv_h * (array[i] - array[i - 1])
    end
    return dStar_array
end

##########
# In-place варианты 
##########

function s!(out, array)
    @inbounds for i in (fict + 1):(N + fict - 1)
        out[i] = 0.5 * (array[i + 1] + array[i])
    end
    return out
end

function d!(out, array)
    inv_h = 1 / h
    @inbounds for i in (fict + 1):(N + fict - 1)
        out[i] = inv_h * (array[i + 1] - array[i])
    end
    return out
end

function sStar!(out, array)
    @inbounds for i in (fict + 1):(N + fict)
        out[i] = 0.5 * (array[i] + array[i - 1])
    end
    return out
end

function dStar!(out, array)
    inv_h = 1 / h
    @inbounds for i in (fict + 1):(N + fict)
        out[i] = inv_h * (array[i] - array[i - 1])
    end
    return out
end


function CreateBuffers()
    Ncentr = N + 2 * fict - 1
    Nface  = N + 2 * fict
    return (
        flux1            = zeros(Nface),
        flux2            = zeros(Nface),
        sStar_rho1       = zeros(Nface),
        sStar_rho2       = zeros(Nface),
        d_rho1_u         = zeros(Ncentr),
        d_rho2_u         = zeros(Ncentr),
        sStar_M          = zeros(Nface),
        theta_mu_diff    = zeros(Ncentr),
        dStar_theta_mu   = zeros(Nface),
        flux_M           = zeros(Nface),
        d_M_Mu           = zeros(Ncentr),
        u_x              = zeros(Nface),
        Jx               = zeros(Nface),
        s_J              = zeros(Ncentr),
        s_u_x            = zeros(Ncentr),
        fluxJ            = zeros(Ncentr),
        d_J_u            = zeros(Nface),
        sStar_rho        = zeros(Nface),
        dStar_mu1        = zeros(Nface),
        dStar_mu2        = zeros(Nface),
        rho_Mu_sum       = zeros(Nface),
        d_E_K_x          = zeros(Nface),
        E_rq             = zeros(Ncentr),
        K_rq             = zeros(Ncentr),
        dStar_E_rq       = zeros(Nface),
        sStar_K_rq       = zeros(Nface),
        Pi_diag          = zeros(Ncentr),
        d_Pi_x           = zeros(Nface),
        dStar_Phi        = zeros(Nface),
        force            = zeros(Nface),
        rho_new_sum      = zeros(Ncentr),
        denom            = zeros(Nface),
        u_x_new          = zeros(Nface),
        rho              = zeros(Ncentr),
        E_jj             = zeros(Ncentr),
        dstar_E_jj       = zeros(Nface),
        u_x_dstar_E      = zeros(Nface),
        conv             = zeros(Ncentr),
        du_dx            = zeros(Ncentr),
        mu_el            = zeros(Ncentr),
        lam_el           = zeros(Ncentr),
        H_arr            = zeros(Ncentr),
        mu1              = zeros(Ncentr),
        mu2              = zeros(Ncentr),
        mu1h             = zeros(Ncentr),
        mu2h             = zeros(Ncentr),
        M_arr            = zeros(Ncentr),
        lap_rho1         = zeros(Ncentr),
        lap_rho2         = zeros(Ncentr),
        dStar_rho1       = zeros(Nface),
        dStar_rho2       = zeros(Nface),
        K_current        = [zeros(3, 3) for _ in 1:Ncentr],
        Pi_el_current    = [zeros(3, 3) for _ in 1:Ncentr],
        Pi_ns_current    = [zeros(3, 3) for _ in 1:Ncentr],
    )
end