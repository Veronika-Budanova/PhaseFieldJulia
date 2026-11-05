#=
Файл, содержащий функции законов сохранения массы, импульса и энергии на одном шаге по времени.
Включен в Main.jl через include.

Методы для матриц - двумерная задача, для векторов - одномерная задача вдоль y (u = 0, d/dx = 0).
ws - рабочие массивы (Workspace.jl), выделенные один раз до цикла по времени.
Порядок арифметических операций совпадает с выделяющей версией, поэтому результат совпадает побитово.
Двумерные методы параллельны (Polyester: @batch и Map2D! из Parallel.jl), одномерные - последовательные.
=#

####################################################################################################
# --- общие для обеих размерностей ---

# --- новый тензор Альманси в ячейке: обновляются xx, yy, zz и xy = yx, ---
# --- компоненты xz, yz, zx, zy схемой не обновляются и берутся из прежнего значения ячейки Eo ---
# --- аргументы Tensor3 перечисляются по столбцам: (1,1), (2,1), (3,1), (1,2), ... ---
@inline AlmansiCell(Eo, e11, e22, e33, e12) = Tensor3(e11,      e12,      Eo[3, 1],
                                                      e12,      e22,      Eo[3, 2],
                                                      Eo[1, 3], Eo[2, 3], e33)

####################################################################################################
# --- двумерная задача ---

# --- компонента (r, q) тензорного поля в отдельный скалярный массив ---
function Component!(out, A::AbstractMatrix{<:AbstractMatrix}, r, q)
    @batch for idx in eachindex(A)
        @inbounds out[idx] = A[idx][r, q]
    end
    return out
end

function ConservMass!(rho1_new, rho2_new, rho1::AbstractMatrix, rho2, U, dt, mu1h, mu2h, M_arr, ws)
    u = U[1]
    v = U[2]

    flux1_x = Map2D!(*, sStarx!(ws.flux1_x, rho1), ws.flux1_x, u)
    flux1_y = Map2D!(*, sStary!(ws.flux1_y, rho1), ws.flux1_y, v)
    flux2_x = Map2D!(*, sStarx!(ws.flux2_x, rho2), ws.flux2_x, u)
    flux2_y = Map2D!(*, sStary!(ws.flux2_y, rho2), ws.flux2_y, v)
    BoundCondVector((flux1_x, flux1_y))
    BoundCondVector((flux2_x, flux2_y))

    d_rho1_u = dx!(ws.d_rho1_u, flux1_x); dy!(ws.tmp_c, flux1_y)
    Map2D!(+, d_rho1_u, d_rho1_u, ws.tmp_c)
    d_rho2_u = dx!(ws.d_rho2_u, flux2_x); dy!(ws.tmp_c, flux2_y)
    Map2D!(+, d_rho2_u, d_rho2_u, ws.tmp_c)

    BoundCondScalar(mu1h)
    BoundCondScalar(mu2h)

    inv_theta     = 1 / theta
    theta_mu_diff = Map2D!((a, b) -> inv_theta * (a - b), ws.theta_mu_diff, mu1h, mu2h)

    flux_M_x = sStarx!(ws.flux_M_x, M_arr); dStarx!(ws.tmp_fx, theta_mu_diff)
    Map2D!(*, flux_M_x, flux_M_x, ws.tmp_fx)
    flux_M_y = sStary!(ws.flux_M_y, M_arr); dStary!(ws.tmp_fy, theta_mu_diff)
    Map2D!(*, flux_M_y, flux_M_y, ws.tmp_fy)
    BoundCondVector((flux_M_x, flux_M_y))

    d_M_Mu = dx!(ws.d_M_Mu, flux_M_x); dy!(ws.tmp_c, flux_M_y)
    Map2D!(+, d_M_Mu, d_M_Mu, ws.tmp_c)

    @batch for idx in eachindex(rho1)
        @inbounds begin
            rho1_new[idx] = rho1[idx] + dt * (-d_rho1_u[idx] + d_M_Mu[idx])
            rho2_new[idx] = rho2[idx] + dt * (-d_rho2_u[idx] - d_M_Mu[idx])
        end
    end
    BoundCondScalar(rho1_new)
    BoundCondScalar(rho2_new)

    return nothing
end

function ConservMomentum!(U_new, rho1_new, rho2_new, rho1::AbstractMatrix, rho2, U, E, Phi, dt,
                          mu1h, mu2h, K_current, Pi_el_current, Pi_ns_current, ws)
    u     = U[1]
    v     = U[2]
    u_new = U_new[1]
    v_new = U_new[2]

    # --- K и Pi_ns посчитаны без фиктивных ячеек, а на гранях у границы области они нужны (периодика) ---
    BoundCondTensor(K_current)
    BoundCondTensor(Pi_ns_current)

    rho = Map2D!(+, ws.rho, rho1, rho2)

    sStarx_rho = sStarx!(ws.sStarx_rho, rho)
    sStary_rho = sStary!(ws.sStary_rho, rho)

    Jx = Map2D!(*, ws.Jx, sStarx_rho, u)
    Jy = Map2D!(*, ws.Jy, sStary_rho, v)

    fluxJ_xx = sx!(ws.fluxJ_xx, Jx); sx!(ws.tmp_c, u)
    Map2D!(*, fluxJ_xx, fluxJ_xx, ws.tmp_c)
    fluxJ_yy = sy!(ws.fluxJ_yy, Jy); sy!(ws.tmp_c, v)
    Map2D!(*, fluxJ_yy, fluxJ_yy, ws.tmp_c)
    BoundCondScalar(fluxJ_xx)
    BoundCondScalar(fluxJ_yy)
    fluxJ_yx = sStarx!(ws.fluxJ_yx, Jy); sStary!(ws.tmp_n, u)
    Map2D!(*, fluxJ_yx, fluxJ_yx, ws.tmp_n)
    fluxJ_xy = sStary!(ws.fluxJ_xy, Jx); sStarx!(ws.tmp_n, v)
    Map2D!(*, fluxJ_xy, fluxJ_xy, ws.tmp_n)

    d_J_u = dStarx!(ws.d_J_u, fluxJ_xx); dy!(ws.tmp_fx, fluxJ_yx)
    Map2D!(+, d_J_u, d_J_u, ws.tmp_fx)
    d_J_v = dStary!(ws.d_J_v, fluxJ_yy); dx!(ws.tmp_fy, fluxJ_xy)
    Map2D!(+, d_J_v, d_J_v, ws.tmp_fy)

    # --- sStar(rho1) * dStar(mu1h) + sStar(rho2) * dStar(mu2h) ---
    rho_Mu_sum_x = sStarx!(ws.rho_Mu_sum_x, rho1); dStarx!(ws.tmp_fx, mu1h)
    Map2D!(*, rho_Mu_sum_x, rho_Mu_sum_x, ws.tmp_fx)
    sStarx!(ws.tmp_fx2, rho2); dStarx!(ws.tmp_fx, mu2h)
    Map2D!((a, b, c) -> a + b * c, rho_Mu_sum_x, rho_Mu_sum_x, ws.tmp_fx2, ws.tmp_fx)
    rho_Mu_sum_y = sStary!(ws.rho_Mu_sum_y, rho1); dStary!(ws.tmp_fy, mu1h)
    Map2D!(*, rho_Mu_sum_y, rho_Mu_sum_y, ws.tmp_fy)
    sStary!(ws.tmp_fy2, rho2); dStary!(ws.tmp_fy, mu2h)
    Map2D!((a, b, c) -> a + b * c, rho_Mu_sum_y, rho_Mu_sum_y, ws.tmp_fy2, ws.tmp_fy)

    d_E_K_x = fill!(ws.d_E_K_x, 0.0)
    d_E_K_y = fill!(ws.d_E_K_y, 0.0)
    for r in 1:3, q in 1:3
        E_rq = Component!(ws.E_rq, E, r, q)
        K_rq = Component!(ws.K_rq, K_current, r, q)
        dStarx!(ws.tmp_fx, E_rq); sStarx!(ws.tmp_fx2, K_rq)
        Map2D!((d, a, b) -> d + a * b, d_E_K_x, d_E_K_x, ws.tmp_fx, ws.tmp_fx2)
        dStary!(ws.tmp_fy, E_rq); sStary!(ws.tmp_fy2, K_rq)
        Map2D!((d, a, b) -> d + a * b, d_E_K_y, d_E_K_y, ws.tmp_fy, ws.tmp_fy2)
    end

    BoundCondTensor(Pi_el_current)

    Pi_xx    = ws.Pi_xx
    Pi_yy    = ws.Pi_yy
    Pi_el_xy = Component!(ws.Pi_el_xy, Pi_el_current, 1, 2)
    @batch for idx in eachindex(E)
        @inbounds begin
            Pi_xx[idx] = Pi_ns_current[idx][1, 1] + Pi_el_current[idx][1, 1]
            Pi_yy[idx] = Pi_ns_current[idx][2, 2] + Pi_el_current[idx][2, 2]
        end
    end

    # --- eta * (dStary(u) + dStarx(v)) + sStarx(sStary(Pi_el_xy)) ---
    Pi_xy = dStary!(ws.Pi_xy, u)
    dStarx!(ws.tmp_n, v)
    sStarx!(ws.tmp_n2, sStary!(ws.tmp_fy, Pi_el_xy))
    Map2D!((a, b, c) -> eta * (a + b) + c, Pi_xy, Pi_xy, ws.tmp_n, ws.tmp_n2)

    d_Pi_x = dStarx!(ws.d_Pi_x, Pi_xx); dy!(ws.tmp_fx, Pi_xy)
    Map2D!(+, d_Pi_x, d_Pi_x, ws.tmp_fx)
    d_Pi_y = dStary!(ws.d_Pi_y, Pi_yy); dx!(ws.tmp_fy, Pi_xy)
    Map2D!(+, d_Pi_y, d_Pi_y, ws.tmp_fy)

    force_x = Map2D!(*, dStarx!(ws.force_x, Phi), ws.force_x, sStarx_rho)
    force_y = Map2D!(*, dStary!(ws.force_y, Phi), ws.force_y, sStary_rho)

    rho_new = Map2D!(+, ws.rho_new, rho1_new, rho2_new)
    denom_x = sStarx!(ws.denom_x, rho_new)
    denom_y = sStary!(ws.denom_y, rho_new)

    @batch for idx in eachindex(u_new)
        @inbounds u_new[idx] = (sStarx_rho[idx] * u[idx] +
                                dt * (-d_J_u[idx] - rho_Mu_sum_x[idx] +
                                      d_E_K_x[idx] + d_Pi_x[idx] + force_x[idx])) / denom_x[idx]
    end
    @batch for idx in eachindex(v_new)
        @inbounds v_new[idx] = (sStary_rho[idx] * v[idx] +
                                dt * (-d_J_v[idx] - rho_Mu_sum_y[idx] +
                                      d_E_K_y[idx] + d_Pi_y[idx] + force_y[idx])) / denom_y[idx]
    end

    return nothing
end

# --- конвективное слагаемое sx(u * dStarx(A)) + sy(v * dStary(A)) ---
function Convection!(conv, A::AbstractMatrix, u, v, ws)
    tmp_fx = dStarx!(ws.tmp_fx, A)
    Map2D!(*, tmp_fx, tmp_fx, u)
    sx!(conv, tmp_fx)
    tmp_fy = dStary!(ws.tmp_fy, A)
    Map2D!(*, tmp_fy, tmp_fy, v)
    sy!(ws.tmp_c, tmp_fy)
    Map2D!(+, conv, conv, ws.tmp_c)
    return conv
end

function ConservAlmansiTensor!(E_new, E::AbstractMatrix{<:AbstractMatrix}, U, dt, ws)
    u = U[1]
    v = U[2]

    E_xx = Component!(ws.E_xx, E, 1, 1)
    E_yy = Component!(ws.E_yy, E, 2, 2)
    E_zz = Component!(ws.E_zz, E, 3, 3)
    E_xy = Component!(ws.E_xy, E, 1, 2)

    du_dx = dx!(ws.du_dx, u)
    dv_dy = dy!(ws.dv_dy, v)
    dv_dx = dx!(ws.dv_dx, sStarx!(ws.tmp_fx, sy!(ws.tmp_c, v)))
    du_dy = dy!(ws.du_dy, sStary!(ws.tmp_fy, sx!(ws.tmp_c, u)))

    conv_xx = Convection!(ws.conv_xx, E_xx, u, v, ws)
    conv_yy = Convection!(ws.conv_yy, E_yy, u, v, ws)
    conv_zz = Convection!(ws.conv_zz, E_zz, u, v, ws)
    conv_xy = Convection!(ws.conv_xy, E_xy, u, v, ws)

    rows = ficty + 1:Ny + ficty - 1
    @batch for j in fictx + 1:Nx + fictx - 1
        @inbounds for i in rows
            exx = E_xx[i, j]
            eyy = E_yy[i, j]
            exy = E_xy[i, j]
            ux  = du_dx[i, j]
            vy  = dv_dy[i, j]
            vx  = dv_dx[i, j]
            uy  = du_dy[i, j]

            e11 = exx + dt * (-conv_xx[i, j] + ux - 2 * exx * ux - 2 * exy * vx)
            e22 = eyy + dt * (-conv_yy[i, j] + vy - 2 * eyy * vy - 2 * exy * uy)
            e33 = E_zz[i, j] + dt * (-conv_zz[i, j])
            e12 = exy + dt * (-conv_xy[i, j] + 0.5 * (vx + uy)
                              - exy * (ux + vy) - exx * uy - eyy * vx)
            E_new[i, j] = AlmansiCell(E_new[i, j], e11, e22, e33, e12)
        end
    end

    return nothing
end

####################################################################################################
# --- одномерная задача (вдоль y): u = 0, d/dx = 0 ---

# --- компонента (r, q) тензорного поля в отдельный скалярный массив ---
function Component!(out, A::AbstractVector{<:AbstractMatrix}, r, q)
    @inbounds for idx in eachindex(A)
        out[idx] = A[idx][r, q]
    end
    return out
end

function ConservMass!(rho1_new, rho2_new, rho1::AbstractVector, rho2, U, dt, mu1h, mu2h, M_arr, ws)
    v = U

    flux1_y = sStary!(ws.flux1_y, rho1); flux1_y .*= v
    flux2_y = sStary!(ws.flux2_y, rho2); flux2_y .*= v
    BoundCondVector(flux1_y)
    BoundCondVector(flux2_y)

    d_rho1_u = dy!(ws.d_rho1_u, flux1_y)
    d_rho2_u = dy!(ws.d_rho2_u, flux2_y)

    BoundCondScalar(mu1h)
    BoundCondScalar(mu2h)

    inv_theta     = 1 / theta
    theta_mu_diff = ws.theta_mu_diff
    @. theta_mu_diff = inv_theta * (mu1h - mu2h)

    flux_M_y = sStary!(ws.flux_M_y, M_arr); flux_M_y .*= dStary!(ws.tmp_fy, theta_mu_diff)
    BoundCondVector(flux_M_y)

    d_M_Mu = dy!(ws.d_M_Mu, flux_M_y)

    @inbounds for idx in eachindex(rho1)
        rho1_new[idx] = rho1[idx] + dt * (-d_rho1_u[idx] + d_M_Mu[idx])
        rho2_new[idx] = rho2[idx] + dt * (-d_rho2_u[idx] - d_M_Mu[idx])
    end
    BoundCondScalar(rho1_new)
    BoundCondScalar(rho2_new)

    return nothing
end

function ConservMomentum!(U_new, rho1_new, rho2_new, rho1::AbstractVector, rho2, U, E, Phi, dt,
                          mu1h, mu2h, K_current, Pi_el_current, Pi_ns_current, ws)
    v     = U
    v_new = U_new

    # --- K и Pi_ns посчитаны без фиктивных ячеек, а на гранях у границы области они нужны (периодика) ---
    BoundCondTensor(K_current)
    BoundCondTensor(Pi_ns_current)

    rho = ws.rho
    @. rho = rho1 + rho2

    sStary_rho = sStary!(ws.sStary_rho, rho)

    Jy = ws.Jy; @. Jy = sStary_rho * v

    fluxJ_yy = sy!(ws.fluxJ_yy, Jy); fluxJ_yy .*= sy!(ws.tmp_c, v)
    BoundCondScalar(fluxJ_yy)

    d_J_v = dStary!(ws.d_J_v, fluxJ_yy)

    # --- sStary(rho1) * dStary(mu1h) + sStary(rho2) * dStary(mu2h) ---
    rho_Mu_sum_y = sStary!(ws.rho_Mu_sum_y, rho1); rho_Mu_sum_y .*= dStary!(ws.tmp_fy, mu1h)
    sStary!(ws.tmp_fy2, rho2); dStary!(ws.tmp_fy, mu2h)
    @. rho_Mu_sum_y = rho_Mu_sum_y + ws.tmp_fy2 * ws.tmp_fy

    d_E_K_y = fill!(ws.d_E_K_y, 0.0)
    for r in 1:3, q in 1:3
        E_rq = Component!(ws.E_rq, E, r, q)
        K_rq = Component!(ws.K_rq, K_current, r, q)
        dStary!(ws.tmp_fy, E_rq); sStary!(ws.tmp_fy2, K_rq)
        @. d_E_K_y = d_E_K_y + ws.tmp_fy * ws.tmp_fy2
    end

    BoundCondTensor(Pi_el_current)

    Pi_yy = ws.Pi_yy
    @inbounds for idx in eachindex(E)
        Pi_yy[idx] = Pi_ns_current[idx][2, 2] + Pi_el_current[idx][2, 2]
    end

    d_Pi_y = dStary!(ws.d_Pi_y, Pi_yy)

    force_y = dStary!(ws.force_y, Phi); force_y .*= sStary_rho

    rho_new = ws.rho_new
    @. rho_new = rho1_new + rho2_new
    denom_y = sStary!(ws.denom_y, rho_new)

    @inbounds for idx in eachindex(v_new)
        v_new[idx] = (sStary_rho[idx] * v[idx] +
                      dt * (-d_J_v[idx] - rho_Mu_sum_y[idx] +
                            d_E_K_y[idx] + d_Pi_y[idx] + force_y[idx])) / denom_y[idx]
    end

    return nothing
end

# --- конвективное слагаемое sy(v * dStary(A)) ---
function Convection!(conv, A::AbstractVector, v, ws)
    tmp_fy = dStary!(ws.tmp_fy, A); tmp_fy .*= v
    sy!(conv, tmp_fy)
    return conv
end

function ConservAlmansiTensor!(E_new, E::AbstractVector{<:AbstractMatrix}, U, dt, ws)
    v = U

    E_xx = Component!(ws.E_xx, E, 1, 1)
    E_yy = Component!(ws.E_yy, E, 2, 2)
    E_zz = Component!(ws.E_zz, E, 3, 3)
    E_xy = Component!(ws.E_xy, E, 1, 2)

    dv_dy = dy!(ws.dv_dy, v)

    conv_xx = Convection!(ws.conv_xx, E_xx, v, ws)
    conv_yy = Convection!(ws.conv_yy, E_yy, v, ws)
    conv_zz = Convection!(ws.conv_zz, E_zz, v, ws)
    conv_xy = Convection!(ws.conv_xy, E_xy, v, ws)

    @inbounds for i in ficty + 1:Ny + ficty - 1
        exx = E_xx[i]
        eyy = E_yy[i]
        exy = E_xy[i]
        vy  = dv_dy[i]

        e11 = exx + dt * (-conv_xx[i])
        e22 = eyy + dt * (-conv_yy[i] + vy - 2 * eyy * vy)
        e33 = E_zz[i] + dt * (-conv_zz[i])
        e12 = exy + dt * (-conv_xy[i] - exy * vy)
        E_new[i] = AlmansiCell(E_new[i], e11, e22, e33, e12)
    end

    return nothing
end
