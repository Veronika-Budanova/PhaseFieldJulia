heaviside(u) = 0.5 * (sign(u) + 1.0)

function ConservMassAndMomentumStep(rho1, rho2, u, E, Phi, dt)
	u_x = [u[j][1] for j in 1:N + 2 * fict]

	
	flux1 = zeros(N + 2 * fict)
	flux2 = zeros(N + 2 * fict)

	sStar_rho1 = sStar(rho1)
	sStar_rho2 = sStar(rho2)

	flux1 = sStar_rho1 .* u_x
	flux2 = sStar_rho2 .* u_x

	flux1 = BoundCondFlux(flux1)
    flux2 = BoundCondFlux(flux2)
	
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
	d_rho1_u = d(flux1)
	d_rho2_u = d(flux2)

	mu_el = MuEl(muelB, rho1, rho2)
	lam_el = LamEl(lamelB, rho1, rho2)

    mu1h = Mu1Hat(E, mu_el, lam_el, rho1, rho2)
    mu2h = Mu2Hat(E, mu_el, lam_el, rho1, rho2)

	mu1h = BoundCondScalar(mu1h)
	mu2h = BoundCondScalar(mu2h)

	M_arr = M(M0, rho1, rho2)
	sStar_M = sStar(M_arr)
	theta_mu1h_mu2h = theta^(-1) .* (mu1h .- mu2h)
	dStar_theta_mu1h_mu2h = dStar(theta_mu1h_mu2h)
	sStar_M_dStar_theta_mu1h_mu2h = sStar_M .* dStar_theta_mu1h_mu2h
	
	sStar_M_dStar_theta_mu1h_mu2h = BoundCondFlux(sStar_M_dStar_theta_mu1h_mu2h)

	d_M_Mu = d(sStar_M_dStar_theta_mu1h_mu2h)
	
	rho1_new = rho1 .+ dt .* (-d_rho1_u .+ d_M_Mu) 
	rho2_new = rho2 .+ dt .* (-d_rho2_u .- d_M_Mu)

	rho1_new = BoundCondScalar(rho1_new)
	rho2_new = BoundCondScalar(rho2_new)

	#ConservMomentum

	rho = rho1 .+ rho2

	fluxJ = zeros(N + 2 * fict)
	Jx = J(rho, u_x)  
	s_J = s(Jx)
	s_u_x = s(u_x)
	fluxJ = s_J .* s_u_x
	fluxJ = BoundCondScalar(fluxJ)
	d_J_u = dStar(fluxJ)

	rho_Mu_sum = sStar(rho1) .* dStar(mu1h) .+
		         sStar(rho2) .* dStar(mu2h)

	K_current = K(E, mu_el, lam_el, rho1, rho2)
	Pi_el_current = PiEl(E, K_current)
	Pi_ns_current = PiNS(u)
    
	d_E_K_x = zeros(N + 2 * fict)
	for r in 1:3, q in 1:3
		E_rq = [E[i][r, q] for i in 1:N + 2 * fict - 1]
		K_rq = [K_current[i][r, q] for i in 1:N + 2 * fict - 1]
		d_E_K_x .+= dStar(E_rq) .* sStar(K_rq)
	end

	d_Pi_x = dStar([Pi_ns_current[i][1, 1] + Pi_el_current[i][1, 1] for i in eachindex(E)])

	sStarrho = sStar(rho)
	dStarPhi = dStar(Phi)
	force = sStarrho .* dStarPhi 

	denom = sStar(rho1_new .+ rho2_new)
	u_x_new = ((sStar(rho) .* u_x) .+ 
			dt .* (-d_J_u .- rho_Mu_sum .+ d_E_K_x .+ d_Pi_x .+ force)) ./ denom

	u_new = [copy(u[i]) for i in 1:N + 2 * fict]
    for i in 1:N + 2 * fict
        u_new[i][1] = u_x_new[i]
    end


	return rho1_new, rho2_new, u_new
end

function ConservAlmansiTensor(E, u, dt)
	u_x = [u[j][1] for j in 1:N + 2 * fict]
	E_new = [zeros(3, 3) for i in 1:N + 2 * fict - 1]

	du_dx = d(u_x)
    
	for j in 1:3
    	E_jj = [E[l][j, j] for l in 1:N + 2 * fict - 1]
        
        dstar_E = dStar(E_jj)           
        conv = s(u_x .* dstar_E)        
        
        for i in fict + 1:N + fict - 1
			if j == 1
            	E_new[i][j,j] = E[i][j,j] + dt * (-conv[i] + du_dx[i] - 2 * E[i][j,j] * du_dx[i])
			else
				E_new[i][j,j] = E[i][j,j] + dt * (-conv[i])
			end
		end
    end

	return E_new
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