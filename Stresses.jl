##########
# П NS, П el
##########
function PiNS(u)
    du_dx = d([u[j][1] for j in 1:N + 2 * fict])
    return [
        begin
            div_u = du_dx[i]
            m = zeros(3, 3)
            m[1, 1] = 2 * eta * du_dx[i] + (zeta - (2 / 3) * eta) * div_u
            m
        end 
		for i in 1:N + 2 * fict - 1
    ]
end

function PiEl(E, K)
	Pi = [zeros(3, 3) for i in 1:N + 2 * fict - 1]  
	for i in fict + 1:N + fict - 1
		for j in 1:3
			for k in 1:3
				sum_r = 0.0
				for r in 1:3
					sum_r += 2 * E[i][j, r] * K[i][r, k]
				end
				Pi[i][j, k] = K[i][j, k] - sum_r
			end
		end
	end
	return Pi
end

##########
# K
##########
function K(E, muel, lamel, rho1, rho2)
    K_ar = [zeros(3, 3) for _ in 1:N + 2*fict - 1]
    for i in fict+1:N+fict-1
        μ = muel[i]
        λ = lamel[i]
        trE = E[i][1,1] + E[i][2,2] + E[i][3,3]
        for j in 1:3
            for k in 1:3
                K_ar[i][j,k] = 2μ * E[i][j,k] + (j == k ? λ * trE : 0.0)
            end
        end
    end
    return K_ar
end

##########
# p
##########
function Pressure(E, muel, lamel, rho1, rho2)
	mu1 = Mu1(E, muel, lamel, rho1, rho2)
	mu2 = Mu2(E, muel, lamel, rho1, rho2)
	psi0 = Psi0(E, muel, lamel, rho1, rho2)
	return rho1 .* mu1 .+ rho2 .* mu2 .- psi0
end