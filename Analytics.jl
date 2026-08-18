#=
Приближенное аналитическое решение задачи падения твердого упругого тела в жидкой среде под  действием внешней массовой силы.

Запуск: julia MakeFigures.jl <config>.jl
Результат: Figures/<config>
=#

####################################################################################################

#=
Новая версия приближенного аналитического решения.

Получено дифференциальное уравнение второго порядка (закон сохранения импульса) для величин x_p=x_p(t) и l_p=l_p(t) (координата
центра масс поршня и его длина соответственно).

На данный момент отсутствует замыкающее соотношение на x_p и l_p.
=#

# --- возвращает параметры задачи одним объектом
function BuildParameters()
    # --- отсчетная геометрия ---
    x0_l   = x0_piston - 0.5 * piston_width   
    x0_r   = x0_piston + 0.5 * piston_width   
    l0_p   = piston_width   

    # --- постоянные физические параметры задачи ---
    kappa = muel + 0.5 * lamel               
    m_l   = rho2A * x0_l                       
    m_r   = rho2A * (L - x0_r)                  
    m_p   = rho1B * l0_p                        
    m_tot = m_l + m_r + m_p

    # --- коэффициенты при производных в балансе импульса ---
    a_in    = 0.5 * (m_l + m_r) + m_p          # при d²x_p
    b_in    = 0.25 * (m_r - m_l)               # при d²l_p
    visc    = 4.0/3.0 * eta + zeta             # (4η/3 + ζ)
    g_force = -g                              # сила направлена в сторону убывания x

    # --- параметры модели одним объектом ---
    return (L=L, x0_piston=x0_piston, v0_piston=v0_piston, x0_l=x0_l, x0_r=x0_r, l0_p=l0_p, kappa=kappa,
            rho1B=rho1B, rho2A=rho2A, Apsi=Apsi, Bpsi=Bpsi,
            m_tot=m_tot, a_in=a_in, b_in=b_in, visc=visc, g_force=g_force)
end
params = BuildParameters()

# --- координаты левой и правой стенки поршня ---
function PistonGeometry(x_p, l_p)
    x_l  = x_p - 0.5 * l_p
    x_r  = x_p + 0.5 * l_p
    return x_l, x_r
end
 
# --- плотности столбов жидкости ---
function ColumnDensities(x_l, x_r, params)
    rho_l = params.rho2A * params.x0_l / x_l
    rho_r = params.rho2A * (params.L - params.x0_r) / (params.L - x_r)
    return rho_l, rho_r
end
 
# --- упругая сила от перепада давлений столбов ---
function SpringForce(rho_l, rho_r, params)
    return params.Apsi * (rho_r^2 - rho_l^2)
end
 
# --- вязкая сила  ---
function ViscousForce(u_l, u_r, x_l, x_r, params)
    return params.visc * (u_r / (params.L - x_r) + u_l / x_l)
end
 
# --- баланс импульса ---
function MomentumBalance(x_p, v_xp, l_p, v_lp, params)
    u_r = v_xp + 0.5 * v_lp
    u_l = v_xp - 0.5 * v_lp
    x_l, x_r     = PistonGeometry(x_p, l_p)
    rho_l, rho_r = ColumnDensities(x_l, x_r, params)
    f_spring     = SpringForce(rho_l, rho_r, params)
    f_visc       = ViscousForce(u_l, u_r, x_l, x_r, params)
    r = params.m_tot * params.g_force - f_spring - f_visc
    return params.a_in, params.b_in, r
end

# --- система дифференциальных уравнений ---
function DiffEquation(dx, x, params, t)
    x_p  = params.x0_piston + x[1] 
    v_xp = x[2]                    
 
    # Замыкание ещё не выведено, в качестве заглушки используется предельный случай жёсткого поршня
    l_p  = params.l0_p
    v_lp = 0.0
    a_in, b_in, r = MomentumBalance(x_p, v_xp, l_p, v_lp, params)
    ac_p = r / a_in
 
    dx[1] = v_xp
    dx[2] = ac_p
end
 
# --- решение ОДУ ---
function SolveAnalytics()
    x0    = [0.0, v0_piston]
    tspan = (0.0, t_max)
    prob  = ODEProblem(DiffEquation, x0, tspan, params)
    return solve(prob, reltol=1e-8, abstol=1e-10)
end
sol = SolveAnalytics()

####################################################################################################

#=
Старая версия приближенного аналитического решения. 

Давала сильное расхождение с численным решением, а также требовала введения калибровочных коэффициентов для регулировки частоты 
колебаний и степени затухания.
=#

#=
rhoA = rho1A + rho2A
rhoB = rho1B + rho2B

x0_left = x0_piston - 0.5 * piston_width
x0_right = L - (x0_piston + 0.5 * piston_width)

p0_left = rhoA * g * x0_left
p0_right = rhoA * g * x0_right

m_piston = rhoB * S * piston_width
m_eff = 1/3 * rhoA * S * (L - piston_width)

c = (4/3 * eta + zeta) * S * (1/x0_left + 1/x0_right)

function HydroForce(x)
    p1 = p0_left * ((x0_left * S) / ((x0_left + x) * S))^gamma
    p2 = p0_right * ((x0_right * S) / ((x0_right - x) * S))^gamma
    return (p1 - p2) * S
end

function DiffEquation(dx, x, p, t)
    dx[1] = x[2]
    dx[2] = 1/(m_coeff * (m_piston + m_eff)) * (-m_piston * g + HydroForce(x[1]) - c_coeff * c * dx[1])
end

x0 = [0.0, v0_piston]
tspan = (0.0, t_max)
prob = ODEProblem(DiffEquation, x0, tspan)
sol = solve(prob)

function BuildDensity(xc_piston)
    n = length(xc_in)
    rho1 = zeros(n)
    rho2 = zeros(n)

    for i in 1:n
        x_c = xc_in[i]
        if bound == "periodic"
            dist = mod(x_c - xc_piston + 0.5 * L, L) - 0.5 * L
        elseif bound == "wall"
            dist = x_c - xc_piston
        end
        if rho_init == "sharp"
            if abs(dist) <= 0.5 * piston_width
                f = 1.0 
            else
                f = 0.0 
            end
        elseif rho_init == "smooth"
            f = 0.5 * (1.0 + tanh(b * (0.5 * piston_width + dist))) *
                0.5 * (1.0 + tanh(b * (0.5 * piston_width - dist)))
        end
        rho1[i] = rho1A + (rho1B - rho1A) * f
        rho2[i] = rho2A + (rho2B - rho2A) * f
    end
    return rho1, rho2
end
=#