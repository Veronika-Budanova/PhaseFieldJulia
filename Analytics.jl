#=
Приближенное аналитическое решение для одномерной задачи и для срезов x = x_slice двумерной задачи.
Поршень одномерен, поэтому решение одинаково для 1D и для любого среза в 2D.

q = (y_p, l_p), dq = (dy_p/dt, dl_p/dt), состояние ОДУ s = [y_p, l_p, dy_p/dt, dl_p/dt].
Участки: нижний столб [0, Y_1], поршень [Y_1, Y_2], верхний столб [Y_2, Ly].

Задаются только M, Y, Rho, U и psi0. Все производные и все интегралы берутся символьно (SymPy), в явном виде,
и один раз переводятся в числовые функции (ToNum). Параметры конфига входят в символьные выражения точными
рациональными числами (Exact), иначе SymPy работает с float-коэффициентами медленно и с ошибками.

M_ik q''_i + C_ijk q'_i q'_j + dΨ/dq_k + D_ik q'_i + F_k = 0

Для MakeFigures.jl: AnalyticState(t), AnalyticProfilesY(yc_arr, y_arr, t).
=#

using SymPy

@syms y z y_p l_p dy_p dl_p Y_b w
q  = [y_p, l_p]
dq = [dy_p, dl_p]

# --- параметр конфига как точное рациональное число ---
Exact(x) = Sym(rationalize(Float64(x)))

# --- производная по времени выражения, зависящего от q(t) ---
d_dt(e) = diff(e, y_p) * dy_p + diff(e, l_p) * dl_p

####################################################################################################

# --- массы нижнего столба, поршня, верхнего столба (сохраняются) ---
function M()
    y0_L = Exact(y0_piston) - Exact(piston_width) / 2
    y0_R = Exact(y0_piston) + Exact(piston_width) / 2
    return y0_L * Exact(rho2A), (y0_R - y0_L) * Exact(rho1B), (Exact(Ly) - y0_R) * Exact(rho2A)
end

# --- границы поршня ---
function Y()
    return [y_p - l_p / 2, y_p + l_p / 2]
end

# --- плотности на трех участках ---
function Rho()
    M_L, M_P, M_R = M()
    Y_1, Y_2 = Y()
    drho_dy = -Exact(g) / (2 * Exact(Apsi))
    return [M_L / Y_1 + drho_dy * (y - Y_1 / 2),
            M_P / (Y_2 - Y_1),
            M_R / (Exact(Ly) - Y_2) + drho_dy * (y - (Y_2 + Exact(Ly)) / 2)]
end

# --- скорость на трех участках ---
function U()
    Y_1, Y_2 = Y()
    dY_1 = d_dt(Y_1)
    dY_2 = d_dt(Y_2)
    rho_L, rho_P, rho_R = Rho()
    return [-integrate(subs(d_dt(rho_L), y => z), (z, 0, y)) / rho_L,
            (rho_P * dY_1 - integrate(subs(d_dt(rho_P), y => z), (z, Y_1, y))) / rho_P,
            integrate(subs(d_dt(rho_R), y => z), (z, y, Exact(Ly))) / rho_R]
end

# --- Psi = ∫ psi0 dy, в поршне E_yy = (1 - (dy0/dy)^2)/2, dy0/dy = l_p0/l_p ---
function Psi()
    rho_L, rho_P, rho_R = Rho()
    E_yy = (1 - (Exact(piston_width) / l_p)^2) / 2
    return IntegrateRegions([Exact(Apsi) * (rho_L - Exact(rho2A))^2,
                             Exact(Bpsi) * (rho_P - Exact(rho1B))^2 + (Exact(muel) + Exact(lamel) / 2) * E_yy^2,
                             Exact(Apsi) * (rho_R - Exact(rho2A))^2])
end

####################################################################################################

# --- ∫_a^b e dy на участке r в явном виде ---
# В столбе все зависит от q только через его подвижную границу Y_b (нижний: y_p = Y_b + l_p/2,
# верхний: y_p = Y_b - l_p/2). После этой замены e - рациональная функция y вида P(y) / (c (y - y_0)^n)
# (знаменатель - степень плотности, линейной по y), она интегрируется подстановкой y = y_0 + w почленно:
# ∫ w^k dw = w^(k+1)/(k+1), ∫ dw/w = log(w). В поршне и при постоянной плотности e - многочлен по y.
function IntegrateY(e, a, b, r)
    if r == 1
        to, back = y_p => Y_b + l_p / 2, Y_b => y_p - l_p / 2
    elseif r == 3
        to, back = y_p => Y_b - l_p / 2, Y_b => y_p + l_p / 2
    else
        to, back = y_p => y_p, Y_b => Y_b
    end
    e, a, b = subs(e, to), subs(a, to), subs(b, to)

    num, den = sympy.fraction(sympy.cancel(sympy.together(e)))
    P_den = sympy.Poly(den, y)
    n = Int(P_den.degree())
    if n == 0
        F = integrate(sympy.expand(num / den), y)
        return subs(subs(F, y => b) - subs(F, y => a), back)
    end

    c   = P_den.all_coeffs()
    y_0 = -c[2] / (n * c[1])
    sympy.cancel(den - c[1] * (y - y_0)^n) == 0 || error("знаменатель не является степенью линейного множителя")

    P_num = sympy.Poly(sympy.expand(subs(num, y => y_0 + w)), w)
    w_a = a - y_0
    w_b = b - y_0
    res = Sym(0)
    for (m, c_m) in enumerate(reverse(P_num.all_coeffs()))
        k = m - 1 - n
        res += c_m * (k == -1 ? log(w_b / w_a) : (w_b^Sym(k + 1) - w_a^Sym(k + 1)) / (k + 1))
    end
    return subs(res / c[1], back)
end

# --- ∫_0^Ly e dy по трем участкам, e = [e_1, e_2, e_3] ---
function IntegrateRegions(e)
    b = [Sym(0), Y()..., Exact(Ly)]
    return sum(IntegrateY(e[r], b[r], b[r + 1], r) for r in 1:3)
end

####################################################################################################

function q_i()
    return q
end

function f_i()
    return [[diff(u, dq[i]) for i in 1:2] for u in U()]
end

function u_i()
    return [[f[i] * dq[i] for i in 1:2] for f in f_i()]
end

function M_ik()
    rho = Rho()
    f   = f_i()
    return [IntegrateRegions([rho[r] * f[r][i] * f[r][k] for r in 1:3]) for i in 1:2, k in 1:2]
end

function C_ijk()
    rho = Rho()
    f   = f_i()
    return [IntegrateRegions([rho[r] * f[r][k] * (diff(f[r][i], q[j]) + f[r][j] * diff(f[r][i], y)) for r in 1:3])
            for i in 1:2, j in 1:2, k in 1:2]
end

function D_ik()
    f = f_i()
    return [(Exact(4)/3*Exact(eta)+Exact(zeta))*IntegrateRegions([diff(f[r][i], y) * diff(f[r][k], y) for r in 1:3]) for i in 1:2, k in 1:2]
end

function F_k()
    rho = Rho()
    f   = f_i()
    return [Exact(g)*IntegrateRegions([rho[r] * f[r][k] for r in 1:3]) for k in 1:2]
end

function dΨdq_k()
    P = Psi()
    return [diff(P, q[k]) for k in 1:2]
end

####################################################################################################

# --- перевод символьного массива в числовую функцию s -> массив чисел Float64 ---
# (точные рациональные коэффициенты с большими числителями lambdify превращает в BigInt/BigFloat)
function ToNum(ex)
    F = [lambdify(e, [y_p, l_p, dy_p, dl_p]) for e in ex]
    return s -> map(Fe -> Float64(Fe(s...)), F)
end

Y_n = ToNum(Y())

# --- числовые функции для правой части ОДУ одного профиля плотности ---
function Numeric()
    return (M_ik   = ToNum(M_ik()),
            C_ijk  = ToNum(C_ijk()),
            D_ik   = ToNum(D_ik()),
            F_k    = ToNum(F_k()),
            dΨdq_k = ToNum(dΨdq_k()))
end

# --- числовые плотность и скорость на трех участках одного профиля (для графиков) ---
function NumericFields()
    return (Rho = [lambdify(r, [y, y_p, l_p, dy_p, dl_p]) for r in Rho()],
            U   = [lambdify(u, [y, y_p, l_p, dy_p, dl_p]) for u in U()])
end

####################################################################################################

# --- p - числовые функции профиля плотности из Numeric ---
function DiffEquation(ds, s, p, t)
    v   = [s[3], s[4]]
    C   = p.C_ijk(s)
    Cvv = [sum(C[i, j, k] * v[i] * v[j] for i in 1:2, j in 1:2) for k in 1:2]
    a   = 1.5*p.M_ik(s) \ (-Cvv - p.dΨdq_k(s) - 6*(p.D_ik(s) * v) - p.F_k(s))
    ds[1] = s[3]
    ds[2] = s[4]
    ds[3] = a[1]
    ds[4] = a[2]
end

s0       = [y0_piston, piston_width, v0_piston, 0.0]
analytic = Dict()
prob = ODEProblem(DiffEquation, s0, (0.0, t_max), Numeric())
analytic = (fields = NumericFields(), sol = solve(prob, reltol=1e-8, abstol=1e-10))


####################################################################################################

# --- для MakeFigures.jl ---
function AnalyticState(t)
    s = analytic.sol(t)
    return (y_p = s[1], l_p = s[2])
end

# --- номер участка, в котором лежит точка yy ---
function Region(yy, s)
    Y_1, Y_2 = Y_n(s)
    return yy <= Y_1 ? 1 : yy >= Y_2 ? 3 : 2
end

# --- ansatz 2: rho1 = 0 в столбах (фаза A), rho2 = 0 в поршне (фаза B) ---
function AnalyticProfilesY(yc_arr, y_arr, t)
    s    = analytic.sol(t)
    fld  = analytic.fields
    reg  = [Region(yy, s) for yy in yc_arr]
    rho  = [fld.Rho[reg[i]](yc_arr[i], s...) for i in eachindex(yc_arr)]
    rho1 = [reg[i] == 2 ? rho[i] : 0.0 for i in eachindex(yc_arr)]
    rho2 = [reg[i] == 2 ? 0.0 : rho[i] for i in eachindex(yc_arr)]
    v    = [fld.U[Region(yy, s)](yy, s...) for yy in y_arr]
    return (rho1 = rho1, rho2 = rho2, v = v)
end