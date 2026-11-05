#=
Главный файл проекта, который отвечает за исполнения всей программы целиком.

Запуск: julia Main.jl <config>.jl   (конфиг ищется в папке Configs)
Размерность задачи задается в конфиге: const dim = 1 (одномерная вдоль y, параметры по x игнорируются) или 2.
Если dim в конфиге не задан, задача считается одномерной.
Двумерная задача распараллелена (Polyester): число потоков задается при запуске, julia -t 8 Main.jl <config>.jl.
Одномерная задача считается последовательно при любом числе потоков.
Результат: численное решение в двоичном коде в файле "data.jls"
В случае необходимости проводить расчеты не с нулевого момента: restart = true в конфиге (с последнего сохраненного чекпоинта)
Логирование: прогресс и проверка отсутствия взрыва решения, сохранения массы, невозрастания полной энергии в файле "log.txt"
Чекпоинт сохраняется только для состояний, прошедших проверки; при рестарте хвост data.jls после чекпоинта отбрасывается
=#

include("Imports.jl")
include("Instruments.jl")
DownloadConfig()

# --- по умолчанию задача одномерная ---
if !@isdefined(dim)
    const dim = 1
end
dim in (1, 2) || error("Недопустимое значение dim = $dim, допустимы 1 и 2")

include("Parallel.jl")
include("Grid.jl")
include("Initialization.jl")
include("BoundaryConditions.jl")
include("Operators.jl")
include("Workspace.jl")
include("Thermodynamics.jl")
include("Stresses.jl")
include("Solver.jl")
include("Verify.jl")
include("Frames.jl")

####################################################################################################
# --- двумерная задача ---

function Main2D(io)
    # --- генерация сетки и срезов для профилей во фреймах ---
    x, xc, y, yc = CreateGrid2D()
    slices       = FrameSlices(x, xc, y, yc)

    # --- продолжаем расчет, только если это запрошено и есть сохраненное состояние ---
    resuming = restart && isfile(joinpath(direct, "checkpoint.jls"))

    if resuming
        # --- восстановление состояния из чекпоинта ---
        chk = LoadCheckpoint()

        rho1 = chk.rho1
        rho2 = chk.rho2
        rho  = chk.rho
        u    = chk.u
        E    = TensorField(chk.E)   # чекпоинт старого формата хранит Matrix{Float64} в ячейках

        t            = chk.t
        step         = chk.step
        next_frame_t = chk.next_frame_t
        frame_count  = chk.frame_count

        energy_hist  = chk.energy_hist
        time_hist    = chk.time_hist

        # --- исходные массы берем из чекпоинта, чтобы проверять сохранение относительно начального состояния ---
        init_mass1   = chk.init_mass1
        init_mass2   = chk.init_mass2

        phi = Phi(yc, rho1)

        # --- дозапись фреймов в существующий файл результатов ---
        # --- фреймы, записанные после чекпоинта (например, фрейм, не прошедший проверку), отбрасываются ---
        data_io = open(joinpath(direct, "data.jls"), "r+")
        if haskey(chk, :data_pos)
            truncate(data_io, chk.data_pos)
        end
        seekend(data_io)

        println(@sprintf("Продолжаю расчет с t = %.5f, шаг %d", t, step))
        flush(io)
    else
        # --- применение начальных условий ---
        u               = InitVelocity2D()
        rho1, rho2, rho = InitDensityChoise(xc, yc)
        E               = InitAlmansiTens2D()

        phi             = Phi(yc, rho1)

        # --- применение граничных условий ---
        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho  = BoundCondScalar(rho)
        u    = BoundCondVector(u)
        E    = BoundCondTensor(E)

        # --- время, шаг по времени, время записи, шаг записи ---
        t            = 0.0
        step         = 0
        next_frame_t = frame_dt
        frame_count  = 0

        # --- массивы для функций времени ---
        energy_hist = []
        time_hist   = []

        # --- начальные массы веществ для проверки сохранения массы ---
        init_mass1 = sum(rho1[ficty + 1:Ny + ficty - 1, fictx + 1:Nx + fictx - 1]) * hx * hy
        init_mass2 = sum(rho2[ficty + 1:Ny + ficty - 1, fictx + 1:Nx + fictx - 1]) * hx * hy

        # --- новый файл для записи результатов расчета в двоичном коде ---
        data_io = open(joinpath(direct, "data.jls"), "w")
        serialize(data_io, FrameHeader(x, xc, y, yc, slices))
    end

    # --- в памяти нужен только последний фрейм предыдущего расчета для сериализации ---
    frames_data = []

    # --- первый фрейм пишем только для нового расчета (при рестарте он уже в data.jls) ---
    if !resuming
        frame_count, mu_el, lam_el = MakeFrame(data_io, slices, rho1, rho2, u, E, phi,
                                               t, step, frame_count, energy_hist, time_hist, frames_data)
        SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                       energy_hist, time_hist, init_mass1, init_mass2, position(data_io))
    end

    # --- время последнего корректного чекпоинта ---
    chk_t  = t
    failed = false

    # --- выделение памяти под массивы, которые перезаписываются на каждом шаге по времени ---
    Ncentr_y = Ny + 2 * ficty - 1
    Ncentr_x = Nx + 2 * fictx - 1

    mu_el         = zeros(Ncentr_y, Ncentr_x)
    lam_el        = zeros(Ncentr_y, Ncentr_x)
    mu1h          = zeros(Ncentr_y, Ncentr_x)
    mu2h          = zeros(Ncentr_y, Ncentr_x)
    M_arr         = zeros(Ncentr_y, Ncentr_x)
    K_current     = fill(zero(Tensor3), Ncentr_y, Ncentr_x)
    Pi_el_current = fill(zero(Tensor3), Ncentr_y, Ncentr_x)
    Pi_ns_current = fill(zero(Tensor3), Ncentr_y, Ncentr_x)

    # --- рабочие массивы операторов и законов сохранения ---
    ws = Workspace2D()

    # --- массивы на следующем шаге по времени ---
    rho1_new = zeros(Ncentr_y, Ncentr_x)
    rho2_new = zeros(Ncentr_y, Ncentr_x)
    u_new    = [zeros(size(u[1])), zeros(size(u[2]))]
    E_new    = fill(zero(Tensor3), Ncentr_y, Ncentr_x)

    # --- реальное время записи предыдущего фрейма (для вывода в лог) ---
    frame_wall_t = time()

    # --- полный цикл по времени ---
    while step < step_max
        # --- проверки на конец расчета ---
        if t >= t_max
            break
        end

        dt = delta_t
        if t + dt > t_max
            dt = t_max - t
        end

        # --- заполнение необходимых массивов ---
        MuEl!(mu_el, rho1, rho2)
        LamEl!(lam_el, rho1, rho2)
        Mu1Hat!(mu1h, E, rho1, rho2, ws)
        Mu2Hat!(mu2h, E, rho1, rho2, ws)
        M!(M_arr, M0, rho1, rho2)
        K!(K_current, E, mu_el, lam_el)
        PiEl!(Pi_el_current, E, K_current)
        PiNS!(Pi_ns_current, u, ws)

        # --- интегрирование по времени законов сохранения и заполнение соответствующих массивов ---
        ConservMass!(rho1_new, rho2_new, rho1, rho2, u, dt, mu1h, mu2h, M_arr, ws)
        ConservMomentum!(u_new, rho1_new, rho2_new, rho1, rho2, u, E, phi, dt,
                         mu1h, mu2h, K_current, Pi_el_current, Pi_ns_current, ws)
        ConservAlmansiTensor!(E_new, E, u, dt, ws)

        rho1, rho1_new = rho1_new, rho1
        rho2, rho2_new = rho2_new, rho2
        u, u_new       = u_new, u
        E, E_new       = E_new, E

        Map2D!(+, rho, rho1, rho2)

        # --- применение граничных условий после интегрирования ---
        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho  = BoundCondScalar(rho)
        u    = BoundCondVector(u)
        E    = BoundCondTensor(E)

        t    += dt
        step += 1

        if t >= next_frame_t - 1e-10
            # --- фрейм пишется и для состояния, не прошедшего проверку (для диагностики) ---
            frame_count, mu_el, lam_el = MakeFrame(data_io, slices, rho1, rho2, u, E, phi,
                                                   t, step, frame_count, energy_hist, time_hist, frames_data)

            next_frame_t += frame_dt

            progress_pct = 100.0 * (t / t_max)
            frame_wall   = time() - frame_wall_t
            frame_wall_t = time()
            println(@sprintf("Прогресс %.1f%% | шаг %d | время %.5f | %.1f с с прошлого фрейма",
                             progress_pct, step, t, frame_wall))
            flush(io)

            # --- копирование лог-файла в соответствующую директорию ---
            cp(log_file, joinpath(direct, log_file), force=true)

            # --- проверка отсутствия взрывов решения, сохранения масс компонентов и невозрастания полной энергии ---
            if !CheckStability(rho1, rho2, u, step, t, E, mu_el, lam_el, phi)
                failed = true
                break
            end
            mass_ok   = CheckMassConservation(rho1, rho2, init_mass1, init_mass2, step, t, mass_rtol)
            energy_ok = CheckEnergyNonIncrease(energy_hist, step, t, energy_rtol)
            if !mass_ok || !energy_ok
                failed = true
                break
            end

            # --- чекпоинт синхронно с фреймом, только после успешных проверок ---
            SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                           energy_hist, time_hist, init_mass1, init_mass2, position(data_io))
            chk_t = t
        end
    end

    # --- состояние последнего посчитанного момента для продолжения (только если проверки не провалены) ---
    if failed
        println(@sprintf("Расчет остановлен проверкой при t = %.5f. Последний корректный чекпоинт: t = %.5f", t, chk_t))
    else
        SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                       energy_hist, time_hist, init_mass1, init_mass2, position(data_io))
    end

    flush(io)
    cp(log_file, joinpath(direct, log_file), force=true)

    close(data_io)
end

####################################################################################################
# --- одномерная задача (вдоль y): u = 0, d/dx = 0, все поля зависят только от y ---

function Main1D(io)
    # --- генерация сетки и внутренних координат для фреймов ---
    y, yc  = CreateGrid1D()
    slices = FrameSlices(y, yc)

    # --- продолжаем расчет, только если это запрошено и есть сохраненное состояние ---
    resuming = restart && isfile(joinpath(direct, "checkpoint.jls"))

    if resuming
        # --- восстановление состояния из чекпоинта ---
        chk = LoadCheckpoint()

        rho1 = chk.rho1
        rho2 = chk.rho2
        rho  = chk.rho
        u    = chk.u
        E    = TensorField(chk.E)   # чекпоинт старого формата хранит Matrix{Float64} в ячейках

        t            = chk.t
        step         = chk.step
        next_frame_t = chk.next_frame_t
        frame_count  = chk.frame_count

        energy_hist  = chk.energy_hist
        time_hist    = chk.time_hist

        # --- исходные массы берем из чекпоинта, чтобы проверять сохранение относительно начального состояния ---
        init_mass1   = chk.init_mass1
        init_mass2   = chk.init_mass2

        phi = Phi(yc, rho1)

        # --- дозапись фреймов в существующий файл результатов ---
        # --- фреймы, записанные после чекпоинта (например, фрейм, не прошедший проверку), отбрасываются ---
        data_io = open(joinpath(direct, "data.jls"), "r+")
        if haskey(chk, :data_pos)
            truncate(data_io, chk.data_pos)
        end
        seekend(data_io)

        println(@sprintf("Продолжаю расчет с t = %.5f, шаг %d", t, step))
        flush(io)
    else
        # --- применение начальных условий ---
        u               = InitVelocity1D()
        rho1, rho2, rho = InitDensityChoise(yc)
        E               = InitAlmansiTens1D()

        phi             = Phi(yc, rho1)

        # --- применение граничных условий ---
        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho  = BoundCondScalar(rho)
        u    = BoundCondVector(u)
        E    = BoundCondTensor(E)

        # --- время, шаг по времени, время записи, шаг записи ---
        t            = 0.0
        step         = 0
        next_frame_t = frame_dt
        frame_count  = 0

        # --- массивы для функций времени ---
        energy_hist = []
        time_hist   = []

        # --- начальные массы веществ для проверки сохранения массы ---
        init_mass1 = sum(rho1[ficty + 1:Ny + ficty - 1]) * hy
        init_mass2 = sum(rho2[ficty + 1:Ny + ficty - 1]) * hy

        # --- новый файл для записи результатов расчета в двоичном коде ---
        data_io = open(joinpath(direct, "data.jls"), "w")
        serialize(data_io, FrameHeader(y, yc, slices))
    end

    # --- в памяти нужен только последний фрейм предыдущего расчета для сериализации ---
    frames_data = []

    # --- первый фрейм пишем только для нового расчета (при рестарте он уже в data.jls) ---
    if !resuming
        frame_count, mu_el, lam_el = MakeFrame(data_io, slices, rho1, rho2, u, E, phi,
                                               t, step, frame_count, energy_hist, time_hist, frames_data)
        SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                       energy_hist, time_hist, init_mass1, init_mass2, position(data_io))
    end

    # --- время последнего корректного чекпоинта ---
    chk_t  = t
    failed = false

    # --- выделение памяти под массивы, которые перезаписываются на каждом шаге по времени ---
    Ncentr_y = Ny + 2 * ficty - 1

    mu_el         = zeros(Ncentr_y)
    lam_el        = zeros(Ncentr_y)
    mu1h          = zeros(Ncentr_y)
    mu2h          = zeros(Ncentr_y)
    M_arr         = zeros(Ncentr_y)
    K_current     = fill(zero(Tensor3), Ncentr_y)
    Pi_el_current = fill(zero(Tensor3), Ncentr_y)
    Pi_ns_current = fill(zero(Tensor3), Ncentr_y)

    # --- рабочие массивы операторов и законов сохранения ---
    ws = Workspace1D()

    # --- массивы на следующем шаге по времени ---
    rho1_new = zeros(Ncentr_y)
    rho2_new = zeros(Ncentr_y)
    u_new    = zeros(size(u))
    E_new    = fill(zero(Tensor3), Ncentr_y)

    # --- реальное время записи предыдущего фрейма (для вывода в лог) ---
    frame_wall_t = time()

    # --- полный цикл по времени ---
    while step < step_max
        # --- проверки на конец расчета ---
        if t >= t_max
            break
        end

        dt = delta_t
        if t + dt > t_max
            dt = t_max - t
        end

        # --- заполнение необходимых массивов ---
        MuEl!(mu_el, rho1, rho2)
        LamEl!(lam_el, rho1, rho2)
        Mu1Hat!(mu1h, E, rho1, rho2, ws)
        Mu2Hat!(mu2h, E, rho1, rho2, ws)
        M!(M_arr, M0, rho1, rho2)
        K!(K_current, E, mu_el, lam_el)
        PiEl!(Pi_el_current, E, K_current)
        PiNS!(Pi_ns_current, u, ws)

        # --- интегрирование по времени законов сохранения и заполнение соответствующих массивов ---
        ConservMass!(rho1_new, rho2_new, rho1, rho2, u, dt, mu1h, mu2h, M_arr, ws)
        ConservMomentum!(u_new, rho1_new, rho2_new, rho1, rho2, u, E, phi, dt,
                         mu1h, mu2h, K_current, Pi_el_current, Pi_ns_current, ws)
        ConservAlmansiTensor!(E_new, E, u, dt, ws)

        rho1, rho1_new = rho1_new, rho1
        rho2, rho2_new = rho2_new, rho2
        u, u_new       = u_new, u
        E, E_new       = E_new, E

        @inbounds for i in eachindex(rho)
            rho[i] = rho1[i] + rho2[i]
        end

        # --- применение граничных условий после интегрирования ---
        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho  = BoundCondScalar(rho)
        u    = BoundCondVector(u)
        E    = BoundCondTensor(E)

        t    += dt
        step += 1

        if t >= next_frame_t - 1e-10
            # --- фрейм пишется и для состояния, не прошедшего проверку (для диагностики) ---
            frame_count, mu_el, lam_el = MakeFrame(data_io, slices, rho1, rho2, u, E, phi,
                                                   t, step, frame_count, energy_hist, time_hist, frames_data)

            next_frame_t += frame_dt

            progress_pct = 100.0 * (t / t_max)
            frame_wall   = time() - frame_wall_t
            frame_wall_t = time()
            println(@sprintf("Прогресс %.1f%% | шаг %d | время %.5f | %.1f с с прошлого фрейма",
                             progress_pct, step, t, frame_wall))
            flush(io)

            # --- копирование лог-файла в соответствующую директорию ---
            cp(log_file, joinpath(direct, log_file), force=true)

            # --- проверка отсутствия взрывов решения, сохранения масс компонентов и невозрастания полной энергии ---
            if !CheckStability(rho1, rho2, u, step, t, E, mu_el, lam_el, phi)
                failed = true
                break
            end
            mass_ok   = CheckMassConservation(rho1, rho2, init_mass1, init_mass2, step, t, mass_rtol)
            energy_ok = CheckEnergyNonIncrease(energy_hist, step, t, energy_rtol)
            if !mass_ok || !energy_ok
                failed = true
                break
            end

            # --- чекпоинт синхронно с фреймом, только после успешных проверок ---
            SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                           energy_hist, time_hist, init_mass1, init_mass2, position(data_io))
            chk_t = t
        end
    end

    # --- состояние последнего посчитанного момента для продолжения (только если проверки не провалены) ---
    if failed
        println(@sprintf("Расчет остановлен проверкой при t = %.5f. Последний корректный чекпоинт: t = %.5f", t, chk_t))
    else
        SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                       energy_hist, time_hist, init_mass1, init_mass2, position(data_io))
    end

    flush(io)
    cp(log_file, joinpath(direct, log_file), force=true)

    close(data_io)
end

####################################################################################################

# --- нужно для того, что бы в процессе расчета лог-файл заполнялся постепенно, а не один раз в самом конце ---
# --- предупреждения @warn тоже пишутся в лог-файл ---
open(log_file, "w") do io
    redirect_stdout(io) do
        with_logger(ConsoleLogger(io)) do
            if dim == 1
                Main1D(io)
            else
                Main2D(io)
            end
        end
    end
end

# --- сообщение о завершении расчета и инструкция по построению графиков ---
println("Расчет завершен.")
println("Данные сохранены в двоичном коде в файл: ", joinpath(direct, "data.jls"))
println("Данные чекпоинта сохранены в двоичном коде в файл: ", joinpath(direct, "checkpoint.jls"))
println("Чтобы продолжить расчет, увеличте t_max, а также измените значение restart на true в конфиг-файле и заново выполните julia Main.jl ", ARGS[1])
println("")
println("Для построения графиков выполните: julia MakeFigures.jl ", ARGS[1])
println("Графики будут сохранены в директорию: ", direct)
