#=
Главный файл проекта, который отвечает за исполнения всей программы целиком.

Запуск: julia Main.jl <config>.jl
Результат: численное решение в двоичном коде в файле "data.jls"
В случае необходимости проводить расчеты не с нулевого момента: restart = true в конфиге (с последнего сохраненного чекпоинта)
Логирование: прогресс и проверка отсутствия взрыва решения, сохранения массы, невозрастания полной энергии в файле "log.txt"
=#

include("Imports.jl")
include("Instruments.jl")
DownloadConfig()
include("Grid.jl")
include("Initialization.jl")
include("BoundaryConditions.jl")
include("Operators.jl")
include("Thermodynamics.jl")
include("Stresses.jl")
include("Solver.jl")
include("Plotting.jl")

function main(io)
    # --- генерация сетки ---
    x, xc = CreateGrid()

    # --- продолжаем расчет, только если это запрошено и есть сохраненное состояние ---
    resuming = restart && isfile(joinpath(direct, "checkpoint.jls"))

    if resuming
        # --- восстановление состояния из чекпоинта ---
        chk = LoadCheckpoint()

        rho1 = chk.rho1
        rho2 = chk.rho2
        rho  = chk.rho
        u    = chk.u
        E    = chk.E

        t            = chk.t
        step         = chk.step
        next_frame_t = chk.next_frame_t
        frame_count  = chk.frame_count

        energy_hist  = chk.energy_hist
        time_hist    = chk.time_hist

        # --- исходные массы берем из чекпоинта, чтобы проверять сохранение относительно начального состояния ---
        init_mass1   = chk.init_mass1
        init_mass2   = chk.init_mass2

        phi = Phi!(zeros(length(rho1)), xc, rho1)

        # --- дозапись фреймов в существующий файл результатов ---
        data_io = open(joinpath(direct, "data.jls"), "a")

        println(@sprintf("Продолжаю расчет с t = %.5f, шаг %d", t, step))
        flush(io)
    else
        # --- применение начальных условий ---
        u               = InitVelocity()
        rho1, rho2, rho = InitDensityChoise(xc)
        E               = InitAlmansiTens()

        phi             = Phi!(zeros(length(rho1)), xc, rho1)

        # --- применение граничных условий ---
        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho  = BoundCondScalar(rho)
        u    = BoundCondVector(u)
        E    = BoundCondScalar(E)

        # --- время, шаг по времени, время записи, шаг записи ---
        t            = 0.0
        step         = 0
        next_frame_t = frame_dt
        frame_count  = 0

        # --- массивы для функций времени ---
        energy_hist = []
        time_hist   = []

        # --- начальные массы веществ для проверки сохранения массы ---
        init_mass1 = sum(rho1[fict + 1:N + fict - 1]) * h
        init_mass2 = sum(rho2[fict + 1:N + fict - 1]) * h

        # --- новый файл для записи результатов расчета в двоичном коде ---
        data_io = open(joinpath(direct, "data.jls"), "w")
        serialize(data_io, (x = x, xc = xc))
    end

    # --- в памяти нужен только последний фрейм для сериализации ---
    frames_data = []

    # --- первый фрейм пишем только для нового расчета (при рестарте он уже в data.jls) ---
    if !resuming
        frame_count, mu_el, lam_el = MakeFrame(data_io, xc, x, rho1, rho2, rho, u, E, phi,
                                               t, step, frame_count, energy_hist, time_hist, frames_data)
        SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                       energy_hist, time_hist, init_mass1, init_mass2)
    end

    # --- создание временных массивов для хранения данных ---
    buf = CreateBuffers()

    # --- выделение памяти под массивы на следующем шаге по времени ---
    rho1_new = zeros(N + 2 * fict - 1)
    rho2_new = zeros(N + 2 * fict - 1)
    u_new    = [zeros(3) for _ in 1:N + 2 * fict]
    E_new    = [zeros(3, 3) for _ in 1:N + 2 * fict - 1]

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
        MuEl!(buf.mu_el, rho1, rho2)
        LamEl!(buf.lam_el, rho1, rho2)
        Mu1Hat!(buf.mu1h, buf, E, rho1, rho2)
        Mu2Hat!(buf.mu2h, buf, E, rho1, rho2)
        M!(buf.M_arr, M0, rho1, rho2)
        K!(buf.K_current, E, buf.mu_el, buf.lam_el)
        PiEl!(buf.Pi_el_current, E, buf.K_current)
        PiNS!(buf.Pi_ns_current, u, buf)

        # --- интегрирование по времени законов сохранения и заполнение соответствующих массивов ---
        ConservMassAndMomentumStep!(rho1_new, rho2_new, u_new, buf,
                                    rho1, rho2, u, E, phi, dt,
                                    buf.mu_el, buf.lam_el, buf.mu1h, buf.mu2h, buf.M_arr,
                                    buf.K_current, buf.Pi_el_current, buf.Pi_ns_current)
        ConservAlmansiTensor!(E_new, buf, E, u, dt)

        rho1, rho1_new = rho1_new, rho1
        rho2, rho2_new = rho2_new, rho2
        u, u_new       = u_new, u
        E, E_new       = E_new, E

        @inbounds for i in 1:N + 2 * fict - 1
            rho[i] = rho1[i] + rho2[i]
        end

        # --- применение граничных условий после интегрирования ---
        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho  = BoundCondScalar(rho)
        u    = BoundCondVector(u)
        E    = BoundCondScalar(E)

        t    += dt
        step += 1

        if t >= next_frame_t - 1e-10
            frame_count, mu_el, lam_el = MakeFrame(data_io, xc, x, rho1, rho2, rho, u, E, phi,
                                                   t, step, frame_count, energy_hist, time_hist, frames_data)

            next_frame_t += frame_dt

            # --- чекпоинт синхронно с фреймом ---
            SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                           energy_hist, time_hist, init_mass1, init_mass2)

            progress_pct = 100.0 * (t / t_max)
            println(@sprintf("Прогресс %.1f%% | шаг %d | время %.5f ", progress_pct, step, t))
            flush(io)

            # --- копирование лог-файла в соответствующую директорию ---
            cp("log.txt", joinpath(direct, "log.txt"), force=true)

            # --- проверка отсутствия взрывов решения, сохранения масс компонентов и невозрастания полной энергии ---
            if !CheckStability(rho1, rho2, u, step, t, E, mu_el, lam_el, phi)
                break
            end
            mass_ok   = CheckMassConservation(rho1, rho2, init_mass1, init_mass2, step, t, mass_rtol)
            energy_ok = CheckEnergyNonIncrease(energy_hist, step, t, energy_rtol)
            if !mass_ok || !energy_ok
                break
            end
        end
    end

    flush(io)
    cp("log.txt", joinpath(direct, "log.txt"), force=true)

    # --- финальный чекпоинт: точное состояние последнего посчитанного момента для продолжения ---
    SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                   energy_hist, time_hist, init_mass1, init_mass2)

    close(data_io)
end

# --- нужно для того, что бы в процессе расчета лог-файл заполнялся постепенно, а не один раз в самом конце ---
open("log.txt", "w") do io
    redirect_stdout(io) do
        main(io)
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