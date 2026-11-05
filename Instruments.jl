#=
Файл со вспомогательными функциями, не связанными непосредственно с численным решением.
Включен в Main.jl через include.
=#

# --- конфиг-файлы лежат в папке Configs, в командной строке передается только имя файла ---
function DownloadConfig()
    global config_file = basename(ARGS[1])
    global config_path = abspath(joinpath("Configs", config_file))
    if !isfile(config_path)
        error("Не найден конфиг-файл: $config_path")
    end
    println("Загружаю конфигурацию: $config_path")
    println("")
    include(config_path)
end

# --- data_pos: позиция конца data.jls на момент чекпоинта (для отбрасывания хвоста при рестарте) ---
function SaveCheckpoint(rho1, rho2, rho, u, E, t, step, frame_count, next_frame_t,
                        energy_hist, time_hist, init_mass1, init_mass2, data_pos)
    checkpoint = (rho1 = rho1, rho2 = rho2, rho = rho, u = u, E = E,
                  t = t, step = step, frame_count = frame_count, next_frame_t = next_frame_t,
                  energy_hist = energy_hist, time_hist = time_hist,
                  init_mass1 = init_mass1, init_mass2 = init_mass2,
                  data_pos = data_pos)
    tmp_path = joinpath(direct, "checkpoint.jls.tmp")
    open(tmp_path, "w") do io
        serialize(io, checkpoint)
    end
    mv(tmp_path, joinpath(direct, "checkpoint.jls"), force = true)
    return nothing
end

function LoadCheckpoint()
    open(joinpath(direct, "checkpoint.jls"), "r") do io
        return deserialize(io)
    end
end