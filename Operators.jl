##########
# s, d, s*, d*
##########

function s(array)

    s_array = zeros(N + 2 * fict - 1)

    for i in (fict + 1):(N + fict - 1)
        s_array[i] = 0.5 * (array[i + 1] + array[i])
    end

    return s_array
end


function d(array)

    d_array = zeros(N + 2 * fict - 1)

    for i in (fict + 1):(N + fict - 1)
        d_array[i] = (1 / h) * (array[i + 1] - array[i])
    end

    return d_array
end


function sStar(array)

    sStar_array = zeros(N + 2 * fict)

    for i in (fict + 1):(N + fict)
        sStar_array[i] = 0.5 * (array[i] + array[i - 1])
    end

    return sStar_array
end


function dStar(array)

    dStar_array = zeros(N + 2 * fict)

    for i in (fict + 1):(N + fict)
        dStar_array[i] = (1 / h) * (array[i] - array[i - 1])
    end

    return dStar_array
end

