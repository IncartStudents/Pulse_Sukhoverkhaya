using Plots
using DataFrames

include("../ECG markup/detector_funcs.jl")

function LeadMarkup(sig, fs)
    # определение полярности и инвертирование
    # level = minimum(sig)+(maximum(sig)+abs(minimum(sig)))/2
    # if level < 0 sig *= -1 end

    # Предподготовка (алгоритм Пана-Томпкинса)
    filtered, delay_1 = lynn_filter(sig, "bandpass") 
    differed, delay_2 = fivepointdiff(filtered)

    # Интегрирование квадрата сигнала в скользящем окне шириной 95 мс (фильтр скользящего среднего)
    # Возведение в квадрат
    sqred = differed.^2
    integrated, delay_3 = movingaverage(sqred, 0.096, fs)

    # Поиск всех реперных точек (изменение направления интегрированного сигнала с возрастания на убывание)
    maxpos = findmax(integrated, 0.2, fs)

    # Поиск точек пересечения нуля дифференцированным сигналом
    zerocrosses = zerocross(differed)

    # Идентификация QRS
    Q, R, S = qrs_points(filtered,zerocrosses,maxpos,fs)

    # Поиск P
    DERFI, delay_3 = lynn_filter(differed, "low")
    P = p_points(DERFI, sig, R, Q, fs)

    # Поиск Т
    # Tud, Tdu, Tod, Tou = t_points(DERFI, R, fs);
    # T = vcat(Tud, Tdu, Tod, Tou);
    T = t_points(DERFI, R, fs)

    # Определение границ волн
    P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te = findwavebounds(differed, DERFI, filtered, P, Q, R, S, T)

    # Коррекция позиций 
    P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te = pos_correct(delay_1, delay_2, delay_3, P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te)

    return P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te
end