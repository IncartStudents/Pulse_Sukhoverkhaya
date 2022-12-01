using DSP
using Statistics

include("../src/readfiles.jl")

struct ZeroCross
    pos::Int64
    type::Int64
end

struct Wave
    b::Int64
    e::Int64
end

# получение из бинарника одного канала и частоты дискретизации
function sigdata(binfile)
    signals, fs, timestart, units = readbin(binfile);
    fs = Int(fs);

    return signals, fs
end

# ФНЧ до 11 (15?) Гц
# y(nT) = 2y(nT - T) - y(nT - 2 T) + x(nT) - 2x(nT- 6T) + x(nT- 12T)
function my_lowpass(signal)

    fsignal = similar(signal)
    for i in 1:lastindex(signal)
        if i < 13
            fsignal[i] = signal[i]
        else
            fsignal[i] = 2*fsignal[i-1] - fsignal[i-2] + signal[i] - 2*signal[i-6] + signal[i-12]
        end
    end

    # aprox = range(fsignal[1], fsignal[length(fsignal)], length(fsignal))
    
    # fsignal = fsignal .- aprox

    delay = 6
    
    return fsignal, delay
end

# ФВЧ от 5 Гц
function my_highpass(signal)

    fsignal = similar(signal)
    for i in 1:lastindex(signal)
        if i < 33
            fsignal[i] = signal[i]
        else
            fsignal[i] = 32*signal[i-16] - (fsignal[i-1] + signal[i] - signal[i-32])
        end
    end

    delay = 16
    
    return fsignal, delay
end

# Полосовой фильтр 5 - 11 (15?) Гц
function lynn_filter(signal, type)

    delay_low = 0; delay_high = 0;

    if type == "low"
        filtered, delay_low = my_lowpass(signal)
    elseif type == "bandpass"
        filtered_low, delay_low = my_lowpass(signal)
        filtered, delay_high = my_highpass(filtered_low)
    end

    aprox = range(filtered[1], filtered[end], length(filtered))
    filtered = filtered .- aprox

    delay = delay_low + delay_high

    return filtered, delay
end

# Дифференцирующий ФНЧ (пятиточечная производная)
# Разностное уравнение: y(nT) = (1/8 T) [-x(nT - 2 T) - 2x(nT - T) + 2x(nT + T) +x(nT+ 2T)]
function fivepointdiff(filtered)

    differed = fill(0.0, length(filtered))

    for i in 3:length(filtered)-2
        differed[i] = (1/8)*(1)*(-filtered[i-2] - 2*filtered[i-1] + 2*filtered[i+1] + filtered[i+2])
    end

    # delay = 2

    delay = 1

    return differed, delay
end

# Фильтр скользящего среднего
function movingaverage(sqred, wsize, fs)

    w_size = wsize*fs |> Int             # ширина окна в отсчетах
    window = sqred[1:w_size]             # (окно) массив w_size последних значений сигнала
    w_cnt = 1
    w_sum = sqred[1]*w_size              # сумма значений в окне

    integrated = fill(0.0, length(sqred))

    for i in 1:lastindex(sqred)

        last = window[w_cnt]              # самое старое значение в окне (будет выкинуто)
        window[w_cnt] = sqred[i]          # запись нового значения в окне (на место самого старого)

        w_sum = w_sum - last + sqred[i]   # сумма значений в окне
        if i>w_size/2 integrated[i] =  w_sum/w_size  end   # среднее в окне

        w_cnt += 1
        if w_cnt > w_size w_cnt = 1 end   # переключение счетчика при необходимости

    end

    delay = w_size/2 |> Int;

    return integrated, delay
end

# Поиск локальных максимумов
function findmax(integrated, radius, fs)

    lastmaxpos = 0
    lastmaxamp = 0
    radius = radius*fs |> Int

    maxpos = Int[]

    SignalLevel1 = 0
    NoiseLevel1 = 0
    Threshold1 = 0

    init = Float64[]

    for i in 1:lastindex(integrated)

        # 2-секундная фаза обучения
        if i/fs < 2
            push!(init, integrated[i])
        elseif i/fs == 2
            push!(init, integrated[i])
            max = maximum(init);
            mn = mean(init);

            SignalLevel1 = 0.875*max;
            NoiseLevel1 = 0.875*mn;
        else
            # Подтверждение предыдущего предполагаемого максимума
            if i-lastmaxpos > radius

                if lastmaxamp < Threshold1
                    NoiseLevel1 = 0.125*lastmaxamp+0.875*NoiseLevel1;
                else
                    push!(maxpos, lastmaxpos)
                
                    lastmaxamp = integrated[i]
                    lastmaxpos = i
                end
            end

            if integrated[i] >= lastmaxamp
                lastmaxamp = integrated[i]
                lastmaxpos = i
            end

        end

        Threshold1 = NoiseLevel1 + 0.25(SignalLevel1 - NoiseLevel1);

    end

    return maxpos
end

# Поиск точек пересечения нуля
function zerocross(differed)
    
    zerocrosses = ZeroCross[]

    for i in 2:lastindex(differed)
        if differed[i] >= 0 && differed[i-1] < 0
            push!(zerocrosses, ZeroCross(i,1))
        elseif differed[i] <= 0 && differed[i-1] > 0
            push!(zerocrosses, ZeroCross(i,0))
        end
    end

    return zerocrosses
end

# Идентификация R, Q, S
function qrs_points(filtered, zerocrosses, maxpos, fs)

    R = Int[]; 
    Q = Int[]; 
    S = Int[];

    Lvl = maximum(filtered)*0.3

    # zc = hcat(zerocrosses.pos, zerocrosses.type);

    for i in 1:lastindex(maxpos)

        # исходим из того, что пик интегрированного всегда правее qrs (НЕ ФАКТ)
        ind = map(x -> x.pos > maxpos[i]-0.5*fs && x.pos < maxpos[i] && x.pos > 1, zerocrosses);
        sample = zerocrosses[ind]
        n = length(sample)

        isq = false;
        isr = false;
        iss = false;

        for j in n:-1:1
            if !isq || !isr || !iss
                if sample[j].type == 0 && filtered[sample[j].pos]>Lvl
                    if !isr
                        push!(R, sample[j].pos)
                        isr = true
                    end
                elseif sample[j].type == 1 
                    if !iss && !isr && abs(filtered[sample[j].pos])>100
                        push!(S, sample[j].pos)
                        iss = true
                    elseif !isq && isr
                        push!(Q, sample[j].pos)
                        isq = true
                    end
                end
            end
        end

        if !isq push!(Q, 0) end
        if !isr push!(R, 0) end
        if !iss push!(S, 0) end

    end

    return Q, R, S
end

# Поиск пика волны P
function p_points(DERFI, sig, R, Q, fs)

    ws = 0.156*fs |> Int64;
    preR = 0.226*fs |> Int64;

    P = Int[];

    k = 0
    for i in R
        k += 1
        if i != 0
            if i-preR > 0 
                b1 = i-preR;
                b2 = i-preR+ws;

                if Q[k] != 0
                    if b2 > Q[k] b2 = Q[k]-1 end
                end

                if mean(sig[b1:b2]) >= 0
                    bend = b2
                    for i in 1:(b2-b1)
                        wmax = argmax(DERFI[b1:bend])
                        if wmax+b1 == bend+1
                            bend = b2 - i
                        else
                            break
                        end
                    end
                    wmin = argmin(DERFI[wmax+b1:b2])

                    wbeg = wmax+b1
                    wend = wmin+wmax+b1

                    isinv = false
                else
                    bend = b2
                    for i in 1:(b2-b1)
                        wmin = argmin(DERFI[b1:bend])
                        if wmin+b1 == bend+1
                            bend = b2 - i
                        else
                            break
                        end
                    end
                    wmax = argmax(DERFI[wmin+b1:b2])

                    wbeg = wmin+b1
                    wend = wmin+wmax+b1

                    isinv = true
                end

                if wend > length(DERFI) wend = length(DERFI) end

                wsign = DERFI[range(wbeg,wend)]

                for j in 2:lastindex(wsign)
                    if isinv
                        if wsign[j] >= 0 && wsign[j-1] < 0
                            push!(P, j+wmax+b1)
                            break
                        end
                    else
                        if wsign[j] <= 0 && wsign[j-1] > 0
                            push!(P, j+wmax+b1)
                            break
                        end
                    end
                end
            end
        end

        if length(P) < k push!(P, 0) end
    end

    return P
end

# Поиск пика волны Т
function t_points(DERFI, R, fs)
    # Рассчет среднего значения R-R-интервалов
    sumRR = 0
    r = R[findall(x -> x!=0, R)]
    if length(r) == 1
        RRav = r[1]
    else
        for i in 2:lastindex(r) sumRR += r[i] - r[i-1] end
        RRav = abs(sumRR/(length(r)-1))
    end

    # Выбор границ окна поиска
    if RRav > 0.7*fs
        bwind = 0.14*fs |> Int64
        ewind = 0.5*fs  |> Int64
    else
        bwind = 0.1*fs |> Int64
        ewind = round(0.6*RRav)  |> Int64
    end

    # Tud = []
    # Tdu = []
    # Tod = []
    # Tou = []
    T = Int[];

    for i in 1:lastindex(R)
        if R[i] != 0 && (bwind+R[i]) < length(DERFI) && (ewind+R[i]) < length(DERFI)

            window = range(bwind+R[i], ewind+R[i])

            dw0 = DERFI[window]

            wmin = argmin(DERFI[window])
            wmax = argmax(DERFI[window])

            if wmax<wmin
                wsig = dw0[range(wmax,wmin)]
                zc = zerocross(wsig)
                zc = map(x -> x.pos, zc)

                if length(zc) != 0
                    if abs(DERFI[wmax])>4*abs(DERFI[wmin])
                        # push!(Tou, zc[1]+wmax+bwind+R[i])
                        push!(T, zc[1]+wmax+bwind+R[i])
                    else
                        # push!(Tud, zc[1]+wmax+bwind+R[i])
                        push!(T, zc[1]+wmax+bwind+R[i])
                    end
                end
            else
                wsig = dw0[range(wmin,wmax)]
                zc = zerocross(wsig)
                zc = map(x -> x.pos, zc)

                if length(zc) != 0
                    mina = minimum(dw0[wmax:length(dw0)])
                    if abs(DERFI[wmax])<4*mina
                        # push!(Tud, zc[1]+wmin+bwind+R[i])
                        push!(T, zc[1]+wmin+bwind+R[i])
                    else
                        if abs(DERFI[wmin])>4*abs(DERFI[wmax])
                            # push!(Tod, zc[1]+wmin+bwind+R[i])
                            push!(T, zc[1]+wmin+bwind+R[i])
                        else
                            # push!(Tdu, zc[1]+wmin+bwind+R[i])
                            push!(T, zc[1]+wmin+bwind+R[i])
                        end
                    end
                end
            end
        end

        if length(T) < i push!(T, 0) end
        
    end

    T = map(x -> x > length(DERFI) ? length(DERFI) : x, T)
    # return Tud, Tdu, Tod, 
    return T
end

# Определение границ волн 
function bounds_coeff(DERFI, pk, dermax, wavetype)

    f = abs(DERFI[pk]*10/DERFI[dermax])

    kb = Float64[]
    ke = Float64[]

    if wavetype == "P"
        kb = 1.35
        ke = 2.0
    elseif wavetype == "Q"
        kb = 1.8
        ke = 1.0
    elseif wavetype == "S"
        kb = 1.0
        if f < 4.0
            ke = 3.0
        elseif f >= 4.0 && f < 4.75
            ke = 8.0
        elseif f >= 4.75 && f < 6.20
            ke = 9.0
        elseif f >= 6.20
            ke = 12.0
        end
    elseif wavetype == "T"
        kb = 2.0
        if f <= 0.13 
            ke = 4.0
        elseif f > 0.13 f < 0.20
            ke = 5.0
        elseif f >= 0.20 f < 0.41
            ke = 6.0
        elseif f >= 0.41
            ke = 7.0
        end
    end

    return kb, ke
end

function wavebounds(P, DERFI, R, wavetype, isinv)

    dermax = find_dermax(DERFI, R)

    # if wavetype == "P" || wavetype == "T"
    #     isinv = false;
    # elseif wavetype == "Q" || wavetype == "S"
    #     isinv = true;
    # else
    #     return "Error"
    # end

    Pwave = Wave[];
    for i in 1:lastindex(P)
        pkb = Float64[];
        pke = Float64[];

        if P[i] != 0
            # идём влево от метки
            for j in P[i]:-1:2
                if isinv
                    if DERFI[j-1] > DERFI[j]
                        pkb = j;
                        break
                    end
                else
                    if DERFI[j-1] < DERFI[j]
                        pkb = j;
                        break
                    end
                end
            end

            # идём вправо от метки
            for j in P[i]+1:lastindex(DERFI)
                if isinv
                    if DERFI[j] < DERFI[j-1]
                        pke = j-1;
                        break
                    end
                else
                    if DERFI[j] > DERFI[j-1]
                        pke = j-1;
                        break
                    end
                end
            end

            # рассчет порога
            if pkb != [] && pke != [] && dermax[i] != 0

                kb, ke = bounds_coeff(DERFI, pke, dermax[i], wavetype)

                THb = DERFI[pkb]/kb;
                THe = DERFI[pke]/ke;

                # идём влево от pkb
                wb = 0;
                for p in pkb:-1:2
                    if isinv
                        if DERFI[p-1] >= THb && DERFI[p] < THb
                            wb = p-1;
                            break
                        end
                    else
                        if DERFI[p-1] <= THb && DERFI[p] > THb
                            wb = p-1;
                            break
                        end
                    end
                end

                # идём вправо от pke
                we = 0;
                for p in pke:lastindex(DERFI)
                    if isinv
                        if DERFI[p] <= THe && DERFI[p-1] > THe
                            we = p;
                            break
                        end
                    else
                        if DERFI[p] >= THe && DERFI[p-1] < THe
                            we = p;
                            break
                        end
                    end
                end

                push!(Pwave, Wave(wb, we))
            end
        end

        if length(Pwave) < i push!(Pwave, Wave(0, 0)) end
    end

    return Pwave
end

function find_dermax(differed,R)
    dermax = Int[]
    k = 0

    for i in R
        k += 1
        for j in i:-1:2
            if differed[j-1] < differed[j]
                push!(dermax, j)
                break
            end
        end
        if length(dermax) < k push!(dermax, 0) end
    end

    return dermax
end

# Коррекция позиций
function pos_correct(delay_1, delay_2, delay_3, P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te)
    dtp = delay_1-delay_2+delay_3;

    T = T .- dtp;
    Tb = Tb .- dtp;
    Te = Te .- dtp;

    P = P .- dtp;
    Pb = Pb .- dtp;
    Pe = Pe .- dtp;

    dqrs = delay_1-delay_2;

    Q = Q .- dqrs;
    R = R .- dqrs;
    S = S .- dqrs;

    Qb = Qb .- dqrs;
    Se = Se .- dqrs;

    P = map((x) -> x < 0 ? 1 : x, P)
    Pb = map((x) -> x < 0 ? 1 : x, Pb)
    Pe = map((x) -> x < 0 ? 1 : x, Pe)
    Q = map((x) -> x < 0 ? 1 : x, Q)
    Qb = map((x) -> x < 0 ? 1 : x, Qb)
    R = map((x) -> x < 0 ? 1 : x, R)
    S = map((x) -> x < 0 ? 1 : x, S)
    Se = map((x) -> x < 0 ? 1 : x, Se)
    T = map((x) -> x < 0 ? 1 : x, T)
    Tb = map((x) -> x < 0 ? 1 : x, Tb)
    Te = map((x) -> x < 0 ? 1 : x, Te)

    return P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te
end

function plot_bounds(Pb,Pe,A)

    x = []
    y = []

    if Pe != []
        for i in 1:lastindex(Pb)
            push!(x, [Pb[i], Pb[i]])
            push!(x, [Pe[i], Pe[i]])
            push!(y, [-A, A])
            push!(y, [-A, A])
        end
    else
        for i in 1:lastindex(Pb)
            push!(x, [Pb[i], Pb[i]])
            push!(y, [-A, A])
        end
    end

    return x, y
end

# Поиск границ волн
function findwavebounds(differed, DERFI, filtered, P, Q, R, S, T)

    p = unique(P); p = filter( x -> x!=0, P)
    q = unique(Q); q = filter( x -> x!=0, Q)
    r = unique(R); r = filter( x -> x!=0, R)
    s = unique(S); s = filter( x -> x!=0, S)
    t = unique(T); t = filter( x -> x!=0, T)

    isinv_p = mean(p)==0 ? false : (mean(filtered[p]) < 0);
    isinv_t = mean(t)==0 ? false : (mean(filtered[t]) < 0);
    isinv_q = mean(q)==0 ? false : (mean(filtered[q]) < 0);
    isinv_s = mean(s)==0 ? false : (mean(filtered[s]) < 0);

    Pwave = wavebounds(P, DERFI, R, "P", isinv_p);
    Twave = wavebounds(T, DERFI, R, "T", isinv_t);

    Qwave = wavebounds(Q, differed, R, "Q", isinv_q);
    Swave = wavebounds(S, differed, R, "S", isinv_s);

    Pb = map((x) -> x.b, Pwave); Pe = map((x) -> x.e, Pwave)
    Tb = map((x) -> x.b, Twave); Te = map((x) -> x.e, Twave)
    Qb = map((x) -> x.b, Qwave); Se = map((x) -> x.e, Swave)

    return P, Pb, Pe, Q, Qb, R, S, Se, T, Tb, Te
end

## отбраковка максимум двух границ и поиск самых ранних onset и поздних end
## g = 6, 6, 6, 10, 12 для Pb, Pe, Qb, Se и Te соответственно
function multilead_bounds(bounds, g, type)
    if length(bounds) <= 2
        if type == "onset" bnd = minimum(bounds)
        elseif type == "end" bnd = maximum(bounds)
        end
    else
        bounds_sorted = sort(bounds)
        if type == "onset"
            min = bounds_sorted[1]
            rej = filter(x -> x < min + g, bounds_sorted)

            if !isempty(rej) && length(rej) <= 2
                ind = findall(x -> !(x in rej), bounds_sorted)
                bounds_sorted = bounds_sorted[ind]
            end

            bnd = minimum(bounds_sorted)

        elseif type == "end"
            max = bounds_sorted[end]
            rej = filter(x -> x > max - g, bounds_sorted)
            if !isempty(rej) && length(rej) <= 2
                ind = findall(x -> !(x in rej), bounds_sorted)
                bounds_sorted = bounds_sorted[ind]
            end

            bnd = maximum(bounds_sorted)
        end
    end

    return bnd
end
