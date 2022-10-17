
include("../src/my_filt.jl")

struct Peak
    min_pos::Int64
    max_pos::Int64
end

struct PDtk
    minmax::Peak
    updown::Peak
end

struct Seg
    bg::Int64
    en::Int64
end

struct Params 
    ipres::Int64
    imx::Int64
    imn0::Int64
    imn::Int64
    Wprev::Float64
    Wnext::Float64
    P1::Float64
    P2::Float64
    Amx::Float64
    Amn::Float64
    Range2::Float64
    speed::Float64
    delta::Float64
    Range::Float64
    speedPrev::Float64
    dSpeed::Float64
    pres::Float64
end

# фильтры + детектор + параметризатор для АД
function calc_ad(seg, Fs)

    # Фильтрация
    fsig_smooth = my_butter(seg, 2, 10, Fs, "low") # сглаживание
    fsig = my_butter(fsig_smooth, 2, 0.3, Fs, "high") # устранение постоянной составляющей

    # детектор
    pk = pkAD(fsig, fs)

    # параметризатор
    events = paramAD(pk, seg, fsig, fs)

    # отбраковка
    bad = discardAD(events, fs)

    return events, bad, fsig_smooth, fsig
end

# тахо для АД
function diff_filter(sig, fs)

    # к-ты
    b = fill(0, trunc(0.040*fs) |> Int)
    b[1] = 1
    b[end]= -1

    a = 1

    # фильтрация
    filtered = DSP.filt(b, a, sig);

    return filtered
end

# пиковый детектор с разрядами для АД
function pkAD(filtered, fs)

    # тахо
    tacho = diff_filter(filtered, fs);

    minDist = trunc(0.250 * fs) |> Int
    minR = 0.1

    N = length(tacho)

    LvP = -Inf;       # уровень предыдущего предполагаемого максимума тахо
    LvR = tacho[1];   # уровень слежения-разряда
    LvZ = 1           # средний уровень пика
    kR = 1/400*(1000/fs);
    kZ = 0.2;

    up = 1;              # предполгаемые детекции вниз и вверх
    mxCnt = 1            # счетчик отсчетов после максимума

    pkUp = 1             # подтверждённый максимум сигнала
    pkDn = 1             # подтверждённый минимум сигнала

    lookDn = false        # false - была детекия вниз

    lookMn = false
    mnVal = filtered[1]  # значение предполагаемого минимума
    mnValPrev = filtered[1]
    mnCnt = 1            # счетчик отсчетов после минимума
    mnCntPrev = 1
    mnPrevPos = 1        # позиция предыдущего минимума

    modeZero = true

    detections = PDtk[]

    # trend = fill(0.0, N)
    # level = fill(0.0, N)

    for i in 1:N

        # trend[i] = LvZ
        # level[i] = LvR

        # по фильтрованному сигналу
        if mnVal >= filtered[i] mnVal = filtered[i]; mnCnt = 0; end # предполагаемый минимум
        if mnValPrev >= filtered[i] mnValPrev = filtered[i]; mnCntPrev = 0; end 

        # по тахо
        if modeZero
            if tacho[i] > LvZ*0.1   # детекция вверх
                up = i
                modeZero = false
            end
        else
            if tacho[i] < LvZ*0.1   # детекция вниз
                if lookDn pkDn = i; lookDn = false end
                modeZero = true
            end
        end

        if tacho[i] > LvR # поиск максимума
            LvR = tacho[i]
            if LvP < LvR
                LvP = LvR
                pkUp = up   # подтверждаем предыдущий пик
                lookDn = true

                mxCnt = 0
                mnCnt = 0
                mnVal = filtered[i]
            end
            if lookMn  
                mnPrevPos = i - mnCntPrev
                lookMn = false
                mnCntPrev = 0
                mnValPrev = filtered[i]
            end
        else # медленный разряд
            LvR = LvR*(1-kR)
            if LvR < minR LvR = minR end
        end

        if mxCnt == minDist # прошло достаточно времени после предыдущего максимума
            if LvP > LvZ && ~lookDn # если пик выше порога по амплитуде и была детекция вниз
                minmax = Peak(mnPrevPos, i - mxCnt) # пик и детекция вверх тахо
                updown = Peak(pkDn, pkUp)           # пики фильтрованного сигнала

                push!(detections, PDtk(minmax, updown))

                mnCntPrev = mnCnt
                mnValPrev = mnVal

                lookMn = true
            end

            err = LvP/2 - LvZ;
            if err > LvZ err = LvZ; end # мин-макс динамика
            if err < -LvZ err = -LvZ; end
            if err > 0 dZ = err*kZ # разная динамика вверх-вниз
            else dZ = err*kZ*4 end

            LvZ = LvZ + dZ
            LvP = 0 #сброс
        end
        
        mxCnt += 1
        mnCnt += 1
        mnCntPrev += 1
    end

    return detections
end

# ------------------------------------------------------------
# Параметризация для АД
function paramAD(pk, raw, filtered, fs)

    L = length(pk)

    events = Params[]

    for i in 2:L-1
        iup_prev = pk[i-1].updown.max_pos
        idown_prev = pk[i-1].updown.min_pos

        iup = pk[i].updown.max_pos
        idown = pk[i].updown.min_pos

        iup_next = pk[i+1].updown.max_pos
        idown_next = pk[i+1].updown.min_pos

        ipres = iup;   # точка привязки давления
        imx = idown;   # макс на исходном
        imn0 = pk[i].minmax.min_pos  # минимум перед
        imn = pk[i+1].minmax.min_pos # минимум после

        Wprev = idown - idown_prev # расстояние от предыдущего (~RR)
        Wnext = idown_next - idown # расстояние до следующего (~RR)
        
        P1 = filtered[iup] # нужно компенсировать задержку до??
        P2 = filtered[idown]

        # размах по мин-максу на фильтрованном
        Amx = filtered[imx]
        Amn = filtered[imn]
        Range2  = Amx - Amn

        # треугольник на давлении
        dpdt = (raw[iup_next] - raw[iup]) / (iup_next - iup)
        speed = dpdt * fs
        
        delta = dpdt * (idown - iup)
        Range = (raw[idown] - raw[iup]) - dpdt*(idown - iup)
        
        dpdt = (raw[iup] - raw[iup_prev]) / (iup - iup_prev)
        speedPrev = dpdt * fs
        
        dSpeed = abs(speed - speedPrev)
        # Range - текущий, Range2 - для сравнения
        
        pres = raw[iup];

        push!(events, Params(ipres, imx, imn0, imn, Wprev, Wnext, P1, P2, Amx, Amn, 
                            Range2, speed, delta, Range, speedPrev, dSpeed, pres))
    end

    return events
end

# Отбраковка для АД
function discardAD(evt, fs)

    N = length(evt)
    bad = fill(0, N)
    bs = Tuple[]

    # время фронта:
    tfront = map(x -> x.imx - x.ipres, evt)

    # скважность
    duty = map((x,y) -> x/y.Wnext, tfront, evt)

    # отношение интервалов до/после:
    kW = map(x -> x.Wprev/x.Wnext, evt)

    # браковка:
    for i in 1:N
        s1 = evt[i].Wprev < 0.25*fs || 1.5*fs < evt[i].Wprev # ширина ДО + рефрактерность (хотя рефрактерность была убрана ранее)
        s2 = tfront[i] < 0.080*fs; # || 0.500*Fs < evt.n1 # ширина фронта % ??? искажения на макс ??? 
        # s3 = evt[i].Range < 1 || 50 < evt[i].Range # абсолютная амплитуда (размах)
        s3 = false
        s4 = 0.8 < duty[i] || duty[i] < 0.10 # ! .10 коэф. скважности ??? искажения на макс ???
        s5 = kW[i] < 0.5 && (evt[i].Wprev < 0.5*fs) # ??? преждевременный
        s6 = kW[i] > 2 && (evt[i].Wnext > 1*fs) # ??? широкий
        s7 = evt[i].speed > 100 || evt[i].speed < -100 # скорость накачки в пульсации
        s8 = evt[i].dSpeed > 20 # ??? широкий

        push!(bs, (s1, s2, s3, s4, s5, s6, s7, s8))
    end

    for j in 1:lastindex(bs)
        for i in lastindex(bs[j]):-1:1
            if bs[j][i]
                bad[j] = bad[j] + 2^(i-1)
            end
        end
    end

    # begs = map((x,y) -> if y != 0 x.imn0 end, evt, bad)
    # begs = begs[begs.!=nothing]

    # ends = map((x,y) -> if y != 0 x.imn end, evt, bad)
    # ends = ends[ends.!=nothing]

    # seg = map((x,y) -> Seg(x,y), begs, ends)

    # new_events = map((x,y,z) -> (tfront = x, duty = y, kW = z), tfront, duty, kW)

    return bad
end
