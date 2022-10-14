using DSP

# Функция поиска валидных пиков в ОДНОМ валидном сегменте АД
# фильтр + детектор + параметризатор

function CALC_AD(seg, Fs)
    # Фильтрация
    fsig_smooth = SmoothFilt(seg, Fs) # сглаживание
    fsig_smooth = fsig_smooth .|> round .|> Int64

    fsig = ConstRemove(fsig_smooth, Fs) # устранение постоянной составляющей
    fsig = fsig*1000 .|> round .|> Int64

    # детектор
    pk = pkAD(fsig, fs)

    # параметризатор
    events = paramAD(pk, seg, fsig, fs)
end

# Сглаживающий фильтр
function SmoothFilt(sig, Fs)

    responsetype = Lowpass(10; fs=Fs)
    designmethod = Butterworth(2)
    fsig = filt(digitalfilter(responsetype, designmethod), sig)

    return fsig
end

# Фильтр постоянной составляющей 
function ConstRemove(sig, Fs)

    responsetype = Highpass(0.3; fs=Fs)
    designmethod = Butterworth(2)
    fsig = filt(digitalfilter(responsetype, designmethod), sig)

    return fsig
end

# тахо
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

struct Peak
    min_pos::Int64
    max_pos::Int64
end

struct PDtk
    minmax::Peak
    updown::Peak
end

# пиковый детектор с разрядами
function pkAD(filtered_x_1000, fs)

    # тахо
    tacho = diff_filter(filtered_x_1000, fs);

    minDist = trunc(0.250 * fs) |> Int
    mm_lsb = 1000
    minR = 0.1 * mm_lsb

    N = length(tacho)

    LvP = -Inf;       # уровень предыдущего предполагаемого максимума тахо
    LvR = tacho[1];   # уровень слежения-разряда
    LvZ = 1 * mm_lsb; # средний уровень пика
    kR = 1/400*(1000/fs);
    kZ = 0.2;

    up = 1;              # предполгаемые детекции вниз и вверх
    mxCnt = 1            # счетчик отсчетов после максимума

    pkUp = 1             # подтверждённый максимум сигнала
    pkDn = 1             # подтверждённый минимум сигнала

    lookDn = false        # false - была детекия вниз

    lookMn = false
    mnVal = filtered_x_1000[1]  # значение предполагаемого минимума
    mnCnt = 1            # счетчик отсчетов после минимума
    mnPrevPos = 1        # позиция предыдущего минимума

    modeZero = true

    detections = PDtk[]

    # trend = fill(0.0, N)
    # level = fill(0.0, N)

    for i in 1:N

        # trend[i] = LvZ
        # level[i] = LvR

        # по фильтрованному сигналу
        if mnVal >= filtered_x_1000[i] mnVal = filtered_x_1000[i]; mnCnt = 0; end # предполагаемый минимум

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
                mnVal = filtered_x_1000[i]
            end
            if lookMn  
                mnPrevPos = i - mnCnt
                lookMn = false
                mnCnt = 0
                mnVal = filtered_x_1000[i]
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
    end

    return detections
end

# коррекция позиций пиков (компенсация задержки)

# тест пикового детектора 
# uppos = map(x -> x.updown.max_pos, detections)
# downpos = map(x -> x.updown.min_pos, detections)

# delay = 5 # компенсация задержки(???????????????)

# correct_up = map(x -> (x-delay) >= 1 ? x-delay : 1, uppos)
# correct_down = map(x -> (x-delay) >= 1 ? x-delay : 1, downpos)

# plot(filtered)

# scatter!(correct_up, filtered[correct_up], markersize = 2, label = "up")
# scatter!(correct_down, filtered[correct_down], markersize = 2, label = "down")

# xlims!(5000, 6000)
# ylims!(-500, 1000)

# ------------------------------------------------------------
struct Params 
    ipres::Float64
    imx::Float64
    imn0::Float64
    imn::Float64
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

struct badset
    s1::Bool
    s2::Bool
    s3::Bool
    s4::Bool
    s5::Bool
    s6::Bool
    s7::Bool
    s8::Bool
end

function discardAD(events, fs)

    seg = Vector[]
    kPress = 1e4

    N = length(events)
    bad = zeros(N, 1)
    bs = badset[]

    # время фронта:
    tfront = map(x -> x.imx - x.ipres, events)

    # скважность
    duty = map((x,y) -> x - y.Wnext, (tfront, events))

    # отношение интервалов до/после:
    kW = map(x -> x.Wprev/x.Wnext, events)

    # браковка:
    for i in 1:N
        s1 = evt[i].Wprev < 0.25*Fs || 1.5*Fs < evt[i].Wprev # ширина ДО + рефрактерность (хотя рефрактерность была убрана ранее)
        s2 = evt[i].n1 < 0.080*Fs; # || 0.500*Fs < evt.n1 # ширина фронта % ??? искажения на макс ??? 
        s3 = evt[i].Range < 0.3*kPres || 5*kPres < evt[i].Range # абсолютная амплитуда (размах)
        s4 = evt[i].duty < 0.12 || 0.8 < evt[i].duty # ! .10 коэф. скважности ??? искажения на макс ???
        s5 = evt[i].kW < 0.5 && (evt[i].Wprev < 0.5*Fs) # ??? преждевременный
        s6 = evt[i].kW > 2 && (evt[i].Wnext > 1*Fs) # ??? широкий
        s7 = evt[i].speed > 10*kPres || evt[i].speed < -10*kPres # скорость накачки в пульсации
        s8 = evt[i].dSpeed > 2*kPres # ??? широкий

        push!(bs, badset(s1, s2, s3, s4, s5, s6, s7, s8))
    end

    return bs
end

