include("../src/my_filt.jl")

struct Edges
    ins1::Int64
    ins2::Int64
    ins3::Int64
    ins4::Int64
    ibeg::Int64
    iend::Int64
    W::Int64
    noise::Float64
    snr::Float64
end

# фильтры + детектор + параметризатор для Тонов
function calc_tone(seg_tone, smoothpres, fs)
    # фильтры
    smoothtone = my_butter(seg_tone, 2, 60, fs, "low") # сглаженный тонов
    ftone = my_butter(smoothtone, 2, 30, fs, "high") # фильтрованный тонов
    fstone = my_butter(abs.(ftone), 2, 10, fs, "low") # огибающая по модулю
    # детектор
    pos = pk_tone(fstone, fs)
    # параметризатор
    edg = find_edges(fstone, pos, fs)
    bad = discard_tone(smoothpres, fstone, pos, edg)

    return pos, bad, fstone
end

# детектор пиков для тонов
function pk_tone(sig, Fs)

    minDist = 0.200*Fs

    LvP = -Inf
    mxCnt = 0
    N = length(sig)
    extr = fill(false, N)

    minR = 90
    LvR = sig[1]  # уровень слежения-разряда
    kR = 1/20/4   # kR = 1/20;

    LvN = 100
    kN = 1/100/4 # kN = 1/100;
    kaN = 2.5

    for i in 1:N
        err = sig[i] - LvN
        LvN += kN*err

        if sig[i] > LvR # поиск максимума
            LvR = sig[i]
            if LvP < LvR
                LvP = LvR
                mxCnt = 0
            end
        else # медленный разряд
            LvR -= LvR*kR
            if LvR < minR LvR = minR end
        end

        if mxCnt == minDist
            if LvP > LvN * kaN # если пик выше по амплитуде, чем шум, и была детекция вниз
                extr[i-mxCnt] = true
            end
            
            LvP = 0; # x(i); # сброс
        end

        mxCnt += 1
    end

    pos = findall(extr)

    return pos
end

# параметризатор для тонов
function find_edges(x, pos, Fs)
    # поиск границ пиков по уровню 0.3
    
    N = length(x)
    Np = length(pos)
    
    Wmax = floor.(0.1*Fs) # макс полуширина, не более
    noiseLen = floor.(0.080*Fs) # длина участка поиска уровня шума
    noiseOffset1 = floor.(0.060*Fs) # fix(0.020*Fs); # отступ от границы ДО
    noiseOffset2 = floor.(0.060*Fs) # fix(0.030*Fs); # отступ от границы ПОСЛЕ

    Acoef = 1/3
    val = x[pos]

    edg = Edges[]
    
    for i=1:Np
        
        Lvl = floor.(val[i]*Acoef)
    
        # поиск ширины
        wBefore = 0
        ibeg = pos[i] |> Int64
        iend = maximum([ibeg-Wmax, 1]) |> Int64
        for k in ibeg:-1:iend
            if x[k] < Lvl break end
            wBefore += 1
        end
        i1 = ibeg-wBefore
        
        wAfter = 0
        ibeg = pos[i] |> Int64
        iend = minimum([pos[i]+Wmax, N]) |> Int64
        for k in ibeg : iend
            if x[k] < Lvl break end
            wAfter += 1
        end
        i2 = ibeg+wAfter
        
        W = wBefore + wAfter
        
        # поиск уровня шума
        ibeg = maximum([1, pos[i]-wBefore-1-noiseOffset1]) |> Int64
        iend = maximum([1, ibeg-noiseLen]) |> Int64
        nsBefore = maximum(x[iend:ibeg])
        ins1 = iend
        ins2 = ibeg
        
        ibeg = minimum([N, pos[i]+wAfter+1+noiseOffset2]) |> Int64
        iend = minimum([N, ibeg+noiseLen]) |> Int64
        nsAfter = maximum(x[ibeg:iend])
        ins3 = ibeg
        ins4 = iend
        
        noise = maximum([nsBefore, nsAfter])
        snr = x[pos[i]]/noise

        push!(edg, Edges(ins1, ins2, ins3, ins4, i1, i2,
                        W, noise, snr))
    end
        
    return edg
end

# браковка для тонов
function discard_tone(smoothpres, fstone, pos, edg)

    bad = fill(0, length(pos))
    pres = smoothpres[pos]

    badset = Tuple[]

    for i in 1:lastindex(edg)
        s1 = pres[i] < 30  # давление в точке меньше 30 мм
        s2 = edg[i].snr < 3.0 # мин сигнал/шум
        s3 = fstone[pos[i]] < 30 # амплитуда тона (0.3 на исходный/1000)
        push!(badset, (s1, s2, s3))
    end

    # браковка
    for i in 1:lastindex(bad)
        for j in 1:lastindex(badset[1])
            if badset[i][j]
                bad[i] = bad[i]+2^j
            end
        end
    end

    return bad
end
