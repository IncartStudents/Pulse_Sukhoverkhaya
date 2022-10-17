using DSP

function calc_tone(seg_tone, seg_pres, fs)
    # фильтры
    smoothpres = floor.(SmoothFilt(seg_pres, fs)) # сглаженный давления
    smoothtone = floor.(smooth_tone(seg_tone, fs)) # сглаженный тонов
    ftone = floor.(my_highpass(smoothtone, fs)) # фильтрованный тонов
    fstone = floor.(SmoothFilt(ftone, fs)) # огибающая по модулю
    # детектор
    pos = pk_tone(fstone, fs)
    # параметризатор
    edg = find_edges(fstone, pos, fs)
    bad = discard_tone(smoothpres, fstone, pos, edg)

    return pos, bad
end

# Сглаживающий фильтр
# Сглаживающий фильтр
function smooth_tone(sig, Fs)

    responsetype = Lowpass(60; fs=Fs)
    designmethod = Butterworth(2)
    fsig = filt(digitalfilter(responsetype, designmethod), sig)

    return fsig
end

function my_highpass(sig, Fs)

    responsetype = Highpass(30; fs=Fs)
    designmethod = Butterworth(2)
    fsig = filt(digitalfilter(responsetype, designmethod), sig)

    return fsig
end

function pk_tone(sig, Fs)

    minAmp = 0
    minDist = 0.200*Fs

    LvP = -Inf
    mxCnt = 0
    N = length(sig)
    extr = fill(false, N)

    minR = 90
    LvZ = 1000    # средний уровень пика
    LvR = sig[1]  # уровень слежения-разряда
    kR = 1/20/4   # kR = 1/20;

    LvN = 100
    kN = 1/100/4 # kN = 1/100;
    kaN = 2.5

    pk = []

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

struct Edges
    ins1::Int64
    ins2::Int64
    ins3::Int64
    ins4::Int64
    ibeg::Int64
    iend::Int64
    W::Int64
    noise::Int64
    snr::Float64
end

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
        nsAfter = maximum(x[ibeg:iend]) |> Int64
        ins3 = ibeg
        ins4 = iend
        
        noise = maximum([nsBefore, nsAfter]) |> Int64
        snr = x[pos[i]]/noise

        push!(edg, Edges(ins1, ins2, ins3, ins4, i1, i2,
                        W, noise, snr))
    end
        
    return edg
end

function discard_tone(smoothpres, fstone, pos, edg)

    kPres = 1e4
    kTone = 1e3

    bad = zeros(length(pos))
    pres = smoothpres[pos]

    badset = Tuple[]

    for i in 1:lastindex(edg)
        s1 = pres[i] < 3*kPres # давление в точке меньше 30 мм
        s2 = edg[i].snr < 3.0 # мин сигнал/шум
        s3 = fstone[pos[i]] < 0.3*kTone # амплитуда тона
        push!(badset, (s1, s2, s3))
    end

    # браковка
    for i in 1:lastindex(bad)
        for j in 1:lastindex(badset[1])
            if badset[i][j]
                bad[i] = j
            end
        end
    end

    return bad
end

function process_seg_tone(Tone, Pres, fs, seg, n, p)

    s1 = Tone[seg[n,1]:seg[n,2]]
    s2 = Pres[seg[n,1]:seg[n,2]]*1000

    pos0, bad = calc_tone(s1, s2, fs)
    pos = pos0[findall(map((x,y) -> x>y, pos0, fill(p[1], length(pos0))))]
    final = pos .+ seg[n,1]

    # псоле параметризации
    notbad = findall(x -> x==0, bad)
    ver = pos0[notbad]
    nb = ver[findall(map((x,y) -> x>y, ver, fill(p[1], length(ver))))]
    paramed = nb .+ seg[n, 1]

    return final, paramed
end