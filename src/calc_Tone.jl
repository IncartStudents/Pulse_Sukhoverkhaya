using DSP

function calc_tone(seg_tone, fs)
    # фильтры
    # smoothpress = floor.(SmoothFilt(seg_pres, fs)) # сглаженный давления
    smoothtone = floor.(smooth_tone(seg_tone, fs)) # сглаженный тонов
    ftone = floor.(my_highpass(smoothtone, fs)) # фильтрованный тонов
    fstone = floor.(SmoothFilt(ftone, fs)) # огибающая по модулю
    # детектор
    pos = pk_tone(fstone, fs)
    # параметризатор

    return pos
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
    
    W = zeros(size(pos))
    noise = W
    i1 = W
    i2 = W
    
    ins1 = W; ins2 = W; ins3 = W; ins4 = W;
    
    for i=1:Np
        
        Lvl = floor.(val(i)*Acoef);
    
        # поиск ширины
        wBefore = 0;
        ibeg = pos[i];
        iend = maximum(ibeg-Wmax, 1);
        for k in ibeg:-1:iend
            if x[k] < Lvl break end
            wBefore += 1;
        end
        i1[i] = ibeg-wBefore;
        
        wAfter = 0;
        ibeg = pos[i];
        iend = minimum(pos[i]+Wmax, N);
        for k in ibeg : iend
            if x[k] < Lvl break end
            wAfter += 1;
        end
        i2[i] = ibeg+wAfter;
        
        W[i] = wBefore + wAfter;
        
        # поиск уровня шума
        ibeg = maximum([1, pos[i]-wBefore-1-noiseOffset1]);
        iend = maximum([1, ibeg-noiseLen]);
        nsBefore = maximum(x[iend:ibeg]);
        ins1[i] = iend;
        ins2[i] = ibeg;
        
        ibeg = minimum([N, pos[i]+wAfter+1+noiseOffset2]);
        iend = minimum([N, ibeg+noiseLen]);
        nsAfter = maximum(x[ibeg:iend]);
        ins3[i] = ibeg;
        ins4[i] = iend;
        
        noise[i] = maximum([nsBefore, nsAfter]);
    end
        
    return pk
end

function process_seg_tone(Tone, fs, seg, n, p)

    s1 = Tone[seg[n,1]:seg[n,2]]

    pos0 = calc_tone(s1, fs)
    pos = findall(map((x,y) -> x>y, pos0, fill(p[1], length(pos0))))

    res = pos0[pos]

    return res
end