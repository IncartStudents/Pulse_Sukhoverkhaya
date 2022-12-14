include("../src/readfiles.jl")     # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/one_seg_calc.jl")

dir = "D:/INCART/Pulse_Data/все базы/КТ 07 АД ЭКГ"  # путь к базе
files = readdir(dir)
allbins = files[findall(x -> split(x, ".")[end] == "bin", files)]

function refADcalc(tonemkp::Vector{ToneEv}, ad::Vector{AD}, Tone::Vector{Float64}, Pres::Vector{Float64}, fs, seg, i::Int64)

    segm = Tone[seg[i].ibeg:seg[i].iend]
    pres = Pres[seg[i].ibeg:seg[i].iend]

    smoothtone = my_butter(segm, 2, 60, fs, Lowpass) # сглаженный тонов
    ftone = my_butter(smoothtone, 2, 30, fs, Highpass) # фильтрованный тонов
    segm = my_butter(abs.(ftone), 2, 10, fs, Lowpass) # огибающая по модулю

    peaks = map(x -> (x.pos - seg[i].ibeg + 1), tonemkp)
    # peaks = map(x -> x.pos, tonemkp) # исп для теста с разметкой гуишного типа

    ipresmax = argmax(pres)
    pumppeaks = filter(x -> x<ipresmax, peaks)
    descpeaks = filter(x -> x>=ipresmax, peaks)
    ipumpmax = argmax(segm[pumppeaks]); apumpmax = maximum(segm[pumppeaks])
    idescmax = argmax(segm[descpeaks]); adescmax = maximum(segm[descpeaks])

    refDADpump = round(Int, pres[pumppeaks[1]]); refSADpump = round(Int, pres[pumppeaks[end]])
    refSADdesc = round(Int, pres[descpeaks[1]]); refDADdesc = round(Int, pres[descpeaks[end]])

    lp = pumppeaks[ipumpmax]
    for j in pumppeaks[ipumpmax:-1:1] 
        if segm[j] <= 0.2*apumpmax 
            refDADpump = abs(j-lp) < 2*fs ? round(Int, pres[j]) : round(Int, pres[lp])
            break 
        end 
        lp = j
    end

    lp = pumppeaks[ipumpmax]
    for j in pumppeaks[ipumpmax:1:end] 
        if segm[j] <= 0.2*apumpmax 
            refSADpump = abs(j-lp) < 2*fs ? round(Int, pres[j]) : round(Int, pres[lp])
            break 
        end 
        lp = j
    end


    if !isempty(ad) && length(ad) >= i
        refSADdesc = ad[i].SAD
        refDADdesc = ad[i].DAD
    else
        lp = descpeaks[idescmax]
        for j in descpeaks[idescmax:-1:1] 
            if segm[j] <= 0.2*adescmax 
                refSADdesc = abs(j-lp) < 2*fs ? round(Int, pres[j]) : round(Int, pres[lp])
                break 
            end 
            lp = j
        end

        lp = descpeaks[idescmax]
        for j in descpeaks[idescmax:1:end] 
            if segm[j] <= 0.2*adescmax 
                refDADdesc = abs(j-lp) < 2*fs ? round(Int, pres[j]) : round(Int, pres[lp])
                break 
            end 
            lp = j
        end
    end

    return (pump = AD(refSADpump, refDADpump), desc = AD(refSADdesc, refDADdesc))
end

for j in 1:length(allbins)

    FILEN = j

    nm = split(allbins[FILEN],".")[1]

    basename = split(dir, "/")[end]
    binfile = "$dir/$nm.bin" # обрабатываемый бинарь (там же рядом должен лежать hdr!)

    signals, fs, timestart, units = readbin(binfile);

    Pres = signals.Pres; # давление
    Tone = signals.Tone; # пульсации

    seg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

    # если есть таблица с реф. давлением на спуске
    ad = try adtablefile = "D:/INCART/Pulse_Data/ad result tables/$basename/$(nm)_table_ad.txt"
            ad = read_ad(adtablefile)
        catch e
            AD[]
        end

    # ######################################################

    Ln = size(seg)[1]

    tonemkp = fill(ToneEv[], Ln)
    presmkp = fill(PresEv[], Ln)
    refAD = fill((pump = AD(0,0), desc = AD(0,0)), Ln)

    for i in 1:Ln
        final_pres, final_tone = one_seg_calc(Pres, Tone, fs, seg, i)

        tonemkp[i] = final_tone
        presmkp[i] = final_pres

        refAD[i] = refADcalc(final_tone, ad, Tone, Pres, fs, seg, i)
    end

    # Сохранение разметки
    ext = [".tone", ".pres"]
    allmkp = [tonemkp, presmkp];

    try readdir("alg markup") catch e mkdir("alg markup") end
    try readdir("alg markup/$basename") catch e mkdir("alg markup/$basename") end

    filename = "alg markup/$basename/$nm"

    for i in 1:lastindex(ext)
        save_markup(allmkp[i], filename, ext[i])
    end

    adfilename = "D:/INCART/Pulse_Data/ref AD/$basename/$nm.ad"
    save_markup(adfilename, refAD)
end

# basenm = split(dir, "/")[end]
# nm = "MB1217211018180349"
# binfile = "D:/INCART/Pulse_Data/все базы/$basenm/$nm"

# signals, fs, timestart, units = readbin(binfile);

# Pres = signals.Pres; # давление
# Tone = signals.Tone; # пульсации

# # если есть таблица с реф. давлением на спуске
# ad = try
#     adtablefile = "D:/INCART/Pulse_Data/ad result tables/$basenm/$(nm)_table_ad.txt"
#     ad = read_ad(adtablefile)
# catch e
#     AD[]
# end

# mkp = ReadRefMkp("formatted alg markup/$basenm/$nm/1/tone.csv")

# bnd = ReadRefMkp("formatted alg markup/$basenm/$nm/1/bounds.csv")

# ad = AD[]
# bounds = Bounds[]
# push!(bounds, bnd.segm)
# mes = refADcalc(mkp, ad, Tone, Pres, fs, bounds, 1)

