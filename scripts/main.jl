include("../src/readfiles.jl")     # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/one_seg_calc.jl")

dir = "D:/INCART/Pulse_Data/все базы/Шумовая база"  # путь к базе
files = readdir(dir)
allbins = files[findall(x -> split(x, ".")[end] == "bin", files)]

function refADcalc(tonemkp::Vector{ToneEv}, ad::Vector{AD}, Tone::Vector{Float64}, Pres::Vector{Float64}, fs, seg, i::Int64)

    segm = Tone[seg[i].ibeg:seg[i].iend]
    pres = Pres[seg[i].ibeg:seg[i].iend]

    smoothtone = my_butter(segm, 2, 60, fs, "low") # сглаженный тонов
    ftone = my_butter(smoothtone, 2, 30, fs, "high") # фильтрованный тонов
    segm = my_butter(abs.(ftone), 2, 10, fs, "low") # огибающая по модулю

    peaks = map(x -> (x.pos - seg[i].ibeg + 1), tonemkp)
    # peaks = filter(x -> x!=0, peaks)

    ipresmax = argmax(pres)
    pumppeaks = filter(x -> x<ipresmax, peaks)
    descpeaks = filter(x -> x>=ipresmax, peaks)
    ipumpmax = argmax(segm[pumppeaks]); apumpmax = maximum(segm[pumppeaks])
    idescmax = argmax(segm[descpeaks]); adescmax = maximum(segm[descpeaks])

    refDADpump = round(pres[pumppeaks[1]]) |> Int64; refSADpump = round(pres[pumppeaks[end]]) |> Int64
    refSADdesc = round(pres[descpeaks[1]]) |> Int64; refDADdesc = round(pres[descpeaks[end]]) |> Int64

    for j in pumppeaks[ipumpmax:-1:1] if segm[j] <= 0.2*apumpmax refDADpump = round(pres[j]) |> Int64; break end end
    for j in pumppeaks[ipumpmax:1:end] if segm[j] <= 0.2*apumpmax refSADpump = round(pres[j]) |> Int64; break end end


    if !isempty(ad) && length(ad) >= i
        refSADdesc = ad[i].SAD
        refDADdesc = ad[i].DAD
    else
        for j in descpeaks[idescmax:-1:1] if segm[j] <= 0.2*adescmax refSADdesc = round(pres[j]) |> Int64; break end end
        for j in descpeaks[idescmax:1:end] if segm[j] <= 0.2*adescmax refDADdesc = round(pres[j]) |> Int64; break end end
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
    ad = AD[]
    try
        adtablefile = "D:/INCART/Pulse_Data/ad result tables/$basename/$(nm)_table_ad.txt"
        ad = read_ad(adtablefile)
    catch e
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

    filename = "alg markup/$basename/$nm"

    for i in 1:lastindex(ext)
        save_markup(allmkp[i], filename, ext[i])
    end

    adfilename = "D:/INCART/Pulse_Data/ref AD/$basename/$nm.ad"
    save_markup(adfilename, refAD)
end