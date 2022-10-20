include("../src/readfiles.jl")     # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/one_seg_calc.jl")

dir = "D:/INCART/Pulse_Data/bin"
files = readdir(dir)
allbins = files[findall(x -> split(x, ".")[end] == "bin", files)]

for j in 1:length(allbins)

    FILEN = j

    nm = split(allbins[FILEN],".")[1]

    binfile = "D:/INCART/Pulse_Data/bin/"*nm*".bin" # обрабатываемый бинарь (там же рядом должен лежать hdr!)

    signals, fs, timestart, units = readbin(binfile);

    Pres = signals.Pres; # давление
    Tone = signals.Tone; # пульсации

    seg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

    # ######################################################

    Ln = size(seg)[1]

    tonemkp = fill(ToneEv[], Ln)
    presmkp = fill(PresEv[], Ln)

    for i in 1:Ln
        final_pres, final_tone = one_seg_calc(Pres, Tone, fs, seg, i)

        tonemkp[i] = final_tone
        presmkp[i] = final_pres
    end

    # Сохранение разметки
    ext = [".tone", ".pres"]
    allmkp = [tonemkp, presmkp];

    filename = "alg markup/"*nm

    for i in 1:lastindex(ext)
        save_markup(allmkp[i], filename, ext[i])
    end
end

