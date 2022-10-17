using Plots

include("../src/readfiles.jl")   # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/one_seg_calc.jl")

# nm = "PX11321102817293"

FILEN = 1

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
    mkp = Vector{Peak}[]
    tonemkp = Vector{Int}[]

    pres_paramed = Vector{Peak}[]
    tone_paramed = Vector{Int64}[]

    ftone = []
    fpres = []
    for i in 1:size(seg)[1]
        final_pres, paramed_pres, final_tone, paramed_tone, fstone, fad = one_seg_calc(Pres, Tone, fs, seg, i)
        push!(ftone, fstone)
        push!(mkp, final_pres)
        push!(tonemkp, final_tone)
        push!(pres_paramed, paramed_pres)
        push!(tone_paramed, paramed_tone)
        push!(fpres, fad)
    end

    # Сохранение разметки (после параметризации & только спуск)
    ext = [".pres", ".tone"]
    allmkp = [pres_paramed, tone_paramed];

    filename = "alg markup/"*nm

    for i in 1:lastindex(ext)
        save_markup(allmkp[i], filename, ext[i])
    end

    # графики 
    for k in 1:size(seg)[1]
        p1 = fpres[k]
        t1 = ftone[k]

        min = map(x -> x.min_pos-seg[k,1]+1, pres_paramed[k])
        max = map(x -> x.max_pos-seg[k,1]+1, pres_paramed[k])

        plot(p1, fmt = :png, legend = false)
        scatter!(min, p1[min], markersize = 2)
        scatter!(max, p1[max], markersize = 2)
        title!("Pres "*string(k))
        savefig("figures/"*nm*"_pres_"*string(k)*".png")

        pks = tone_paramed[k] .- seg[k,1] 
        plot(t1, fmt = :png, legend = false)
        scatter!(pks, t1[pks], markersize = 2)
        title!("Tone")
        savefig("figures/"*nm*"_tone_"*string(k)*".png")
    end
end

###### графики

# plot(p1, fmt = :png, legend = false)
# title!("Pres raw (from alg)")
# # savefig("figures/pres_raw_from_alg.png")

# min = map(x -> x.min_pos-seg[2,1]+1, mkp[2])
# max = map(x -> x.max_pos-seg[2,1]+1, mkp[2])

# scatter!(min, p1[min], markersize = 2)
# scatter!(max, p1[max], markersize = 2)
# title!("Alg markup")
# savefig("figures/pres_alg_mkp.png")

# plot(t1, fmt = :png, legend = false)
# title!("Tone raw (from alg)")
# # savefig("figures/tone_raw_from_alg.png")

# pks = tonemkp[2] .- seg[2,1] 

# scatter!(pks, t1[pks], markersize = 2)
# title!("Alg markup")
# # savefig("figures/tone_alg_mkp.png")

# # после параметризации
# min = map(x -> x.min_pos-seg[2,1]+1, pres_paramed[2])
# max = map(x -> x.max_pos-seg[2,1]+1, pres_paramed[2])

# plot(p1, fmt = :png, legend = false)
# scatter!(min, p1[min], markersize = 2)
# scatter!(max, p1[max], markersize = 2)
# title!("Pres")
# savefig("figures/"*nm*".png")

# pks = tone_paramed[2] .- seg[2,1] 
# plot(t1, fmt = :png, legend = false)
# scatter!(pks, t1[pks], markersize = 2)
# title!("Tone")
# savefig("figures/"*nm*".png")

# ######################################################

### скаттерограммы признаков
# wprev = map(x -> x.Wprev, events)
# scatter(wprev[1:end-1], wprev[2:end], legend = false)

# tfrnt = map(x -> x.tfront, new_events)
# scatter(tfrnt[1:end-1], tfrnt[2:end], legend = false)

# rng = map(x -> x.Range, events)
# scatter(rng[1:end-1], rng[2:end], legend = false)

# dty = map(x -> x.duty, new_events)
# scatter(dty[1:end-1], dty[2:end], legend = false)

# kw = map(x -> x.kW, new_events)
# scatter(kw[1:end-1], kw[2:end], legend = false)