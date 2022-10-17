using Plots

include("../src/readfiles.jl")   # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 
include("../src/calc_AD.jl")     # добавления файла с функциями для работы с АД
include("../src/calc_tone.jl")   # добавления файла с функциями для работы с тонами

binfile = "D:/INCART/Pulse_Data/bin/PX11321102817293.bin" # обрабатываемый бинарь (там же рядом должен лежать hdr!)

signals, fs, timestart, units = readbin(binfile);

Pres = signals.Pres; # давление
Tone = signals.Tone; # пульсации

seg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

# ######################################################
mkp = Vector{Peak}[]
tonemkp = Vector{Int}[]
for i in 1:size(seg)[1]
    mkpi, p = process_seg(Pres, fs, seg, i)
    push!(mkp, mkpi)

    pos = process_seg_tone(Tone, fs, seg, i, p)
    push!(tonemkp, pos)
end

# Сохранение разметки (только спуск)
ext = [".pres", ".tone"]
allmkp = [mkp, tonemkp]

filename = "alg markup/"*split(split(binfile, "/")[end], ".")[end-1]

for i in 1:lastindex(ext)
    save_markup(allmkp[i], filename, ext[i])
end

###### графики
p1 = Pres[seg[2,1]:seg[2,2]]
t1 = Tone[seg[2,1]:seg[2,2]]

plot(p1, fmt = :png, legend = false)
title!("Pres raw (from alg)")
savefig("figures/pres_raw_from_alg.png")

min = map(x -> x.min_pos-seg[2,1]+1, mkp[2])
max = map(x -> x.max_pos-seg[2,1]+1, mkp[2])

scatter!(min, p1[min], markersize = 2)
scatter!(max, p1[max], markersize = 2)
title!("Alg markup")
savefig("figures/pres_alg_mkp.png")

plot(t1, fmt = :png, legend = false)
title!("Tone raw (from alg)")
savefig("figures/tone_raw_from_alg.png")

pks = tonemkp[2] .- seg[2,1] .+ 1

scatter!(pks, t1[pks], markersize = 2)
title!("Alg markup")
savefig("figures/tone_alg_mkp.png")
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