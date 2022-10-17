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