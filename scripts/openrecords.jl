using Plots

include("../src/readfiles.jl")   # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями

binfile = "D:/INCART/Pulse_Data/bin/PX11321102817293.bin"
pls = "C:/ktiflg/export/PX11321102817293.pls"
ton = "C:/ktiflg/export/PX11321102817293.ton"

signals, fs, timestart, units = readbin(binfile);

ECG = signals.LR    # экг
Pres = signals.Pres # давление
Tone = signals.Tone # пульсации

# парсинг разметки
pres_ref = pls_ton_parse(pls)
tone_ref = pls_ton_parse(ton)

# получение валидных фрагментов из разметки
valid_pres_ref = markup_seg(pres_ref)
valid_tone_ref = markup_seg(tone_ref)

n = 2

segmP = Pres[valid_pres_ref[n][1]:valid_pres_ref[n][2]]
segmT = Tone[valid_tone_ref[n][1]:valid_tone_ref[n][2]]

presPb = map(x -> x.bg-valid_pres_ref[n][1]+1, pres_ref[n])
presPe = map(x -> x.en-valid_pres_ref[n][1]+1, pres_ref[n][1:end-1])

tonePb = map(x -> x.bg-valid_tone_ref[n][1]+1, tone_ref[n])
# tonePe = map(x -> x.en-valid_tone_ref[n][1]+1, tone_ref[n][1:end-1])

plot(segmP*1000)
scatter!(presPb, segmP[presPb]*1000, markersize = 2)
scatter!(presPe, segmP[presPe]*1000, markersize = 2, fmt = :png)
xlims!(4000, 8000)
ylims!(100000, 120000)

plot(segmT)
scatter!(tonePb, segmT[tonePb], markersize = 2, fmt = :png)
# scatter!(tonePe, segmT[tonePe], markersize = 2)
xlims!(4000, 6000)

savefig("tone_auto_markup.png")