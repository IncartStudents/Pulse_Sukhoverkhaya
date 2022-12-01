using Plots

include("../src/readfiles.jl")   # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями

# nm = "PX113211018174339"
nm = "PX113211028173836"

binfile = "D:/INCART/Pulse_Data/bin/"*nm*".bin"
pls = "C:/ktiflg/export/"*nm*".pls"
ton = "C:/ktiflg/export/"*nm*".ton"

signals, fs, timestart, units = readbin(binfile);

ECG = signals[1]    # экг
Pres = signals.Pres # давление
Tone = signals.Tone # пульсации

# парсинг разметки
pres_ref = pls_ton_parse(pls);
tone_ref = pls_ton_parse(ton);

# получение валидных фрагментов из разметки
valid_pres_ref = markup_seg(pres_ref)
valid_tone_ref = markup_seg(tone_ref)

n = 2

segmP = Pres[valid_pres_ref[n][1]:valid_pres_ref[n][2]]
segmT = Tone[valid_tone_ref[n][1]:valid_tone_ref[n][2]]
segmECG = ECG[valid_tone_ref[n][1]:valid_tone_ref[n][2]]

presPb = map(x -> x.bg-valid_pres_ref[n][1]+1, pres_ref[n])
presPe = map(x -> x.en-valid_pres_ref[n][1]+1, pres_ref[n][1:end-1])

tonePb = map(x -> x.bg-valid_tone_ref[n][1]+1, tone_ref[n])
# tonePe = map(x -> x.en-valid_tone_ref[n][1]+1, tone_ref[n][1:end-1])

plot(segmP*1000, fmt = :png, legend = false)
title!("Pres raw")
# savefig("figures/pres_raw.png")

scatter!(presPb, segmP[presPb]*1000, markersize = 2)
scatter!(presPe, segmP[presPe]*1000, markersize = 2, fmt = :png, legend = false)
title!("Pres ref markup")
# savefig("figures/pres_ref_mkp.png")

plot(segmT, fmt = :png, legend = false)
title!("Tone raw")
# savefig("figures/tone_raw.png")

scatter!(tonePb, segmT[tonePb], markersize = 2, fmt = :png)
title!("Tone ref markup")
# savefig("figures/tone_ref_markup.png")
