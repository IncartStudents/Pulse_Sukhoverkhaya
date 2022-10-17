include("../src/my_filt.jl")     # ФВЧ и ФНЧ баттерворда с задаваемой частотой
include("../src/calc_AD.jl")     # для работы с АД
include("../src/calc_Tone.jl")   # для работы с тонами
include("../src/help_func.jl")   # вспомогательные функции

function one_seg_calc(Pres, Tone, fs, seg, n)

    # выделение сегмента
    pres = Pres[seg[n,1]:seg[n,2]]
    tone = Tone[seg[n,1]:seg[n,2]]

    # обработка АД
    events_pres, bad_ad, ad_smooth, fad = calc_ad(pres, fs)
    # обработка тонов
    pos_tone, bad_tone, fstone = calc_tone(tone, ad_smooth, fs)

    # момент начала стравливания воздуха
    presstr = maximum(pres)              # амплитуда вершины треугольника давления
    p = findall(x -> x == presstr, pres) # позиция вершины треугольника давления

    # 1.1. РАЗМЕТКА АД ДО ПАРАМЕТРИЗАЦИИ (ТОЛЬКО СПУСК)
    min = map(x -> x.imn, events_pres)
    max = map(x -> x.imx, events_pres)

    # выбор меток правее вершины треугольника давления (спуск)
    down_pres = map((x,y) -> x>y ? true : false, max, fill(p[1], length(max)))
    ind_pres = findall(down_pres)

    # запись разметки в структуру
    final_pres = map((x,y) -> Peak(x+seg[n,1], y+seg[n,1]), min[ind_pres], max[ind_pres])

    # 1.1. РАЗМЕТКА ТОНОВ ДО ПАРАМЕТРИЗАЦИИ (ТОЛЬКО СПУСК)
    pos_fin = pos_tone[findall(map((x,y) -> x>y, pos_tone, fill(p[1], length(pos_tone))))]
    final_tone = pos_fin .+ seg[n,1]

    # 2.1 РАЗМЕТКА АД ПОСЛЕ ПАРАМЕТРИЗАЦИИ
    notbad = findall(x -> x==0, bad_ad)
    ver = max[notbad]
    down = map((x,y) -> x>y ? true : false, ver, fill(p[1], length(ver)))
    ind = findall(down)

    paramed_pres = map((x,y) -> Peak(x+seg[n,1], y+seg[n,1]), min[notbad][ind], max[notbad][ind])

    # 2.1 РАЗМЕТКА ТОНОВ ПОСЛЕ ПАРАМЕТРИЗАЦИИ
    notbad = findall(x -> x==0, bad_tone)
    ver = pos_tone[notbad]
    nb = ver[findall(map((x,y) -> x>y, ver, fill(p[1], length(ver))))]
    paramed_tone = nb .+ seg[n, 1]

    return final_pres, paramed_pres, final_tone, paramed_tone, fstone, fad
end