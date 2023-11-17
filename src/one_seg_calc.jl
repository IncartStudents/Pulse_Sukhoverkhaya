include("../src/my_filt.jl")     # ФВЧ и ФНЧ баттерворда с задаваемой частотой
include("../src/calc_AD.jl")     # для работы с АД
include("../src/calc_Tone.jl")   # для работы с тонами
include("../src/help_func.jl")   # вспомогательные функции

function one_seg_calc(Pres, Tone, fs, seg, n)

    # выделение сегмента
    pres = Pres[seg[n].ibeg:seg[n].iend]
    tone = Tone[seg[n].ibeg:seg[n].iend]

    # Обработка АД
    events_pres, bad_ad, ad_smooth = calc_ad(pres, fs)

    # РАЗМЕТКА АД
    # сборка в вектор структур
    final_pres = map((x,y) -> PresEv(x.imx.+seg[n].ibeg, x.imn.+seg[n].ibeg, y), events_pres, bad_ad)

    # Обработка тонов
    pos_tone, bad_tone = calc_tone(tone, ad_smooth, fs)

    # РАЗМЕТКА ТОНОВ
    # сборка в вектор структур
    final_tone = map((x,y) -> ToneEv(x, y), pos_tone.+seg[n].ibeg, bad_tone)

    return final_pres, final_tone
end