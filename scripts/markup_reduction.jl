using DataFrames

# разметка первого тестового алгоритма в гуишную
include("../src/help_func.jl")
include("../src/readfiles.jl") 
include("../src/refmkpguifunctions.jl") 

# пути
rawpath = "D:/INCART/Pulse_Data/все базы"
mkppath = "alg markup"
adpath = "D:/INCART/Pulse_Data/ref AD"
basenm = "Шумовая база"
# filename = "PX11321102817293"

redalgmarkuppath = "formatted alg markup"

allbasefiles = readdir("$rawpath/$basenm")
allbins = filter(x -> split(x, ".")[end]=="bin", allbasefiles)
allfnames = map(x -> split(x, ".")[1], allbins)

for filename in allfnames
    # чтение исходного сигнала
    signals, fs, _, _ = readbin("$rawpath/$basenm/$filename")
    Tone = signals.Tone  # пульсации
    Pres = signals.Pres  # давление

    # получение валидных сегментов (тем же способом, что и в алгоритме - в противном случае вытянуть из самой разметки)
    vseg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

    Pres_mkp = test_markup_parse("$mkppath/$basenm/$filename.pres")
    Tone_mkp = test_markup_parse("$mkppath/$basenm/$filename.tone")

    # перетягиваем разметку из алгоритма в гуишную с поправкой на позицию начала сегмента, так как в гуишной
    # начало каждого сегмента принимается за ноль (но хранится информация о положении начала сегмента во всем сигнале)
    Pres_guimkp = map((x,d) -> map(y -> PresGuiMkp(y.ibeg-d.ibeg, y.iend-d.ibeg, 1), x), Pres_mkp, vseg)
    Tone_guimkp = map((x,d) -> map(y -> ToneGuiMkp(y.pos-d.ibeg, 1), x), Tone_mkp, vseg)

    # получение референтных границ АД и рабочей зоны на накачке и на спуске (ПО АМПЛИТУДЕ)
    ad_alg0 = read_alg_ad("$adpath/$basenm/$filename.ad")

    # границы рабочей зоны (пик треугольника давления - 10 мм, минимум спуска + 10 мм)
    # на спуске (ПО АМПЛИТУДЕ)
    wzbounds = map(x -> get_workzone_bounds(Pres[x.ibeg:x.iend]), vseg)

    # запись разметки в структуру файлов
    try readdir(redalgmarkuppath)
    catch e mkdir(redalgmarkuppath) end
    try readdir("$redalgmarkuppath/$basenm")
    catch e mkdir("$redalgmarkuppath/$basenm") end
    try readdir("$redalgmarkuppath/$basenm/$filename")
    catch e mkdir("$redalgmarkuppath/$basenm/$filename") end

    for i in 1:lastindex(vseg)

        ad_alg = (pump = Bounds(ad_alg0[i].pump.SAD, ad_alg0[i].pump.DAD), desc = Bounds(ad_alg0[i].desc.SAD, ad_alg0[i].desc.DAD))
        x = vseg[i]
        # получение позиций границ рабочей зоны и АД
        pump_adbounds = get_bounds(Pres[x.ibeg:x.iend], ad_alg.pump, true)
        desc_adbounds = get_bounds(Pres[x.ibeg:x.iend], ad_alg.desc, false)
        iadbounds = (pump = pump_adbounds, desc = desc_adbounds)

        pump_wzbounds = get_bounds(Pres[x.ibeg:x.iend], wzbounds[i].pump, true)
        desc_wzbounds = get_bounds(Pres[x.ibeg:x.iend], wzbounds[i].desc, false)
        iwzbounds = (pump = pump_wzbounds, desc = desc_wzbounds)

        fold = "$redalgmarkuppath/$basenm/$filename/$i"
        try readdir(fold)
        catch e mkdir(fold) end    
        SaveRefMarkup("$fold/pres.csv", Pres_guimkp[i])
        SaveRefMarkup("$fold/tone.csv", Tone_guimkp[i])
        SaveRefMarkup("$fold/bounds.csv", Pres[x.ibeg:x.iend], x, iadbounds, iwzbounds)
    end
end