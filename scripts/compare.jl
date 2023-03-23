# используется скрипт для сравнения от Юлии Александровны - переписать свой по тому же алгоритму
include("../src/beat2beat.jl") 
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 
include("../src/readfiles.jl") 

dir = "D:/INCART/Pulse_Data/all bases/KT 07 AD ECG"
files = readdir(dir)
allbins = files[findall(x -> split(x, ".")[end] == "bin", files)]

for i in 1:lastindex(allbins)
    FILEN = i

    nm = split(allbins[FILEN],".")[1]

    binfile = "D:/INCART/Pulse_Data/bin/$nm.bin"
    signals, fs, timestart, units = readbin(binfile);

    ref_pres_filename = "C:/ktiflg/export/$nm.pls"
    alg_pres_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/$nm.pres"

    ref_tone_filename = "C:/ktiflg/export/$nm.ton"
    alg_tone_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/$nm.tone"

    # парсинг разметки
    pres_ref = pls_ton_parse(ref_pres_filename);
    pres_alg = test_markup_parse(alg_pres_filename);

    tone_ref = pls_ton_parse(ref_tone_filename);
    tone_alg = test_markup_parse(alg_tone_filename);

    # ПО ДАВЛЕНИЮ сравниваем позиции максимумов
    pres_compare_result = compare_result(pres_ref, pres_alg, 0.15*fs, "pres")

    # ПО ТОНАМ сравниваем позиции пиков (певый столбец в обеих разметках
    tone_compare_result = compare_result(tone_ref, tone_alg, 0.15*fs, "tone")

    filename = "alg markup/"*nm
    ext = [".tonecomp",".prescomp"]
    compresults = [tone_compare_result, pres_compare_result]
    for i in 1:lastindex(ext)
        save_compare_results(compresults[i], filename, ext[i])
    end
end

