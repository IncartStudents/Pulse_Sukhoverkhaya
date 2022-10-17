# используется скрипт для сравнения от Юлии Александровны - переписать свой по тому же алгоритму
include("../src/beat2beat.jl") 
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 

ref_pres_filename = "C:/ktiflg/export/PX11321102817293.pls"
alg_pres_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/PX11321102817293.pres"

ref_tone_filename = "C:/ktiflg/export/PX11321102817293.ton"
alg_tone_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/PX11321102817293.tone"

# парсинг разметки
pres_ref = pls_ton_parse(ref_pres_filename)
pres_alg = pls_ton_parse(alg_pres_filename)

tone_ref = pls_ton_parse(ref_tone_filename)
tone_alg = pls_ton_parse(alg_tone_filename)

# по давлению сравниваем позиции минимумов
presTP = Vector{Int}[]
presFP = Vector{Int}[]
for i in 1:lastindex(pres_ref)
    fcr = map(x -> x.bg, pres_ref[i])
    fc = map(x -> x.bg, pres_alg[i])

    outcomp = calc_indexpairs(fcr, fc, 0.5*fs)

    push!(presTP, findall(x -> x[1] != -1 && x[2] != -1, outcomp))
    push!(presFP, findall(x -> x[1] == -1 && x[2] != -1, outcomp))
end

# по тонам сравниваем позиции пиков (певый столбец в обеих разметках)
toneTP = Vector{Int}[]
toneFP = Vector{Int}[]
for i in 1:lastindex(tone_ref)
    fcr = map(x -> x.bg, tone_ref[i])
    fc = map(x -> x.bg, tone_alg[i])

    outcomp = calc_indexpairs(fcr, fc, 0.5*fs)

    push!(toneTP, findall(x -> x[1] != -1 && x[2] != -1, outcomp))
    push!(toneFP, findall(x -> x[1] == -1 && x[2] != -1, outcomp))
end

presTP
presFP
toneTP
toneFP
