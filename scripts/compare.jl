# используется скрипт для сравнения от Юлии Александровны - переписать свой по тому же алгоритму
include("../src/beat2beat.jl") 
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 
include("../src/readfiles.jl") 

# nm = "PX11321102817293"
FILEN = 5

dir = "D:/INCART/Pulse_Data/bin"
files = readdir(dir)
allbins = files[findall(x -> split(x, ".")[end] == "bin", files)]

nm = split(allbins[FILEN],".")[1]

binfile = "D:/INCART/Pulse_Data/bin/"*nm*".bin"

ref_pres_filename = "C:/ktiflg/export/"*nm*".pls"
alg_pres_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/"*nm*".pres"

ref_tone_filename = "C:/ktiflg/export/"*nm*".ton"
alg_tone_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/"*nm*".tone"

# парсинг разметки
pres_ref = pls_ton_parse(ref_pres_filename);
pres_alg = pls_ton_parse(alg_pres_filename);

tone_ref = pls_ton_parse(ref_tone_filename);
tone_alg = pls_ton_parse(alg_tone_filename);

# по давлению сравниваем позиции минимумов
presTP = Vector{Int}[]
presFP = Vector{Int}[]
presFN = Vector{Int}[]
for i in 1:minimum([lastindex(pres_ref),lastindex(pres_alg)])
    fcr = map(x -> x.bg, pres_ref[i])
    fc = map(x -> x.bg, pres_alg[i])

    outcomp = calc_indexpairs(fcr, fc, 0.5*fs)

    push!(presTP, findall(x -> x[1] != -1 && x[2] != -1, outcomp))
    push!(presFP, findall(x -> x[1] == -1 && x[2] != -1, outcomp))
    push!(presFN, findall(x -> x[1] != -1 && x[2] == -1, outcomp))
end

# по тонам сравниваем позиции пиков (певый столбец в обеих разметках)
toneTP = Vector{Int}[]
toneFP = Vector{Int}[]
toneFN = Vector{Int}[]
for i in 1:minimum([lastindex(tone_ref),lastindex(tone_alg)])

    # только тоны в реф. разметке пишутся и на накачке, и на спуске - остальное и в реф.,
    # и в тестовой пишется и сравнивается только на спуске. чтобы не тянуть значения точек
    # начала спуска для каждого сегмента в файл разметки, тут костыль: первое значение из файла
    # разметки считаем точкой начала спуска, и точки из реф. берём только правее неё.
    fc = map(x -> x.bg, tone_alg[i])
    p = fc[1]
    down_fcr = map(x -> if x.bg>p x.bg end, tone_ref[i])
    fcr = down_fcr[down_fcr.!=nothing]

    outcomp = calc_indexpairs(fcr, fc, 0.5*fs)

    push!(toneTP, findall(x -> x[1] != -1 && x[2] != -1, outcomp))
    push!(toneFP, findall(x -> x[1] == -1 && x[2] != -1, outcomp))
    push!(toneFN, findall(x -> x[1] != -1 && x[2] == -1, outcomp))
end

# presTP
# presFP
# toneTP
# toneFP

function calc_tfpn(tfpn)
    n = 0
    for i in tfpn
        n += length(i)
    end

    return n
end

NpresTP = calc_tfpn(presTP)
NpresFP = calc_tfpn(presFP)
NpresFN = calc_tfpn(presFN)

SensPress = NpresTP/(NpresTP+NpresFN)

NtoneTP = calc_tfpn(toneTP)
NtoneFP = calc_tfpn(toneFP)
NtoneFN = calc_tfpn(toneFN) # КТResult3 выдаёт странную разметку по тонам, поэтому пока тут много FN

SensTone = NtoneTP/(NtoneTP+NtoneFN)