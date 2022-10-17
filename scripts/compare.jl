# используется скрипт для сравнения от Юлии Александровны - переписать свой по тому же алгоритму
include("../src/beat2beat.jl") 
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 

ref_pres_filename = "C:/ktiflg/export/PX11321102817293.pls"
alg_pres_filename = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/PX11321102817293.pres"

# парсинг разметки
pres_ref = pls_ton_parse(ref_pres_filename)
pres_alg = pls_ton_parse(alg_pres_filename)

TP = Vector{Int}[]
FP = Vector{Int}[]

# по давлению сравниваем позиции минимумов
for i in 1:lastindex(pres_ref)
    fcr = map(x -> x.bg, pres_ref[i])
    fc = map(x -> x.bg, pres_alg[i])

    outcomp = calc_indexpairs(fcr, fc, 0.5*fs)

    push!(TP, findall(x -> x[1] != -1 && x[2] != -1, outcomp))
    push!(FP, findall(x -> x[1] == -1 && x[2] != -1, outcomp))
end

outcomp = calc_indexpairs(tonePb, pos, 0.5*fs)
TP = findall(x -> x[1] != -1 && x[2] != -1, outcomp)
FP = findall(x -> x[1] == -1 && x[2] != -1, outcomp)
