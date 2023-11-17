# Сравнение тестовой и референтной разметок гуишного формата
using DataFrames
using CSV

include("../src/beat2beat.jl") 
include("../src/refmkpguifunctions.jl")
include("../src/readfiles.jl") 

mutable struct Metrics
    filename::String # имя файла
    meas::Int64      # номер измерения

    TP::Int64
    FP::Int64
    TN::Int64
    FN::Int64
    Ac::Float64
    Se::Float64
    Sp::Float64
    PPV::Float64
    NPV::Float64
end

function getmetrics(compareres::Vector{Tuple{Int64, Int64}}, rbad::Vector{Int64}, tbad::Vector{Int64})

    TP = 0; TN = 0; FP = 0; FN = 0; Ac = 0.0; Se = 0.0; Sp = 0.0; PPV = 0.0; NPV = 0.0


    # на этом этапе причина браковки в тесте не интересует, поэтому метки типа из теста переводим в булевые
    tbad = map(x -> x == 0 ? false : true, tbad)

    for (ref, test) in compareres
        if ref != -1 && test != -1      # есть и в реф, и в тест (предполагаемый TP)
            if rbad[ref] != 2 && !tbad[test] TP += 1; # знач/незнач в реф и не збраковано в тесте
            elseif rbad[ref] == 2 && !tbad[test] FP += 1; # шум в реф и не забраковано в тесте
            elseif rbad[ref] != 1 && tbad[test] TN += 1; # не знач или шум в реф и забраковано в тесте
            elseif rbad[ref] == 1 && tbad[test] FN += 1; # знач в реф и забраковано в тесте
            end
        elseif ref == -1 && test == -1  # такого не должно случаться, но на всякий проверим
            println("оба -1!!!")
        elseif ref == -1 && test != -1  # нету в реф, есть в тест (предполагаемый FP)
            if !tbad[test] FP += 1 else TN += 1 end # если в тесте не забраковано, то FP, в противном случае - TN
        elseif ref != -1 && test == -1  # есть в реф, нету в тест (предполагаемый FN)
            if rbad[ref] == 1 FP += 1 else TN += 1 end
        end
    end

    Ac = (TP+TN)/(TP+TN+FP+FN); Ac = round(Ac*100, digits = 2)

    if (TP+FN) != 0 Se = TP/(TP+FN); Se = round(Se*100, digits = 2) end
    if (TN+FP) != 0 Sp = TN/(TN+FP); Sp = round(Sp*100, digits = 2) end

    if (TP+FP) != 0 PPV = TP/(TP+FP); PPV = round(PPV*100, digits = 2) end
    if (TN+FN) != 0 NPV = TN/(TN+FN); NPV = round(NPV*100, digits = 2) end

    return Metrics("", 0, TP, FP, TN, FN, Ac, Se, Sp, PPV, NPV)
end

function getmetrics(metrics::Vector{Metrics}) # среднее по всем метрикам

    sTP = 0; sTN = 0; sFP = 0; sFN = 0; Ac = 0.0; Se = 0.0; Sp = 0.0; PPV = 0.0; NPV = 0.0

    TP = map(x -> x.TP, metrics)
    TN = map(x -> x.TN, metrics)
    FP = map(x -> x.FP, metrics)
    FN = map(x -> x.FN, metrics)

    sTP = sum(TP)
    sTN = sum(TN)
    sFP = sum(FP)
    sFN = sum(FN)

    Ac = (sTP+sTN)/(sTP+sTN+sFP+sFN); Ac = round(Ac*100, digits = 2)

    if (sTP+sFN) != 0 Se = sTP/(sTP+sFN); Se = round(Se*100, digits = 2) end
    if (sTN+sFP) != 0 Sp = sTN/(sTN+sFP); Sp = round(Sp*100, digits = 2) end

    if (sTP+sFP) != 0 PPV = sTP/(sTP+sFP); PPV = round(PPV*100, digits = 2) end
    if (sTN+sFN) != 0 NPV = sTN/(sTN+sFN); NPV = round(NPV*100, digits = 2) end

    return Metrics("", 0, sTP, sFP, sTN, sFN, Ac, Se, Sp, PPV, NPV)
end

function isinbounds(pos::Int64, bounds::Bounds) # находится ли точка внутри границ (случай двух границ)
    ibeg = minimum([bounds.ibeg, bounds.iend])
    iend = maximum([bounds.ibeg, bounds.iend])
    return pos >= ibeg && pos <= iend
end

function isinbounds(pos::Int64, bounds::Vector{Bounds}) # находится ли точка хотя бы в одном из наборов по 2 границы
    for bnd in bounds
        if isinbounds(pos, bnd) return true end
    end
    return false
end

basename = "polifunctional"                       # название базы
testdir = "formatted alg markup/$basename"      # путь к тестовой разметке
refdir = "ref markup/$basename"                 # путь к референтной разметке
bindir = "D:/INCART/Pulse_Data/all bases/$basename" # путь к бинарям

# сравниваем тестовую с реферетной, а не наоборот, поэтому отталкиваемся от именн фалов в тестовой
filenames = readdir(testdir) # названия папок, соответсвующие названиям файлов, со вложенными папками по каждому измерению

Stat = Metrics[] # статистики по всем файлам и измерениям в них

# цикл по всем папкам (соответсвующим именам исходных файлов)
for i in 1:lastindex(filenames)
# i = 1 # для отладки

    measures = readdir("$testdir/$(filenames[i])") # все измерения (папка по каждому) в папке файла с тестовой разметкой

    imetrics = Metrics[]  # вектор метрик по кажому измерению в файле
    # цикл по всем папкам (соответсвующим номеру измерения). Внутри папки стандартные файлы: bounds.csv, pres.csv, tone.csv
    for j in 1:lastindex(measures)
    # j = 3
        # читааем разметку по событиям тонов из тест и реф
        # при этом, если файла (каких-то из вложенных папок, где он лежит) нету в реф, то просто пропускаем и идем дальше
        ref_toneEv = try ReadRefMkp("$refdir/$(filenames[i])/$(measures[j])/tone.csv") catch e continue end
        # ref_toneEv = ReadRefMkp("$refdir/$(filenames[i])/$(measures[j])/tone.csv")
        test_toneEv = ReadRefMkp("$testdir/$(filenames[i])/$(measures[j])/tone.csv")
        ref_bounds = ReadRefMkp("$testdir/$(filenames[i])/$(measures[j])/bounds.csv")

        # определения радиуса сравнения нужно знать частоту дискретизации, поэтому зачитываем хедер соответсвующего файла
        fs = readhdr("$bindir/$(filenames[i]).hdr")[2]

        # сравнение разметок
        rcomp = 0.2 # радиус поиска пары в секундах

        # отправляем на сранение только позиции пиков внутри рабочей зоны
        ref_inbounds = filter(x -> isinbounds(x.pos, [ref_bounds.iwz.pump, ref_bounds.iwz.desc]), ref_toneEv)
        test_inbounds = filter(x -> isinbounds(x.pos, [ref_bounds.iwz.pump, ref_bounds.iwz.desc]), test_toneEv)

        refpos = map(x -> x.pos, ref_inbounds); refbad = map(x -> x.type, ref_inbounds)
        testpos = map(x -> x.pos, test_inbounds); testbad = map(x -> x.type, test_inbounds)

        outcomp = calc_indexpairs(refpos, testpos, floor(rcomp*fs))

        metrics = getmetrics(outcomp, refbad, testbad)
        metrics.filename = filenames[i]; metrics.meas = j
        push!(imetrics, metrics)
        push!(Stat, metrics)
    end

    # среднее по всем измерениям в файле
    fmetrics = getmetrics(imetrics)
    fmetrics.filename = "Total $(filenames[i])"
    push!(Stat, fmetrics) 

end

total = getmetrics(Stat)
total.filename = "Total"
push!(Stat, total)

res = DataFrame(Stat)

CSV.write("results/tonestats_poli.csv", res, delim = ";")