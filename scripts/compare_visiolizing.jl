using Plots
using DataFrames
using CSV

include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 
include("../src/readfiles.jl") 
include("../src/my_filt.jl")   

struct Stats
    filename::String # имя файла
    meas::Int64      # номер измерения
    TP::Int64
    FP::Int64
    FN::Int64
    Se::Float64
    PVP::Float64
end

struct FNCauses
    filename::String # имя файла
    meas::Int64      # номер измерения
    FN::Int64
    causes::Vector{Int}
end

function figures_n_stats(allbins, sigtype)
    if sigtype != "tone" && sigtype != "pres" return "Invalid signal type!" end

    statistics = fill(Stats[], lastindex(allbins))
    causesstat = fill(FNCauses[], lastindex(allbins))

    for FILEN in 1:lastindex(allbins)

        nm = split(allbins[FILEN],".")[1]

        binfile = "D:/INCART/Pulse_Data/bin/$nm.bin"
        signals, fs, _, _ = readbin(binfile)

        Tone = signals.Tone # пульсации
        Pres = signals.Pres # давление

        # получение валидных сегментов
        vseg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

        # парсинг разметки с результатами сравнения
        tone_compare = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/$nm.$(sigtype)comp"
        outcomp = parse_compare_results(tone_compare)

        # референтные границы САД-ДАД, внутри которых считаем статистики
        adtablefile = "D:/INCART/Pulse_Data/ad result tables/$(nm)_table_ad.txt"
        ad = read_ad(adtablefile)

        # Se = fill(0.0, lastindex(outcomp))
        # PVP = fill(0.0, lastindex(outcomp))

        statistics[FILEN] = fill(Stats("", 0, 0, 0, 0, 0.0, 0.0), length(outcomp))
        causesstat[FILEN] = fill(FNCauses("", 0, 0, Int[]), length(outcomp))

        for h in 1:lastindex(outcomp)
            # # позиции реф разметки 
            # ref_pos = map(x -> x.iref != -1 ? x.iref : NaN, outcomp[h])
            # ref_pos = ref_pos[findall(x -> !isnan(x), ref_pos)]

            # # позиции тест разметки (без отбракованных)
            # test_pos = map(x -> x.itest != -1 && x.bad == 0 ? x.itest : NaN , outcomp[h])
            # test_pos = test_pos[findall(x -> !isnan(x), test_pos)]

            # TP позиции
            TPpos = map(x -> x.code == "ee" ? x.itest : NaN, outcomp[h])
            TPpos = TPpos[findall(x -> !isnan(x), TPpos)]

            # FP позиции
            FPpos = map(x -> x.code == "oe" ? x.itest : NaN, outcomp[h])
            FPpos = FPpos[findall(x -> !isnan(x), FPpos)]

            # FN позиции
            FNpos = map(x -> if x.code == "eo" x.iref elseif x.code == "eb" x.itest else NaN end, outcomp[h])
            FNpos = FNpos[findall(x -> !isnan(x), FNpos)]

            # Забракованные
            badpos = map(x -> x.bad != 0 ? x.itest : NaN, outcomp[h])
            badpos = badpos[findall(x -> !isnan(x), badpos)]

            # расшифровка причин браковки
            # bad1 = map(x -> my_bitand(x.bad,2^0)!=0 ? x.itest : NaN, outcomp[h])
            # bad1 = bad1[findall(x -> !isnan(x), bad1)]
            # bad2 = map(x -> my_bitand(x.bad,2^1)!=0 ? x.itest : NaN, outcomp[h])
            # bad2 = bad2[findall(x -> !isnan(x), bad2)]
            # bad3 = map(x -> my_bitand(x.bad,2^2)!=0 ? x.itest : NaN, outcomp[h])
            # bad3 = bad3[findall(x -> !isnan(x), bad3)]
            # bad4 = map(x -> my_bitand(x.bad,2^3)!=0 ? x.itest : NaN, outcomp[h])
            # bad4 = bad4[findall(x -> !isnan(x), bad4)]

            # статистики БЕЗ УЧЕТА ОТБРАКОВАННЫХ СОБЫТИЙ (в случае с давлением, кроме нахождения слева от точки начала спуска)
            # считаем только внутри референтных границ САД-ДАД

            bounds = get_ad_bounds(Pres[vseg[h].ibeg:vseg[h].iend], ad[h])
            if bounds.isad == 0 || bounds.idad == 0 bsad = vseg[h].ibeg; bdad = vseg[h].iend
            else bsad = bounds.isad+vseg[h].ibeg; bdad = bounds.idad+vseg[h].ibeg end

            TP = length(findall(x -> x.code == "ee" && x.itest >= bsad && x.itest <= bdad, outcomp[h]))
            FP = length(findall(x -> x.code == "oe" && x.itest >= bsad && x.itest <= bdad, outcomp[h]))
            FN = length(findall(x -> (x.code == "eo" || x.code == "eb") && x.itest >= bsad && x.itest <= bdad, outcomp[h]))

            Se = round((TP/(TP+FN))*100, digits=2)
            PVP = round((TP/(TP+FP))*100, digits=2)
            statistics[FILEN][h] = Stats(nm, h, TP, FP, FN, Se, PVP)

            # поиск причины или набора причин отбраковки для FN внутри реф САД-ДАД (если они есть)
            css = map(x -> (x.code == "eo" || x.code == "eb") && x.itest >= bsad && x.itest <= bdad
                            && x.bad != 0 ? errorcode(x.bad, 10) : -1, outcomp[h])
            css = css[findall(x -> x != -1, css)]
            allcss = Int64[]
            for i in css allcss = vcat(allcss, i) end

            causesstat[FILEN][h] = FNCauses(nm, h, FN, allcss)

            # рисование
            if sigtype == "tone"
                seg = Tone[vseg[h].ibeg:vseg[h].iend]

                smoothtone = my_butter(seg, 2, 60, fs, Lowpass) # сглаженный тонов
                ftone = my_butter(smoothtone, 2, 30, fs, Highpass) # фильтрованный тонов
                fstone = my_butter(abs.(ftone), 2, 10, fs, Lowpass) # огибающая по модулю
            elseif sigtype == "pres"
                seg = Pres[vseg[h].ibeg:vseg[h].iend]

                fsig_smooth = my_butter(seg, 2, 10, fs, Lowpass) # сглаживание
                fstone = my_butter(fsig_smooth, 2, 0.3, fs, Highpass) # устранение постоянной составляющей
            end

            plot(fstone, label = "$sigtype")
            # scatter!(ref_pos.-vseg[h].ibeg.+1, fstone[ref_pos.-vseg[h].ibeg.+1], markersize = 2, label = "ref")
            # scatter!(test_pos.-vseg[h].ibeg.+1, fstone[test_pos.-vseg[h].ibeg.+1], markersize = 2, label = "test")

            scatter!(TPpos.-vseg[h].ibeg.+1, fstone[TPpos.-vseg[h].ibeg.+1], markersize = 3, label = "TP", color = :green)
            scatter!(FPpos.-vseg[h].ibeg.+1, fstone[FPpos.-vseg[h].ibeg.+1], markersize = 3, label = "FP", color = :blue)
            scatter!(FNpos.-vseg[h].ibeg.+1, fstone[FNpos.-vseg[h].ibeg.+1], markersize = 3, label = "FN", color = :yellow)
            scatter!(badpos.-vseg[h].ibeg.+1, fstone[badpos.-vseg[h].ibeg.+1], markersize = 1, label = "bad", color = :red)
            # scatter!(bad1.-vseg[h].ibeg.+1, fstone[bad1.-vseg[h].ibeg.+1].+100, markersize = 3, markershape = :star, label = "s1", color = :black)
            # scatter!(bad2.-vseg[h].ibeg.+1, fstone[bad2.-vseg[h].ibeg.+1].+200, markersize = 3, markershape = :star, label = "s2", color = :green)
            # scatter!(bad3.-vseg[h].ibeg.+3, fstone[bad3.-vseg[h].ibeg.+1].+300, markersize = 3, markershape = :star, label = "s3", color = :pink)
            # scatter!(bad4.-vseg[h].ibeg.+1, fstone[bad4.-vseg[h].ibeg.+1].+400, markersize = 3, markershape = :star, label = "s4", color = :blue)
            plot!([[bounds.isad, bounds.isad],[bounds.idad, bounds.idad]],[[minimum(fstone), maximum(fstone)],[minimum(fstone), maximum(fstone)]], label = ["SAD" "DAD"], color = :black)
            title!("Se = $(Se), PVP = $(PVP)")

            savefig("compare figures/$sigtype/$nm $h.png")
        end
    end

    return statistics, causesstat
end

function savestatistics(statistics, tone_causesstat, sigtype)
    if sigtype != "tone" && sigtype != "pres" return "Invalid signal type!" end

    allTP = allFP = allFN = 0
    allSemean = allPVPmean = 0.0

    # основная статистика
    allstats = Stats[] 
    k = 0
    for i in statistics
        sumTP = sum(map(x -> x.TP, i)) 
        sumFP = sum(map(x -> x.FP, i)) 
        sumFN = sum(map(x -> x.FN, i)) 
        meanSe = round(mean(map(x -> x.Se, i)), digits=2)
        meanPVP = round(mean(map(x -> x.PVP, i)), digits=2)

        allstats = vcat(vcat(allstats, i), Stats("Sum TP, FP, FN.\n Mean Se, PVP.", 0, sumTP, sumFP, sumFN, meanSe, meanPVP))

        allTP += sumTP
        allFP += sumFP
        allFN += sumFN
        allSemean += meanSe
        allPVPmean += meanPVP
        k += 1
    end

    allstats = vcat(allstats, Stats("All Sum TP, FP, FN.\n All Mean Se, PVP.", 0, 
                                    allTP, allFP, allFN, round(allSemean/k, digits=2), round(allPVPmean/k, digits=2)))

    # статистика по причинам ошибок (FN в частности)
    tone_causes = Int[]
    for i in tone_causesstat for j in i tone_causes = vcat(tone_causes, j.causes) end end
    tone_causes = sort(unique(tone_causes))

    numofcauses = Vector{Vector{Int}}[]
    headerres = NTuple[]
    FNtotal = 0
    causestotal = fill(0, length(tone_causes))

    for i in tone_causesstat   # по каждому файлу
    FNs = 0
    cntrec = fill(0, length(tone_causes))
    for k in i             # по кадому измерению
        cnt = fill(0, length(tone_causes))
        for p in k.causes
            if !isempty(p)
                icauses = findall(x -> x==p, tone_causes)[1]
                if !isempty(icauses)
                    for j in tone_causes[icauses]
                        ind = findall(x -> x==j[1], tone_causes)[1]
                        cnt[ind] += 1
                    end
                end
            end
        end
        FNs += k.FN
        cntrec += cnt
        headerres = [headerres..., (filename = k.filename, meas = k.meas, FN = k.FN)]
        numofcauses = [numofcauses..., Tuple(cnt)]
    end
    FNtotal += FNs
    causestotal += cntrec
    headerres = [headerres..., (filename = "Sum", meas = 0, FN = FNs)]
    numofcauses = [numofcauses..., Tuple(cntrec)]
    end
    headerres = [headerres..., (filename = "Total", meas = 0, FN = FNtotal)]
    numofcauses = [numofcauses..., Tuple(causestotal)]

    df = headerres |> DataFrame
    ndf = numofcauses |> DataFrame
    currnames = names(ndf)
    for i in 1:lastindex(currnames) rename!(ndf, currnames[i] => Symbol.(tone_causes[i])) end


    stats_df = allstats |> DataFrame
    causes_df = hcat(df, ndf)

    CSV.write("results/$(sigtype)_statistics.csv", stats_df, delim = ";")
    CSV.write("results/$(sigtype)_causes.csv", causes_df, delim = ";")
end

dir = "D:/INCART/Pulse_Data/bin"
files = readdir(dir)
allbins = files[findall(x -> split(x, ".")[end] == "bin", files)]

tone_statistics, tone_causesstat = figures_n_stats(allbins, "tone")
pres_statistics, pres_causesstat = figures_n_stats(allbins, "pres")

savestatistics(tone_statistics, tone_causesstat, "tone")
savestatistics(pres_statistics, pres_causesstat, "pres")

