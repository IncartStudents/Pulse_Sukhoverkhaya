using CImGui
using CImGui: ImVec2
using ImPlot
using CSV
using DataFrames
using Gtk
using FileIO
using Images
using Plots
using ImPlot.LibCImGui: ImGuiCond_Always, ImPlotAxisFlags_NoGridLines

include("../src/Renderer.jl")
using .Renderer

include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 
include("../src/readfiles.jl") 
include("../src/my_filt.jl")  

mutable struct PlotElements
    signal::Vector{Float64}
    iSAD::Int64
    iDAD::Int64
    TPpoints::Vector{Int64}
    FPpoints::Vector{Int64}
    FNpoints::Vector{Int64}
    BADpoints::Vector{Int64}
end

struct Data
    tone_stat::DataFrame
    pres_stat::DataFrame
    tone_causes::DataFrame
    pres_causes::DataFrame
end

mutable struct Globals
    data::Data
    tab_item::Int
    selected_item::Int
    needdisplay::Bool
    plotdata::PlotElements

    function Globals()
        df = DataFrame(""=>0)
        data = Data(df, df, df, df)
        tab_item = 1
        selected_item = 0
        needdisplay = false
        plotdata = PlotElements(Float64[], 0, 0, Int[], Int[], Int[], Int[])

        new(data, tab_item, selected_item, needdisplay, plotdata)
    end
end

function SelectFolderButton(v::Globals) # Кнопка загрузки папки с результатами
    if CImGui.Button("Загрузить папку с результатами сравнения")
        fname = open_dialog_native("Выберите папку", action = GtkFileChooserAction.SELECT_FOLDER);
        files = readdir(fname)
        pres_stats = tone_stats = pres_causes = tone_causes = DataFrame(""=>0)
        for i in files
            if length(split(i, "pres")) != 1   # если файл по давлению
                if length(split(i, "statistics")) != 1 # если файл с основными статистиками
                    pres_stats = CSV.read("$fname/$i", DataFrame)
                elseif length(split(i, "causes")) != 1 # если файл с причинами FN
                    pres_causes = CSV.read("$fname/$i", DataFrame)
                end
            elseif length(split(i, "tone")) != 1 # если файл по пульсациям
                if length(split(i, "statistics")) != 1 # если файл с основными статистиками
                    tone_stats = CSV.read("$fname/$i", DataFrame)
                elseif length(split(i, "causes")) != 1 # если файл с причинами FN
                    tone_causes = CSV.read("$fname/$i", DataFrame)
                end
            end
        end
        v.data = Data(tone_stats, pres_stats, tone_causes, pres_causes)
    end
end

function SelectTable(v::Globals) # Вкладки с выбором таблицы
    CImGui.NewLine()
    if CImGui.BeginTabBar("Выбор таблицы", CImGui.ImGuiTabBarFlags_None)
        if CImGui.BeginTabItem("Давление") v.tab_item = 1 
        CImGui.EndTabItem() end
        if CImGui.BeginTabItem("Тоны") v.tab_item = 2 
        CImGui.EndTabItem() end
        CImGui.EndTabBar()
    end
end

function ResultTables(s::Renderer.GR, v::Globals) # Таблицы с результатами
    CImGui.SetNextWindowPos(ImVec2(0, 0))
    CImGui.SetNextWindowSize(ImVec2(s.w, s.h/2))
    CImGui.Begin("Таблица результатов")
        SelectFolderButton(v) # Кнопка загрузки папки с результатами
        FNCausesTable(s, v) # Таблица с причинами FN для конкретной записи
        SelectTable(v) # Вкладки с выбором таблицы

        # Таблицы
        if v.tab_item == 1 stats = v.data.pres_stat # Если выбрана вкладка Давление
        elseif v.tab_item == 2 stats = v.data.tone_stat end # Если выбрана вкладка Тоны

        row,col=size(stats)
        nms=names(stats)

        CImGui.Columns(col, "Заголовки", false)
        CImGui.Separator()
        for i in nms
            CImGui.Text(i)
            CImGui.NextColumn()
        end
        CImGui.Columns(1)
        CImGui.Separator()

        CImGui.BeginChild("##ScrollingRegion", ImVec2(s.w, CImGui.GetFontSize() * 13), false, CImGui.ImGuiWindowFlags_HorizontalScrollbar);
        CImGui.Columns(col,"Метрики")
        for i in 1:row
            for j in nms
                if j == nms[1] && isempty(findall("Sum", string(stats[i,j])))
                    isselected = (i == v.selected_item)
                    CImGui.PushID(i)
                    if CImGui.Selectable(string(stats[i,j]), isselected, CImGui.ImGuiSelectableFlags_SpanAllColumns) 
                        v.selected_item = i 
                        v.needdisplay = true  
                        GeneratePlotData(split(string(stats[i,j]), " ")[1], v) # Генерация данных для графика с разметкой
                    end
                    CImGui.PopID()
                    CImGui.NextColumn()
                else CImGui.Text(string(stats[i,j])); CImGui.NextColumn() end
            end
            CImGui.Separator()
        end
        CImGui.Columns(1)
        CImGui.Separator()
        CImGui.EndChild()

    CImGui.End()
end

function FNCausesTable(s::Renderer.GR, v::Globals) # Таблица с причинами FN для конкретной записи
    CImGui.NewLine()
    CImGui.SameLine(s.w*0.5)
    CImGui.BeginChild("##ScrollingRegion2", ImVec2(s.w*0.5, CImGui.GetFontSize()*4.5), false, CImGui.ImGuiWindowFlags_HorizontalScrollbar)
        if v.tab_item == 1 causes = v.data.pres_causes # Если выбрана вкладка Давление
        elseif v.tab_item == 2 causes = v.data.tone_causes end # Если выбрана вкладка Тоны

        if v.selected_item != 0
            CImGui.Text("Причины отбраковки, повлиявшие на количество FN:")

            _,col=size(causes)
            nms=names(causes)[3:end]

            CImGui.Columns(col-2,"Причины FN")
            CImGui.Separator()
            for i in nms
                CImGui.Text(i)
                CImGui.NextColumn()
            end
            CImGui.Separator()
            for j in nms
                CImGui.Text(string(causes[v.selected_item, j]))
                CImGui.NextColumn()
            end
            CImGui.Columns(1)
            CImGui.Separator()
        end
    CImGui.EndChild()
end

function GeneratePlotData(name, v::Globals) # Генерация данных для графика с разметкой

    if v.tab_item == 1 sigtype = "pres" elseif v.tab_item == 2 sigtype = "tone" end

    binfile = "D:/INCART/Pulse_Data/bin/$name.bin"
    signals, fs, _, _ = readbin(binfile)

    Tone = signals.Tone # пульсации
    Pres = signals.Pres # давление

    # получение валидных сегментов
    vseg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

    # парсинг разметки с результатами сравнения
    tone_compare = "D:/ИНКАРТ/Pulse_Sukhoverkhaya/alg markup/$name.$(sigtype)comp"
    outcomp = parse_compare_results(tone_compare)

    # референтные границы САД-ДАД, внутри которых считаем статистики
    adtablefile = "D:/INCART/Pulse_Data/ad result tables/$(name)_table_ad.txt"
    ad = read_ad(adtablefile)

    if sigtype == "pres" h = v.data.pres_stat[v.selected_item, "meas"]
    elseif sigtype == "tone" h = v.data.tone_stat[v.selected_item, "meas"] end

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

    bounds = get_ad_bounds(Pres[vseg[h].ibeg:vseg[h].iend], ad[h])
    if bounds.isad == 0 || bounds.idad == 0 bsad = vseg[h].ibeg; bdad = vseg[h].iend
    else bsad = bounds.isad+vseg[h].ibeg; bdad = bounds.idad+vseg[h].ibeg end

    if sigtype == "tone"
        seg = Tone[vseg[h].ibeg:vseg[h].iend]

        smoothtone = my_butter(seg, 2, 60, fs, "low") # сглаженный тонов
        ftone = my_butter(smoothtone, 2, 30, fs, "high") # фильтрованный тонов
        fstone = my_butter(abs.(ftone), 2, 10, fs, "low") # огибающая по модулю
    elseif sigtype == "pres"
        seg = Pres[vseg[h].ibeg:vseg[h].iend]

        fsig_smooth = my_butter(seg, 2, 10, fs, "low") # сглаживание
        fstone = my_butter(fsig_smooth, 2, 0.3, fs, "high") # устранение постоянной составляющей
    end

    v.plotdata = PlotElements(fstone, bounds.isad, bounds.idad, 
                                TPpos.-vseg[h].ibeg.+1, FPpos.-vseg[h].ibeg.+1, 
                                FNpos.-vseg[h].ibeg.+1, badpos.-vseg[h].ibeg.+1)
end

function MarkupPlots(s::Renderer.GR, v::Globals) # Графики с разметкой
    CImGui.SetNextWindowPos(ImVec2(0, s.h/2)) 
    CImGui.SetNextWindowSize(ImVec2(s.w, s.h/2))
    CImGui.Begin("Разметка")

        # Пока что просто подгружаем ранее построенные Графики
        # img = load("compare figures/pres/PX11321102817293 1.png")
        # CImGui.Image(pointer_from_objref(img), ImVec2(300,300))

        # if v.needdisplay 
        #     if v.tab_item == 1 
        #         sfold = "pres"
        #         fname = v.data.pres_stat[v.selected_item,"filename"] 
        #         meas = v.data.pres_stat[v.selected_item,"meas"] 
        #     elseif v.tab_item == 2 
        #         sfold = "tone" 
        #         fname = v.data.pres_stat[v.selected_item,"filename"] 
        #         meas = v.data.pres_stat[v.selected_item,"meas"] 
        #     end
        #     img = load("compare figures/$sfold/$fname $meas.png")
        #     display(img)
        #     v.needdisplay = false 
        # end

        # Отрисовка
        if !isempty(v.plotdata.signal)
            ymin = minimum(v.plotdata.signal)
            ymax = maximum(v.plotdata.signal)
            ImPlot.SetNextPlotLimits(0, length(v.plotdata.signal), ymin, ymax, ImGuiCond_Always)
            if ImPlot.BeginPlot("Разметка", C_NULL, C_NULL, ImVec2(s.w, s.h/2.5), flags = ImPlotAxisFlags_NoGridLines)
                ImPlot.PlotLine(v.plotdata.signal)
                ImPlot.SetNextMarkerStyle(ImPlotMarker_Circle, 5, ImPlot.GetColormapColor(1))
                ImPlot.PlotScatter(v.plotdata.TPpoints, v.plotdata.signal[v.plotdata.TPpoints])
                ImPlot.SetNextMarkerStyle(ImPlotMarker_Circle, 5, ImPlot.GetColormapColor(2))
                ImPlot.PlotScatter(v.plotdata.FPpoints, v.plotdata.signal[v.plotdata.FPpoints])
                ImPlot.SetNextMarkerStyle(ImPlotMarker_Circle, 5, ImPlot.GetColormapColor(3))
                ImPlot.PlotScatter(v.plotdata.FNpoints, v.plotdata.signal[v.plotdata.FNpoints])
                ImPlot.PlotLine([v.plotdata.iSAD, v.plotdata.iSAD], [ymin, ymax])
                ImPlot.PlotLine([v.plotdata.iDAD, v.plotdata.iDAD], [ymin, ymax])
            end
        end

    CImGui.End()
end

function ui(s::Renderer.GR, v::Globals)  # Сборка
    CImGui.StyleColorsLight()
    ResultTables(s::Renderer.GR, v::Globals) # Таблицы с результатами
    MarkupPlots(s::Renderer.GR, v::Globals) # Графики с разметкой
end

function show_gui() # Main
    state = Globals();
    size = Renderer.GR();
    Renderer.render(
        ()->ui(size, state),
        width=1700,
        height=1600,
        title="",
        v = size
    )
end

show_gui();