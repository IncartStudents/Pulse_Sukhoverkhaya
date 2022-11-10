using CImGui
using CImGui: ImVec2, ImVec4
using ImPlot
using CSV
using DataFrames
using Gtk
using FileIO
using Images
using Plots
using ImPlot.LibCImGui: ImGuiCond_Once, ImGuiCond_Always, ImPlotAxisFlags_NoGridLines

include("../src/Renderer.jl")
using .Renderer

include("../src/help_func.jl")
include("../src/readfiles.jl") 
include("../src/my_filt.jl")  
include("../src/refmkpguifunctions.jl") 

# include(joinpath(pathof(ImPlot), "..", "..", "demo", "implot_demo.jl"))
# show_demo()


struct Signal    # Нужные каналы и частота дискретизации, вытянутые из бинаря + границы валидных сегментов
    ECG::Vector{Float64}
    Pres::Vector{Float64}
    Tone::Vector{Float64}
    fs::Int64
end

mutable struct PlotBounds   # Данные границ сегмента, рабочей зоны и АД ( рз и ад - на накачке и на спуске)
    AD::NamedTuple{(:pump, :desc), NTuple{2, Bounds}}
    workreg::NamedTuple{(:pump, :desc), NTuple{2, Bounds}}
    segbounds::Bounds
end

struct Mkp       # Разметки тонов и пульсаций (каждая в векторе своих структур)
    Pres::Vector{PresGuiMkp}
    Tone::Vector{ToneGuiMkp}
    bounds::PlotBounds
end

struct Area      # Границы нарисованной прямогольной зоны + на каком графике нарисована
    begpos::ImVec2
    endpos::ImVec2
    whichplot::String
end 

mutable struct PlotElem    # Данные для построения графика по подному сигналу
    sig::Vector{Float64}

    ibegs::Vector{Int64}
    iends::Vector{Int64}
    type::Vector{Int64} # 0 - незначимые (отсутствие метки), 1 - значимые, 2 - шумовые события

    limits::ImPlotLimits
end

mutable struct PlotData     # Данные для построения графиков по всем сигналам
    ECG::PlotElem
    rawPres::PlotElem
    Pres::PlotElem
    Tone::PlotElem
end

mutable struct PlotLinks    # Для синхронизации масштабов графиков
    xmin::Ref
    xmax::Ref
    ymin::Ref
    ymax::Ref
    linkx::Bool
    linky::Bool
end

mutable struct PlotID       # ID каждого графика
    ECG::Int64
    rawPres::Int64
    Pres::Int64
    Tone::Int64
end

mutable struct Globals      # "Глобальные" переменные (главная структура)
    filename::String
    signal::Signal
    markup::Mkp
    plotbounds::PlotBounds
    dataforplotting::PlotData
    plotlinks::PlotLinks
    combo_item::Int64
    typecombo_item::Int64
    mode::Int64
    selected_peaks::Vector{Int64}
    area::Area
    pt0::ImVec2
    plotsid::PlotID
    allfiles::Vector{String}
    tabledata::Vector{Tuple}
    fold::String
    selecteditem::Int64
    allbases::Vector{String}
    selectedbase::Int64
    isguistarted::Bool

    function Globals()
        filename = ""
        signal = Signal(Float64[], Float64[], Float64[], 0)
        bnds0 = Bounds(0,0)
        plotbounds = PlotBounds((pump = bnds0, desc = bnds0), (pump = bnds0, desc = bnds0), bnds0)
        markup = Mkp(PresGuiMkp[], ToneGuiMkp[], plotbounds)
        tup = PlotElem(Float64[], Int64[], Int64[], Int64[], ImPlotLimits(ImPlotRange(0.0,0.0), ImPlotRange(0.0,0.0)))
        dataforplotting = PlotData(tup, tup, tup, tup)
        plotlinks = PlotLinks(Ref(0.0), Ref(1.0), Ref(0.0), Ref(1.0), true, false)
        combo_item = 1
        typecombo_item = 1
        mode = 0
        selected_peaks = Int64[]
        area = Area(ImVec2(0.0,0.0), ImVec2(0.0,0.0), "")
        pt0 = ImVec2(0.0,0.0)
        plotsid = PlotID(-2,-1,1,2)
        allfiles = String[]
        tabledata = Tuple[]
        fold = ""
        selecteditem = 1
        allbases = String[]
        selectedbase = 0
        isguistarted = false

        new(filename, signal, markup, plotbounds, dataforplotting, 
            plotlinks, combo_item, typecombo_item, mode, 
            selected_peaks, area, pt0, plotsid, allfiles, tabledata, fold,
            selecteditem, allbases, selectedbase, isguistarted)
    end
end

function GeneratePlotData(v::Globals)   # Генерация новых данных для построения графиков

    # референтные значения САД-ДАД и рабочей зоны + границы сегмента
    ad = v.markup.bounds.AD             # по амплитуде
    wz = v.markup.bounds.workreg        # по амплитуде
    vseg = v.markup.bounds.segbounds    # по отсчетам

    # поиск границ САД-ДАД и рабочей зоны
    pump_adbounds = get_bounds(v.signal.Pres[vseg.ibeg:vseg.iend], ad.pump, true)
    desc_adbounds = get_bounds(v.signal.Pres[vseg.ibeg:vseg.iend], ad.desc, false)

    pump_wzbounds = get_bounds(v.signal.Pres[vseg.ibeg:vseg.iend], wz.pump, true)
    desc_wzbounds = get_bounds(v.signal.Pres[vseg.ibeg:vseg.iend], wz.desc, false)

    v.plotbounds = PlotBounds((pump = pump_adbounds, desc = desc_adbounds), (pump = pump_wzbounds, desc = desc_wzbounds), vseg)

    # ЭКГ
    ECG = v.signal.ECG[vseg.ibeg:vseg.iend]

    # тоны
    seg = v.signal.Tone[vseg.ibeg:vseg.iend]

    smoothtone = my_butter(seg, 2, 60, v.signal.fs, "low")          # сглаженный тонов
    ftone = my_butter(smoothtone, 2, 30, v.signal.fs, "high")       # фильтрованный тонов
    tone_sig = my_butter(abs.(ftone), 2, 10, v.signal.fs, "low")    # огибающая по модулю

    tone_peaks = map(x -> x.pos, v.markup.Tone)

    # пульсации
    seg = v.signal.Pres[vseg.ibeg:vseg.iend]

    fsig_smooth = my_butter(seg, 2, 10, v.signal.fs, "low") # сглаживание
    pres_sig = my_butter(fsig_smooth, 2, 0.3, v.signal.fs, "high") # устранение постоянной составляющей

    pres_begs = map(x -> x.ibeg, v.markup.Pres)
    pres_ends = map(x -> x.iend, v.markup.Pres)

    # расставление типов меток по умолчанию (все, что за границами рабочей зоны - в незначимые, остальные - (пока) в значимые)
    tone_peaks_type = map(x -> x < pump_wzbounds.ibeg || (x > pump_wzbounds.iend && x < desc_wzbounds.ibeg) || x > desc_wzbounds.iend ? 0 : 1, tone_peaks)
    pres_peaks_type = map(x -> x < pump_wzbounds.ibeg || (x > pump_wzbounds.iend && x < desc_wzbounds.ibeg) || x > desc_wzbounds.iend ? 0 : 1, pres_begs)

    ECGtup = PlotElem(ECG, Int64[], Int64[], Int64[], ImPlotLimits(ImPlotRange(0.0, length(ECG) |> Float64),ImPlotRange(minimum(ECG), maximum(ECG))))
    RawPrestup = PlotElem(seg, Int64[], Int64[], Int64[], ImPlotLimits(ImPlotRange(0.0, length(seg) |> Float64),ImPlotRange(minimum(seg), maximum(seg))))
    Prestup = PlotElem(pres_sig, pres_begs, pres_ends, pres_peaks_type, ImPlotLimits(ImPlotRange(0.0, length(pres_sig) |> Float64),ImPlotRange(minimum(pres_sig), maximum(pres_sig))))
    Tonetup = PlotElem(tone_sig, tone_peaks, tone_peaks, tone_peaks_type, ImPlotLimits(ImPlotRange(0.0, length(tone_sig) |> Float64),ImPlotRange(minimum(tone_sig), maximum(tone_sig))))

    v.dataforplotting = PlotData(ECGtup, RawPrestup, Prestup, Tonetup)
end

function ReadData(v::Globals) # Чтение данных из нового выбранного файла

    # чтение сигнала из бинаря
    signals, fs, _, _ = readbin(v.fold*"/"*v.allfiles[v.combo_item])

    ECG = signals[1]     # ЭКГ
    Tone = signals.Tone  # пульсации
    Pres = signals.Pres  # давление

    # чтение разметки: если есть в реф - из реф, если нету - из тест
    # если есть - зачитываем сразу 
    srcdir = "formatted alg markup"
    dstdir = "ref markup"
    basename = split(v.fold, "/")[end]
    name = v.allfiles[v.combo_item]
    measure = v.selecteditem

    dir = try readdir("$dstdir/$basename/$name/$measure"); dstdir catch e srcdir end

    # зачитывание разметки из реф
    Pres_mkp = ReadRefMkp("$dir/$basename/$name/$measure/pres.csv")
    Tone_mkp = ReadRefMkp("$dir/$basename/$name/$measure/tone.csv")
    bnds = ReadRefMkp("$dir/$basename/$name/$measure/bounds.csv")

    v.markup = Mkp(Pres_mkp, Tone_mkp, PlotBounds(bnds.ad, bnds.wz, bnds.segm))
    v.signal = Signal(ECG, Pres, Tone, fs)

    GeneratePlotData(v)
end

function SaveRefMarkupButton(v::Globals) # Кнопка сохранения исправленной референтной разметки тонов и пульсаций + границ сегмента, рабочей зоны и АД для текущего измерения выбранного файла
    if !isempty(v.allfiles)
        CImGui.SameLine(CImGui.GetWindowContentRegionWidth()-250)
        if CImGui.Button("Сохранить разметку")
            dirname0 = "ref markup"
            try readdir(dirname0)
            catch e mkdir(dirname0) end
            basename = split(v.fold, "/")[end]
            dirname0 = "ref markup/$basename"
            try readdir(dirname0)
            catch e mkdir(dirname0) end
            dirname0 = "ref markup/$basename/$(v.allfiles[v.combo_item])"
            try readdir(dirname0)
            catch e mkdir(dirname0) end
            dirname = "ref markup/$basename/$(v.allfiles[v.combo_item])/$(v.tabledata[v.selecteditem][1][2])"
            try readdir(dirname)
            catch e mkdir(dirname) end

            Pres = v.dataforplotting.rawPres.sig

            pres_markup = map((x,y,z) -> PresGuiMkp(x, y, z), v.dataforplotting.Pres.ibegs, v.dataforplotting.Pres.iends, v.dataforplotting.Pres.type)
            tone_markup = map((x,y) -> ToneGuiMkp(x, y), v.dataforplotting.Tone.ibegs, v.dataforplotting.Tone.type)
            ad = (pump = AD(round(Int, Pres[v.plotbounds.AD.pump.iend]), round(Int, Pres[v.plotbounds.AD.pump.ibeg])), 
                    desc = AD(round(Int, Pres[v.plotbounds.AD.desc.ibeg]), round(Int, Pres[v.plotbounds.AD.desc.iend])))
            wz = (pump = Bounds(round(Int, Pres[v.plotbounds.workreg.pump.iend]), round(Int, Pres[v.plotbounds.workreg.pump.ibeg])), 
                    desc = Bounds(round(Int, Pres[v.plotbounds.workreg.desc.ibeg]), round(Int, Pres[v.plotbounds.workreg.desc.iend])))
            segbounds = v.plotbounds.segbounds

            SaveRefMarkup("$dirname/pres.csv", pres_markup) 
            SaveRefMarkup("$dirname/tone.csv", tone_markup) 
            SaveRefMarkup("$dirname/bounds.csv", segbounds, ad, wz) 

            MakeTableData(v)
        end
    end
end

function FilenamesTable(v::Globals)
    if !isempty(v.allfiles)
        CImGui.NewLine()
        CImGui.BeginChild("##filenames_header", ImVec2(CImGui.GetWindowContentRegionWidth(), CImGui.GetTextLineHeightWithSpacing()*1.3))
        CImGui.Columns(1, "Заголовк")
        CImGui.Separator()
        CImGui.TextColored(ImVec4(0.45, 0.7, 0.80, 1.00), "Имена файлов базы")
        CImGui.Columns(1)
        CImGui.Separator()
        CImGui.EndChild()

        CImGui.BeginChild("##filenames_scrollingregion", (0, 500))
        CImGui.Columns(1, "Имена файлов")
        for i in 1:lastindex(v.allfiles)
            CImGui.PushID(i)
            if CImGui.Selectable(v.allfiles[i], i == v.combo_item)
                v.combo_item = i
                ReadData(v)
                MakeTableData(v) 
                GeneratePlotData(v)
                ChangePlotsID(v)
            end
            CImGui.PopID()
            CImGui.NextColumn()
            CImGui.Separator()
        end
        CImGui.Columns(1)
        CImGui.EndChild()

    end
end

function MakeTableData(v::Globals) # чтение реф границ АД для всех измерений в файле
    dir0 = "ref markup"
    dir = "formatted alg markup"
    basename = split(v.fold, "/")[end]
    name = v.allfiles[v.combo_item]
    allmesfiles = readdir("$dir/$basename/$name")
    v.tabledata = fill(("Номер измерения" => 0, "САД реф." => 0, "ДАД реф." => 0), length(allmesfiles))
    for i in 1:lastindex(allmesfiles)
        bnds = try ReadRefMkp("$dir0/$basename/$name/$(allmesfiles[i])/bounds.csv")
                catch e ReadRefMkp("$dir/$basename/$name/$(allmesfiles[i])/bounds.csv") end
        v.tabledata[i] = ("Номер измерения" => i, "САД реф." => bnds.ad.desc.ibeg, "ДАД реф." => bnds.ad.desc.iend)
    end
end

function MeasuresTable(v::Globals)   # Таблица с границами АД для всех измерений выбранного файла
    if !isempty(v.allfiles)
        if isempty(v.tabledata) MakeTableData(v) end
        CImGui.NewLine()
        names = map(x -> x[1], v.tabledata[1])
        col = length(names)
        CImGui.BeginChild("##header", ImVec2(CImGui.GetWindowContentRegionWidth(), CImGui.GetTextLineHeightWithSpacing()*1.3))
        CImGui.Columns(col, "Заголовки")
        CImGui.Separator()
        for i in names
            CImGui.TextColored(ImVec4(0.45, 0.7, 0.80, 1.00), i)
            CImGui.NextColumn()
        end
        CImGui.Columns(1)
        CImGui.Separator()
        CImGui.EndChild()

        CImGui.BeginChild("##scrollingregion", (0, 300))
        CImGui.Columns(col, "Измерения и границы АД")
        for i in 1:lastindex(v.tabledata)
            for j in v.tabledata[i]
                CImGui.PushID(i)
                if CImGui.Selectable(string(j[2]), i == v.selecteditem, CImGui.ImGuiSelectableFlags_SpanAllColumns)
                    v.selecteditem = i 
                    ReadData(v)
                    GeneratePlotData(v)
                    ChangePlotsID(v)
                end
                CImGui.PopID()
                CImGui.NextColumn()
            end
            CImGui.Separator()
        end
        CImGui.Columns(1)
        CImGui.EndChild()
    end
end

function MenuWindow(v::Globals)   # Окно меню
    CImGui.Begin("Меню")
        BasesTable(v)
        FilenamesTable(v)
        MeasuresTable(v)
    CImGui.End()
end

function PlotLine(id, args...; label = false) # Функция рисования сигнала либо вертикальных границ
    CImGui.PushID(id)
    if length(args) == 1
        if label == false ImPlot.PlotLine(args[1])
        else ImPlot.PlotLine(args[1], label_id = label) end
    elseif length(args) == 2
        if label == false ImPlot.PlotLine(args[1], args[2])
        else ImPlot.PlotLine(args[1], args[2], label_id = label) end
    end
    CImGui.PopID()
end

function Scatter(id, x, y, markerstyle, markersize, label) # Функция нанесения разметки
    CImGui.PushID(id)
    ImPlot.SetNextMarkerStyle(markerstyle, markersize)
    if label == false ImPlot.PlotScatter(x, y)
    else ImPlot.PlotScatter(x, y, label_id = label) end
    CImGui.PopID()
end

function isin(point::ImVec2, searchbox::NamedTuple{(:xmin, :xmax, :ymin, :ymax), NTuple{4, Float64}}) # Находится ли точка внутри указанной прямоугольной зоны
    return (point.x > searchbox.xmin && point.x < searchbox.xmax && point.y > searchbox.ymin && point.y < searchbox.ymax)
end

function InsideArea(v::Globals, whichplot)  # Действия с точками внутри выбранной прямоугольной зоны
    if !CImGui.IsMouseDown(0) && v.area.begpos != v.area.endpos
        if whichplot == "tone"
            allbegs = v.dataforplotting.Tone.ibegs
            allends = v.dataforplotting.Tone.iends
            alltypes = v.dataforplotting.Tone.type
            sig = v.dataforplotting.Tone.sig
        elseif whichplot == "pulse"
            allbegs = v.dataforplotting.Pres.ibegs
            allends = v.dataforplotting.Pres.iends
            alltypes = v.dataforplotting.Pres.type
            sig = v.dataforplotting.Pres.sig
        end
        
        xmin = minimum([v.area.begpos.x, v.area.endpos.x]) |> Float64
        xmax = maximum([v.area.begpos.x, v.area.endpos.x]) |> Float64
        ymin = minimum([v.area.begpos.y, v.area.endpos.y]) |> Float64
        ymax = maximum([v.area.begpos.y, v.area.endpos.y]) |> Float64

        searchbox = (xmin = xmin, xmax = xmax, 
                        ymin = ymin, ymax = ymax)

        if v.mode == 3 # ретипизация меток в области
            newalltypes = map((x,y) -> isin(ImVec2(x,sig[x]), searchbox) ? (v.typecombo_item-1) : y, allbegs, alltypes)
            if whichplot == "tone" v.dataforplotting.Tone.type = newalltypes
            elseif whichplot == "pulse" v.dataforplotting.Pres.type = newalltypes end
        elseif v.mode == 4  # удаление меток в области
            allin = map(x -> isin(ImVec2(x,sig[x]), searchbox) ? true : false, allbegs)
            ind = findall(allin)
            allbegs = deleteat!(allbegs, ind); allends = deleteat!(allends, ind); 
            alltypes = deleteat!(alltypes, ind)
            if whichplot == "tone"
                v.dataforplotting.Tone.ibegs = allbegs
                v.dataforplotting.Tone.iends = allends
                v.dataforplotting.Tone.type = alltypes
            elseif whichplot == "pulse"
                v.dataforplotting.Pres.ibegs = allbegs
                v.dataforplotting.Pres.iends = allends
                v.dataforplotting.Pres.type = alltypes
            end
        end

    end
end

function MouseClick(v::Globals, whichplot)  # Ответ на определенные клики мышью по графику
    if ImPlot.IsPlotHovered() && CImGui.IsMouseClicked(0) && unsafe_load(CImGui.GetIO().KeyCtrl)  # точка (Ctrl + клик левой кнопкой) 
        pt = ImPlot.GetPlotMousePos()
        if whichplot == "tone"
            allbegs = v.dataforplotting.Tone.ibegs
            allends = v.dataforplotting.Tone.iends
            alltypes = v.dataforplotting.Tone.type
            sig = v.dataforplotting.Tone.sig
        elseif whichplot == "pulse"
            allbegs = v.dataforplotting.Pres.ibegs
            allends = v.dataforplotting.Pres.iends
            alltypes = v.dataforplotting.Pres.type
            sig = v.dataforplotting.Pres.sig
        end

        # чтобы реже случались промахи, ищем в окрестностях позиции курсора во время клика максимум сигнала
        mp = round(Int, pt.x)
        sr = 100
        newpeak = argmax(sig[mp-sr:mp+sr])[1] + (mp-sr)

        if v.mode == 3     # добавить или ретипизировать метку
            dscale = 100.0
            searchbox = (xmin = newpeak-dscale, xmax = newpeak+dscale, 
                        ymin = sig[newpeak]-dscale, ymax = sig[newpeak]+dscale)
            retyped = false
            for i in allbegs
                if isin(ImVec2(i, sig[i]), searchbox) # если выбрана существующая метка - ретипизируем
                    ind = findfirst(x -> x==i, allbegs)
                    alltypes[ind] = v.typecombo_item - 1
                    if whichplot == "tone"
                        v.dataforplotting.Tone.ibegs = allbegs
                        v.dataforplotting.Tone.iends = allends
                        v.dataforplotting.Tone.type = alltypes
                    elseif whichplot == "pulse"
                        v.dataforplotting.Pres.ibegs = allbegs
                        v.dataforplotting.Pres.iends = allends
                        v.dataforplotting.Pres.type = alltypes
                    end
                    retyped = true
                end
            end
            if !retyped && whichplot == "tone" # в противном случае устанавливаем новую метку
                push!(allbegs, newpeak); push!(alltypes, v.typecombo_item-1)
                v.dataforplotting.Tone.ibegs = allbegs
                v.dataforplotting.Tone.iends = allbegs
                v.dataforplotting.Tone.type = alltypes
            end
        elseif v.mode == 4 # удалить метку
            dscale = 100.0
            searchbox = (xmin = pt.x-dscale, xmax = pt.x+dscale, 
                        ymin = pt.y-dscale, ymax = pt.y+dscale)
            for i in allbegs
                if isin(ImVec2(i, sig[i]), searchbox)
                    ind = findfirst(x -> x==i, allbegs)
                    if whichplot == "tone"
                        deleteat!(allbegs, ind); deleteat!(alltypes, ind)
                        v.dataforplotting.Tone.ibegs = allbegs
                        v.dataforplotting.Tone.iends = allbegs
                        v.dataforplotting.Tone.type = alltypes
                    elseif whichplot == "pulse"  ### НЕПРАВИЛЬНО РАБОТАЕТ
                        deleteat!(allbegs, ind)
                        deleteat!(allends, ind)
                        deleteat!(alltypes, ind)
                        v.dataforplotting.Pres.ibegs = allbegs
                        v.dataforplotting.Pres.iends = allends
                        v.dataforplotting.Pres.type = alltypes
                    end
                end
            end
        end
    end

    # при нажатии кнопки мыши, сначала разово true становится IsMouseClicked и затем, при удержании, IsMouseDown,
    # поэтому в момент начала удержания (клика) берём позицию мыши
    if ImPlot.IsPlotHovered() && CImGui.IsMouseClicked(0) && unsafe_load(CImGui.GetIO().KeyAlt)
        ipt = ImPlot.GetPlotMousePos()
        v.pt0 = ImVec2(ipt.x, ipt.y)
        v.area = Area(ImVec2(0.0,0.0), ImVec2(0.0,0.0), whichplot)
    end
    if ImPlot.IsPlotHovered() && CImGui.IsMouseDown(0) && unsafe_load(CImGui.GetIO().KeyAlt) # область
        pt = ImPlot.GetPlotMousePos()
        v.area = Area(ImVec2(v.pt0.x, v.pt0.y), ImVec2(pt.x, pt.y), whichplot)
    end

    if ImPlot.IsPlotHovered() && CImGui.IsMouseDown(0)  # захват (для передвижения границ) (курсор наведен на любой график, левая клавиша мыши зажата)
        pt = ImPlot.GetPlotMousePos()
        r = 0.015 * (v.plotbounds.segbounds.iend-v.plotbounds.segbounds.ibeg)
        if v.mode == 1 # двигаем границы рабочей зоны
            leftpump = v.plotbounds.workreg.pump.ibeg; rightpump = v.plotbounds.workreg.pump.iend
            leftdesc = v.plotbounds.workreg.desc.ibeg; rightdesc = v.plotbounds.workreg.desc.iend
            if abs(pt.x-leftpump) <= r # двигаем левую границу на накачке
                leftpump = round(Int, pt.x)
            elseif abs(pt.x-leftdesc) <= r # двигаем левую границу на спуске
                leftdesc = round(Int, pt.x) 
            elseif abs(pt.x-rightpump) <= r #двигаем правую границу на накачке
                rightpump = round(Int, pt.x) 
            elseif abs(pt.x-rightdesc) <= r #двигаем правую границу на спуске
                rightdesc = round(Int, pt.x) 
            end
            v.plotbounds.workreg = (pump = Bounds(leftpump, rightpump), desc = Bounds(leftdesc, rightdesc))
        elseif v.mode == 2 # двигаем границы реф АД
            leftpump = v.plotbounds.AD.pump.ibeg; rightpump = v.plotbounds.AD.pump.iend
            leftdesc = v.plotbounds.AD.desc.ibeg; rightdesc = v.plotbounds.AD.desc.iend
            if abs(pt.x-leftpump) <= r          # двигаем левую границу на накачке
                leftpump = round(Int, pt.x)
            elseif abs(pt.x-leftdesc) <= r      # двигаем левую границу на спуске
                leftdesc = round(Int, pt.x)
            elseif abs(pt.x-rightpump) <= r     #двигаем правую границу на накачке
                rightpump = round(Int, pt.x)
            elseif abs(pt.x-rightdesc) <= r     #двигаем правую границу на спуске
                rightdesc = round(Int, pt.x)
            end
            v.plotbounds.AD = (pump = Bounds(leftpump, rightpump), desc = Bounds(leftdesc, rightdesc))
        end

        return whichplot
    end
end

function ModeRadioButton(v::Globals)   # Радио-кнопки выбора режима работы
    CImGui.SameLine(200)
    CImGui.RadioButton("Изменить границы рабочей зоны", v.mode == 1) && (v.mode == 1 ? v.mode = 0 : v.mode = 1;); CImGui.SameLine()
    rbwidth = CImGui.CalcItemWidth()
    CImGui.RadioButton("Изменить границы АД", v.mode == 2) && (v.mode == 2 ? v.mode = 0 : v.mode = 2;); CImGui.SameLine()
    CImGui.RadioButton("Добавить или ретипизировать метку", v.mode == 3) && (v.mode == 3 ? v.mode = 0 : v.mode = 3;); CImGui.SameLine()
    CImGui.RadioButton("Удалить метку", v.mode == 4) && (v.mode == 4 ? v.mode = 0 : v.mode = 4;)

    return rbwidth
end

function BoundsInput(id::String, pos::Int64, sig, bound::Int64, ispump::Bool) # Поле ввода и отображения значения
    CImGui.SameLine(pos)
    str01 = round(Int, sig[bound])
    str1 = string(str01)
    CImGui.SetNextItemWidth(70)
    CImGui.InputText(id, str1, 4)
    str1 = filter(x -> x!='\0', str1)
    newbound = bound
    if parse(Int64, str1) != str01 
        for i in 2:lastindex(sig)
            if ispump 
                if sig[i] > parse(Int64, str1) && sig[i-1] <= parse(Int64, str1) newbound = i-1; return newbound end
            else 
                if sig[i] <= parse(Int64, str1) && sig[i-1] > parse(Int64, str1) newbound = i; return newbound end
            end
        end
    end

    return newbound
end

function BoundsFields(v::Globals) # Поля ввода и отображения значений границ рабочей зоны и АД
    CImGui.NewLine()

    sig = v.dataforplotting.rawPres.sig

    idadpump = v.plotbounds.AD.pump.ibeg
    isadpump = v.plotbounds.AD.pump.iend
    isaddesc = v.plotbounds.AD.desc.ibeg
    idaddesc = v.plotbounds.AD.desc.iend

    iwendpump = v.plotbounds.workreg.pump.ibeg
    iwbegpump = v.plotbounds.workreg.pump.iend
    iwbegdesc = v.plotbounds.workreg.desc.ibeg
    iwenddesc = v.plotbounds.workreg.desc.iend

    CImGui.SameLine(670)
    CImGui.TextUnformatted("Накачка:")
    isadpump = BoundsInput("##SADpump", 775, sig, isadpump, true)
    idadpump = BoundsInput("##DADpump", 850, sig, idadpump, true)


    CImGui.SameLine(240)
    CImGui.TextUnformatted("Накачка:")
    iwbegpump = BoundsInput("##WBEGpump", 345, sig, iwbegpump, true)
    iwendpump = BoundsInput("##WENDpump", 420, sig, iwendpump, true)

    CImGui.NewLine()
    CImGui.SameLine(670)
    CImGui.TextUnformatted("Спуск:")
    isaddesc = BoundsInput("##SADdesc", 775, sig, isaddesc, false)
    idaddesc = BoundsInput("##DADdesc", 850, sig, idaddesc, false)

    CImGui.SameLine(240)
    CImGui.TextUnformatted("Спуск:")
    iwbegdesc = BoundsInput("##WBEGdesc", 345, sig, iwbegdesc, false)
    iwenddesc = BoundsInput("##WENDdesc", 420, sig, iwenddesc, false)

    v.plotbounds.AD = (pump = Bounds(idadpump, isadpump), desc = Bounds(isaddesc, idaddesc))
    v.plotbounds.workreg = (pump = Bounds(iwendpump, iwbegpump), desc = Bounds(iwbegdesc, iwenddesc))
end

function ChangeTypeCombo(v::Globals) # Выпадающий список типов меток
    CImGui.SameLine(970)
    CImGui.SetNextItemWidth(420)
    types = ["Незначимая","Значимая","Шум"]
    if CImGui.BeginCombo("##type", types[v.typecombo_item])
        for i in 1:lastindex(types)
            if CImGui.Selectable(string(types[i]), types[i] == types[v.typecombo_item]) 
                v.typecombo_item = i 
            end
        end
        CImGui.EndCombo()
    end
end

function Info() # Справка
    CImGui.TextDisabled("Функции")
    if CImGui.IsItemHovered()
        CImGui.BeginTooltip()
        CImGui.PushTextWrapPos(CImGui.GetFontSize() * 100.0)
        CImGui.BulletText("Ctrl + клик ЛКМ - ТОЧКА")
        CImGui.BulletText("Alt + удержание ЛКМ - ОБЛАСТЬ")
        CImGui.BulletText("Alt + клик ЛКМ - УДАЛЕНИЕ ОБЛАСТИ")
        CImGui.BulletText("Удержание ЛКМ - ПЕРЕТАСКИВАНИЕ ГРАНИЦ")
        CImGui.Text("\n!Выполняемые действия с точкой или областью зависят от выбранного режима.")
        CImGui.PopTextWrapPos()
        CImGui.EndTooltip()
    end
end

function ChangePlotsID(v::Globals) # Изменение айдишников каждого графика (чтобы масштабы перестраивались только тогда, котогда нужно)
    v.plotsid.ECG -= 2
    v.plotsid.rawPres -= 1
    v.plotsid.Pres += 1
    v.plotsid.Tone += 2
end

function FigureWindow(v::Globals) # Окно графиков с разметкой
    CImGui.Begin("Разметка")

    if !isempty(v.allfiles)
        Info()
        ModeRadioButton(v)
        BoundsFields(v)
        ChangeTypeCombo(v)
        CImGui.NewLine()
        CImGui.BulletText("Двойной щелчок ЛКМ по графику для восстановления корректного масштаба.")
        SaveRefMarkupButton(v)

        ImPlot.PushColormap(ImPlotColormap_Deep)

        CImGui.PushID(v.plotsid.ECG)
            sig = v.dataforplotting.ECG.sig
            isadpump = v.plotbounds.AD.pump.ibeg; idadpump = v.plotbounds.AD.pump.iend
            wbegpump = v.plotbounds.workreg.pump.ibeg; wendpump = v.plotbounds.workreg.pump.iend
            isaddesc = v.plotbounds.AD.desc.ibeg; idaddesc = v.plotbounds.AD.desc.iend
            wbegdesc = v.plotbounds.workreg.desc.ibeg; wenddesc = v.plotbounds.workreg.desc.iend
            ymin = v.dataforplotting.ECG.limits.Y.Min; ymin = ymin > 0 ? ymin*0.8 : ymin*1.1
            ymax = v.dataforplotting.ECG.limits.Y.Max; ymax = ymax > 0 ? ymax*1.1 : ymax*0.8
            xmin = v.dataforplotting.ECG.limits.X.Min; xmax = v.dataforplotting.ECG.limits.X.Max

            flags = v.mode == 1 || v.mode == 2 ? ImGuiCond_Always : ImGuiCond_Once

            ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
            v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
            C_NULL, C_NULL, C_NULL, C_NULL)
            ImPlot.SetNextPlotLimits(xmin, xmax, ymin, ymax, flags)
            if ImPlot.BeginPlot("ЭКГ", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/5),
                                x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
                PlotLine(1, sig)
                PlotLine(2, [isadpump, isadpump], [ymin, ymax])
                PlotLine(2, [idadpump, idadpump], [ymin, ymax])
                PlotLine(3, [isaddesc, isaddesc], [ymin, ymax])
                PlotLine(3, [idaddesc, idaddesc], [ymin, ymax])
                PlotLine(4, [wbegpump, wbegpump], [ymin, ymax])
                PlotLine(4, [wendpump, wendpump], [ymin, ymax])
                PlotLine(5, [wbegdesc, wbegdesc], [ymin, ymax])
                PlotLine(5, [wenddesc, wenddesc], [ymin, ymax])

                MouseClick(v, "ecg")

                ImPlot.EndPlot()
            end
        CImGui.PopID()

        CImGui.PushID(v.plotsid.rawPres)
            sig = v.dataforplotting.rawPres.sig
            isadpump = v.plotbounds.AD.pump.ibeg; idadpump = v.plotbounds.AD.pump.iend
            wbegpump = v.plotbounds.workreg.pump.ibeg; wendpump = v.plotbounds.workreg.pump.iend
            isaddesc = v.plotbounds.AD.desc.ibeg; idaddesc = v.plotbounds.AD.desc.iend
            wbegdesc = v.plotbounds.workreg.desc.ibeg; wenddesc = v.plotbounds.workreg.desc.iend
            ymin = v.dataforplotting.rawPres.limits.Y.Min; ymax = v.dataforplotting.rawPres.limits.Y.Max
            xmin = v.dataforplotting.rawPres.limits.X.Min; xmax = v.dataforplotting.rawPres.limits.X.Max

            flags = v.mode == 1 || v.mode == 2 || unsafe_load(CImGui.GetIO().KeyAlt) ? ImGuiCond_Always : ImGuiCond_Once

            ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
            v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
            C_NULL, C_NULL, C_NULL, C_NULL)
            ImPlot.SetNextPlotLimits(xmin, xmax, ymin*0.8, ymax*1.1, flags)
            if ImPlot.BeginPlot("Давление", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/5),
                                x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
                PlotLine(1, sig)
                PlotLine(2, [isadpump, isadpump], [ymin*0.8, ymax*1.1], label = "САД и ДАД на накачке")
                PlotLine(2, [idadpump, idadpump], [ymin*0.8, ymax*1.1], label = "САД и ДАД на накачке")
                PlotLine(3, [isaddesc, isaddesc], [ymin*0.8, ymax*1.1], label = "САД и ДАД на спуске")
                PlotLine(3, [idaddesc, idaddesc], [ymin*0.8, ymax*1.1], label = "САД и ДАД на спуске")
                PlotLine(4, [wbegpump, wbegpump], [ymin*0.8, ymax*1.1], label = "Рабочая зона на накачке")
                PlotLine(4, [wendpump, wendpump], [ymin*0.8, ymax*1.1], label = "Рабочая зона на накачке")
                PlotLine(5, [wbegdesc, wbegdesc], [ymin*0.8, ymax*1.1], label = "Рабочая зона на спуске")
                PlotLine(5, [wenddesc, wenddesc], [ymin*0.8, ymax*1.1], label = "Рабочая зона на спуске")


                MouseClick(v, "pres")

                ImPlot.EndPlot()
            end
        CImGui.PopID()

        CImGui.PushID(v.plotsid.Pres)
            sig = v.dataforplotting.Pres.sig
            isadpump = v.plotbounds.AD.pump.ibeg; idadpump = v.plotbounds.AD.pump.iend
            wbegpump = v.plotbounds.workreg.pump.ibeg; wendpump = v.plotbounds.workreg.pump.iend
            isaddesc = v.plotbounds.AD.desc.ibeg; idaddesc = v.plotbounds.AD.desc.iend
            wbegdesc = v.plotbounds.workreg.desc.ibeg; wenddesc = v.plotbounds.workreg.desc.iend
            ymin = v.dataforplotting.Pres.limits.Y.Min; ymax = v.dataforplotting.Pres.limits.Y.Max
            xmin = v.dataforplotting.Pres.limits.X.Min; xmax = v.dataforplotting.Pres.limits.X.Max

            type = v.dataforplotting.Pres.type; ibegs = v.dataforplotting.Pres.ibegs; iends = v.dataforplotting.Pres.iends
            begs0 = map((x,y) -> x == 0 ? y : -1, type, ibegs); begs0 = filter(x -> x!=-1, begs0)
            ends0 = map((x,y) -> x == 0 ? y : -1, type, iends); ends0 = filter(x -> x!=-1, ends0)

            begs1 = map((x,y) -> x == 1 ? y : -1, type, ibegs); begs1 = filter(x -> x!=-1, begs1)
            ends1 = map((x,y) -> x == 1 ? y : -1, type, iends); ends1 = filter(x -> x!=-1, ends1)

            begs2 = map((x,y) -> x == 2 ? y : -1, type, ibegs); begs2 = filter(x -> x!=-1, begs2)
            ends2 = map((x,y) -> x == 2 ? y : -1, type, iends); ends2 = filter(x -> x!=-1, ends2)

            flags = v.mode == 1 || v.mode == 2 || unsafe_load(CImGui.GetIO().KeyAlt) ? ImGuiCond_Always : ImGuiCond_Once

            ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
            v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
            C_NULL, C_NULL, C_NULL, C_NULL)
            ImPlot.SetNextPlotLimits(xmin, xmax, ymin*1.1, ymax*1.1, flags)
            if ImPlot.BeginPlot("Пульсации", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/5),
                                x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
                PlotLine(1, sig)
                PlotLine(2, [isadpump, isadpump], [ymin, ymax]*1.1)
                PlotLine(2, [idadpump, idadpump], [ymin, ymax]*1.1)
                PlotLine(3, [isaddesc, isaddesc], [ymin, ymax]*1.1)
                PlotLine(3, [idaddesc, idaddesc], [ymin, ymax]*1.1)
                PlotLine(4, [wbegpump, wbegpump], [ymin, ymax]*1.1)
                PlotLine(4, [wendpump, wendpump], [ymin, ymax]*1.1)
                PlotLine(5, [wbegdesc, wbegdesc], [ymin, ymax]*1.1)
                PlotLine(5, [wenddesc, wenddesc], [ymin, ymax]*1.1)

                if !isempty(begs1)
                    Scatter(6, begs1, sig[begs1], ImPlotMarker_Circle, 5, "Значимые")
                    Scatter(6, ends1, sig[ends1], ImPlotMarker_Square, 5, "Значимые")
                end

                if !isempty(begs0)
                    Scatter(7, begs0, sig[begs0], ImPlotMarker_Circle, 5, "Незначимые")
                    Scatter(7, ends0, sig[ends0], ImPlotMarker_Square, 5, "Незначимые")
                end

                if !isempty(begs2)
                    Scatter(8, begs2, sig[begs2], ImPlotMarker_Circle, 5, "Шум")
                    Scatter(8, ends2, sig[ends2], ImPlotMarker_Square, 5, "Шум")
                end

                MouseClick(v, "pulse")

                if v.area.whichplot == "pulse"
                    InsideArea(v, "pulse")
                    if v.area.begpos != v.area.endpos
                        PlotLine(9, [v.area.begpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.begpos.y])
                        PlotLine(9, [v.area.begpos.x, v.area.endpos.x], [v.area.endpos.y, v.area.endpos.y])
                        PlotLine(9, [v.area.begpos.x, v.area.begpos.x], [v.area.begpos.y, v.area.endpos.y])
                        PlotLine(9, [v.area.endpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.endpos.y])
                    end
                end

                ImPlot.EndPlot()
            end
        CImGui.PopID()

        CImGui.PushID(v.plotsid.Tone)
            sig = v.dataforplotting.Tone.sig
            isadpump = v.plotbounds.AD.pump.ibeg; idadpump = v.plotbounds.AD.pump.iend
            wbegpump = v.plotbounds.workreg.pump.ibeg; wendpump = v.plotbounds.workreg.pump.iend
            isaddesc = v.plotbounds.AD.desc.ibeg; idaddesc = v.plotbounds.AD.desc.iend
            wbegdesc = v.plotbounds.workreg.desc.ibeg; wenddesc = v.plotbounds.workreg.desc.iend
            ymin = v.dataforplotting.Tone.limits.Y.Min; ymax = v.dataforplotting.Tone.limits.Y.Max
            xmin = v.dataforplotting.Tone.limits.X.Min; xmax = v.dataforplotting.Tone.limits.X.Max

            type = v.dataforplotting.Tone.type; peaks = v.dataforplotting.Tone.ibegs
            peaks0 = map((x,y) -> x == 0 ? y : -1, type, peaks); peaks0 = filter(x -> x!=-1, peaks0)
            peaks1 = map((x,y) -> x == 1 ? y : -1, type, peaks); peaks1 = filter(x -> x!=-1, peaks1)
            peaks2 = map((x,y) -> x == 2 ? y : -1, type, peaks); peaks2 = filter(x -> x!=-1, peaks2)

            flags = v.mode == 1 || v.mode == 2 || unsafe_load(CImGui.GetIO().KeyAlt) ? ImGuiCond_Always : ImGuiCond_Once

            ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
            v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
            C_NULL, C_NULL, C_NULL, C_NULL)
            ImPlot.SetNextPlotLimits(xmin, xmax, ymin*1.1, ymax*1.1, flags)
            if ImPlot.BeginPlot("Тоны", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/5),
                                x_flags = ImPlotAxisFlags_NoGridLines, y_flags = ImPlotAxisFlags_NoGridLines)
                PlotLine(1, sig)
                PlotLine(2, [isadpump, isadpump], [ymin, ymax]*1.1)
                PlotLine(2, [idadpump, idadpump], [ymin, ymax]*1.1)
                PlotLine(3, [isaddesc, isaddesc], [ymin, ymax]*1.1)
                PlotLine(3, [idaddesc, idaddesc], [ymin, ymax]*1.1)
                PlotLine(4, [wbegpump, wbegpump], [ymin, ymax]*1.1)
                PlotLine(4, [wendpump, wendpump], [ymin, ymax]*1.1)
                PlotLine(5, [wbegdesc, wbegdesc], [ymin, ymax]*1.1)
                PlotLine(5, [wenddesc, wenddesc], [ymin, ymax]*1.1)

                if !isempty(peaks1)
                    Scatter(6, peaks1, sig[peaks1], ImPlotMarker_Circle, 5, "Значимые")
                end

                if !isempty(peaks0)
                    Scatter(7, peaks0, sig[peaks0], ImPlotMarker_Circle, 5, "Незначимые")
                end

                if !isempty(peaks2)
                    Scatter(8, peaks2, sig[peaks2], ImPlotMarker_Circle, 5, "Шум")
                end

                if !isempty(v.selected_peaks) 
                    Scatter(9, v.selected_peaks, sig[v.selected_peaks], ImPlotMarker_Circle, 5, false)
                end

                MouseClick(v, "tone")

                if v.area.whichplot == "tone"
                    InsideArea(v, "tone")
                    if v.area.begpos != v.area.endpos
                        PlotLine(10, [v.area.begpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.begpos.y])
                        PlotLine(10, [v.area.begpos.x, v.area.endpos.x], [v.area.endpos.y, v.area.endpos.y])
                        PlotLine(10, [v.area.begpos.x, v.area.begpos.x], [v.area.begpos.y, v.area.endpos.y])
                        PlotLine(10, [v.area.endpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.endpos.y])
                    end
                end

                ImPlot.EndPlot()
            end
        CImGui.PopID()

        ImPlot.PopColormap()
        
        CImGui.End()
    end
end

function ReadBase(v::Globals)
    folder = v.fold
    listoffiles = readdir(folder)
    allbins = map(x -> split(x,".")[end] == "bin" ? split(x,".")[1] : "", listoffiles)
    v.allfiles = filter(x -> !isempty(x), allbins)
    ReadData(v)
end

function BasesTable(v::Globals)
    if !isempty(v.allbases)
        CImGui.NewLine()
        CImGui.BeginChild("##bases_header", ImVec2(CImGui.GetWindowContentRegionWidth(), CImGui.GetTextLineHeightWithSpacing()*1.3))
        CImGui.Columns(1, "##bh")
        CImGui.Separator()
        CImGui.TextColored(ImVec4(0.45, 0.7, 0.80, 1.00), "Имена баз")
        CImGui.Columns(1)
        CImGui.Separator()
        CImGui.EndChild()

        allbases = map(x -> split(x, "/")[end], v.allbases)
        CImGui.BeginChild("##bases_scrollingregion", (0, 100))
        CImGui.Columns(1, "##fn")
        for i in 1:lastindex(allbases)
            CImGui.PushID(i)
            if CImGui.Selectable(allbases[i], i == v.selectedbase)
                v.selectedbase = i
                v.fold = v.allbases[i]
                ReadBase(v)
                GeneratePlotData(v)
                ChangePlotsID(v)
            end
            CImGui.PopID()
            CImGui.NextColumn()
            CImGui.Separator()
        end
        CImGui.Columns(1)
        CImGui.EndChild()

    end
end

function LoadAllBases(v::Globals)
    if !v.isguistarted
        dir = "D:/INCART/Pulse_Data/все базы"
        allfolds = readdir(dir)
        v.allbases = map(x -> "$dir/$x", allfolds)
        v.isguistarted = true
    end
end

function ui(v::Globals) # Главное окно программы (Окно меню + окно графиков с разметкой)
    LoadAllBases(v)
    MenuWindow(v)
    FigureWindow(v)
end

function show_gui() # Main
    state = Globals();
    Renderer.render(
        ()->ui(state),
        width = 1800,
        height = 1600,
        title = "",
        v = Renderer.GR()
    )
end

show_gui();
