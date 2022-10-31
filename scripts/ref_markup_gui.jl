using CImGui
using CImGui: ImVec2
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

# include(joinpath(pathof(ImPlot), "..", "..", "demo", "implot_demo.jl"))
# show_demo()

struct Signal
    ECG::Vector{Float64}
    Pres::Vector{Float64}
    Tone::Vector{Float64}
    fs::Int64
    validsegs::Vector{NamedTuple{(:ibeg, :iend), Tuple{Int64, Int64}}}
end

struct Mkp
    Pres::Vector{Vector{PresEv}}
    Tone::Vector{Vector{ToneEv}}
end

struct Bounds
    ibeg::Int64
    iend::Int64
end

mutable struct PlotElem
    sig::Vector{Float64}

    ibegs::Vector{Int64}
    iends::Vector{Int64}

    ylims::NamedTuple{(:min, :max), Tuple{Float64, Float64}}

    length::Int64
end

mutable struct PlotData
    ECG::PlotElem
    rawPres::PlotElem
    Pres::PlotElem
    Tone::PlotElem
end

mutable struct PlotBounds
    AD::Bounds
    workreg::Bounds
end

mutable struct PlotID
    ecg::Int64
    raw_pres::Int64
    pres::Int64
    tone::Int64
end

mutable struct PlotLinks
    xmin::Ref
    xmax::Ref
    ymin::Ref
    ymax::Ref
    linkx::Bool
    linky::Bool
end

mutable struct Globals
    filename::String
    signal::Signal
    markup::Mkp
    plotbounds::PlotBounds
    dataforplotting::PlotData
    plotid::PlotID
    plotlinks::PlotLinks
    combo_item::Int64
    mode::Int64
    selected_peaks::Vector{Int64}

    function Globals()
        filename = ""
        signal = Signal(Float64[], Float64[], Float64[], 0, NamedTuple{(:ibeg, :iend), Tuple{Int64, Int64}}[])
        markup = Mkp(Vector{PresEv}[], Vector{ToneEv}[])
        tup = PlotElem(Float64[], Int64[], Int64[], (min = 0.0, max = 0.0), 0)
        plotbounds = PlotBounds(Bounds(0,0), Bounds(0,0))
        dataforplotting = PlotData(tup, tup, tup, tup)
        plotid = PlotID(-2, -1, 1, 2)
        plotlinks = PlotLinks(Ref(0.0), Ref(1.0), Ref(0.0), Ref(1.0), true, false)
        combo_item = 1
        mode = 0
        selected_peaks = Int64[]

        new(filename, signal, markup, plotbounds, dataforplotting, plotid, plotlinks, combo_item, mode, selected_peaks)
    end
end

function GeneratePlotData(v::Globals)
    # референтные границы САД-ДАД, внутри которых считаем статистики
    adtablefile = "D:/INCART/Pulse_Data/ad result tables/$(v.filename)_table_ad.txt"
    ad = read_ad(adtablefile)

    if length(ad) < v.combo_item 
        bsad = 0; bdad = v.signal.validsegs[v.combo_item].iend - v.signal.validsegs[v.combo_item].ibeg
    else
        bounds = get_ad_bounds(v.signal.Pres[v.signal.validsegs[v.combo_item].ibeg:v.signal.validsegs[v.combo_item].iend], ad[v.combo_item])
        if bounds.isad == 0 || bounds.idad == 0 bsad = 0; bdad = v.signal.validsegs[v.combo_item].iend - v.signal.validsegs[v.combo_item].ibeg
        else bsad = bounds.isad; bdad = bounds.idad end
    end

    # ЭКГ
    ECG = v.signal.ECG[v.signal.validsegs[v.combo_item].ibeg:v.signal.validsegs[v.combo_item].iend]

    # тоны
    seg = v.signal.Tone[v.signal.validsegs[v.combo_item].ibeg:v.signal.validsegs[v.combo_item].iend]

    smoothtone = my_butter(seg, 2, 60, v.signal.fs, "low") # сглаженный тонов
    ftone = my_butter(smoothtone, 2, 30, v.signal.fs, "high") # фильтрованный тонов
    tone_sig = my_butter(abs.(ftone), 2, 10, v.signal.fs, "low") # огибающая по модулю

    tone_peaks = map(x -> x.pos-v.signal.validsegs[v.combo_item].ibeg+1, v.markup.Tone[v.combo_item])

    # пульсации
    seg = v.signal.Pres[v.signal.validsegs[v.combo_item].ibeg:v.signal.validsegs[v.combo_item].iend]

    fsig_smooth = my_butter(seg, 2, 10, v.signal.fs, "low") # сглаживание
    pres_sig = my_butter(fsig_smooth, 2, 0.3, v.signal.fs, "high") # устранение постоянной составляющей

    pres_begs = map(x -> x.ibeg-v.signal.validsegs[v.combo_item].ibeg+1, v.markup.Pres[v.combo_item])
    pres_ends = map(x -> x.iend-v.signal.validsegs[v.combo_item].ibeg+1, v.markup.Pres[v.combo_item])

    # границы рабочей зоны (пик треугольника давления - 10 мм, минимум спуска + 10 мм)
    lvlbeg = maximum(seg) - 10; lvlend = seg[end] + 10
    wbeg = 0; wend = 0
    for i in 2:lastindex(seg) if seg[i] <= lvlbeg && seg[i-1] > lvlbeg wbeg = i
                                elseif seg[i] < lvlend && seg[i-1] >= lvlend wend = i end end

    ad_bounds = Bounds(bsad, bdad)
    workreg = Bounds(wbeg, wend)
    ECGtup = PlotElem(ECG, Int64[], Int64[], (min = minimum(ECG), max = maximum(ECG)), length(ECG))
    RawPrestup = PlotElem(seg, Int64[], Int64[], (min = minimum(seg), max = maximum(seg)), length(seg))
    Prestup = PlotElem(pres_sig, pres_begs, pres_ends, (min = minimum(pres_sig), max = maximum(pres_sig)), length(pres_sig))
    Tonetup = PlotElem(tone_sig, tone_peaks, tone_peaks, (min = minimum(tone_sig), max = maximum(tone_sig)), length(tone_sig))

    v.dataforplotting = PlotData(ECGtup, RawPrestup, Prestup, Tonetup)
    v.plotbounds = PlotBounds(ad_bounds, workreg)
end

function ChangePlotsID(v::Globals)
    v.plotid.ecg -= 2
    v.plotid.raw_pres -= 2
    v.plotid.pres += 2
    v.plotid.tone += 2
end

function MeasuresCombo(v::Globals)
    if !isempty(v.filename)
        if CImGui.BeginCombo("##meas", string(v.combo_item))
            for i in 1:length(v.markup.Pres)
                if CImGui.Selectable(string(i), i == v.combo_item) 
                    v.combo_item = i 
                    GeneratePlotData(v)
                    ChangePlotsID(v)
                end
            end
            CImGui.EndCombo()
        end
    end
end

function LoadButtons(v::Globals)
    # Выбирается бинарь, из той же папки подтягивается соответсвующий hdr, 
    # из папки results (пока что) подтягивается соответствующий файл разметки
    if CImGui.Button("Загрузить запись")
        fname = open_dialog_native("Выберите bin-file");

        # чтение данных из файла
        signals, fs, _, _ = readbin(fname)

        ECG = signals.LR     # ЭКГ
        Tone = signals.Tone  # пульсации
        Pres = signals.Pres  # давление

        # получение валидных сегментов
        vseg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)

        v.signal = Signal(ECG, Pres, Tone, fs, vseg)

        # парсинг тестовой разметки
        v.filename = split(split(fname, ".")[end-1], "\\")[end]
        Pres_mkp = test_markup_parse("alg markup/$(v.filename).pres")
        Tone_mkp = test_markup_parse("alg markup/$(v.filename).tone")

        v.markup = Mkp(Pres_mkp, Tone_mkp)

        GeneratePlotData(v)
    end
end

function MenuWindow(v::Globals)
    # CImGui.SetNextWindowPos(ImVec2(0,0))
    # CImGui.SetNextWindowSize(ImVec2(s.w, s.h/2))
    CImGui.Begin("Меню")
        LoadButtons(v)
        MeasuresCombo(v) # Комбо-бокс с выбором измерения
    CImGui.End()
end

function PlotLine(id, args...; label = false)
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

function Scatter(id, x, y, markerstyle, markersize, label)
    CImGui.PushID(id)
    ImPlot.SetNextMarkerStyle(markerstyle, markersize)
    if label == false ImPlot.PlotScatter(x, y)
    else ImPlot.PlotScatter(x, y, label_id = label) end
    CImGui.PopID()
end

function isin(point::ImVec2, searchbox::NamedTuple{(:xmin, :xmax, :ymin, :ymax), NTuple{4, Float64}})
    return (point.x > searchbox.xmin && point.x < searchbox.xmax && point.y > searchbox.ymin && point.y < searchbox.ymax)
end

function MouseClick(v::Globals, ymax, whichplot)
    if ImPlot.IsPlotHovered() && CImGui.IsMouseClicked(0) && unsafe_load(CImGui.GetIO().KeyCtrl)  # Для выбора (перед удалением) или добавления точки
        pt = ImPlot.GetPlotMousePos()
        allpeaks = v.dataforplotting.Tone.ibegs
        sig = v.dataforplotting.Tone.sig
        if v.mode == 3     # добавить точку
            # чтобы реже случались промахи, ищем в окрестностях позиции курсора во время клика максимум сигнала
            mp = round(pt.x) |> Int64
            sr = 100
            newpeak = argmax(sig[mp-sr:mp+sr]) + (mp-sr)
            push!(allpeaks, newpeak); sort(allpeaks)
            v.dataforplotting.Tone.ibegs = allpeaks
            v.dataforplotting.Tone.iends = allpeaks
        elseif v.mode == 4 # выбрать перед удалением
            dscale = 0.02*ymax
            searchbox = (xmin = pt.x-dscale, xmax = pt.x+dscale, 
                        ymin = pt.y-dscale, ymax = pt.y+dscale)
            for i in allpeaks
                if isin(ImVec2(i, sig[i]), searchbox) 
                    if !isempty(findall(x -> x==i, v.selected_peaks)) # если пик уже был выбран - отменяем выбор
                        v.selected_peaks = filter(x -> x!=i, v.selected_peaks)
                    else                                             # в противном случае выбираем пик
                        push!(v.selected_peaks, i) 
                    end
                end
            end
        end
    end
    if ImPlot.IsPlotHovered() && CImGui.IsMouseDown(0)  # Для передвижения границ (курсор наведен на любой график, левая клавиша мыши зажата)
        pt = ImPlot.GetPlotMousePos()
        if v.mode == 1 # двигаем границы рабочей зоны
            left = v.plotbounds.workreg.ibeg; right = v.plotbounds.workreg.iend
            if abs(pt.x-left) <= 500.0 # двигаем левую границу
                newleft = round(pt.x) |> Int64
                v.plotbounds.workreg = Bounds(newleft, right)
            elseif abs(pt.x-right) <= 500.0 #двигаем правую границу
                newright = round(pt.x) |> Int64
                v.plotbounds.workreg = Bounds(left, newright)
            end
        elseif v.mode == 2 # двигаем границы реф АД
            left = v.plotbounds.AD.ibeg; right = v.plotbounds.AD.iend
            if abs(pt.x-left) <= 500.0 # двигаем левую границу
                newleft = round(pt.x) |> Int64
                v.plotbounds.AD = Bounds(newleft, right)
            elseif abs(pt.x-right) <= 500.0 #двигаем правую границу
                newright = round(pt.x) |> Int64
                v.plotbounds.AD = Bounds(left, newright)
            end
        end
    end
end

function ModeRadioButton(v::Globals)
    CImGui.SameLine()
    CImGui.RadioButton("Изменить границы рабочей зоны", v.mode == 1) && (v.mode == 1 ? v.mode = 0 : v.mode = 1;); CImGui.SameLine()
    rbwidth = CImGui.CalcItemWidth()
    CImGui.RadioButton("Изменить границы АД", v.mode == 2) && (v.mode == 2 ? v.mode = 0 : v.mode = 2;); CImGui.SameLine()
    CImGui.RadioButton("Добавить метку", v.mode == 3) && (v.mode == 3 ? v.mode = 0 : v.mode = 3;); CImGui.SameLine()
    CImGui.RadioButton("Выбрать метку", v.mode == 4) && (v.mode == 4 ? (v.mode = 0; v.selected_peaks = Int64[]) : v.mode = 4;);

    return rbwidth
end

function BoundsFields(v::Globals)
    CImGui.NewLine(); CImGui.SameLine(400)
    sig = v.dataforplotting.rawPres.sig
    isad = v.plotbounds.AD.ibeg
    idad = v.plotbounds.AD.iend

    SAD = sig[isad]
    DAD = sig[idad]
    
    str01 = round(SAD) |> Int64
    str1 = string(str01)
    CImGui.SetNextItemWidth(70)
    CImGui.InputText("##SAD", str1, length(str1)+1)
    str1 = filter(x -> x!='\0', str1)
    if parse(Int64, str1) != str01 
        for i in 2:lastindex(sig)
            if sig[i] <= parse(Int64, str1) && sig[i-1] > parse(Int64, str1) isad = i; break end
        end
        v.plotbounds.AD = Bounds(isad, idad)
    end

    CImGui.SameLine()
    str2 = string(round(DAD) |> Int64)
    CImGui.SetNextItemWidth(70)
    CImGui.InputText("##DAD", str2, length(str2)+1)
end

function DeleteMark(v::Globals)
    if v.mode == 4
        CImGui.SameLine()
        if CImGui.SmallButton("удалить")
            allpeaks = v.dataforplotting.Tone.ibegs
            for i in v.selected_peaks allpeaks = filter(x -> x!=i, allpeaks) end
            v.dataforplotting.Tone.ibegs = allpeaks
            v.dataforplotting.Tone.iends = allpeaks
            v.selected_peaks = Int64[]
        end
    end
end

function FigureWindow(v::Globals)

    CImGui.Begin("Разметка")

    if !isempty(v.filename)
        if CImGui.Button("Вернуть исходный масштаб")
            ChangePlotsID(v)
        end
        ModeRadioButton(v)
        DeleteMark(v)
        BoundsFields(v)

    ImPlot.PushColormap(ImPlotColormap_Deep)

    CImGui.PushID(v.plotid.ecg)
        sig = v.dataforplotting.ECG.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        ymin = v.dataforplotting.ECG.ylims.min; ymax = v.dataforplotting.ECG.ylims.max
        xmax = v.dataforplotting.ECG.length

        flags = v.mode == 1 || v.mode == 2 ? ImGuiCond_Always : ImGuiCond_Once

        ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
        v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
        C_NULL, C_NULL, C_NULL, C_NULL)
        ImPlot.SetNextPlotLimits(0, xmax, ymin, ymax, flags)
        if ImPlot.BeginPlot("ЭКГ", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin, ymax]*1.1)
            PlotLine(2, [idad, idad], [ymin, ymax]*1.1)

            MouseClick(v, ymax, "ecg")

            ImPlot.EndPlot()
        end
    CImGui.PopID()

    CImGui.PushID(v.plotid.raw_pres)
        sig = v.dataforplotting.rawPres.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        wbeg = v.plotbounds.workreg.ibeg; wend = v.plotbounds.workreg.iend
        ymin = v.dataforplotting.rawPres.ylims.min; ymax = v.dataforplotting.rawPres.ylims.max
        xmax = v.dataforplotting.rawPres.length

        flags = v.mode == 1 || v.mode == 2 ? ImGuiCond_Always : ImGuiCond_Once

        ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
        v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
        C_NULL, C_NULL, C_NULL, C_NULL)
        ImPlot.SetNextPlotLimits(0, xmax, ymin*0.8, ymax*1.1, flags)
        if ImPlot.BeginPlot("Давление", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin*0.8, ymax*1.1], label = "САД и ДАД")
            PlotLine(2, [idad, idad], [ymin*0.8, ymax*1.1], label = "САД и ДАД")
            PlotLine(3, [wbeg, wbeg], [ymin*0.8, ymax*1.1], label = "Рабочая зона")
            PlotLine(3, [wend, wend], [ymin*0.8, ymax*1.1], label = "Рабочая зона")

            MouseClick(v, ymax, "pres")

            ImPlot.EndPlot()
        end
    CImGui.PopID()

    CImGui.PushID(v.plotid.pres)
        sig = v.dataforplotting.Pres.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        wbeg = v.plotbounds.workreg.ibeg; wend = v.plotbounds.workreg.iend
        ymin = v.dataforplotting.Pres.ylims.min; ymax = v.dataforplotting.Pres.ylims.max
        xmax = v.dataforplotting.Pres.length
        begs = v.dataforplotting.Pres.ibegs; ends = v.dataforplotting.Pres.iends

        flags = v.mode == 1 || v.mode == 2 ? ImGuiCond_Always : ImGuiCond_Once

        ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
        v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
        C_NULL, C_NULL, C_NULL, C_NULL)
        ImPlot.SetNextPlotLimits(0, xmax, ymin*1.1, ymax*1.1, flags)
        if ImPlot.BeginPlot("Пульсации", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin, ymax]*1.1)
            PlotLine(2, [idad, idad], [ymin, ymax]*1.1)
            PlotLine(3, [wbeg, wbeg], [ymin, ymax]*1.1)
            PlotLine(3, [wend, wend], [ymin, ymax]*1.1)
            Scatter(4, begs, sig[begs], ImPlotMarker_Circle, 5, "begs")
            Scatter(5, ends, sig[ends], ImPlotMarker_Circle, 5, "ends")

            MouseClick(v, ymax, "pulse")

            ImPlot.EndPlot()
        end
    CImGui.PopID()

    CImGui.PushID(v.plotid.tone)
        sig = v.dataforplotting.Tone.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        wbeg = v.plotbounds.workreg.ibeg; wend = v.plotbounds.workreg.iend
        ymin = v.dataforplotting.Tone.ylims.min; ymax = v.dataforplotting.Tone.ylims.max
        xmax = v.dataforplotting.Tone.length; peaks = v.dataforplotting.Tone.ibegs

        flags = v.mode == 1 || v.mode == 2 ? ImGuiCond_Always : ImGuiCond_Once

        ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
        v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
        C_NULL, C_NULL, C_NULL, C_NULL)
        ImPlot.SetNextPlotLimits(0, xmax, ymin*1.1, ymax*1.1, flags)
        if ImPlot.BeginPlot("Тоны", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin, ymax]*1.1)
            PlotLine(2, [idad, idad], [ymin, ymax]*1.1)
            PlotLine(3, [wbeg, wbeg], [ymin, ymax]*1.1)
            PlotLine(3, [wend, wend], [ymin, ymax]*1.1)
            Scatter(4, peaks, sig[peaks], ImPlotMarker_Circle, 5, false)

            if !isempty(v.selected_peaks) 
                Scatter(5, v.selected_peaks, sig[v.selected_peaks], ImPlotMarker_Circle, 5, false)
            end

            MouseClick(v, ymax, "tone")

            ImPlot.EndPlot()
        end
    CImGui.PopID()

    ImPlot.PopColormap()
    
    CImGui.End()
    end
end

function ui(v::Globals)
    MenuWindow(v)
    FigureWindow(v)
end

function show_gui() # Main
    state = Globals();
    Renderer.render(
        ()->ui(state),
        width = 1700,
        height = 1600,
        title = "",
        v = Renderer.GR()
    )
end

show_gui();