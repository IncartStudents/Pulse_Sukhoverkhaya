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

struct Bounds
    ibeg::Int64
    iend::Int64
end

struct Signal
    ECG::Vector{Float64}
    Pres::Vector{Float64}
    Tone::Vector{Float64}
    fs::Int64
    validsegs::Vector{Bounds}
end

struct Mkp
    Pres::Vector{Vector{PresEv}}
    Tone::Vector{Vector{ToneEv}}
end

struct Area
    begpos::ImVec2
    endpos::ImVec2
    whichplot::String
end

mutable struct PlotElem
    sig::Vector{Float64}

    ibegs::Vector{Int64}
    iends::Vector{Int64}
    type::Vector{Int64} # 0 - незначимые (отсутствие метки), 1 - значимые, 2 - шумовые события

    limits::ImPlotLimits
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
    segbounds::Bounds
end

mutable struct PlotLinks
    xmin::Ref
    xmax::Ref
    ymin::Ref
    ymax::Ref
    linkx::Bool
    linky::Bool
end

mutable struct PlotID
    ECG::Int64
    rawPres::Int64
    Pres::Int64
    Tone::Int64
end

mutable struct Globals
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

    function Globals()
        filename = ""
        signal = Signal(Float64[], Float64[], Float64[], 0, Bounds[])
        markup = Mkp(Vector{PresEv}[], Vector{ToneEv}[])
        tup = PlotElem(Float64[], Int64[], Int64[], Int64[], ImPlotLimits(ImPlotRange(0.0,0.0), ImPlotRange(0.0,0.0)))
        plotbounds = PlotBounds(Bounds(0,0), Bounds(0,0), Bounds(0,0))
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

        new(filename, signal, markup, plotbounds, dataforplotting, 
            plotlinks, combo_item, typecombo_item, mode, 
            selected_peaks, area, pt0, plotsid, allfiles, tabledata, fold,
            selecteditem)
    end
end

function GeneratePlotData(v::Globals)
    # референтные границы САД-ДАД, внутри которых считаем статистики
    # adtablefile = "D:/INCART/Pulse_Data/ad result tables/$(v.allfiles[v.combo_item])_table_ad.txt"
    # ad = read_ad(adtablefile)

    ad = map(x -> AD(x[2][2],x[3][2]), v.tabledata)

    if length(ad) < v.selecteditem 
        bsad = 1; bdad = v.signal.validsegs[v.selecteditem].iend - v.signal.validsegs[v.selecteditem].ibeg
    else
        bounds = get_ad_bounds(v.signal.Pres[v.signal.validsegs[v.selecteditem].ibeg:v.signal.validsegs[v.selecteditem].iend], ad[v.selecteditem])
        if bounds.isad == 0 || bounds.idad == 0 bsad = 1; bdad = v.signal.validsegs[v.selecteditem].iend - v.signal.validsegs[v.selecteditem].ibeg
        else bsad = bounds.isad; bdad = bounds.idad end
    end

    # ЭКГ
    ECG = v.signal.ECG[v.signal.validsegs[v.selecteditem].ibeg:v.signal.validsegs[v.selecteditem].iend]

    # тоны
    seg = v.signal.Tone[v.signal.validsegs[v.selecteditem].ibeg:v.signal.validsegs[v.selecteditem].iend]

    smoothtone = my_butter(seg, 2, 60, v.signal.fs, "low") # сглаженный тонов
    ftone = my_butter(smoothtone, 2, 30, v.signal.fs, "high") # фильтрованный тонов
    tone_sig = my_butter(abs.(ftone), 2, 10, v.signal.fs, "low") # огибающая по модулю

    tone_peaks = map(x -> x.pos-v.signal.validsegs[v.selecteditem].ibeg+1, v.markup.Tone[v.selecteditem])

    # пульсации
    seg = v.signal.Pres[v.signal.validsegs[v.selecteditem].ibeg:v.signal.validsegs[v.selecteditem].iend]

    fsig_smooth = my_butter(seg, 2, 10, v.signal.fs, "low") # сглаживание
    pres_sig = my_butter(fsig_smooth, 2, 0.3, v.signal.fs, "high") # устранение постоянной составляющей

    pres_begs = map(x -> x.ibeg-v.signal.validsegs[v.selecteditem].ibeg+1, v.markup.Pres[v.selecteditem])
    pres_ends = map(x -> x.iend-v.signal.validsegs[v.selecteditem].ibeg+1, v.markup.Pres[v.selecteditem])

    # границы рабочей зоны (пик треугольника давления - 10 мм, минимум спуска + 10 мм)
    lvlbeg = maximum(seg) - 30; lvlend = seg[end] >= 30 ? seg[end] : 30
    wbeg = 1; wend = length(seg)
    for i in 2:lastindex(seg) if seg[i] <= lvlbeg && seg[i-1] > lvlbeg wbeg = i
                                elseif seg[i] < lvlend && seg[i-1] >= lvlend wend = i end end

    # расставление типов меток по умолчанию (все, что за границами рабочей зоны - в незначимые, остальные - (пока) в значимые)
    tone_peaks_type = map(x -> x < wbeg || x > wend ? 0 : 1, tone_peaks)
    pres_peaks_type = map(x -> x < wbeg || x > wend ? 0 : 1, pres_begs)

    ad_bounds = Bounds(bsad, bdad)
    workreg = Bounds(wbeg, wend)
    ECGtup = PlotElem(ECG, Int64[], Int64[], Int64[], ImPlotLimits(ImPlotRange(0.0, length(ECG) |> Float64),ImPlotRange(minimum(ECG), maximum(ECG))))
    RawPrestup = PlotElem(seg, Int64[], Int64[], Int64[], ImPlotLimits(ImPlotRange(0.0, length(seg) |> Float64),ImPlotRange(minimum(seg), maximum(seg))))
    Prestup = PlotElem(pres_sig, pres_begs, pres_ends, pres_peaks_type, ImPlotLimits(ImPlotRange(0.0, length(pres_sig) |> Float64),ImPlotRange(minimum(pres_sig), maximum(pres_sig))))
    Tonetup = PlotElem(tone_sig, tone_peaks, tone_peaks, tone_peaks_type, ImPlotLimits(ImPlotRange(0.0, length(tone_sig) |> Float64),ImPlotRange(minimum(tone_sig), maximum(tone_sig))))

    v.dataforplotting = PlotData(ECGtup, RawPrestup, Prestup, Tonetup)
    v.plotbounds = PlotBounds(ad_bounds, workreg, v.signal.validsegs[v.selecteditem])
end

function MeasuresCombo(v::Globals)
    if !isempty(v.allfiles)
        if CImGui.BeginCombo("##meas", string(v.combo_item))
            for i in 1:length(v.markup.Pres)
                if CImGui.Selectable(string(i), i == v.combo_item) 
                    v.combo_item = i 
                    GeneratePlotData(v)
                end
            end
            ChangePlotsID(v)
            CImGui.EndCombo()
        end
    end
end

function ReadData(fname, v::Globals)
    # чтение данных из файла
    signals, fs, _, _ = readbin(fname)

    ECG = signals[1]     # ЭКГ
    Tone = signals.Tone  # пульсации
    Pres = signals.Pres  # давление

    # получение валидных сегментов
    vseg = get_valid_segments(Pres, Tone, 15, -1e7, 30*fs)
    vseg = map(x -> Bounds(x.ibeg, x.iend), vseg)

    v.signal = Signal(ECG, Pres, Tone, fs, vseg)

    # парсинг тестовой разметки
    Pres_mkp = test_markup_parse("alg markup/$(v.allfiles[v.combo_item]).pres")
    Tone_mkp = test_markup_parse("alg markup/$(v.allfiles[v.combo_item]).tone")

    # референтные границы САД-ДАД, внутри которых считаем статистики
    adtablefile = "D:/INCART/Pulse_Data/ad result tables/$(v.allfiles[v.combo_item])_table_ad.txt"
    ad0 = read_ad(adtablefile)

    n = length(Pres_mkp)
    ad = fill(AD(0,0), maximum([n,length(ad0)]))
    ad[1:length(ad0)] = ad0

    v.markup = Mkp(Pres_mkp, Tone_mkp)

    v.tabledata = map((y,z) -> ("Номер записи" => y, "САД реф." => z.SAD, "ДАД реф." => z.DAD), 
                                                                        range(1,n), ad)

    GeneratePlotData(v)
end

function LoadButtons(v::Globals)
    # Выбирается бинарь, из той же папки подтягивается соответсвующий hdr, 
    # из папки results (пока что) подтягивается соответствующий файл разметки
    if CImGui.Button("Загрузить базу")
        folder = open_dialog_native("Выберите папку", action = GtkFileChooserAction.SELECT_FOLDER)
        v.fold = folder
        listoffiles = readdir(folder)
        allbins = map(x -> split(x,".")[end] == "bin" ? split(x,".")[1] : "", listoffiles)
        v.allfiles = filter(x -> !isempty(x), allbins)
        ReadData(folder*"/"*v.allfiles[v.combo_item], v)
        # v.allfiles[v.combo_item] = split(split(fname, ".")[end-1], "\\")[end]
    end
end

function SaveRefMarkup(filename::String, ext::String, markup::PlotElem)
    open(filename*ext, "w") do io

        if ext == ".pres" write(io, "beg   end   type")
        elseif ext == ".tone" write(io, "pos   type") end

        for j in 1:lastindex(markup.ibegs)
            if ext == ".pres"
                write(io, "\n$(markup.ibegs[j])   $(markup.iends[j])   $(markup.type[j])")
            elseif ext == ".tone"
                write(io, "\n$(markup.ibegs[j])   $(markup.type[j])")
            end
        end

    end
end

function SaveRefMarkup(filename::String, ext::String, markup::PlotBounds)
    open(filename*ext, "w") do io
        write(io, "var   beg   end")
        write(io, "\nSEG   $(markup.segbounds.ibeg)   $(markup.segbounds.iend)")
        write(io, "\nAD   $(markup.AD.ibeg)   $(markup.AD.iend)")
        write(io, "\nWZ   $(markup.workreg.ibeg)   $(markup.workreg.iend)")
    end
end

function SaveRefMarkupButton(v::Globals)
    if !isempty(v.allfiles)
        CImGui.SameLine()
        if CImGui.Button("Сохранить разметку")
            dirname = "ref markup/$(v.allfiles[v.combo_item])/measure $(v.combo_item)"
            try readdir(dirname)
            catch e mkdir(dirname) end

            filename = "$dirname/$(v.combo_item)"
            ext = [".pres", ".tone", ".bounds"]
            markup = [v.dataforplotting.Pres, v.dataforplotting.Tone, v.plotbounds]

            for i in 1:lastindex(ext) SaveRefMarkup(filename, ext[i], markup[i]) end
        end
    end
end

function FilenamesCombo(v::Globals)
    CImGui.SameLine(200)
    if !isempty(v.allfiles)
        CImGui.SetNextItemWidth(500)
        if CImGui.BeginCombo("Имена файлов базы", v.allfiles[v.combo_item])
            for i in 1:length(v.allfiles)
                if CImGui.Selectable(v.allfiles[i], i == v.combo_item) 
                    v.combo_item = i 
                    ReadData(v.fold*"/"*v.allfiles[v.combo_item], v)
                    GeneratePlotData(v)
                end
            end
            ChangePlotsID(v)
            CImGui.EndCombo()
        end
    end
end

function MeasuresTable(v::Globals)
    if !isempty(v.allfiles)
        CImGui.NewLine()
        names = map(x -> x[1], v.tabledata[1])
        col = length(names)
        CImGui.Columns(col, "Записи и измерения")
        CImGui.Separator()
        for i in names
            CImGui.Text(i)
            CImGui.NextColumn()
        end
        CImGui.Separator()
        for i in 1:lastindex(v.tabledata)
            for j in v.tabledata[i]
                CImGui.PushID(i)
                if CImGui.Selectable(string(j[2]), i == v.selecteditem, CImGui.ImGuiSelectableFlags_SpanAllColumns)
                    v.selecteditem = i 
                    GeneratePlotData(v)
                end
                CImGui.PopID()
                CImGui.NextColumn()
            end
            CImGui.Separator()
        end
        CImGui.Columns(1)
    end
end

function MenuWindow(v::Globals)
    # CImGui.SetNextWindowPos(ImVec2(0,0))
    # CImGui.SetNextWindowSize(ImVec2(s.w, s.h/2))
    CImGui.Begin("Меню")
        LoadButtons(v)
        FilenamesCombo(v)
        # MeasuresCombo(v) # Комбо-бокс с выбором измерения
        SaveRefMarkupButton(v)
        MeasuresTable(v)
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

function InsideArea(v::Globals, whichplot)
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

function MouseClick(v::Globals, ymax, whichplot)
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
        mp = round(pt.x) |> Int64
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

        return whichplot
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

        return whichplot
    end

    if ImPlot.IsPlotHovered() && CImGui.IsMouseDown(0)  # захват (для передвижения границ) (курсор наведен на любой график, левая клавиша мыши зажата)
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

        return whichplot
    end
end

function ModeRadioButton(v::Globals)
    CImGui.SameLine(200)
    CImGui.RadioButton("Изменить границы рабочей зоны", v.mode == 1) && (v.mode == 1 ? v.mode = 0 : v.mode = 1;); CImGui.SameLine()
    rbwidth = CImGui.CalcItemWidth()
    CImGui.RadioButton("Изменить границы АД", v.mode == 2) && (v.mode == 2 ? v.mode = 0 : v.mode = 2;); CImGui.SameLine()
    CImGui.RadioButton("Добавить или ретипизировать метку", v.mode == 3) && (v.mode == 3 ? v.mode = 0 : v.mode = 3;); CImGui.SameLine()
    CImGui.RadioButton("Удалить метку", v.mode == 4) && (v.mode == 4 ? v.mode = 0 : v.mode = 4;)

    return rbwidth
end

function BoundsInput(id::String, pos::Int64, sig, bound::Int64)
    CImGui.SameLine(pos)
    str01 = round(sig[bound]) |> Int64
    str1 = string(str01)
    CImGui.SetNextItemWidth(70)
    CImGui.InputText(id, str1, 4)
    str1 = filter(x -> x!='\0', str1)
    newbound = bound
    if parse(Int64, str1) != str01 
        for i in 2:lastindex(sig)
            if sig[i] <= parse(Int64, str1) && sig[i-1] > parse(Int64, str1) newbound = i; break end
        end
    end

    return newbound
end

function BoundsFields(v::Globals)
    CImGui.NewLine()

    sig = v.dataforplotting.rawPres.sig

    isad = v.plotbounds.AD.ibeg
    idad = v.plotbounds.AD.iend

    iwbeg = v.plotbounds.workreg.ibeg
    iwend = v.plotbounds.workreg.iend

    isad = BoundsInput("##SAD", 670, sig, isad)
    idad = BoundsInput("##DAD", 745, sig, idad)

    iwbeg = BoundsInput("##WBEG", 240, sig, iwbeg)
    iwend = BoundsInput("##WEND", 315, sig, iwend)

    v.plotbounds.AD = Bounds(isad, idad)
    v.plotbounds.workreg = Bounds(iwbeg, iwend)
end

function ChangeTypeCombo(v::Globals)
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

function Info()
    CImGui.TextDisabled("Справка")
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

function ChangePlotsID(v::Globals)
    v.plotsid.ECG -= 2
    v.plotsid.rawPres -= 1
    v.plotsid.Pres += 1
    v.plotsid.Tone += 2
end

function FigureWindow(v::Globals)

    CImGui.Begin("Разметка")

    if !isempty(v.allfiles)
        Info()
        ModeRadioButton(v)
        BoundsFields(v)
        ChangeTypeCombo(v)

    ImPlot.PushColormap(ImPlotColormap_Deep)

    CImGui.PushID(v.plotsid.ECG)
        sig = v.dataforplotting.ECG.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        ymin = v.dataforplotting.ECG.limits.Y.Min; ymax = v.dataforplotting.ECG.limits.Y.Max
        xmin = v.dataforplotting.ECG.limits.X.Min; xmax = v.dataforplotting.ECG.limits.X.Max

        flags = v.mode == 1 || v.mode == 2 ? ImGuiCond_Always : ImGuiCond_Once

        ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
        v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
        C_NULL, C_NULL, C_NULL, C_NULL)
        ImPlot.SetNextPlotLimits(xmin, xmax, ymin, ymax, flags)
        if ImPlot.BeginPlot("ЭКГ", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin, ymax]*1.1)
            PlotLine(2, [idad, idad], [ymin, ymax]*1.1)

            MouseClick(v, ymax, "ecg")

            ImPlot.EndPlot()
        end
    CImGui.PopID()

    CImGui.PushID(v.plotsid.rawPres)
        sig = v.dataforplotting.rawPres.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        wbeg = v.plotbounds.workreg.ibeg; wend = v.plotbounds.workreg.iend
        ymin = v.dataforplotting.rawPres.limits.Y.Min; ymax = v.dataforplotting.rawPres.limits.Y.Max
        xmin = v.dataforplotting.rawPres.limits.X.Min; xmax = v.dataforplotting.rawPres.limits.X.Max

        flags = v.mode == 1 || v.mode == 2 || unsafe_load(CImGui.GetIO().KeyAlt) ? ImGuiCond_Always : ImGuiCond_Once

        ImPlot.LinkNextPlotLimits(v.plotlinks.linkx ? v.plotlinks.xmin : C_NULL, v.plotlinks.linkx ? v.plotlinks.xmax : C_NULL,
        v.plotlinks.linky ? v.plotlinks.ymin : C_NULL, v.plotlinks.linky ? v.plotlinks.ymax : C_NULL,
        C_NULL, C_NULL, C_NULL, C_NULL)
        ImPlot.SetNextPlotLimits(xmin, xmax, ymin*0.8, ymax*1.1, flags)
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

    CImGui.PushID(v.plotsid.Pres)
        sig = v.dataforplotting.Pres.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        wbeg = v.plotbounds.workreg.ibeg; wend = v.plotbounds.workreg.iend
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
        if ImPlot.BeginPlot("Пульсации", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoDecorations, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin, ymax]*1.1)
            PlotLine(2, [idad, idad], [ymin, ymax]*1.1)
            PlotLine(3, [wbeg, wbeg], [ymin, ymax]*1.1)
            PlotLine(3, [wend, wend], [ymin, ymax]*1.1)

            if !isempty(begs1)
                Scatter(4, begs1, sig[begs1], ImPlotMarker_Circle, 5, "Значимые")
                Scatter(4, ends1, sig[ends1], ImPlotMarker_Square, 5, "Значимые")
            end

            if !isempty(begs0)
                Scatter(5, begs0, sig[begs0], ImPlotMarker_Circle, 5, "Незначимые")
                Scatter(5, ends0, sig[ends0], ImPlotMarker_Square, 5, "Незначимые")
            end

            if !isempty(begs2)
                Scatter(6, begs2, sig[begs2], ImPlotMarker_Circle, 5, "Шум")
                Scatter(6, ends2, sig[ends2], ImPlotMarker_Square, 5, "Шум")
            end

            whichplot = MouseClick(v, ymax, "pulse")

            if v.area.whichplot == "pulse"
                InsideArea(v, "pulse")
                if v.area.begpos != v.area.endpos
                    PlotLine(8, [v.area.begpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.begpos.y])
                    PlotLine(8, [v.area.begpos.x, v.area.endpos.x], [v.area.endpos.y, v.area.endpos.y])
                    PlotLine(8, [v.area.begpos.x, v.area.begpos.x], [v.area.begpos.y, v.area.endpos.y])
                    PlotLine(8, [v.area.endpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.endpos.y])
                end
            end

            ImPlot.EndPlot()
        end
    CImGui.PopID()

    CImGui.PushID(v.plotsid.Tone)
        sig = v.dataforplotting.Tone.sig
        isad = v.plotbounds.AD.ibeg; idad = v.plotbounds.AD.iend
        wbeg = v.plotbounds.workreg.ibeg; wend = v.plotbounds.workreg.iend
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
        if ImPlot.BeginPlot("Тоны", C_NULL, C_NULL, ImVec2(CImGui.GetWindowContentRegionMax().x, CImGui.GetWindowContentRegionMax().y/4.5),
                            x_flags = ImPlotAxisFlags_NoGridLines, y_flags = ImPlotAxisFlags_NoGridLines)
            PlotLine(1, sig)
            PlotLine(2, [isad, isad], [ymin, ymax]*1.1)
            PlotLine(2, [idad, idad], [ymin, ymax]*1.1)
            PlotLine(3, [wbeg, wbeg], [ymin, ymax]*1.1)
            PlotLine(3, [wend, wend], [ymin, ymax]*1.1)

            if !isempty(peaks1)
                Scatter(4, peaks1, sig[peaks1], ImPlotMarker_Circle, 5, "Значимые")
            end

            if !isempty(peaks0)
                Scatter(5, peaks0, sig[peaks0], ImPlotMarker_Circle, 5, "Незначимые")
            end

            if !isempty(peaks2)
                Scatter(6, peaks2, sig[peaks2], ImPlotMarker_Circle, 5, "Шум")
            end

            if !isempty(v.selected_peaks) 
                Scatter(7, v.selected_peaks, sig[v.selected_peaks], ImPlotMarker_Circle, 5, false)
            end

            whichplot = MouseClick(v, ymax, "tone")

            if v.area.whichplot == "tone"
                InsideArea(v, "tone")
                if v.area.begpos != v.area.endpos
                    PlotLine(8, [v.area.begpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.begpos.y])
                    PlotLine(8, [v.area.begpos.x, v.area.endpos.x], [v.area.endpos.y, v.area.endpos.y])
                    PlotLine(8, [v.area.begpos.x, v.area.begpos.x], [v.area.begpos.y, v.area.endpos.y])
                    PlotLine(8, [v.area.endpos.x, v.area.endpos.x], [v.area.begpos.y, v.area.endpos.y])
                end
            end

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
        width = 1800,
        height = 1600,
        title = "",
        v = Renderer.GR()
    )
end

show_gui();
