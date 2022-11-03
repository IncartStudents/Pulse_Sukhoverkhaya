using CImGui
using CImGui: ImVec2
using ImPlot
using CSV
using DataFrames
using Gtk
using FileIO
using Images
using Plots
using ImPlot.LibCImGui: ImGuiCond_Once, ImPlotAxisFlags_NoGridLines

include("../src/Renderer.jl")
using .Renderer

include("../src/help_func.jl")
include("../src/readfiles.jl") 
include("../src/my_filt.jl")  

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

struct PlotElem
    sig::Vector{Float64}
    peaks::Vector{Any}
    iSAD::Int64
    iDAD::Int64
    ylims::NamedTuple{(:min, :max), Tuple{Float64, Float64}}
    length::Int64
end

mutable struct PlotID
    pres::Int64
    tone::Int64
end

mutable struct Globals

end

function GeneratePlotData(v::Globals)

    if v.tab_item == 1 sigtype = "pres" elseif v.tab_item == 2 sigtype = "tone" end

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

    if sigtype == "tone"
        seg = v.signal.Tone[v.signal.validsegs[v.combo_item].ibeg:v.signal.validsegs[v.combo_item].iend]

        smoothtone = my_butter(seg, 2, 60, v.signal.fs, "low") # сглаженный тонов
        ftone = my_butter(smoothtone, 2, 30, v.signal.fs, "high") # фильтрованный тонов
        sig = my_butter(abs.(ftone), 2, 10, v.signal.fs, "low") # огибающая по модулю

        peaks = map(x -> x.pos-v.signal.validsegs[v.combo_item].ibeg+1, v.markup.Tone[v.combo_item])
    elseif sigtype == "pres"
        seg = v.signal.Pres[v.signal.validsegs[v.combo_item].ibeg:v.signal.validsegs[v.combo_item].iend]

        fsig_smooth = my_butter(seg, 2, 10, v.signal.fs, "low") # сглаживание
        sig = my_butter(fsig_smooth, 2, 0.3, v.signal.fs, "high") # устранение постоянной составляющей

        peaks = map(x -> (begs = x.ibeg-v.signal.validsegs[v.combo_item].ibeg+1, 
                        ends = x.iend-v.signal.validsegs[v.combo_item].ibeg+1), v.markup.Pres[v.combo_item])
    end

    v.dataforplotting = PlotElem(sig, peaks, bsad, bdad, (min = minimum(sig), max = maximum(sig)), length(sig))
end

function LoadButton( v::Globals)
    # Выбирается бинарь, из той же папки подтягивается соответсвующий hdr, 
    # из папки results (пока что) подтягивается соответствующий файл разметки
    if CImGui.Button("Загрузить запись")
        fname = open_dialog_native("Выберите bin-file");

        # чтение данных из файла
        signals, fs, _, _ = readbin(fname)

        ECG = signals[1]     # ЭКГ
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

end

function ui(v::Globals)
    MenuWindow(v)
    FigureWindow(v)
end
