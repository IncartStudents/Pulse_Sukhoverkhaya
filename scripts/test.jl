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
include("../ECG markup/onelead.jl")


function ui() # Главное окно программы (Окно меню + окно графиков с разметкой)
    CImGui.Begin("A")
    CImGui.End()

    CImGui.Begin("B")
    CImGui.End()
end

function show_gui() # Main
    Renderer.render(
        ()->ui(),
        width = 1800,
        height = 1600,
        title = "",
        v = Renderer.GR()
    )
end

show_gui();
