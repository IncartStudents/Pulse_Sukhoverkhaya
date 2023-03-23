using CImGui
using CImGui: ImVec2, ImVec4
using ImPlot
using ImPlot.LibCImGui: ImGuiCond_Once, ImGuiCond_Always, ImPlotAxisFlags_NoGridLines

include("../src/Renderer.jl")
using .Renderer

mutable struct Globals
    a::Int

    function Globals()
        a = 0
        new(a)
    end
end



function ui(v::Globals) # Главное окно программы (Окно меню + окно графиков с разметкой)
    CImGui.Begin("test")
        CImGui.InputInt("1", v.a)
    CImGui.End()
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