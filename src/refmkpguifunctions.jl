struct PresGuiMkp
    ibeg::Int64
    iend::Int64
    type::Int64
end

struct ToneGuiMkp
    pos::Int64
    type::Int64
end

struct Bounds    # Границы рабочей зоны или АД
    ibeg::Int64
    iend::Int64
end


function get_workzone_bounds(seg::Vector{Float64})
    lvlbeg_desc = round(maximum(seg) - 10) |> Int64; lvlend_desc = seg[end] >= 10 ? Int64(round(seg[end]))+10 : 10
    lvlbeg_pump = round(maximum(seg) - 10) |> Int64; lvlend_pump = seg[1] >= 10 ? Int64(round(seg[1]))+10 : 10
    
    return (pump = Bounds(lvlbeg_pump, lvlend_pump), desc = Bounds(lvlbeg_desc, lvlend_desc))
end

function SaveRefMarkup(filename::String, markup::Vector{PresGuiMkp}) # Сохранение исправленной референтной разметки
    open(filename, "w") do io

        write(io, "beg   end   type")

        for j in 1:lastindex(markup)
            write(io, "\n$(markup[j].ibeg)   $(markup[j].iend)   $(markup[j].type)")
        end

    end
end

function SaveRefMarkup(filename::String, markup::Vector{ToneGuiMkp}) # Сохранение исправленной референтной разметки
    open(filename, "w") do io

        write(io, "pos   type")

        for j in 1:lastindex(markup)
            write(io, "\n$(markup[j].pos)   $(markup[j].type)")
        end

    end
end

function SaveRefMarkup(filename::String, Presseg::Vector{Float64}, segbounds::Bounds, ad::NamedTuple{(:pump, :desc), Tuple{Bounds, Bounds}}, wr::NamedTuple{(:pump, :desc), Tuple{Bounds, Bounds}}) # Сохранение границ сегмента, рабочей зоны и АД
    
    ibeg = segbounds.ibeg
    iend = segbounds.iend
    
    isadpump = ad.pump.ibeg; Asadpump = round(Int, Presseg[isadpump])
    idadpump = ad.pump.iend; Adadpump = round(Int, Presseg[idadpump])
    isaddesc = ad.desc.ibeg; Asaddesc = round(Int, Presseg[isaddesc])
    idaddesc = ad.desc.iend; Adaddesc = round(Int, Presseg[idaddesc])

    iwbegpump = wr.pump.ibeg; Awbegpump = round(Int, Presseg[iwbegpump])
    iwendpump = wr.pump.iend; Awendpump = round(Int, Presseg[iwendpump])
    iwbegdesc = wr.desc.ibeg; Awbegdesc = round(Int, Presseg[iwbegdesc])
    iwenddesc = wr.desc.iend; Awenddesc = round(Int, Presseg[iwenddesc])

    open(filename, "w") do io
        write(io, "       segmbegpos   segmendpos   iSAD   ADamp   iDAD   DADamp   iWRh   WRh   iWRl   WRl")
        write(io, "\npump   $ibeg   $iend   $isadpump   $Asadpump   $idadpump   $Adadpump   $iwendpump   $Awendpump   $iwbegpump   $Awbegpump")
        write(io, "\ndesc   $ibeg   $iend   $isaddesc   $Asaddesc   $idaddesc   $Adaddesc   $iwbegdesc   $Awbegdesc   $iwenddesc   $Awenddesc")
    end

end

function ReadRefMkp(filename::String)

    key = split(split(filename, ".")[end-1], "/")[end] |> Symbol

    if key == :bounds 
        bnds = (pump = Bounds(0,0), desc = Bounds(0,0)) 
        markup = (ad = bnds, wz = bnds, segm = Bounds(0,0))
        iad = fill(Bounds(0,0), 2); iwz = fill(Bounds(0,0), 2)
        ad = fill(Bounds(0,0), 2); wz = fill(Bounds(0,0), 2)
        sb = Bounds(0,0)
    elseif key == :pres markup = PresGuiMkp[]
    elseif key == :tone markup = ToneGuiMkp[]
    else return "Invalid input!" end

    open(filename) do file # Открываем файл
        needwrite = false
        k = 0
        while !eof(file)
            line = readline(file)
            if needwrite
                vars = split(line, "   ")
                if key == :bounds vars = vars[2:end] end
                vars = map(x -> parse(Int64, x), vars)

                if key == :pres push!(markup, PresGuiMkp(vars[1], vars[2], vars[3]))
                elseif key == :tone push!(markup, ToneGuiMkp(vars[1], vars[2]))
                elseif key == :bounds 
                    sb = Bounds(vars[1], vars[2])
                    iad[k] = Bounds(vars[3], vars[5])
                    iwz[k] = Bounds(vars[7], vars[9])
                    ad[k] = Bounds(vars[4], vars[6])
                    wz[k] = Bounds(vars[8], vars[10])
                end
            end
            needwrite = true
            k += 1
        end
        if key == :bounds
            markup = (iad = (pump = iad[1], desc = iad[2]), iwz = (pump = iwz[1], desc = iwz[2]), 
                        ad = (pump = ad[1], desc = ad[2]), wz = (pump = wz[1], desc = wz[2]), 
                        segm = sb)
        end
    end

    return markup
end

# ReadRefMkp("formatted alg markup/bin/PX11321102817293/1/bounds.csv")