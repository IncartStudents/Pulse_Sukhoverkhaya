# Поиск границ сегментов (валидных) по набору меток (1 или 0) айдишников отсчетов

# Вход: вектор булевых значений (меток валидности каждого отсчета)
# Выход: вектор границ сегментов, внутри которых заключены валидные отсчеты (окруженные с обеих сторон невалидными)
function IDtoSegBounds(ID::Vector{Bool})

    begs = fill(false, length(ID))
    ends = fill(false, length(ID))
    bounds = Vector{Int}[]

    if ID[1] begs[1] = true end
    for i in 2:lastindex(ID)
        if ID[i-1] && !ID[i]
            ends[i-1] = true
        elseif !ID[i-1] && ID[i]
            begs[i] = true
        end
    end
    if ID[end-1] ends[end] = true end

    bounds = hcat(findall(begs), findall(ends))

    return bounds
end

### test
# ID = [1,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0,1,0,0,1,1,1,1] |> Vector{Bool}
# findall(ID)
# begs, ends = IDtoSegBounds(ID)

##----------------------------------------------------##

# парсинг разметки из файлов
struct Markup
    bg::Int64
    en::Int64
    amp::Int64
end

function pls_ton_parse(filename)

    markup = Vector{Markup}[]

    open(filename) do file # Открываем файл
        
        one_valid_markup = Markup[]
        new_valid = false

        while !eof(file)

            line = readline(file)

            if length(split(line,"#")) != 1
                if new_valid
                    push!(markup, one_valid_markup)
                    one_valid_markup = Markup[]
                end
                new_valid = true
            else
                if new_valid 
                    values = split(line, "	")
                    vl = map(x -> parse(Int, x), values)
                    push!(one_valid_markup, Markup(vl[1], vl[2], vl[3]))
                end
            end

        end

        push!(markup, one_valid_markup)
    end

    return markup
end

function markup_seg(mkp)
    valid_seg = Vector[]
    for i in mkp
        # если важно сохранить количество сегментов (чтобы соотнести со вторым сигналом)
        # if isempty(i)
        #     push!(valid_seg, Int[])
        # else
        #     push!(valid_seg, [i[1].bg, i[end].bg])
        # end
        if !isempty(i)
            push!(valid_seg, [i[1].bg, i[end].bg])
        end
    end

    return valid_seg
end