using DSP
using Statistics

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
                    if length(vl) == 3
                        push!(one_valid_markup, Markup(vl[1], vl[2], vl[3]))
                    elseif length(vl) == 2
                        push!(one_valid_markup, Markup(vl[1], vl[2], 0))
                    elseif length(vl) == 1
                        push!(one_valid_markup, Markup(vl[1], 0, 0))
                    end
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

# побитовое И
function my_bitand(a::Int, b::Int)

    a_bit = bitstring(a)
    b_bit = bitstring(b)

    comp = map((x,y) -> x==y=='1' ? true : false, a_bit, b_bit)
    pow = length(a_bit) .- findall(comp)

    res = 0
    for i in pow res += 2^i end

    return res
end

# получение валидных сегментов
function get_valid_segments(Pres, Tone, pres_min_amp, tone_min_amp, min_seg_len)
    ivalid = fill(false, length(Pres))
    for i in 1:lastindex(ivalid) ivalid[i] = Pres[i] > pres_min_amp && Tone[i] > tone_min_amp end # выбор валидных фрагментов

    seg = IDtoSegBounds(ivalid) # получение границ валидных сегментов

    MinSegLen = 30*fs   # минимальная допустимая длительность сегментов

    len = seg[:,2] - seg[:,1]       # расчет длительностей сегментов
    seg = seg[len .> min_seg_len, :]  # получение сегментов длиннее установленного порога
end

### Сохранение разметки
function save_markup(markup, filename, ext)
    open(filename*ext, "w") do io

        if ext == ".pres"
            write(io, "beg	end")
        elseif ext == ".tone"
            write(io, "peak")
        end

        for j in 0:lastindex(markup)-1
            write(io, "\n"*"#$j")
            for i in markup[j+1]
                if ext == ".pres"
                    write(io, "\n"*string(i.min_pos)*"	"*string(i.max_pos))
                elseif ext == ".tone"
                    write(io, "\n"*string(i))
                end
            end
        end
    end
end

# инициализация состояния фильтра (эмпирический подбор)
function filterInitState(b, a)
    # функция эмпирически определяет начальное состояние фильтра
    
    # возможны два варианта: ВЧ - начинается с нуля, НЧ - начинается с единицы
    
    # чтобы использовать фильтр, необходимо умножить первый отсчет сигнала
    # x(1) на его начальные состояния z: 
    # z = x(1)*z;
    
    x1 = 1; zero = 0;
    x = ones(1,10);
    order = length(b)-1;
    while length(a)-1 < order
        a = [a; 0];
    end
    if a[1] != 1
        b = b./a[1];
        a = a./a[1];
    end
    z0 = zeros(order, 1);
    z1 = zeros(order, 1);
    for k in 1:order # z - суммы в звеньях задержки (Direct Form II Transposed)
       z0[k] = sum(x1*b[1+k:end]- zero*a[1+k:end]);
       z1[k] = sum(x1*b[1+k:end]- x1*a[1+k:end]);
    end

    y0 = filt(b,a,x,z0);
    y1 = filt(b,a,x,z1); 

    if (1-y1[end]) > y0[end]
        z = z0; # ВЧ / полосовые фильтры - инициализировать с нуля
    else
        z = z1; # НЧ фильтр - инициализировать с единицы * x(1) 
    end

    return z
end