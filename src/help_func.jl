using DSP
using Statistics

struct Markup
    bg::Int64
    en::Int64
    amp::Int64
end

struct PresEv   # каждое событие (пара) на сигнале давления в разметке
    ibeg::Int64 # максимум (начало события)
    iend::Int64 # минимум (конец события)
    bad::Int64  # число, в бинарной записи которого зашифрованы причины браковки
                # (каждая из 8 причин + нахождение слева от точки спуска + 
                # + отсутствие браковки. каждой причине соответвует определённая
                # степень двойки)
end

struct ToneEv  # каждое событие на сигнале тонов в разметке
    pos::Int64 # позиция пика
    bad::Int64 # см комментарии к одноименному полю структуры PresEv
end

struct Comp   # результат сравнения по каждой точке
    iref::Int64   # соответсвующие позиции референтной
    itest::Int64  # и тестовой меток (-1 в поле, если нет соответствия)
    bad::Int64    # для тестовых меток - информация о причине отбрабовки (см errors.txt)
    code::String # составной код ошибки (см errors.txt)
end

# Поиск границ сегментов (валидных) по набору меток (1 или 0) айдишников отсчетов

# Вход: вектор булевых значений (меток валидности каждого отсчета)
# Выход: вектор границ сегментов, внутри которых заключены валидные отсчеты (окруженные с обеих сторон невалидными)
function IDtoSegBounds(ID::Vector{Bool})

    begs = fill(false, length(ID))
    ends = fill(false, length(ID))

    if ID[1] begs[1] = true end
    for i in 2:lastindex(ID)
        if ID[i-1] && !ID[i]
            ends[i-1] = true
        elseif !ID[i-1] && ID[i]
            begs[i] = true
        end
    end
    if ID[end-1] ends[end] = true end

    bounds = map((x,y) -> (ibeg = x, iend = y), findall(begs), findall(ends))

    return bounds
end
##----------------------------------------------------##

# парсинг референтной разметки
function pls_ton_parse(filename)

    markup = Vector{Markup}[]

    open(filename) do file # Открываем файл
        
        one_valid_markup = Markup[]
        new_valid = false

        while !eof(file)

            line = readline(file)

            if length(split(line,"#")) != 1
                if new_valid && one_valid_markup != Markup[]
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

# парсинг тестовой разметки (pres и tone - каждый в вектор своих структур)
function test_markup_parse(filename)
    
    ispres = istone = false
    if split(filename, ".")[end] == "pres"
        ispres = true
        markup = Vector{PresEv}[]
    elseif split(filename, ".")[end] == "tone"
        istone = true
        markup = Vector{ToneEv}[]
    end

    open(filename) do file # Открываем файл
        
        if ispres one_valid_markup = PresEv[]
        elseif istone one_valid_markup = ToneEv[] end

        new_valid = false

        while !eof(file)

            line = readline(file)

            if length(split(line,"#")) != 1
                if new_valid && one_valid_markup != Markup[]
                    push!(markup, one_valid_markup)
                    if ispres one_valid_markup = PresEv[]
                    elseif istone one_valid_markup = ToneEv[] end
                end
                new_valid = true
            else
                if new_valid 
                    values = split(line, "   ")
                    vl = map(x -> parse(Int, x), values)
                    if length(vl) == 3
                        push!(one_valid_markup, PresEv(vl[1], vl[2], vl[3]))
                    elseif length(vl) == 2
                        push!(one_valid_markup, ToneEv(vl[1], vl[2]))
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

# получение по одному коду браковки кодов составляющих его причин
function errorcode(error, max_reasons_number)
    reasons = fill(false, max_reasons_number+1)
    # на всякий случай проверяем до max_reasons_number а не max_reasons_number - 1
    for i in 0:max_reasons_number
        reasons[i+1] = my_bitand(error, 2^i) != 0
    end

    reasons = findall(reasons).-1

    return reasons
end

# получение валидных сегментов
function get_valid_segments(Pres, Tone, pres_min_amp, tone_min_amp, min_seg_len)
    ivalid = fill(false, length(Pres))
    for i in 1:lastindex(ivalid) ivalid[i] = Pres[i] > pres_min_amp && Tone[i] > tone_min_amp end # выбор валидных фрагментов

    seg = IDtoSegBounds(ivalid) # получение границ валидных сегментов

    len = map(x -> x.iend - x.ibeg, seg)    # расчет длительностей сегментов
    valseg = map((x,y) -> x>min_seg_len ? y : false, len, seg)  # получение сегментов длиннее установленного порога
    valseg = valseg[findall(x -> x!=false, valseg)]

    return valseg
end

struct AD
    SAD::Int64
    DAD::Int64
end

### Сохранение разметки
function save_markup(markup, filename, ext)
    open(filename*ext, "w") do io

        if ext == ".pres"
            write(io, "beg   end   bad")
        elseif ext == ".tone"
            write(io, "pos   bad")
        end

        for j in 0:lastindex(markup)-1
            write(io, "\n"*"#$j")
            for i in markup[j+1]
                if ext == ".pres"
                    write(io, "\n$(i.ibeg)   $(i.iend)   $(i.bad)")
                elseif ext == ".tone"
                    write(io, "\n$(i.pos)   $(i.bad)")
                end
            end
        end
    end
end

function save_markup(filename, ad::Vector{NamedTuple{(:pump, :desc), NTuple{2, AD}}})
    open(filename, "w") do io
        write(io, "SADpump DADpump SADdesc DADdesc")
        for i in ad
            write(io, "\n$(i.pump.SAD)   $(i.pump.DAD)   $(i.desc.SAD)   $(i.desc.DAD)")
        end
    end
end

function read_ad(adtablefile)

    ad = AD[]

    open(adtablefile) do file # Открываем файл
        
        new_valid = false

        while !eof(file)

            line = readline(file)

            if length(split(line,"#")) != 1
                new_valid = true
            else
                if new_valid 
                    values = split(line, "\t")
                    # values = values[findall(x -> !isempty(x) && x!= "\t", values)]
                    vl = map(x -> parse(Int, x), values[6:7])
                    push!(ad, AD(vl[1], vl[2]))
                end
            end

        end
    end

    return ad
end

function read_alg_ad(adtablefile)
    ad = NamedTuple{(:pump, :desc), NTuple{2, AD}}[]
    needwrite = false
    open(adtablefile) do file # Открываем файл
        while !eof(file)
            line = readline(file)
            if needwrite
                values = split(line, "   ")
                values = map(x -> parse(Int, x), values)
                push!(ad, (pump = AD(values[1], values[2]), desc = AD(values[3], values[4])))
            end
            needwrite = true
        end
    end

    return ad
end

function get_ad_bounds(segm::Vector{Float64}, ad::AD, ispump::Bool)
    a = b = 0
    for i in 2:lastindex(segm)
        if ispump
            if segm[i] >= ad.DAD && segm[i-1] < ad.DAD b = i
            elseif segm[i] > ad.SAD && segm[i-1] <= ad.SAD a = i-1; break end
        else
            if segm[i] <= ad.SAD && segm[i-1] > ad.SAD a = i
            elseif segm[i] < ad.DAD && segm[i-1] >= ad.DAD b = i-1; break end
        end
    end

    if !ispump println("$(ad.SAD)   $(ad.DAD)") end

    bounds = (isad=a, idad=b)

    return bounds
end

function parse_compare_results(filename)

    markup = Vector{Comp}[]

    open(filename) do file # Открываем файл
        
        one_valid_markup = Comp[]
        new_valid = false

        while !eof(file)

            line = readline(file)

            if length(split(line,"#")) != 1
                if new_valid && one_valid_markup != Comp[]
                    push!(markup, one_valid_markup)
                    one_valid_markup = Comp[]
                end
                new_valid = true
            else
                if new_valid 
                    values = split(line, " ")
                    values = values[findall(x -> !isempty(x) && x!= " ", values)]
                    vl = map(x -> parse(Int, x), values[1:end-1])
                    if length(vl) == 3
                        push!(one_valid_markup, Comp(vl[1], vl[2], vl[3], values[end]))
                    # elseif length(vl) == 2
                    #     push!(one_valid_markup, Markup(vl[1], vl[2], 0))
                    # elseif length(vl) == 1
                    #     push!(one_valid_markup, Markup(vl[1], 0, 0))
                    end
                end
            end

        end

        push!(markup, one_valid_markup)
    end

    return markup
end

function save_compare_results(compare_result, filename, ext)
    open(filename*ext, "w") do io

        w = 6

        write(io, "$(rpad("ref", w))   $(rpad("test", w))   $(rpad("bad", w))   $(rpad("code", w))")

        for j in 0:lastindex(compare_result)-1
            write(io, "\n"*"#$j")
            for i in compare_result[j+1]
                ref = rpad(i.iref, w); test = rpad(i.itest, w)
                bad = rpad(i.bad, w); code = rpad(i.code, w)
                write(io, "\n$ref   $test   $bad   $code")
            end
        end
    end
end

function compare_result(ref, alg, rcomp, sigtype)
    # rcomp = 0.15*fs
    # если не совпадает число циклов измерения в разметках, берем меньшее
    N = minimum([lastindex(ref),lastindex(alg)])

    compare_result = fill(Comp[],N)

    for i in 1:N

        if sigtype == "tone"
            fc = map(x -> x.pos, alg[i])
            fcr = map(x -> x.bg, ref[i])
        elseif sigtype == "pres"
            fcr = map(x -> x.en, ref[i])
            fc = map(x -> x.ibeg, alg[i])
        else
            return
        end

        outcomp = calc_indexpairs(fcr, fc, rcomp)

        ln = length(outcomp)
        compare_result[i] = fill(Comp(0,0,0,""), ln)
        for j in 1:ln
            refpoint = outcomp[j][1]
            testpoint = outcomp[j][2]
            bad = 0; iref = itest = -1
            if refpoint != -1
                if sigtype == "tone" iref = ref[i][refpoint].bg
                elseif sigtype == "pres" iref = ref[i][refpoint].en end
            end
            if testpoint != -1
                if sigtype == "tone" itest = alg[i][testpoint].pos
                elseif sigtype == "pres" itest = alg[i][testpoint].ibeg end
                bad = alg[i][testpoint].bad
            end

            coderef = codetest = ' '
            if iref == -1 coderef = 'o' else coderef = 'e' end
            if itest == -1 codetest = 'o' else 
            if bad == 0 codetest = 'e' else codetest = 'b' end end

            compare_result[i][j] = Comp(iref,itest,bad,coderef*codetest)
        end
    end

    return compare_result
end

# инициализация состояния фильтра (эмпирический подбор)
# function filterInitState(b, a)
#     # функция эмпирически определяет начальное состояние фильтра
    
#     # возможны два варианта: ВЧ - начинается с нуля, НЧ - начинается с единицы
    
#     # чтобы использовать фильтр, необходимо умножить первый отсчет сигнала
#     # x(1) на его начальные состояния z: 
#     # z = x(1)*z;
    
#     x1 = 1; zero = 0;
#     x = ones(1,10);
#     order = length(b)-1;
#     while length(a)-1 < order
#         a = [a; 0];
#     end
#     if a[1] != 1
#         b = b./a[1];
#         a = a./a[1];
#     end
#     z0 = zeros(order, 1);
#     z1 = zeros(order, 1);
#     for k in 1:order # z - суммы в звеньях задержки (Direct Form II Transposed)
#        z0[k] = sum(x1*b[1+k:end]- zero*a[1+k:end]);
#        z1[k] = sum(x1*b[1+k:end]- x1*a[1+k:end]);
#     end

#     y0 = filt(b,a,x,z0);
#     y1 = filt(b,a,x,z1); 

#     if (1-y1[end]) > y0[end]
#         z = z0; # ВЧ / полосовые фильтры - инициализировать с нуля
#     else
#         z = z1; # НЧ фильтр - инициализировать с единицы * x(1) 
#     end

#     return z
# end