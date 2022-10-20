using Dates, CSV, Tables, DataFrames

"""
словарь с кодировками типов, которые используются в hdr-файлах
"""
const string2datatype = Dict{String, DataType}(
    "int8"    => Int8,
    "uint8"   => UInt8,
    "int16"   => Int16,
    "uint16"  => UInt16,
    "int32"   => Int32,
    "uint32"  => UInt32,
    "int64"   => Int64,
    "uint64"  => UInt64,
    "float"   => Float32,
    "float32" => Float32,
    "double"  => Float64,
    "float64" => Float64,
    "string" => String
)

"""
чтение hdr-файла заголовка
"""
function readhdr(filepath::AbstractString)
    open(filepath, "r") do io
        out = readhdr(io)
    end
end

function readhdr(io::IO)

    lines = [rstrip(readline(io)) for _ in 1:5] #, enc"windows-1251") # read and decode from windows-1251 to UTF-8 string

    delim = (' ', '\t')
    ln = split(lines[1], delim)
    num_ch, fs, lsb = parse(Int, ln[1]), parse(Float64, ln[2]), parse(Float64, ln[3])
    type = Int32
    if (length(ln) > 3) # optional field
        type = string2datatype[ln[4]]
    end

    ln = split(lines[2], delim)
    ibeg, iend = parse(Int, ln[1]), parse(Int, ln[2])
    timestart = parse(DateTime, ln[3])

    names = String.(split(lines[3], delim))
    lsbs = parse.(Float64, split(lines[4], delim))
    units = String.(split(lines[5], delim))

    if num_ch != length(names) # фикс, если в начале указано неверное кол-во каналов
        num_ch = length(names)
    end

    # опционально - читаем следующую строку c типами, если она подходит под описание
    line = rstrip(readline(io))
    t_ = String.(split(line, delim))
    types = if length(t_) == num_ch && t_[1]!="tohead"
        map(x->string2datatype[x], t_)
    else
        fill(type, num_ch)
    end

    length(names) == length(lsbs) == length(units) || error("разное количество полей")

    return num_ch, fs, ibeg, iend, timestart, names, lsbs, units, type, types
end


"""
чтение bin-файла с каналами, рядом должен лежать hdr-файл
"""
function readbin(filepath::AbstractString, range::Union{Nothing, UnitRange{Int}} = nothing)
    # защита от дурака
    fpath, ext = splitext(filepath)
    hdrpath = fpath * ".hdr"
    binpath = fpath * ".bin"

    num_ch, fs, _, _, timestart, names, lsbs, units, type, _ = readhdr(hdrpath)

    offset = (range !== nothing) ? range.start - 1 : 0
    
    elsize = num_ch * sizeof(type)
    byteoffset = offset * elsize # 0-based
    maxlen = (filesize(binpath) - byteoffset) ÷ elsize # 0-based
    len = (range !== nothing) ? min(maxlen, length(range)) : maxlen

    if len <= 0
        data = Matrix{type}(undef, num_ch, 0)
    else
        data = Matrix{type}(undef, num_ch, len)
        open(binpath, "r") do io
            seek(io, byteoffset)
            read!(io, data)
        end
    end

    channels = [data[ch, :] .* lsbs[ch] for ch in 1:num_ch] |> Tuple # matrix -> vector of channel vectors
    sym_names = Symbol.(names) |> Tuple # column names: String -> Symbol 
    
    named_channels = NamedTuple{sym_names}(channels)
    return named_channels, fs, timestart, units
end

# чтение комбинированного текстового формата, состоящего из хедера и данных
function readevt(filepath::String)
    open(filepath, "r") do io
        # complex code to read some_data
        num_ch, fs, ibeg, iend, timestart, names, lsbs, units, type, types = readhdr(io)

        # будем читать строки до тех пор, пока не стречим разделитель хедера и данных
        while readline(io) !== "#" 
        end
        #@show line = readline(io)
        #@show line !== "#"

        sym_names = Symbol.(names)
        #@show sym_names
        pairs = map(zip(sym_names, types)) do (k, v)
            k => v
        end

        data = CSV.File(io; 
            header = sym_names, 
            types = Dict(pairs), 
            delim = '\t',
        ) |> columntable # named tuple of columns

        return data, fs, timestart, units
    end
end


#############
##(ADDED)####
function readmarkup(fname::String)
    markup = CSV.read(fname, DataFrame)
    p = markup.Peaks  
    
    return p
end

