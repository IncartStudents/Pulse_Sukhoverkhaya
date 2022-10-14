using Plots

include("../src/readfiles.jl")   # добавление файла с функциями чтения bin- и hdr-файлов
include("../src/help_func.jl")   # добавление файла со вспомогательными функциями 
include("../src/calc_AD.jl")     # добавления файла с функциями для работы с АД

binfile = "D:/INCART/Pulse_Data/bin/PX113211027165348.bin" # обрабатываемый бинарь (там же рядом должен лежать hdr!)

signals, fs, timestart, units = readbin(binfile);

Pres = signals.Pres; # давление
Tone = signals.Tone; # пульсации

ivalid = fill(false, length(Pres))
for i in 1:lastindex(ivalid) ivalid[i] = Pres[i] > 15 && Tone[i] > -1e7 end; # выбор валидных фрагментов

seg = IDtoSegBounds(ivalid) # получение границ валидных сегментов

MinSegLen = 30*fs   # минимальная допустимая длительность сегментов

len = seg[:,2] - seg[:,1]       # расчет длительностей сегментов
seg = seg[len .> MinSegLen, :]  # получение сегментов длиннее установленного порога

n = 2

s1 = Pres[seg[n,1]:seg[n,2]]

sfilt = SmoothFilt(s1, fs)
filtered  = ConstRemove(sfilt, fs)
filtered *= 1000

plot(s1)
plot!(sfilt)
plot!(filtered)

#  здесь должен быть вызов функций обработки каждого сегмента каждого сигнала в отдельности в отдельности
pk = pkAD(filtered, fs)
events = paramAD(pk, s1, filtered, fs)

bs = discardAD(events, fs)

