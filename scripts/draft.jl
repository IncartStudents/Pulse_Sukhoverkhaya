using ImPlot

include("../src/my_filt.jl")  
include("../src/readfiles.jl") 

signals, fs, _, _ = readbin("D:/INCART/Pulse_Data/все базы/Шумовая база/Кустовская_23-01-21_11-29-20_.bin")
segm = signals[1][75211:75211+7000]
fsegm = my_butter(segm, 2, 20, fs, "low")
plot(segm)
plot!(fsegm[2000:end])