using DSP

function my_butter(sig, order, freq, fs, Ftype::Type{<:FilterType})
    responsetype = Ftype(freq; fs=fs)
    designmethod = Butterworth(order)
    fsig = filt(digitalfilter(responsetype, designmethod), sig)

    return fsig
end