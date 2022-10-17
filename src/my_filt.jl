using DSP

function my_butter(sig, order, freq, fs, type)
    if type == "low" responsetype = Lowpass(freq; fs=fs)
    elseif type == "high" responsetype = Highpass(freq; fs=fs)
    end
    designmethod = Butterworth(order)
    fsig = filt(digitalfilter(responsetype, designmethod), sig)

    return fsig
end