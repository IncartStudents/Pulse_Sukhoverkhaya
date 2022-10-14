
# примеры
map(+, [1, 2, 3], [10, 20, 30])
map((x) -> x * 2, [1, 2, 3])

function BoundsCoeff(f, wavetype)
    kb = []
    ke = []
    if wavetype == "P"
        kb = 1.35
        ke = 2.0
    elseif wavetype == "Q"
        kb = 1.8
        ke = 1
    elseif wavetype == "S"
        kb = 1
        if f < 4
            ke = 3
        elseif f >= 4 && f < 4.75
            ke = 8
        elseif f >= 4.75 && f < 6.20
            ke = 9
        elseif f >= 6.20
            ke = 12
        end
    elseif wavetype == "T"
        kb = 2
        if f <= 0.13 
            ke = 4
        elseif f > 0.13 f < 0.20
            ke = 5
        elseif f >= 0.20 f < 0.41
            ke = 6
        elseif f >= 0.41
            ke = 7
        end
    end

    return kb, ke
end

kb,ke = map(BoundsCoeff, [1,2,6,4,8], ["P","Q","R","S","T"])

function demo_fcn(in)

    out1 = in(1);
    out2 = in(2);

    return out1, out2
end

a = [1,2,3,4,5,6,7,8,9,10]
sample = map(x -> x>5 && x<10 && x>1, a)
a = a[sample]