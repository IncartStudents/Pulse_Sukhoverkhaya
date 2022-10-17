# функция сравнения двух разметок
#находим индексы пар для двух столбцов
# алгоритм из п. 201.12.1.101.2.3.2
# https://docs.cntd.ru/document/1200157909 
function calc_indexpairs(
    times1::Vector, # позиции в первой разметке
    times2::Vector, # позиции во второй разметке
    t_radius #радиус поиска соседа (в точках. t_radius = t_sec*Fs)
 )

    N1 = length(times1)
    N2 = length(times2)

    indexpairs = Vector{NTuple{2,Int}}(undef, N1 + N2) # для хранения индексов пар
  
    pair_ind = (-1, -1)
    Np = 0 # N1 + N2
    i1 = 1
    i2 = 1

    T = isempty(times1) ? 0 : first(times1) #если пустой times1 - T=0, иначе первый элемент times1
    t = isempty(times2) ? 0 : first(times2) 

    stop1 = false
    stop2 = false
    # работаем по алгоритму 
    while !(stop1 && stop2)
        if t <= T   # найдено раньше эталона
            if i2 < N2
                t_ = times2[i2 + 1]
            else
                t_ = t + Inf
            end
            delta = abs(T - t)
            delta_ = abs(T - t_)
            if delta < delta_ && delta < t_radius
                pair_ind = (i1, i2)
                if i1 < N1
                    i1 = i1 + 1
                else
                    stop1 = true
                end
                T = times1[i1]
            else
                if stop2
                    pair_ind = (-1, -1)
                else
                    pair_ind = (-1, i2)
                end
            end
            if i2 < N2
                i2 = i2 + 1
            else
                stop2 = true
            end
            t = t_ # == time2[i2]

        else # T < t
            if i1 < N1
                T_ = times1[i1 + 1]
            else
                T_ = T + Inf #3*t150ms;
            end
            delta = abs(T - t)
            delta_ = abs(T_ - t)
            if delta < delta_ && delta < t_radius
                # T-t paired
                pair_ind = (i1, i2);
                if i2 < N2
                    i2 = i2 + 1
                else
                    stop2 = true
                end
                t = times2[i2]
            else
                # T-O/X extra detecton
                if stop1
                    pair_ind = (-1, -1)
                else
                    pair_ind = (i1, -1)
                end
            end
            if i1 < N1
                i1 = i1 + 1
            else
                stop1 = true
            end
            T = T_; # == times1[i1]
        end
        if pair_ind[1] == -1 && pair_ind[2] == -1
            continue;
        end
        Np = Np + 1

        indexpairs[Np] = pair_ind # уже аллоцировали заранее
    end

    resize!(indexpairs, Np)

    return indexpairs
end

