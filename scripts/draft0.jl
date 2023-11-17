using Plots

x = range(0, 100, step=1)
y = fill(6480000, length(x))

plot(x,y, label="Затраты")
plot!(x, x*100000, label="Выручка")
xlabel!("Объём, шт")
ylabel!("Доход/издержки, руб")
xlims!()