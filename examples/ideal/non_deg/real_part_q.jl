include("common.jl")

bbar_arr = 10.0 .^ (-4:0)
ombar_arr = 10.0 .^ (-2:1)
qbar_arr = 10.0 .^ range(-7, 5, 100)

PLOTS = []
for bbar in bbar_arr
    for ombar in ombar_arr
        @show bbar
        @show ombar
        P = plot(
            xlab = L"\bar q",
            xscale = :log10,
            ylab = L"\Re \chi_0(\bar \omega, \bar q, \bar \beta)",
            palette = :seaborn_colorblind,
            title = "\$\\bar\\omega = $(@sprintf("%01.0e", ombar)) \$, \$\\bar\\beta = $(@sprintf("%01.0e", bbar)) \$",
            titlefont = font(11, "Computer Modern"),
            guidefont = font(9, "Computer Modern"),
        )
        #gen_rf = -QEDprobing.real_lindhard_nonzero_temperature.(ombar,qbar_arr,bbar)
        gen_rf = -_real_lindhard_new.(ombar, qbar_arr, bbar)
        #plot!(P,qbar_arr,gen_rf,lw=5,lab=@sprintf("%01.0e",ombar),c=:orange)
        plot!(P, qbar_arr, gen_rf, lw = 5, lab = "", c = :orange)


        non_deg_rf = -_real_lindhard_non_deg.(ombar, qbar_arr, bbar)
        plot!(P, qbar_arr, non_deg_rf, ls = :dot, c = :blue, lab = "")

        push!(PLOTS, P)
    end
end

Pall = plot(
    PLOTS...,
    plot_title = "Ideal Electron Gas, no screening (general vs. non. deg.)",
    layout = (length(bbar_arr), length(ombar_arr)),
    size = (1300, 1200),
)

file_path = joinpath(PLOTDIR, "real_part_q.pdf")
savefig(Pall, file_path)
@info "saved in $file_path"
