include("common.jl")

bbar_arr = 10.0 .^ (-4:0)
qbar_arr = 10.0 .^ (-2:2)

ombar_arr = 10.0 .^ range(-10, 10, 300)

PLOTS = []
for bbar in bbar_arr
    for qbar in qbar_arr
        @show bbar
        @show qbar
        P = plot(
            xlab = L"\bar \omega",
            xscale = :log10,
            ylab = L"\Re \chi_0(\bar \omega, \bar q, \bar \beta)",
            palette = :seaborn_colorblind,
            title = "\$\\bar q= $(@sprintf("%01.0e", qbar)) \$, \$\\bar\\beta = $(@sprintf("%01.0e", bbar)) \$",
            titlefont = font(11, "Computer Modern"),
            guidefont = font(9, "Computer Modern"),
        )
        #gen_rf = -QEDprobing.real_lindhard_nonzero_temperature.(ombar,qbar_arr,bbar)
        gen_rf = -_real_lindhard_new.(ombar_arr, qbar, bbar)
        #plot!(P,qbar_arr,gen_rf,lw=5,lab=@sprintf("%01.0e",ombar),c=:orange)
        plot!(P, ombar_arr, gen_rf, lw = 5, lab = "", c = :orange)


        non_deg_rf = -_real_lindhard_non_deg.(ombar_arr, qbar, bbar)
        plot!(P, ombar_arr, non_deg_rf, ls = :dot, c = :blue, lab = "")

        push!(PLOTS, P)
    end
end

Pall = plot(
    PLOTS...,
    plot_title = "Ideal Electron Gas, no screening (general vs. non. deg.)",
    layout = (length(bbar_arr), length(qbar_arr)),
    size = (1400, 1200),
)

file_path = joinpath(PLOTDIR, "real_part_om.pdf")
savefig(Pall, file_path)
@info "saved in $file_path"
