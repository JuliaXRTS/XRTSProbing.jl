include("common.jl")

beta_arr = 10.0 .^ range(-6, 4, length = 300)
c = general_chem_pot_check.(beta_arr)
c_non_deg = non_deg_chem_pot_check.(beta_arr)

P = plot(
    xlab = L"\bar\beta",
    xscale = :log10,
    ylab = L"\vert 3 \int_0^\infty x^2 f(x^2,\bar\beta)\,\mathrm{d}x - 1\vert",
    yscale = :log10,
    ylim = (1.0e-10, 1.1),
    legend = :topleft,
    plot_title = "Ideal Electron gas, no screening",
    size = (800, 500),
    bottom_margin = 10px,
    left_margin = 20px,
)

plot!(P, beta_arr, c, lab = "Ichimura (2018)", lw = 3, c = :orange)

plot!(P, beta_arr, c_non_deg, lab = "non-deg.", c = :black, ls = :dash)

file_path = joinpath(PLOTDIR, "chemical_potential_test.pdf")
savefig(P, file_path)
@info "saved in $file_path"
