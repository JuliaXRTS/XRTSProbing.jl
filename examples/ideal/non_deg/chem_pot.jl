include("common.jl")


#beta_arr = range(0.1,3,length=300)
beta_arr = 10.0 .^ range(-6, 0, 300)
c = _chemical_potential_normalized.(beta_arr)
c_non_deg = _non_deg_chem_pot.(beta_arr)

P = plot(
    xlab = L"\bar\beta",
    xscale = :log10,
    ylab = L"\bar\mu(\bar\beta)",
    legend = :bottomright,
    #title = "Ideal Electron gas, no screening"
)

plot!(P, beta_arr, c, lab = "Ichimura (2018)", lw = 3, c = :orange)

plot!(P, beta_arr, c_non_deg, lab = "non-deg.", c = :black, ls = :dash)


beta_arr = 10.0 .^ range(0, 6, 300)

c = _chemical_potential_normalized.(beta_arr)
c_non_deg = _non_deg_chem_pot.(beta_arr)

P2 = plot(
    xlab = L"\bar\beta",
    xscale = :log10,
    ylab = L"\bar\mu(\bar\beta)",
    legend = :bottomright,
    #title = "Ideal Electron gas, no screening"
)

plot!(P2, beta_arr, c, lab = "Ichimura (2018)", lw = 3, c = :orange)

plot!(P2, beta_arr, c_non_deg, lab = "non-deg.", c = :black, ls = :dash)


Pall = plot(
    P,
    P2,
    plot_title = "Ideal Electron gas, no screening",
    size = (800, 500),
    bottom_margin = 10px,
)
file_path = joinpath(PLOTDIR, "chemical_potential.pdf")
savefig(Pall, file_path)
@info "saved in $file_path"
