include("common.jl")

beta_arr = 10.0 .^ (-5:-1)

P = plot(
    xscale = :log10,
    xlab = L"\nu",
    ylab = L"g(\nu,\bar\beta)",
    legend = :topleft,
    legend_title = L"\bar\beta",
    palette = :seaborn_colorblind,
    #title = "\$-g(\\nu,\\bar\\beta)\$ (general vs. non-deg)"
)

for beta in reverse(beta_arr)
    #nu_arr = range(1e-3,10000,103)
    nu_arr = 10.0 .^ range(-4, 4, 300)

    gen_res = -general_integal.(nu_arr, beta)
    non_deg_res = -non_deg_integal.(nu_arr, beta)


    plot!(P, nu_arr, gen_res, lw = 3, lab = @sprintf("%01.0e", beta))
    plot!(P, nu_arr, non_deg_res, ls = :dash, c = :black, lab = "")
end


beta_arr = 10.0 .^ (-5:1)

P2 = plot(
    xscale = :log10,
    xlab = L"\nu",
    yscale = :log10,
    ylab = "relative error",
    legend = :outerright,
    legend_title = L"\bar\beta",
    palette = :seaborn_colorblind,
)

for beta in reverse(beta_arr)
    #nu_arr = range(1e-3,10000,103)
    nu_arr = 10.0 .^ range(-4, 4, 103)

    gen_res = general_integal.(nu_arr, beta)
    non_deg_res = non_deg_integal.(nu_arr, beta)

    rel_err = @. abs(gen_res - non_deg_res) / abs(gen_res)


    plot!(P2, nu_arr, rel_err, lw = 3, lab = @sprintf("%01.0e", beta))
end

l = @layout [a{0.45w} b]
Pall = plot(
    P,
    P2,
    plot_title = "\$-g(\\nu,\\bar\\beta)\$ (general vs. non-deg)",
    layout = l,
    size = (800, 400),
    bottom_margin = 10px,
)

file_path = joinpath(PLOTDIR, "g_func_compare.pdf")
savefig(Pall, file_path)
@info "saved in $file_path"
