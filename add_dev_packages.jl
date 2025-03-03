using Pkg: Pkg

# QEDprocesses with perturbative Compton
#Pkg.add(; url="https://github.com/szabo137/QEDprocesses.jl", rev="dev-perturbative-compton")

# QEDcore with FlatPhaseSpaceLayout

Pkg.add(; url = "https://github.com/szabo137/QEDbase.jl", rev = "extent-psl-api")
#Pkg.add(; path = "/Users/uwe/Work/jlProjects/QEDjl-project/forks/QEDbase.jl")

# QEDcore with FlatPhaseSpaceLayout
Pkg.add(; url = "https://github.com/QEDjl-project/QEDcore.jl", rev = "dev")
