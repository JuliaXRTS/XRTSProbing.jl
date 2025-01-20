using Pkg: Pkg

# QEDprocesses with perturbative Compton
#Pkg.add(; url="https://github.com/szabo137/QEDprocesses.jl", rev="dev-perturbative-compton")

# QEDcore with FlatPhaseSpaceLayout
Pkg.add(; url = "https://github.com/QEDjl-project/QEDbase.jl", rev = "dev")

# QEDcore with FlatPhaseSpaceLayout
Pkg.add(; url = "https://github.com/szabo137/QEDcore.jl", rev = "add-rambo")
