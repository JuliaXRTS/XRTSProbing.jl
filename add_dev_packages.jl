using Pkg: Pkg

# QEDprocesses with perturbative Compton
#Pkg.add(; url="https://github.com/szabo137/QEDprocesses.jl", rev="dev-perturbative-compton")

# QEDcore with FlatPhaseSpaceLayout
Pkg.add(; url = "https://github.com/szabo137/QEDbase.jl", rev = "refac-test-implementation")

# QEDcore with FlatPhaseSpaceLayout
Pkg.add(; url = "https://github.com/szabo137/QEDcore.jl", rev = "add-rambo")
