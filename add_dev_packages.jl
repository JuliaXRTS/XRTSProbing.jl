using Pkg: Pkg

# QEDprocesses with perturbative Compton
#Pkg.add(; url="https://github.com/szabo137/QEDprocesses.jl", rev="dev-perturbative-compton")

Pkg.add(; url = "https://github.com/QEDjl-project/QEDbase.jl", rev = "dev") #Pkg.add(; path = "/Users/uwe/Work/jlProjects/QEDjl-project/forks/QEDbase.jl")
Pkg.add(; url = "https://github.com/QEDjl-project/QEDcore.jl", rev = "dev")
Pkg.add(; url = "https://github.com/QEDjl-project/QEDprocesses.jl", rev = "dev")
