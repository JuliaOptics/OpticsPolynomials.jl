using Documenter, OpticsPolynomials

makedocs(modules = [OpticsPolynomials], doctest = false, sitename = "OpticsPolynomials")
deploydocs(repo = "github.com/JuliaOptics/OpticsPolynomials.jl.git", devbranch="main")
