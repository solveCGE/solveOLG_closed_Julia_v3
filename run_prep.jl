using Dates
using Statistics
using NLsolve
using Printf
using Distributed
using Plots
using Parameters

# load some handy auxiliary functions
include("loadfunctions.jl")

# control center
const tend            = 300   # number of periods
const nag             = 100   # number of age groups (nag = x => max age = x-1)
const budget_bal      = 3     # budget closing instrument (1.. tauW, 2.. tauF, 3.. tauC, 4.. taul, 5.. tauprof, 6.. cG) to fulfill given debt path
const genplots        = false # generate plots

# initializes all variables globally
include("initdata.jl")

# load main functions
include("algo.jl")
include("hh.jl")
include("firm.jl")