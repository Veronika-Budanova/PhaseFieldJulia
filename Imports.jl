#=
Файл, содержащий подключение внешних пакетов.
Включен в Main.jl через include.
=#

using PyCall

plt = pyimport("matplotlib.pyplot")
patches = pyimport("matplotlib.patches")
plt.switch_backend("Agg")
animation = pyimport("matplotlib.animation")

using LinearAlgebra
using ForwardDiff
using Printf
using Random
using ProgressBars
using Profile
using Serialization
using Logging
using DifferentialEquations
using Polyester
using StaticArrays
using SymPy
