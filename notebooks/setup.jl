import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


using DelimitedFiles
using LinearAlgebra
using Printf
using SparseArrays

using BSON
using HDF5
using Interpolations
using OrdinaryDiffEq, DiffEqCallbacks
using Postprocessing
using UnPack
using RecursiveArrayTools
using Roots
using StaticArrays
using SummationByPartsOperators
import SymPy; sp=SymPy

import FFTW; FFTW.set_num_threads(1)

using PyCall, LaTeXStrings; import PyPlot; plt=PyPlot
inset_locator = pyimport("mpl_toolkits.axes_grid.inset_locator")

cycler = pyimport("cycler").cycler
line_cycler   = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                 cycler(marker=["4", "2", "3", "1", "+", "x", "."]))

plt.rc("axes", prop_cycle=line_cycler)
plt.rc("text", usetex=true)
plt.rc("text.latex", preamble="\\usepackage{newpxtext}\\usepackage{newpxmath}\\usepackage{commath}\\usepackage{mathtools}")
plt.rc("font", family="serif", size=18.)
plt.rc("savefig", dpi=100)
plt.rc("legend", loc="best", fontsize="medium", fancybox=true, framealpha=0.5)
plt.rc("lines", linewidth=2.5, markersize=10, markeredgewidth=2.5)


function relaxation!(integrator)
    told = integrator.tprev
    uold = integrator.uprev
    tnew = integrator.t
    unew = integrator.u

    γ = one(tnew)
    terminate_integration = false
    γlo = one(γ)/2
    γhi = 3*one(γ)/2
    energy_old = relaxation_functional(uold, integrator.p)
    sign_condition = (relaxation_functional(γlo, unew, uold, integrator.p)-energy_old) * (relaxation_functional(γhi, unew, uold, integrator.p)-energy_old)
    if iszero(sign_condition)
        # everything is fine
        γ = one(tnew)
    elseif sign_condition > 0
        terminate_integration = true
        @warn "Terminating integration because no solution γ can be found."
    else
        γ = find_zero(g -> relaxation_functional(g, unew, uold, integrator.p)-energy_old, (γlo, γhi), Roots.AlefeldPotraShi())
    end
    if γ < eps(typeof(γ))
        terminate_integration = true
        @warn "Terminating integration because γ=$γ is too small."
    end

#     println(γ)
    @. unew = uold + γ * (unew - uold)
    DiffEqBase.set_u!(integrator, unew)
    if !(tnew ≈ first(integrator.opts.tstops))
        tγ = told + γ * (tnew - told)
        DiffEqBase.set_t!(integrator, tγ)
    end

    if terminate_integration
        terminate!(integrator)
    end

    nothing
end
