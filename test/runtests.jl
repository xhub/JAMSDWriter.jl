using JAMSDWriter, Compat, JuMP
using Base.Test

include("nl_convert.jl")
include("nl_linearity.jl")
include("nl_write.jl")

#solver = JuMP.UnsetSolver()
solvers = Any[]
push!(solvers, JAMSDSolver())

examples_path = joinpath(dirname(dirname(@__FILE__)), "examples")
for example in ["jump_nltrig.jl", "jump_nlexpr.jl"]
    include(joinpath(examples_path, example))
end
for example in [
                "jump_minlp.jl",
                "jump_nonlinearbinary.jl",
                "jump_no_obj.jl",
                "jump_const_obj.jl",
                "jump_pruning.jl"
                ]
    include(joinpath(examples_path, example))
end

# need to support n-ary max and min operators
#for example in ["jump_maxmin.jl"]
#    include(joinpath(examples_path, example))
#end

# The following list is obtained though
#    grep "solvers = Any" JuMP/test/solvers.jl
lp_solvers = Any[]
ip_solvers = Any[]
ip_dual_solvers = Any[]
semi_solvers = Any[]
sos_solvers = Any[]
lazy_solvers, lazy_soc_solvers, lazylocal_solvers, cut_solvers, cutlocal_solvers, heur_solvers, info_solvers = Any[], Any[], Any[], Any[], Any[], Any[], Any[]
quad_solvers = Any[]
quad_soc_solvers = Any[]
rsoc_solvers = Any[]
nlp_solvers = Any[]
convex_nlp_solvers = Any[]
minlp_solvers = Any[]
sdp_solvers = Any[]

push!(lp_solvers, JAMSDSolver())
push!(lp_solvers, JAMSDSolver("cplex"))
#push!(lp_solvers, JAMSDSolver("gurobi"))
push!(lp_solvers, JAMSDSolver("mosek"))
push!(lp_solvers, JAMSDSolver("xpress"))

push!(quad_solvers, JAMSDSolver())
push!(quad_solvers, JAMSDSolver("cplex"))
#push!(quad_solvers, JAMSDSolver("gurobi"))
push!(quad_solvers, JAMSDSolver("mosek"))
push!(lp_solvers, JAMSDSolver("xpress"))

quad_mip_solvers = copy(quad_solvers)

# Investigate if we can solve those
#soc_solvers = copy(quad_solvers)
soc_solvers = Any[]

push!(convex_nlp_solvers, JAMSDSolver())
push!(convex_nlp_solvers, JAMSDSolver("conopt"))
# buggy on my system
#push!(convex_nlp_solvers, JAMSDSolver("ipopt"))
push!(convex_nlp_solvers, JAMSDSolver("minos"))
push!(convex_nlp_solvers, JAMSDSolver("snopt"))
push!(convex_nlp_solvers, JAMSDSolver("mosek"))

push!(nlp_solvers, JAMSDSolver("baron"))
push!(nlp_solvers, JAMSDSolver("knitro"))

append!(convex_nlp_solvers, nlp_solvers)

push!(minlp_solvers, JAMSDSolver())



include(Pkg.dir("JuMP","test","qcqpmodel.jl"))

#
ipt = false

include("jump_test_nonlinear.jl")
