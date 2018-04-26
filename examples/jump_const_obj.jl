using JuMP, Base.Test, JAMSDWriter

# Example with no objective (#50)

if !isdefined(:solver); solver = JAMSDWriter.JAMSDSolver(); end

@testset "example: jump_no_obj" begin
    m = Model(solver=solver)
    @variable(m, 0 <= yp <= 1, Int)
    @variable(m, 0 <= l <= 1000.0)
    @variable(m, 0 <= f <= 1000.0)
    @NLconstraint(m, .087 * l >= f ^ 2)
    @constraint(m, l <= yp * 1000.0)
    @objective(m, Min, 5)

    @test solve(m) == :Optimal
    @test getobjectivevalue(m) == 5
end
