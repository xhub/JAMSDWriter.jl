@testset "nl_linearity" begin
    @testset "check simplification of formulae" begin
        # First term is true, we should choose `then`
        expr = :(ifelse(1 > 0, x[1] ^ 2, x[2] + 1))
        lin_expr = JAMSDWriter.LinearityExpr(expr)
        @test JAMSDWriter.convert_formula(lin_expr.c) == :(x[1] ^ 2)
        @test lin_expr.linearity == :nonlinear

        # First term is false, we should choose `else`
        expr = :(ifelse(1 < 0, x[1] ^ 2, x[2] + 1))
        lin_expr = JAMSDWriter.LinearityExpr(expr)
        @test JAMSDWriter.convert_formula(lin_expr.c) == :(x[2] + 1)
        @test lin_expr.linearity == :linear

        # First term isn't constant, we can't simplify
        expr = :(ifelse(1 > x[1], x[1] ^ 2, x[2] + 1))
        lin_expr = JAMSDWriter.LinearityExpr(expr)
        @test JAMSDWriter.convert_formula(lin_expr.c) == expr
        @test lin_expr.linearity == :nonlinear
    end
end
