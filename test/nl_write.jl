facts("[nl_write] Temp file handling") do
    # Turn on debug mode so files persist
    old_debug = JAMSDWriter.CONFIG[:debug]
    JAMSDWriter.setdebug(true)

#    filename = "test"
#    filepath = joinpath(JAMSDWriter.solverdata_dir, filename)
#    JAMSDWriter.clean_solverdata()

#    context("all temp files deleted successfully") do
#        @fact length(readdir(JAMSDWriter.solverdata_dir)) --> 1
#    end

    solver = JAMSDWriter.JAMSDSolver()
    m = Model(solver=solver)
    @variable(m, x >= 0)
    @objective(m, Min, x)
    solve(m)

    # Reset debug mode and clean up
    JAMSDWriter.setdebug(old_debug)
    JAMSDWriter.clean_solverdata()
end
