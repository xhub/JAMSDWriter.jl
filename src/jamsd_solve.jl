function jamsd_set_modeltype(m::JAMSDMathProgModel)
    discrete = any((m.vartypes .== :Int) + (m.vartypes .== :Bin) .> 0)
    if discrete
        if m.model_type == qcp
            m.model_type = miqcp
        elseif m.model_type == nlp
             m.model_type = minlp
        end
    end
    jamsd_set_modeltype(m.jamsd_ctx, m.model_type)
    res = ccall((:ctx_setobjsense, "libjamsd"), Cint, (Ptr{context}, Cint), m.jamsd_ctx, sense_to_jamsd[m.sense])
    res != 0 && error("return code $res from JAMSD")
end

function jamsd_set_modeltype(ctx::Ptr{context}, model_type)
    res = ccall((:ctx_setmodeltype, "libjamsd"), Cint, (Ptr{context}, Cint), ctx, model_type)
    res != 0 && error("return code $res from JAMSD")
end

function jamsd_solve(ctx::Ptr{context}, ctx_dest::Ptr{context}, solver_name::String, emp::Ptr{Void}=C_NULL)
    res = ccall((:model_compress, "libjamsd"), Cint, (Ptr{context}, Ptr{context}, Ptr{Void}, Ptr{Void}), ctx, ctx_dest, emp, C_NULL)
    res != 0 && error("return code $res from JAMSD")
    res = ccall((:ctx_exportmodel, "libjamsd"), Cint, (Ptr{context}, Ptr{context}), ctx, ctx_dest)
    res != 0 && error("return code $res from JAMSD")

    JAMSDWriter.emp_hack(emp)
    ccall((:ctx_writemodel, "libjamsd"), Cint, (Ptr{JAMSDWriter.context}, Cstring), ctx_dest, "validation.gms")
    ccall((:ctx_setsolverstr, "libjamsd"), Cint, (Ptr{JAMSDWriter.context}, Cstring), ctx_dest, "CONVERTD")
    ccall((:ctx_callsolver, "libjamsd"), Cint, (Ptr{context},), ctx_dest)
    # switch back to the default solver

    ccall((:ctx_setsolverstr, "libjamsd"), Cint, (Ptr{JAMSDWriter.context}, Cstring), ctx_dest, solver_name)

    if emp != C_NULL
        res = ccall((:empinfo_solve, "libjamsd"), Cint, (Ptr{empinfo},), emp)
        return res
    else
        return ccall((:ctx_callsolver, "libjamsd"), Cint, (Ptr{context},), ctx_dest)
    end
end

function jamsd_setup_gams()
    ctx = ccall((:ctx_alloc, "libjamsd"), Ptr{context}, (Cuint,), 0)

    gamscntr_template = readstring(joinpath(solverdata_dir, "gamscntr.dat"))
    if isdir("gams_dir")
        rm("gams_dir", recursive=true)
    end
    mkdir("gams_dir")
    gamscntr_file = open("gams_dir/gamscntr.dat", "w")
    cur_dir = pwd() * "/gams_dir"

    println(gamscntr_file, replace(gamscntr_template, r"@@SUB@@", cur_dir))
    close(gamscntr_file)

    res = ccall((:gams_set_gamscntr, "libjamsd"), Cint, (Ptr{context}, Cstring), ctx, pwd() * "/gams_dir/gamscntr.dat")
    res != 0 && error("return code $res from JAMSD")

    gamsdir = split(gamscntr_template, "\n")[29]

    res = ccall((:gams_set_gamsdir, "libjamsd"), Cint, (Ptr{context}, Cstring), ctx, gamsdir)
    res != 0 && error("return code $res from JAMSD")

    ENV["PATH"] *= ":" * gamsdir

    return ctx
end
