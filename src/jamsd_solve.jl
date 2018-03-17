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
    res = ccall((:ctx_setobjsense, jamsd_libname), Cint, (Ptr{context}, Cint), m.jamsd_ctx, sense_to_jamsd[m.sense])
    res != 0 && error("return code $res from JAMSD")
end

function jamsd_set_modeltype(ctx::Ptr{context}, model_type)
    res = ccall((:ctx_setmodeltype, jamsd_libname), Cint, (Ptr{context}, Cint), ctx, model_type)
    res != 0 && error("return code $res from JAMSD")
end

function jamsd_solve(ctx::Ptr{context}, ctx_dest::Ptr{context}, solver_name::String, emp::Ptr{empinfo}=Ptr{empinfo}(C_NULL))

    CONFIG[:solver_log] && hack_solver_log()

    if emp != C_NULL
        res = ccall((:empinfo_transform, jamsd_libname), Cint, (Ptr{empinfo}, Ptr{context}), emp, ctx_dest)
        res != 0 && error("return code $res from JAMSD")

        hack_exportempinfo(ctx_dest, ctx, emp)

        if CONFIG[:export_gms]
            ccall((:ctx_writemodel, jamsd_libname), Cint, (Ptr{context}, Cstring), ctx_dest, "validation.gms")
            ccall((:ctx_setsolverstr, jamsd_libname), Cint, (Ptr{context}, Cstring), ctx_dest, "CONVERTD")
            ccall((:ctx_callsolver, jamsd_libname), Cint, (Ptr{context},), ctx_dest)
        end

        ccall((:ctx_setsolverstr, jamsd_libname), Cint, (Ptr{JAMSDWriter.context}, Cstring), ctx_dest, solver_name)

        res = ccall((:empinfo_solve, jamsd_libname), Cint, (Ptr{empinfo},), emp)
        return res
    else
        res = ccall((:model_compress, jamsd_libname), Cint, (Ptr{context}, Ptr{context}, Ptr{empinfo}, Ptr{Void}), ctx, ctx_dest, emp, C_NULL)
        res != 0 && error("return code $res from JAMSD")
        res = ccall((:ctx_exportmodel, jamsd_libname), Cint, (Ptr{context}, Ptr{context}), ctx, ctx_dest)
        res != 0 && error("return code $res from JAMSD")

        JAMSDWriter.emp_hack(emp)

        if CONFIG[:export_gms]
            ccall((:ctx_writemodel, jamsd_libname), Cint, (Ptr{context}, Cstring), ctx_dest, "validation.gms")
            ccall((:ctx_setsolverstr, jamsd_libname), Cint, (Ptr{context}, Cstring), ctx_dest, "CONVERTD")
            ccall((:ctx_callsolver, jamsd_libname), Cint, (Ptr{context},), ctx_dest)
        end

        # switch back to the default solver
        ccall((:ctx_setsolverstr, jamsd_libname), Cint, (Ptr{JAMSDWriter.context}, Cstring), ctx_dest, solver_name)

        return ccall((:ctx_callsolver, jamsd_libname), Cint, (Ptr{context},), ctx_dest)
    end
end

function jamsd_setup_gams()
    ctx = ccall((:ctx_alloc, jamsd_libname), Ptr{context}, (Cuint,), 0)

    gamscntr_template_file = joinpath(solverdata_dir, "gamscntr.dat")

    if !isfile(gamscntr_template_file)
        jamsd_init_gams_solverdata()
    end
    if !isfile(gamscntr_template_file)
        error("Could not create template GAMS control file! Make sure that GAMS is properly installed and available via the system path")
    end

    gamscntr_template = readstring(gamscntr_template_file)

    gams_dir = "gams_dir"
    cur_dir = joinpath(pwd(), "gams_dir")
    if isdir(gams_dir)
        try
            rm("gams_dir", recursive=true, force=true)
            mkdir(gams_dir)
        catch
            run(`cmd /C RMDIR /s /q gams_dir`)
            gams_dir = mktempdir(pwd())
            cur_dir = gams_dir
        end
    else
        mkdir(gams_dir)
    end

    gamscntr_file = open(joinpath(gams_dir, "gamscntr.dat"), "w")

    println(gamscntr_file, replace(gamscntr_template, r"@@SUB@@", cur_dir))
    close(gamscntr_file)

    # we need an empty Matrixfile
    touch(joinpath(cur_dir, "gamsmatr.dat"))

    res = ccall((:gams_set_gamscntr, jamsd_libname), Cint, (Ptr{context}, Cstring), ctx, joinpath(cur_dir, "gamscntr.dat"))
    res != 0 && error("return code $res from JAMSD")

    # hm bad hack
    gamsdir = split(gamscntr_template, term_str)[29]

    CONFIG[:debug] && println("DEBUG: gamsdir is ``$gamsdir''")

    res = ccall((:gams_set_gamsdir, jamsd_libname), Cint, (Ptr{context}, Cstring), ctx, gamsdir)
    res != 0 && error("return code $res from JAMSD")

    ENV["PATH"] *= ":" * gamsdir

    return ctx
end

function jamsd_init_gams_solverdata()
    if !isdir(solverdata_dir)
        error("No directory named $solverdata_dir. ")
    end

    gamscntr = joinpath(solverdata_dir, "gamscntr.dat")
    if isfile(gamscntr)
        return
    end

    substr = joinpath(solverdata_dir, "tototututata")
    gms_file = joinpath(solverdata_dir, "dummy.gms")

    rm(substr, force=true, recursive=true)
    mkdir(substr)

    path = ENV["PATH"]
    println("$path")
    run(`gams $gms_file scrdir=$substr lo=0`)

    out_gamscntr = open(gamscntr, "w")
    input = readstring(open(joinpath(substr, "gamscntr.dat")))
    input = replace(input, pwd(), "@@SUB@@")
    println(out_gamscntr, replace(input, substr, "@@SUB@@"))

    close(out_gamscntr)

    rm(substr, recursive=true)
end
