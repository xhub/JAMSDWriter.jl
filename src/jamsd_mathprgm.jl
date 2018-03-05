function jamsd_set_modeltype(m::JAMSDMathProgModel, idx)
    discrete = any((m.vartypes .== :Int) + (m.vartypes .== :Bin) .> 0)
    if discrete
        if m.model_type == qcp
            m.model_type = miqcp
        elseif m.model_type == nlp
             m.model_type = minlp
        end
    end
    res = ccall((:ctx_setmodeltype, jamsd_libname), Cint, (Ptr{context}, Cint), m.jamsd_ctx, m.model_type)
    res != 0 && error("return code $res from JAMSD")
    res = ccall((:ctx_setobjsense, jamsd_libname), Cint, (Ptr{context}, Cint), m.jamsd_ctx, sense_to_jamsd[m.sense])
    res != 0 && error("return code $res from JAMSD")
end

function jamsd_declare_mathprgm(mp, ctx::Ptr{context}, id::Int)
    m = mp.emp.m.internalModel.inner

    jamsd_mp = emp_mp_alloc(ctx, id)
    if mp.objequ > 0 || mp.objvar > 0
        typ = 0
        # mode is AGENT_OPT2MCP
        mod = 1

        emp_mp_start(jamsd_mp, typ, mod)
        emp_mp_objdir(jamsd_mp, mp.sense)
        emp_mp_objequ(jamsd_mp, mp.objequ-1)

        if mp.objvar > 0
            emp_mp_objvar(jamsd_mp, m.v_index_map[mp.objvar])
        else
            emp_mp_objvar(jamsd_mp, mp.objvar-1)
        end

        for eidx in mp.equs
            if eidx != mp.objequ
                emp_mp_constraint(jamsd_mp, eidx-1)
            end
        end

        for vidx in mp.vars
            emp_mp_var(jamsd_mp, m.v_index_map[vidx])
        end
    else
        typ = 2
        mod = 2
        emp_mp_start(jamsd_mp, typ, mod)
        VIvars = keys(mp.matching)
        equ_seen = Vector{Int}(length(mp.equs))
        sidx = 1
        for vidx in mp.vars
            if vidx in VIvars
                eidx = mp.matching[vidx]
                emp_mp_vipair(jamsd_mp, eidx-1, m.v_index_map[vidx])
                equ_seen[sidx] = eidx
                sidx += 1
                CONFIG[:debug] && println("DEBUG: eqn $eidx perp x[$vidx]")
            else
                # QVI case
                emp_mp_var(jamsd_mp, vidx-1)
            end
        end

        for idx in mp.equs
            if !(idx in equ_seen)
                CONFIG[:debug] && println("DEBUG: eqn $eidx is a constraint")
                emp_mp_constraint(jamsd_mp, idx-1)
            end
        end
    end

    return jamsd_mp
end
