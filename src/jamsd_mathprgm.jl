function jamsd_set_modeltype(m::JAMSDMathProgModel, idx)
    discrete = any((m.vartypes .== :Int) + (m.vartypes .== :Bin) .> 0)
    if discrete
        if m.model_type == qcp
            m.model_type = miqcp
        elseif m.model_type == nlp
             m.model_type = minlp
        end
    end
    res = ccall((:ctx_setmodeltype, "libjamsd"), Cint, (Ptr{context}, Cint), m.jamsd_ctx, m.model_type)
    res != 0 && error("return code $res from JAMSD")
    res = ccall((:ctx_setobjsense, "libjamsd"), Cint, (Ptr{context}, Cint), m.jamsd_ctx, sense_to_jamsd[m.sense])
    res != 0 && error("return code $res from JAMSD")
end

function jamsd_declare_mathprgm(mp, ctx::Ptr{context}, id::Int)
    jamsd_mp = emp_mp_alloc(ctx, id)
    if mp.objequ > 0 || mp.objvar > 0
        typ = 0
        # mode is AGENT_OPT2MCP
        mod = 1
        emp_mp_start(jamsd_mp, typ, mod)
        emp_mp_objdir(jamsd_mp, mp.sense)
        emp_mp_objequ(jamsd_mp, mp.objequ-1)
        emp_mp_objvar(jamsd_mp, mp.objvar-1)
        for idx in mp.equs
            if idx != mp.objequ
                emp_mp_constraint(jamsd_mp, idx-1)
            end
        end
        for idx in mp.vars
            emp_mp_var(jamsd_mp, idx-1)
        end
    else
        typ = 2
        mod = 2
        emp_mp_start(jamsd_mp, typ, mod)
        VIvars = keys(mp.matching)
        equ_seen = Vector{Bool}(length(mp.equs))
        fill!(equ_seen, false)
        for vidx in mp.vars
            if vidx in VIvars
                eidx = mp.matching[vidx]
                emp_mp_vipair(jamsd_mp, eidx-1, vidx-1)
                equ_seen[eidx] = true;
                CONFIG[:debug] && println("DEBUG: eqn $eidx perp x[$vidx]")
            else
                # QVI case
                emp_mp_var(jamsd_mp, vidx-1)
            end
        end

        for idx in mp.equs
            if !equ_seen[idx]
                CONFIG[:debug] && println("DEBUG: eqn $eidx is a constraint")
                emp_mp_constraint(jamsd_mp, idx-1)
            end
        end
    end

    return jamsd_mp
end
