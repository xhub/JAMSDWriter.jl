mutable struct LinearityExpr
    c
    linearity::Symbol
    coeff::Float64
end

LinearityExpr(c, linearity) = LinearityExpr(c, linearity, 1)

eval(c::LinearityExpr) = eval(get_expr(c))
Base.print(io::IO, c::LinearityExpr) = print(io::IO, "($(c.c),$(c.linearity))")
Base.show(io::IO, c::LinearityExpr) = print(io::IO, "($(c.c),$(c.linearity))")

LinearityExpr(c) = throw("Couldn't determine linearity of expression:\n$c")
LinearityExpr(c::Real) = LinearityExpr(c, :const)
LinearityExpr(c::LinearityExpr) = c  # Already processed; ignore!
function LinearityExpr(c::Expr)
    if c.head == :call
        for i = 2:length(c.args)
            c.args[i] = LinearityExpr(c.args[i])
        end

        args = c.args[2:end]
        if c.args[1] in [:+, :-]
            if check_for_linearity(:nonlinear, args)
                linearity = :nonlinear
            elseif check_for_linearity(:linear, args)
                linearity = :linear
            else
                linearity = :const
            end
        elseif c.args[1] == :*
            if check_for_linearity(:nonlinear, args) ||
               length(filter((c) -> check_linearity(:linear, c), args)) > 1
                linearity = :nonlinear
            elseif check_for_linearity(:linear, args)
                linearity = :linear
            else
                linearity = :const
            end
        elseif c.args[1] == :/
            if c.args[3].linearity != :const
                linearity = :nonlinear
            else
                linearity = c.args[2].linearity
            end
        elseif c.args[1] == :^
            if c.args[3].linearity != :const
                linearity = :nonlinear
            elseif c.args[2].linearity == :linear
                if c.args[3].c == 1
                    linearity = :linear
                else
                    linearity = :nonlinear
                end
            else
                linearity = c.args[2].linearity
            end
        elseif c.args[1] == :ifelse
            if c.args[2].linearity == :const
                # We know which branch to take already - simplify
                c.args[2] = pull_up_constants(c.args[2])
                if c.args[2].c != 0
                    linearity = c.args[3].linearity
                    c = c.args[3].c
                else
                    linearity = c.args[4].linearity
                    c = c.args[4].c
                end
            else
                linearity = :nonlinear
            end
        else
            if check_for_linearity(:linear, args) ||
               check_for_linearity(:nonlinear, args)
                linearity = :nonlinear
            else
                linearity = :const
            end
        end

    elseif c.head == :ref
        linearity = :linear
    elseif c.head == :comparison
        indices = 1:2:length(c.args)
        for i in indices
            c.args[i] = LinearityExpr(c.args[i])
        end
        if check_for_linearity(:linear, c.args[indices]) ||
           check_for_linearity(:nonlinear, c.args[indices])
            linearity = :nonlinear
        else
            linearity = :const
        end
    elseif c.head in [:&&, :||]
        for i = 1:length(c.args)
            c.args[i] = LinearityExpr(c.args[i])
        end
        if check_for_linearity(:linear, c.args) ||
           check_for_linearity(:nonlinear, c.args)
            linearity = :nonlinear
        else
            linearity = :const
        end
    end
    return LinearityExpr(c, linearity)
end

check_linearity(linearity::Symbol, c::LinearityExpr) = c.linearity == linearity

function check_for_linearity(linearity::Symbol, args::Array)
    return !isempty(filter((c) -> check_linearity(linearity, c), args))
end

get_expr(c) = c
get_expr(c::LinearityExpr) = get_expr(c.c)
function get_expr(c::Expr)
    map!(get_expr, c.args, c.args)
    return c
end

pull_up_constants(c) = c
function pull_up_constants(c::LinearityExpr)
    if c.linearity == :const
        c.c = eval(c)
    elseif isa(c.c, Expr)
        map!(pull_up_constants, c.c.args, c.c.args)
    end
    return c
end

function prune_linear_terms!(c::LinearityExpr, lin_constr::Dict{Int, Float64},
                             constant::Float64=0.0, negative_tree::Bool=false)
    if c.linearity != :nonlinear
        constant = add_linear_tree!(c, lin_constr, constant, negative_tree)
        c = LinearityExpr(:(0), :const)
        return true, c, constant
    else
        expr = c.c
        if expr.head == :call
            if expr.args[1] == :+
                n = length(expr.args)
                pruned = falses(n - 1)
                for i in 2:n
                    pruned[i - 1], expr.args[i], constant = prune_linear_terms!(
                        expr.args[i], lin_constr, constant, negative_tree)
                end
                if sum(.!pruned) > 1
                    inds = vcat([1], collect(2:n)[.!pruned])
                    c.c.args = expr.args[inds]
                else
                    c = expr.args[findfirst(.!pruned) + 1]
                end
            elseif expr.args[1] == :-
                if length(expr.args) == 3
                    pruned_first, expr.args[2], constant = prune_linear_terms!(
                        expr.args[2], lin_constr, constant, negative_tree)
                    pruned_second, expr.args[3], constant = prune_linear_terms!(
                        expr.args[3], lin_constr, constant, !negative_tree)
                    if pruned_first
                        new_expr = Expr(:call, :-, expr.args[3])
                        c = LinearityExpr(new_expr, :nonlinear)
                    elseif pruned_second
                        c = expr.args[2]
                    end
                end
            end
        end
        return false, c, constant
    end
end

function add_linear_tree!(c::LinearityExpr, lin_constr::Dict{Int, Float64},
                          constant::Float64=0.0, negative_tree::Bool=false)
    c = collate_linear_terms(c)
    negative_tree && negate(c)
    constant = add_tree_to_constr!(c, lin_constr, constant)
    return constant
end

function collate_linear_terms(c::LinearityExpr)
    if isa(c.c, Expr) && c.c.head == :call
        for i in 2:length(c.c.args)
            c.c.args[i] = collate_linear_terms(c.c.args[i])
        end
    end
    if c.linearity == :linear && isa(c.c, Expr) && c.c.head == :call
        func = c.c.args[1]
        if func == :-
            if length(c.c.args) == 2
                c = negate(c.c.args[2])
            else
                c.c.args[3] = negate(c.c.args[3])
                c.c.args[1] = :+
            end
        elseif func == :*
            # This can be n-ary multiplication, but there is only one variable
            # First, find index of variable in the args
            x = findfirst(ex -> ex.linearity == :linear, c.c.args[2:end]) + 1
            # Multiply variable term by each of the constants in the args
            x_expr = c.c.args[x]
            for i in setdiff(2:length(c.c.args), x)
                x_expr = multiply(x_expr, c.c.args[i].c)
            end
            c = x_expr
        elseif func == :/
            c = multiply(c.c.args[2], 1 / c.c.args[3].c)
        elseif func == :^
            c = c.c.args[2]
        end
    end
    return c
end

negate(c::LinearityExpr) = multiply(c::LinearityExpr, -1)
function multiply(c::LinearityExpr, a::Real)
    if isa(c.c, Expr) && c.c.head == :call
        @assert c.c.args[1] == :+
        map!(arg -> multiply(arg, a), c.c.args[2:end], c.c.args[2:end])
    elseif c.linearity == :const
        c.c *= a
    else
        c.coeff *= a
    end
    return c
end

function add_tree_to_constr!(c::LinearityExpr, lin_constr::Dict{Int, Float64},
                             constant::Float64)
    if isa(c.c, Expr) && c.c.head == :call
        @assert c.c.args[1] == :+
        for i in 2:length(c.c.args)
            constant = add_tree_to_constr!(c.c.args[i], lin_constr, constant)
        end
    elseif c.linearity == :const
        constant += c.c
    else
        lin_constr[c.c.args[2]] += c.coeff
    end
    return constant
end
