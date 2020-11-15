struct ERI end

function compute_integrals(::ERI, basis_functions)
    T = Float64

    # Build list of unique primitives
    primitives = [(ζ, P, λ) for cgto in basis_functions
                  for (c, ζ, P, λ) in cgto]
    primitives = unique(sort(primitives, by=reverse))

    # TODO compute product of contraction and normalisation here

    # Construct mapping primitives -> cgto
    cgtos = [Tuple{Int, T}[] for _ = 1:length(primitives)]
    for (ibas, cgto) in enumerate(basis_functions)
        for (c, ζ, P, λ) in cgto
            idx = findfirst(isequal((ζ, P, λ)), primitives)
            @assert !isnothing(idx)

            coefficient = c * compute_normalisation(ζ, λ)
            push!(cgtos[idx], (ibas, coefficient))
        end
    end

    # Compute pairs of primitives
    pairs = []
    for (a, prim_a) in enumerate(primitives)
        for b = 1:a
            prim_b = primitives[b]

            # Compute ζab, Pab and Kab
            ζa, Pa, λa = prim_a
            ζb, Pb, λb = prim_b
            @assert λa ≥ λb  # Otherwise recursion not implemented

            wζab = ζa * ζb / (ζa + ζb)
            ζab  = ζa + ζb  # ζ / η
            Pab  = (ζa * Pa + ζb * Pb) / (ζa + ζb)  # P / Q
            Kab  = sqrt(2) * π^(5/4) / (ζa + ζb) * exp(-wζab * sum(abs2, Pa - Pb))

            # Collect contraction coefficients
            target_cgtos = []
            for (ibas_a, coeff_a)  in cgtos[a], (ibas_b, coeff_b) in cgtos[b]
                # TODO Could do some screening here

                # The (a, b) double loop only covers the lower triangle, since
                # <a |Op b> = <b | Op a>. But to get the integral of the contracted
                # cGTO containing the primitives a and b, both the terms <b | Op a>
                # and <a |Op b> are needed in the sum over (product of pairs of) 
                # contraction coefficients. To account for this the contraction
                # coefficient product is augmented by a factor of two in such cases.
                symfac = (a != b && ibas_a == ibas_b) ? 2 : 1
                push!(target_cgtos, (ibas_a, ibas_b, symfac * coeff_a * coeff_b))
            end

            if !isempty(target_cgtos)
                push!(pairs, (ζab, Pab, Kab, target_cgtos, (λa, λb), (Pa, Pb)))
            end
        end
    end

    # TODO This is super horrible ... probably need a better data structure
    # Total number of cartesian cgto basis functions
    # fun[1][4] == λ of the first contracted GTO, assumed to be the same in all
    # primitives (which is true, but should have a better interface)
    n_s = count(fun -> fun[1][4] == 0, basis_functions)  # Total number of s functions
    n_p = count(fun -> fun[1][4] == 1, basis_functions)  # Total number of p functions
    n_bas = sum(basis_functions) do fun
        c, ζ, P, λ = fun[1]  # λ of first primitive == λ of full cgto
        @assert λ ≤ 1
        3^λ
    end

    n_p > 0 && error("wrong here!")
    # FIXME This uses the (incorrect) assumption that p functions
    #       always come *after* the s functions ... this is only true for the
    #       primitives!
    function cgto_range(λ, ibas)
        λ == 0 && return ibas:ibas
        inp = (ibas - n_s - 1)
        λ == 1 && return n_s .+ 3inp .+ (ibas:ibas+3)
    end

    out = zeros(T, (n_bas, n_bas, n_bas, n_bas))
    for (ab, pair_ab) in enumerate(pairs)
        for cd in 1:ab
            pair_cd = pairs[cd]

            # Compute stuff we need about the four centres to drive the recursion
            ζab, Pab, Kab, target_ab, (λa, λb), centres_ab = pair_ab
            ζcd, Pcd, Kcd, target_cd, (λc, λd), centres_cd = pair_cd
            @assert (λa, λb) ≥ (λc, λd)  # Otherwise recursion not implemented
            ζabcd = ζab + ζcd  # ζ + η
            wζabcd = ζab * ζcd / (ζab + ζcd)  # ρ
            Pabcd = (ζab * Pab + ζcd * Pcd) / (ζab + ζcd)  # W
            λabcd = (λa, λb, λc, λd)  # AM on each centre
            centres = (centres_ab..., centres_cd...) # A, B, C, D
            quad = (ζabcd, Pabcd, wζabcd, (ζab, ζcd),
                    (Pab, Pcd), (Kab, Kcd), λabcd, centres)

            # TODO This obviously needs to be generalised for higher AM functions
            #      Is it true that the symmetrised loop (0,0,0,0) to λabcd will always cut it here?
            interm = Dict{NTuple{4, Int}, Array{T, 5}}()  # Intermediates
            for λtest in ((0, 0, 0, 0), (1, 0, 0, 0), (1, 0, 1, 0),
                          (1, 1, 0, 0), (1, 1, 1, 0), (1, 1, 1, 1))
                all(λabcd .≥ λtest) || break
                interm[λtest] = compute_intermediate(ERI(), Val(λtest), quad, interm)
            end

            coeffs = interm[λabcd]
            m_max = size(coeffs, 5) - 1
            auxvec = auxilary_vector.(Ref(ERI()), 0:m_max, Ref(quad))
            auxvec = reshape(auxvec, 1, 1, 1, 1, m_max + 1)
            res = dropdims(sum(coeffs .* auxvec, dims=5), dims=5)

            # Store results
            # TODO Probably later want to call some callback or otherwise
            #      consume the data
            for (ibas_a, ibas_b, cab) in target_ab
                for (ibas_c, ibas_d, ccd) in target_cd
                    rnge_a = cgto_range(λa, ibas_a)
                    rnge_b = cgto_range(λb, ibas_b)
                    rnge_c = cgto_range(λc, ibas_c)
                    rnge_d = cgto_range(λd, ibas_d)

                    # Same as with forming the pairs, the loop over quads takes
                    # symmetry into account and thus again suppresses some of the
                    # terms in the contraction sum, for which a factor of two is needed.
                    symfac = ((ab != cd) && (ibas_a, ibas_b) == (ibas_c, ibas_d)) ? 2 : 1

                    # Multiply in factors (1/2) which will be undone when undoing
                    # the symmetry reduction below
                    symfac *= (  (ibas_a == ibas_b ? 1/2 : 1)
                               * (ibas_c == ibas_d ? 1/2 : 1)
                               * ((ibas_a, ibas_b) == (ibas_c, ibas_d) ? 1/2 : 1))


                    out[rnge_a, rnge_b, rnge_c, rnge_d] .+= symfac .* cab .* ccd .* res
                end  # target_cd
            end  # target_ab
        end
    end

    # Undo the symmetry reduction
    out .= out .+ permutedims(out, [3, 4, 1, 2])
    out .= out .+ permutedims(out, [2, 1, 3, 4])
    out .= out .+ permutedims(out, [1, 2, 4, 3])

    out
end

function compute_normalisation(ζ::T, λ) where T <: Real
    dfacs = sqrt(prod(@. doublefactorial(max(1, 2 * λ - 1))))
    (2ζ / T(π))^(3 / T(4)) * (4ζ)^(sum(λ) / T(2)) / T(dfacs)
end

# TODO The formulas for the p functions
#      have not been checked properly and might be wrong

# Compute (ps, ss) intermediate
function compute_intermediate(::ERI, ::Val{(1, 0, 0, 0)}, quad, interm::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    #                                     Note for future generalisation:
    Pa    = reshape(Pa,    3, 1, 1, 1)  # These have the shape of the centre
    Pab   = reshape(Pab,   3, 1, 1, 1)  # to be lowered and the data comes
    Pabcd = reshape(Pabcd, 3, 1, 1, 1)  # also from the centre to be lowered

    #           i  j  k  l  m
    out = zeros(3, 1, 1, 1, 2)
    # A bit annoying here is that out[:, :, :, :, 1] is actually m=0
    @. begin
        # (psss) coefficient for m=0
        out[:, :, :, :, 1:1] = (Pab   - Pa ) * interm[(0, 0, 0, 0)]
        # (psss) coefficient for m=1
        out[:, :, :, :, 2:2] = (Pabcd - Pab) * interm[(0, 0, 0, 0)]
    end
    return out
end

# Compute (ps, ps) intermediate
function compute_intermediate(::ERI, ::Val{(1, 0, 1, 0)}, quad, interm::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pc    = reshape(Pc,    1, 1, 3, 1)
    Pcd   = reshape(Pcd,   1, 1, 3, 1)
    Pabcd = reshape(Pabcd, 1, 1, 3, 1)

    #           i  j  k  l  m
    out = zeros(3, 1, 3, 1, 3)
    @. begin
        out[:, :, :, :, 1:2] += (Pcd   - Pc ) * interm[(1, 0, 0, 0)]
        # (psss)^(0) has a coefficient for m=0 and m=1, but we want (psss)^(1),
        # so just shift the m-range where we set the values one up.
        out[:, :, :, :, 2:3] += (Pabcd - Pcd) * interm[(1, 0, 0, 0)]

        for i in 1:3
            out[i, :, i, :, 2] += interm[(0, 0, 0, 0)] / 2ζabcd
        end
    end
    out
end

# Compute (pp, ss) intermediate
function compute_intermediate(::ERI, ::Val{(1, 1, 0, 0)}, quad, interm::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pb    = reshape(Pb,    1, 3, 1, 1)
    Pab   = reshape(Pab,   1, 3, 1, 1)
    Pabcd = reshape(Pabcd, 1, 3, 1, 1)

    #           i  j  k  l  m
    out = zeros(3, 3, 1, 1, 3)
    @. begin
        out[:, :, :, :, 1:2] += (Pab   - Pb ) * interm[(1, 0, 0, 0)]
        out[:, :, :, :, 2:3] += (Pabcd - Pab) * interm[(1, 0, 0, 0)]

        for i in 1:3
            out[i, i, :, :, 1] += interm[(0, 0, 0, 0)] / 2ζab
            out[i, i, :, :, 2] -= interm[(0, 0, 0, 0)] * wζabcd / ζab
        end
    end
    out
end

# Compute (pp, ps) intermediate
function compute_intermediate(::ERI, ::Val{(1, 1, 1, 0)}, quad, interm::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pc    = reshape(Pc,    1, 1, 3, 1)
    Pcd   = reshape(Pcd,   1, 1, 3, 1)
    Pabcd = reshape(Pabcd, 1, 1, 3, 1)

    #           i  j  k  l  m
    out = zeros(3, 3, 3, 1, 4)
    @. begin
        out[:, :, :, :, 1:3] += (Pcd   - Pc ) * interm[(1, 1, 0, 0)]
        out[:, :, :, :, 2:4] += (Pabcd - Pcd) * interm[(1, 1, 0, 0)]

        for i in 1:3
            out[i, :, i, :, 2:3] += permutedims(interm[(1, 0, 0, 0)], [2, 1, 3, 4, 5]) / 2ζabcd
            out[:, i, i, :, 2:3] += interm[(1, 0, 0, 0)] / 2ζabcd
        end
    end
    out
end

# Compute (pp, pp) intermediate
function compute_intermediate(::ERI, ::Val{(1, 1, 1, 1)}, quad, interm::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pd    = reshape(Pd,    1, 1, 1, 3)
    Pcd   = reshape(Pcd,   1, 1, 1, 3)
    Pabcd = reshape(Pabcd, 1, 1, 1, 3)

    #           i  j  k  l  m
    out = zeros(3, 3, 3, 3, 5)
    @. begin
        out[:, :, :, :, 1:4] += (Pcd   - Pd ) * interm[(1, 1, 1, 0)]
        out[:, :, :, :, 2:5] += (Pabcd - Pcd) * interm[(1, 1, 1, 0)]

        for i in 1:3
            out[i, :, :, i, 2:4] += permutedims(interm[(1, 0, 1, 0)], [2, 1, 3, 4, 5]) / 2ζabcd
            out[:, i, :, i, 2:4] += interm[(1, 0, 1, 0)] / 2ζabcd
            out[:, :, i, i, 1:3] += interm[(1, 1, 0, 0)] / 2ζcd
            out[:, :, i, i, 2:5] -= interm[(1, 1, 0, 0)] * wζabcd / ζcd
        end
    end
    out
end

# The auxiliary integral
function compute_intermediate(::ERI, ::Val{(0, 0, 0, 0)}, quad, ::Dict)
    ζabcd = quad[1]  # ζ + η
    Kab, Kcd = quad[6]
    reshape([Kab * Kcd / sqrt(ζabcd)], 1, 1, 1, 1, 1)
end

function auxilary_vector(::ERI, m, quad)
    wζabcd = quad[3]  # ρ
    Pab, Pcd = quad[5]
    boys(m, wζabcd * sum(abs2, Pab - Pcd))
end
