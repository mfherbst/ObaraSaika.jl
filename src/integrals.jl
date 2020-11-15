struct ERI end

function compute_integrals(::ERI, basis_functions)
    T = Float64

    # Build list of unique primitives
    primitives = [(ζ, P, λ) for cgto in basis_functions
                  for (c, ζ, P, λ) in cgto]
    primitives = unique(sort(primitives))

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
    n_s = findall(fun -> fun[1][4] == 0, basis_functions)  # Total number of s functions
    n_p = findall(fun -> fun[1][4] == 1, basis_functions)  # Total number of p functions
    n_bas = sum(basis_functions) do fun
        c, ζ, P, λ = fun[1]  # λ of first primitive == λ of full cgto
        @assert λ ≤ 1
        3^λ
    end

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

            results = Dict{NTuple{4, Int}, Array{T, 5}}()

            results[(0, 0, 0, 0)] = compute_intermediate(ERI(), Val((0, 0, 0, 0)),
                                                         quad, results)

            # First do the s integrals
            # Then do the (ps | ss) integrals (if needed for current batch)
            # Then do the (ps | ps) integrals (if needed for current batch)
            # Then do the (pp | ss) integrals (if needed for current batch)
            # Then do the (pp | ps) integrals (if needed for current batch)
            # Then do the (pp | pp) integrals (if needed for current batch)

            coeffs = results[(0, 0, 0, 0)]
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

# Compute (ps, ss) intermediate
function compute_intermediate(::ERI, ::Val{(1, 0, 0, 0)}, quad, results::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pa    = reshape(Pa,    3, 1, 1, 1)
    Pab   = reshape(Pab,   3, 1, 1, 1)
    Pabcd = reshape(Pabcd, 3, 1, 1, 1)

    #           i  j  k  l  m
    out = zeros(3, 1, 1, 1, 2)
    @. begin
        out[:, :, :, :, 1] = (Pab   - Pa ) * results[(0, 0, 0, 0)]
        out[:, :, :, :, 2] = (Pabcd - Pab) * results[(0, 0, 0, 0)]
    end
    return out
end

# Compute (ps, ps) intermediate
function compute_intermediate(::ERI, ::Val{(1, 0, 1, 0)}, quad, results::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pc    = reshape(Pc,    1, 1, 3, 1)
    Pcd   = reshape(Pcd,   1, 1, 3, 1)
    Pabcd = reshape(Pabcd, 1, 1, 3, 1)

    #           i  j  k  l  m
    out = zeros(3, 1, 3, 1, 3)
    @. begin
        out[:, :, :, :, 1:2] += (Pcd   - Pc ) * results[(1, 0, 0, 0)]
        out[:, :, :, :, 2:3] += (Pabcd - Pcd) * results[(1, 0, 0, 0)]

        for i in 1:3
            out[i, :, i, :, 2] += results[(0, 0, 0, 0)] / 2ζabcd
        end
    end
    out
end

# Compute (pp, ss) intermediate
function compute_intermediate(::ERI, ::Val{(1, 1, 0, 0)}, quad, results::Dict)
    ζabcd, Pabcd, wζabcd, (ζab, ζcd), (Pab, Pcd), (Kab, Kcd), λabcd, centres = quad
    Pa, Pb, Pc, Pd = centres

    Pb    = reshape(Pb,    1, 3, 1, 1)
    Pab   = reshape(Pab,   1, 3, 1, 1)
    Pabcd = reshape(Pabcd, 1, 3, 1, 1)

    #           i  j  k  l  m
    out = zeros(3, 3, 1, 1, 3)
    @. begin
        out[:, :, :, :, 1:2] += (Pab   - Pb ) * results[(1, 0, 0, 0)]
        out[:, :, :, :, 2:3] += (Pabcd - Pab) * results[(1, 0, 0, 0)]

        for i in 1:3
            out[i, :, i, :, 1] += results[(0, 0, 0, 0)] / 2ζabcd
            out[i, :, i, :, 2] -= results[(0, 0, 0, 0)] * wζabcd / ζabcd
        end
    end
    out
end

# ppps
# pppp



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
