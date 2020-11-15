struct ERI end
# using Infiltrator

function compute_integrals(::ERI, basis_functions)
    T = Float64
    n_bas = length(basis_functions)

    # Build list of unique primitives
    primitives = [(ζ, P, am) for cgto in basis_functions
                  for (c, ζ, P, am) in cgto]
    primitives = unique(sort(primitives))

    # TODO compute product of contraction and normalisation here

    # Construct mapping primitives -> cgto
    cgtos = [Tuple{Int, T}[] for _ = 1:length(primitives)]
    for (ibas, cgto) in enumerate(basis_functions)
        for (c, ζ, P, am) in cgto
            idx = findfirst(isequal((ζ, P, am)), primitives)
            @assert sum(am) == 0  # Only s functions for now
            @assert !isnothing(idx)

            coefficient = c * compute_normalisation(ζ, am)
            push!(cgtos[idx], (ibas, coefficient))
        end
    end

    # Compute pairs of primitives
    pairs = []
    for (a, prim_a) in enumerate(primitives)
        for b = 1:length(primitives)  # use symmetry and start at a
            prim_b = primitives[b]

            # Compute ζab, Pab and Kab
            ζa, Pa, ama = prim_a
            ζb, Pb, amb = prim_b

            wζab = ζa * ζb / (ζa + ζb)
            ζab  = ζa + ζb
            Pab  = (ζa * Pa + ζb * Pb) / (ζa + ζb)
            Kab  = sqrt(2) * π^(5/4) / (ζa + ζb) * exp(-wζab * sum(abs2, Pa - Pb))

            # Collect contraction coefficients
            target_cgtos = []
            for (ibas_a, coeff_a)  in cgtos[a], (ibas_b, coeff_b) in cgtos[b]
                # TODO Could do some screening here
                push!(target_cgtos, (ibas_a, ibas_b, coeff_a * coeff_b))
            end

            if !isempty(target_cgtos)
                push!(pairs, (ζab, Pab, Kab, target_cgtos))
            end
        end
    end

    values = zeros(T, length(pairs), length(pairs))
    for (ab, pair_ab) in enumerate(pairs)
        for cd in 1:length(pairs)   # use symmetry and start at ab
            pair_cd = pairs[cd]

            m = 0
            coeffs = auxilary_coefficients(ERI(), pair_ab, pair_cd)
            vec = auxilary_vector(ERI(), m, pair_ab, pair_cd)

            values[ab, cd] = coeffs * vec
        end
    end

    # @infiltrate

    out = zeros(T, (n_bas, n_bas, n_bas, n_bas))
    for (ab, pair_ab) in enumerate(pairs), (cd, pair_cd) in enumerate(pairs)
        for (ibas_a, ibas_b, cab) in pair_ab[4], (ibas_c, ibas_d, ccd) in pair_cd[4]
            out[ibas_a, ibas_b, ibas_c, ibas_d] += values[ab, cd] * cab * ccd
        end
    end
    out
end


function compute_normalisation(ζ::T, am) where T <: Real
    dfacs = sqrt(prod(@. doublefactorial(max(1, 2 * am - 1))))
    (2ζ / T(π))^(3 / T(4)) * (4ζ)^(sum(am) / T(2)) / T(dfacs)
end

function auxilary_coefficients(::ERI, pair_ab, pair_cd)
    ζab, Pab, Kab = pair_ab
    ζcd, Pcd, Kcd = pair_cd
    Kab * Kcd / sqrt(ζab + ζcd)
end

function auxilary_vector(::ERI, m, pair_ab, pair_cd)
    ζab, Pab, Kab = pair_ab
    ζcd, Pcd, Kcd = pair_cd
    wζabcd = ζab * ζcd / (ζab + ζcd)
    boys(m, wζabcd * sum(abs2, Pab - Pcd))
end
