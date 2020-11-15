using ObaraSaika
using Test
using PyCall
using QuadGK

import ObaraSaika: ERI, compute_integrals, boys, compute_normalisation

@testset "ObaraSaika.jl" begin
    @testset "Normalisation" begin
        f(x) = exp(-1.5 * abs(x)^2)
        integral = 4π * quadgk(x -> f(x) * f(x) * x * x, 0, 100)[1]
        @test compute_normalisation(1.5, [0]) ≈ 1 / sqrt(integral)
    end

    @testset "Boys" begin
        boystest(m, t) = quadgk(x -> x^(2m) * exp(-t * x^2), 0, 1)[1]
        @test boystest(0, 1.5)   ≈ boys(0, 1.5)
        @test boystest(0, 0.375) ≈ boys(0, 0.375)
        @test boystest(1, 1.5)   ≈ boys(1, 1.5)
        @test boystest(1, 0.375) ≈ boys(1, 0.375)
    end

    @testset "H 1s" begin
        bas = pyimport("pyscf.gto").parse("H S\n1.5 1.0\nH S\n0.9 1.0")
        mol  = pyimport("pyscf.gto").M(basis=bas, atom="H 0 0 0", spin=1, unit="Bohr")
        ref = mol.intor("int2e")

        #                    c    zeta  P         am
        basis_functions = [[(1.0, 1.5, [0, 0, 0], 0)],
                           [(1.0, 0.9, [0, 0, 0], 0)]]
        res = compute_integrals(ERI(), basis_functions)
        @test res ≈ ref
    end

    @testset "H₂ 1s" begin
        bas = pyimport("pyscf.gto").parse("H S\n1.0 1.0\n")
        mol  = pyimport("pyscf.gto").M(basis=bas, atom="H 0 0 0; H 0 1 0", unit="Bohr")
        ref = mol.intor("int2e")

        #                    c    zeta  P         am
        basis_functions = [[(1.0, 1.0, [0, 0, 0], 0)],
                           [(1.0, 1.0, [0, 1, 0], 0)]]
        res = compute_integrals(ERI(), basis_functions)
        @test res ≈ ref
    end

    @testset "H sto-3g" begin
        mol  = pyimport("pyscf.gto").M(basis="sto-3g", atom="H 0 0 0", spin=1, unit="Bohr")
        ref = mol.intor("int2e")

        #                    c           zeta         P         am
        basis_functions = [[(0.15432897, 3.42525091, [0, 0, 0], 0),
                            (0.53532814, 0.62391373, [0, 0, 0], 0),
                            (0.44463454, 0.16885540, [0, 0, 0], 0)]]
        res = compute_integrals(ERI(), basis_functions)
        @test res ≈ ref atol=1e-7  # TODO That's not a great agreement
    end

    @testset "H₂ sto-3g" begin
        mol  = pyimport("pyscf.gto").M(basis="sto-3g", atom="H 0 0 0; H 0 1 0", unit="Bohr")
        ref = mol.intor("int2e")

        #                    c           zeta         P         am
        basis_functions = [[(0.15432897, 3.42525091, [0, 0, 0], 0),
                            (0.53532814, 0.62391373, [0, 0, 0], 0),
                            (0.44463454, 0.16885540, [0, 0, 0], 0)],
                           [(0.15432897, 3.42525091, [0, 1, 0], 0),
                            (0.53532814, 0.62391373, [0, 1, 0], 0),
                            (0.44463454, 0.16885540, [0, 1, 0], 0)]]
        res = compute_integrals(ERI(), basis_functions)
        @test res ≈ ref atol=1e-7  # TODO That's not a great agreement
    end
end
