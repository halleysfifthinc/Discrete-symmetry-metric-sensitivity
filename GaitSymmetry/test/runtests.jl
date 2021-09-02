using GaitSymmetry
using Test, InteractiveUtils

@testset "GaitSymmetry.jl" begin
    @testset "inverses" begin
        x = abs(randn()) + 1
        y = abs(randn()) + 1

        for F in subtypes(SymmetryFunc)
            s = F(x,y)
            @test F(inv(F)(s)...) ≈ s
        end

        @test (inv(Plo05)(0.6931471805599453; scaled=false) .≈ (1., 2.)) |> prod
        @test_throws DomainError inv(Plo05)(-10.)
    end

    @testset "conversions" begin
        x = abs(randn()) + 1
        y = abs(randn()) + 1

        for F1 in subtypes(SymmetryFunc), F2 in subtypes(SymmetryFunc)
            s1 = F1(x,y)
            # Plo05 and Allen2011 is non-directional, so any conversion to another
            # metric must only consider magnitude
            if F1 ∈ (Plo05,)
                # Sel86 needs inverses (not magnitude) to remove consideration of
                # directionality
                if F2 === Sel86
                    @test (convert(F2, F1)(s1) ≈ inv(F2(x,y)) ||
                           convert(F2, F1)(s1) ≈ F2(x,y))
                else
                    @test abs(convert(F2, F1)(s1)) ≈ abs(F2(x,y))
                end
            else
                @test convert(F2, F1)(s1) ≈ F2(x,y)
            end
        end
    end
end
