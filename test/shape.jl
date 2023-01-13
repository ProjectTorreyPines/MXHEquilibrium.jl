R0 = 6.2
Z0 = 0.0
a = 2.0
ϵ = a/R0
κ = 1.7
δ = 0.33
δ₀ = asin(δ)

S = MillerShape(R0,Z0,ϵ,κ,δ)

@testset "Shape Tests" begin

    @testset "Promote and Convert Tests" begin
        S32 = MillerShape{Float32}(S)
        @test typeof(S32) == MillerShape{Float32}
        SS, SS32, v = promote(S, S32, 1.0f0)
        @test typeof(SS) == typeof(S)
        @test typeof(SS32) == typeof(S)
        @test typeof(v) == Float64
        @test convert(MShape{Float64}, 1.f0) == 1e0
    end

    @testset "Point Tests" begin
        @test (true,true) == (S(0.0)./R0  .≈ (1 + ϵ, 0.0))
        @test (true,true) == (S(pi)./R0   .≈ (1 - ϵ, 0.0))
        @test (true,true) == (S(pi/2)./R0 .≈ (1 - δ*ϵ, κ*ϵ))
    end

    @testset "Boundary Tests" begin
        @test round(circumference(S)/R0,digits=2) == 2.79
    end

    @testset "Curvature Tests" begin
        @test curvature(S, 0.0)*R0  ≈ (1 + δ₀)^2/(ϵ*κ^2)
        @test curvature(S, pi)*R0   ≈ (1 - δ₀)^2/(ϵ*κ^2)
        @test curvature(S, pi/2)*R0 ≈ κ/(ϵ*cos(δ₀)^2)
    end

    @testset "MXH Fitting" begin
        R0 = 1.6667475890195798
        Z0 = -0.14918543837718545
        ϵ = 0.3278211449811159
        κ = 1.8166930389369045
        c0 = -0.0847654065659985
        c = [-0.12517240241115193, 0.04846721164851622, 0.022137477337172622,
             -0.02377160268474721, -0.008327873921997632, 0.004168909440627845,
             -0.001275800378739416, -0.005074242254891684, -0.004038080173331539,
             0.0008125282874464121]
        s = [0.49443563266437024, 0.0640865222160445, -0.06090895885121784,
             0.006162625943678958, 0.017797802404712248, -0.0021482577173381435,
             -0.004270022610390773, -0.0023722092647870804, 0.0009513715077412459,
             -0.0009515964309653494]

        MXH = MXHShape(R0,Z0,ϵ,κ,c0,c,s)
        bdry = plasma_boundary(MXH,N=500)
        fMXH = fit(bdry,MXHShape(10))
        for f in fieldnames(MXHShape)
            (f in (:c, :s)) && continue
            @test getfield(MXH,f) ≈ getfield(fMXH,f) rtol = 1e-1
        end
        @test isapprox(MXH.c, fMXH.c, rtol=1e-1)
        @test isapprox(MXH.s, fMXH.s, rtol=1e-1)
    end

    @testset "Integration" begin
        g = readg((@__DIR__)*"/g150219.03200")
        M = efit(g,clockwise_phi=false)
        efit_bdry = plasma_boundary(M)
        mxh = fit(efit_bdry,MXHShape(5))
        mxh_bdry = plasma_boundary(mxh,N=1000)
        @test circumference(mxh,rtol=1e-6) ≈ circumference(mxh_bdry) atol=1e-3
        @test area(mxh,rtol=1e-6) ≈ area(mxh_bdry) atol=1e-3
        @test surface_area(mxh,rtol=1e-6) ≈ surface_area(mxh_bdry) atol=1e-3
        @test volume(mxh,rtol=1e-6) ≈ volume(mxh_bdry,dx=0.001,dy=0.001) atol=1e-3
    end
end
