# ITER parameters taken from
# “One size fits all” analytic solutions to the Grad–Shafranov equation
# Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818
R0 = 6.2
Z0 = 0.0
a = 2.0
B0 = 5.3
ϵ = a / R0
κ = 1.7
δ = 0.33
MS = MillerShape(R0, Z0, ϵ, κ, δ)
alpha = -0.155
qstar = 1.57
β_t = 5.02

S1 = solovev(B0, MS, alpha, qstar; B0_dir=1, Ip_dir=1)
S2 = solovev(B0, MS, alpha, qstar; B0_dir=-1, Ip_dir=1)
S3 = solovev(B0, MS, alpha, qstar; B0_dir=1, Ip_dir=-1)
S4 = solovev(B0, MS, alpha, qstar; B0_dir=-1, Ip_dir=-1)

btip1 = (1, 1)
btip2 = (-1, 1)
btip3 = (1, -1)
btip4 = (-1, -1)

cc1 = 3
cc2 = 3
cc3 = 3
cc4 = 3

r1 = magnetic_axis(S1) .+ (0.1, 0.0)
r2 = magnetic_axis(S2) .+ (0.1, 0.0)
r3 = magnetic_axis(S3) .+ (0.1, 0.0)
r4 = magnetic_axis(S4) .+ (0.1, 0.0)

test_data = ((cc1, S1, r1, btip1), (cc2, S2, r2, btip2), (cc3, S3, r3, btip3), (cc4, S4, r4, btip4))

@testset "Solov'ev Tests" begin

    @testset "Promote & Convert Tests" begin
        S1F32 = convert(SolovevEquilibrium{Float32,7,Float32}, S1)
        @test typeof(S1F32) == SolovevEquilibrium{Float32,7,Float32}
        SS, SS1, v = promote(S1, S1F32, 1.0f0)
        @test typeof(SS) == typeof(S1)
        @test typeof(SS1) == typeof(S1)
        @test typeof(v) == Float64
        @test convert(typeof(S1), 1) == one(Float64)
    end

    @testset "Beta Tests" begin
        @test round(beta_t(S1); digits=2) == β_t
        @test round(beta_p(S1); digits=2) == round(beta_t(S1) * qstar^2 / ϵ^2; digits=2)
    end

    @testset "Bt-Ip Signs" begin
        for (cc, S, r, btip) in test_data
            Bt, Ip = (sign(Bfield(S, r...)[2]), sign(Jfield(S, r...)[2]))
            @test (Bt, Ip) == btip
        end
    end

    @testset verbose = true "CurlB == μ₀J" begin
        for (cc, S, r, btip) in test_data
            @test curlB(S, r...) ≈ mu0 * Jfield(S, r...) rtol = 0.01
            @test MXHEquilibrium.curlB_autodiff(S, r...) ≈ mu0 * Jfield(S, r...) rtol = 0.01
        end
    end

    @testset verbose = true "3D tests" begin
        for (cc, S, r, btip) in test_data
            @test Bfield(S, r...) ≈ Bfield(S, r[1], zero(r[1]), r[2]) rtol = 0.01
            @test Jfield(S, r...) ≈ Jfield(S, r[1], zero(r[1]), r[2]) rtol = 0.01

            @test gradB(S, r...) ≈ gradB(S, r[1], zero(r[1]), r[2]) rtol = 0.01
            @test gradB(S, r[1], zero(r[1]), r[2]) ≈ MXHEquilibrium.gradB_autodiff(S, r[1], zero(r[1]), r[2]) rtol = 0.01

            @test curlB(S, r...) ≈ curlB(S, r[1], zero(r[1]), r[2]) rtol = 0.01
            @test curlB(S, r[1], zero(r[1]), r[2]) ≈ MXHEquilibrium.curlB_autodiff(S, r[1], zero(r[1]), r[2]) rtol = 0.01
        end
    end

    @testset verbose = true "COCOS Consistancy" begin
        for (cc, S, r, btip) in test_data
            psi = S(r...)
            sigma_F = sign(poloidal_current(S, psi))
            sigma_pprime = sign(pressure_gradient(S, psi))
            sigma_q = sign(safety_factor(S, psi))
            sigma_dpsi = -sign(psi) # Solov'ev psi at boundary zero
            @test check_cocos(btip[1], btip[2], sigma_F, sigma_pprime, sigma_q, sigma_dpsi, cc; verbose=true)
        end
    end

    @testset "ForwardDiff" begin
        function opti(alpha)
            S = solovev(B0, MS, alpha, qstar; B0_dir=1, Ip_dir=1)
            return (beta_p(S) - 2)^2
        end

        @test (ForwardDiff.derivative(opti, alpha) !== nothing)
    end

    beta_t_of(; R0=R0, Z0=Z0, ϵ=ϵ, κ=κ, δ=δ, B0=B0, alpha=alpha, qstar=qstar) =
        solovev(B0, MillerShape(R0, Z0, ϵ, κ, δ), alpha, qstar; B0_dir=1, Ip_dir=1).beta_t

    @testset "ForwardDiff.derivative ϵ" begin
        @test isfinite(ForwardDiff.derivative(v -> beta_t_of(; ϵ=v), ϵ))
    end

    @testset "ForwardDiff.derivative κ" begin
        @test isfinite(ForwardDiff.derivative(v -> beta_t_of(; κ=v), κ))
    end

    @testset "ForwardDiff.derivative δ" begin
        @test isfinite(ForwardDiff.derivative(v -> beta_t_of(; δ=v), δ))
    end

    @testset "ForwardDiff.derivative R0" begin
        @test isfinite(ForwardDiff.derivative(v -> beta_t_of(; R0=v), R0))
    end

    @testset "ForwardDiff.derivative Z0" begin
        @test ForwardDiff.derivative(v -> beta_t_of(; Z0=v), Z0) ≈ 0.0 atol = 1e-5
    end

    @testset "ForwardDiff.derivative B0" begin
        @test ForwardDiff.derivative(v -> beta_t_of(; B0=v), B0) ≈ 0.0 atol = 1e-6
    end

    @testset "ForwardDiff.derivative alpha" begin
        @test isfinite(ForwardDiff.derivative(v -> beta_t_of(; alpha=v), alpha))
    end

    @testset "ForwardDiff.derivative qstar" begin
        @test isfinite(ForwardDiff.derivative(v -> beta_t_of(; qstar=v), qstar))
    end

end
