
dir = @__DIR__
g1 = readg(dir * "/g150219.03200") # (+Bt, +Ip)
g2 = readg(dir * "/g133221.01150") # (-Bt, +Ip)
g3 = readg(dir * "/g159177.02700") # (+Bt, -Ip)
g4 = readg(dir * "/g153298.04400") # (-Bt, -Ip)

btip1 = (1, 1)
btip2 = (-1, 1)
btip3 = (1, -1)
btip4 = (-1, -1)

# COCOS for the GEQDSK files
cc1 = 1
cc2 = 5
cc3 = 3
cc4 = 7

M1 = efit(g1; clockwise_phi=false)
M2 = efit(g2; clockwise_phi=false)
M3 = efit(g3; clockwise_phi=false)
M4 = efit(g4; clockwise_phi=false)

r1 = (M1.axis[1] + 0.1, M1.axis[2])
r2 = (M2.axis[1] + 0.1, M2.axis[2])
r3 = (M3.axis[1] + 0.1, M3.axis[2])
r4 = (M4.axis[1] + 0.1, M4.axis[2])

test_data = ((cc1, g1, M1, r1, btip1), (cc2, g2, M2, r2, btip2), (cc3, g3, M3, r3, btip3), (cc4, g4, M4, r4, btip4))

@testset "EFIT Tests" begin

    @testset "Bt-Ip Signs" begin
        for (cc, g, M, r, btip) in test_data
            @test (sign(g.bcentr), sign(g.current)) == btip
            Bt, Ip = (sign(Bfield(M, r...)[2]), sign(Jfield(M, r...)[2]))
            @test (Bt, Ip) == btip
        end
    end

    @testset verbose = true "CurlB == μ₀J" begin
        for (cc, g, M, r, btip) in test_data
            @test curlB(M, r...) ≈ mu0 * Jfield(M, r...) rtol = 0.01
        end
    end

    @testset verbose = true "3D tests" begin
        for (cc, g, M, r, btip) in test_data
            @test Bfield(M, r...) ≈ Bfield(M, r[1], zero(r[1]), r[2]) rtol = 0.01
            @test Jfield(M, r...) ≈ Jfield(M, r[1], zero(r[1]), r[2]) rtol = 0.01

            @test gradB(M, r...) ≈ gradB(M, r[1], zero(r[1]), r[2]) rtol = 0.01
            @test gradB(M, r[1], zero(r[1]), r[2]) ≈ MXHEquilibrium.gradB_autodiff(M, r[1], zero(r[1]), r[2]) rtol = 0.01

            @test curlB(M, r...) ≈ curlB(M, r[1], zero(r[1]), r[2]) rtol = 0.01
            @test curlB(M, r[1], zero(r[1]), r[2]) ≈ MXHEquilibrium.curlB_autodiff(M, r[1], zero(r[1]), r[2]) rtol = 0.01
        end
    end

    @testset verbose = true "COCOS Identify" begin
        for (cc, g, M, r, btip) in test_data
            @test identify_cocos(g; clockwise_phi=false) == (cc,)
        end
    end

    @testset verbose = true "COCOS Consistancy" begin
        for (cc, g, M, r, btip) in test_data
            @test check_cocos(g, cc)
        end
    end

    @testset verbose = true "COCOS Transformation" begin
        for (cc_in, g, M, r, btip) in test_data
            for cc in 1:8
                _cc = cocos(cc)
                g_new = transform_cocos(g, cc_in, _cc)
                @test identify_cocos(g_new; clockwise_phi=_cc.sigma_RpZ < 0) == (cc,)
                @test check_cocos(g_new, cc)
                @test curlB(M, r...) ≈ mu0 * Jfield(M, r...) rtol = 0.01
            end
        end
    end

end
