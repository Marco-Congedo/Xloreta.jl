using Xloreta, Test

tol = 1e-12

# number of electrodes, data samples, voxels
Ne, Ns, Nv=20, 200, 3000

# fake leadfield in common average reference
K = ℌ(Ne)*randn(Ne, Nv)

# fake data
X=randn(Ne, Ns)

# sample covariance matrix of the fake data
C=(1/Ns)*(X*X')

# fake weights for weighted minimum norm solutions
weights=abs.(randn(Nv))

TsLor1 = sLORETA(K, 1)     # model-driven sLORETA with α=1
TsLor2 = sLORETA(K, 0)    # model-driven sLORETA with α=0
TsLor3 = sLORETA(K, 1, C)  # data-driven sLORETA with α=1

TeLor1 = eLORETA(K, 1)     # model-driven eLORETA with α=1
TeLor1 = eLORETA(K, 1)    # model-driven eLORETA with α=1
TeLor1 = eLORETA(K, 1, C)  # data-driven eLORETA with α=1

@testset "Testing Localization error by Point-Spread functions..." begin
    @test psfLocError(K, TsLor1) == 0
    @test psfLocError(K, TsLor2) == 0
    @test psfLocError(K, TsLor3) == 0

    @test psfLocError(K, TeLor1) == 0
    @test psfLocError(K, TeLor1) == 0
    @test psfLocError(K, TeLor1) == 0
end;