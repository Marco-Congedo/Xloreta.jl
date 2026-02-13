using Xloreta, Test, LinearAlgebra

tol = 1e-12

# number of electrodes, data samples, voxels
Ne, Ns, Nv=20, 200, 3000

# fake leadfield in common average reference
K = ℌ(Ne)*randn(Ne, Nv)

# fake data
X=randn(Ne, Ns)

# sample covariance matrix of the fake data
C=Symmetric((1/Ns)*(X*X'))

# fake weights for weighted minimum norm solutions
weights=abs.(randn(Nv))

# run only as the localization error is not zero
Tmn1 = minNorm(K, 1)    # unweighted model-driven min norm with α=1
Tmn2 = minNorm(K, 10)   # unweighted model-driven min norm with α=10
Tmn3 = minNorm(K, 1; W=weights) # weighted model-driven min norm with α=1
Tmn4 = minNorm(K, 1, C) # data-driven min norm with α=1


TsLor1 = sLORETA(K, 1)     # model-driven sLORETA with α=1
TsLor2 = sLORETA(K, 0)    # model-driven sLORETA with α=0
TsLor3 = sLORETA(K, 1, C)  # data-driven sLORETA with α=1

@testset "Testing Localization error of sLORETA..." begin
    @test psfLocError(K, TsLor1) == 0
    @test psfLocError(K, TsLor2) == 0
    @test psfLocError(K, TsLor3) == 0
end;

TeLor1 = eLORETA(K, 1; verbose=false)     # model-driven eLORETA with α=1
TeLor1 = eLORETA(K, 1; verbose=false)    # model-driven eLORETA with α=1
TeLor1 = eLORETA(K, 1, C; verbose=false)  # data-driven eLORETA with α=1

@testset "Testing Localization error of eLORETA..." begin
    @test psfLocError(K, TeLor1) == 0
    @test psfLocError(K, TeLor1) == 0
    @test psfLocError(K, TeLor1) == 0
end;