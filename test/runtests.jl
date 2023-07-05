using DataDrivenAcoustics
using Test


function test2d(datapm)
    x1 = transfercoef(datapm, nothing, AcousticReceiver(50.0, -5.0))
    x2 = transfercoef(datapm, nothing, AcousticReceiver(50.0, -10.0))
    x3 = transfercoef(datapm, nothing, AcousticReceiver(50.0, -15.0))
    x = transfercoef(datapm, nothing, [AcousticReceiver(50.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x


    x = transfercoef(datapm, nothing, AcousticReceiverGrid2D(50.0, 0.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] == x
    x = transfercoef(datapm, nothing, AcousticReceiverGrid2D(50.0, 10.0, 3, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (3, 3)
    @test [x1, x2, x3] == x[1,:]

    x1 = transmissionloss(datapm, nothing, AcousticReceiver(50.0, -5.0))
    x2 = transmissionloss(datapm, nothing, AcousticReceiver(50.0, -10.0))
    x3 = transmissionloss(datapm, nothing, AcousticReceiver(50.0, -15.0))
    x = transmissionloss(datapm, nothing, [AcousticReceiver(50.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x
    x = transmissionloss(datapm, nothing, AcousticReceiverGrid2D(50.0, 0.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] == x
    x = transmissionloss(datapm, nothing, AcousticReceiverGrid2D(50.0, 10.0, 3, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (3, 3)
    @test [x1, x2, x3] == x[1,:]
end

function test3d(datapm)
    x1 = transfercoef(datapm, nothing, AcousticReceiver(50.0, 0.0, -5.0))
    x2 = transfercoef(datapm, nothing, AcousticReceiver(50.0,0.0, -10.0))
    x3 = transfercoef(datapm, nothing, AcousticReceiver(50.0, 0.0, -15.0))
    x = transfercoef(datapm, nothing, [AcousticReceiver(50.0, 0.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x


    x = transfercoef(datapm, nothing, AcousticReceiverGrid3D(50.0, 0.0, 1, 0.0, 1.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] == x
    x = transfercoef(datapm, nothing, AcousticReceiverGrid3D(50.0, 10.0, 3, 0.0, 1.0, 2, -5.0, -5.0, 3))
    @test x isa AbstractArray
    @test size(x) == (3, 2, 3)
    @test [x1, x2, x3] == x[1, 1,:]


    x1 = transmissionloss(datapm, nothing, AcousticReceiver(50.0, 0.0,  -5.0))
    x2 = transmissionloss(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0))
    x3 = transmissionloss(datapm, nothing, AcousticReceiver(50.0, 0.0, -15.0))
    x = transmissionloss(datapm, nothing, [AcousticReceiver(50.0, 0.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x
    x = transmissionloss(datapm, nothing, AcousticReceiverGrid3D(50.0, 0.0, 1, 0.0, 1.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] == x
    x = transmissionloss(datapm, nothing, AcousticReceiverGrid3D(50.0, 10.0, 3, 0.0, 1.0, 2, -5.0, -5.0, 3))
    @test x isa AbstractArray
    @test size(x) == (3, 2, 3)
    @test [x1, x2, x3] == x[1,1,:]
end


@test RayBasis2D in models()
@test RayBasis2DCurv in models()
@test RayBasis3D in models()
@test RayBasis3DRCNN in models()
@test GPR in models()


env = UnderwaterEnvironment()
pm = PekerisRayModel(env, 7)

Random.seed!(1)

txpos = [0.0, -5.0]
rxpos = rand(2, 500) .* [80.0, -20.0] .+ [1.0, 0.0]
tloss = Array{Float32}(undef, 1, size(rxpos)[2])
for i in 1 : 1 : size(rxpos)[2]
    tloss[1, i] = Float32(transmissionloss(pm, AcousticSource(txpos[1], txpos[2], 1000.0), AcousticReceiver(rxpos[1,i], rxpos[2,i]); mode=:coherent))
end
dataenv = DataDrivenUnderwaterEnvironment(rxpos, tloss; frequency = 1000.0, soundspeed = 1540.0);



datapm = RayBasis2D(dataenv; inilearnrate = 0.005, seed = true)
@test datapm isa RayBasis2D
x = transfercoef(datapm, nothing, AcousticReceiver(50.0, -10.0)) 
@test x isa Complex
@test imag(x) ≈ -0.0089012 atol=0.0002
@test abs(x) ≈ 0.00921678 atol=0.0002
y = transmissionloss(datapm, nothing, AcousticReceiver(50.0, -10.0))
@test -10 * log10(abs2(x)) ≈ y atol=0.1
test2d(datapm)


arr = arrivals(datapm, nothing, AcousticReceiver(50.0, -10.0))
@test arr isa AbstractVector{<:RayArrival}
@test length(arr) == 44
@test all([(amp2db(abs(arr[j].phasor)) < amp2db(abs(arr[j-1].phasor))) for j ∈ 2:44])
@test abs(arr[1].phasor) ≈ 0.035303 atol=0.001
@test real(arr[2].phasor) < 0.0
@test imag(arr[2].phasor) ≈ 0.022897 atol = 0.001


datapm = RayBasis2DCurv(dataenv; inilearnrate = 0.005, seed = true)
@test datapm isa RayBasis2DCurv
x = transfercoef(datapm, nothing, AcousticReceiver(50.0, -10.0)) 
@test x isa Complex
@test imag(x) ≈ -0.012282 atol=0.0002
@test abs(x) ≈ 0.029002 atol=0.0002
y = transmissionloss(datapm, nothing, AcousticReceiver(50.0, -10.0))
@test -10 * log10(abs2(x)) ≈ y atol=0.1
test2d(datapm)

arr = arrivals(datapm, nothing, AcousticReceiver(50.0, -10.0))
@test arr isa AbstractVector{<:RayArrival}
@test length(arr) == 56
@test all([(amp2db(abs(arr[j].phasor)) < amp2db(abs(arr[j-1].phasor))) for j ∈ 2:56])
@test abs(arr[1].phasor) ≈ 0.019757 atol = 0.001
@test real(arr[2].phasor) ≈ -0.0177007 atol = 0.001
@test imag(arr[2].phasor) ≈ 0.0073227 atol = 0.001


kern = Mat(1/2, 0.0, 0.0)
datapm = GPR(dataenv, kern; logObsNoise = -5.0, seed = true, ratioₜ = 1.0)
@test datapm isa GPR
test2d(datapm)
x = transfercoef(datapm, nothing, AcousticReceiver(50.0, -10.0)) 
@test abs(x) ≈ 0.024228 atol=0.0002
y = transmissionloss(datapm, nothing, AcousticReceiver(50.0, -10.0))
@test -10 * log10(abs2(x)) ≈ y atol=0.1



Random.seed!(1)
txpos = [0.0, 0.0, -5.0]
rxpos = rand(3, 500) .* [100.0, 0.0, -20.0] .+ [1.0, 0.0, 0.0];
tloss = Array{Float32}(undef, 1, size(rxpos)[2])
for i in 1 : 1 : size(rxpos)[2]
    tloss[1, i] = Float32(transmissionloss(pm, AcousticSource(txpos[1], txpos[2], txpos[3], 1000.0), AcousticReceiver(rxpos[1,i], rxpos[2,i], rxpos[3,i]); mode=:coherent))
end
dataenv = DataDrivenUnderwaterEnvironment(rxpos, tloss; frequency = 1000.0, soundspeed = 1540.0);
datapm = RayBasis3D(dataenv; inilearnrate = 0.005, seed  = true)
@test datapm isa RayBasis3D
test3d(datapm)

x = transfercoef(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0)) 
@test x isa Complex
@test imag(x) ≈ -0.028361 atol=0.0002
@test abs(x) ≈ 0.028391 atol=0.0002
y = transmissionloss(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0))
@test -10 * log10(abs2(x)) ≈ y atol=0.1

arr = arrivals(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0))
@test arr isa AbstractVector{<:RayArrival}
@test length(arr) == 52
@test all([(amp2db(abs(arr[j].phasor)) < amp2db(abs(arr[j-1].phasor))) for j ∈ 2:52])
@test abs(arr[1].phasor) ≈ 0.0197353 atol=0.001
@test real(arr[2].phasor) ≈ -0.011336 atol = 0.001
@test imag(arr[2].phasor) ≈ -0.016114 atol = 0.001



Random.seed!(1)

RCNN = Chain(  
    x -> (x ./ 0.5f0 .* π .- 0.5f0) .* 2.0f0, #normalization of incident angle
    Dense(1, 30, sigmoid),
    Dense(30, 50, sigmoid),  
    Dense(50, 2),
)
dataenv = DataDrivenUnderwaterEnvironment(rxpos, tloss; frequency = 1000.0, soundspeed = 1540.0, waterdepth = 20.0, tx = AcousticSource(0.0, 0.0, -5.0, 1000.0))
datapm = RayBasis3DRCNN(dataenv, RCNN; seed = true, inilearnrate = 0.05, ncount = 500)
@test datapm isa RayBasis3DRCNN
x = transfercoef(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0)) 
@test x isa Complex
@test imag(x) ≈ -0.0017830 atol=0.0002
@test abs(x) ≈ 0.029992 atol=0.0002
y = transmissionloss(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0))
@test -10 * log10(abs2(x)) ≈ y atol=0.


arr = arrivals(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0))

@test arr isa AbstractVector{<:RayArrival}
@test length(arr) == 6
@test all([(amp2db(abs(arr[j].phasor)) < amp2db(abs(arr[j-1].phasor))) for j ∈ 2:6])
@test abs(arr[1].phasor) ≈ 0.019894 atol=0.001
@test real(arr[2].phasor)≈ 0.015285 atol = 0.001
@test imag(arr[2].phasor) ≈ -0.0115367 atol = 0.001



kern = Mat(1/2, [0.0, 0.0, 0.0], 0.0)
datapm = GPR(dataenv, kern; logObsNoise = -5.0, seed = true, ratioₜ = 1.0)
@test datapm isa GPR
test3d(datapm)
x = transfercoef(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0)) 
@test abs(x) ≈ 0.023884 atol=0.0002
y = transmissionloss(datapm, nothing, AcousticReceiver(50.0, 0.0, -10.0))
@test -10 * log10(abs2(x)) ≈ y atol=0.1






