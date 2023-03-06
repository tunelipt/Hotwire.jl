# Testing probe3d stuff

let

    k² = (0.0, 0.0, 0.0)
    h² = (1.0, 1.0, 1.0)

    U = collect(1.0:0.5:10.0)
    T0 = 20.0+273.15
    E1 = sqrt.(2.0 .+ 0.7 .* U .^ 0.42)
    E2 = sqrt.(2.2 .+ 0.65 .* U .^ 0.42)
    E3 = sqrt.(1.95 .+ 0.75 .* U .^ 0.41)
    cal = (CalibrCurve(E1, U, KingFit(E1, U), T0),
           CalibrCurve(E2, U, KingFit(E2, U), T0),
           CalibrCurve(E3, U, KingFit(E3, U), T0))

    wires = (Wire(5.0, 5.0*1.8, T0),
             Wire(5.2, 5.2*1.8, T0),
             Wire(5.3, 5.3*1.8, T0))
    
    probe = Probe3d(wires, cal, k², h²)
    for i in 1:length(U)
        Uⁱ = probe(E1[i], E2[i], E3[i], T0)
        @test Uⁱ[1] ≈ U[i] rtol=1e-4
    end


    k² = (0.03, 0.04, 0.05)
    h² = (1.05, 1.08, 1.02)

    probe = Probe3d(wires, cal, k², h²)
    for i in 1:length(U)
        Uⁱ = probe(E1[i], E2[i], E3[i], T0)
        @test Uⁱ[1] ≈ U[i] rtol=1e-4
    end

end
