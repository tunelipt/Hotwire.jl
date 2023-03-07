# Testing probe3d stuff

let

    k² = (0.0, 0.0, 0.0)
    h² = (1.0, 1.0, 1.0)

    T0 = 20.0+273.15
    E = collect(1.0:20.0)
    trivialcal = CalibrCurve(E, E, x->x, 293.15)
    wire1 = Wire(3.5, 3.5*1.8, T0)
    wire2 = Wire(3.7, 3.7*1.8, T0)
    wire3 = Wire(3.6, 3.6*1.8, T0)
    wires = (wire1, wire2, wire3)

    cal = (trivialcal, trivialcal, trivialcal)
    
    probe = Probe3d(wires, cal, k², h²)

    Ux, Uy, Uz = probe(1.0, 1.0, 1.0, T0)
    @test Ux ≈ 1.0
    @test Uy ≈ 0.0 atol=1e-5
    @test Uz ≈ 0.0 atol=1e-5
    

    k² = (0.03, 0.04, 0.05)
    h² = (1.05, 1.08, 1.02)
    probe = Probe3d(wires, cal, k², h²)

    Ux, Uy, Uz = probe(1.0, 1.0, 1.0, T0)
    @test Ux ≈ 1.0
    @test Uy ≈ 0.0 atol=1e-5
    @test Uz ≈ 0.0 atol=1e-5
    
end
