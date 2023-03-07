# Testing probe3d stuff

let

    k² = (0.0, 0.0, 0.0)
    h² = (1.0, 1.0, 1.0)

    U = collect(1.0:0.5:10.0)
    T0 = 20.0+273.15
    trivialcal(E,T) = E
    wire1 = Wire(3.5, 3.5*1.8)
    wire2 = Wire(3.7, 3.7*1.8)
    wire3 = Wire(3.6, 3.6*1.8)
    wires = (wire1, wire2, wire3)

    cal = (trivialcal, trivialcal, trivialcal)
    
    probe = Probe3d(wires, cal, k², h²)
    for i in 1:length(U)
        Uⁱ = probe(3, 3, 3, T0)
        println(U¹)
        #@test Uⁱ[1] ≈ U[i] rtol=1e-4
    end


    k² = (0.03, 0.04, 0.05)
    h² = (1.05, 1.08, 1.02)
    
    #probe = Probe3d(wires, cal, k², h²)
    #for i in 1:length(U)
    #    Uⁱ = probe(E1[i], E2[i], E3[i], T0)
    #    @test Uⁱ[1] ≈ U[i] rtol=1e-4
    #end

end
