using CurveFit
using Polynomials
using Statistics

struct CalibrCurve{CTA<:Hotwire.AbstractThermalAnemometer}
    "Sensor information"
    sensor::CTA
    "Polynomial fitting the calibration curve"
    fit::Polynomial{Float64}
    "Reference temperature of the calibration"
    T0::Float64
    "Parameters for King's law extrapolation for low speeds"
    king::NTuple{3, Float64}
    "Minimum output voltage with nonzero calibration speed"
    Emin::Float64
end

"""
    cal = CalibrCurve(sensor, V, E, temp, N)

Create a calibration for a thermal anemometer

Parameters:
 * `sensor` an anemometer object (<:AbstractThermalAnemometer)
 * `V` a vector with calibration velocities
 * `E` a vector with calibration voltages
 * `temp` a vector with calibration temperatures or a single calibration temperature
 * `N` degree of polynomial fit

This method will correct all voltages to the mean calibration temperature. Then it will fit the nonzero calibration velocities using a polynomial fit of degree N. 

The speeds between zero and the first nonzero velocity are then interpolated using King's law.

# Examples
```jldoctest
julia> U = [0.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 18.0];

julia> E = [0.979, 1.333, 1.385, 1.423, 1.48, 1.521, 1.551, 1.575, 1.624];

julia> temp = [20.02, 19.75, 19.79, 19.73, 19.67, 19.75, 19.85, 19.85, 20.12];

julia> Rw = Thermistor(5e3, 3500, 25.0);

julia> cta = CTASensor(Rw, 100.0, 20/120);

julia> cal = CalibrCurve(cta, U, E, temp, 5)
CalibrCurve{CTASensor{Thermistor}}(CTASensor{Thermistor}(Thermistor(5000.0, 3500.0, 298.15), 100.0, 0.16666666666666666, 1.0), Polynomial(2526.816258773263 - 3916.1946161935334*x - 1547.7951632370953*x^2 + 6109.915344234665*x^3 - 3944.5331575193045*x^4 + 814.5086726564189*x^5), 19.836666666666666, (0.9897386382393427, 0.5558472575990494, 0.5), 1.332625512625792)

julia> cal.(E, temp) - U  # Check the calibration
9-element Vector{Float64}:
  0.0
  0.00033274230130109217
 -0.002686962091977829
  0.005428070340446567
 -0.007450571094204328
  0.006889739326171096
 -0.0019912770775185606
 -0.0007380705280706223
  0.00021632881207622745

julia> 

```
"""
function CalibrCurve(sensor::CTA, V, E, temp, N=5) where{CTA<:Hotwire.AbstractThermalAnemometer}
    
    nz = V .> 0  # Points with velocity above zero 
    Tm = mean(temp)
   
    # Correct output to reference temperature
    Ec = tempcorr.(Ref(sensor), E, temp, Tm)
        
    V1 = V[nz]
    E1 = Ec[nz]

    # Fit the nonzero velocities to a polynomial of degree N
    fit = curve_fit(Polynomial, E1, V1, N)

    # Now if the speed is lower than te minimum calibration speed, attempt to
    # extrapolate using King's law: E^2 = A + B + U^n
    # For now, n=0.5
    Ea = mean(Ec[.!nz])
    Eb = minimum(E1)
    Vb = fit.(Eb)
    
    A = sqrt(Ea)
    B = (Eb^2 - A) / sqrt(Vb)
    
    return CalibrCurve(sensor, fit, Tm, (A, B, 0.5), Eb)
    
end


"""
    U = cal(E [,temp])

Apply the calibration to obtain wind speed.

Parameters:
 * `E` Output voltage
 * `temp` fluid temperature. If not provided assume calibration temperature
 * Returns velocity

# Examples
```jldoctest
julia> cal(1.5)
6.925627034695708

julia> cal(1.5, 15)
5.902143445069172

julia> cal(1.5, 20)
6.9642311432613395

julia> cal(1.5, 25)
8.309847269677324

```

"""
function (cal::CalibrCurve)(E, temp)
    Ec = tempcorr(cal.sensor, E, temp, cal.T0)
    
    if Ec < cal.Emin
        A = cal.king[1]
        B = cal.king[2]
        n = cal.king[3]
        if Ec < sqrt(A)
            U = zero(E)
        else
            U = ( (Ec^2 - A) / B ) ^ (1/n)
        end
    else
        U = cal.fit.(Ec)
    end
    
    return U
end

(cal::CalibrCurve)(E) = cal(E, cal.T0)
