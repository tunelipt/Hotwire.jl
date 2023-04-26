using Statistics


abstract type AbstractCalibr end
abstract type AbstractCalibr1d end

struct CalibrCurve{T,Fit} <: AbstractCalibr1d
    "Anemometer calibration output"
    E::Vector{T}
    "Calibration Velocity"
    U::Vector{T}
    "Polynomial fitting the calibration curve"
    fit::Fit
    "Reference temperature of the calibration"
    T0::T
end
Base.broadcastable(cal::CalibrCurve) = Ref(cal)

reftemp(cal::CalibrCurve) = cal.T0

"""
    cal = calibr_curve(sensor, V, E, temp, N)

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

julia> temp = [20.02, 19.75, 19.79, 19.73, 19.67, 19.75, 19.85, 19.85, 20.12] .+ 273.15;

julia> Rw = Thermistor(5e3, 3500, 25.0+273.15);

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
function calibr_curve(sensor::CTA, V, E, temp, fitfun::Function) where{CTA<:Hotwire.AbstractThermalAnemometer}
    
    nz = V .> 0  # Points with velocity above zero 
    Tm = mean(temp)
   
    # Correct output to reference temperature
    Ec = tempcorr.(Ref(sensor), E, temp, Tm)
    Uc = [u for u in V]
    
    # Fit the nonzero velocities to a polynomial of degree N
    fit = fitfun(Ec, Uc) 

    return CalibrCurve(Ec, Uc, fit, Tm) 
end

calibr_curve(sensor, V, E, temp, N=5) =
    calibr_curve(sensor, V, E, temp, make_poly_fitfun(N))



"""
    `U = cal(sensor, E [,temp])`

Apply the calibration to obtain wind speed.

Parameters:
 * `sensor` Structure representing the anemometer
 * `E` Output voltage
 * `temp` fluid temperature. If not provided assume calibration temperature
 * Returns velocity

# Examples
```jldoctest
julia> cal(w, 1.5)
6.925627034695708

julia> cal(w, 1.5, 15)
5.902143445069172

julia> cal(w, 1.5, 20)
6.9642311432613395

julia> cal(w, 1.5, 25)
8.309847269677324

```

"""
function (cal::CalibrCurve)(sensor, E, temp)
    Ec = tempcorr(sensor, E, temp, cal.T0)
    
    return cal.fit(Ec)
end

(cal::CalibrCurve)(sensor, E) = cal(sensor, E, cal.T0)
