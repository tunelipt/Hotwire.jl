

"""
    tempcorr(w::CTASensor, temp)

Corrects the anemometer output due to differences in flow temperature.

 * `w` `CTASensor` object
 * `temp` Temperature

This correction assumes that only the temperature changed and not fluid properties. Another function that takes that into account will be developped later.

The correction works by using the definition of convective heat transfer:

``Q = h A (Tw - Too)``

where ``Tw`` is the operating temperature and ``Too`` is the flow temperature. Imagine that for the *same* flow speed the flow temperature changes to ``Too'``. Now, the assuming the the convection coefficient (``h``) remains the same, the heat output will change:

``Q' = h A (Tw - Too')``

Since the anemometer output is proportial to the voltage accross the sensor element, dividing the above equations and setting

``Q = α²E²/Rw``

we get the following equation for output correction:

``E' = E × ( (Tw - Too') / (Tw - Too) )^(1/2)``

If between the resistance element and the flow is an insulator (very common in thermistors and hot film probes), the correction equation remains the same event though this insulator should be taken into account:

``Q = h A (Tw - Too - Q⋅β)``

where ``Q⋅β`` corrects the sensor temperature (``Tw``) to the outer surface mean temperature.

"""
function tempcorr(w::CTASensor, E, temp)
    Tw = optemperature(w)
    Too = temperature(w.R)
    return E * sqrt( (Tw - Too) / (Tw - temp) )
end



"""
    tempcorr(w::CCASensor, temp)

Corrects the anemometer output due to differences in flow temperature.

 * `w` `CCASensor` object
 * `temp` Temperature

This correction assumes that only the temperature changed and not fluid properties. Another function that takes that into account will be developped later.

The correction works by using the definition of convective heat transfer:

``Q = h A (Tw - Too)``

where ``Tw`` is the operating temperature and ``Too`` is the flow temperature. Imagine that for the *same* flow speed the flow temperature changes to ``Too'``. Now, the assuming the the convection coefficient (``h``) remains the same, the heat output will change:

``Q' = h A (Tw - Too')``

Since the anemometer output is proportial to the voltage accross the sensor element, dividing the above equations and setting

``Q = α²E²/Rw``

we get the following equation for output correction:

``E' = E × ( (Tw - Too') / (Tw - Too) )^(1/2)``

If between the resistance element and the flow is an insulator (very common in thermistors and hot film probes), the correction equation remains the same event though this insulator should be taken into account:

``Q = h A (Tw - Too - Q⋅β)``

where ``Q⋅β`` corrects the sensor temperature (``Tw``) to the outer surface mean temperature.

"""
function tempcorr(w::CCASensor{Resistor}, E, temp)

    α = w.R.α
    To = temperature(w.R)
    return E / (1.0 + α*(temp - To))
end

function funroot(f, x1, x2, eps=1e-7, maxiter=100)
    y1 = f(x1)
    y2 = f(x2)
    if abs(y1) < eps return x1 end
    if abs(y2) < eps return x2 end
        
    err = 0.0
    
    for i in 1:maxiter
        a = (y2-y1) / (x2 - x1)
        b = y1 - a*x1
        xn = -b/a
        err = abs(xn - x2)
        if err < eps
            return xn
        end
        x1 = x2
        x2 = xn
        y1 = y2
        y2 = f(x2)
        if abs(y2) < eps
            return x2
        end
        
    end

    error("funroot failed to converge after $maxiter iterations and with $err error!")
    return x2
end


function tempcorr(w::CCASensor, E, temp, eps=1e-7, maxiter=100)

    R = w.R
    g = gain(w)
    I = current(w)
    Tw = temperature(R, E/(g*I))
    To = temperature(w.R)
    ΔT = Tw - temp

    f = Eo -> E*(temperature(R, Eo/(g*I)) - To) - Eo * ΔT
    return funroot(f, E, 1.01*E, eps, maxiter)
    #error("tempcorr(w::CCASensor, $E, $temp failed to converge in $maxiter iterations with $err error!")
    return Eo
end



