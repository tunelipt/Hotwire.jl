# Operating conditions correction

The output of a thermal anemometer depends onde the electronic characteristics of the sensor but on the fluid conditions as well. Temperature, pressure and composition will affect the measurements if they are different from the calibration conditions. Changes in the properties of the fluid will mean changes in anemometer output and the result should be corrected accordingly.

## The basic idea of output correction

The convection from a heated element is modeled using [Newton's law of cooling](https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling):

``\dot{Q} = h\cdot A\cdot \left(T_s - T_a\right)``

In this equation,

 * ``\dot{Q}`` is the heat being dissipated
 * ``h`` is the convection coefficient that depends on geometry, fluid properties and *fluid velocity*
 * ``A`` is the surface area that is loosing heat
 * ``T_s`` is the surface temperature
 * ``T_a`` is the ambient temperature

What couples the electric part of anemometer with the heat transfer is the Joule heating, such that

``\dot{Q} = E^2 \cdot I_w = \frac{E^2}{R_w}``

where ``E`` is the voltage across the heating element, ``I_w`` is the electric current and ``R_w`` is the electrical resistance that changes with temperature.

A calibration consists in measuring the output voltage ``E_c`` for different fluid velocities ``U_c`` given an electrical circuit, a calibration fluid (usually air) at a given pressure ``P_c`` and temperature ``T_c``.

If, when using the anemometer, any of these parameters are changed, the output ``E`` of the anemometer also changes and the idea here is to find an equivalent output in calibration conditions.

In what follows, we will assume that the circuit is a constant temperature anemometer.


## Correction parameters

The most common type of correction is due to change in fluid temperature. The pressure might also change. Another, less common situation in the case of anemometers but important when dealing the flow meters, is changing the working fluid. The calibration could be carried out using air but the measurements are done in helium for instance.

The electrical configuration of the circuit might also change. In the case of a CTA anemometer this might correspond to a change in operating resistance ``R_w``. In the case of constant current anemometers this corresponds to a change in current ``I_w``.

The geometry of the sensor will always be considered unchanging and dynamic effects are not considered here.

## Simplest correction: constant fluid properties

If the thermodynamic and transport properties of the fluid do not change, the output of the anemometer in calibration conditions

``E_c^2 = h(U_c) \cdot A \cdot R_{wc} \cdot (T_{wc} - T_{ac})``

here we assume that the surface temperature is the mean temperature of the resistor element. Since the properties do not change during measurement, the sensor resistance ``R_w`` and the corresponding sensor temperature might be different. The room temperature might also be different. Therefore the output _for the same fluid velocity_ is:

``E^2 = h(U_c) \cdot A \cdot R_w \cdot (T_wc - T_a)``

under these conditions, dividing both equations, we get:

``E_c = E \cdot \sqrt{ \frac{R_{wc}}{R_w}\cdot\frac{T_{wc}-T_{ac}}{T_w - T_c} }``


## Changes in fluid properties

If the fluid properties change, either due to changes in room temperature and pressure or changes in fluid composition, the correction is a little more complex.

In heat transfer, it is common to use the Nusselt number to calculate the convection coefficient. The Nusselt number (``Nu``) depends on Reynolds number (``Re``) and Prandtl number (``Pr``):

``Nu = \frac{h\cdot L}{k} = Nu(Re,Pr)``

where

 * ``L`` is a characteristic dimension of the sensor
 * ``k`` is the heat conductivity of the fluid

Searching the literature for correlations ``Nu(Re,Pr)`` for different geometries, it is very common to find

``Nu(Re,Pr) = f(Re) \cdot Pr^n``

where ``n`` is an exponent, usually close to 1/3, 0.3, 0.4 or 0.25. Fixing this value ``n``, we get:

``h = f(Re) Pr^n \times \frac{k}{L}``

so that

``E^2 = f(Re) Pr^n \times \frac{k}{L} \cdot R_w \cdot (T_w - T_a)``

During calibration

``E_c^2 = f(Re_c) Pr_c^n \times \frac{k_c}{L} \cdot R_{wc} \cdot (T_{wc} - T_{ac})``

Given ``E``to calculate the corresponding calibration value ``E_c``, _we assume_ that ``Re_c = Re`` and divide one equation by the other:

``E_c^2 = E^2 \times \frac{Pr_c^n}{Pr^n} \cdot \frac{k_c}{k} \cdot \frac{R_{wc}}{R_w} \cdot \frac{T_{wc} - T_{ac}}{T_w - T_a}``

from the calibration curve, we calculate ``U_c``. To obtain the fluid velocity, remember that we assumed that ``Re = Re_c``:


``U = U_c \times \frac{\rho_c}{\rho} \cdot \frac{\mu}{\mu_c}``



## Different sensor configurations

Up to this point, it was assumed that the resistor temperature ``T_w`` is the same as the surface temperature ``T_s`` that is exposed to the fluid. This is not always the case! There might be an protection, the electrical leads might be large and also transfer heat. Each case should be analyzed individually. This package assumes that the format for the heat transfer is:

   ``E^2 = Y(Re, T_a, P_a, x, \ldots) \cdot R_w\cdot (T_w - T_a)``

In this equation ``\ldots`` is any set of parameters used to characterize the model. As an example we will consider a sensor that is a thermistos with wire leads and a glass capsule around it.

## Glass bead thermistors

As an example of correction, we will model the heat transfer from MF-58 glass bead thermistors. These thermistors are not really appropriate for velocity measurements but they are extremely cheap and easily obtainable. What makes these thermistors bad for anemometry is that they have a large glass insulating shell around the thermistor element and the electrical leads are large.

The model developed here, assumes that the outer surface has a mean temperature ``T_s`` that is related to the thermistor temperature by the following relation:

``T_s = T_w - \beta \cdot \dot{Q}``

where
 * ``T_s`` is the mean outside surface temperature of the glass bead
 * ``T_w`` is the temperature of the internal, actual thermistor element
 * ``\beta`` is a coefficient that takes into account the geometry and thermal conductivity of the glass bead
 * ``\dot{Q}`` is the Joule heating of the internal thermistor element

The heat loss through the leads should also be taken into consideration. Each lead, will be modeled as a fin whose base temperature equals to the mean surface temperature ``T_s``. The fin is assumed to be very long (_infinitely long_). The temperature distribution of each lead has an expoential character:

``T(x) - T_a = (T_s - T_a)\cdot\exp(-m\cdot x)``

where

``m^2 = \frac{h_f\cdot P_f}{k_f\cdot A_f}``

with

 * ``h_f`` is the convection coefficient of the leads
 * ``P_f`` is the perimeter of the lead
 * ``k_f`` is the thermal conductivity of the lead material
 * ``A_f`` is the area of the cross section of the lead

With these hypotheses, the flow loss of the sensor is

``\dot{Q} = \dot{Q_1} + \dot{Q_2} = hA(T_s - T_a) + N\cdot k_fA_f\sqrt{\frac{h_fP_f}{k_fA_f}} (T_s - T_a)``

In this equation, ``N`` is the number of leads leaving the sensor. It will usually be 2 or 1 (when one of the leads is not exposed).


Combining the equations above, we have

``E^2 = \frac{X}{1 + \beta\cdot X} (T_w - T_s)``

where

``X = hA + N\cdot k_fA_f\sqrt{\frac{h_fP_f}{k_fA_f}}``

We still have one problem, we do not know either `h` or `h_f`. We will postulate that

``h_f = \gamma h``

To estimate the value of ``gamma``, the idea is to use most common forms of expression for ``Nu = Nu(Re,Pr)``:

``Nu = a \cdot Re^q \cdot Pr^n``

For the body of the sensor, we have a characteristic length ``D`` and for the leads, this is ``d``. Using the above approximation we get

``\gamma \approx \left( \frac{d}{D} \right)^q``

The parameter ``q`` is usually is of the order of 0.3-0.5. A value will have to be selected.

Again we use the approximation for ``h``

``h = f(Re) \times \frac{k}{D} \times Pr^n``

Now, we have:

``X = hA + N\sqrt{\gamma h k_f A_f P} = f(Re)\frac{k\cdot Pr^n}{D}\cdot A + N\sqrt{\frac{f(Re)\cdot k \cdot Pr^n\gamma k_f A_f P_f}{D}} = c_1 \cdot f(Re) \phi + c_2 \sqrt{ f(Re)\phi}``

where
 * ``\phi = k\cdot Pr^n``
 * ``c_1 = A/D``
 * ``c_2 = N\sqrt{\gamma k_f A_f P_f / D}``




And the expression for ``Y`` is:

``Y = \frac{X}{1+\beta X}``

### Correction procedure

The objective here is to determine, _for the same Reynolds number_ ``Re\equiv Re_c``, the corresponding calibration output voltage ``E_c``. Since ``R_w``, ``T_w`` is known from electronic configuration, knowing the room temperature, the room temperature ``T_a`` is measured, ``Y`` is calculated as

``Y = \frac{E^2}{R_w(T_w - T_a)}``

Now we can calculate ``X`` and from that ``f(Re)`` is obtained. Using calibration conditions, ``X_c`` and ``Y_c`` are calculated and we finally obtain the corrected output voltage corresponding to calibration conditions ``E_c``. Using the calibration curve, ``U_c`` can be calculated and as above, 


``U = U_c \times \frac{\rho_c}{\rho} \cdot \frac{\mu}{\mu_c}``


Now we have an algorithm that can be used to correct anemometer output for different operating fluid, temperature and pressure.



### Correction for constant current anemometer (CCA)

When using the constant temperature anemometer, ``R_w`` and therefore ``T_w`` where fixed. This simplifies _a lot_ the calculations. In the case of constant temperature anemometer things are trickier. Since the current is maintained constant, as velocity increases, temperature decreases. But the output voltage is directly related to the constant temperature ``I_w`` and therefore we can determine the sensor resistance and temperature:

``R_w = \frac{E}{I_w}``

The idea here is to select a reference temperature and from the calibration curve and use the same procedure used for the CTA sensor to correct the value to this reference temperature (and resistance). The oposite procedure 


