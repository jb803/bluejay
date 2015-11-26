Rough Overall Design
========

Aim is to provide an alternative to the Pro98 M2505 

M2505 Features
------

The M2505 has the following properties: 

* Impulse:	7450Ns
* F_avg: 	2491N
* m_prop: 	3.34kg
* m_tot:	6.25kg

BlueJay Features
------

From the varPo.m script we find that at a low chamber pressure of 10bar
we can produce a feasible rocket (i.e. TWR>5) with:

* Impulse: 	7450Ns
* F_avg:	1kN
* m_prop:	5.21kg

This assumes that the combustion gases have properties with values typical of those seen
in the exit of gas turbine combustors.

The mass fraction of this vehcicle is ~25%

This requires:
* m_dot:	0.6989kgs^-1

Assuming a stoichiometric mixture:

C3H7OH + 9N2O -> N2 + 3CO2 + 4H20

This implies Y_f = 0.1316 and hence m_dot_f = 0.092 kgs^-1

Fuel Pump
------

As Ro_C3H7OH = 786kgm^-3 then we get a volumetric flow rate of Q=1.17e-4 m^3s^-1

In order to prevent blowback due to combustion instabilities we assume: p_f = 1.5*p0_chamber

This requires a head of 195m

For a 40mm rotor ( roughly the smallest size we can easily make ) a specific diameter of 24.5. According to the Cordier diagram this is at the very
bottom of a radial pump. Increasig the volumetric flow rate may help. For the time being (25/11/15) a turbopump will be designed to this specification and we wil have to up F_avg if nec.


Oxidizer Pump
-------

If nitrous were treated as a cryogen then it may be preferable to use another pump for the oxidiser.

The density of liquid nitrous is 1230 kgm^-3

Again, assuming stoichiometry, this gives a volumetric flow rate of 4.93e-4 m^3s^-1

This results in a specific diameter of 10.6 - much more comfortably within the region dominated by radial pumps.

Plan
------
The plan is to:

1. Embark on matlab scripts to non-dimensionally design the fuel pump
2. Implement combustion modelling to calculate the actual temperature of the mix
3. Refine the rough design outline
