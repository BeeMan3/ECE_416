Diode I-V Sweep

** This is a template to get you started in the right direction
** It won't necessarily have everything you need to extract your 
** I-V curves but should get you a long way.

.options post

* FILL IN THE PARAMETERS below to match the values you calculated 
*IS - Reverse saturation current (A)
*CJ0 - Zero-bias Junction Capacitance (F/m2)
*VJ - Built-in voltage (V)
*Rs - Parasitic resistance/Series resistance (Ohms)

.MODEL diode1 D (IS=5.7332E-27 CJO=7.6621E-16 VJ=3.4328E-16 RS=1E-2
* LEAVE THESE PARAMETERS ALONE
                 +level=1 M=0.5 N=1.2)

* This part defines your spice circuit
Vsrc n1 0 DC 0
D1 n1 0 diode1

.DC Vsrc -1.5 1.5 .01

* 516 students: don't forget your transient analysis!
*Vsrc1 n1 0 PULSE()
*.tran 1ns 2us

.PRINT I(D1)
.PROBE I(D1)



.end
