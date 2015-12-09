# pyProp

Blade Element Momentum theory program to calculate small-scale low-Re rotor characteristics.

I.Evdokimov, "Blade element theory in the UAV multirotor Blade Optimization", 
CAE Conference Pacengo del Garda (Verona) - Italy, 19 - 20 October 2015
http://www.caeconference.com/speakers/evdokimov.html

==================================================================================================

HOW TO USE:

==================================================================================================

In shell:
	python pyBladeCalc.py APCFOIL.csv apcsf_10x7_geom.txt

APCFOIL.csv - file with characteristics of airfoil.
The data should be formatted as it shows in the first line:	
	RE;ALPHA;CX;CY;MZ
Chatacteristics of the test blade were obtained using XFOIL programm. For airfoil very similar to
Eppler E63.

apcsf_10x7_geom.txt - the standard file with geometrical characteristics of the blade downloaded 
from database http://m-selig.ae.illinois.edu/props/volume-1/propDB-volume-1.html 

==================================================================================================

RESULTS:

==================================================================================================

