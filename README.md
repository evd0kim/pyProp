# pyProp

Blade Element Momentum theory program to calculate small-scale low-Re rotor characteristics.

I.Evdokimov, "Blade element theory in the UAV multirotor Blade Optimization", 
CAE Conference Pacengo del Garda (Verona) - Italy, 19 - 20 October 2015
http://www.caeconference.com/speakers/evdokimov.html

PYTHON MODULES

STANDARD: csv, sys

SCIENTIFIC: numpy, scipy, pylab, matplotlib

==================================================================================================

HOW TO USE:

==================================================================================================

Launch in shell:

	python pyBladeCalc.py APCFOIL.csv apcsf_10x7_geom.txt

APCFOIL.csv - file with characteristics of the airfoil.
The data should be formatted as it shows in the first line:	
	RE;ALPHA;CX;CY;MZ
Chatacteristics of the test blade were obtained using XFOIL programm. For airfoil very similar to
Eppler E63.

apcsf_10x7_geom.txt - the standard file with geometrical characteristics of the blade downloaded 
from database [1]

Launch in python shell:

1. Importing sys - module:

	import sys
	
2. Assign promt parameters to list:

   	sys.argv = ['', 'APCFOIL.csv', 'apcsf_10x7_geom.txt']
	
3. Launching script

   	execfile('pyBladeCalc.py')

==================================================================================================

RESULTS:

==================================================================================================

The output file format of the pyBladeCalc.py is like files with experimental data in database [1].
It can be successfully plotted using GNUPLOT or similar software.

The result files of the pyBladeCalcStatic.py is different. The static data is obtained through the
different rotation frequencies. In according with [1] the user has to compare only Cp and Ct in
static mode, but [2] allows to calculate efficiency at hover so the output file includes this 
characteristic.

The script pyBladeCalc.py also produces *.png files with resulting angles and Cp, Ct and 
efficiency characteristics. 

The pyBladeCalcStatic.py code icludes commented pieces with printing commands and plots
resulting angles using matplotlib package.

==================================================================================================

REFERENCES:

==================================================================================================

1. http://m-selig.ae.illinois.edu/props/volume-1/propDB-volume-1.html 
2. (Russian) Jur'ev B. N.: “Ajerodinamicheskij raschet vertoletov”, Oborongnz, 1956
