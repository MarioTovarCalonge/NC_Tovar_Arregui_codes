System requirements:

-Python3 (checked with python version 3.6.8)
-Python3 scipy, matplotlib and numpy packages

Once python and its dependencies are available, no explicit installation, or special hardware is needed.


##########################################################################################

Executing the following command from this folder:

python3 IGRA_positive.py

will produce an example of a simulation of a IGRA positive vaccine clinical trial formally equivalent to the one published in the main text. 

-Standard output will be written, by default in this folder.

##########################################################################################

Expected output: as a final result of the script, the algorithm will create two images and several files:

-Igra_positive_VEdis.png show the balance curves between vaccine parameters for this kind of trial for different balances of initial 
population in states F and L.

-curve_L()_F().txt include all the points of the above curves for initial contitions of L() and F().
for instance, curve_L95.00_F5.00.txt contains the dataset for ploting L=95%, F=5% curve of Igra_positive images.

