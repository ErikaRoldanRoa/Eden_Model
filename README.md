## Topology and Geometry of the Eden Model

DEVELOPERS: <br />
Fedor Manin (manin@math.ucsb.edu) <br />
Anna Krymova (anna.krymova@tum.de) <br />
Erika Roldan (erika.roldan@ma.tum.de) <br />
Benjamin Schweinhart (bschweinhart@albany.edu) <br />



DATE: December 15, 2020

LICENSE: GNU GENERAL PUBLIC LICENSE (see license.txt)

## Overview 
This package includes software to run simulations of the Eden Growth Model in Z^d in dimensions 2-5, and analyze the topology (Betti numbers and persistent homology) and local geometry of the structure. The software is also able to read a .txt file with the time at which tiles should be added. This allows analysis of simulations from other stochastic models. Note: the cells in a .txt file should be in chronological order. For graphical representation, the program can create a picture of a two-dimensional growth model, and can output a .txt file for 3-dimensional growth modles which can be inputed to MAYA to produce an interactive 3-dimensional image.

GUDHI is used to compute homology and persistent homology in 3D, 4D, and 5D. If you use this functionality, make sure to cite this library.

To represent the topology and local geometry of the Eden growth model, the software can build plots showing the following:
* the frequencies of the changes in Betti numbers (Figure 6.1 in the paper),
* the distributionvolumes of top dimensional "holes" (the finite components of the complement of the structure due to Alexander duality, Figure 6.5 in the paper),
* the growth of the Betti numbers and the perimeter (Figure 6.2 in the paper),
* in two and three dimensions, the frequencies of top dimensional holes with specific shapes with 3 and 4 cells (Table 6.4 in the paper) 

The folder *2d/files* contains the data of the 2D simulations that were used in the paper TOPOLOGY AND LOCAL GEOMETRY OF THE EDEN MODEL https://arxiv.org/pdf/2005.12349.pdf.

All plots and graphs are saved in the project folder.

## Acknowledgments
Fedor Manin supported in part by NSF DMS-2001042. <br />
Erika Roldan was supported in part by NSF-DMS #1352386 and NSF-DMS #1812028 during 2018-2019. <br />
This project received funding from the European Union’s Horizon 2020 research and innovation program under the
Marie Skłodowska-Curie grant agreement No. 754462.


## Citations 

If you use this code, cite the paper TOPOLOGY AND LOCAL GEOMETRY OF THE EDEN MODEL. If you use the computations of homology and persistent homology for 3D-5D, cite the GUDHI package.

```
@article{manin2020topology,
   title={Topology and local geometry of the Eden model},
   author={Manin, Fedor and Roldan, Erika and Schweinhart, Benjamin},
   journal={arXiv preprint arXiv:2005.12349},
   year={2020}
}
```

## Dependencies:

Python 3.8.

GUDHI. http://gudhi.gforge.inria.fr/


## Installation

To install this package with conda run one of the following:
```
conda install -c conda-forge gudhi
```
```
conda install -c conda-forge/label/cf201901 gudhi
```
```
conda install -c conda-forge/label/cf202003 gudhi
```

## Usage
Download the whole folder and run the file *main.py*.<br />
At first, the system will ask you to enter the dimension:
```
Please, enter the desired dimension of your model (from 2 to 5): 
```
Then, you have to specify if you have a file with a pre-defined model or not:
```
Do you have a file with a model? 
0 -- you want to generate a new model 
1 -- you have a file
```
In the case of а 2d or а 3d model, on the next step you decide if you want a picture:
```
Do you want a picture of your model? (with a large model it can take time)  
0 -- no 
1 -- yes
```
In case you want to read the cubical complex from a file, the program asks you to specify the file format:
```
What is the format of the file? 
0 -- list of tuples 
1 -- Perseus
```
And then you should give a name of a file:
```
Name of the file (for example, filename.txt):
```
If you chose to generate a new model, then now the system asks you to enter the size of the model:
```
How many tiles would you like in your model?
```
After that, the modeling and analysis is complete. 
Before creating barcodes, the program will ask you to specify the minimum length of the interval for every barcode (depending on the dimension of the model).
It is done in order to omit short-lived homology groups, i.e. short intervals.

When all calculations are finished, you will see the sentence:
```
WE ARE DONE! CHECK THE FOLDER!
```
Now, you are welcome to check the results in the corresponding folder. 



