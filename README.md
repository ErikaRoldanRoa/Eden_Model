# Topology and local geometry of the Eden Cell Growth Model

DEVELOPERS: 
Fedor Manin (manin@math.ucsb.edu)
Anna Krymova (anna.krymova@tum.de)
Erika Roldan (erika.roldan@ma.tum.de)
Benjamin Schweinhart (bschweinhart@albany.edu)

DATE: December 15, 2020

LICENSE: GNU GENERAL PUBLIC LICENSE (see license.txt)

# Overview 
This library allow you to run simulations of the Eden Growth Model for dimensions 2-5 and analyze the topology (Betti numbers and persistent homology) and local geometry of the structure. The software is also able to read a txt file with the time at which tiles should be added. This allows to analize simulations comming from other stochastic models comming from First Passage Percolation. For the 2-dimensional model, the program can create a picture of it. For the 3-dimensional model, a txt file is created: the user has to copy-paste its content to MAYA to get the 3-dimensional model of the Eden Model.

For computing homology and persistent homology in 3D, 4D and 5D we use the library GUDHI, make sure to cite this library if you are using homology or persistent homology in this dimensions.

The software builds:
* the plot of the frequencies of Betti numbers, 
* the plot of the frequency of the volume of top dimensional "holes" (the finite components of the complement of the structure due to Alexander duality),
* the plot of the growth rates of Betti numbers and the perimeter,
* for 2d and 3d model, the plots of the frequency of the number of top dimensional holes for specific shapes with 3 and 4 cells

The folder ** contains the data of the 2D simulations that were used in the paper TOPOLOGY AND LOCAL GEOMETRY OF THE EDEN MODEL.

All the plots and graphs are saved in the folder of the project.

# Acknoledgments
Fedor Manin supported in part by NSF DMS-2001042.
Erika Roldan was supported in part by NSF-DMS #1352386 and NSF-DMS #1812028 during 2018-2019. 
This project received funding from the European Union’s Horizon 2020 research and innovation program under the
Marie Skłodowska-Curie grant agreement No. 754462.


# Citations 

If you use this code, cite the paper. If you use the computations of persistent homology cite the GUDHI package.


DEPENDENCIES:

Python 3.8.

GUDHI. http://gudhi.gforge.inria.fr/

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

What things you need to install the software and how to install them

```
gudhi
```

### Installing

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

## Running the tests

At first, the system will ask you to enter the dimension:
```
Please, enter the desired dimension of your model (from 2 to 5): 
```
Then, you have to specify if you have a file with a pre-defined model or not:
```
Do you have a file with a model? 
 1 -- you have a file 
 0 -- you want to generate a new model
```
In the case of а 2d or а 3d model, on the next step you decide if you want a picture:
```
Do you want a picture of your model? (with a large model it can take time) 
 1 -- yes 
 0 -- no
```
If in the previous step you chose to generate a new model, then now the system asks you to enter the size of the model:
```
How many tiles would you like in your model?
```
After that, all the modeling and analysis is done. 
Before creating barcodes, the program will ask you to specify the minimum length of the interval for every barcode (depending on the dimension of the model).
It is done in order to omit short-lived homology groups, i.e. short intervals.

When all calculations are finished, you will see the sentence:
```
WE ARE DONE! CHECK THE FOLDER!
```
And now you are welcome to check the results in the corresponding folder. 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

