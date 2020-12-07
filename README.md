

DEVELOPERS:
DATE: November 21, 2020

LICENSE:

OVERVIEW:
This software provides 

# Eden_Model
This software provides methods to build and analyze 2d, 3d, 4d, and 5d Eden Models. The software is able to read a pre-defined model from a txt file, as well as create a model based on the number of tiles, specified by the user. For the 2-dimensional model, the program can create a picture of it. For the 3-dimensional model, a txt file is created: the user has to copy-paste its content to MAYA to get the 3-dimensional model of the Eden Model. 
The software builds:
* the plot of the frequencies of Betti numbers, 
* the plot of the Frequency of the Number of Holes,
* the plot of the growth rates of Betti numbers and the perimeter
* for 2d and 3d model, the plots of the Frequency of the Number of Holes for specific hole shapes.
All the plots and graphs are saved in the folder of the project.
Gudhi package is used to build the persistent barcodes.

## Authors


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

