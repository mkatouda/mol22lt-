# mol2tolt.py


## Overview
Prepare Moltemplate LT file (*.lt) from antechamber output mol2 file (*.mol2)
or acpype output Gromacs gro file (*.gro) and itp file (*.itp). 
mol2tolt.py works with Python 3.


## Usage example
### Mol2 -> Moltemplate
#### Input files
Mol2 file: rr-dihydroxyoctane.mol2  

Note: Atom type of mol2 file (*.mol2) should be assigned by antechamber.  

#### Command
```
$ python ./gaff-itp2moltemp.py rr-dihydroxyoctane.mol2
```

#### Output file
Moltemplate LT file: rr-dihydroxyoctane.lt  
Moltemplate LT file with removing first atom: rr-dihydroxyoctane-a.lt  
Moltemplate LT file with removing last atom:  rr-dihydroxyoctane-b.lt  
Moltemplate LT file with removing first and last atoms: rr-dihydroxyoctane-c.lt  

### Gromacs -> Moltemplate
#### Input files
Gromacs gro file: rr-dihydroxyoctane_GMX.gro  
Gromacs itp file: rr-dihydroxyoctane_GMX.itp    

#### Command
```
$ python ./gaff-itp2moltemp.py rr-dihydroxyoctane_GMX.gro
```

#### Output file
Moltemplate LT file: rr-dihydroxyoctane.lt  
Moltemplate LT file with removing first atom: rr-dihydroxyoctane-a.lt  
Moltemplate LT file with removing last atom:  rr-dihydroxyoctane-b.lt  
Moltemplate LT file with removing first and last atoms: rr-dihydroxyoctane-c.lt  

## License
The source code is licensed MIT. The website content is licensed CC BY 4.0,see LICENSE.
