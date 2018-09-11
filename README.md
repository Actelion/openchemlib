# OpenChemLib
# Open Source Chemistry Library

## Fractal dimension
The fractal dimension program calculates the fractal dimension for molecules. The main class is
com.actelion.research.chem.properties.fractaldimension.FractalDimensionMoleculeMain.java

Input is a file with a list of SMILES strings. Output is a tab separated file with a header line. One line is produced
for each SMILES string in the result file.

### SMILES code for five molecules
C(C(C1)C2)C3CC2CC1C3 <br/>
O=C(C[C@@H]1OCC=C2[C@H](C3)[C@@H]1[C@H]1[C@@]4(CC5)[C@H]3N5C2)N1c1c4cccc1 <br/>
Cc(ccc(Cl)c1)c1N(CC1)CCN1C(C(C1)COc(cc2)c1cc2OC)=O <br/>
CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCC(NCCO)=O <br/>
OCC(C(C(C1O)O)O)OC1O <br/>

### Java command line example to run the calculation on Linux
java -server -Xmx1g -classpath openChemLib.jar com.actelion.research.chem.properties.fractaldimension.FractalDimensionMoleculeMain $*

