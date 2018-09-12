# Fractal dimension calculation
The fractal dimension program calculates the fractal dimension for molecules. The main class is
com.actelion.research.chem.properties.fractaldimension.FractalDimensionMoleculeMain.java

Input is a file with a list of SMILES strings. Output is a tab separated file with a header line. One line is produced
for each SMILES string in the result file.

## Java command line example to run the calculation on Linux
java -server -Xmx1g -classpath openChemLib.jar com.actelion.research.chem.properties.fractaldimension.FractalDimensionMoleculeMain -i /home/user/data/adamantane.smi -w /home/user/tmp 

