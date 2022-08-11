## OpenChemLib
*OpenChemLib* is Java based framework providing cheminformatics core functionality and user interface components. Its main focus is on organics chemistry and small molecules. It is built around a StereoMolecule class, which represents a molecule using atom and bond tables, provides atom neighbours, ring and aromaticity information, and supports MDL's concept of enhanced stereo representation. Additional classes provide, 2D-depiction, descriptor calculation, molecular similarity and substructure search, reaction search, property prediction, conformer generation, support for molfile and SMILES formats, energy minimization, ligand-protein interactions, and more. *OpenChemLib's idcode* represents molecules, fragments or reactions as canonical, very compact string that includes stereo and query features.
Different to other cheminformatics frameworks, *OpenChemLib* also provides user interface components that allow to easily embed chemical functionality into Java applications, e.g. to display or edit chemical structures or reactions.

### Dependencies
*OpenChemLib* requires JRE 8 or newer including JavaFX. Otherwise, there are no dependencies.

### How to download the project
```bash
git clone https://github.com/Actelion/openchemlib.git
```

### Build the project
To build the project with maven run the following from within the project directory:
```bash
./mvnw package
```
To build the project using the JDK only (Mac, Linux) run this from within the project directory:
```
./buildOpenChemLib
```

### Folder 'examples'
Contains examples for working with the *OpenChemLib* library.

### Logo
![logo](logo.png)
