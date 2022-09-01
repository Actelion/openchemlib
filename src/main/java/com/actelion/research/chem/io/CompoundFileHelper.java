/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.chem.io;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.MolfileParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public abstract class CompoundFileHelper {
	public static final int cFileTypeMask = 0x007FFFFF;
	public static final int cFileTypeDataWarrior = 0x00000001;
	public static final int cFileTypeDataWarriorTemplate = 0x00000002;
	public static final int cFileTypeDataWarriorQuery = 0x00000004;
	public static final int cFileTypeDataWarriorMacro = 0x00000008;
	public static final int cFileTypeTextTabDelimited = 0x00000010;
    public static final int cFileTypeTextCommaSeparated = 0x00000020;
	public static final int cFileTypeTextSemicolonSeparated = 0x00000040;
	public static final int cFileTypeTextVLineSeparated = 0x00000080;
	public static final int cFileTypeTextAnyCSV = cFileTypeTextCommaSeparated | cFileTypeTextSemicolonSeparated | cFileTypeTextVLineSeparated;
    public static final int cFileTypeTextAny = cFileTypeTextTabDelimited | cFileTypeTextAnyCSV;
	public static final int cFileTypeSDV3 = 0x00000100;
    public static final int cFileTypeSDV2 = 0x00000200;
    public static final int cFileTypeSD = cFileTypeSDV3 | cFileTypeSDV2;
	public static final int cFileTypeRXN = 0x00000400;
	public static final int cFileTypeSOM = 0x00000800;
	public static final int cFileTypeJPG = 0x00001000;
	public static final int cFileTypeGIF = 0x00002000;
	public static final int cFileTypePNG = 0x00004000;
	public static final int cFileTypeSVG = 0x00008000;
	public static final int cFileTypePictureFile = cFileTypeJPG | cFileTypeGIF | cFileTypePNG | cFileTypeSVG;
    public static final int cFileTypeRDV3 = 0x00010000;
    public static final int cFileTypeRDV2 = 0x00020000;
    public static final int cFileTypeRD = cFileTypeRDV3 | cFileTypeRDV2;
	public static final int cFileTypeMOL = 0x00040000;
	public static final int cFileTypeMOL2 = 0x00080000;
	public static final int cFileTypePDB = 0x00100000;
	public static final int cFileTypeMMTF = 0x00200000;
	public static final int cFileTypeProtein = cFileTypePDB | cFileTypeMMTF;
	public static final int cFileTypeSDGZ = 0x00400000;
    public static final int cFileTypeUnknown = -1;
	public static final int cFileTypeDirectory = -2;

	public static final int cFileTypeCompoundFiles =
			  CompoundFileHelper.cFileTypeMOL
			| CompoundFileHelper.cFileTypeMOL2
			| CompoundFileHelper.cFileTypeSD
			| CompoundFileHelper.cFileTypeDataWarrior;

	// explicitly supported compression format (SD-files only)
	public static final String cGZipExtention = ".gz";

	public static final int cFileTypeDataWarriorCompatibleData = cFileTypeDataWarrior | cFileTypeTextAny | cFileTypeRD | cFileTypeSD | cFileTypeSDGZ;
	public static final int cFileTypeDataWarriorTemplateContaining = cFileTypeDataWarrior | cFileTypeDataWarriorQuery | cFileTypeDataWarriorTemplate;

	private static File sCurrentDirectory;
	private int mRecordCount,mErrorCount;

	public abstract String selectOption(String message, String title, String[] option);
	public abstract File selectFileToOpen(String title, int filetypes);
	public abstract String selectFileToSave(String title, int filetype, String newFileName);
	public abstract void showMessage(String message);

	public static File getCurrentDirectory() {
		return sCurrentDirectory;
		}

	public static void setCurrentDirectory(File d) {
		sCurrentDirectory = d;
		}

	public ArrayList<StereoMolecule> readStructuresFromFile(boolean readIdentifier) {
        File file = selectFileToOpen("Please select a compound file", cFileTypeCompoundFiles);
        return readStructuresFromFile(file, readIdentifier);
	    }

	public ArrayList<String> readIDCodesFromFile() {
        File file = selectFileToOpen("Please select a compound file", cFileTypeCompoundFiles);
        return readIDCodesFromFile(file);
	    }

	public ArrayList<StereoMolecule> readStructuresFromFile(File file, boolean readIdentifier) {
        if (file == null)
            return null;

        ArrayList<StereoMolecule> moleculeList = new ArrayList<StereoMolecule>();
        readChemObjectsFromFile(file, moleculeList, null, null, readIdentifier, false);

        return moleculeList;
	    }

	/**
	 * Reads all compounds as idcode list from the given file.
	 * @param file MOL-, mol2-, SD- or DataWarrior file
	 * @return
	 */
	public ArrayList<String> readIDCodesFromFile(File file) {
        if (file == null)
            return null;

        ArrayList<String> idcodeList = new ArrayList<String>();
        readChemObjectsFromFile(file, null, idcodeList, null, false, false);

        return idcodeList;
	    }

	/**
	 * Reads all compounds as idcode list with identifiers from the given file.
	 * Therefore, it asks for an identifier column.
	 * @param file if null the user is asked for a file
	 * @param readIDCoords if true, then the id-coords are SPACE delimited attached to the idcode
	 * @return list of String[2] with idcode (index 0) and molecule name (index 1)
	 */
	public ArrayList<String[]> readIDCodesWithNamesFromFile(File file, boolean readIDCoords) {
		if (file == null)
			file = selectFileToOpen("Please select substance file",
									CompoundFileHelper.cFileTypeMOL
								  | CompoundFileHelper.cFileTypeSD
								  | CompoundFileHelper.cFileTypeDataWarrior);

		if (file == null)
			return null;

		ArrayList<String[]> idcodeWithIDList = new ArrayList<String[]>();
		readChemObjectsFromFile(file, null, null, idcodeWithIDList, false, readIDCoords);

		return idcodeWithIDList;
		}

	private void readChemObjectsFromFile(File file,
                                         ArrayList<StereoMolecule> moleculeList,
	                                     ArrayList<String> idcodeList,
										 ArrayList<String[]> idcodeWithIDList,
	                                     boolean readIdentifier, boolean readIDCoords) {
	    mRecordCount = 0;
	    mErrorCount = 0;
	    String filename = file.getName();
	    int index = filename.indexOf('.');
	    String extention = (index == -1) ? "" : filename.substring(index).toLowerCase();

	    if (extention.equals(".mol")
		 || extention.equals(".mol2")) {
		    StereoMolecule mol = null;
	    	if (extention.equals(".mol"))
			    mol = new MolfileParser().getCompactMolecule(file);
	    	else
			    try { mol = new Mol2FileParser().load(filename); } catch (Exception e) { e.printStackTrace(); }

		    if (mol != null && mol.getAllAtoms() != 0) {
			    if (moleculeList != null)
			        moleculeList.add(mol);
			    if (idcodeList != null || idcodeWithIDList != null) {
			        Canonizer canonizer = new Canonizer(mol);
			        String idcode = canonizer.getIDCode();
				    String coords = canonizer.getEncodedCoordinates();
			        if (idcode != null && coords.length() != 0 && readIDCoords)
			            idcode = idcode+" "+coords;
			        if (idcodeList != null)
			            idcodeList.add(idcode);
			        if (idcodeWithIDList != null) {
					    String[] idcodeWithID = new String[2];
					    idcodeWithID[0] = idcode;
					    idcodeWithID[1] = mol.getName();
					    idcodeWithIDList.add(idcodeWithID);
					    }
				    }
			    }
	    	return;
		    }

	    CompoundFileParser parser = (extention.equals(".sdf")) ?
	                                           new SDFileParser(file)
	                              : (extention.equals(".dwar")) ?
	                                           new DWARFileParser(file)
	                              : (extention.equals(".ode")) ?
	                                           new ODEFileParser(file) : null;

	    // If we create molecules,
	    // then we might set the name field with the proper identifier
	    int indexOfID = -1;
	    if (idcodeWithIDList != null || readIdentifier) {
	        String[] fieldNames = parser.getFieldNames();
	        if (fieldNames != null && fieldNames.length != 0) {
	            String id = selectOption("Select compound name or identifier", filename, fieldNames);
	            if (id != null)
	            	for (int i=0; i<fieldNames.length; i++)
	            		if (fieldNames[i].equals(id))
	            			{ indexOfID = i; break; }
	            }
	        if (parser instanceof SDFileParser)
	            parser = new SDFileParser(file, fieldNames);
	        }

	    while (parser.next()) {
	        mRecordCount++;
	        boolean isError = false;
	        if (moleculeList != null) {
	            StereoMolecule mol = parser.getMolecule();
	            if (mol != null) {
	                if (indexOfID != -1)
	                    mol.setName(parser.getFieldData(indexOfID));
                    moleculeList.add(mol);
                    }
	            else {
	                isError = true;
	                }
	            }

	        if (idcodeList != null || idcodeWithIDList != null) {
	            String idcode = parser.getIDCode();
	            if (idcode != null) {
					if (readIDCoords) {
						String coords = parser.getCoordinates();
						if (coords != null)
							idcode = idcode.concat(" ").concat(coords);
						}
					}

				String id = null;
				if (idcodeWithIDList != null) {
					id = parser.getMoleculeName();	// default
					if (indexOfID != -1)
						id = parser.getFieldData(indexOfID);
					}

				if (idcode != null) {
					if (idcodeList != null)
						idcodeList.add(idcode);
					if (idcodeWithIDList != null) {
						String[] idcodeWithID = new String[2];
						idcodeWithID[0] = idcode;
						idcodeWithID[1] = id;
						idcodeWithIDList.add(idcodeWithID);
						}
					}
	            else
                    isError = true;
    	        }

	        if (isError)
	            mErrorCount++;
    	    }
    	}

    public int getRecordCount() {
        return mRecordCount;
        }

	public int getErrorCount() {
	    return mErrorCount;
	    }
	
	public static CompoundFileFilter createFileFilter(int filetypes, boolean isSaving) {
		if (filetypes == cFileTypeDirectory)
			return new CompoundFileFilter() {
				@Override
				public boolean accept(File f) {
					return f.isDirectory();
				}
			};

		CompoundFileFilter filter = new CompoundFileFilter();
		if ((filetypes & cFileTypeDataWarrior) != 0) {
            filter.addExtension("dwar");
            if (!isSaving)
                filter.addExtension("ode");  // old extention
			filter.addDescription("DataWarrior data files");
			}
		if ((filetypes & cFileTypeDataWarriorTemplate) != 0) {
            filter.addExtension("dwat");
            if (!isSaving)
                filter.addExtension("odt");  // old extention
			filter.addDescription("DataWarrior template files");
			}
		if ((filetypes & cFileTypeDataWarriorQuery) != 0) {
            filter.addExtension("dwaq");
            if (!isSaving)
                filter.addExtension("odq");  // old extention
			filter.addDescription("DataWarrior query files");
			}
		if ((filetypes & cFileTypeDataWarriorMacro) != 0) {
            filter.addExtension("dwam");
			filter.addDescription("DataWarrior macro files");
			}
		if ((filetypes & cFileTypeTextTabDelimited) != 0) {
			filter.addExtension("tsv");
			filter.addExtension("txt");
			filter.addDescription("TAB delimited text files");
			}
        if ((filetypes & cFileTypeTextAnyCSV) != 0) {
            filter.addExtension("csv");
            filter.addDescription("Comma [,;|] separated text files");
            }
		if ((filetypes & cFileTypeRXN) != 0) {
			filter.addExtension("rxn");
			filter.addDescription("MDL reaction files");
			}
		if ((filetypes & cFileTypeSD) != 0) {
			filter.addExtension("sdf");
			filter.addDescription("MDL SD-files");
			}
		if ((filetypes & cFileTypeSDGZ) != 0) {
			filter.addExtension("sdf.gz");
			filter.addDescription("gzipped MDL SD-files)");
			}
		if ((filetypes & cFileTypeRD) != 0) {
			filter.addExtension("rdf");
			filter.addDescription("MDL RD-files");
			}
		if ((filetypes & cFileTypeSOM) != 0) {
            filter.addExtension("dwas");
            if (!isSaving)
                filter.addExtension("som");  // old extention
			filter.addDescription("DataWarrior self organized map");
			}
		if ((filetypes & cFileTypeJPG) != 0) {
			filter.addExtension("jpg");
			filter.addExtension("jpeg");
			filter.addDescription("JPEG image files");
			}
		if ((filetypes & cFileTypeGIF) != 0) {
			filter.addExtension("gif");
			filter.addDescription("GIF image files");
			}
		if ((filetypes & cFileTypePNG) != 0) {
			filter.addExtension("png");
			filter.addDescription("PNG image files");
			}
		if ((filetypes & cFileTypeSVG) != 0) {
			filter.addExtension("svg");
			filter.addDescription("scalable vector graphics files");
			}

        if (filetypes == cFileTypeDataWarriorCompatibleData) {
            filter.setDescription("DataWarrior compatible files");
            }
        if (filetypes == cFileTypeDataWarriorTemplateContaining) {
            filter.setDescription("Files containing a DataWarrior template");
            }
        if (filetypes == cFileTypePictureFile) {
            filter.setDescription("Image files");
            }
		if ((filetypes & cFileTypePDB) != 0) {
			filter.addExtension("pdb");
			filter.addDescription("Protein Data Bank files");
			}
		if ((filetypes & cFileTypeMMTF) != 0) {
			filter.addExtension("mmtf");
			filter.addDescription("Binary Protein Data Bank files");
			}
		if ((filetypes & cFileTypeMOL) != 0) {
			filter.addExtension("mol");
			filter.addDescription("MDL Molfiles");
			}
		if ((filetypes & cFileTypeMOL2) != 0) {
			filter.addExtension("mol2");
			filter.addDescription("Tripos Mol2 files");
			}

		return filter;
		}

	/**
	 * Return the extension portion of the file's name.
	 * Known joint extensions (currently only ".sdf.gz") are returned where they exist.
	 * @return extension without the dot
	 */
	public static String getExtension(File file) {
		int index = -1;
		String filename = (file == null) ? null : file.getName();
		if (filename != null)
			index = getExtensionIndex(filename);

		return (index == -1) ? null : filename.substring(index+1).toLowerCase();
		}

	/**
	 * Returns the index of the dot separating filename from its extension.
	 * Known joint extensions (currently only ".sdf.gz") are recognized.
	 * @param filename
	 * @return index of extension dot or -1
	 */
	private static int getExtensionIndex(String filename) {
		int i = filename.lastIndexOf('.');

		if (i>0 && i<filename.length()-1
		 && filename.substring(i).equalsIgnoreCase(CompoundFileHelper.cGZipExtention))
			i = filename.lastIndexOf('.', i-1);

		return (i>0 && i<filename.length()-1) ? i : -1;
		}

	/**
	 * Provided that fileName has a leading file path, then path and separator are removed.
	 * Provided that fileName has a recognized extension, then the extension is removed.
	 * @param filePath file name with or without complete path and with or without extension
	 * @return naked file name without leading path and extension
	 */
	public static String removePathAndExtension(String filePath) {
		int i1 = filePath.lastIndexOf(File.separatorChar);
		int i2 = (getFileType(filePath) != cFileTypeUnknown) ? getExtensionIndex(filePath) : -1;
		if (i1 == -1)
			return (i2 == -1) ? filePath : filePath.substring(0, i2);
		else
			return (i2 == -1 || i2 < i1) ? filePath.substring(i1+1) : filePath.substring(i1+1, i2);
		}

	/**
	 * Provided that filePath has a recognized extension, then the extension is removed.
	 * @param filePath file name with or without complete path and with or without extension
	 * @return file name or path without extension
	 */
	public static String removeExtension(String filePath) {
		int i = (getFileType(filePath) != cFileTypeUnknown) ?
				filePath.lastIndexOf('.') : -1;
		return (i == -1) ? filePath : filePath.substring(0, i);
		}

	/**
	 * Note: If
	 * @param filename
	 * @return one or multiple filtetypes that matching the extension of the given filename
	 */
	public static int getFileType(String filename) {
        int index = getExtensionIndex(filename);

        if (index == -1)
            return cFileTypeUnknown;

        String extension = filename.substring(index).toLowerCase();
        if (extension.equals(".dwar") || extension.equals(".ode"))
            return cFileTypeDataWarrior;
        if (extension.equals(".dwat") || extension.equals(".odt"))
            return cFileTypeDataWarriorTemplate;
        if (extension.equals(".dwaq") || extension.equals(".odq"))
            return cFileTypeDataWarriorQuery;
        if (extension.equals(".dwas") || extension.equals(".som"))
            return cFileTypeSOM;
        if (extension.equals(".dwam"))
            return cFileTypeDataWarriorMacro;
        if (extension.equals(".txt") || extension.equals(".tsv"))
            return cFileTypeTextTabDelimited;
        if (extension.equals(".csv"))
            return cFileTypeTextAnyCSV;
        if (extension.equals(".sdf"))
            return cFileTypeSD;
		if (extension.equals(".sdf.gz"))
			return cFileTypeSDGZ;
        if (extension.equals(".rdf"))
            return cFileTypeRD;
        if (extension.equals(".rxn"))
            return cFileTypeRXN;
        if (extension.equals(".jpg") || extension.equals(".jpeg"))
            return cFileTypeJPG;
		if (extension.equals(".gif"))
			return cFileTypeGIF;
        if (extension.equals(".png"))
            return cFileTypePNG;
        if (extension.equals(".svg"))
            return cFileTypeSVG;
		if (extension.equals(".mol"))
			return cFileTypeMOL;
		if (extension.equals(".mol2"))
			return cFileTypeMOL2;
		if (extension.equals(".pdb"))
			return cFileTypePDB;
		if (extension.equals(".mmtf"))
			return cFileTypeMMTF;

        return cFileTypeUnknown;
        }

    /**
     * @param fileTypes
     * @return list of all extensions that are covered by fileTypes
     */
    public ArrayList<String> getExtensionList(int fileTypes) {
    	ArrayList<String> list = new ArrayList<String>();
    	int type = 0x00000001;
    	while ((type & cFileTypeMask) != 0) {
    		if ((type & fileTypes) != 0)
    			for (String extension:getExtensions(type))
    			    if (!list.contains(extension))
    				    list.add(extension);

    		type <<= 1;
    		}
    	return list;
    	}

	/**
	 * @param filetype
	 * @return preferred file extension including the dot
	 */
	public static String getExtension(int filetype) {
    	String[] extensions = getExtensions(filetype);
    	return extensions.length == 0 ? "" : extensions[0];
		}

	/**
	 * @param filetype
	 * @return file extensions including the dot
	 */
    public static String[] getExtensions(int filetype) {
		ArrayList<String> extensions = new ArrayList<>();
		switch (filetype) {
		case cFileTypeDataWarrior:
			extensions.add(".dwar");
			break;
		case cFileTypeDataWarriorQuery:
			extensions.add(".dwaq");
			break;
		case cFileTypeDataWarriorTemplate:
			extensions.add(".dwat");
			break;
		case cFileTypeDataWarriorMacro:
			extensions.add(".dwam");
			break;
		case cFileTypeTextTabDelimited:
			extensions.add(".txt");
			extensions.add(".tsv");
			break;
		case cFileTypeTextAnyCSV:
        case cFileTypeTextCommaSeparated:
		case cFileTypeTextSemicolonSeparated:
		case cFileTypeTextVLineSeparated:
            extensions.add(".csv");
            break;
		case cFileTypeSD:
        case cFileTypeSDV2:
        case cFileTypeSDV3:
			extensions.add(".sdf");
			break;
		case cFileTypeRD:
        case cFileTypeRDV2:
        case cFileTypeRDV3:
			extensions.add(".rdf");
			break;
		case cFileTypeRXN:
			extensions.add(".rxn");
			break;
		case cFileTypeSOM:
			extensions.add(".dwas");
			break;
		case cFileTypeJPG:
			extensions.add(".jpeg");
			extensions.add(".jpg");
			break;
		case cFileTypeGIF:
			extensions.add(".gif");
			break;
		case cFileTypePNG:
			extensions.add(".png");
			break;
		case cFileTypeSVG:
			extensions.add(".svg");
			break;
		case cFileTypeMOL:
			extensions.add(".mol");
			break;
		case cFileTypeMOL2:
			extensions.add(".mol2");
			break;
		case cFileTypePDB:
			extensions.add(".pdb");
			break;
		case cFileTypeMMTF:
			extensions.add(".mmtf");
			break;
		case cFileTypeSDGZ:
			extensions.add(".sdf.gz");
			break;
			}
		return extensions.toArray(new String[0]);
		}

	public void saveRXNFile(Reaction rxn) {
		String fileName = selectFileToSave("Select reaction file", cFileTypeRXN, "Untitled Reaction");
		if (fileName != null) {
			String extension = ".rxn";
			int dotIndex = fileName.lastIndexOf('.');
			int slashIndex = fileName.lastIndexOf(File.separator);
			if (dotIndex == -1
			 || dotIndex < slashIndex)
				fileName = fileName.concat(extension);
		    else if (!fileName.substring(dotIndex).equalsIgnoreCase(extension)) {
				showMessage("uncompatible file name extension.");
			    return;
				}

			try {
				BufferedWriter theWriter = new BufferedWriter(new FileWriter(new File(fileName)));
				new RXNFileCreator(rxn).writeRXNfile(theWriter);
				theWriter.close();
				}
			catch (IOException e) {
				showMessage("IOException: "+e);
				}
			}
		}
	}
