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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.StringReader;

import com.actelion.research.chem.MolfileParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.io.BOMSkipper;


public class RXNFileParser
{
    private static final String RXN_V3_COUNTS_LINE = "M  V30 COUNTS";
    private static final String V30_BEGIN_REACTANT = "M  V30 BEGIN REACTANT";
    private static final String V30_BEGIN_PRODUCT = "M  V30 BEGIN PRODUCT";
    private static final String RXN_MAGIC = "$RXN";
    private static final String RXN_V3_MAGIC = "$RXN V3000";
    private static final String END_CTAB = "M  V30 END CTAB";
    private static final String END_MOL_TAG = "M  END";
    private static final String MOL_MAGIC = "$MOL";
    private static final String AUX_MOLFILE_HEADER = "\nActelion Java MolfileCreator 2.0\n\n  0  0  0  0  0  0              0 V3000\n";
    private static final int DEFAULT_CAPACITY = 32768;

    public RXNFileParser()
	{
	}

	public Reaction getReaction(String buffer) throws Exception
	{
		Reaction theReaction = new Reaction();
		BufferedReader theReader = new BufferedReader(new StringReader(buffer));
		parse(theReaction, theReader);

		return theReaction;
	}

	public Reaction getReaction(File file) throws Exception
	{
		Reaction theReaction = new Reaction();
		BufferedReader theReader = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"));
        BOMSkipper.skip(theReader);
		parse(theReaction, theReader);

		return theReaction;
	}

	public boolean parse(Reaction theReaction, String buffer)
		throws Exception
	{
		BufferedReader theReader = new BufferedReader(new StringReader(buffer));

		return parse(theReaction, theReader);
	}

	public boolean parse(Reaction theReaction, File file)
		throws Exception
	{
		BufferedReader theReader = new BufferedReader(new InputStreamReader(new FileInputStream(file), "UTF-8"));
        BOMSkipper.skip(theReader);

		return parse(theReaction, theReader);
	}

    private boolean parse(Reaction theReaction, BufferedReader theReader) throws Exception
    {
        String theLine = theReader.readLine();
        boolean ok = false;
        if ((theLine == null) || !theLine.startsWith(RXN_MAGIC)) {
            throw new Exception("'$RXN' tag not found");
        }
        if (theLine.equals(RXN_V3_MAGIC)) {
            ok = parseV3(theReaction,theReader);
        } else {
            ok = parseV2(theReaction,theReader);
        }
        return ok;
    }

    // First line is already parsed
    private boolean parseV3(Reaction theReaction, BufferedReader theReader) throws Exception
    {
        String theLine = null;
        for (int i = 0; i < 3; i++) {
            theLine = theReader.readLine();
        }

        // preferrably decode the idcode based encoding from the comment line, if present
        if (theLine.startsWith(RXNFileCreator.RXN_CODE_TAG)) {
            String encoding = theLine.substring(RXNFileCreator.RXN_CODE_TAG.length());
            if (ReactionEncoder.decode(encoding, true, theReaction) != null)
            	return true;
        }

        theLine = theReader.readLine();
        MolfileParser molParser = new MolfileParser();
        if (theLine != null && theLine.startsWith(RXN_V3_COUNTS_LINE)) {
            String t = theLine.substring(13).trim();
            String p[] = t.split(" ");
            int reactantCount = Integer.parseInt(p[0]);
            int productCount = Integer.parseInt(p[1]);
            if (reactantCount > 0) {
                theLine = theReader.readLine();
                if (V30_BEGIN_REACTANT.equals(theLine)) {
                    for (int i = 0; i < reactantCount; i++) {
                        StereoMolecule molecule = new StereoMolecule();
                        StringBuffer molfile = new StringBuffer(DEFAULT_CAPACITY);
                        molfile.append(AUX_MOLFILE_HEADER);
                        do {
                            theLine = theReader.readLine();
                            molfile.append(theLine);
                            molfile.append("\n");
                        } while ((theLine != null) && !theLine.startsWith(END_CTAB));
                        molParser.parse(molecule,molfile);
                        theReaction.addReactant(molecule);
                    }
                }
                
                theLine = theReader.readLine(); //  end reactant line
            }
            if (productCount > 0) {
                theLine = theReader.readLine();
                if (V30_BEGIN_PRODUCT.equals(theLine)) {
                    for (int i = 0; i < productCount; i++) {
                        StereoMolecule molecule = new StereoMolecule();
                        StringBuffer molfile = new StringBuffer(DEFAULT_CAPACITY);
                        molfile.append(AUX_MOLFILE_HEADER);
                        do {
                            theLine = theReader.readLine();
                            molfile.append(theLine);
                            molfile.append("\n");
                        } while ((theLine != null) && !theLine.startsWith(END_CTAB));
                        molParser.parse(molecule,molfile);
                        theReaction.addProduct(molecule);
                    }
                    theLine = theReader.readLine(); //  end product line
                }
            }
            return true;
        }
        return false;
    }

	/* This is the old approach, RXN specs were missing so we invented one
    private boolean parseOld(Reaction theReaction, BufferedReader theReader)
		throws Exception
	{
		String theLine = theReader.readLine();

		if ((theLine == null) || !theLine.startsWith(RXN_MAGIC)) {
			throw new Exception("'$RXN' tag not found");
		}

		for (int i = 0; i < 4; i++) {
			theLine = theReader.readLine();
		}
		int reactantCount = Integer.parseInt(theLine.substring(0, 3).trim());
		int productCount = Integer.parseInt(theLine.substring(3, 6).trim());

		MolfileParser molParser = new MolfileParser();

		for (int i = 0; i < reactantCount; i++) {
			theLine = theReader.readLine();

			if ((theLine == null) || !theLine.startsWith(MOL_MAGIC)) {
				throw new Exception("'$MOL' tag not found");
			}

            StereoMolecule reactant = new StereoMolecule();
			StringBuffer molfile = new StringBuffer(DEFAULT_CAPACITY);

			do {
				theLine = theReader.readLine();
				molfile.append(theLine);
				molfile.append("\n");
			} while ((theLine != null) && !theLine.startsWith(END_MOL_TAG));

			if (theLine == null) {
				throw new Exception("'M  END' not found");
			}

			molParser.parse(reactant, molfile);

			theReaction.addReactant(reactant);
		}

		for (int i = 0; i < productCount; i++) {
			theLine = theReader.readLine();

			if ((theLine == null) || !theLine.startsWith(MOL_MAGIC)) {
				throw new Exception("'$MOL' tag not found");
			}

            StereoMolecule product = new StereoMolecule();
			StringBuffer molfile = new StringBuffer(DEFAULT_CAPACITY);

			do {
				theLine = theReader.readLine();
				molfile.append(theLine);
				molfile.append("\n");
			} while ((theLine != null) && !theLine.startsWith(END_MOL_TAG));

			if (theLine == null) {
				throw new Exception("'M  END' not found");
			}

			molParser.parse(product, molfile);

			theReaction.addProduct(product);
		}

		theReader.close();

		return true;

		//return theReaction;
	}	*/

    /**
     * @param theReaction
     * @param theReader
     * @return
     * @throws Exception
     */
    // First line is already parsed
    private boolean parseV2(Reaction theReaction, BufferedReader theReader) throws Exception
    {
        String theLine = null;
        for (int i = 0; i < 3; i++) {
            theLine = theReader.readLine();
        }

		// preferrably decode the idcode based encoding from the comment line, if present
		if (theLine.startsWith(RXNFileCreator.RXN_CODE_TAG)) {
			String encoding = theLine.substring(RXNFileCreator.RXN_CODE_TAG.length());
			if (ReactionEncoder.decode(encoding, true, theReaction) != null)
				return true;
		}

		theLine = theReader.readLine();
		int reactantCount = Integer.parseInt(theLine.substring(0, 3).trim());
        int productCount = Integer.parseInt(theLine.substring(3, 6).trim());

        MolfileParser molParser = new MolfileParser();

        for (int i = 0; i < reactantCount; i++) {
            theLine = theReader.readLine();

            if ((theLine == null) || !theLine.startsWith(MOL_MAGIC)) {
                throw new Exception("'$MOL' tag not found");
            }

            StereoMolecule reactant = new StereoMolecule();
            StringBuffer molfile = new StringBuffer(DEFAULT_CAPACITY);

            do {
                theLine = theReader.readLine();
                molfile.append(theLine);
                molfile.append("\n");
            } while ((theLine != null) && !theLine.startsWith(END_MOL_TAG));

            if (theLine == null) {
                throw new Exception("'M  END' not found");
            }

            molParser.parse(reactant, molfile);

            theReaction.addReactant(reactant);
        }

        for (int i = 0; i < productCount; i++) {
            theLine = theReader.readLine();

            if ((theLine == null) || !theLine.startsWith(MOL_MAGIC)) {
                throw new Exception("'$MOL' tag not found");
            }

            StereoMolecule product = new StereoMolecule();
            StringBuffer molfile = new StringBuffer(DEFAULT_CAPACITY);

            do {
                theLine = theReader.readLine();
                molfile.append(theLine);
                molfile.append("\n");
            } while ((theLine != null) && !theLine.startsWith(END_MOL_TAG));

            if (theLine == null) {
                throw new Exception("'M  END' not found");
            }

            molParser.parse(product, molfile);

            theReaction.addProduct(product);
        }

        theReader.close();

        return true;

        //return theReaction;
    }

}
