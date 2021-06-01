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

/*
 * @history
 * 27.01.2005    cxr     added scaling before writing individual molfiles!
 * 16.08.2010    cxr     adapted to be consistent with symyx ct-file format 2010
 */

package com.actelion.research.chem.io;

import com.actelion.research.chem.MolfileV3Creator;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;

public class RXNFileV3Creator
{
	public static final String RXN_CODE_TAG = RXNFileCreator.RXN_CODE_TAG;
    private static final String NL = System.lineSeparator();

    private StringBuffer rxnbuffer = null;

    public RXNFileV3Creator(Reaction rxn)
    {
        this(rxn, null);
    }

	public RXNFileV3Creator(Reaction r, String programName) {
        Reaction rxn = new Reaction(r);
        try {
            StringWriter theWriter = new StringWriter();
            theWriter.write("$RXN V3000"+NL);
            theWriter.write(programName != null ? programName : "");
            theWriter.write(NL+NL);
			theWriter.write(RXN_CODE_TAG+ReactionEncoder.encode(r, true,
                    ReactionEncoder.INCLUDE_MAPPING | ReactionEncoder.INCLUDE_COORDS | ReactionEncoder.INCLUDE_CATALYSTS));
			theWriter.write(NL);
            int rcnt = rxn.getReactants();
            int pcnt = rxn.getProducts();
            theWriter.write(String.format("M  V30 COUNTS %d %d"+NL,rcnt,pcnt));

			double scalingFactor = getScalingFactor(rxn);

            if (rcnt > 0) {
                theWriter.write("M  V30 BEGIN REACTANT"+NL);
                for (int i=0; i<rxn.getReactants(); i++) {
                theWriter.write(MolfileV3Creator.writeCTAB(rxn.getReactant(i), scalingFactor));
                }
                theWriter.write("M  V30 END REACTANT"+NL);
            }
            if (pcnt > 0) {
                theWriter.write("M  V30 BEGIN PRODUCT"+NL);
                for (int i=0; i<rxn.getProducts(); i++) {
                    theWriter.write(MolfileV3Creator.writeCTAB(rxn.getProduct(i), scalingFactor));
                }
                theWriter.write("M  V30 END PRODUCT"+NL);
            }
            theWriter.write("M  END"+NL);
            rxnbuffer = theWriter.getBuffer();
            theWriter.close();
        } catch (Exception e) {
            System.err.println("Error in RXNFileCreator: " + e);
        }
    }

    private double getScalingFactor(Reaction rxn) {
        double avbl = 0;
        int bondCount = 0;
        for (int m=0; m<rxn.getMolecules(); m++) {
            avbl += rxn.getMolecule(m).getAverageBondLength() * rxn.getMolecule(m).getAllBonds();
            bondCount += rxn.getMolecule(m).getAllBonds();
        }

        if (bondCount != 0)
            return bondCount / avbl;

        return 1.0;
    }

    public String getRXNfile()
    {
        return rxnbuffer != null ? rxnbuffer.toString() : null;
    }


    public void writeRXNfile(Writer theWriter) throws IOException {
        if (rxnbuffer == null)
            throw new IOException("NULL RXNFileBuffer!");
        theWriter.write(rxnbuffer.toString());
    }
}
