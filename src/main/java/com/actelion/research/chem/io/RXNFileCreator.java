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

import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;

public class RXNFileCreator {
	public static final String RXN_CODE_TAG = "OCL_RXN_V1.0:";
	private static final String NL = System.lineSeparator();

    private StringBuffer rxnbuffer = null;

    public RXNFileCreator(Reaction rxn)
    {
        this(rxn, null);
    }

    public RXNFileCreator(Reaction r, String programName) {
        Reaction rxn = new Reaction(r);
        try {
            StringWriter theWriter = new StringWriter();
            theWriter.write("$RXN"+NL);
            theWriter.write(programName != null ? programName : "");
			theWriter.write(NL+NL);
			theWriter.write(RXN_CODE_TAG+ ReactionEncoder.encode(r, true,
                    ReactionEncoder.INCLUDE_MAPPING
                        | ReactionEncoder.INCLUDE_COORDS
                        | ReactionEncoder.INCLUDE_CATALYSTS
                        | ReactionEncoder.RETAIN_REACTANT_AND_PRODUCT_ORDER));
			theWriter.write(NL);
            theWriter.write("  "+rxn.getReactants()+"  "+rxn.getProducts() + NL);

            double scale = getScalingFactor(rxn);

            for (int i=0; i<rxn.getMolecules(); i++) {
                theWriter.write("$MOL"+NL);
                new MolfileCreator(rxn.getMolecule(i),true, scale, null).writeMolfile(theWriter);
            }
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
            return (double)bondCount / avbl;

        return 1.0;
    }

    public String getRXNfile() {
        return rxnbuffer != null ? rxnbuffer.toString() : null;
    }


    public void writeRXNfile(Writer theWriter) throws IOException {
        if (rxnbuffer == null)
            throw new IOException("NULL RXNFileBuffer!");
        theWriter.write(rxnbuffer.toString());
    }
}
