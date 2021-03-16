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

package com.actelion.research.gui.dnd;

import com.actelion.research.chem.*;
import com.actelion.research.chem.dnd.ChemistryFlavors;

import java.awt.datatransfer.*;
import java.io.IOException;

public class MoleculeTransferable implements Transferable,ClipboardOwner {
    protected StereoMolecule mMol;

    public MoleculeTransferable(StereoMolecule mol) {
            mMol = mol;
            }

    @Override
    public synchronized DataFlavor[] getTransferDataFlavors() {
            return ChemistryFlavors.MOLECULE_FLAVORS;
            }

    @Override
    public boolean isDataFlavorSupported( DataFlavor flavor ) {
        for (DataFlavor f:ChemistryFlavors.MOLECULE_FLAVORS)
            if (f.equals(flavor))
                return true;

        return false;
        }

    @Override
    public synchronized Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException,IOException {
        if (flavor.equals(ChemistryFlavors.DF_SERIALIZED_MOLECULE)) {
            return new StereoMolecule(mMol);
        } else if (flavor.equals(ChemistryFlavors.DF_MDLMOLFILEV3)) {
            return new MolfileV3Creator(mMol).getMolfile();
        } else if (flavor.equals(ChemistryFlavors.DF_MDLMOLFILE)) {
            return new MolfileCreator(mMol).getMolfile();
        } else if (flavor.equals(ChemistryFlavors.DF_SMILES)) {
            return new IsomericSmilesCreator(mMol).getSmiles();
        } else if (flavor.equals(ChemistryFlavors.DF_IDCODE) || flavor.equals(DataFlavor.stringFlavor)) {
            final Canonizer canonizer = new Canonizer(mMol);
            return String.format("%s %s", canonizer.getIDCode(), canonizer.getEncodedCoordinates(true));
        } else
            throw new UnsupportedFlavorException(flavor);
    }

    public String toString()
    {
        return "MoleculeTransferable";
    }

    public void lostOwnership(Clipboard clipboard, Transferable contents) {
    }
}

