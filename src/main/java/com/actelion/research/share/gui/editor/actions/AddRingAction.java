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

package com.actelion.research.share.gui.editor.actions;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Project:
 * User: rufenec
 * Date: 2/1/13
 * Time: 4:03 PM
 */
public class AddRingAction extends BondHighlightAction
{
    Model model;
    int ringSize = 0;
    boolean aromatic = false;

    public AddRingAction(Model m, int ringSize, boolean aromatic)
    {
        super(m);
        this.model = m;
        this.ringSize = ringSize;
        this.aromatic = aromatic;
    }

    public boolean onMouseDown(IMouseEvent evt)
    {
        return false;
    }

    public boolean onMouseUp(IMouseEvent evt)
    {
        model.pushUndo();
        StereoMolecule mol = model.getMolecule();
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        boolean ok = false;
        if (mol != null)
        {
            int atom = getAtomAt(mol,pt);
            int bond = getBondAt(mol,pt);
            if (atom != -1) {
                ok = mol.addRing((float) pt.getX(), (float) pt.getY(), ringSize, aromatic, Molecule.getDefaultAverageBondLength());
                model.setSelectedBond(-1);
            } else if (bond != -1) {
                ok = mol.addRing((float) pt.getX(), (float) pt.getY(), ringSize, aromatic, Molecule.getDefaultAverageBondLength());
                model.setSelectedAtom(-1);
            } else {
                ok = mol.addRing((float) pt.getX(), (float) pt.getY(), ringSize, aromatic, Molecule.getDefaultAverageBondLength());
                if (model.isReaction())
                    model.needsLayout(true);
            }
        } else {
            mol = new StereoMolecule();
            ok= mol.addRing((float) pt.getX(), (float) pt.getY(), ringSize, aromatic, Molecule.getDefaultAverageBondLength());
            model.setValue(mol, true);
        }
        if (ok)
            mol.ensureHelperArrays(Molecule.cHelperNeighbours);
        return ok;
    }

    @Override
    protected boolean onDrag(GenericPoint pt)
    {
        return true;
    }
}
