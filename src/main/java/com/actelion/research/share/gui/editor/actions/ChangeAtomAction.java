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

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Project:
 * User: rufenec
 * Date: 2/1/13
 * Time: 4:13 PM
 */
public class ChangeAtomAction extends AtomHighlightAction {

    int theAtomNo = 6;

    public ChangeAtomAction(Model model, int atomNo) {
        super(model);
        theAtomNo = atomNo;
    }


    @Override
    public boolean onMouseUp(IMouseEvent evt) {
        model.pushUndo();
        int theAtom = model.getSelectedAtom();
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        StereoMolecule mol = model.getMolecule();
        if (theAtom != -1) {
            mol.setAtomicNo(theAtom, theAtomNo);
        } else {
            int atom = mol.addAtom((float) pt.getX(), (float) pt.getY());
            mol.setAtomicNo(atom, theAtomNo);
        }
        setHighlightAtom(mol, -1);
        return true;
    }

    @Override
    boolean trackHighLight(GenericPoint pt) {
        int lastAtom = model.getSelectedAtom();
        boolean ok = super.trackHighLight(pt);
        int theAtom = model.getSelectedAtom();
        return ok || lastAtom != theAtom;
    }
}
