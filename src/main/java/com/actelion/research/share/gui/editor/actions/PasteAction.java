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

import com.actelion.research.share.gui.editor.Model;

import java.awt.*;


/**
 * Project:
 * User: rufenec
 * Date: 5/16/13
 * Time: 3:39 PM
 */
public abstract class PasteAction extends CommandAction
{
    java.awt.Dimension bounds;
//    ClipboardHandler mClipboardHandler =new ClipboardHandler();

    public PasteAction(Model model,java.awt.Dimension bounds)
    {
        super(model);
        this.bounds = bounds;
    }

    @Override
    public void onCommand()
    {
//        ClipboardHandler handler =new ClipboardHandler();
//        StereoMolecule mol = handler.pasteMolecule();
//        if (mol != null) {
////                model.addMolecule(mol, new Rectangle2D(0, 0,bounds.getWidth(), bounds.getHeight()));
//        }
        paste();
    }

    public Dimension getBounds()
    {
        return bounds;
    }

    private void paste()
    {
        if ((model.getMode() & Model.MODE_REACTION) != 0 && pasteReaction())
            return;

        pasteMolecule();
    }

    public abstract boolean pasteMolecule();
    public abstract boolean pasteReaction();


/*
    private boolean pasteReaction()
    {
        boolean ret = false;
//        if (mClipboardHandler != null) {
//            Reaction rxn = mClipboardHandler.pasteReaction();
//            if (rxn != null) {
//                StereoMolecule mMol = model.getMolecule();
//                for (int i = 0; i < rxn.getMolecules(); i++) {
//                    rxn.getMolecule(i).setFragment(mMol.isFragment());
//                }
//                model.setReaction(rxn);
//                ret = true;
//            }
//        }
        return ret;
    }
*/

/*
    private boolean pasteMolecule()
    {
        boolean ret = false;
//        if (mClipboardHandler != null) {
//            StereoMolecule mol = mClipboardHandler.pasteMolecule();
//            if (mol != null && mol.getAllAtoms() != 0) {
//                StereoMolecule mMol = model.getMolecule();
//                if (mMol.getAllAtoms() == 0) {
//                    boolean isFragment = mMol.isFragment();
//                    mol.copyMolecule(mMol);
//                    mMol.setFragment(isFragment);
//                    model.notifyChange();
////                    moleculeChanged(true);
//                } else {
//                    int avbl = (int) mMol.getAverageBondLength();
//                    Depictor d = new Depictor(mol);
//                    System.err.println("Implement pasteMolecule()");
////                    d.updateCoords(this.getGraphics(), new Rectangle2D.Float(0, 0,
////                            this.getWidth(),
////                            this.getHeight()),
////                        AbstractDepictor.cModeInflateToMaxAVBL + avbl
////                    );
//                    int originalAtoms = mMol.getAllAtoms();
//                    mMol.addMolecule(mol);
//                    for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
//                        mMol.setAtomSelection(atom, atom >= originalAtoms);
//                    }
//                    model.notifyChange();
////                    moleculeChanged(true);
//                }
//            }
//            ret = true;
//        }
        return ret;
    }
*/

}


