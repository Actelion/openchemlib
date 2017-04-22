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
import com.actelion.research.share.gui.editor.Model;

/**
 * Project:
 * User: rufenec
 * Date: 5/16/13
 * Time: 3:31 PM
 */
public abstract class CopyAction extends CommandAction
{
    public CopyAction(Model model)
    {
        super(model);
    }

    @Override
    public void onCommand()
    {
//        ClipboardHandler handler = new ClipboardHandler();
//        handler.copyMolecule(model.getHighlightedMolecule());
        copy();
    }


    private void copy()
    {
        int mMode = model.getMode();

        boolean isReaction = ((mMode & Model.MODE_REACTION) != 0);
        boolean selectionFound = false;
        boolean isBothSideSelection = false;
        boolean isOnProductSide = false;

        StereoMolecule mMol = model.getMolecule();
        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            if (mMol.isSelectedAtom(atom)) {
                if (!selectionFound) {
                    selectionFound = true;
                    if (!isReaction) {
                        break;
                    }
                    isOnProductSide = model.isOnProductSide(mMol.getAtomX(atom), mMol.getAtomY(atom));
                } else {
                    if (isOnProductSide != model.isOnProductSide(mMol.getAtomX(atom), mMol.getAtomY(atom))) {
                        isBothSideSelection = true;
                        break;
                    }
                }
            }
        }

        if (isReaction) {
            if (isBothSideSelection) {
                copyReaction(true);
            } else if (selectionFound) {
                copyMolecule(true);
            } else {
                copyReaction(false);
            }
        } else {
            copyMolecule(selectionFound);
        }
    }

    public abstract boolean copyReaction(boolean selectionOnly);
    public abstract boolean copyMolecule(boolean selectionOnly);

//    private boolean copyReaction(boolean selectionOnly)
//    {
////       ClipboardHandler mClipboardHandler = new ClipboardHandler();
////        Reaction rx = selectionOnly ? getSelectedReaction() : model.getReaction();
////        if (rx != null && mClipboardHandler != null) {
////            return mClipboardHandler.copyReaction(rx);
////        }
////
//        return false;
//    }

//    private Reaction getSelectedReaction()
//    {
//        Reaction rxn = new Reaction();
//        StereoMolecule[] mFragment = model.getFragments();
//        int mReactantCount = model.getReactantCount();
//        for (int i = 0; i < mFragment.length; i++) {
//            StereoMolecule selectedMol = getSelectedCopy(mFragment[i]);
//            if (selectedMol != null) {
//                if (i < mReactantCount) {
//                    rxn.addReactant(selectedMol);
//                } else {
//                    rxn.addProduct(selectedMol);
//                }
//            }
//        }
//        return rxn;
//    }

//    private boolean copyMolecule(boolean selectionOnly)
//    {
////        ClipboardHandler mClipboardHandler = new ClipboardHandler();
////        StereoMolecule mMol = model.getMolecule();
////        if (mMol.getAllAtoms() != 0 && mClipboardHandler != null) {
////            return mClipboardHandler.copyMolecule(selectionOnly ? getSelectedCopy(mMol) : mMol);
////        }
////
//        return false;
//    }

/*    private StereoMolecule getSelectedCopy(StereoMolecule sourceMol)
     {
         int atomCount = 0;
         for (int atom = 0; atom < sourceMol.getAllAtoms(); atom++) {
             if (sourceMol.isSelectedAtom(atom)) {
                 atomCount++;
             }
         }

         if (atomCount == 0) {
             return null;
         }

         int bondCount = 0;
         for (int bond = 0; bond < sourceMol.getAllBonds(); bond++) {
             if (sourceMol.isSelectedBond(bond)) {
                 bondCount++;
             }
         }

         boolean[] includeAtom = new boolean[sourceMol.getAllAtoms()];
         for (int atom = 0; atom < sourceMol.getAllAtoms(); atom++) {
             includeAtom[atom] = sourceMol.isSelectedAtom(atom);
         }

         StereoMolecule destMol = new StereoMolecule(atomCount, bondCount);
         sourceMol.copyMoleculeByAtoms(destMol, includeAtom, false, null);
         return destMol;
     }*/
}
