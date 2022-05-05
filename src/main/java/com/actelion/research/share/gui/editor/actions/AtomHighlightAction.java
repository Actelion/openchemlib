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
import com.actelion.research.chem.NamedSubstituents;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.share.gui.DialogResult;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.dialogs.IAtomQueryFeaturesDialog;
import com.actelion.research.share.gui.editor.geom.GeomFactory;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IKeyEvent;
import com.actelion.research.share.gui.editor.io.IMouseEvent;


/**
 * Project:
 * User: rufenec
 * Date: 2/1/13
 * Time: 4:13 PM
 */
public abstract class AtomHighlightAction extends DrawAction
{

    AtomHighlightAction(Model model)
    {
        super(model);
    }

    protected GenericPoint lastHightlightPoint = null;

    boolean trackHighLight(GenericPoint pt)
    {
        StereoMolecule mol = model.getMolecule();
        int currentAtom = model.getSelectedAtom();
        int atom = findAtom(mol, pt);
        // Update at least when current selected atom and new selected atom differ
        boolean ok = atom != -1;
        lastHightlightPoint = pt;
        String keyStrokes = model.getKeyStrokeBuffer().toString();
        if (currentAtom != -1 && atom != currentAtom) {
            if (keyStrokes.length() > 0) {
                StereoMolecule currentMol = model.getMolecule();
                expandAtomKeyStrokes(currentMol, currentAtom, keyStrokes);
                model.getKeyStrokeBuffer().setLength(0);
            }
        }
        if (mol != null) {
            setHighlightAtom(mol, atom);
        } else {
            if (model.getSelectedAtom() != -1) {
                model.getKeyStrokeBuffer().setLength(0);
                ok = true;
            }
        }
        atom = model.getSelectedAtom();
        ok = ok | atom != -1;
        return ok;
    }

    void setHighlightAtom(StereoMolecule mol, int atom)
    {
        model.setSelectedAtom(atom);
    }


    @Override
    public boolean onMouseDown(IMouseEvent evt)
    {
        return false;
    }


    @Override
    public boolean onMouseMove(IMouseEvent evt, boolean drag)
    {
        if (!drag) {
            GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
            boolean ok = trackHighLight(pt);
            return ok;
        } else {
            return false;
        }

    }

    @Override
    public boolean paint(IDrawContext _ctx)
    {
        boolean ok = false;
        int theAtom = model.getSelectedAtom();
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            if (theAtom != -1) {
                drawAtomHighlight(_ctx, mol, theAtom);
                ok = true;
            }
        }
        return ok;
    }

    @Override
    public boolean onKeyPressed(IKeyEvent evt)
    {
        GeomFactory factory = model.getGeomFactory();
        if (evt.getCode().equals(factory.getDeleteKey())) {
            int theAtom = model.getSelectedAtom();
            StereoMolecule mol = model.getMolecule();
            if (theAtom != -1) {
                mol.deleteAtom(theAtom);
                setHighlightAtom(mol, -1);
                return true;
            } else {
                boolean update = false;
                StereoMolecule m = model.getMolecule();
                if (m.deleteSelectedAtoms()) {
                    update = true;
                }
                return update;
            }
        } else if (handleCharacter(evt)) {
            return true;
        }
        return false;
    }

    private boolean handleCharacter(IKeyEvent evt)
    {
        int theAtom = model.getSelectedAtom();
        StereoMolecule mol = model.getMolecule();
        StringBuilder keyStrokeBuffer = model.getKeyStrokeBuffer();
        if (mol == null) {
            return false;
        }
        int newRadical;
        String code = evt.getText();
        if (theAtom != -1 /*&& !evt.isShiftDown()*/) {
            char c = code != null && code.length() > 0 ? code.charAt(0) : 0;
            boolean isFirst = (keyStrokeBuffer.length() == 0);
            if (isFirst) {
                switch (c) {
                    case '+':
                    case '-':
                        return mol.changeAtomCharge(theAtom, c == '+');

                    case '.':
                        newRadical = (mol.getAtomRadical(theAtom) == Molecule.cAtomRadicalStateD) ? 0 : Molecule.cAtomRadicalStateD;
                        mol.setAtomRadical(theAtom, newRadical);
                        return true;
                    case ':':
                        newRadical = (mol.getAtomRadical(theAtom) == Molecule.cAtomRadicalStateT) ? Molecule.cAtomRadicalStateS
                            : (mol.getAtomRadical(theAtom) == Molecule.cAtomRadicalStateS) ? 0 : Molecule.cAtomRadicalStateT;
                        mol.setAtomRadical(theAtom, newRadical);
                        return true;

                    case 'q':
                        return mol.isFragment() ? showAtomQFDialog(theAtom) : false;
                    case '?':
                        return mol.changeAtom(theAtom, 0, 0, -1, 0);

                    default:
                        break;
                }
            }
            return handleCharacter(mol, theAtom, evt);
        } else {
            return handleCharsNonSelected(code);
        }
    }

    private boolean handleCharsNonSelected(String code)
    {

        if ("h".equals(code)) {
            model.flip(true);
            return true;
        } else if ("v".equals(code)) {
            model.flip(false);
            return true;
        }
        return false;
    }

    private void expandAtomKeyStrokes(StereoMolecule mol, int highliteAtom, String keyStrokes)
    {
        int atomicNo = Molecule.getAtomicNoFromLabel(keyStrokes);
        if (atomicNo != 0) {
            if (mol.changeAtom(highliteAtom, atomicNo, 0, -1, 0)) {
                model.changed();
                return;
            }
        }

        StereoMolecule substituent = NamedSubstituents.getSubstituent(keyStrokes);
        if (substituent != null) {

            // Copy the the fragment containing the attachment point into a new molecule.
            // Then attach the substituent, create new atom coordinates for the substituent,
            // while retaining coordinates of the fragment.
            StereoMolecule fragment = new StereoMolecule();
            fragment.addFragment(mol, highliteAtom, null);
            double sourceAVBL = fragment.getAverageBondLength();
            int firstAtomInFragment = fragment.getAllAtoms();
            for (int atom = 0; atom < fragment.getAllAtoms(); atom++) {
                fragment.setAtomMarker(atom, true);
            }
            fragment.addSubstituent(substituent, 0);
            new CoordinateInventor(CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS).invent(fragment);

            double dx = mol.getAtomX(highliteAtom) - sourceAVBL * fragment.getAtomX(0);
            double dy = mol.getAtomY(highliteAtom) - sourceAVBL * fragment.getAtomY(0);

            // Attach the substituent to the complete molecule and take coodinates from the
            // previously created fragment-substituent species.
            int firstAtomInMol = mol.getAllAtoms();
            mol.addSubstituent(substituent, highliteAtom);
            int substituentAtoms = mol.getAllAtoms() - firstAtomInMol;
            for (int i = 0; i < substituentAtoms; i++) {
                mol.setAtomX(firstAtomInMol + i, sourceAVBL * fragment.getAtomX(firstAtomInFragment + i) + dx);
                mol.setAtomY(firstAtomInMol + i, sourceAVBL * fragment.getAtomY(firstAtomInFragment + i) + dy);
            }
            mol.setStereoBondsFromParity();
            model.changed();
        }
    }

    private boolean handleCharacter(StereoMolecule mol, int theAtom, IKeyEvent evt)
    {
        GeomFactory factory = model.getGeomFactory();
        StringBuilder keyStrokeBuffer = model.getKeyStrokeBuffer();
        boolean isFirst = (keyStrokeBuffer.length() == 0);
        String code = evt.getText();
        char c = code != null && code.length() > 0 ? code.charAt(0) : 0;
        if (evt.getCode().equals(factory.getDeleteKey())) {
            if (theAtom != -1) {
                mol.deleteAtom(theAtom);
                setHighlightAtom(mol, -1);
                return true;
            } else {
                if (mol.deleteSelectedAtoms()) {
                    return true;
                }
            }
        } else if (!isFirst && evt.getCode().equals(factory.getEscapeKey())) {
            keyStrokeBuffer.setLength(0);
            return true;
        } else if (!isFirst && evt.getCode().equals(factory.getBackSpaceKey())) {
            keyStrokeBuffer.setLength(keyStrokeBuffer.length() - 1);
            return true;
        } else if (evt.getCode().equals(factory.getEnterKey())) {
            expandAtomKeyStrokes(mol, theAtom, keyStrokeBuffer.toString());
            keyStrokeBuffer.setLength(0);
            return true;
        } else if (isFirst && c > 48 && c <= 57) {
            addChain(mol, theAtom, c-47);
            return true;
        } else if ((c > 48 && c <= 57) || (c >= 65 && c <= 90) || (c >= 97 && c <= 122) || (c >= 48 && c <= 57) || (c == '-')) {
            keyStrokeBuffer.append(c);
            return true;
        }
        return false;
    }

    private void addChain(StereoMolecule mol, int theAtom, int chainAtoms)
    {
        if (mol.getFreeValence(theAtom) > 0) {
            int atom1 = theAtom;
            int hydrogenCount = mol.getAllAtoms() - mol.getAtoms();
            for (int i = 1; i < chainAtoms; i++) {
                GenericPoint pt = suggestNewX2AndY2(atom1);
                int atom2 = mol.addAtom(pt.getX(), pt.getY());
                if (atom2 == -1) {
                    break;
                }
                mol.addBond(atom1, atom2);
                atom1 = atom2 - hydrogenCount;    // new atom was added behind all hydrogens and travels now to the front
                mol.ensureHelperArrays(Molecule.cHelperNeighbours);
            }
        }
    }


    private boolean showAtomQFDialog(int atom)
    {
        GeomFactory factory = model.getGeomFactory();
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            boolean showReactionHints = (model.getMode() & Model.MODE_REACTION) != 0;
            IAtomQueryFeaturesDialog dlg = factory.createAtomQueryFeatureDialog(mol, atom, showReactionHints);
            return dlg.doModalAt(lastHightlightPoint.getX(), lastHightlightPoint.getY()) == DialogResult.IDOK;
        }
        return false;
    }

    public int findAtom(StereoMolecule mol, GenericPoint pt)
    {
        int foundAtom = -1;
        double pickx = pt.getX();
        double picky = pt.getY();
        double avbl = mol.getAverageBondLength();
        double foundDistanceSquare = Float.MAX_VALUE;
        double maxDistanceSquare = avbl * avbl / (avbl / 3);
        int mAllAtoms = mol.getAllAtoms();
        for (int atom = 0; atom < mAllAtoms; atom++) {
            double x = mol.getAtomX(atom);
            double y = mol.getAtomY(atom);
            double distanceSquare = (pickx - x) * (pickx - x) + (picky - y) * (picky - y);
            if (distanceSquare < maxDistanceSquare && distanceSquare < foundDistanceSquare) {
                foundDistanceSquare = distanceSquare;
                foundAtom = atom;
            }
        }
        return foundAtom;
    }


}
