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
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Project:
 * User: rufenec
 * Date: 1/28/13
 * Time: 1:07 PM
 */
public abstract class BondBaseAction extends BondHighlightAction
{

    protected BondBaseAction(Model model)
    {
        super(model);
    }

    public abstract int getBondType();
    public void onAddBond(int srcAtom, int targetAtom)
    {
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            int bondType = getBondType();
            if (bondType == Molecule.cBondTypeSingle)
                bondType = mol.suggestBondType(srcAtom, targetAtom);
            mol.addBond(srcAtom, targetAtom, bondType);
            mol.ensureHelperArrays(Molecule.cHelperNeighbours);
        } else {
        }
    }

    public void onChangeBond(int bond)
    {
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            mol.changeBond(bond, Molecule.cBondTypeIncreaseOrder);
            mol.ensureHelperArrays(Molecule.cHelperNeighbours);
        }
    }

    public boolean onMouseUp(IMouseEvent evt)
    {

        boolean ok = true;
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        model.pushUndo();
        int sourceAtom = getAtomAt(origin);
        int selectedAtom = model.getSelectedAtom();
        model.setSelectedBond(-1);
        StereoMolecule mol = model.getMoleculeAt(origin, true);
        model.setSelectedAtom(sourceAtom);
        if (!dragging) {
            if (mol != null && sourceAtom != -1) {
                if (mol.getAllConnAtoms(sourceAtom) != Model.MAX_CONNATOMS) {
                    GenericPoint targetPoint = suggestNewX2AndY2(sourceAtom);
                    int stopAtom = mol.findAtom((float) targetPoint.getX(), (float) targetPoint.getY());
                    if (stopAtom != -1) {
                        int bondType = getBondType();
                        if (bondType == Molecule.cBondTypeSingle)
                            bondType = mol.suggestBondType(sourceAtom, stopAtom);
                        mol.addOrChangeBond(sourceAtom, stopAtom, bondType);
                    } else {
                        int targetAtom = mol.addAtom((float) targetPoint.getX(), (float) targetPoint.getY(), 0.0f);
                        onAddBond(sourceAtom, targetAtom);
                    }
                    ok = true;
                }
            } else if (mol != null) {
                int bond = getBondAt(mol, pt);
                if (bond != -1) {
                    onChangeBond(bond);
                } else {
                    sourceAtom = mol.addAtom((float) pt.getX(), (float) pt.getY());
                    GenericPoint targetPoint = suggestNewX2AndY2(sourceAtom);
                    int targetAtom = mol.addAtom((float) targetPoint.getX(), (float) targetPoint.getY(), 0.0f);
                    onAddBond(sourceAtom, targetAtom);
                }
                ok = true;
            } else  {
                mol = model.getMolecule();
                sourceAtom = mol.addAtom((float) evt.getX(), (float) evt.getY());
                GenericPoint p = suggestNewX2AndY2(sourceAtom);
                int t = mol.addAtom((float) p.getX() /*+ mol.getAverageBondLength()*/, (float) p.getY());
                onAddBond(sourceAtom, t);
                // This creates a new Fragment, so make sure scheme gets layouted (if in RXN mode)
                if (model.isReaction())
                    model.needsLayout(true);
                ok = true;
            }
        } else { // dragging
            if (mol != null) {
                if (sourceAtom != -1) {
                    int targetAtom = selectedAtom;
                    if (targetAtom == -1) {
                        double dx = origin.getX() - pt.getX();
                        double dy = origin.getY() - pt.getY();
                        GenericPoint targetPoint = pt;
                        if (dx * dx + dy * dy < Model.MIN_BOND_LENGTH_SQUARE) {
                            targetPoint = suggestNewX2AndY2(sourceAtom);
                        }
                        targetAtom = mol.addAtom((float) targetPoint.getX(), (float) targetPoint.getY(), 0.0f);
                    }
                    StereoMolecule tm = model.getMoleculeAt(pt, true);
                    if (mol == tm) {
                        onAddBond(sourceAtom, targetAtom);
                    } else if (tm != null) {
                        mol.addMolecule(tm);
                        targetAtom = mol.findAtom((float) pt.getX(), (float) pt.getY());
                        model.deleteMolecule(tm);
                        onAddBond(sourceAtom, targetAtom);
                    }
                    ok = true;
                } else { // mol == null
                    int startAtom = mol.addAtom((float) origin.getX(), (float) origin.getY());
                    int endAtom = mol.addAtom((float) pt.getX(), (float) pt.getY());
                    onAddBond(startAtom, endAtom);
                    if (model.isReaction())
                        model.needsLayout(true);
                    ok = true;
                }
            } else {
                mol = model.getMolecule();
                int startAtom = mol.addAtom((float) origin.getX(), (float) origin.getY());
                int endAtom = mol.addAtom((float) pt.getX(), (float) pt.getY());
                onAddBond(startAtom, endAtom);
                ok = true;
            }
        }
        dragging = false;
        return ok;

    }


    public boolean paint(IDrawContext _ctx)
    {
        boolean ok = super.paint(_ctx);
        if (dragging) {
            drawBondLine(_ctx);
            ok = true;
        }

        return ok;

    }

    private void drawBondLine(IDrawContext ctx)
    {
        GenericPoint point = origin;
        if (point != null && last != null) {
            int atom = getAtomAt(point);
            StereoMolecule mol = model.getMoleculeAt(point, true);
            if (mol != null && atom != -1) {
                point = new GenericPoint(mol.getAtomX(atom), mol.getAtomY(atom));
            }
            ctx.save();
            ctx.drawLine(point.getX(), point.getY(), last.getX(), last.getY());
            ctx.restore();
        }
    }


}
