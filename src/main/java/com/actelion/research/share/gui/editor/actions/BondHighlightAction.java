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
import com.actelion.research.share.gui.DialogResult;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.dialogs.IBondQueryFeaturesDialog;
import com.actelion.research.share.gui.editor.geom.GeomFactory;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IKeyEvent;
import com.actelion.research.share.gui.editor.io.IMouseEvent;


/**
 * Project:
 * User: rufenec
 * Date: 1/28/13
 * Time: 1:07 PM
 */
public abstract class BondHighlightAction extends AtomHighlightAction
{
    GenericPoint origin = null;
    GenericPoint last = null;
    boolean dragging = false;

    public BondHighlightAction(Model model)
    {
        super(model);
    }


    @Override
    boolean trackHighLight(GenericPoint pt)
    {
        int lastAtom = model.getSelectedAtom();
        int lastBond = model.getSelectedBond();

        lastHightlightPoint = pt;
        StereoMolecule mol = model.getMoleculeAt(pt, true);
        if (super.trackHighLight(pt)) {
            setHighlightBond(null, -1);
            return true;
        }

        int bond = getBondAt(mol, pt);
        if (bond >= 0) {
            setHighlightBond(mol, bond);
            setHighlightAtom(mol, -1);
            return true;
        }
        boolean update = lastAtom != -1 || lastBond != bond ;
        setHighlightBond(null, -1);
        setHighlightAtom(null, -1);
        return update;
    }

    void setHighlightBond(StereoMolecule mol, int bond)
    {
        model.setSelectedBond(bond);
    }


    public boolean onMouseDown(IMouseEvent evt)
    {
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        origin = pt;
        return false;
    }


    public boolean onMouseMove(IMouseEvent evt, boolean drag)
    {
        dragging = drag;
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        if (!drag) {
            return trackHighLight(pt);
        } else {
            return onDrag(pt);
        }
    }

    protected boolean onDrag(GenericPoint pt)
    {
        double dx = Math.abs(pt.getX() - origin.getX());
        double dy = Math.abs(pt.getY() - origin.getY());
        if (dx > 5 || dy > 5) {
            trackHighLight(pt);
            last = pt;
        } else {
            last = null;
        }
        return true;
    }


    public boolean paint(IDrawContext _ctx)
    {
        StereoMolecule mol = model.getMolecule();
        boolean ok = false;
        if (mol != null) {
            int theAtom = model.getSelectedAtom();
            int theBond = model.getSelectedBond();
            if (theBond != -1) {
                drawBondHighlight(_ctx, mol, theBond);
                ok = true;
            } else if (theAtom != -1) {
                return super.paint(_ctx);
            }
        }
        return ok;

    }

    @Override
    public boolean onKeyPressed(IKeyEvent evt)
    {
        int theBond = model.getSelectedBond();
        GeomFactory factory = model.getGeomFactory();
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            if (evt.getCode().equals(factory.getDeleteKey())) {
                if (theBond != -1) {
                    mol.deleteBondAndSurrounding(theBond);
                    setHighlightBond(mol, -1);
                    return true;
                }
            } else {
                if (handleCharacter(evt.getText())) {
                    return true;
                }
            }
        }
        return super.onKeyPressed(evt);
    }

    private boolean handleCharacter(String code)
    {
        StereoMolecule mol = model.getMolecule();
        int theBond = model.getSelectedBond();

        if (mol != null && theBond != -1) {
            if (code != null && code.length() > 0) {
                char c = code.charAt(0);
                switch (c) {
                    case 'q':
                        return mol.isFragment() ? showBondQFDialog(theBond) : false;

                    case '4':
                    case '5':
                    case '6':
                    case '7':
                        return mol.addRingToBond(theBond, c - '0', false, Molecule.getDefaultAverageBondLength());

                    case 'b':
                        return mol.addRingToBond(theBond, 6, true, Molecule.getDefaultAverageBondLength());

                    case '1':
                        return mol.changeBond(theBond, Molecule.cBondTypeSingle);

                    case '2':
                        return mol.changeBond(theBond, Molecule.cBondTypeDouble);

                    case '3':
                        return mol.changeBond(theBond, Molecule.cBondTypeTriple);

                    case 'u':
                        return mol.changeBond(theBond, Molecule.cBondTypeUp);

                    case 'd':
                        return mol.changeBond(theBond, Molecule.cBondTypeDown);

                    case 'c':
                        return mol.changeBond(theBond, Molecule.cBondTypeCross);

                    case 'm':
                        return mol.changeBond(theBond, Molecule.cBondTypeMetalLigand);
                }
            }
        }
        return false;
    }

    private boolean showBondQFDialog(int bond)
    {
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            GeomFactory factory = model.getGeomFactory();
            IBondQueryFeaturesDialog dlg = factory.createBondFeaturesDialog(/*new BondQueryFeaturesDialog(*/mol, bond);
            return dlg.doModalAt(lastHightlightPoint.getX(),lastHightlightPoint.getY()) == DialogResult.IDOK;
        }
        return false;
    }

}

