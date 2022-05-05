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
import com.actelion.research.share.gui.editor.geom.ICursor;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Project:
 * User: rufenec
 * Date: 3/26/13
 * Time: 3:42 PM
 */
public class NewChainAction extends BondHighlightAction {

    private int sourceAtom = -1;
    private int numChainAtoms = 0;
    private double[] mChainAtomX = null;
    private double[] mChainAtomY = null;
    private int[] mChainAtom = null;

    public NewChainAction(Model model) {
        super(model);
    }

    public boolean onMouseDown(IMouseEvent evt) {
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        StereoMolecule mol = model.getMolecule();
        boolean update = false;
        origin = pt;
        sourceAtom = findAtom(mol, pt);
        if (sourceAtom != -1) {
            //if (mol != null)
            {
                if (mol.getAllConnAtoms(sourceAtom) == Model.MAX_CONNATOMS) {
                    return false;
                }
                origin = new GenericPoint(mol.getAtomX(sourceAtom), mol.getAtomY(sourceAtom));
                update = true;
                numChainAtoms = 0;
                mChainAtomX = null;
                mChainAtomY = null;
                mChainAtom = null;
            }
        } else {
            origin = new GenericPoint(evt.getX(), evt.getY());
            update = true;
            numChainAtoms = 0;
            mChainAtomX = null;
            mChainAtomY = null;
            mChainAtom = null;
        }

        return update;
    }

    public boolean onMouseUp(IMouseEvent evt) {
        boolean ok = false;
        model.pushUndo();
        GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
        StereoMolecule mol = model.getMolecule();
        if (numChainAtoms == 0) {
            mol = model.getMoleculeAt(pt, false);
            if (mol != null) {
                int atom = model.getSelectedAtom();
                if (atom != -1) {
                    addSingleBondAtAtom(mol, atom);
                }
            }

        } else if (numChainAtoms > 0) {
            if (sourceAtom == -1) {
                sourceAtom = mol.addAtom((float) origin.getX(), (float) origin.getY());
            }

            if (mChainAtom[0] == -1) {
                mChainAtom[0] = mol.addAtom((float) mChainAtomX[0], (float) mChainAtomY[0]);
            }

            if (mChainAtom[0] != -1) {
                mol.addBond(sourceAtom, mChainAtom[0]);
            }
            if(model.isReaction())
                model.needsLayout(true);
        }

        if (numChainAtoms > 1) {
            for (int i = 1; i < numChainAtoms; i++) {
                if (mChainAtom[i] == -1) {
                    mChainAtom[i] = mol.addAtom((float) mChainAtomX[i], (float) mChainAtomY[i]);
                }
                if (mChainAtom[i] != -1) {
                    mol.addBond(mChainAtom[i - 1], mChainAtom[i]);
                }
            }
            if(model.isReaction())
                model.needsLayout(true);
        }
        highlightAtom(mol, -1);
        ok = true;
        dragging = false;
        return ok;

    }

    private void addSingleBondAtAtom(StereoMolecule mol, int atom) {
        GenericPoint p = suggestNewX2AndY2(atom);
        int targetAtom = mol.findAtom((float) p.getX(), (float) p.getY());
        if (targetAtom != -1) {
            mol.addOrChangeBond(atom, targetAtom, mol.suggestBondType(atom, targetAtom));
        } else {
            targetAtom = mol.addAtom((float) p.getX(), (float) p.getY(), 0.0f);
            mol.addBond(atom, targetAtom, mol.suggestBondType(atom, targetAtom));
            mol.ensureHelperArrays(Molecule.cHelperNeighbours);
        }
    }


    @Override
    protected boolean onDrag(GenericPoint pt) {
        StereoMolecule mol = model.getMolecule();
        boolean repaintNeeded = false;
        if (mol != null) {
            double lastX, lastY;
            if (numChainAtoms > 0) {
                lastX = mChainAtomX[numChainAtoms - 1];
                lastY = mChainAtomY[numChainAtoms - 1];
            } else {
                lastX = 0.0;
                lastY = 0.0;
            }
            double avbl = mol.getAverageBondLength();
            double s0 = avbl;//.floor();
            double s1 = (0.866 * avbl);//.floor();
            double s2 = (0.5 * avbl);//.floor();
            double dx = pt.getX() - origin.getX();
            double dy = pt.getY() - origin.getY();
            double a = 1.0;// sqrt(avbl/2*avbl/2);
            double b = 1.0;//sqrt(avbl/2*avbl/2);


            if (Math.abs(dy) > Math.abs(dx)) {
                numChainAtoms = (int) (2 * Math.abs(dy) / (s0 + s2));
                if ((int) Math.abs(dy) % (int) (s0 + s2) > s0) {
                    numChainAtoms++;
                }
                mChainAtomX = new double[numChainAtoms];
                mChainAtomY = new double[numChainAtoms];
                if (pt.getX() < origin.getX()) {
                    b = -b;
                }
                if (pt.getY() < origin.getY()) {
                    a = -a;
                }
                if (numChainAtoms > 0) {
                    mChainAtomX[0] = origin.getX() + s1 * b;
                    mChainAtomY[0] = origin.getY() + s2 * a;
                    for (int i = 1; i < numChainAtoms; i++) {
                        if ((i & 1) == 0) {
                            mChainAtomX[i] = mChainAtomX[i - 1] + s0 * b;
                            mChainAtomY[i] = mChainAtomY[i - 1] + s2 * a;
                        } else {
                            mChainAtomX[i] = mChainAtomX[i - 1];
                            mChainAtomY[i] = mChainAtomY[i - 1] + s0 * a;
                        }
                    }
                }
            } else {
                numChainAtoms = (int) (Math.abs(dx) / s1);
                mChainAtomX = new double[numChainAtoms];
                mChainAtomY = new double[numChainAtoms];
                if (pt.getX() < origin.getX()) {
                    s1 = -s1;
                }
                if (pt.getY() < origin.getY()) {
                    s2 = -s2;
                }
                for (int i = 0; i < numChainAtoms; i++) {
                    mChainAtomX[i] = origin.getX() + (i + 1) * s1;
                    mChainAtomY[i] = origin.getY();
                    if ((i & 1) == 0) {
                        mChainAtomY[i] += s2;
                    }
                }
            }
            if (numChainAtoms > 0) {
                mChainAtom = new int[numChainAtoms];
                for (int i = 0; i < numChainAtoms; i++) {
                    mChainAtom[i] = mol.findAtom((float) mChainAtomX[i], (float) mChainAtomY[i]);
                    if (mChainAtom[i] != -1) {
                        mChainAtomX[i] = mol.getAtomX(mChainAtom[i]);
                        mChainAtomY[i] = mol.getAtomY(mChainAtom[i]);
                    }
                }
                if (mChainAtomX[numChainAtoms - 1] != lastX
                        || mChainAtomY[numChainAtoms - 1] != lastY) {
                    repaintNeeded = true;
                }
            } else if (lastX != 0 || lastY != 0) {
                repaintNeeded = true;
            }
        }
        return repaintNeeded;

    }

    @Override
    public boolean paint(IDrawContext ctx) {
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            if (!dragging) {
                super.paint(ctx);
            } else {
                int theAtom = model.getSelectedAtom();
                drawChain(ctx, theAtom != -1 ? new GenericPoint(mol.getAtomX(theAtom), mol.getAtomY(theAtom)) : origin);
            }
        }
        return false;
    }

    private void drawChain(IDrawContext ctx, GenericPoint pt) {
        if (numChainAtoms > 0) {
            drawLine(ctx, pt.getX(), pt.getY(), mChainAtomX[0], mChainAtomY[0]);
        }
        if (numChainAtoms > 1) {
            for (int i = 1; i < numChainAtoms; i++) {
                drawLine(ctx, mChainAtomX[i - 1], mChainAtomY[i - 1], mChainAtomX[i], mChainAtomY[i]);
            }
        }
    }


    void drawLine(IDrawContext _ctx, double x1, double y1, double x2, double y2) {
        _ctx.drawLine(x1, y1, x2, y2);
    }

    @Override
    public int getCursor() {
        return ICursor.TOOL_CHAINCURSOR;
    }

}
