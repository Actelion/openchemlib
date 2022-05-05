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
import com.actelion.research.share.gui.DrawConfig;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.geom.IColor;
import com.actelion.research.share.gui.editor.geom.ICursor;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IKeyEvent;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Basic class which handles all actions which interact with the drawing surface
 */
public abstract class DrawAction implements Action
{
    public static final double HIGHLIGHT_ATOM_RADIUS = 5;
    public static final int MAX_CONNATOMS = 8;
    public static final String UNKNOWN = "<unknown>";
    public static final int KEYSTROKEFONTSIZE = 24;

    protected Model model;

    public DrawAction(Model m)
    {
        model = m;
    }


    @Override
    public void onCommand()
    {
        // NoOp
    }

    @Override
    public final boolean isCommand()
    {
        return false;
    }

    @Override
    public boolean onKeyPressed(IKeyEvent evt)
    {
        return false;
    }

    @Override
    public boolean onKeyReleased(IKeyEvent evt)
    {
        return false;
    }

    @Override
    public int getCursor()
    {
        return ICursor.DEFAULT;
    }

    @Override
    public boolean onDoubleClick(IMouseEvent ev)
    {
        return false;
    }

    public void onActionLeave()
    {

    }

    public void onActionEnter()
    {

    }

    /**
     * returns the atom at the current point
     *
     * @param pt
     * @return atom no or -1 if no atom there
     */
    int getAtomAt(GenericPoint pt)
    {
        StereoMolecule mol = model.getMoleculeAt(pt, false);
        return getAtomAt(mol, pt);
    }



    int getBondAt(GenericPoint pt)
    {
        StereoMolecule mol = model.getMoleculeAt(pt, true);
        return getBondAt(mol, pt);
    }

    int getBondAt(StereoMolecule mol, GenericPoint pt)
    {
        if (mol != null) {
            return mol.findBond((float) pt.getX(), (float) pt.getY());
        }
        return -1;
    }



    int getAtomAt(StereoMolecule mol, GenericPoint pt)
    {
        if (mol != null) {
            return mol.findAtom((float) pt.getX(), (float) pt.getY());
        }
        return -1;
    }


    protected void drawBondHighlight(IDrawContext ctx, StereoMolecule mol, int theBond)
    {
        int width = (int) (0.32f * mol.getAverageBondLength());
        if (width < 5)
            width = (int)HIGHLIGHT_ATOM_RADIUS;

        double x1 = mol.getAtomX(mol.getBondAtom(0, theBond));
        double y1 = mol.getAtomY(mol.getBondAtom(0, theBond));
        double x2 = mol.getAtomX(mol.getBondAtom(1, theBond));
        double y2 = mol.getAtomY(mol.getBondAtom(1, theBond));

        DrawConfig cfg = model.getGeomFactory().getDrawConfig();

        ctx.save();
        ctx.setLineWidth(width);
        ctx.setStroke(cfg.getHighLightColor());
//        ctx.setStroke(IColor.BLUE);
        ctx.drawLine(x1, y1, x2, y2);
        ctx.restore();
    }


    protected void highlightAtom(StereoMolecule mol, int atom)
    {
//        model.setSelectedMolecule(mol);
        model.setSelectedAtom(atom);
    }


    protected void drawAtomHighlight(IDrawContext ctx, StereoMolecule mol, int theAtom)
    {
        drawAtomHighlightElement(ctx, mol, theAtom);
        if (model.getKeyStrokeBuffer().length() > 0) {
            drawAtomKeyStrokes(ctx, mol, theAtom);
        }
    }

    private void drawAtomHighlightElement(IDrawContext ctx, StereoMolecule mol, int theAtom)
    {
        int radius = (int) (0.32f * mol.getAverageBondLength());
        if (radius < 5)
            radius = (int)HIGHLIGHT_ATOM_RADIUS;

        DrawConfig cfg = model.getGeomFactory().getDrawConfig();
        GenericPoint highlightPoint = new GenericPoint(mol.getAtomX(theAtom), mol.getAtomY(theAtom));
        ctx.save();
        ctx.setFill(cfg.getHighLightColor());
        ctx.fillElipse(
                highlightPoint.getX() - radius, highlightPoint.getY() - radius ,
                2 * radius, 2 * radius
        );
        ctx.restore();
    }


    protected void drawAtomKeyStrokes(IDrawContext ctx, StereoMolecule mol, int theAtom)
    {

        String s = model.getKeyStrokeBuffer().toString();
        int validity = model.getAtomKeyStrokeValidity(s);
        GenericPoint highlightPoint = new GenericPoint(mol.getAtomX(theAtom), mol.getAtomY(theAtom));
        DrawConfig cfg = model.getGeomFactory().getDrawConfig();

        ctx.save();
        ctx.setFill((
                validity == Model.KEY_IS_ATOM_LABEL) ? cfg.getForegroundColor()
                : (validity == Model.KEY_IS_SUBSTITUENT) ? IColor.BLUE
                : (validity == Model.KEY_IS_VALID_START) ? IColor.GRAY : IColor.RED);

        if (validity == Model.KEY_IS_INVALID) {
            s = s + UNKNOWN;
        }
        ctx.setFont(ctx.getFont(), KEYSTROKEFONTSIZE, false);
        ctx.fillText(s, highlightPoint.getX(), highlightPoint.getY());
        ctx.restore();
    }

    protected GenericPoint suggestNewX2AndY2(int atom)
    {
        StereoMolecule mol = model.getMolecule();// .getSelectedMolecule();
        mol.ensureHelperArrays(Molecule.cHelperNeighbours);

        double newAngle = Math.PI * 2 / 3;
        if (atom != -1) {
            double[] angle = new double[DrawAction.MAX_CONNATOMS + 1];
            //angle[0] = Math.PI * 3 / 4;
            for (int i = 0; i < mol.getAllConnAtoms(atom); i++) {
                angle[i] = mol.getBondAngle(atom, mol.getConnAtom(atom, i));
            }

            if (mol.getAllConnAtoms(atom) == 0) {
                newAngle = Math.PI * 2 /3 ;
            } else if (mol.getAllConnAtoms(atom) == 1) {
                if (angle[0] < -Math.PI * 5 / 6) {
                    newAngle = Math.PI / 3;
                } else if (angle[0] < -Math.PI / 2) {
                    newAngle = Math.PI * 2 / 3;
                } else if (angle[0] < -Math.PI / 6) {
                    newAngle = Math.PI / 3;
                } else if (angle[0] < 0.0) {
                    newAngle = Math.PI * 2 / 3;
                } else if (angle[0] < Math.PI / 6) {
                    newAngle = -Math.PI * 2 / 3;
                } else if (angle[0] < Math.PI / 2) {
                    newAngle = -Math.PI / 3;
                } else if (angle[0] < Math.PI * 5 / 6) {
                    newAngle = -Math.PI * 2 / 3;
                } else {
                    newAngle = -Math.PI / 3;
                }
            } else {
                for (int i = mol.getAllConnAtoms(atom) - 1; i > 0; i--) {  // bubble sort
                    for (int j = 0; j < i; j++) {
                        if (angle[j] > angle[j + 1]) {
                            double temp = angle[j];
                            angle[j] = angle[j + 1];
                            angle[j + 1] = temp;
                        }
                    }
                }
                angle[mol.getAllConnAtoms(atom)] = angle[0] + Math.PI * 2;

                int largestNo = 0;
                double largestDiff = 0.0;
                for (int i = 0; i < mol.getAllConnAtoms(atom); i++) {
                    double angleDiff = angle[i + 1] - angle[i];
                    if (largestDiff < angleDiff) {
                        largestDiff = angleDiff;
                        largestNo = i;
                    }
                }
                newAngle = (angle[largestNo] + angle[largestNo + 1]) / 2;
            }
        }
        double avbl = mol.getAverageBondLength();
        return new GenericPoint(
                mol.getAtomX(atom) + avbl * Math.sin(newAngle),
                mol.getAtomY(atom) + avbl * Math.cos(newAngle)
        );
    }

}
