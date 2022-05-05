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

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.geom.GeomFactory;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.io.IKeyEvent;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

/**
 * Project:
 * User: rufenec
 * Date: 5/22/13
 * Time: 4:00 PM
 */
public class AtomMapAction extends AtomHighlightAction
{

    private GenericPoint firstPoint = null;
    private GenericPoint lastPoint = null;
    private int secondAtom = -1;

    public AtomMapAction(Model model)
    {
        super(model);
    }

    public  void onActionEnter()
    {
        int mode = model.getDisplayMode();
        model.setDisplayMode(mode | AbstractDepictor.cDModeShowMapping| AbstractDepictor.cDModeSuppressCIPParity);
    }

    public void onActionLeave()
    {
        int mode = model.getDisplayMode();
        if ((mode & AbstractDepictor.cDModeShowMapping) != 0) {
            mode &= ~(AbstractDepictor.cDModeShowMapping | AbstractDepictor.cDModeSuppressCIPParity);
            model.setDisplayMode(mode);
        }
    }

//    @Override
//    public boolean onMouseDown(ACTMouseEvent evt)
//    {
//        origin = new GenericPoint(evt.getX(),evt.getY());
//        return super.onMouseDown(evt);
//    }


    @Override
    public boolean onKeyPressed(IKeyEvent evt)
    {
        GeomFactory builder = model.getGeomFactory();
        if (evt.getCode().equals(builder.getDeleteKey())) {
            StereoMolecule mMol = model.getMolecule();
            boolean found = false;
            for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                if (mMol.getAtomMapNo(atom) != 0) {
                    mMol.setAtomMapNo(atom, 0, false);
                    found = true;
                }
            }
            return found;
        }
        return super.onKeyPressed(evt);
    }

    @Override
    public boolean onMouseMove(IMouseEvent evt, boolean drag)
    {
        firstPoint = lastPoint = null;
        if (model.isReaction()) {
            StereoMolecule mol = model.getMolecule();
            if (!drag) {
                GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
                secondAtom = -1;
                if(trackHighLight(pt)) {
                    int mCurrentHiliteAtom = model.getSelectedAtom();
                    if (mCurrentHiliteAtom != -1) {
                        int mapNo = mol.getAtomMapNo(mCurrentHiliteAtom);
                        if (mapNo != 0) {
                            for (int atom = 0; atom < mol.getAtoms(); atom++) {
                                if (atom != mCurrentHiliteAtom
                                    && mol.getAtomMapNo(atom) == mapNo) {
                                    secondAtom = atom;
                                    break;
                                }
                            }
                        }
                    }
                    return true;
                }
            } else {
                int atom = model.getSelectedAtom();
                if (mol != null && atom != -1) {
                    GenericPoint pt = new GenericPoint(evt.getX(), evt.getY());
                    firstPoint = new GenericPoint(mol.getAtomX(atom), mol.getAtomY(atom));
                    lastPoint = pt;
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    boolean trackHighLight(GenericPoint pt) {
        int lastAtom = model.getSelectedAtom();
        boolean ok = super.trackHighLight(pt);
        int theAtom = model.getSelectedAtom();
        return ok || lastAtom != theAtom;
    }

    @Override
    public boolean onMouseUp(IMouseEvent ev)
    {
        int mode = model.getDisplayMode();
        if ((mode & AbstractDepictor.cDModeShowMapping) == 0) {
            mode |= AbstractDepictor.cDModeShowMapping;
            //System.out.println("Display Mode " + mode);
            model.setDisplayMode(mode);
        }
        int atom = model.getSelectedAtom();
        if (atom != -1) {
            model.mapReaction(atom,firstPoint,lastPoint);
            //assistedMap(atom);
        }
        model.setSelectedAtom(-1);
        firstPoint = lastPoint = null;
        return true;
    }

//    private void assistedMap(int atom)
//    {
//        StereoMolecule mol = model.getSelectedMolecule();
//        int freeMapNo = model.getNextMapNo();
//        if (mol != null) {
//            StereoMolecule source = model.getFragmentAt(firstPoint, false);
//            StereoMolecule target = model.getFragmentAt(lastPoint, false);
//            if (target != null && target != source) {
//                int dest = mol.findAtom((int) lastPoint.getX(), (int) lastPoint.getY());
//                if (dest != -1) {
//                    mol.setAtomMapNo(atom, freeMapNo, false);
//                    mol.setAtomMapNo(dest, freeMapNo, false);
//                }
//                model.tryAutoMapReaction();
//            }
//        }
//    }
//

    @Override
    public boolean paint(IDrawContext ctx)
    {
        boolean ok = false;
        ctx.save();
        if (model.isReaction()) {
            StereoMolecule mol = model.getMolecule();
            if (firstPoint != null && lastPoint != null) {
                StereoMolecule source = model.getFragmentAt(firstPoint, false);
                StereoMolecule target = model.getFragmentAt(lastPoint, false);
                if (target != null && target != source) {
                    int theAtom = mol.findAtom((float)lastPoint.getX(),(float)lastPoint.getY());
                    if (theAtom != -1)
                        drawAtomHighlight(ctx, mol, theAtom);
                }
//                ctx.setStroke(IColor.RED);
                GeomFactory builder = model.getGeomFactory();
                ctx.setStroke(builder.getMapToolColor());
                ctx.drawLine(firstPoint.getX(), firstPoint.getY(), lastPoint.getX(), lastPoint.getY());
            } else if( secondAtom != -1) {
                drawAtomHighlight(ctx, mol, secondAtom);
            }
        }
        ok = super.paint(ctx);
        ctx.restore();
        return ok;
    }

}
