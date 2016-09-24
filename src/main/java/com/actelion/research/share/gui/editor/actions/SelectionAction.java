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
import com.actelion.research.share.gui.DialogResult;
import com.actelion.research.share.gui.editor.Model;
import com.actelion.research.share.gui.editor.chem.IDrawingObject;
import com.actelion.research.share.gui.editor.dialogs.IAtomQueryFeaturesDialog;
import com.actelion.research.share.gui.editor.dialogs.IBondQueryFeaturesDialog;
import com.actelion.research.share.gui.editor.geom.GeomFactory;
import com.actelion.research.share.gui.editor.geom.ICursor;
import com.actelion.research.share.gui.editor.geom.IDrawContext;
import com.actelion.research.share.gui.editor.geom.IPolygon;
import com.actelion.research.share.gui.editor.io.IKeyEvent;
import com.actelion.research.share.gui.editor.io.IMouseEvent;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.List;


/**
 * Project:
 * User: rufenec
 * Date: 1/24/13
 * Time: 5:57 PM
 */
public class SelectionAction extends BondHighlightAction//DrawAction
{

    private volatile IPolygon polygon ;
    private GeomFactory factory ;
    int atom = -1;
    int bond = -1;
    boolean shift = false, duplicate = false,rectangular = false;

    public SelectionAction(Model model) {
        super(model);
        factory = model.getGeomFactory();
        polygon = factory.createPolygon();

    }

    @Override
    public boolean onKeyPressed(IKeyEvent evt) {
        shift = evt.isShiftDown();
        rectangular = evt.isAltDown();
        return super.onKeyPressed(evt);
    }

    @Override
    public boolean onKeyReleased(IKeyEvent evt) {
        shift = false;
        rectangular = false;
        duplicate = false;

        return super.onKeyReleased(evt);
    }

    @Override
    public boolean onMouseDown(IMouseEvent evt) {
        java.awt.geom.Point2D pt = new Point2D.Double(evt.getX(), evt.getY());
        StereoMolecule mol = model.getMoleculeAt(pt, true);

        polygon = factory.createPolygon();
        polygon.add(pt);

        duplicate = false;

        last = origin = new Point2D.Double(pt.getX(), pt.getY());
        atom = getAtomAt(mol, origin);
        bond = getBondAt(mol, origin);
        if (atom != -1) {
            if (!mol.isSelectedAtom(atom)) {
                if (!shift)
                    deselectAllAtoms();
                mol.setAtomSelection(atom, true);
            }
        } else if (bond != -1) {
            if (!mol.isSelectedBond(bond)) {
                if (!shift)
                    deselectAllAtoms();
                int a1 = mol.getBondAtom(0, bond);
                int a2 = mol.getBondAtom(1, bond);
                mol.setAtomSelection(a1, true);
                mol.setAtomSelection(a2, true);
            }
        }
        return false;
    }

    @Override
    public boolean onMouseUp(IMouseEvent ev) {
        polygon = factory.createPolygon();
        atom = bond = -1;
        origin = last = null;
        duplicate = false;
        return true;
    }

    @Override
    public boolean onMouseMove(IMouseEvent evt, boolean drag) {
        boolean ok = false;
        java.awt.geom.Point2D pt = new Point2D.Double(evt.getX(), evt.getY());
        if (drag) {
            double dx = last.getX() - pt.getX();
            double dy = last.getY() - pt.getY();
            if (!shift || duplicate) {
                ok = moveAtomsAndBonds(dx, dy, model.getSelectedDrawingObject() != null);
                if (ok) {
                    moveSelectedDrawItems(dx, dy);
                } else if (model.getSelectedDrawingObject() != null) {
                    model.getSelectedDrawingObject().move((float) -dx, (float) -dy);
                    ok = true;
                }
                if (!ok) {
                    ok = selectItems(pt);
                }
            } else if (shift) {
                duplicate = true;
                duplicateSelected();
            }
        } else {
            ok = trackHighLight(pt);
        }

        last = pt;
        return ok;
    }

    @Override
    boolean trackHighLight(java.awt.geom.Point2D pt) {
        boolean selected = false;
        IDrawingObject lastSelected = model.getSelectedDrawingObject();
        java.util.List<IDrawingObject> drawables = model.getDrawingObjects();
        model.setSelectedDrawingObject(null);

        for (IDrawingObject d : drawables) {
            if (d.getBoundingRect().contains(pt.getX(), pt.getY())) {
                model.setSelectedDrawingObject(d);
                selected = true;
                break;
            }
        }
        boolean ok = selected || lastSelected != null || super.trackHighLight(pt);
        return ok;

    }



    @Override
    public boolean onDoubleClick(IMouseEvent evt) {
//        StereoMolecule mol = model.getSelectedMolecule();
        java.awt.geom.Point2D pt = new Point2D.Double(evt.getX(), evt.getY());
        StereoMolecule mol = model.getMoleculeAt(pt, true);
        if (mol != null) {
            int atom = mol.findAtom((float) pt.getX(), (float) pt.getY());
            int bond = mol.findBond((float) pt.getX(), (float) pt.getY());
            boolean mShiftIsDown = evt.isShiftDown();
            int mMode = model.getMode();

            if (mol.isFragment()) {
                if (atom != -1) {
                    return showAtomQFDialog(atom);
                } else if (bond != -1) {
                    return showBondQFDialog(bond);
                }
//             else if (mCurrentHiliteObject != null) {
//                if (!mShiftIsDown) {
//                    for (int i = 0; i < mol.getAllAtoms(); i++)
//                        mol.setAtomSelection(i, false);
//                    for (int i = 0; i < mDrawingObjectList.size(); i++)
//                        ((AbstractDrawingObject) mDrawingObjectList.get(i)).setSelected(false);
//                }
//
//                mCurrentHiliteObject.setSelected(true);
//                update(UPDATE_REDRAW);
//            }
            } else {
                int rootAtom = -1;
                if (atom != -1) {
                    rootAtom = atom;
                } else if (bond != -1) {
                    rootAtom = mol.getBondAtom(0, bond);
                }

                if (rootAtom != -1 /*|| mCurrentHiliteObject != null*/) {
                    if (!mShiftIsDown) {
                        deselectAllAtoms();
                        for (int i = 0; i < mol.getAllAtoms(); i++) {
                            mol.setAtomSelection(i, true);
                        }
//                    if (mDrawingObjectList != null)
//                        for (AbstractDrawingObject drawingObject : mDrawingObjectList)
//                            drawingObject.setSelected(false);
                    }

//                if (rootAtom != -1) {
//                    model.setSelectedFragment(rootAtom);
//                } else {
////                    mCurrentHiliteObject.setSelected(true);
//                }
                    return true;
                }
            }
        } else {
            IDrawingObject drawingObject = model.getSelectedDrawingObject();
            if (drawingObject != null) {
                //System.out.printf("Doubleclick on seleted object\n");
                if (evt.isShiftDown()) {
                    deselectAllAtoms();
                    deselectAllDrawingObjects();
//                    model.setSelectedDrawingObject(null);
                }
                drawingObject.setSelected(true);
                model.setSelectedDrawingObject(drawingObject);
            }

        }
        return false;
    }


    @Override
    public int getCursor() {
        int ha = model.getSelectedAtom();
        int hb = model.getSelectedBond();
        StereoMolecule mol = model.getMolecule();

        if (shift && rectangular)
            return ICursor.TOOL_SELECTRECTPLUSCURSOR;

        if (shift)
            return ICursor.TOOL_LASSOPLUSCURSOR;

        if (rectangular)
            return ICursor.TOOL_SELECTRECTCURSOR;

        if (ha != -1 && mol.isSelectedAtom(ha))
            return shift ? ICursor.TOOL_HANDPLUSCURSOR : ICursor.TOOL_HANDCURSOR;

        if (hb != -1 && mol.isSelectedBond(hb))
            return shift ? ICursor.TOOL_HANDPLUSCURSOR : ICursor.TOOL_HANDCURSOR;

        if (ha != -1 || hb != -1)
            return ICursor.TOOL_POINTERCURSOR;
        return ICursor.TOOL_LASSOCURSOR;
    }


    @Override
    public boolean paint(IDrawContext ctx) {
        ctx.save();
        super.paint(ctx);
        ctx.setStroke(factory.getSelectionColor());
        if (rectangular && origin != null && last != null) {
            drawDashedRect(ctx);
            return true;
        } else if (polygon != null && polygon.size() > 1) {
            drawPolygon(ctx);
            return true;
        } else {
            ctx.setStroke(factory.getHighLightColor());
            ctx.setFill(factory.getHighLightColor());
            IDrawingObject obj =  model.getSelectedDrawingObject();
            if(obj != null) {
                obj.draw(ctx,null);
            }
        }
        ctx.restore();
        return false;
    }

    private void drawDashedRect(IDrawContext ctx) {
        java.awt.geom.Rectangle2D rc = makeRect(origin, last);
        if (rc.getWidth() > 5 && rc.getHeight() > 5) {
            drawDashedRect(ctx, rc.getMinX(), rc.getMinY(), rc.getWidth(), rc.getHeight(), new int[]{
                    5,
                    2
            });
        }
    }

    private void drawPolygon(IDrawContext ctx) {
        ctx.drawPolygon(polygon);
    }

    private void drawDashedLine(IDrawContext context,
                                double srcx, double srcy,
                                double targetx, double targety,
                                int[] dashPattern) {
        context.drawDashedLine(srcx, srcy, targetx, targety, dashPattern);
    }

    private void drawDashedRect(IDrawContext ctx, double x, double y, double w, double h, int[] pattern) {
        drawDashedLine(ctx, x, y, x + w, y, pattern);
        drawDashedLine(ctx, x, y, x, y + h, pattern);
        drawDashedLine(ctx, x, y + h, x + w, y + h, pattern);
        drawDashedLine(ctx, x + w, y + h, x + w, y, pattern);
    }

    private boolean moveAtomsAndBonds(double dx, double dy, boolean force) {
        boolean ok = false;
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            if (mol != null && atom != -1 || force) {
                translateAtoms(mol, dx, dy, true);
                ok = true;
            } else if (mol != null && bond != -1 || force) {
                translateBonds(mol, dx, dy, true);
                ok = true;
            }
        }
        return ok;
    }

    private boolean selectItems(java.awt.geom.Point2D pt) {
        boolean ok = false;
        if (rectangular) {
            selectRectanglarRegion(null);
            ok = true;
        } else {
            if (selectPolygonRegion(null, pt)) {
                ok = true;
            }
        }
        return ok;
    }

    private boolean moveSelectedDrawItems(double dx, double dy) {
        boolean ok = false;
        for (IDrawingObject selectedOne : model.getDrawingObjects()) {
            if (selectedOne.isSelected()) {
                selectedOne.move((float) -dx, (float) -dy);
                ok = true;
            }
        }
        return ok;
    }

    private void duplicateSelected() {
        StereoMolecule mol = model.getMolecule();

        int originalAtoms = mol.getAllAtoms();
        int originalBonds = mol.getAllBonds();
        int[] atomMap = new int[mol.getAllAtoms()];
        int esrGroupCountAND = mol.renumberESRGroups(Molecule.cESRTypeAnd);
        int esrGroupCountOR = mol.renumberESRGroups(Molecule.cESRTypeOr);
        for (int atom = 0; atom < originalAtoms; atom++) {
            if (mol.isSelectedAtom(atom)) {
                int newAtom = mol.getAllAtoms();
                atomMap[atom] = newAtom;
                mol.copyAtom(mol, atom, esrGroupCountAND, esrGroupCountOR);
            }
        }
        for (int bond = 0; bond < originalBonds; bond++) {
            if (mol.isSelectedBond(bond)) {
                mol.copyBond(mol, bond, esrGroupCountAND, esrGroupCountOR, atomMap, false);
            }
        }
        for (int atom = 0; atom < originalAtoms; atom++) {
            mol.setAtomSelection(atom, false);
        }
        for (int atom = originalAtoms; atom < mol.getAllAtoms(); atom++) {
            mol.setAtomMapNo(atom, 0, false);
        }
    }


    private boolean selectPolygonRegion(StereoMolecule m, java.awt.geom.Point2D pt) {
        if (polygon.size() > 1 && Math.abs(pt.getX() - polygon.get(polygon.size() - 1).getX()) < 10
                && Math.abs(pt.getY() - polygon.get(polygon.size() - 1).getY()) < 10) {
            return false;
        }
        if (origin == null) {
            throw new RuntimeException("NUll DOWN Point!");
        }
        polygon.remove(origin);
        polygon.add(pt);
        polygon.add(origin);

        deselectAllAtoms();
        deselectAllDrawingObjects();
        selectFromPolygonRegion();
        return true;
    }

    private void selectFromPolygonRegion() {
        StereoMolecule mol = model.getMolecule();
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            boolean isSelected = polygon.contains(mol.getAtomX(i), mol.getAtomY(i));
            mol.setAtomSelection(i, isSelected);
//            model.setSelectedMolecule(mol);
        }

        List<IDrawingObject> drawables = model.getDrawingObjects();
        for (IDrawingObject d : drawables) {
            Rectangle2D r = d.getBoundingRect();
            if (polygon.contains(r.getCenterX(), r.getCenterY()))
                d.setSelected(true);
        }
    }

    private void deselectAllAtoms() {
        StereoMolecule mol = model.getMolecule();
        deselectAtoms(mol);
    }

    private void selectRectanglarRegion(StereoMolecule mol) {
        java.awt.geom.Rectangle2D rc = makeRect(origin, last);
        boolean selected = false;
        if (mol != null) {
            selectAtomsInRectangle(mol, rc);
            selected = true;
        } else {
            StereoMolecule m = model.getMolecule();
            deselectAtoms(m);
            java.awt.geom.Rectangle2D bounds = factory.getBoundingRect(m);
            if (bounds != null && bounds.intersects(rc.getX(), rc.getY(), rc.getWidth(), rc.getHeight())) {
                selectRectanglarRegion(m);
                //break;
            }
        }
        if (!selected) {
            selectDrawingObjectsInRectangle(rc);
        }
    }

    private void deselectAtoms(StereoMolecule mol) {
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            mol.setAtomSelection(i, false);
        }
    }

    private void selectDrawingObjectsInRectangle(java.awt.geom.Rectangle2D rc) {
        for (IDrawingObject dw : model.getDrawingObjects()) {
            dw.setSelected(false);
            java.awt.geom.Rectangle2D r = dw.getBoundingRect();
            if (rc.contains(r.getCenterX(), r.getCenterY()))
                dw.setSelected(true);
        }
    }

    private void selectAtomsInRectangle(StereoMolecule mol, java.awt.geom.Rectangle2D rc) {
        for (int i = 0; i < mol.getAllAtoms(); i++) {
            boolean isSelected = rc.contains(mol.getAtomX(i), mol.getAtomY(i));
            mol.setAtomSelection(i, isSelected);
        }
//        model.setSelectedMolecule(mol);
    }

    private void translateBonds(StereoMolecule mol, double dx, double dy, boolean selected) {
        int a1 = mol.getBondAtom(0, bond);
        int a2 = mol.getBondAtom(1, bond);
        if (selected) {
            translateAtoms(mol, dx, dy, true);
        } else {
            translateAtom(mol, a1, dx, dy);
            translateAtom(mol, a2, dx, dy);
        }
    }

    private void translateAtom(StereoMolecule mol, int a, double dx, double dy) {
        mol.setAtomX(a, mol.getAtomX(a) - dx);
        mol.setAtomY(a, mol.getAtomY(a) - dy);
    }

    private void translateAtoms(StereoMolecule mol, double dx, double dy, boolean allSelected) {
        if (allSelected) {
            for (int i = 0; i < mol.getAllAtoms(); i++) {
                if (mol.isSelectedAtom(i)) {
                    translateAtom(mol, i, dx, dy);
                }
            }
        } else {
            translateAtom(mol, atom, dx, dy);
        }
    }

    private java.awt.geom.Rectangle2D makeRect(java.awt.geom.Point2D origin, java.awt.geom.Point2D pt) {
        double x = Math.min(origin.getX(), pt.getX());
        double y = Math.min(origin.getY(), pt.getY());
        double w = Math.abs(origin.getX() - pt.getX());
        double h = Math.abs(origin.getY() - pt.getY());
        return new java.awt.geom.Rectangle2D.Double(x, y, w, h);
    }


    private void deselectAllDrawingObjects() {
        for (IDrawingObject d : model.getDrawingObjects()) {
            d.setSelected(false);
        }
    }


    private boolean showAtomQFDialog(int atom) {
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            IAtomQueryFeaturesDialog dlg = factory.createAtomQueryFeatureDialog(/* new AtomQueryFeaturesDialog*/mol, atom);
            return dlg.doModalAt(lastHightlightPoint.getX(), lastHightlightPoint.getY()) == DialogResult.IDOK;
        }
        return false;
    }

    private boolean showBondQFDialog(int bond) {
        StereoMolecule mol = model.getMolecule();
        if (mol != null) {
            IBondQueryFeaturesDialog dlg = factory.createBondFeaturesDialog( /*new BondQueryFeaturesDialog(*/mol, bond);
            return dlg.doModalAt(lastHightlightPoint.getX(), lastHightlightPoint.getY()) == DialogResult.IDOK;
        }
        return false;
    }

}

