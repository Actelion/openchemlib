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

package com.actelion.research.share.gui.editor;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.reaction.IReactionMapper;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.gui.generic.GenericPoint;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.share.gui.Arrow;
import com.actelion.research.share.gui.ChemistryGeometryHelper;
import com.actelion.research.share.gui.editor.chem.AbstractExtendedDepictor;
import com.actelion.research.share.gui.editor.chem.IDrawingObject;
import com.actelion.research.share.gui.editor.geom.GeomFactory;
import com.actelion.research.share.gui.editor.listeners.IChangeListener;
import com.actelion.research.share.gui.editor.listeners.IValidationListener;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Project:
 * User: rufenec
 * Date: 1/24/13
 * Time: 5:02 PM
 */
public abstract class Model
{
    public interface AtomHighlightCallback
    {
        void onHighlight(int atom, boolean selected);
    }

    public interface BondHighlightCallback
    {
        void onHighlight(int atom, boolean selected);
    }

    //    protected static GeomFactory builder = GeomFactory.getGeomFactory() ;
    protected GeomFactory geomFactory;
    public static final int KEY_IS_ATOM_LABEL = 1;
    public static final int KEY_IS_SUBSTITUENT = 2;
    public static final int KEY_IS_VALID_START = 3;
    public static final int KEY_IS_INVALID = 4;

    //    private StereoMolecule selectedMolecule;
    public static final int MODE_MULTIPLE_FRAGMENTS = 1;
    public static final int MODE_MARKUSH_STRUCTURE = 2;
    public static final int MODE_REACTION = 4;
    public static final int MODE_DRAWING_OBJECTS = 8;

    public static final int MAX_CONNATOMS = 8;
    public static final int MIN_BOND_LENGTH_SQUARE = 100;

    private static final float FRAGMENT_MAX_CLICK_DISTANCE = 24.0f;
    private static final float FRAGMENT_GROUPING_DISTANCE = 1.2f;    // in average bond lengths
    private static final float FRAGMENT_CLEANUP_DISTANCE = 1.5f;    // in average bond lengths
    private static final float DEFAULT_ARROW_LENGTH = 0.08f;        // relative to panel width

    public static final int FAKE_ATOM_NO = 100;

/*
    protected static final int UPDATE_NONE = 0;
    protected static final int UPDATE_REDRAW = 1;
    // redraw molecules and drawing objects with their current coordinates
    protected static final int UPDATE_CHECK_VIEW = 2;
    // redraw with on-the-fly coordinate transformation only if current coords do not fit within view area
    // (new coords are generated for one draw() only; the original coords are not changed)
    protected static final int UPDATE_CHECK_COORDS = 3;
    // redraw with in-place coordinate transformation only if current coords do not fit within view area
    // (the original atom and object coords are replaced by the new ones)
    protected static final int UPDATE_SCALE_COORDS = 4;
    // redraw with in-place coordinate transformation; molecules and objects are scaled to fill
    // the view unless the maximum average bond length reaches the optimum.
    protected static final int UPDATE_SCALE_COORDS_USE_FRAGMENTS = 5;
    // as UPDATE_SCALE_COORDS but uses fragments from mFragment rather than creating a
    // fresh mFragment list from mMol. Used for setting a reaction or fragment list from outside.
    protected static final int UPDATE_INVENT_COORDS = 6;
    // redraw with in-place coordinate transformation; first all molecules' coordinates
    // are generated from scratch, then molecules and objects are scaled to fill
    // the view unless the maximum average bond length reaches the optimum.
*/

    public static final int MAX_UNDO_SIZE = 5;
    //    private List<DrawItems> _undoList = new ArrayList<DrawItems>();
    private List<StereoMolecule> _undoList = new ArrayList<StereoMolecule>();
    private int selectedESRType = 0;

    private int selectedAtom = -1;
    private int selectedBond = -1;
    private int displayMode = 0;
//    private int reactantIndex = 0;
    private int mReactantCount;
    private int[] mFragmentNo;
    private boolean mAtomColorSupported;
    private GenericPoint arrowPos = new GenericPoint(0,0);

    private List<IValidationListener> validationListeners = new ArrayList<IValidationListener>();
    private List<IChangeListener> changeListeners = new ArrayList<IChangeListener>();
    private boolean needslayout = false;
    private int mMode = 0;

    private Dimension displaySize = new Dimension(0, 0);

    private StereoMolecule mMol = new StereoMolecule();        // molecule being modified directly by the drawing editor
    private StringBuilder mAtomKeyStrokeBuffer = new StringBuilder();
    private List<IDrawingObject> mDrawingObjectList;//, mUndoDrawingObjectList;
    private StereoMolecule[] mFragment;    // in case of MODE_MULTIPLE_FRAGMENTS contains valid stereo fragments
    private IDrawingObject selectedDrawingObject;

    private IReactionMapper mMapper;
    private ImageProvider imageProvider;


    public Model(GeomFactory factory, int mode)
    {
        this.geomFactory = factory;
        mDrawingObjectList = new ArrayList<IDrawingObject>();
        mMode = mode;
        if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
            mMode |= (MODE_MULTIPLE_FRAGMENTS);
        }

        if ((mMode & (MODE_DRAWING_OBJECTS | MODE_REACTION)) != 0) {

        }
        if (isReactionMode()) {
            arrowPos = new GenericPoint(0,0);
            Arrow arrow = new Arrow(factory.getDrawConfig(), 0, 0, 0, 0);
            mDrawingObjectList.add(arrow);
        }
    }

    public GeomFactory getGeomFactory()
    {
        return geomFactory;
    }

/*
    public void resizeReaction(Dimension os, Dimension ns)
    {
        Reaction rxn = getReaction(false);
        double offsetx = ns.getWidth() - os.getWidth();
        double offsety = ns.getHeight() - os.getHeight();
        double scale;

        if (Math.abs(offsetx) > Math.abs(offsety)) {
            scale = ns.getHeight() / os.getHeight();
        } else {
            scale = ns.getWidth() / os.getWidth();

        }
        ChemistryHelper.transformReaction(rxn, offsetx, offsety, scale);
        setValue(rxn);

    }
*/

    public void cleanReaction(boolean cleanAll)
    {
        Reaction reaction = getReaction();
//        System.out.printf("cleanreaction %s\n",reaction.getReactants());
        Dimension dim = getDisplaySize();
        double w = dim.getWidth();
        double h = dim.getHeight();
        double width = w / 5;

        if (w > 0 && h > 0) {
            IDrawingObject arrow = getDrawingObjects().get(0);// new Arrow(mx , my, dx,20);
            arrow.setRect((float) (0.5f * w), (float) (0.5f * h), (float) (0.5f * .16 * w), 20);
            arrowPos = new GenericPoint((0.5f * w),  (0.5f * h));
            mMode = MODE_MULTIPLE_FRAGMENTS;
            if (cleanAll)
                cleanupCoordinates(true, true);
            ChemistryGeometryHelper.scaleInto(reaction, 0, 0, dim.getWidth(), dim.getHeight(), width);
            setValue(reaction);
        }
    }


/*
    public void setSelectedMolecule(StereoMolecule selectedMolecule)
    {

        if (selectedMolecule == null) {
//            Exception e = new Exception("Passed NULL selected Molecule");
//            e.printStackTrace();

        }
//        this.selectedMolecule = selectedMolecule;
    }
*/

//    public StereoMolecule getSelectedMolecule()
//    {
//
//        //return selectedMolecule;
//        return mMol;
//    }

    public StereoMolecule getSelectedCopy(StereoMolecule sourceMol)
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
    }


    public void scale(float dx, float dy)
    {
        mMol.scaleCoords(Math.min(dx, dy));
    }

    public Reaction getSelectedReaction()
    {
        Reaction rxn = new Reaction();
        for (int i = 0; i < mFragment.length; i++) {
            StereoMolecule selectedMol = getSelectedCopy(mFragment[i]);
            if (selectedMol != null) {
                if (i < mReactantCount) {
                    rxn.addReactant(selectedMol);
                } else {
                    rxn.addProduct(selectedMol);
                }
            }
        }
        return rxn;
    }


    private int findFragment(float x, float y)
    {
        int fragment = -1;
        double minDistance = Float.MAX_VALUE;
        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            double dx = x - mMol.getAtomX(atom);
            double dy = y - mMol.getAtomY(atom);
            double distance = Math.sqrt(dx * dx + dy * dy);
            if (distance < FRAGMENT_MAX_CLICK_DISTANCE
                && minDistance > distance) {
                minDistance = distance;
                fragment = mFragmentNo[atom];
            }
        }

        return fragment;
    }

    public StereoMolecule[] getFragments()
    {
        return mFragment;
    }

    public void setFragments(StereoMolecule[] fragment)
    {
        mMol.deleteMolecule();
        mFragment = fragment;
        for (int i = 0; i < fragment.length; i++) {
            mMol.addMolecule(mFragment[i]);
        }
        pushUndo();

        mFragmentNo = new int[mMol.getAllAtoms()];
        for (int atom = 0, f = 0; f < mFragment.length; f++) {
            for (int j = 0; j < mFragment[f].getAllAtoms(); j++) {
                mFragmentNo[atom++] = f;
            }
        }

        mMode = MODE_MULTIPLE_FRAGMENTS;
        notifyChange();
//          fireEvent(new DrawAreaEvent(this, DrawAreaEvent.TYPE_MOLECULE_CHANGED, false));

//          update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
    }

    public Reaction getReaction()
    {
        if ((mMode & MODE_REACTION) == 0) {
            return null;
        }
        Reaction rxn = new Reaction();
        syncFragments();
        for (int i = 0; i < mFragment.length; i++) {
            if (i < mReactantCount) {
                rxn.addReactant(mFragment[i]);
            } else {
                rxn.addProduct(mFragment[i]);
            }
        }
//        for (int i = 0; i < mFragment.length; i++) {
//            for (int j = 0; j < mFragment[i].getAllAtoms(); j++) {
//                System.out.printf("getReaction Fragment[%s] atom[%s]=%s\n",i,j,mFragment[i].getAtomicNo(j));
//            }
//        }
//
        return rxn;
    }


    public void setReaction(Reaction rxn)
    {

        Dimension displaySize = getDisplaySize();
//        System.out.printf("SetReaction %s\n", displaySize);

        mMol = new StereoMolecule();
        mFragment = new StereoMolecule[rxn.getMolecules()];
        mReactantCount = rxn.getReactants();
        boolean isFragment = false;
        for (int i = 0; i < rxn.getMolecules(); i++) {
            isFragment |= rxn.getMolecule(i).isFragment();
            StereoMolecule molecule = rxn.getMolecule(i);
            GenericRectangle boundingRect = ChemistryGeometryHelper.getBoundingRect(molecule);
//            if(i < mReactantCount) {
//                arrowPos = new GenericPoint(boundingRect.getX()+boundingRect.getWidth(),boundingRect.getY()+boundingRect.getHeight()/2);
//            }
            mFragment[i] = molecule;
            mMol.addMolecule(mFragment[i]);
        }

//        System.out.printf("SetReaction %s mols = %d\n", rxn.getReactants(),rxn.getMolecules());

//        ChemistryHelper.scaleInto(rxn,0,0,(double)displaySize.width,(double)displaySize.height);

        mMol.setFragment(isFragment);

        mFragmentNo = new int[mMol.getAllAtoms()];
        for (int atom = 0, f = 0; f < mFragment.length; f++) {
            for (int j = 0; j < mFragment[f].getAllAtoms(); j++) {
                mFragmentNo[atom++] = f;
            }
        }

        try {
            mMol.validate();
        } catch (Exception e) {
            System.out.println("WARNING:" + e);
//            e.printStackTrace();
        }
//        fireEvent(new DrawAreaEvent(this, DrawAreaEvent.TYPE_MOLECULE_CHANGED, false));
        mMode = MODE_MULTIPLE_FRAGMENTS | MODE_REACTION;

//        for (int i = 0; i < mFragment.length; i++) {
//            for (int j = 0; j < mFragment[i].getAllAtoms(); j++) {
//                System.out.printf("setReaction Fragment[%s] atom[%s]=%s\n",i,j,mFragment[i].getAtomicNo(j));
//            }
//        }

        notifyChange();
//        update(UPDATE_SCALE_COORDS_USE_FRAGMENTS);
    }

    public MarkushStructure getMarkushStructure()
    {
        if ((mMode & MODE_MARKUSH_STRUCTURE) == 0) {
            return null;
        }

        MarkushStructure markush = new MarkushStructure();
        for (int i = 0; i < mFragment.length; i++) {
            if (i < mReactantCount) {
                markush.addCore(mFragment[i]);
            } else {
                markush.addRGroup(mFragment[i]);
            }
        }
        return markush;
    }

    public void setMarkushStructure(MarkushStructure markush)
    {
        mMol.deleteMolecule();
        mFragment = new StereoMolecule[markush.getCoreCount() + markush.getRGroupCount()];
        mReactantCount = markush.getCoreCount();
        boolean isFragment = false;
        for (int i = 0; i < markush.getCoreCount() + markush.getRGroupCount(); i++) {
            mFragment[i] = (i < markush.getCoreCount()) ? markush.getCoreStructure(i)
                : markush.getRGroup(i - markush.getCoreCount());
            isFragment |= mFragment[i].isFragment();
            mMol.addMolecule(mFragment[i]);
        }
        mMol.setFragment(isFragment);
        pushUndo();

        mFragmentNo = new int[mMol.getAllAtoms()];
        for (int atom = 0, f = 0; f < mFragment.length; f++) {
            for (int j = 0; j < mFragment[f].getAllAtoms(); j++) {
                mFragmentNo[atom++] = f;
            }
        }

//        fireEvent(new DrawAreaEvent(this, DrawAreaEvent.TYPE_MOLECULE_CHANGED, false));

        mMode = MODE_MULTIPLE_FRAGMENTS | MODE_MARKUSH_STRUCTURE;
        notifyChange();
    }

    public void setDisplayMode(int dMode)
    {
        displayMode = dMode;
        notifyChange();
    }


    public int getMode()
    {
        return mMode;
    }


    public StereoMolecule getMolecule()
    {
        return mMol;
    }


    public boolean isAtomColorSupported()
    {
        return mAtomColorSupported;
    }

    public void setAtomColorSupported(boolean acs)
    {
        mAtomColorSupported = acs;
    }

    protected void cleanupCoordinates(boolean multifragment, boolean invent)
    {
        int selectedAtomCount = 0;
        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            if (mMol.isSelectedAtom(atom)) {
                selectedAtomCount++;
            }
        }
        boolean selectedOnly = (selectedAtomCount != 0 && selectedAtomCount != mMol.getAllAtoms());

        if (!multifragment) {
            AbstractDepictor depictor = createDepictor(getMolecule());
            cleanupMoleculeCoordinates(depictor, invent, selectedOnly);
        } else {
            AbstractExtendedDepictor depictor = createExtendedDepictor();
            cleanupMultiFragmentCoordinates(depictor, selectedOnly, invent);
        }

        if (selectedOnly)
            mMol.removeAtomMarkers();
    }


    private void cleanupMoleculeCoordinates(AbstractDepictor depictor, boolean invent, boolean selectedOnly)
    {
        //if (mUpdateMode == UPDATE_INVENT_COORDS)
        if (invent) {
            if (selectedOnly) {
                for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                    mMol.setAtomMarker(atom, !mMol.isSelectedAtom(atom));
                }
            }

            new CoordinateInventor(selectedOnly ? CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mMol);
//			mMol.setStereoBondsFromParity(); not needed anymore
        }

        DepictorTransformation dt = depictor.simpleValidateView(new GenericRectangle(0, 0, this.getWidth(), this.getHeight()), AbstractDepictor.cModeInflateToMaxAVBL);
        if (dt != null)
            dt.applyTo(mMol);

    }


    private float getHeight()
    {
        return (float) displaySize.getHeight();
    }

    private float getWidth()
    {
        return (float) displaySize.getWidth();
    }

/*
    private void cleanupMultiFragmentCoordinatesEx(AbstractExtendedDepictor depictor, boolean selectedOnly, boolean invent)
    {
        //if (selectedOnly && mUpdateMode == UPDATE_INVENT_COORDS)
        if (invent) {
            int[] fragmentAtom = new int[mFragment.length];
            for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                int fragment = mFragmentNo[atom];
                mFragment[fragment].setAtomMarker(fragmentAtom[fragment], !mMol.isSelectedAtom(atom));
                fragmentAtom[fragment]++;
            }
        }

        java.awt.geom.Rectangle2D.Double[] boundingRect = new java.awt.geom.Rectangle2D.Double[mFragment.length];
//		float fragmentWidth = 0.0f;
        for (int fragment = 0; fragment < mFragment.length; fragment++) {
            //if (mUpdateMode == UPDATE_INVENT_COORDS)
            if (invent) {
                new CoordinateInventor(selectedOnly ?
                    CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mFragment[fragment]);
                mFragment[fragment].setStereoBondsFromParity();
            }
            AbstractDepictor d = createDepictor(mFragment[fragment]);// new Depictor(mFragment[fragment]);
//            depictor.updateCoords(null, null, AbstractDepictor.cModeInflateToMaxAVBL);
            depictor.updateCoords(null, null, AbstractDepictor.cModeInflateToMaxAVBL);

//            boundingRect[fragment] = d.getBoundingRect();
//			fragmentWidth += boundingRect[fragment].width;
        }

        double spacing = FRAGMENT_CLEANUP_DISTANCE * AbstractDepictor.cOptAvBondLen;
        double avbl = mMol.getAverageBondLength();
//        float arrowWidth = ((mMode & MODE_REACTION) == 0) ?
//            0f
//            : (mUpdateMode == UPDATE_SCALE_COORDS_USE_FRAGMENTS) ?
//            DEFAULT_ARROW_LENGTH * getWidth()
//            : ((IArrow) mDrawingObjectList.get(0)).getLength() * AbstractDepictor.cOptAvBondLen / avbl;

        double rawX = 0.5 * spacing;
//        for (int fragment = 0; fragment <= mFragment.length; fragment++) {
//            if ((mMode & MODE_REACTION) != 0 && fragment == mReactantCount) {
//                ((IArrow) mDrawingObjectList.get(0)).setCoordinates(
//                    rawX - spacing / 2, getHeight() / 2, rawX - spacing / 2 + arrowWidth, getHeight() / 2);
//                rawX += arrowWidth;
//            }
//
//            if (fragment == mFragment.length) {
//                break;
//            }
//
//            float dx = rawX - boundingRect[fragment].x;
//            float dy = 0.5f * (getHeight() - boundingRect[fragment].height)
//                - boundingRect[fragment].y;
//            mFragment[fragment].translateCoords(dx, dy);
//
//            rawX += spacing + boundingRect[fragment].width;
//        }

//        depictor.updateCoords(null, new java.awt.geom.Rectangle2D.Double(0, 0, getWidth(), getHeight()),
//            AbstractDepictor.cModeInflateToMaxAVBL);

        int[] fragmentAtom = new int[mFragment.length];
        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            int fragment = mFragmentNo[atom];
            mMol.setAtomX(atom, mFragment[fragment].getAtomX(fragmentAtom[fragment]));
            mMol.setAtomY(atom, mFragment[fragment].getAtomY(fragmentAtom[fragment]));

            fragmentAtom[fragment]++;
        }

        mMol.setStereoBondsFromParity();
    }

*/
    private int maxUpdateMode()
    {
        return AbstractDepictor.cModeInflateToMaxAVBL /*+ HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen)*/;
    }


    private void cleanupMultiFragmentCoordinates(AbstractExtendedDepictor depictor, boolean selectedOnly, boolean invent)
    {
        if (selectedOnly && invent) {
            int[] fragmentAtom = new int[mFragment.length];
            for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                int fragment = mFragmentNo[atom];
                mFragment[fragment].setAtomMarker(fragmentAtom[fragment], !mMol.isSelectedAtom(atom));
                fragmentAtom[fragment]++;
            }
        }

        GenericRectangle[] boundingRect = new GenericRectangle[mFragment.length];
        //		float fragmentWidth = 0.0f;
        for (int fragment = 0; fragment < mFragment.length; fragment++) {
            if (invent) {
                new CoordinateInventor(selectedOnly ? CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS : 0).invent(mFragment[fragment]);
//              mFragment[fragment].setStereoBondsFromParity(); not needed anymore
            }
            AbstractDepictor d = createDepictor(mFragment[fragment]);
            d.updateCoords(null, null, AbstractDepictor.cModeInflateToMaxAVBL);
            boundingRect[fragment] = d.getBoundingRect();
        }

        double spacing = FRAGMENT_CLEANUP_DISTANCE * AbstractDepictor.cOptAvBondLen;
        double avbl = mMol.getAverageBondLength();
        double arrowWidth = isReaction() ? DEFAULT_ARROW_LENGTH * getWidth() : 0;
//            0f : true ?
//            : (mMode == UPDATE_SCALE_COORDS_USE_FRAGMENTS) ?
//            DEFAULT_ARROW_LENGTH * getWidth()
//            : mDrawingObjectList.get(0).getBoundingRect().getWidth() * AbstractDepictor.cOptAvBondLen / avbl;

        double rawX = 0.5 * spacing;
        for (int fragment = 0; fragment <= mFragment.length; fragment++) {
            if (isReaction() && fragment == mReactantCount) {
                mDrawingObjectList.get(0).setRect(
                        (float) (rawX - spacing / 20),
                        getHeight() / 2,
                        (float) (/*rawX - spacing / 2 + */arrowWidth),
                        getHeight() / 2);
                rawX += arrowWidth;
            }

            if (fragment == mFragment.length) {
                break;
            }

            double dx = rawX - boundingRect[fragment].x;
            double dy = 0.5 * (getHeight() - boundingRect[fragment].height)
                - boundingRect[fragment].y;
            mFragment[fragment].translateCoords(dx, dy);

            rawX += spacing + boundingRect[fragment].width;
        }

        depictor.updateCoords(null, new GenericRectangle(0, 0, getWidth(), getHeight()), maxUpdateMode());

        int[] fragmentAtom = new int[mFragment.length];
        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            int fragment = mFragmentNo[atom];
            mMol.setAtomX(atom, mFragment[fragment].getAtomX(fragmentAtom[fragment]));
            mMol.setAtomY(atom, mFragment[fragment].getAtomY(fragmentAtom[fragment]));
            fragmentAtom[fragment]++;
        }

        mMol.setStereoBondsFromParity();
    }

    public void analyzeFragmentMembership()
    {
        mMol.ensureHelperArrays(Molecule.cHelperParities);

        int[] fragmentNo = new int[mMol.getAllAtoms()];
        int fragments = mMol.getFragmentNumbers(fragmentNo, false, true);

        fragments = joinCloseFragments(fragmentNo, fragments);
        sortFragmentsByPosition(fragmentNo, fragments);
        mFragmentNo = fragmentNo;

        mFragment = mMol.getFragments(fragmentNo, fragments);
    }

    private void syncFragments()
    {
        mMol.ensureHelperArrays(Molecule.cHelperParities);
        int[] fragmentNo = new int[mMol.getAllAtoms()];
        int fragments = mMol.getFragmentNumbers(fragmentNo, false, true);
        mFragment = mMol.getFragments(fragmentNo, fragments);
        fragments = joinCloseFragments(fragmentNo, fragments);

        sortFragmentsByPosition(fragmentNo, fragments);
        mFragmentNo = fragmentNo;
        mFragment = mMol.getFragments(fragmentNo, fragments);

        for (StereoMolecule m : mFragment) {
            m.ensureHelperArrays(Molecule.cHelperParities);
        }
    }

    private int joinCloseFragments(int[] fragmentNo, int fragments)
    {
        if (fragments < 2) {
            return fragments;
        }

        boolean[][] mergeFragments = new boolean[fragments][];
        for (int i = 1; i < fragments; i++) {
            mergeFragments[i] = new boolean[i];
        }

        double avbl = mMol.getAverageBondLength();
        for (int atom1 = 1; atom1 < mMol.getAllAtoms(); atom1++) {
            for (int atom2 = 0; atom2 < atom1; atom2++) {
                double dx = mMol.getAtomX(atom2) - mMol.getAtomX(atom1);
                double dy = mMol.getAtomY(atom2) - mMol.getAtomY(atom1);
                double distance = Math.sqrt(dx * dx + dy * dy);
                if (distance < FRAGMENT_GROUPING_DISTANCE * avbl) {
                    int fragment1 = fragmentNo[atom1];
                    int fragment2 = fragmentNo[atom2];
                    if (fragment1 != fragment2) {
                        if (fragment1 > fragment2) {
                            mergeFragments[fragment1][fragment2] = true;
                        } else {
                            mergeFragments[fragment2][fragment1] = true;
                        }
                    }
                }
            }
        }

        int[] newFragmentIndex = new int[fragments];
        for (int fragment = 0; fragment < fragments; fragment++) {
            newFragmentIndex[fragment] = fragment;
        }

        int mergeCount = 0;
        for (int i = 1; i < fragments; i++) {
            for (int j = 0; j < i; j++) {
                if (mergeFragments[i][j]) {
                    int index1 = newFragmentIndex[i];
                    int index2 = newFragmentIndex[j];
                    if (index1 != index2) {
                        mergeCount++;
                        int minIndex = Math.min(index1, index2);
                        int maxIndex = Math.max(index1, index2);
                        for (int k = 0; k < fragments; k++) {
                            if (newFragmentIndex[k] == maxIndex) {
                                newFragmentIndex[k] = minIndex;
                            } else if (newFragmentIndex[k] > maxIndex) {
                                newFragmentIndex[k]--;
                            }
                        }
                    }
                }
            }
        }

        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            fragmentNo[atom] = newFragmentIndex[fragmentNo[atom]];
        }

        return fragments - mergeCount;
    }

    private void sortFragmentsByPosition(int[] fragmentNo, int fragments)
    {
        int[][] fragmentDescriptor = new int[fragments][((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) ? 2 : 1];
        for (int fragment = 0; fragment < fragments; fragment++) {
            fragmentDescriptor[fragment][0] = fragment;
        }
        Point[] fragmentCOG = calculateFragmentCenterOfGravity(fragmentNo, fragments);

        if (isReactionMode()) {
            mReactantCount = 0;
//            Arrow arrow = ((mMode & MODE_REACTION) != 0) ? (Arrow) mDrawingObjectList.get(0) : null;
//            System.out.printf("Arrow placement %s\n",arrow.getBoundingRect());
            for (int fragment = 0; fragment < fragments; fragment++) {
                fragmentDescriptor[fragment][1] =
                    isOnProductSide(fragmentCOG[fragment].x,fragmentCOG[fragment].y)
//                    arrow.isOnProductSide(fragmentCOG[fragment].x,fragmentCOG[fragment].y)
                        ? 1 : 0;
                if (fragmentDescriptor[fragment][1] == 0) {
                    mReactantCount++;
                }
            }
        } else if ((mMode & MODE_MARKUSH_STRUCTURE) != 0) {
            mReactantCount = fragments;
            for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                if (mMol.getAtomicNo(atom) == 0 && fragmentDescriptor[fragmentNo[atom]][1] == 0) {
                    fragmentDescriptor[fragmentNo[atom]][1] = 1;
                    mReactantCount--;
                }
            }
        }

        final Point[] cog = fragmentCOG;
        Arrays.sort(fragmentDescriptor, new Comparator<int[]>()
        {
            public int compare(int[] fragmentDescriptor1, int[] fragmentDescriptor2)
            {
                if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
                    if (fragmentDescriptor1[1] != fragmentDescriptor2[1]) {
                        return (fragmentDescriptor1[1] == 0) ? -1 : 1;
                    }
                }

                return cog[fragmentDescriptor1[0]].x < cog[fragmentDescriptor2[0]].x ? -1 : 1;

//                return (cog[fragmentDescriptor1[0]].x
//                    + cog[fragmentDescriptor1[0]].y
//                    < cog[fragmentDescriptor2[0]].x
//                    + cog[fragmentDescriptor2[0]].y) ? -1 : 1;
            }
        });

        int[] newFragmentIndex = new int[fragments];
        Point[] centerOfGravity = new Point[fragments];
        for (int fragment = 0; fragment < fragments; fragment++) {
            int oldIndex = fragmentDescriptor[fragment][0];
            newFragmentIndex[oldIndex] = fragment;
            centerOfGravity[fragment] = fragmentCOG[oldIndex];
        }

        fragmentCOG = centerOfGravity;
        for (int atom1 = 0; atom1 < mMol.getAllAtoms(); atom1++) {
            fragmentNo[atom1] = newFragmentIndex[fragmentNo[atom1]];
        }

//        System.out.printf("SortFragments Count %s\n",mReactantCount);

    }

    private boolean isReactionMode()
    {
        return (mMode & MODE_REACTION) != 0;
    }

    public boolean isOnProductSide(double x, double y)
    {
//        Arrow arrow = ((mMode & MODE_REACTION) != 0) ? (Arrow) mDrawingObjectList.get(0) : null;
//        if (arrow != null) {
//            System.out.printf("Arrow determines product side");
//            return arrow.isOnProductSide((float) x, (float) y);
//        }
        if (isReactionMode()) {
            return x > arrowPos.getX();
        }
        return x > getDisplaySize().getWidth() / 2;
    }

    private Point[] calculateFragmentCenterOfGravity(int[] fragmentNo, int fragments)
    {
        Point[] fragmentCOG = new Point[fragments];
        int[] fragmentAtoms = new int[fragments];
        for (int fragment = 0; fragment < fragments; fragment++) {
            fragmentCOG[fragment] = new Point(0, 0);
        }
        for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
            fragmentCOG[fragmentNo[atom]].x += mMol.getAtomX(atom);
            fragmentCOG[fragmentNo[atom]].y += mMol.getAtomY(atom);
            fragmentAtoms[fragmentNo[atom]]++;
        }
        for (int fragment = 0; fragment < fragments; fragment++) {
            fragmentCOG[fragment].x /= fragmentAtoms[fragment];
            fragmentCOG[fragment].y /= fragmentAtoms[fragment];
        }
        return fragmentCOG;
    }

    public void setMapper(IReactionMapper mapper)
    {
        mMapper = mapper;
    }

    public void mapReaction(int atom, GenericPoint left, GenericPoint right)
    {
        StereoMolecule mol = getMolecule();// getSelectedMolecule();
        if (mol != null && left != null && right != null) {
            int freeMapNo = getNextMapNo();
            StereoMolecule source = getFragmentAt(left, false);
            StereoMolecule target = getFragmentAt(right, false);
            boolean b1 = isOnProductSide(left.getX(), left.getY());
            boolean b2 = isOnProductSide(right.getX(), right.getY());

            if (target != null && target != source && b1 ^ b2) {
                int dest = mol.findAtom((int) right.getX(), (int) right.getY());
                if (dest != -1) {
                    mol.setAtomMapNo(atom, freeMapNo, false);
                    mol.setAtomMapNo(dest, freeMapNo, false);
                }
                if (mMapper != null)
                    tryAutoMapReaction();
            }
        }
    }


    /**
     * *********************************************************************************************************
     */
    public int getNextMapNo()
    {
        int freeMapNo = 1;
//        for (StereoMolecule mol : drawElement.molecules)
        {
            StereoMolecule mol = mMol;
            for (int i = 0; i < mol.getAtoms(); i++) {
                freeMapNo = Math.max(mol.getAtomMapNo(i) + 1, freeMapNo);
            }
        }
        return freeMapNo;
    }

    public void popUndo()
    {
        if (_undoList.size() > 0) {
            setValue(_undoList.get(_undoList.size() - 1), false);
            _undoList.remove(_undoList.size() - 1);
        }
    }

    public void pushUndo()
    {
//        _undoList.add(new StereoMolecule(getMol()));
        _undoList.add(new StereoMolecule(mMol));
        if (_undoList.size() > MAX_UNDO_SIZE) {
            _undoList.remove(0);
        }
//        mUndoDrawingObjectList = (mDrawingObjectList == null) ?
//            null : new ArrayList<IDrawingObject>(mDrawingObjectList);
    }


    public int getESRType()
    {
        return selectedESRType;
    }

    public void setESRType(int type)
    {
        selectedESRType = type;
        notifyChange();
    }

    public void addValidationListener(IValidationListener l)
    {
        if (!validationListeners.contains(l)) {
            validationListeners.add(l);
        }
    }

    public void removeValidationListener(IValidationListener l)
    {
        if (validationListeners.contains(l)) {
            validationListeners.remove(l);
        }
    }

    public void addChangeListener(IChangeListener l)
    {
        if (!changeListeners.contains(l)) {
            changeListeners.add(l);
        }
    }

    public void removeChangeListener(IChangeListener l)
    {
        if (changeListeners.contains(l)) {
            changeListeners.remove(l);
        }
    }


    public void setDisplaySize(Dimension displaySize)
    {
        if (isReactionMode()) {
            if (displaySize.getWidth() != 0 && displaySize.getHeight() != 0) {
                if (displaySize.getWidth() != this.displaySize.getWidth() || displaySize.getHeight() != this.displaySize.getHeight()) {
                    double dx = displaySize.getWidth()/this.displaySize.getWidth();
                    double dy = displaySize.getHeight()/this.displaySize.getHeight();
                    double scale = Math.min(dx,dy);
                    scale(scale);
                }
            }
//
        }
        this.displaySize = displaySize;
    }

    private void scale(double scale)
    {
        if (!Double.isInfinite(scale)) {
            if (scale != 1 && scale > 0) {
//                System.out.printf("Scale %f\n",scale);
                AbstractDepictor d = createDepictor(mMol);
                DepictorTransformation dt = d.simpleValidateView(new GenericRectangle(0, 0, this.getWidth(), this.getHeight()),
                    AbstractDepictor.cModeInflateToMaxAVBL + (int)mMol.getAverageBondLength());
                if (dt != null)
                    dt.applyTo(mMol);
            }
        }
    }

    public Dimension getDisplaySize()
    {
        return displaySize;
    }

    public void deleteMolecule(StereoMolecule tm)
    {
        System.err.println("DeleteMolecule needs to be implemented????");
//        drawElement.removeMolecule(tm);
    }


    public final void setValue(StereoMolecule value, boolean b)
    {
        needsLayout(b);
        mMol = value;
        notifyChange();
    }

    public void setValue(Reaction rxn)
    {
        setReaction(rxn);
    }


    public void changed()
    {
        notifyChange();
    }

    public void valueInvalidated()
    {
        for (IValidationListener l : validationListeners) {
            l.valueInvalidated();
        }
    }

    private void notifyChange()
    {
        for (IChangeListener l : changeListeners) {
            l.onChange();
        }
    }


    public StereoMolecule getMoleculeAt(GenericPoint pt, boolean includeBond)
    {
        StereoMolecule mol = mMol;
        {
            if (mol.findAtom(pt.getX(), pt.getY()) != -1) {
                return mol;
            }
/*
            for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
                GenericPoint ap = new GenericPoint(mol.getAtomX(atom), mol.getAtomY(atom));
                if (Math.abs(ap.distance(pt)) < 5) {
//                    System.out.println("getMoleculeAt Atom\n");
                    return mol;
                }
            }
*/
            if (includeBond) {
                for (int i = 0; i < mol.getAllBonds(); i++) {
                    int source = mol.getBondAtom(0, i);
                    int target = mol.getBondAtom(1, i);
                    java.awt.geom.Line2D line =
                        new java.awt.geom.Line2D.Double(mol.getAtomX(source), mol.getAtomY(source), mol.getAtomX(target), mol.getAtomY(target));
                    double dist = line.ptSegDist(pt.getX(), pt.getY());
                    if (dist < 5) {
/*
                        if (Math.abs(new GenericPoint(mol.getAtomX(source), mol.getAtomY(source)).distance(pt)) < 5 &&
                                (Math.abs(new GenericPoint(mol.getAtomX(target), mol.getAtomY(target)).distance(pt)) < 5)) {
                            return mol;
                        }
*/
                        return mol;
                    }
                }
            }
        }
        return null;
    }

    private int getFragmentByAtom(int atom)
    {
        if (atom >= 0 && atom < getMolecule().getAllAtoms()) {
            int idx = mFragmentNo[atom];
            if (idx >= 0 && idx < mFragmentNo.length)
                return idx;
        }
        return -1;
    }

    public void selectFragmentByAtom(int rootAtom)
    {
        int fragment = getFragmentByAtom(rootAtom);
        for (int i = 0; fragment != -1 && i < mMol.getAllAtoms(); i++) {
            if (mFragmentNo[i] == fragment) {
                mMol.setAtomSelection(i, true);
            }
        }
    }


    private boolean isPointOnAtomOrBond(StereoMolecule mol, GenericPoint pt, boolean includeBond)
    {
        for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
            GenericPoint ap = new GenericPoint(mol.getAtomX(atom), mol.getAtomY(atom));
            if (Math.abs(ap.distance(pt)) < 5) {
                return true;
            }
        }
        if (includeBond) {
            for (int i = 0; i < mol.getAllBonds(); i++) {
                int source = mol.getBondAtom(0, i);
                int target = mol.getBondAtom(1, i);
                java.awt.geom.Line2D line =
                    new java.awt.geom.Line2D.Double(
                        mol.getAtomX(source), mol.getAtomY(source),
                        mol.getAtomX(target), mol.getAtomY(target));
                double dist = line.ptSegDist(pt.getX(), pt.getY());
                if (dist < 5) {
/*
                    if (Math.abs(new GenericPoint(mol.getAtomX(source), mol.getAtomY(source)).distance(pt)) < 5 &&
                            (Math.abs(new GenericPoint(mol.getAtomX(target), mol.getAtomY(target)).distance(pt)) < 5)) {
                        return true;
                    }
*/
                    return true;
                }
            }
        }
        return false;
    }

    public StereoMolecule getFragmentAt(GenericPoint pt, boolean includeBond)
    {
        for (StereoMolecule mol : getFragments()) {
            if (isPointOnAtomOrBond(mol, pt, includeBond))
                return mol;
        }
        return null;
    }


    public static int rowFromESRType(int type)
    {
        switch (type) {
            case Molecule.cESRTypeAbs:
                return 0;
            case Molecule.cESRTypeOr:
                return 1;
            case Molecule.cESRTypeAnd:
                return 2;
        }
        return 0;
    }

    public static int esrTypeFromRow(int row)
    {
        switch (row) {
            case 0:
                return Molecule.cESRTypeAbs;
            case 1:
                return Molecule.cESRTypeOr;
            case 2:
                return Molecule.cESRTypeAnd;
        }
        return Molecule.cESRTypeAbs;

    }

    public int getSelectedAtom()
    {
        return selectedAtom;
    }

    public void setSelectedAtom(int theAtom)
    {
        if (selectedAtom != theAtom) {
            if (atomHighlightCallback != null) {
                atomHighlightCallback.onHighlight(theAtom != -1 ? theAtom : selectedAtom, theAtom != -1);
            }
        }
        this.selectedAtom = theAtom;
    }

    public int getSelectedBond()
    {
        return selectedBond;
    }

    public void setSelectedBond(int theBond)
    {
        if (selectedBond != theBond) {
            if (bondHighlightCallback != null) {
                bondHighlightCallback.onHighlight(theBond != -1 ? theBond : selectedBond, theBond != -1);
            }
        }
        this.selectedBond = theBond;
    }


    public final void setMode(int mode)
    {
        this.mMode = mode;
        if ((mMode & (MODE_REACTION | MODE_MARKUSH_STRUCTURE)) != 0) {
            mMode |= (MODE_MULTIPLE_FRAGMENTS);
        }

//        if ((mMode & (MODE_DRAWING_OBJECTS | MODE_REACTION)) != 0) {
//            drawElement.drawingObjects = new ArrayList<DrawingObject>();
//            drawElement.drawingObjects.add(new JFXReactionArrow());
//        }
    }


    public List<IDrawingObject> getDrawingObjects()
    {
        return mDrawingObjectList;
    }

    public void addDrawingObject(IDrawingObject o)
    {
        if (!mDrawingObjectList.contains(o)) {
            pushUndo();
            mDrawingObjectList.add(o);
//            drawElement.drawingObjects.add(o);
        }
    }

    public boolean isReaction()
    {
        return isReactionMode();
    }

    public boolean isFragment()
    {
        boolean fragement = false;
        fragement = mMol.isFragment();
        return fragement;
    }

    public void setFragment(boolean fragment)
    {
        mMol.setFragment(fragment);
        notifyChange();
    }

    public void setNewMolecule()
    {
        StereoMolecule mol = new StereoMolecule();
        mol.setFragment(isFragment());
        setValue(mol, true);
    }


    public void needsLayout(boolean set)
    {
        needslayout = set;
//        if (set)
//            System.out.printf("NeedsLayout %s\n",set);
    }

    public boolean needsLayout()
    {
        return needslayout;
    }

    public int getDisplayMode()
    {
        return displayMode;
    }

    GenericPoint calculateCenter(StereoMolecule r)
    {
        float x = 0;
        float y = 0;
        int atoms = r.getAllAtoms();
        for (int atom = 0; atom < atoms; atom++) {
            x += r.getAtomX(atom);
            y += r.getAtomY(atom);
        }
        return new GenericPoint(x / atoms, y / atoms);
    }


    public String getIDCode()
    {
        if (!isReaction()) {
            StereoMolecule mol = getMolecule();//getSelectedMolecule();
            if (mol != null && mMol.getAllAtoms() > 0) {
                Canonizer can = new Canonizer(mol);
                return (can.getIDCode() + " " + can.getEncodedCoordinates());
            }
        } else {
            Reaction rxn = getReaction();
            String idc = ReactionEncoder.encode(rxn, true, ReactionEncoder.INCLUDE_DEFAULT);
//            if (idc != null)
//                ReactionEncoder.decode(idc, ReactionEncoder.INCLUDE_DEFAULT, null);
            return idc;
        }
        return null;
    }

    public StringBuilder getKeyStrokeBuffer()
    {
        return mAtomKeyStrokeBuffer;
    }


    public int getAtomKeyStrokeValidity(String s)
    {
        if (Molecule.getAtomicNoFromLabel(s) != 0)
            return KEY_IS_ATOM_LABEL;
        if (NamedSubstituents.getSubstituentIDCode(s) != null)
            return KEY_IS_SUBSTITUENT;
        if (isValidAtomKeyStrokeStart(s))
            return KEY_IS_VALID_START;
        return KEY_IS_INVALID;
    }

    /**
     * @param s
     * @return true if s is either a valid atom symbol or a valid substituent name
     */
    private boolean isValidAtomKeyStroke(String s)
    {
        return Molecule.getAtomicNoFromLabel(s) != 0
            || NamedSubstituents.getSubstituentIDCode(s) != null;
    }

    /**
     * @param s
     * @return true if adding one or more chars may still create a valid key stroke sequence
     */
    private boolean isValidAtomKeyStrokeStart(String s)
    {
        if (s.length() < 3)
            for (int i = 1; i < Molecule.cAtomLabel.length; i++) {
                if (Molecule.cAtomLabel[i].startsWith(s))
                    return true;
            }

        return NamedSubstituents.isValidSubstituentNameStart(s);
    }

    public int getMarkushCount()
    {
        return 0;
    }

    public void tryAutoMapReaction()
    {
        SSSearcher sss = new MySSSearcher();

        Reaction rxn = getReaction();//new Reaction(reaction);

        // Mark the manually mapped atoms, so we may re-assign them later
        for (int i = 0; i < rxn.getMolecules(); i++) {
            StereoMolecule mol = rxn.getMolecule(i);
            for (int a = 0; a < mol.getAtoms(); a++) {
                if (mol.getAtomMapNo(a) > 0) {
                    mol.setAtomicNo(a, FAKE_ATOM_NO + mol.getAtomMapNo(a));
                }
            }
        }
        rxn = mMapper.mapReaction(rxn, sss);
        if (rxn != null) {
            int offset = 0;
            // Sync the display molecule with the reaction fragments
            {
                int[] fragmentAtom = new int[mFragment.length];
                for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                    int fragment = mFragmentNo[atom];
                    if (mFragment[fragment].getAtomicNo(fragmentAtom[fragment]) > FAKE_ATOM_NO) {
                        mMol.setAtomMapNo(atom, mFragment[fragment].getAtomicNo(fragmentAtom[fragment]) - FAKE_ATOM_NO, false);
                        offset = Math.max(mMol.getAtomMapNo(atom), offset);
                    }
                    fragmentAtom[fragment]++;
                }
            }
            {
                int[] fragmentAtom = new int[mFragment.length];
                for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                    int fragment = mFragmentNo[atom];
                    if (mFragment[fragment].getAtomMapNo(fragmentAtom[fragment]) > 0 && (mFragment[fragment].getAtomicNo(fragmentAtom[fragment]) <= FAKE_ATOM_NO)) {
                        mMol.setAtomMapNo(atom, mFragment[fragment].getAtomMapNo(fragmentAtom[fragment]) + offset, true);
                    }
                    fragmentAtom[fragment]++;
                }
            }
        }
        syncFragments();
    }


    public String getMolFile(boolean v3)
    {
        if (v3)
            return new MolfileV3Creator(mMol).getMolfile();
        else
            return new MolfileCreator(mMol).getMolfile();
    }

    public void setMolFile(String molFile)
    {
        try {
            MolfileParser p = new MolfileParser();
            StereoMolecule mol = new StereoMolecule();
            p.parse(mol, molFile);
            setValue(mol, true);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public String getSmiles()
    {
        return new IsomericSmilesCreator(mMol).getSmiles();
    }

    public void setSmiles(String smiles)
    {
        try {
            SmilesParser p = new SmilesParser();
            StereoMolecule mol = new StereoMolecule();
            p.parse(mol, smiles);
            setValue(mol, true);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int getReactantCount()
    {
//        System.out.printf("getReactantCount Count %s\n",mReactantCount);
        return mReactantCount;
    }


    static class MySSSearcher extends SSSearcher
    {
        @Override
        public boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom)
        {
            if (mMolecule.getAtomicNo(moleculeAtom) == mFragment.getAtomicNo(fragmentAtom))
                if (mMolecule.isAromaticAtom(moleculeAtom) || mFragment.isAromaticAtom(fragmentAtom))
                    return true;
            return super.areAtomsSimilar(moleculeAtom, fragmentAtom);
        }

        @Override
        public boolean areBondsSimilar(int moleculeBond, int fragmentBond)
        {
            if (mMolecule.isAromaticBond(moleculeBond) || mMolecule.isDelocalizedBond(moleculeBond) ||
                mFragment.isAromaticBond(fragmentBond) || mFragment.isDelocalizedBond(fragmentBond)
                )
                return true;
            return super.areBondsSimilar(moleculeBond, fragmentBond);
            //return true;
        }
    }

    private GenericPoint calculateCenterOfGravity()
    {
        int atoms = mMol.getAllAtoms();
        double sumx = 0;
        double sumy = 0;
        for (int atom = 0; atom < atoms; atom++) {
            sumx += mMol.getAtomX(atom);
            sumy += mMol.getAtomY(atom);
        }
        return atoms > 0 ? new GenericPoint(sumx / atoms, sumy / atoms) : null;
    }

    public void flip(boolean horiz)
    {
        GenericPoint pt = calculateCenterOfGravity();
        if (pt != null) {
            // center
            moveCoords((float) -pt.getX(), (float) -pt.getY());
            if (horiz) {
                scaleCoords(-1, 1);
            } else {
                scaleCoords(1, -1);
            }
            moveCoords((float) pt.getX(), (float) pt.getY());

            // invert stereo bonds
            for (int bond = 0; bond < mMol.getAllBonds(); bond++) {
                if (mMol.getBondType(bond) == Molecule.cBondTypeUp)
                    mMol.setBondType(bond, Molecule.cBondTypeDown);
                else if (mMol.getBondType(bond) == Molecule.cBondTypeDown)
                    mMol.setBondType(bond, Molecule.cBondTypeUp);
            }
        }
    }

    private void scaleCoords(float scalex, float scaley)
    {
        int atoms = mMol.getAllAtoms();
        for (int atom = 0; atom < atoms; atom++) {
            mMol.setAtomX(atom, mMol.getAtomX(atom) * scalex);
            mMol.setAtomY(atom, mMol.getAtomY(atom) * scaley);
        }
    }

    private void moveCoords(float cx, float cy)
    {
        int atoms = mMol.getAllAtoms();
        for (int atom = 0; atom < atoms; atom++) {
            mMol.setAtomX(atom, mMol.getAtomX(atom) + cx);
            mMol.setAtomY(atom, mMol.getAtomY(atom) + cy);
        }
    }


    private AtomHighlightCallback atomHighlightCallback = null;
    private BondHighlightCallback bondHighlightCallback = null;

    public void registerAtomHighlightCallback(AtomHighlightCallback cb)
    {
        atomHighlightCallback = cb;
    }

    public void registerBondHighlightCallback(BondHighlightCallback cb)
    {
        bondHighlightCallback = cb;
    }


    public void addMolecule(StereoMolecule mol,double x ,double y)
    {
        if (mol != null && mol.getAllAtoms() != 0) {
            if (mMol.getAllAtoms() == 0) {
                int avbl = 0;
                boolean isFragment = mMol.isFragment();
                scaleIntoView(mol, avbl,0,0);
                mol.copyMolecule(mMol);
                mMol.setFragment(isFragment);
                notifyChange();
            } else {
                int avbl = (int) mMol.getAverageBondLength();
                scaleIntoView(mol, avbl,x,y);
                int originalAtoms = mMol.getAllAtoms();
                boolean isFragment = mMol.isFragment();
                mMol.addMolecule(mol);
                for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
                    mMol.setAtomSelection(atom, atom >= originalAtoms);
                }
                mMol.setFragment(isFragment);
                notifyChange();
            }
        }

    }

    private void scaleIntoView(StereoMolecule mol, int avbl,double cx,double cy)
    {
        AbstractDepictor d = createDepictor(mol);
        DepictorTransformation dt = d.simpleValidateView(new GenericRectangle(0, 0, this.getWidth(), this.getHeight()), AbstractDepictor.cModeInflateToMaxAVBL + avbl);

        if (dt != null) {
            dt.move(cx,cy);
            System.out.printf("Transform %s %s\n",dt.getOffsetX(),cx);
            dt.applyTo(mol);
        }
    }

    public IDrawingObject getSelectedDrawingObject()
    {
        return selectedDrawingObject;
    }

    public void setSelectedDrawingObject(IDrawingObject sel)
    {
        if (selectedDrawingObject != sel) {
            setSelectedAtom(-1);
            setSelectedBond(-1);
        }
        this.selectedDrawingObject = sel;
//        for (IDrawingObject d : getDrawingObjects())
//            d.setSelected(false);
        if (sel != null)
            this.selectedDrawingObject.setSelected(true);
    }

/*    public void cleanupCoordinates()
    {
        AbstractDepictor depictor = createDepictor(getMolecule());
        DepictorTransformation dt = depictor.simpleValidateView(new Rectangle2D.Double(0, 0, this.getWidth(), this.getHeight()), AbstractDepictor.cModeInflateToMaxAVBL);
        if (dt != null)
            dt.applyTo(mMol);

    }*/

    public void cleanMolecule(boolean invent, boolean selectedOnly)
    {
        cleanupCoordinates(false, invent);
        valueInvalidated();
    }

    //    public abstract void cleanMolecule(boolean selectedOnly);
    protected abstract AbstractExtendedDepictor createExtendedDepictor();

    protected abstract AbstractDepictor createDepictor(StereoMolecule stereoMolecule);

    public abstract void analyzeReaction();

    public abstract boolean copyMolecule(boolean selected);

    public abstract boolean copyReaction(boolean selected);

    public abstract StereoMolecule pasteMolecule(double cx,double cy);

    public abstract Reaction pasteReaction(double cx,double cy);

    public ImageProvider getImageProvider()
    {
        return imageProvider;
    }


    public void setImageProvider(ImageProvider p)
    {
        imageProvider = p;
    }


/*
    private void cleanupMultiFragmentCoordinates(ExtendedDepictor depictor)
   	{
   		Rectangle2D.Double[] boundingRect = new Rectangle2D.Double[mFragment.length];
   //		float fragmentWidth = 0.0f;
   		for (int fragment = 0; fragment < mFragment.length; fragment++) {
   			Depictor d = new Depictor(mFragment[fragment]);
   			d.updateCoords(null, null, AbstractDepictor.cModeInflateToMaxAVBL);
   			boundingRect[fragment] = d.getBoundingRect();
   //			fragmentWidth += boundingRect[fragment].width;
   		}

   		double spacing = FRAGMENT_CLEANUP_DISTANCE * AbstractDepictor.cOptAvBondLen;
   		double avbl = mMol.getAverageBondLength();
   		double arrowWidth = DEFAULT_ARROW_LENGTH * getWidth();

   		double rawX = 0.5 * spacing;
   		for (int fragment = 0; fragment <= mFragment.length; fragment++) {
   			if (isReactionMode() && fragment == mReactantCount) {
   				 mDrawingObjectList.get(0).setRect(
                    (float)(rawX - spacing / 2),
                    getHeight() / 2,
                    (float)(rawX - spacing / 2 + arrowWidth),
                    getHeight() / 2);
   				rawX += arrowWidth;
   			}

   			if (fragment == mFragment.length) {
   				break;
   			}

   			double dx = rawX - boundingRect[fragment].x;
   			double dy = 0.5 * (getHeight() - boundingRect[fragment].height)
   				- boundingRect[fragment].y;
   			mFragment[fragment].translateCoords(dx, dy);

   			rawX += spacing + boundingRect[fragment].width;
   		}

   		depictor.updateCoords(null, new Rectangle2D.Double(0, 0, getWidth(), getHeight()), maxUpdateMode());

   		int[] fragmentAtom = new int[mFragment.length];
   		for (int atom = 0; atom < mMol.getAllAtoms(); atom++) {
   			int fragment = mFragmentNo[atom];
   			mMol.setAtomX(atom, mFragment[fragment].getAtomX(fragmentAtom[fragment]));
   			mMol.setAtomY(atom, mFragment[fragment].getAtomY(fragmentAtom[fragment]));
   			fragmentAtom[fragment]++;
   		}

//   		mMol.setStereoBondsFromParity();
   	}
*/


}
