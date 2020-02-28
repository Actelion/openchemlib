/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import java.util.ArrayList;
import java.util.Collection;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

public abstract class DefaultCompoundCollectionModel<T> implements CompoundCollectionModel<T> {
    private ArrayList<T> mCompoundList;
    private ArrayList<CompoundCollectionListener> mListenerList;
    private boolean mEnableEvents;

    public DefaultCompoundCollectionModel() {
        mCompoundList = new ArrayList<T>();
        mListenerList = new ArrayList<>();
        mEnableEvents = true;
        }

    public void addCompound(T compound) {
        int size = mCompoundList.size();
        mCompoundList.add(compound);
        informListeners(size, size);
        }

    public void addCompound(int index, T compound) {
        int size = mCompoundList.size();
        mCompoundList.add(index, compound);
        informListeners(index, size);
        }

    public void addCompoundList(Collection<T> list) {
    	if (list.size() != 0) {
    		for (T compound:list)
    			mCompoundList.add(compound);
	        informListeners(mCompoundList.size() - list.size(), list.size() - 1);
    		}
        }

    public void setCompound(int index, T compound) {
        mCompoundList.set(index, compound);
        informListeners(index, index);
        }

    public void setCompoundList(Collection<T> list) {
        int size = mCompoundList.size() - 1;
        mCompoundList.clear();
		for (T compound:list)
			mCompoundList.add(compound);
        informListeners(0, Math.max(size, list.size() - 1));
        }

    public void addCompoundCollectionListener(CompoundCollectionListener l) {
        mListenerList.add(l);
        }

    public T getCompound(int index) {
        return mCompoundList.get(index);
        }

    public void removeCompoundCollectionListener(CompoundCollectionListener l) {
        mListenerList.remove(l);
        }

    public void clear() {
        int size = mCompoundList.size() - 1;
        mCompoundList.clear();
        informListeners(0, size);
        }

    public void addMoleculeList(Collection<StereoMolecule> list) {
    	mEnableEvents = false;
    	int oldSize = mCompoundList.size();
    	for (StereoMolecule mol:list)
    		addMolecule(mCompoundList.size(), mol);
    	mEnableEvents = true;
        informListeners(oldSize, mCompoundList.size() - 1);
    	}

    abstract public void addMolecule(int index, StereoMolecule mol);
    abstract public void setMolecule(int index, StereoMolecule mol);
    abstract public StereoMolecule getMolecule(int index);

    public StereoMolecule getMoleculeForDisplay(int index) {
        return getMolecule(index);
        }

    public int getSize() {
        return mCompoundList.size();
        }

    public void remove(int index) {
        mCompoundList.remove(index);
        informListeners(index, mCompoundList.size());
        }

    private void informListeners(int fromIndex, int toIndex) {
    	if (mEnableEvents)
    		for (CompoundCollectionListener l:mListenerList)
    			l.collectionUpdated(fromIndex, toIndex);
        }

    /**
     * This version of the DefaultCompoundCollectionModel collects molecules as StereoMolecules.
     * It is the preferred model when the number of handled molecules
     * is limited and when specific molecule features beyond idcode, atom coordinates, and molecule
	 * name must not get lost.
     */
    public static class Molecule extends DefaultCompoundCollectionModel<StereoMolecule> {
        public StereoMolecule getMolecule(int index) {
            return getCompound(index);
            }

        public void setMolecule(int index, StereoMolecule mol) {
            setCompound(index, mol);
            }

        public void addMolecule(int index, StereoMolecule mol) {
            addCompound(index, mol);
            }
        }

    /**
     * This version of the DefaultCompoundCollectionModel collects
     * molecules as IDCodes. The molecule access functions result
     * in appropriate conversion.
     * It is the preferred model when the number of handled molecules
     * is potentially high and when atom selection or atom colors need to be retained.
     */
    public static class IDCode extends DefaultCompoundCollectionModel<String> {
        public StereoMolecule getMolecule(int index) {
            return new IDCodeParser(true).getCompactMolecule(getCompound(index));
            }

        public void setMolecule(int index, StereoMolecule mol) {
        	Canonizer c = new Canonizer(mol);
            setCompound(index, c.getIDCode()+" "+c.getEncodedCoordinates());
            }

        public void addMolecule(int index, StereoMolecule mol) {
        	Canonizer c = new Canonizer(mol);
            addCompound(index, c.getIDCode()+" "+c.getEncodedCoordinates());
            }
        }

    /**
     * This version of the DefaultCompoundCollectionModel collects
     * molecules as String[2] with idcodes & idcoords (index 0) and molecule name (index 1).
	 * The molecule access functions result in appropriate conversion.
     * It is the preferred model when the number of handled molecules
     * is potentially high and when the molecule name/ID needs to be retained.
     */
    public static class IDCodeWithName extends DefaultCompoundCollectionModel<String[]> {
        public StereoMolecule getMolecule(int index) {
			StereoMolecule mol = new IDCodeParser().getCompactMolecule(getCompound(index)[0]);
			mol.setName(getCompound(index)[1]);
            return mol;
        }

        public void setMolecule(int index, StereoMolecule mol) {
			String[] idcodeWithName = new String[2];
            Canonizer c = new Canonizer(mol);
			idcodeWithName[0] = c.getIDCode().concat(" ").concat(c.getEncodedCoordinates());
			idcodeWithName[1] = mol.getName();
            setCompound(index, idcodeWithName);
        }

        public void addMolecule(int index, StereoMolecule mol) {
			String[] idcodeWithName = new String[2];
            Canonizer c = new Canonizer(mol);
			idcodeWithName[0] = c.getIDCode().concat(" ").concat(c.getEncodedCoordinates());
			idcodeWithName[1] = mol.getName();
            addCompound(index, idcodeWithName);
        }
    }

    /**
     * This version of the DefaultCompoundCollectionModel collects
     * IDCodes and StereoMolecules without conversion into their native type.
     * It is the preferred model if the molecule source type is not predictable.
     */
    public static class Native extends DefaultCompoundCollectionModel<Object> {
        public StereoMolecule getMolecule(int index) {
            Object o = getCompound(index);
            return (o instanceof StereoMolecule) ? (StereoMolecule)o
                 : (o instanceof String) ? new IDCodeParser(true).getCompactMolecule((String)o)
                 : null;
            }

        public void setMolecule(int index, StereoMolecule mol) {
            setCompound(index, mol);
            }

        public void addMolecule(int index, StereoMolecule mol) {
            addCompound(index, mol);
            }
        }
    }
