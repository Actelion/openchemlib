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

import java.util.Collection;

import com.actelion.research.chem.StereoMolecule;

public interface CompoundCollectionModel<T> {
    // A CompoundCollectionModel maintains a list of opaque objects which somehow
    // contain chemical structures. The default implementation expects these
    // objects to be idcodes or StereoMolecules. addMolecule(),setMolecule() and
    // getMolecule() must provide a direct access to StereoMolecules, no matter what
    // kind of object is used in the list.
    public void addCompoundCollectionListener(CompoundCollectionListener l);
    public void removeCompoundCollectionListener(CompoundCollectionListener l);
    public void clear();
    public void remove(int index);
    public int getSize();

    // Methods using StereoMolecule
    public void addMolecule(int index, StereoMolecule mol);
    public void addMoleculeList(Collection<StereoMolecule> list);
    public void setMolecule(int index, StereoMolecule mol);
    public StereoMolecule getMolecule(int index);
    public StereoMolecule getMoleculeForDisplay(int index);

    // Methods using the native object
    public void addCompound(T compound);
    public void addCompound(int index, T compound);
    public void addCompoundList(Collection<T> list);
    public T getCompound(int i);
    public void setCompound(int index, T compound);
    public void setCompoundList(Collection<T> list);
    }
