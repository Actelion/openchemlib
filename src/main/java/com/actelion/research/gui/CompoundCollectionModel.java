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
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.chem.StereoMolecule;

import java.util.Collection;

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
