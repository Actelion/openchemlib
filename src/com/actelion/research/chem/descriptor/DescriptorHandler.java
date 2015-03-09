/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
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

package com.actelion.research.chem.descriptor;

/**
 * DescriptorHandler is the base Interface for any molecular descriptor
 * handling StereoMolecules or Reactions
 */

public interface DescriptorHandler<T, U extends Object> extends ISimilarityCalculator<T> {
    static final String FAILED_STRING = "Calculation Failed";
    static final byte[] FAILED_BYTES = FAILED_STRING.getBytes();
    public abstract DescriptorInfo getInfo();
    public abstract String getVersion();
	public abstract String encode(T o);
    public abstract T decode(String s);
    public abstract T decode(byte[] bytes);
    public abstract T createDescriptor(U chemObject);
    public abstract boolean calculationFailed(T o);
    public abstract DescriptorHandler<T,U> getDeepCopy();
    }
