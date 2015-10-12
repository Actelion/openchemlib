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

public class DescriptorHandlerStandard2DFactory implements DescriptorConstants,DescriptorHandlerFactory {
	protected static DescriptorHandlerStandard2DFactory sFactory;

	public static DescriptorHandlerFactory getFactory() {
		if (sFactory == null) {
			synchronized(DescriptorHandlerStandard2DFactory.class) {
				if (sFactory == null)
					sFactory = new DescriptorHandlerStandard2DFactory();
				}
			}
		return sFactory;
		}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public DescriptorHandler getDefaultDescriptorHandler(String shortName) {

		if (DESCRIPTOR_FFP512.shortName.equals(shortName))
			return DescriptorHandlerFFP512.getDefaultInstance();
		if (DESCRIPTOR_PFP512.shortName.equals(shortName))
			return DescriptorHandlerPFP512.getDefaultInstance();
		if (DESCRIPTOR_HashedCFp.shortName.equals(shortName))
			return DescriptorHandlerHashedCFp.getDefaultInstance();
		if (DESCRIPTOR_SkeletonSpheres.shortName.equals(shortName))
			return DescriptorHandlerSkeletonSpheres.getDefaultInstance();
		if (DESCRIPTOR_OrganicFunctionalGroups.shortName.equals(shortName))
			return DescriptorHandlerFunctionalGroups.getDefaultInstance();

        if (DESCRIPTOR_ReactionIndex.shortName.equals(shortName))
            return DescriptorHandlerReactionIndex.getDefaultInstance();

		return null;
		}

    @SuppressWarnings({ "unchecked", "rawtypes" })
	public DescriptorHandler create(String shortName) {

		if (DESCRIPTOR_FFP512.shortName.equals(shortName))
			return new DescriptorHandlerFFP512();
		if (DESCRIPTOR_PFP512.shortName.equals(shortName))
			return new DescriptorHandlerPFP512();
		if (DESCRIPTOR_HashedCFp.shortName.equals(shortName))
			return new DescriptorHandlerHashedCFp();
		if (DESCRIPTOR_SkeletonSpheres.shortName.equals(shortName))
			return new DescriptorHandlerSkeletonSpheres();
		if (DESCRIPTOR_OrganicFunctionalGroups.shortName.equals(shortName))
			return new DescriptorHandlerFunctionalGroups();

        if (DESCRIPTOR_ReactionIndex.shortName.equals(shortName))
            return new DescriptorHandlerReactionIndex();

		return null;
		}
	}
