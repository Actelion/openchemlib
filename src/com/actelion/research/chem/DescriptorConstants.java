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

public interface DescriptorConstants {
    public static final int DESCRIPTOR_TYPE_UNKNOWN = -1;
    public static final int DESCRIPTOR_TYPE_MOLECULE = 1;
    public static final int DESCRIPTOR_TYPE_REACTION = 2;

    public static final DescriptorInfo DESCRIPTOR_FFP512 = 
                            new DescriptorInfo("FragmentFingerprint512",
                                               "FragFp",
                                               DESCRIPTOR_TYPE_MOLECULE,
                                               true,
                                               true,
                                               false);
    public static final DescriptorInfo DESCRIPTOR_PFP512 = 
                            new DescriptorInfo("PathFingerprint512",
                                               "PathFp",
                                               DESCRIPTOR_TYPE_MOLECULE,
                                               true,
                                               true,
                                               false);
    public static final DescriptorInfo DESCRIPTOR_HashedCFp = 
                            new DescriptorInfo("HashedSphericalFingerprint512",
                                               "SphereFp",
                                               DESCRIPTOR_TYPE_MOLECULE,
                                               true,
                                               true,
                                               true);	// for the creation of up/down bonds
    public static final DescriptorInfo DESCRIPTOR_SkeletonSpheres = 
                                               new DescriptorInfo("HashedSkeletonSphereCount1024",
                                               "SkelSpheres",
                                               DESCRIPTOR_TYPE_MOLECULE,
                                               false,
                                               true,
                                               true);	// for the creation of up/down bonds
    public static final DescriptorInfo DESCRIPTOR_OrganicFunctionalGroups = 
    										   new DescriptorInfo("FunctionalGroupTreeCount1024",
    										   "OrgFunctions",
    										   DESCRIPTOR_TYPE_MOLECULE,
    										   false,
    										   false,
    										   true);	// for the creation of up/down bonds
    public static final DescriptorInfo DESCRIPTOR_CenteredSkeletonFragments = 
										        new DescriptorInfo("CenteredSkeletonFragments",
										        "CentSkelFrags",
										        DESCRIPTOR_TYPE_MOLECULE,
										        false,
										        true,
										        true);	// for the creation of up/down bonds
    public static final DescriptorInfo DESCRIPTOR_TopoPPHistDist = 
                            new DescriptorInfo("TopologicalPharmacophoreHistograms",
                                               "TopPPHist",
                                               DESCRIPTOR_TYPE_MOLECULE,
                                               false,
                                               false,
                                               false);
    public static final DescriptorInfo DESCRIPTOR_Flexophore = 
                            new DescriptorInfo("Flexophore",
                                               "Flexophore",
                                               DESCRIPTOR_TYPE_MOLECULE,
                                               false,
                                               false,
                                               false);
    public static final DescriptorInfo DESCRIPTOR_Flexophore_HighRes = 
        					new DescriptorInfo("FlexophoreHighResolution",
        									   "FlexophoreHighRes",
        									   DESCRIPTOR_TYPE_MOLECULE,
        									   false,
        									   false,
        									   false);
    public static final DescriptorInfo DESCRIPTOR_ReactionIndex = 
                            new DescriptorInfo("ReactionIndex",
                                               "RxnIdx",
                                               DESCRIPTOR_TYPE_REACTION,
                                               false,
                                               false,
                                               false);
    public static final DescriptorInfo DESCRIPTOR_IntegerVector = 
    						new DescriptorInfo("IntegerVector",
    											"IntVec",
    											DESCRIPTOR_TYPE_UNKNOWN,
    											false,
    											false,
    											false);
   
    public static final DescriptorInfo DESCRIPTOR_MAX_COMMON_SUBSTRUCT = 
        					new DescriptorInfo("MaximumCommonSubstructure",
        										"MCS",
        										DESCRIPTOR_TYPE_MOLECULE,
        										false,
        										true,
        										false);

    public static final DescriptorInfo DESCRIPTOR_SUBSTRUCT_QUERY_IN_BASE = 
							new DescriptorInfo("SubStructureQueryInBase",
								"SSSQinB",
								DESCRIPTOR_TYPE_MOLECULE,
								false,
								false, // ??? TODO check
								false);
    
    public static final DescriptorInfo DESCRIPTOR_FULL_FRAGMENT_SET = 
            new DescriptorInfo("FullFragmentSet",
                               "FullFragSet",
                               DescriptorConstants.DESCRIPTOR_TYPE_MOLECULE,
                               true,
                               true,
                               false);

    public static final DescriptorInfo[] DESCRIPTOR_LIST = {
                                                DESCRIPTOR_FFP512,
                                                DESCRIPTOR_PFP512,
                                                DESCRIPTOR_HashedCFp,
                                                DESCRIPTOR_SkeletonSpheres,
                                                DESCRIPTOR_OrganicFunctionalGroups,
                                                DESCRIPTOR_Flexophore
                                                };

    public static final DescriptorInfo[] DESCRIPTOR_EXTENDED_LIST = {
                                                DESCRIPTOR_FFP512,
                                                DESCRIPTOR_PFP512,
                                                DESCRIPTOR_HashedCFp,
                                                DESCRIPTOR_SkeletonSpheres,
                                                DESCRIPTOR_CenteredSkeletonFragments,
                                                DESCRIPTOR_FULL_FRAGMENT_SET,
                                                DESCRIPTOR_MAX_COMMON_SUBSTRUCT,
                                                DESCRIPTOR_SUBSTRUCT_QUERY_IN_BASE,
                                                DESCRIPTOR_TopoPPHistDist,
                                                DESCRIPTOR_OrganicFunctionalGroups,
                                                DESCRIPTOR_Flexophore,
                                                DESCRIPTOR_Flexophore_HighRes,
                                                DESCRIPTOR_ReactionIndex,
                                                DESCRIPTOR_IntegerVector,
                                                DESCRIPTOR_FULL_FRAGMENT_SET
                                                };
    }

