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

package com.actelion.research.chem.descriptor;

import com.actelion.research.chem.SSSearcherWithIndex;

public interface DescriptorConstants {
	int DESCRIPTOR_TYPE_UNKNOWN = -1;
	int DESCRIPTOR_TYPE_MOLECULE = 1;
	int DESCRIPTOR_TYPE_REACTION = 2;

	DescriptorInfo DESCRIPTOR_FFP512 =
                            new DescriptorInfo("FragmentFingerprint512",
												"FragFp",
												SSSearcherWithIndex.cIndexVersion,
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												false);
	DescriptorInfo DESCRIPTOR_SSSFP =
							new DescriptorInfo("SSSPathFingerprint",
												"SSSPathFp",
												"1.0",
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												false);
	DescriptorInfo DESCRIPTOR_PFP512 =
                            new DescriptorInfo("PathFingerprint512",
												"PathFp",
												"1.1",
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												false);
	DescriptorInfo DESCRIPTOR_HashedCFp =
                            new DescriptorInfo("HashedSphericalFingerprint512",
												"SphereFp",
												"2.1",
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												true);	// for the creation of up/down bonds
	DescriptorInfo DESCRIPTOR_SkeletonSpheres =
							new DescriptorInfo("HashedSkeletonSphereCount1024",
												"SkelSpheres",
												"1.1",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												true,
												true,
												true);	// for the creation of up/down bonds
	DescriptorInfo DESCRIPTOR_ALLFRAG =
							new DescriptorInfo("AllFragments1024",
												"AllFragFp",
												"1.0",
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												false);
	DescriptorInfo DESCRIPTOR_OrganicFunctionalGroups =
												new DescriptorInfo("FunctionalGroupTreeCount1024",
												"OrgFunctions",
												"1.0",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												false,
												false,
												true);	// for the creation of up/down bonds
    DescriptorInfo DESCRIPTOR_CenteredSkeletonFragments =
												new DescriptorInfo("CenteredSkeletonFragments",
												"CentSkelFrags",
												"1.0",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												false,
												true,
												true);	// for the creation of up/down bonds
    DescriptorInfo DESCRIPTOR_TopoPPHistDist =
                            new DescriptorInfo("TopologicalPharmacophoreHistograms",
												"TopPPHist",
												"version",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												false,
												false,
												false);
    DescriptorInfo DESCRIPTOR_Flexophore =
							new DescriptorInfo("Flexophore",
												"Flexophore",
												"5.0",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												false,
												false,
												false);

	DescriptorInfo DESCRIPTOR_ShapeAlign =
        					new DescriptorInfo("PharmacophoreEnhancedShapeAlignment",
        									   "PheSA",
												"2.1",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												false,
												false,
												false);

    DescriptorInfo DESCRIPTOR_ShapeAlignSingleConf =
        					new DescriptorInfo("PharmacophoreEnhancedShapeAlignmentSingleConfQuery",
        									   "PheSASingleConf",
												"2.1",
												DESCRIPTOR_TYPE_MOLECULE,
        									   false,
												false,
												false,
												false);

	DescriptorInfo DESCRIPTOR_ReactionFP =
							new DescriptorInfo("ReactionFingerprint",
												"RxnFP",
												DescriptorHandlerReactionFP.cVersion,
												DESCRIPTOR_TYPE_REACTION,
												false,	// it is a binary, but similarity calculations have a bias on reaction center bits
												false,
												false,
												false);
    DescriptorInfo DESCRIPTOR_IntegerVector =
    						new DescriptorInfo("IntegerVector",
    											"IntVec",
												"1.0",
    											DESCRIPTOR_TYPE_UNKNOWN,
    											false,
												true,
    											false,
    											false);
   
    DescriptorInfo DESCRIPTOR_MAX_COMMON_SUBSTRUCT =
        					new DescriptorInfo("MaximumCommonSubstructure",
        										"Structure",
												"1.0",
        										DESCRIPTOR_TYPE_MOLECULE,
        										false,
												false,
        										true,
        										false);

    DescriptorInfo DESCRIPTOR_SUBSTRUCT_QUERY_IN_BASE =
							new DescriptorInfo("SubStructureQueryInBase",
												"SSSQinB",
												"1.0",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												false,
												false, // ??? TODO check
												false);
    
    DescriptorInfo DESCRIPTOR_FULL_FRAGMENT_SET =
							new DescriptorInfo("FullFragmentSet",
												"FullFragSet",
												"1.0",
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												false);
    
    DescriptorInfo DESCRIPTOR_PhysicoChemicalProperties =
							new DescriptorInfo("DescriptorPhysicoChemicalProperties",
												"PhysChem",
												"version",
												DESCRIPTOR_TYPE_MOLECULE,
												false,
												true,
												false,
												false);

    DescriptorInfo DESCRIPTOR_BINARY_SKELETONSPHERES =
							new DescriptorInfo("BinarySkeletonSpheres",
												"BinSkelSpheres",
												"10052017",
												DESCRIPTOR_TYPE_MOLECULE,
												true,
												true,
												true,
												false);
    DescriptorInfo DESCRIPTOR_PTREE =
			new DescriptorInfo("PharmacophoreTree",
							   "PTree",
								"1.0",
								DESCRIPTOR_TYPE_MOLECULE,
								false,
								false,
								false,
								false);


    DescriptorInfo[] DESCRIPTOR_LIST = {
                                                DESCRIPTOR_FFP512,
                                                DESCRIPTOR_PFP512,
											    DESCRIPTOR_ALLFRAG,
                                                DESCRIPTOR_HashedCFp,
                                                DESCRIPTOR_SkeletonSpheres,
                                                DESCRIPTOR_OrganicFunctionalGroups,
                                                DESCRIPTOR_Flexophore,
												DESCRIPTOR_ReactionFP
                                                };

    DescriptorInfo[] DESCRIPTOR_EXTENDED_LIST = {
                                                DESCRIPTOR_FFP512,
                                                DESCRIPTOR_PFP512,
											    DESCRIPTOR_SSSFP,
											    DESCRIPTOR_ALLFRAG,
                                                DESCRIPTOR_HashedCFp,
                                                DESCRIPTOR_SkeletonSpheres,
												DESCRIPTOR_BINARY_SKELETONSPHERES,
                                                DESCRIPTOR_CenteredSkeletonFragments,
                                                DESCRIPTOR_MAX_COMMON_SUBSTRUCT,
                                                DESCRIPTOR_SUBSTRUCT_QUERY_IN_BASE,
                                                DESCRIPTOR_TopoPPHistDist,
                                                DESCRIPTOR_OrganicFunctionalGroups,
                                                DESCRIPTOR_Flexophore,
												DESCRIPTOR_ShapeAlign,
												DESCRIPTOR_ReactionFP,
                                                DESCRIPTOR_IntegerVector,
                                                DESCRIPTOR_FULL_FRAGMENT_SET,
                                                DESCRIPTOR_PhysicoChemicalProperties,
                                                DESCRIPTOR_BINARY_SKELETONSPHERES
                                                };
    }

