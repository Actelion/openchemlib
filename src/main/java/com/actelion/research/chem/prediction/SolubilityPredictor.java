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

package com.actelion.research.chem.prediction;

import com.actelion.research.chem.AtomTypeCalculator;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.text.DecimalFormat;
import java.text.NumberFormat;


public class SolubilityPredictor {
	public static final float cSolubilityUnknown = -999;

	private static final long[] cAtomType = {
		           0x80002L,           0x80004L,           0x80009L,           0x8000dL,
		           0x90002L,           0x90003L,           0x90004L,           0x90005L,
		           0x90008L,           0x90009L,           0x9000dL,           0x9000eL,
		           0xa8002L,           0xa8003L,           0xa8004L,           0xa8009L,
		           0xa800dL,           0xc8002L,           0xc8003L,           0xc8004L,
		           0xc8005L,           0xc8008L,           0xc8009L,          0x110002L,
		          0x110003L,          0x110004L,          0x110008L,          0x128004L,
		          0x148004L,          0x148008L,          0x190002L,          0x190003L,
		         0x1090002L,         0x1090003L,         0x1090004L,         0x1090005L,
		         0x1090008L,         0x1090009L,         0x109000dL,         0x109000eL,
		         0x10a8002L,        0x40080004L,        0x48080002L,        0x48080004L,
		        0x48090002L,        0x48090003L,        0x48090004L,        0x48090008L,
		        0x4809000cL,        0x54080002L,        0x54090002L,        0x54090003L,
		        0x54090004L,        0x54090008L,        0x540a8002L,        0x540a8008L,
		        0x540a800cL,        0x64090002L,        0x64090003L,        0x64090004L,
		        0x64090008L,        0x6409000cL,        0x640a8002L,        0x640a8003L,
		        0x640a8004L,        0x640c8002L,        0x640c8003L,        0x640c8004L,
		        0x640c8008L,        0x88090002L,        0x88090003L,        0x880a8002L,
		        0x880a8003L,        0x880c8002L,        0x880c8003L,        0x94090002L,
		        0x94090003L,        0x940a8002L,        0x940a8003L,        0xa4128002L,
		        0xc8090002L,        0xd4090002L,        0xd40a8002L,        0xd40c8002L,
		       0x809010002L,       0x809010003L,       0x809010004L,       0x809010008L,
		       0x815010002L,       0x815010003L,       0x815010004L,       0x815010008L,
		       0x815028002L,       0x815028003L,       0x815028004L,       0x825010002L,
		       0x825010003L,       0x825028002L,       0x848090002L,       0x848090003L,
		       0x848090004L,       0x848090008L,       0x8480a8002L,       0x8480a8003L,
		       0x8480a8004L,       0x8480c8002L,       0x8480c8003L,       0x8480c8004L,
		       0x8480c8008L,       0x848110002L,       0x848110003L,       0x848128002L,
		       0x848128003L,       0x84812800eL,       0x848190002L,       0x849090002L,
		       0x849090003L,       0x849090004L,       0x849090008L,       0x854080002L,
		       0x854090002L,       0x8540a8002L,       0x8540c8002L,       0x855090002L,
		     0x24048090002L,     0x24048090003L,     0x2a048090002L,     0x2a048090003L,
		     0x2a054090002L,     0x2a054090003L,     0x2a0540a8002L,     0x32048090002L,
		     0x32048090003L,     0x32054090002L,     0x320540a8002L,     0x32064090002L,
		     0x320640c8002L,     0x320640c8003L,     0x44048090002L,     0x44054090002L,
		     0x44054090003L,     0x440540a8002L,     0x44064090002L,     0x440640a8002L,
		     0x440640c8002L,     0x4a048090002L,     0x4a048090008L,     0x4a054090002L,
		     0x4a054090003L,     0x4a0540a8002L,     0x4a0540a8003L,     0x4a0540a8007L,
		     0x4a0540a8008L,     0x4a064090002L,     0x4a064090008L,     0x4a0640a8002L,
		     0x4a0640c8002L,     0x52054090002L,     0x520540a8002L,     0x520640a8002L,
		     0x520640c8002L,    0x404808090002L,    0x404808090003L,    0x4048080a8002L,
		    0x4048080c8002L,    0x404809010002L,    0x404809010003L,    0x40a808090002L,
		    0x40a808090003L,    0x40a8080a8002L,    0x40a8080c8002L,    0x40a8080c8003L,
		    0x40a809010002L,    0x40a809010003L,    0x40a814090002L,    0x40a814090003L,
		    0x40a8140a8002L,    0x40a8140c8002L,    0x40a815010002L,    0x40a815028002L,
		    0x412808090002L,    0x4128080a8002L,    0x4128080c8002L,    0x412809010002L,
		    0x412814090002L,    0x4128140a8002L,    0x4128140c8002L,    0x412815010002L,
		    0x424048090002L,    0x424048090003L,    0x424054090002L,    0x424054090003L,
		    0x4240540a8002L,    0x424064090002L,    0x424064090003L,    0x4240640a8002L,
		    0x4240640c8002L,    0x4240640c8003L,    0x424088090002L,    0x4240880a8002L,
		    0x4240880c8002L,    0x4240880c8008L,    0x424094090002L,    0x424094090008L,
		    0x4240940a8002L,    0x4240940a8007L,    0x4240940c8002L,    0x42409412800eL,
		    0x424809010002L,    0x424815010002L,    0x424815010003L,    0x424815028002L,
		    0x424815028003L,    0x424848090002L,    0x424848090003L,    0x4248480a8002L,
		    0x4248480a8003L,    0x424848110002L,    0x424848128002L,    0x424849090002L,
		    0x42a048090002L,    0x42a054090002L,    0x42a088090002L,    0x42a094090002L,
		    0x42a0940a8002L,  0x12024048090002L,  0x15024048090002L,  0x15024048090003L,
		  0x1502a048090002L,  0x1502a054090002L,  0x1502a0540a8002L,  0x19024048090002L,
		  0x1902a048090002L,  0x1902a054090002L,  0x1902a0540a8002L,  0x19032048090002L,
		  0x19032054090002L,  0x190320540a8002L,  0x19032064090002L,  0x190320640a8002L,
		  0x190320640c8002L,  0x25024048090007L,  0x2502a048090007L,  0x2502a054090007L,
		  0x2502a0540a8007L,  0x25032054090007L,  0x250320540a8007L,  0x250320640a8007L,
		  0x250320640c8007L,  0x2504a048090008L,  0x2504a054090008L,  0x2504a0540a8008L,
		  0x2504a064090008L,  0x2504a0640a8008L,  0x2902a054090007L, 0x212024048090002L,
		 0x21202a048090002L, 0x21202a054090002L, 0x21202a0540a8002L, 0x212032048090002L,
		 0x2120320640c8002L, 0x21204a0540a8007L, 0x21204a094090008L, 0x21204a0940a8008L,
		 0x21204a0940c8008L, 0x212424048090002L, 0x212424054090002L, 0x2124240540a8002L,
		 0x212424848090002L, 0x2124248480a8002L, 0x2124248480c8002L };

	private static final float[] cIncrement = {
		   -0.19f,    1.27f,  -0.701f,   2.691f,  -0.227f,    0.03f,   0.106f,  -0.476f,
		  -0.447f,  -0.191f,  -0.333f,   0.086f,   0.247f,  -0.062f,   0.016f,   0.387f,
		   0.235f,  -0.432f,  -0.903f,    0.39f,   0.581f,   4.524f,  -0.635f,   -0.28f,
		    0.77f,   -0.05f,   1.087f,   0.192f,   0.196f,   -0.52f,   0.542f,   0.363f,
		   0.792f,   0.592f,   0.964f,   0.354f,  -0.685f,  -0.315f,  -0.413f,  -0.595f,
		    0.22f,  -1.432f,  -2.254f,    0.44f,   -0.27f,  -0.133f,  -0.269f,   0.267f,
		   0.572f,  -0.568f,   0.174f,  -0.185f,  -0.235f,   -0.17f,  -0.181f,  -0.342f,
		  -0.348f,  -0.437f,  -0.804f,  -0.412f,  -0.215f,  -0.625f,  -0.831f,   0.497f,
		  -0.431f,  -1.331f,   0.507f,  -0.632f,  -0.599f,  -0.156f,  -0.353f,  -0.164f,
		  -0.441f,  -0.497f,   -1.06f,   0.115f,  -0.225f,  -0.154f,  -0.031f,  -1.574f,
		  -1.093f,  -0.738f,   -0.45f,  -0.556f,  -0.181f,   2.384f,    1.75f,  -1.666f,
		  -1.066f,   1.327f,   0.803f,  -1.505f,  -2.537f,   -0.17f,   0.149f,   0.521f,
		   2.905f,  -0.252f,    0.86f,   0.361f,   0.403f,   0.005f,   1.146f,   0.936f,
		    -0.3f,   0.209f,  -0.583f,  -0.024f,   -0.01f,   0.611f,   0.486f,   0.862f,
		  -0.035f,  -0.596f,   1.161f,   1.647f,   0.844f,   0.125f,   0.142f,  -0.171f,
		   0.442f,   0.088f,   3.066f,   1.652f,  -0.203f,  -0.018f,   0.023f,   0.073f,
		   0.254f,   0.554f,   0.595f,  -0.406f,  -0.637f,  -0.174f,  -0.101f,  -0.543f,
		  -2.406f,  -3.292f,  -0.053f,  -0.193f,    1.85f,  -1.261f,  -0.656f,   -0.73f,
		  -0.938f,   0.128f,   1.154f,   0.242f,  -0.529f,  -0.278f,  -0.802f,   0.912f,
		  -1.381f,   0.463f,   1.074f,  -0.628f,  -0.962f,  -1.832f,  -1.499f,  -2.116f,
		  -2.207f,  -1.317f,   2.501f,  -0.849f,  -0.602f,  -0.622f,   2.122f,  -2.226f,
		   0.913f,  -2.259f,   -1.25f,   1.394f,  -1.402f,   2.073f,  -2.957f,   0.291f,
		  -3.476f,  -2.727f,  -3.132f,   -2.12f,  -0.725f,  -0.297f,   0.083f,   0.347f,
		  -1.425f,   -1.66f,  -1.282f,  -1.265f,   0.719f,   0.138f,   1.302f,   0.859f,
		   1.359f,   0.659f,   -0.94f,     0.9f,   0.319f,  -2.571f,   1.109f,   0.972f,
		   1.653f,   2.602f,   0.729f,   1.066f,   1.067f,  -0.311f,   0.031f,  -0.203f,
		  -0.681f,  -1.258f,    1.07f,  -3.096f,  -0.228f,   1.933f,   0.119f,   2.108f,
		   0.113f,   1.628f,   1.308f,   3.336f,   0.754f,  -0.465f,  -0.397f,   0.077f,
		  -0.479f,  -0.153f,   0.141f,   2.135f,   0.234f,   0.461f,    0.67f,  -0.361f,
		  -1.039f,  -0.483f,   0.137f,  -0.768f,  -0.511f,   3.424f,  -0.855f,  -0.585f,
		  -1.567f,   3.343f,    1.84f,   0.389f,   1.122f,    1.63f,   1.335f,   0.366f,
		  -0.557f,   0.432f,   0.204f,   0.882f,   0.466f,  -0.458f,   0.404f,   0.657f,
		   1.115f,   1.976f,   1.786f,  -0.036f,   -1.05f,   1.045f,   0.044f,   1.033f,
		   -1.08f,   2.539f,   2.235f,    2.29f,   3.121f,   3.932f,    2.75f };


	public SolubilityPredictor() {
		}


	public float assessSolubility(StereoMolecule mol) {
		float logS = -0.530f;

		mol.normalizeAmbiguousBonds();
		mol.ensureHelperArrays(Molecule.cHelperRings);

		for (int atom=0; atom<mol.getAtoms(); atom++) {
			long type = -1;
			try {
				type = AtomTypeCalculator.getAtomType(mol, atom, AtomTypeCalculator.cPropertiesForSolubility);
				}
			catch (Exception e) {}

			for (int i =0; i<cIncrement.length; i++) {
				if (cAtomType[i] == type) {
					logS += cIncrement[i];
					break;
					}
				}
			}

		return logS;
		}


	public ParameterizedStringList getDetail(StereoMolecule mol) {
		ParameterizedStringList detail = new ParameterizedStringList();
		detail.add("Solubility values are estimated applying an atom-type based increment system.",
							ParameterizedStringList.cStringTypeText);
		detail.add("Atom-types are 64-bit numbers describing atoms and their near surrounding.",
							ParameterizedStringList.cStringTypeText);
		detail.add("Recognized atom types and their contributions are:",
							ParameterizedStringList.cStringTypeText);
		detail.add("Base value = -0.530",
							ParameterizedStringList.cStringTypeText);

		int count[] = new int[cIncrement.length];

		if (mol != null) {
			mol.normalizeAmbiguousBonds();
			mol.ensureHelperArrays(Molecule.cHelperRings);

			for (int atom=0; atom<mol.getAtoms(); atom++) {
				long type = -1;
				try {
					type = AtomTypeCalculator.getAtomType(mol, atom, AtomTypeCalculator.cPropertiesForSolubility);
					}
				catch (Exception e) {}

				for (int i =0; i<cIncrement.length; i++) {
					if (cAtomType[i] == type) {
						count[i]++;
						break;
						}
					}
				}
			}
		NumberFormat formatter = new DecimalFormat("#0.000");

		for (int i=0; i<cIncrement.length; i++)
			if (count[i] != 0)
				detail.add(""+count[i]+" * "+
									formatter.format(cIncrement[i])
									+"   AtomType: 0x"
								  +Long.toHexString(cAtomType[i]),
								  ParameterizedStringList.cStringTypeText);

		return detail;
		}
	}
