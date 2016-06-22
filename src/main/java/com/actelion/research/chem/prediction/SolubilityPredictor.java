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
import com.actelion.research.chem.StereoMolecule;

import java.text.DecimalFormat;
import java.text.NumberFormat;


public class SolubilityPredictor {
	public static final float cSolubilityUnknown = -999;

	private static final float[] cIncrement = {
	 -0.190f,  1.270f, -0.701f,  2.691f, -0.227f,  0.030f,  0.106f, -0.476f,
	 -0.447f, -0.191f, -0.333f,  0.086f,  0.247f, -0.062f,  0.016f,  0.387f,
	  0.235f, -0.432f, -0.903f,  0.390f,  0.581f,  4.524f, -0.635f,  0.792f,
	  0.592f,  0.964f,  0.354f, -0.685f, -0.315f, -0.413f, -0.595f,  0.220f,
	 -0.280f,  0.770f, -0.050f,  1.087f,  0.192f,  0.196f, -0.520f,  0.542f,
	  0.363f, -0.181f,  2.384f,  1.750f, -1.666f, -1.066f,  1.327f,  0.803f,
	 -1.505f, -2.537f, -0.170f,  0.149f,  0.521f,  2.905f, -0.252f, -1.432f,
	 -2.254f,  0.440f, -0.270f, -0.133f, -0.269f,  0.267f,  0.572f, -0.568f,
	  0.174f, -0.185f, -0.235f, -0.170f, -0.181f, -0.342f, -0.348f, -0.437f,
	 -0.804f, -0.412f, -0.215f, -0.625f, -0.831f,  0.497f, -0.431f, -1.331f,
	  0.507f, -0.632f, -0.599f,  0.860f,  0.361f,  0.403f,  0.005f,  1.146f,
	  0.936f, -0.300f,  0.209f, -0.583f, -0.024f, -0.010f,  1.647f,  0.844f,
	  0.125f,  0.142f, -0.171f,  0.442f,  0.088f,  3.066f,  1.652f, -0.156f,
	 -0.353f, -0.164f, -0.441f, -0.497f, -1.060f,  0.611f,  0.486f,  0.115f,
	 -0.225f, -0.154f, -0.031f,  0.862f, -0.035f, -0.596f, -1.574f, -1.093f,
	  1.161f, -0.738f, -0.450f, -0.556f, -0.622f,  2.122f, -1.402f,  2.073f,
	 -3.132f, -2.120f,  0.347f, -1.265f, -1.317f,  2.501f, -2.226f,  0.913f,
	 -2.957f,  0.291f, -0.725f, -1.425f, -0.203f, -0.018f, -0.849f, -2.259f,
	 -3.476f, -0.297f, -1.660f,  0.023f,  0.073f,  0.254f,  0.554f,  0.595f,
	 -0.602f, -1.250f,  1.394f, -2.727f,  0.083f, -1.282f, -0.406f, -0.637f,
	 -0.174f, -0.101f, -0.543f, -2.406f, -3.292f, -0.681f, -1.258f,  1.070f,
	 -3.096f, -0.228f,  0.719f,  0.138f,  1.302f,  0.859f,  1.359f,  0.659f,
	 -0.940f,  0.900f,  0.319f, -2.571f,  1.933f,  0.119f,  2.108f,  0.113f,
	  3.336f,  0.754f, -0.465f, -0.053f, -0.193f,  1.850f, -1.261f, -0.656f,
	 -0.730f, -0.938f,  1.109f,  0.972f,  1.653f,  2.602f,  1.628f, -0.397f,
	  0.128f,  1.154f,  0.242f, -0.529f, -0.278f, -0.802f,  0.912f, -1.381f,
	  0.463f,  1.074f, -0.628f, -0.962f,  0.729f,  1.066f,  1.067f, -0.311f,
	  0.031f,  1.308f,  0.077f, -0.479f, -0.203f, -1.832f, -1.499f, -2.116f,
	 -2.207f, -0.153f,  0.141f,  2.135f,  0.234f,  0.461f,  0.670f, -0.361f,
	 -1.039f, -0.483f,  0.137f, -0.768f, -0.511f,  3.424f, -0.855f, -0.585f,
	 -1.567f,  0.657f,  1.115f,  1.976f,  1.786f, -0.036f, -1.050f,  2.539f,
	  2.235f,  2.290f,  3.121f,  3.932f,  2.750f,  3.343f,  1.840f,  0.389f,
	  1.122f,  1.630f,  1.335f,  0.366f, -0.557f,  1.045f,  0.432f,  0.204f,
	  0.882f,  0.466f, -0.458f,  0.044f,  1.033f, -1.080f,  0.404f };

	private static final long[] cAtomType = {
			262146L,			262148L,			262153L,			262157L,
			264194L,			264195L,			264196L,			264197L,
			264200L,			264201L,			264205L,			264206L,
			267266L,			267267L,			267268L,			267273L,
			267277L,			271362L,			271363L,			271364L,
			271365L,			271368L,			271369L,			395266L,
			395267L,			395268L,			395269L,			395272L,
			395273L,			395277L,			395278L,			398338L,
			526338L,			526339L,			526340L,			526344L,
			529412L,			533508L,			533512L,			788482L,
			788483L,		 136448002L,		 136448003L,		 136448004L,
		 136448008L,		 139593730L,		 139593731L,		 139593732L,
		 139593736L,		 139596802L,		 139596803L,		 139596804L,
		 143788034L,		 143788035L,		 143791106L,		 268697604L,
		 270794754L,		 270794756L,		 270796802L,		 270796803L,
		 270796804L,		 270796808L,		 270796812L,		 273940482L,
		 273942530L,		 273942531L,		 273942532L,		 273942536L,
		 273945602L,		 273945608L,		 273945612L,		 278136834L,
		 278136835L,		 278136836L,		 278136840L,		 278136844L,
		 278139906L,		 278139907L,		 278139908L,		 278144002L,
		 278144003L,		 278144004L,		 278144008L,		 405014530L,
		 405014531L,		 405014532L,		 405014536L,		 405017602L,
		 405017603L,		 405017604L,		 405021698L,		 405021699L,
		 405021700L,		 405021704L,		 405145602L,		 405145603L,
		 405145604L,		 405145608L,		 408158210L,		 408160258L,
		 408163330L,		 408167426L,		 408291330L,		 539232258L,
		 539232259L,		 539235330L,		 539235331L,		 539239426L,
		 539239427L,		 539363330L,		 539363331L,		 542377986L,
		 542377987L,		 542381058L,		 542381059L,		 542509058L,
		 542509059L,		 542509070L,		 546837506L,		 807667714L,
		 807798786L,		 810813442L,		 810816514L,		 810820610L,
	  139722885122L,	  139722885123L,	  142944110594L,	  142944110595L,
	  142947256322L,	  142947259394L,	  147239077890L,	  147242223618L,
	  277161838594L,	  277161838595L,	  277164984322L,	  277164984323L,
	  277164987394L,	  277164987395L,	  277169178626L,	  277169181698L,
	  277296187394L,	  277296187395L,	  280383064066L,	  280386209794L,
	  280386212866L,	  280390404098L,	  280390407170L,	  280517412866L,
	  280517412867L,	  280520558594L,	  280520558595L,	  280520561666L,
	  284678031362L,	  284681177090L,	  284681177091L,	  284681180162L,
	  284685371394L,	  284685374466L,	  284812380162L,	  284812380163L,
	  284815525890L,	  284815528962L,	  284819720194L,	  284819727362L,
	  284819727363L,	  414600792066L,	  414603937794L,	  414603937795L,
	  414603940866L,	  414603940867L,	  414735140866L,	  414735140867L,
	  414738286594L,	  414738286595L,	  414738289666L,	  414742480898L,
	  414742480899L,	  414742483970L,	  414742488066L,	  414742488067L,
	  414869358594L,	  414869358595L,	  414869361666L,	  414869361667L,
	  414869489666L,	  417956366338L,	  417959512066L,	  552174094338L,
	  552177240066L,	  552177240067L,	  552177243138L,	  552181434370L,
	  552181437442L,	  552181441538L,	  552308312066L,	  552308315138L,
	  552308319234L,	  552308319240L,	  552308443138L,	  552311457794L,
	  555395319810L,	  555395319816L,	  555398465538L,	  555398465539L,
	  555398468610L,	  555398468611L,	  555398468615L,	  555398468616L,
	  555402659842L,	  555402659848L,	  555402662914L,	  555402667010L,
	  555529537538L,	  555529537544L,	  555529540610L,	  555529540615L,
	  555529544706L,	  555529668610L,	  555532683266L,	  555532686338L,
	  555667032078L,	  559693432834L,	  559693435906L,	  559697630210L,
	  559697634306L,   283951296153602L,   287249831036930L,   287249831036931L,
   287253052262402L,   287253055408130L,   287253055411202L,   291647877548034L,
   291651098773506L,   291651101919234L,   291651101922306L,   291655393740802L,
   291655396886530L,   291655396889602L,   291655401080834L,   291655401083906L,
   291655401088002L,   424688784508930L,   424692005734402L,   424692008880130L,
   424692008883202L,   424696300701698L,   424696308048898L,   424826223462402L,
   424826226608130L,   424826226611202L,   424826357680130L,   424826357683202L,
   424826357687298L,   568724807747591L,   568728028973063L,   568728032118791L,
   568728032121863L,   568732327086087L,   568732327089159L,   568732331283463L,
   568732331287559L,   568862249849863L,   569002906880008L,   569002910025736L,
   569002910028808L,   569002914220040L,   569002914223112L,   569003041097736L,
   569003041100808L,   569003041104904L,   573126078629895L };


	public SolubilityPredictor() {
		}


	public float assessSolubility(StereoMolecule mol) {
		float logS = -0.530f;

		for (int atom=0; atom<mol.getAtoms(); atom++) {
			long type = -1;
			try {
				type = AtomTypeCalculator.getAtomType(mol, atom,
											 AtomTypeCalculator.cPropertiesForSolubility);
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
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				long type = -1;
				try {
					type = AtomTypeCalculator.getAtomType(mol, atom,
												 AtomTypeCalculator.cPropertiesForSolubility);
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
