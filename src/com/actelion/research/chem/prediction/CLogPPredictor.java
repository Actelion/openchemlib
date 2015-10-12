/*
* Copyright (c) 1997 - 2015
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

import java.util.TreeMap;

import com.actelion.research.chem.AtomTypeCalculator;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.SortedList;

public class CLogPPredictor {
	private static final int ATOM_TYPE_MODE = AtomTypeCalculator.cPropertiesForCLogPCharges;

	protected static final long[] ATOM_TYPE = {
	           0x40002L,           0x40004L,           0x40802L,           0x40803L,
	           0x40804L,           0x40805L,           0x40808L,           0x40809L,
	           0x4080dL,           0x4080eL,           0x41402L,           0x41403L,
	           0x41404L,           0x41409L,           0x4140eL,           0x42402L,
	           0x42403L,           0x42404L,           0x42405L,           0x4240eL,
	           0x60802L,           0x60803L,           0x60804L,           0x60805L,
	           0x60808L,           0x60809L,           0x6080dL,           0x6080eL,
	           0x61402L,           0x61404L,           0x61409L,           0x80004L,
	           0x80802L,           0x80803L,           0x80804L,           0x80808L,
	           0x81404L,           0x82404L,           0x82408L,           0xc0802L,
	           0xc0803L,         0x8220842L,         0x8220843L,         0x8220844L,
	         0x8220848L,         0x8520842L,         0x8520843L,         0x8520844L,
	         0x8520848L,         0x8521442L,         0x8521443L,         0x8521444L,
	         0x8521448L,         0x8920842L,         0x8920843L,         0x8921442L,
	         0x8921443L,         0x8921444L,        0x10240002L,        0x10240008L,
	        0x10240048L,        0x10240802L,        0x10240803L,        0x10240804L,
	        0x10240808L,        0x1024080cL,        0x10240842L,        0x10240843L,
	        0x10240844L,        0x10240848L,        0x10540802L,        0x10540803L,
	        0x10540804L,        0x10540808L,        0x10540842L,        0x10540843L,
	        0x10540844L,        0x10540848L,        0x10541402L,        0x10541408L,
	        0x10541442L,        0x10940004L,        0x10940802L,        0x10940803L,
	        0x10940804L,        0x10940808L,        0x10940842L,        0x10940843L,
	        0x10940844L,        0x10940848L,        0x10941402L,        0x10941403L,
	        0x10941404L,        0x10941442L,        0x10942402L,        0x10942408L,
	        0x1094240aL,        0x10942442L,        0x18240802L,        0x18240803L,
	        0x18240804L,        0x18240808L,        0x1824080cL,        0x18240842L,
	        0x18240843L,        0x18240844L,        0x18240848L,        0x1824084cL,
	        0x18241402L,        0x18241403L,        0x18241404L,        0x18241442L,
	        0x18241448L,        0x18242402L,        0x18242403L,        0x18242404L,
	        0x18242408L,        0x18242442L,        0x18242443L,        0x18242444L,
	        0x18242448L,        0x18260802L,        0x18260803L,        0x18260804L,
	        0x18260808L,        0x18260842L,        0x18260843L,        0x18260844L,
	        0x18260848L,        0x18540802L,        0x18540804L,        0x18540842L,
	        0x18541402L,        0x18560802L,        0x18560803L,        0x20240802L,
	        0x20240803L,        0x20240842L,        0x20240843L,        0x20241402L,
	        0x20241403L,        0x20241442L,        0x20241443L,        0x20241444L,
	        0x20242402L,        0x20242403L,        0x20242442L,        0x20242443L,
	        0x20260802L,        0x20260803L,        0x20260842L,        0x20260843L,
	        0x20261402L,        0x20540802L,        0x20540803L,        0x20540842L,
	        0x20540843L,        0x20540844L,        0x20540848L,        0x20541402L,
	        0x20541403L,        0x20541442L,        0x20541443L,        0x20560802L,
	        0x20560803L,        0x20560843L,        0x20561402L,        0x20981402L,
	        0x30240802L,        0x30260802L,        0x30540802L,        0x30541402L,
	        0x30542402L,        0x30560802L,      0x2088220842L,      0x2088220843L,
	      0x2148220842L,      0x2148220843L,      0x2148520842L,      0x2148520848L,
	      0x2148521442L,      0x2248220842L,      0x2248520842L,      0x2248920842L,
	      0x4088220842L,      0x4088220843L,      0x4088520842L,      0x4088520843L,
	      0x4088521442L,      0x4088521443L,      0x4088920842L,      0x4088921442L,
	      0x4090240043L,      0x4090240802L,      0x4090240803L,      0x4090240842L,
	      0x4090240843L,      0x4148220842L,      0x4148220843L,      0x4148520842L,
	      0x4148520843L,      0x4148521442L,      0x4148920842L,      0x4148921442L,
	      0x4150240802L,      0x4150240803L,      0x4150240842L,      0x4150240843L,
	      0x4150540802L,      0x4150540803L,      0x4150540842L,      0x4150541402L,
	      0x4150541442L,      0x4248220842L,      0x4248520842L,      0x4248521442L,
	      0x4248920842L,      0x4248921442L,      0x4250240802L,      0x4250240803L,
	      0x4250240842L,      0x4250240843L,      0x4250540802L,      0x4250540803L,
	      0x4250540842L,      0x4250541402L,      0x4250940802L,      0x4250942402L,
	      0x6088220842L,      0x6088220843L,      0x6088520842L,      0x6088520843L,
	      0x6088521442L,      0x6088521443L,      0x6088920842L,      0x6088921442L,
	      0x6090240802L,      0x6090240803L,      0x6090240842L,      0x6090240843L,
	      0x6090540802L,      0x6090540803L,      0x6090540842L,      0x6090540843L,
	      0x6090541402L,      0x6090541442L,      0x6090940802L,      0x6090940803L,
	      0x6090940842L,      0x6090940843L,      0x6090941403L,      0x6090941442L,
	      0x6090942403L,      0x6098240802L,      0x6098240803L,      0x6098240842L,
	      0x6098240843L,      0x6098241402L,      0x6098241442L,      0x6098260802L,
	      0x6148220842L,      0x6148520842L,      0x6148521442L,      0x6150240802L,
	      0x6150240842L,      0x6150540802L,      0x6150540842L,      0x6150940842L,
	      0x6158240802L,      0x6158260802L,      0x8090240802L,      0x8090240842L,
	      0x8090240843L,      0x8090540802L,      0x8090540842L,      0x8090540843L,
	      0x8090541402L,      0x8090541442L,      0x8090940802L,      0x8090940842L,
	      0x8090941442L,      0x8090942402L,      0x8090942442L,      0x8098240802L,
	      0x8098240842L,      0x8098241402L,      0x8098241442L,      0x8098242402L,
	      0x8098242408L,      0x8098242442L,      0x8098260802L,      0x8098260842L,
	      0x8098560842L,      0x8150240802L,      0x8150240808L,      0x8150240842L,
	      0x8150240843L,      0x8150240848L,      0x8150540802L,      0x8150540842L,
	      0x8150541402L,      0x8150541407L,      0x8150541442L,      0x8150541448L,
	      0x8150940802L,      0x8150940842L,      0x8150941402L,      0x8150941442L,
	      0x8150942442L,      0x8158240802L,      0x8158240808L,      0x8158240842L,
	      0x8158240843L,      0x8158240848L,      0x8158241402L,      0x8158241408L,
	      0x8158241442L,      0x8158242402L,      0x8158242442L,      0x8158260802L,
	      0x8158260808L,      0x8158260842L,      0x8158540802L,      0x8158541402L,
	      0x8158560802L,      0x8250240842L,      0x8250540802L,      0x8250540842L,
	      0x8250541402L,      0x8250541442L,      0x8250940842L,      0x8250941402L,
	      0x8250941442L,      0x8258241402L,      0x8258241442L,      0x8258541402L,
	   0x1024090240043L,   0x1024090240801L,   0x1024090240802L,   0x1024090240842L,
	   0x1054090240802L,   0x1054090240842L,   0x1054150240802L,   0x1054150240842L,
	   0x1054150540802L,   0x1054150540842L,   0x1054150541402L,   0x1094088220843L,
	   0x1094090240802L,   0x1094090240803L,   0x1094090240842L,   0x1094090240843L,
	   0x1094150240802L,   0x1094150540802L,   0x1094150541402L,   0x1094250240802L,
	   0x1094250240842L,   0x1094250540802L,   0x1094250940802L,   0x1094250941402L,
	   0x1094250942402L,   0x1824090240802L,   0x1824090240842L,   0x1824150240802L,
	   0x1824150240842L,   0x1824150540802L,   0x1824150540842L,   0x1824150541402L,
	   0x1824250942402L,   0x1826090240802L,   0x1826090240842L,   0x1826090540802L,
	   0x1826090540842L,   0x1826098241402L,   0x1826098241442L,   0x1854090240802L,
	   0x1854090240842L,   0x2024250240803L,   0x2054150240807L,   0x2054150540807L,
	   0x2054150541407L,   0x2054150541447L,   0x2054250540807L,   0x2054250541407L,
	   0x2054250941407L,   0x205425094144aL,   0x2056090541407L,   0x2058150240808L,
	   0x2058150240848L,   0x2058150540808L,   0x2058150541408L,   0x2058150541448L,
	   0x2058158240808L,   0x2058158240848L,   0x2058158241408L,   0x2058158241448L,
	   0x2058158260808L,   0x2058158260848L,   0x2094150540807L,   0x2094150541407L,
	   0x2094250541407L,   0x2096090541407L,   0x4000000040803L,   0x4000000040804L,
	   0x4000000040808L,   0x4000000041404L,   0x4000000042404L,   0x4000000060804L,
	   0x4000000061404L,   0x4000000081403L,   0x4000008220843L,   0x4000010240803L,
	   0x4000010540803L,   0x4000020580803L,   0x4000020581403L,   0x4000030241403L,
	   0x4004088220843L,   0x4004088520843L,   0x4004090240803L,   0x4004090240843L,
	   0x4004148220843L,   0x4004148520843L,   0x4006088520843L,   0x4008090240803L,
	   0x4008090540803L,   0x4008090540843L,   0x4008098241403L,   0x4008150540803L,
	   0x4008150541403L,   0x4008158241403L,   0x5024090240803L,   0x5024090240843L,
	   0x5054090240803L,   0x5054090240843L,   0x5824090240803L,   0x8000000040803L,
	   0x8000008220843L,   0x8000010240803L,   0x8000010240843L,   0x8000020242443L,
	   0x8004090240803L,   0x8004090240843L };

	protected static final float[] INCREMENT = {
	  0.6967f,  0.0000f,  0.4886f, -0.4727f, -0.0749f,  0.6262f,  0.2735f,  0.5700f,
	  0.7010f,  0.9534f, -0.2809f, -0.8260f, -0.1785f, -1.6203f, -1.0960f,  0.1395f,
	 -0.2975f, -1.2908f,  1.0162f, -1.3825f,  0.5111f, -0.4357f, -0.1041f,  0.3424f,
	 -0.0615f,  0.6035f,  0.7227f,  0.4346f, -0.3310f, -0.4980f, -1.4915f,  0.3317f,
	  0.4292f, -0.5824f, -0.1834f,  0.1306f, -0.5015f, -0.5258f,  0.4244f, -0.1610f,
	 -0.2778f,  0.2766f,  0.3593f,  0.7715f,  0.3150f, -0.2652f, -0.0965f,  0.4202f,
	  0.1871f, -0.3684f, -0.0778f,  0.8943f,  0.3694f,  0.2879f,  0.4489f, -0.2601f,
	  0.4771f,  0.1923f,  0.4597f,  0.3384f,  0.6633f,  0.4544f,  0.1597f,  0.6339f,
	  0.3504f,  0.0449f,  0.3420f,  0.2611f,  0.4046f,  0.5219f, -0.3632f, -0.4108f,
	  0.3057f, -0.1456f, -0.2713f, -0.5193f,  0.4526f,  0.5539f, -0.7070f, -0.4881f,
	 -0.4100f,  0.0000f,  0.1479f,  0.3448f,  0.4298f,  0.5579f, -0.1265f, -0.0425f,
	  0.0767f,  0.6635f, -0.3812f, -0.8368f,  1.0287f, -0.1021f,  0.3587f, -0.5945f,
	  0.1692f, -0.1218f,  0.4381f,  0.1695f,  0.4525f,  0.3352f,  0.1583f,  0.4036f,
	 -0.0480f,  0.5023f, -0.2649f,  0.7691f, -0.3552f,  1.0301f, -0.1141f, -0.5932f,
	  0.1749f,  0.1313f, -0.1804f,  0.3994f,  0.2291f,  0.3169f,  0.3599f, -0.0039f,
	 -0.2956f,  0.4907f,  0.3540f,  0.2192f,  0.1565f,  0.6935f,  0.3618f,  0.6735f,
	  0.5778f, -0.5636f,  0.5569f,  0.3038f, -0.3276f, -0.4659f,  0.3818f,  0.3283f,
	  0.2239f,  0.2043f,  0.0590f, -0.4835f,  0.6165f, -0.4011f,  0.5578f, -0.2164f,
	 -0.0175f,  0.2981f,  0.1100f,  0.2715f,  0.4409f, -0.1609f,  0.3775f, -0.1346f,
	 -0.6992f, -0.4670f,  0.1566f,  0.0468f, -0.1321f,  1.3686f,  0.0000f, -0.4116f,
	  1.0186f, -0.3935f,  0.5223f,  0.2839f,  0.5129f,  0.1266f,  0.0103f,  1.5193f,
	  0.2705f,  0.4294f,  0.0120f, -0.3397f,  0.1483f,  0.2806f,  0.3206f,  0.5662f,
	 -0.0987f, -0.1005f, -0.3576f,  0.0961f, -0.6401f,  0.1921f, -0.1533f, -0.4170f,
	  0.1094f,  0.8231f, -0.3784f,  0.4032f, -0.6461f,  0.8035f,  0.2029f, -0.3745f,
	  0.3317f,  0.1841f,  0.7071f,  0.1227f,  0.7949f,  0.0350f,  0.3818f, -0.1554f,
	  0.3785f, -0.2405f,  0.2359f,  0.3463f, -0.4925f, -0.0929f, -0.4352f, -0.2207f,
	 -0.9960f, -0.7238f, -0.5469f, -1.2939f, -0.0136f,  0.2791f, -0.1653f, -0.1238f,
	  0.4951f,  0.2899f,  0.0657f,  0.7189f,  0.0570f,  0.6619f, -0.6381f, -0.8073f,
	  0.2355f,  0.3048f, -0.0199f, -0.0752f,  0.2764f,  0.8011f, -0.1744f,  0.1581f,
	 -0.3848f,  0.5993f,  0.5268f, -0.0417f,  0.3770f,  0.6998f,  0.5940f,  0.5912f,
	 -0.5571f,  0.0238f, -0.2475f,  0.0307f, -0.3875f, -0.7437f,  0.5144f,  0.0057f,
	  0.7655f,  0.1720f, -2.5624f, -0.3066f,  0.3647f,  0.4733f, -0.3401f, -0.1450f,
	  0.7088f, -0.1318f,  0.0426f, -0.1203f, -0.3624f,  0.5358f, -0.3701f, -0.5648f,
	 -0.1972f, -0.8769f, -0.3675f, -0.2004f,  0.1336f, -0.1699f,  0.4461f,  0.1559f,
	  1.1168f,  0.2365f, -0.2206f,  0.4480f, -0.4053f, -0.1361f,  0.2199f,  0.0536f,
	 -0.0210f,  0.6985f,  0.9643f,  0.1727f, -0.0329f, -0.1893f,  0.0702f,  0.1496f,
	 -1.3825f,  0.4146f, -0.5028f,  0.3832f,  0.9545f, -0.4152f, -1.0369f, -0.1830f,
	  0.5883f, -0.2918f, -0.5294f, -0.6541f, -0.1906f, -0.8484f, -0.3457f,  0.9541f,
	 -0.7924f, -0.6020f,  0.0800f, -0.2596f,  0.8382f, -0.2668f, -0.1106f,  0.0362f,
	 -0.3189f, -0.7278f, -0.0894f, -0.2277f, -0.2394f, -0.2962f,  0.7776f, -0.0118f,
	 -0.4358f,  0.3749f, -0.6070f, -0.1857f,  0.1139f, -0.4416f, -0.3704f, -0.7487f,
	 -0.1079f, -0.2992f, -0.3277f,  0.0251f, -0.9188f,  0.2945f, -0.2234f,  0.3468f,
	  0.3317f,  0.2891f,  0.2613f, -0.0344f, -0.6005f, -0.6258f, -0.5434f, -0.7712f,
	 -0.9057f, -0.1668f, -0.9905f, -1.4915f, -0.0372f, -1.1638f,  0.1262f, -0.5248f,
	 -0.1538f, -0.3682f,  0.3249f,  0.0650f,  0.0511f, -0.4607f,  0.2231f,  0.2822f,
	  0.1397f,  0.2833f, -0.1226f, -0.4592f, -0.3435f, -0.6654f, -0.5056f, -0.8631f,
	  0.1536f, -0.4051f,  0.0891f, -0.6972f, -0.4699f, -0.6774f, -0.0622f, -0.9300f,
	  0.1337f, -0.4938f,  0.3948f, -0.4075f, -0.6411f, -0.0091f, -0.1333f, -0.5192f,
	 -0.1661f,  0.3317f, -0.6427f, -0.0707f,  0.4806f,  0.3828f,  0.2229f,  0.6160f,
	 -0.0884f, -0.0471f,  0.1106f,  0.3821f,  0.0922f,  0.0806f,  0.3371f,  0.1884f,
	  0.1381f, -0.2392f, -3.6435f, -2.1509f,  0.4975f, -0.3620f, -2.5384f, -1.6821f,
	 -0.3651f,  0.6262f, -1.6185f, -1.0650f,  0.8374f,  0.3685f,  0.2577f,  0.3791f,
	 -3.2333f, -1.7948f, -0.6592f, -1.3148f, -0.7380f,  0.0534f, -1.7552f, -1.8039f,
	 -1.1340f, -0.5653f, -1.2454f,  0.9473f,  1.4231f,  1.0112f, -1.9498f, -2.0249f,
	 -1.2349f,  0.3280f, -3.9189f, -2.1995f,  0.1889f, -1.2314f, -1.8023f, -0.2995f,
	 -0.4067f, -0.1316f };
	
	public static final float cCLogPUnknown = -999f;
	private static SortedList<Long> sSortedTypeList;

	public CLogPPredictor() {
		if (sSortedTypeList == null) {
			synchronized(CLogPPredictor.class) {
				if (sSortedTypeList == null) {
					sSortedTypeList = new SortedList<Long>();
					for (Long l:ATOM_TYPE)
						sSortedTypeList.add(l);
					}
				}
			}
		}

	/**
	 * Ambiguous bonds are normalized. 
	 * @param mol
	 * @return
	 */
	public float assessCLogP(StereoMolecule mol) {
		float cLogP = 0.0f;

		mol.normalizeAmbiguousBonds();
		mol.ensureHelperArrays(Molecule.cHelperRings);

		for (int atom=0; atom<mol.getAtoms(); atom++) {
			try {
				int index = sSortedTypeList.getIndex(AtomTypeCalculator.getAtomType(mol, atom, ATOM_TYPE_MODE));
				if (index != -1)
					cLogP += INCREMENT[index];
				}
			catch (Exception e) {}	// unsupported atom type exceptions are tolerable
			}
		
		return cLogP;
		}

	public ParameterizedStringList getDetail(StereoMolecule mol) {
		ParameterizedStringList detail = new ParameterizedStringList();
		detail.add("cLogP Values are estimated applying an atom-type based increment system.",
							ParameterizedStringList.cStringTypeText);
		detail.add("Atom-types are 64-bit numbers describing atoms and their near surrounding.",
							ParameterizedStringList.cStringTypeText);
		detail.add("Recognized atom types and their contributions are:",ParameterizedStringList.cStringTypeText);

		mol.normalizeAmbiguousBonds();
		mol.ensureHelperArrays(Molecule.cHelperRings);

		if (mol != null) {
			int errorCount = 0;
			TreeMap<Long,Integer> countMap = new TreeMap<Long,Integer>();
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				try {
					long atomType = AtomTypeCalculator.getAtomType(mol, atom, ATOM_TYPE_MODE);
					Integer typeCount = countMap.get(new Long(atomType));
					if (typeCount == null)
						countMap.put(new Long(atomType), new Integer(1));
					else
						countMap.put(new Long(atomType), new Integer(typeCount.intValue()+1));
					}
				catch (Exception e) {
					errorCount++;
					}
				}

			if (errorCount != 0)
				detail.add("Warning: "+errorCount + " atom type(s) could not be determined.", ParameterizedStringList.cStringTypeText);

			for (Long type:countMap.keySet()) {
				if (sSortedTypeList.contains(type))
					detail.add(countMap.get(type) + " * "+ INCREMENT[sSortedTypeList.getIndex(type)] + " AtomType: 0x" + Long.toHexString(type),ParameterizedStringList.cStringTypeText);
				else
					detail.add("Warning: For atom type 0x"+Long.toHexString(type)+" ("+countMap.get(type)+" times found) is no increment available.", ParameterizedStringList.cStringTypeText);
				}
			}
		
		return detail;
		}
	}