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
import com.actelion.research.util.SortedList;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.TreeMap;

public class CLogPPredictor {
	private static final int ATOM_TYPE_MODE = AtomTypeCalculator.cPropertiesForCLogPCharges;

	protected static final long[] ATOM_TYPE = {
		           0x80002L,           0x80004L,           0x90002L,           0x90003L,
		           0x90004L,           0x90005L,           0x90008L,           0x90009L,
		           0x9000dL,           0x9000eL,           0x92003L,           0x92004L,
		           0x92008L,           0x94003L,           0xa8002L,           0xa8003L,
		           0xa8004L,           0xa8009L,           0xa800eL,           0xaa004L,
		           0xc8002L,           0xc8003L,           0xc8004L,           0xc8005L,
		           0xc800eL,           0xca004L,          0x100004L,          0x110002L,
		          0x110003L,          0x110004L,          0x110008L,          0x128004L,
		          0x12a003L,          0x148004L,          0x148008L,          0x190002L,
		          0x190003L,         0x1090002L,         0x1090003L,         0x1090004L,
		         0x1090005L,         0x1090008L,         0x1090009L,         0x109000dL,
		         0x109000eL,         0x1092004L,         0x10a8002L,         0x10a8004L,
		         0x10a8009L,         0x10aa004L,        0x48080002L,        0x48080008L,
		        0x48080048L,        0x48090002L,        0x48090003L,        0x48090004L,
		        0x48090008L,        0x4809000cL,        0x48090042L,        0x48090043L,
		        0x48090044L,        0x48090048L,        0x48092003L,        0x48094003L,
		        0x48094043L,        0x54090002L,        0x54090003L,        0x54090004L,
		        0x54090008L,        0x54090042L,        0x54090043L,        0x54090044L,
		        0x54090048L,        0x54092003L,        0x540a8002L,        0x540a8008L,
		        0x540a8042L,        0x64080004L,        0x64090002L,        0x64090003L,
		        0x64090004L,        0x64090008L,        0x64090042L,        0x64090043L,
		        0x64090044L,        0x64090048L,        0x640a8002L,        0x640a8003L,
		        0x640a8004L,        0x640a8042L,        0x640c8002L,        0x640c8008L,
		        0x640c800aL,        0x640c8042L,        0x88090002L,        0x88090003L,
		        0x88090042L,        0x88090043L,        0x880a8002L,        0x880a8003L,
		        0x880a8042L,        0x880a8043L,        0x880a8044L,        0x880c8002L,
		        0x880c8003L,        0x880c8042L,        0x880c8043L,        0x880cc043L,
		        0x94090002L,        0x94090003L,        0x94090042L,        0x94090043L,
		        0x94090044L,        0x94090048L,        0x940a8002L,        0x940a8003L,
		        0x940a8042L,        0x940a8043L,        0x94112003L,        0x9412a003L,
		        0xa4128002L,        0xc8090002L,        0xc80aa003L,        0xd4090002L,
		        0xd40a8002L,        0xd40c8002L,       0x809010042L,       0x809010043L,
		       0x809010044L,       0x809010048L,       0x809012043L,       0x809014043L,
		       0x815010042L,       0x815010043L,       0x815010044L,       0x815010048L,
		       0x815028042L,       0x815028043L,       0x815028044L,       0x815028048L,
		       0x825010042L,       0x825010043L,       0x825028042L,       0x825028043L,
		       0x825028044L,       0x848090002L,       0x848090003L,       0x848090004L,
		       0x848090008L,       0x84809000cL,       0x848090042L,       0x848090043L,
		       0x848090044L,       0x848090048L,       0x84809004cL,       0x8480a8002L,
		       0x8480a8003L,       0x8480a8004L,       0x8480a8042L,       0x8480a8048L,
		       0x8480c8002L,       0x8480c8003L,       0x8480c8004L,       0x8480c8008L,
		       0x8480c8042L,       0x8480c8043L,       0x8480c8044L,       0x8480c8048L,
		       0x848110002L,       0x848110003L,       0x848110042L,       0x848110043L,
		       0x848128002L,       0x848128003L,       0x848128043L,       0x848190002L,
		       0x8481a8002L,       0x849090002L,       0x849090003L,       0x849090004L,
		       0x849090008L,       0x849090042L,       0x849090043L,       0x849090044L,
		       0x849090048L,       0x854090002L,       0x854090004L,       0x854090042L,
		       0x8540a8002L,       0x854110002L,       0x854128002L,       0x855090002L,
		       0x855090003L,     0x24048080043L,     0x24048090002L,     0x24048090003L,
		     0x24048090042L,     0x24048090043L,     0x24048092003L,     0x24048092043L,
		     0x24048094003L,     0x24048094043L,     0x2a048090002L,     0x2a048090003L,
		     0x2a048090042L,     0x2a048090043L,     0x2a054090002L,     0x2a054090003L,
		     0x2a054090042L,     0x2a0540a8002L,     0x2a0540a8042L,     0x32048090002L,
		     0x32048090003L,     0x32048090042L,     0x32048090043L,     0x32054090002L,
		     0x32054090003L,     0x32054090042L,     0x320540a8002L,     0x32064090002L,
		     0x320640c8002L,     0x44048090002L,     0x44048090042L,     0x44048090043L,
		     0x44048092003L,     0x44054090002L,     0x44054090042L,     0x44054090043L,
		     0x44054092003L,     0x44054092043L,     0x440540a8002L,     0x440540a8042L,
		     0x44064090002L,     0x44064090042L,     0x440640a8042L,     0x440640c8002L,
		     0x440640c8042L,     0x4a048090002L,     0x4a048090008L,     0x4a048090042L,
		     0x4a048090043L,     0x4a048090048L,     0x4a054090002L,     0x4a054090042L,
		     0x4a054092003L,     0x4a0540a8002L,     0x4a0540a8007L,     0x4a0540a8042L,
		     0x4a0540a8048L,     0x4a0540aa003L,     0x4a064090002L,     0x4a064090042L,
		     0x4a0640a8002L,     0x4a0640a8042L,     0x4a0640c8042L,     0x52048090042L,
		     0x52054090002L,     0x52054090042L,     0x520540a8002L,     0x520540a8042L,
		     0x52064090042L,     0x520640a8002L,     0x520640a8042L,    0x404808090042L,
		    0x404808090043L,    0x404808092043L,    0x4048080a8042L,    0x4048080a8043L,
		    0x4048080aa043L,    0x4048080c8042L,    0x404809010042L,    0x404809010043L,
		    0x40a808090042L,    0x40a808090043L,    0x40a808092043L,    0x40a8080a8042L,
		    0x40a8080a8043L,    0x40a8080aa043L,    0x40a8080c8042L,    0x40a809010042L,
		    0x40a809010043L,    0x40a814090042L,    0x40a814090043L,    0x40a8140a8042L,
		    0x40a8140c8042L,    0x40a815010042L,    0x40a815010048L,    0x40a815028042L,
		    0x412808090042L,    0x4128080a8042L,    0x4128080c8042L,    0x412809010042L,
		    0x412814090042L,    0x4128140a8042L,    0x4128140c8042L,    0x412815010042L,
		    0x412825010042L,    0x424048090002L,    0x424048090003L,    0x424048090042L,
		    0x424048090043L,    0x424054090002L,    0x424054090003L,    0x424054090042L,
		    0x424054090043L,    0x4240540a8002L,    0x4240540a8042L,    0x424064090002L,
		    0x424064090003L,    0x424064090042L,    0x424064090043L,    0x4240640a8003L,
		    0x4240640a8042L,    0x4240640c8003L,    0x424088090002L,    0x424088090042L,
		    0x4240880a8002L,    0x4240880a8042L,    0x4240880aa003L,    0x4240880c8002L,
		    0x4240880c8008L,    0x4240880c8042L,    0x424094090002L,    0x424094090008L,
		    0x424094090042L,    0x424094090043L,    0x424094090048L,    0x4240940a8002L,
		    0x4240940a8008L,    0x4240940a8042L,    0x4240940aa003L,    0x4240940c8002L,
		    0x4240940c8042L,    0x4240a40a8002L,    0x4240a40a8042L,    0x424809010042L,
		    0x424809010043L,    0x424815010042L,    0x424815010043L,    0x424815012043L,
		    0x424815028042L,    0x424815028043L,    0x424825010042L,    0x424825028042L,
		    0x424848090002L,    0x424848090003L,    0x424848090042L,    0x424848090043L,
		    0x4248480a8002L,    0x4248480a8042L,    0x424848110002L,    0x424848110042L,
		    0x424848128002L,    0x424848128008L,    0x424848128042L,    0x424849090002L,
		    0x42a048090002L,    0x42a048090042L,    0x42a054090002L,    0x42a054090042L,
		    0x42a064090042L,    0x42a094090002L,    0x42a0940a8002L,    0x42a0a40a8002L,
		    0x42a809010042L,    0x42a815010042L,    0x42a815028042L,    0x42a848090002L,
		    0x42a848110042L,    0x42a848128002L,    0x42a849090002L,  0x12024048080043L,
		  0x12024048090001L,  0x12024048090002L,  0x12024048090042L,  0x12024048092003L,
		  0x12024048092043L,  0x15024048090002L,  0x15024048090042L,  0x15024048092003L,
		  0x15024048092043L,  0x1502a048090002L,  0x1502a048090042L,  0x1502a054090002L,
		  0x1502a054090042L,  0x1502a0540a8002L,  0x19024048090002L,  0x19024048090003L,
		  0x19024048090042L,  0x19024048090043L,  0x1902a048090002L,  0x1902a054090002L,
		  0x1902a0540a8002L,  0x19032048090002L,  0x19032048090042L,  0x19032054090002L,
		  0x19032064090002L,  0x190320640a8002L,  0x190320640c8002L,  0x22032048090003L,
		  0x2502a048090007L,  0x2502a054090007L,  0x2502a0540a8007L,  0x2502a0540a8047L,
		  0x25032054090007L,  0x250320540a8007L,  0x250320640a8007L,  0x250320640a804aL,
		  0x2504a048090008L,  0x2504a048090048L,  0x2504a054090008L,  0x2504a0540a8008L,
		  0x2504a0540a8048L,  0x2902a054090007L,  0x2902a0540a8007L,  0x290320540a8007L,
		 0x202404064090043L, 0x212024048090002L, 0x212024048090042L, 0x212024048092003L,
		 0x21202a048090002L, 0x21202a048090042L, 0x21202a054090002L, 0x21202a054090042L,
		 0x21202a0540a8002L, 0x2120320640c8002L, 0x21204a0540a8007L, 0x21204a094090008L,
		 0x21204a094090048L, 0x21204a0940a8008L, 0x21204a0940a8048L, 0x2120520540a8007L,
		 0x212424048090002L, 0x212424048090042L, 0x212424054090002L, 0x212424054090042L,
		 0x212424094128008L, 0x212424094128048L, 0x2124248480a8002L, 0x2124248480a8042L,
		 0x215024048090002L, 0x215024048090042L };

	protected static final float[] INCREMENT = {
		  0.6967f,     0.0f,  0.4886f, -0.4727f, -0.0749f,  0.6262f,  0.2735f,    0.57f,
		   0.701f,  0.9534f, -3.6435f, -2.1509f,  0.4975f, -2.1995f, -0.2809f,  -0.826f,
		 -0.1785f, -1.6203f,  -1.096f,  -0.362f,  0.1395f, -0.2975f, -1.2908f,  1.0162f,
		 -1.3825f, -2.5384f,  0.3317f,  0.4292f, -0.5824f, -0.1834f,  0.1306f, -0.5015f,
		  0.6262f, -0.5258f,  0.4244f,  -0.161f, -0.2778f,  0.5111f, -0.4357f, -0.1041f,
		  0.3424f, -0.0615f,  0.6035f,  0.7227f,  0.4346f, -1.6821f,  -0.331f,  -0.498f,
		 -1.4915f, -0.3651f,  0.4597f,  0.3384f,  0.6633f,  0.4544f,  0.1597f,  0.6339f,
		  0.3504f,  0.0449f,   0.342f,  0.2611f,  0.4046f,  0.5219f,  -1.065f, -1.2314f,
		 -1.8023f, -0.3632f, -0.4108f,  0.3057f, -0.1456f, -0.2713f, -0.5193f,  0.4526f,
		  0.5539f,  0.8374f,  -0.707f, -0.4881f,   -0.41f,     0.0f,  0.1479f,  0.3448f,
		  0.4298f,  0.5579f, -0.1265f, -0.0425f,  0.0767f,  0.6635f, -0.3812f, -0.8368f,
		  1.0287f, -0.1021f,  0.3587f, -0.5945f,  0.1692f, -0.1218f,  0.3283f,  0.2239f,
		  0.2043f,   0.059f, -0.4835f,  0.6165f, -0.4011f,  0.5578f, -0.2164f, -0.0175f,
		  0.2981f,    0.11f,  0.2715f, -0.2995f,  -0.467f,  0.1566f,  0.0468f, -0.1321f,
		  1.3686f,     0.0f, -0.4116f,  1.0186f, -0.3935f,  0.5223f,  0.3685f,  0.2577f,
		  1.5193f,  0.2705f,  0.3791f,   0.012f, -0.3397f,  0.1483f,  0.2766f,  0.3593f,
		  0.7715f,   0.315f, -1.6185f,  0.1889f, -0.2652f, -0.0965f,  0.4202f,  0.1871f,
		 -0.3684f, -0.0778f,  0.8943f,  0.3694f,  0.2879f,  0.4489f, -0.2601f,  0.4771f,
		  0.1923f,  0.4381f,  0.1695f,  0.4525f,  0.3352f,  0.1583f,  0.4036f,  -0.048f,
		  0.5023f, -0.2649f,  0.7691f, -0.3552f,  1.0301f, -0.1141f, -0.5932f,  0.1749f,
		  0.1313f, -0.1804f,  0.3994f,  0.2291f,  0.3169f,  0.3599f, -0.0039f, -0.2956f,
		  0.4409f, -0.1609f,  0.3775f, -0.1346f,  0.2839f,  0.5129f,  0.1266f,  0.4294f,
		  0.2806f,  0.4907f,   0.354f,  0.2192f,  0.1565f,  0.6935f,  0.3618f,  0.6735f,
		  0.5778f, -0.5636f,  0.5569f,  0.3038f, -0.3276f, -0.6992f,  0.0103f, -0.4659f,
		  0.3818f,  0.3317f,  0.1841f,  0.7071f,  0.1227f,  0.7949f, -0.6592f, -1.3148f,
		 -0.4067f, -0.1316f, -0.4925f, -0.0929f, -0.4352f, -0.2207f,  -0.996f, -0.7238f,
		 -0.5469f, -1.2939f, -0.0136f,  0.0657f,  0.7189f,   0.057f,  0.6619f, -0.6381f,
		 -0.8073f,  0.2355f,  0.3048f, -0.0199f, -0.0752f,  0.4461f,  0.1559f,  1.1168f,
		 -1.8039f,  0.2365f, -0.2206f,   0.448f,  -1.134f, -0.5653f, -0.4053f, -0.1361f,
		  0.2199f,  0.0536f,  -0.021f,  0.6985f,  0.9643f, -0.4152f, -1.0369f,  -0.183f,
		  0.5883f, -0.2918f, -0.5294f, -0.6541f,  0.9473f, -0.1906f, -0.8484f, -0.3457f,
		  0.9541f,  1.4231f, -0.7924f,  -0.602f,    0.08f, -0.2596f,  0.8382f, -0.4416f,
		 -0.3704f, -0.7487f, -0.1079f, -0.2992f, -0.3277f,  0.0251f, -0.9188f,  0.1094f,
		  0.8231f, -3.2333f,   0.035f,  0.3818f,  -0.738f,  0.2791f,  0.3206f,  0.5662f,
		 -0.3784f,  0.4032f, -1.7948f, -0.1554f,  0.3785f,  0.0534f, -0.1653f, -0.0987f,
		 -0.1005f, -0.6461f,  0.8035f, -0.2405f, -0.1238f, -0.3576f,  0.0961f, -0.6401f,
		  0.2029f,  0.2359f,  0.4951f,  0.1921f, -0.3745f,  0.3463f,  0.2899f, -0.1533f,
		  -0.417f,   0.377f,  0.6998f,   0.594f,  0.5912f, -0.5571f,  0.0238f, -0.2475f,
		  0.0307f, -0.3875f, -0.7437f,  0.5144f,  0.0057f,  0.7655f,   0.172f, -2.5624f,
		 -0.3066f,  0.3647f,  0.1727f, -0.0329f, -0.1893f,  0.0702f, -1.2454f,  0.1496f,
		 -1.3825f,  0.4146f, -0.2668f, -0.1106f,  0.0362f, -0.3189f, -0.7278f, -0.0894f,
		 -0.2277f, -0.2394f,  1.0112f, -0.2962f,  0.7776f,  0.2945f, -0.2234f,  0.2764f,
		  0.8011f, -0.1744f,  0.1581f, -1.7552f, -0.3848f,  0.5993f,  0.5268f, -0.0417f,
		  0.4733f, -0.3401f,  -0.145f,  0.7088f, -0.1318f,  0.0426f, -0.5028f,  0.3832f,
		 -0.0118f, -0.4358f,  0.3749f, -0.1203f, -0.5648f, -0.1972f, -0.8769f, -0.3675f,
		 -0.2004f,  -0.607f, -0.1857f,  0.3468f, -0.3624f,  0.5358f, -0.3701f,  0.1336f,
		  0.9545f,  0.1139f, -0.1699f,  0.3317f,  0.2891f,  0.2613f, -0.0344f, -1.9498f,
		 -2.0249f, -0.6005f, -0.6258f, -1.2349f,   0.328f, -0.5434f, -0.7712f, -0.9057f,
		 -0.1668f, -0.9905f, -0.0372f, -1.1638f,  0.1262f, -0.5248f, -0.1538f, -0.3682f,
		  0.3249f,   0.065f,  0.0511f, -0.4607f,  0.2231f,  0.2822f,  0.1397f, -0.4938f,
		  0.3948f, -0.4075f, -0.6411f, -0.0091f, -0.1333f, -0.5192f, -0.1661f,  0.3317f,
		 -0.0707f,  0.4806f,  0.3828f,  0.2229f,   0.616f,  0.3371f,  0.1884f,  0.1381f,
		 -1.4915f,  0.2833f, -0.1226f, -3.9189f, -0.4592f, -0.3435f, -0.6654f, -0.5056f,
		 -0.8631f,  0.1536f, -0.6427f, -0.0884f, -0.0471f,  0.1106f,  0.3821f, -0.2392f,
		 -0.4051f,  0.0891f, -0.6972f, -0.4699f,  0.0922f,  0.0806f, -0.6774f, -0.0622f,
		   -0.93f,  0.1337f };


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
				long type = AtomTypeCalculator.getAtomType(mol, atom, ATOM_TYPE_MODE);
				int index = sSortedTypeList.getIndex(type);
				if (index != -1)
					cLogP += INCREMENT[index];
//				else
//					System.out.println("type not found for atom "+atom+":"+Long.toHexString(type)+" "+AtomTypeCalculator.getTypeString(type, ATOM_TYPE_MODE));
				}
			catch (Exception e) {}	// unsupported atom type exceptions are tolerable
			}
		
		return cLogP;
		}

	/**
	 * Normalizes ambiguous bonds and assigns cLogP increments to every atom
	 * based on its enhanced atom type.
	 * @param mol
	 * @param increment not smaller than non-H atom count of mol
	 * @return
	 */
	public void getCLogPIncrements(StereoMolecule mol, float[] increment) {
		mol.normalizeAmbiguousBonds();
		mol.ensureHelperArrays(Molecule.cHelperRings);

		for (int atom=0; atom<mol.getAtoms(); atom++) {
			try {
				int index = sSortedTypeList.getIndex(AtomTypeCalculator.getAtomType(mol, atom, ATOM_TYPE_MODE));
				if (index != -1)
					increment[atom] = INCREMENT[index];
			}
			catch (Exception e) {}	// unsupported atom type exceptions are tolerable
		}
	}

	public ParameterizedStringList getDetail(StereoMolecule mol) {
		ParameterizedStringList detail = new ParameterizedStringList();
		detail.add("cLogP Values are estimated applying an atom-type based increment system.",
							ParameterizedStringList.cStringTypeText);
		detail.add("Atom-types are 64-bit numbers describing atoms and their near surrounding.",
							ParameterizedStringList.cStringTypeText);
		detail.add("Recognized atom types and their contributions are:",ParameterizedStringList.cStringTypeText);

		if (mol != null) {
			mol.normalizeAmbiguousBonds();
			mol.ensureHelperArrays(Molecule.cHelperRings);

			int errorCount = 0;
			TreeMap<Long,Integer> countMap = new TreeMap<Long,Integer>();
			NumberFormat formatter = new DecimalFormat("#0.000");
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
					detail.add(countMap.get(type) + " * "+
							formatter.format(INCREMENT[sSortedTypeList.getIndex(type)]) + " AtomType: 0x" + Long.toHexString(type),ParameterizedStringList.cStringTypeText);
				else
					detail.add("Warning: For atom type 0x"+Long.toHexString(type)+" ("+countMap.get(type)+" times found) is no increment available.", ParameterizedStringList.cStringTypeText);
				}
			}
		
		return detail;
		}
	}