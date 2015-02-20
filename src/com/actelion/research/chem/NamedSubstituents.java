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

package com.actelion.research.chem;

import java.util.TreeMap;

public class NamedSubstituents {
	private static final String[][] SUBSTITUENT_LIST = {
			{ "Ac", "gCaHA`AIf`@" },
			{ "Alloc", "gNph@l@ILzuR@@" },
			{ "Allyl", "gC`HL@IVt@@" },
			{ "Bn;Bzl;Benzyl", "daD@`F@DjUZxHH@@" },
			{ "Boc", "daxD`@S@AIgijj@@" },
			{ "BOM;BzOM", "deTH`@H@Re[TYj`@@@" },
			{ "Bs", "dmtDpAdLS`aPAIe]jf`@e`@@" },
			{ "Bt", "dew@`@aJ@DiY]paej`@@@" },
			{ "Btm", "did@P@BJ@Die_ahH@@@" },
			{ "Bu;n-Bu", "gJPHB@IRuP@" },
			{ "Bz;Benzoyl", "didH`@p@RYm^Eh@@@@" },
			{ "Bzh", "dg|@`N@LdbRbtJUB]aAP@@@@" },
			{ "Cbz", "dmtD`@S@AIgYVUZh@@@@" },
			{ "Cy", "gOpHL@IToWUU@@" },
			{ "cyclobutyl", "gKPHL@IThuT@@" },
			{ "cycloheptyl", "daD@`L@DjWVzjj`@" },
			{ "cyclooctyl", "did@`L@DjWWajjj@@" },
			{ "cyclopentyl", "gFpHL@ITimUP@" },
			{ "cyclopropyl", "gBPHL@Qxjh@" },
			{ "DEAE", "daz@`@x@RiUjj`@" },
			{ "DEIPS", "diD@P@\\B@DjfVjj`@" },
			{ "DMIPS", "gNpD@xD@RjZjh@" },
			{ "DMPM", "dcLD`@kPCIEMDdcttDDT@@" },
			{ "DMPS", "deT@P@\\B@LddTjPsU@@@@" },
			{ "DMTr", "fak@b@@Mt@ISZ{SMjo{NQKfm@AU@@@E@@@@" },
			{ "DNP", "dkmB`hdDt~@HeNfS{HihheCAUhBHX@@" },
			{ "DNS;Dan", "fhi`a@KPP@HH@YIHYheEhYKQgKP@@QP@@@" },
			{ "DPIPS", "fdyAA@H@\\B@FRRIQSQIHzp_Qjh@h@@@@@" },
			{ "DPTBS", "fleAA@H@\\B@FRRIQSRIIWNbEMU@EP@@@@@" },
			{ "DTBMS", "dmT@P@\\B@Djffjjjh@@" },
			{ "Et", "eMBD@ch@" },
			{ "Fmoc", "fde@b@@Hp@IL{LrjxeVCzKUT@@@P@@@" },
			{ "i-Am", "gGPHJ@YIDZj@@" },
			{ "i-Bu", "gJPHB@ITuP@" },
			{ "Im", "gFtHAj@IRnKSP@" },
			{ "i-Pr", "gC`HL@Qz`@" },
			{ "MDIPS", "diD@P@\\B@DjfZjj`@" },
			{ "MDPS", "foA@A@@NA@CIIEEBdeeVLzj@@@@@@" },
			{ "Me", "eFBH@c@@" },
			{ "MEM", "gNphAR@IRoUT@@" },
			{ "Mes", "deT@`J@DjY{[`bB`@@" },
			{ "MMTr", "ffcAB@@Z@Dim]ifuWYrI\\uh@Jh@@@@@@" },
			{ "MOM", "gCaHA`AJZ`@" },
			{ "MPM;PMB", "deTH`@d@Rfuunh@J@@" },
			{ "Ms", "gJPdH`DD@cuh@" },
			{ "MTM", "gC`D@DX@Rfh@" },
			{ "m-Tolyl", "daD@`N@DjWjXHB@@" },
			{ "N3", "gClHaE`@RnReX@" },
			{ "n-Am;Am", "gGPHJ@IRmU@@" },
			{ "neo-Am", "gGPHJ@IUMU@@" },
			{ "nitro;NO2", "gChhhE`BRnRYh@" },
			{ "Np", "deVDaHAI@HeNR[e_aZ@B@@" },
			{ "n-Pr;Pr", "gC`HL@IST@@" },
			{ "o-Tolyl", "daD@`J@DjYvxH`@@" },
			{ "Ph;Phenyl", "gOpHL@IToVD@@@" },
			{ "Pht", "dcLL`@RU@Dfyed]ZBA`@@" },
			{ "Piv;Pv", "gNqHA`AIffj`@" },
			{ "PMBM", "dcLD`@T`AJUm]FZh@J@@" },
			{ "PNB", "dcNLaHAEt@bTyInUvxV`@f@@" },
			{ "Poc", "didD`@S@AIgexVjj`@" },
			{ "PPi", "diDFsHSB[`|J|A@Lxn{lddqdZih@@" },
			{ "p-Tolyl", "daD@`N@DjWzXHB@@" },
			{ "s-Am", "gGPHL@YIDZj@@" },
			{ "s-Bu;s-Butyl", "gJPHL@ITuP@" },
			{ "SEM", "diDHPFApD@rRQUJjj`@" },
			{ "SES", "dedDpHP@``AgCIICeHmUT@@" },
			{ "t-Am", "gGPHB@IUMU@@" },
			{ "TBDMS;TBS", "dax@P@\\B@Djfjjh@@" },
			{ "TBDPS", "fdy@A@@NA@CIIEEEIde]XOhuPAT@@@@@" },
			{ "TBMPS", "dg\\HPHApH@rRQJJPjg]UAT@@@" },
			{ "t-Bu;t-Butyl", "gJPHB@Q}T@@" },
			{ "TDS", "ded@P@\\B@LddTeeUUP@@" },
			{ "Tf", "daxDhHP@``BiAiCiCIICHmU@@" },
			{ "TFA", "gNqBJIARFdF@YEHYUL@@" },
			{ "Thexyl", "gNpHB@IUMUT@@" },
			{ "THF", "gFqH@PAJYujj@@" },
			{ "THP", "gOqH@PAJYZzjh@" },
			{ "TIPS", "dmT@P@\\B@DjfYjjjh@@" },
			{ "TMS", "gJPD@xD@czh@" },
			{ "Tos;Ts", "dmtDPHP@``CIICLeaeZ@B@@" },
			{ "Troc", "diDDHJxHaHcH`PCHiBeJjf@@" },
			{ "Trt", "fbm@B@A@FRQIRKQPiIZdoIcdHJ`@@@@@@" },
			{ "Xyl", "did@`J@DjYynBHH@@" },
			};
	private static TreeMap<String,String> sIDCodeMap = null;

	private static void createMap() {
		sIDCodeMap = new TreeMap<String,String>();
		for (String[] sp:SUBSTITUENT_LIST) {
			String[] key = sp[0].split(";");
			for (String k:key)
				sIDCodeMap.put(normalize(k), sp[1]);
			}
		}

	private static String normalize(String s) {
		return s.toLowerCase().replace("-", "");
		}

	public static String getSubstituentIDCode(String name) {	
		if (sIDCodeMap == null)
			createMap();

		return sIDCodeMap.get(normalize(name));
		}

	public static StereoMolecule getSubstituent(String name) {	
		String idcode = getSubstituentIDCode(name);
		return (idcode == null) ? null : new IDCodeParser(false).getCompactMolecule(idcode);
		}

	/**
	 * @param s
	 * @return whether appending some characters to s may lead to a known substituent name
	 */
	public static boolean isValidSubstituentNameStart(String s) {	
		if (sIDCodeMap == null)
			createMap();

		s = normalize(s);
		for (String n:sIDCodeMap.keySet())
			if (n.startsWith(s))
				return true;
		return false;
		}
	}
