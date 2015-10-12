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

import java.util.ArrayList;

public class SSSearcherWithIndex {
	public static final String cIndexVersion = "1.2.1";

	public static final String[] cKeyIDCode = {
        "QM@HzAmdqjF@",
        "RF@Q``",
        "qC`@ISTAlQE`",
        "`J@H",
        "QM@HzAmdqbF@",
        "qC`@ISTAlQEhqPp@",
        "sJP@DiZhAmQEb",
        "RF@QPvR@",
        "QM@HzA@",
        "qC`@ISTAlQEhpPp@",
        "qC`@Qz`MbHl",
        "sJP@DiZhAmQEcFZF@",
        "RFPDXH",
        "qC`@IVtAlQE`",
        "QM@HvAmdqfF@",
        "sGP@DiVj`FsDVM@",
        "`L@H",
        "sJP@DizhAmQEcFBF@",
        "sJP@DjvhAmQEb",
        "sFp@DiTt@@AlqEcP",
        "sGP@LdbMU@MfHlZ",
        "QMHAIhD",
        "QM@HzAy@",
        "sJP@DkVhAmQEb",
        "sNp@DiUjj@[\\QXu`",
        "sJP@DiZhAmQEcFBF@",
        "sGP@DjVj`FsDVM@",
        "RFPDTH",
        "RG@DXOH@",
        "sGP@Divj`FsDVMcAC@",
        "sGP@Dj}j`FsDVM@",
        "qC`@Qz`MbHmFRF@",
        "sNp@LdbJjj@[\\QXu`",
        "QMHAIhGe@",
        "QM@HzAyd`",
        "QM`AIhD",
        "qC`@ISTA@",
        "sGP@DkUj`FsDVM@",
        "qC`@IVtAlQEhqPp@",
        "sNp@DiUjj@[\\QXuqea`@",
        "KAx@@IRjuUPAlHPfES\\",
        "QM`BN`P",
        "sJP@DjZhAmQEcFJF@",
        "Hid@@DjU^nBBH@FtaBXUMp`",
        "sNp@Diujj@[\\QXuq`a`@",
        "sJP@DjvhAmQEcFZF@",
        "sJP@DjZhAmQEcFFF@",
        "sOp@DjWkB@@FwDVM\\YhX@",
        "sNp@Dj}Zj@[\\QXu`",
        "sNp@DiWjj@[\\QXuq`a`@",
        "sOp@DjWkB@@D",
        "KAx@@ITouUPAlHPfES\\",
        "KAx@@YIDTjjh@vDHSBin@",
        "sNp@DkUZj@[\\QXu`",
        "RFPDXOH@",
        "QM`BN`^L`",
        "qC`@ISTAy@",
        "sGP@LdbMU@MfHl[FVF@",
        "qCb@AIZ`H",
        "KAx@@IRjuUPAlHPfES]FFa`@",
        "KAx@@ITnuUPAlHPfES\\",
        "HiD@@DiUVjj`AmHPfES\\H",
        "sNp@DjUjj@[\\QXu`",
        "sJP@DkVhAmQEcFJF@",
        "sGP@DjVj`FsDVMcCC@",
        "qC`@Qz`MbHmFBF@",
        "sJP@DkfhAmQEb",
        "qC`@IVtAlQEhsPp@",
        "sGP@Djuj`FsDVM@",
        "sGP@Dj}j`FsDVMcMC@",
        "sJP@DiZhA@",
        "KAx@@ISjuUPAlHPfES]F@a`@",
        "sJP@DjZhAmQEcFRF@",
        "KAx@@IRnuUPAlHPfES]F@a`@",
        "HiD@@DjWvjj`AmHPfES\\H",
        "QMHAIhGd@",
        "sNp@DiUjj@[\\QXuq`a`@",
        "KAx@@IVjmUPAlHPfES\\",
        "sGP@DjVj`FsDVMcMC@",
        "QM`AIhGe@",
        "HiD@@LdbJRjjh@[RDIaTwB",
        "qCp@AIZ`H",
        "sGP@LdbMU@MfHl[FFF@",
        "QMDARVA@",
        "sNp@LdbJjj@[\\QXuqba`@",
        "sNp@LdbJjj@[\\QXuqca`@",
        "sGP@Dkej`FsDVM@",
        "qCb@AIZ`OI@",
        "HaD@@DjUZxHH@AlHPfES]FLa`@",
        "sGP@DkYj`FsDVM@",
        "qCb@AIV`H",
        "sNp@LdbJjj@[\\QXuqea`@",
        "sGP@DkUj`FsDVMcEC@",
        "sFp@DiTt@@Axa@",
        "Hmt@@DjU_ZxHHj@AmhPfES\\Lj",
        "QM`BN`^P",
        "qCb@AIZ`OH`",
        "sFp@DiTt@@AxaP",
        "sGP@Djuj`FsDVMcEC@",
        "sGP@Djuj`FsDVMcIC@",
        "sGP@DkUj`FsDVMcKC@",
        "sJP@DkfhAmQEcFRF@",
        "sGP@DjVj`FsDVMcIC@",
        "HaD@@DjUZxHH@AlHPfES]FFa`@",
        "qC`@IRtDVqDV@",
        "sNp@Dj}Zj@[\\QXuqfa`@",
        "KAx@@ITnuUPAlHPfES]FFa`@",
        "HiD@@DkUUjj`AmHPfES\\H",
        "sJQ@@dkU@H",
        "qC`@Qz`H",
        "KAx@@IUkmUPAlHPfES\\",
        "KAx@@ITouUPAlHPfES]FJa`@",
        "sJP@H~j@[TQX`",
        "sGP@DjZj`FsDVM@",
        "sJP@DkVhAmQEcFFF@",
        "sJX@@eKU@H",
        "sJP@DizhAy@",
        "QMHAIhGbP",
        "KAx@@ITouUPAlHPfES]FNa`@",
        "HaD@@DjUZxHD@AlHPfES\\",
        "HaD@@DjUZxHH@A@",
        "sNp@LdbJjj@[\\QXuqaa`@",
        "Hed@@LdbRQUUUP@vTHSBinFP",
        "KAx@@ITouUPAlHPfES]FLa`@",
        "sNp@DkUZj@[\\QXuqba`@",
        "KAx@@ITjuUPAlHPfES]FNa`@",
        "KAx@@YIDTjjh@vDHSBincGPp@",
        "HaD@@DjYvxH`@AlHPfES]FLa`@",
        "RF@QP`",
        "qCb@AIj`H",
        "sNp@DjUjj@[\\QXuqaa`@",
        "sNp@DkVZj@[\\QXu`",
        "KAx@@YIDUJjh@vDHSBin@",
        "sGP@DkYj`FsDVMcIC@",
        "sGP@DjVj`FsDVMcAC@",
        "sGP@DiVj`D",
        "sJP@DkVhAmQEcFZF@",
        "sNp@LdbLjj@[\\QXu`",
        "QM@HvAmdqbF@",
        "HaD@@DjWjXHB@AlHPfES\\",
        "sNp@DjwZj@[\\QXuqba`@",
        "sNp@LdbJjj@[\\QXuqda`@",
        "sFp@DiTt@@Axa`",
        "HiD@@Djuujj`AmHPfES\\H",
        "sNp@DkUZj@[\\QXuqca`@",
        "sJP@DiZhAy@",
        "KAx@@YIDTjjh@vDHSBincCPp@",
        "KAx@@IWNmUPAlHPfES\\",
        "KAx@@IVkMUPAlHPfES\\",
        "sJQ@@dju@H",
        "qCb@AIZ`OH@",
        "qC`@ISTAxa@",
        "sNp@DjyZj@[\\QXu`",
        "Hid@@DjUfaBB`@FtaBXUMp`",
        "HiD@@DiUVjj`AmHPfES\\LXBF@",
        "KAx@@IUjmUPAlHPfES\\",
        "HiD@@DjWvjj`AmHPfES\\LXjF@",
        "sJP@DjVhAmQEb",
        "qCb@AIV`OH`",
        "HiD@@LdbJRjjh@[RDIaTwCFDa`@",
        "KAx@@YIDTjjh@vDHSBinc@Pp@",
        "sNp@DjUjj@[\\QXuqda`@",
        "qC`@Qz`OED",
        "sJP@DkfhAmQEcFZF@",
        "KAx@@YIDbjjh@vDHSBincDPp@",
        "sGP@Djyj`FsDVMcMC@",
        "KAx@@IVrmUPAlHPfES\\",
        "qCp@AIZ`OI@",
        "sJX@@dkU@H",
        "sJQ@@dkU@OH`",
        "sNp@Di]ZjBBvxbqk@",
        "Hkl@@DjU_Uk``bj`@[VDIaTwCJzX",
        "sGP@DjZj`FsDVMcEC@",
        "Hid@@DjU^nBBH@FtaBXUMpqcHX@",
        "sNp@DkeZj@[\\QXu`",
        "sNp@DjYjj@[\\QXuqca`@",
        "sGQ@@djuT@`",
        "HiD@@LdbJTjjh@[RDIaTwB",
        "sOp@DjWkB@@Gd`",
        "HeT@@LdbbRKBDQD@CYPaLJfxY@",
        "qCr@XIKTA@",
        "HiD@@DjW^jj`AmHPfES\\LXJF@",
        "HeT@@DjU]k``b`@[JDIaTwCH",
        "sGP@Djuj`FsDVMcCC@",
        "`IH`B",
        "sOp@DjWkB@@GdX",
        "sJQ@@eKU@H",
        "KAx@@YIDUJjh@vDHSBincBPp@",
        "sJX@@eKU@OH@",
        "KAx@@YIDTjjh@vDHSBincAPp@",
        "sOq@@drm\\@@@`",
        "KAx@@IUkMUPAlHPfES\\",
        "qCp@AIj`H",
        "Hed@@DjUUjjj@FraBXUMpr",
        "sGX@@eJuT@`",
        "sGP@DkUj`FsDVMcCC@",
        "HiD@@Dj}Ujj`AmHPfES\\LXrF@",
        "KAx@@ITouUPAlHPfES]FHa`@",
        "Hed@@DjWujjj@FraBXUMpsFIa`@",
        "sGP@DiUj``mfHlZ",
        "sFp@DiTvjhAlqEcP",
        "Hid@@DjU^nBBH@FtaBXUMpq`XX@",
        "sJP@DkVdAmQEb",
        "qCp@AIZ`OH`",
        "QMhDRVA@",
        "qC`@ISJAlQE`",
        "qCp@BOTAyhl",
        "sJX@@eOU@ODB",
        "sFp@DiTt@@AyaB",
        "sGP@DkUj`FsDVMcMC@",
        "Hid@@DjYUaBH`@FtaBXUMpqcHX@",
        "qC`@Qz`OH@",
        "HiD@@DjUVjj`AmHPfES\\LXZF@",
        "sJP@H~j@[TQXqda`@",
        "sJX@@eKU@OI@",
        "sNp@Djejj@[\\QXu`",
        "sJQ@@dsU@H",
        "sJQ@@dkU@OI`",
        "KAx@@YIMDVjh@vDHSBin@",
        "Hid@@DjU^nBBD@FtaBXUMp`",
        "sNp@DkgZj@[\\QXuqca`@",
        "qC`@IRtDVqDVcEC@",
        "Hed@@LdbRQeUUP@vTHSBinFP",
        "sNp@DiUjj@P",
        "qC`@IRtDT",
        "sNp@DkYZj@[\\QXuqca`@",
        "KAx@@IUkmUPAlHPfES]FDa`@",
        "KAx@@IVjmUPAlHPfES]FNa`@",
        "sOx@@drm\\@@@`",
        "KAx@@ITjuUPAlHPfES]FBa`@",
        "QMDARVAyH",
        "sJP`@dfvhA@",
        "HeT@@DjU_k``b`@[JDIaTwCLXfF@",
        "KAx@@IToUUPAlHPfES]FJa`@",
        "sGP@DkYj`FsDVMcEC@",
        "qCb@AIZ`ODH",
        "`I@`B",
        "KAx@@IUzmUPAlHPfES]FFa`@",
        "sNp@DkfZj@[\\QXu`",
        "KAx@@ITnuUPAlHPfES]F@a`@",
        "HiD@@LddURjjh@[RDIaTwB",
        "sNp@Dj~Zj@[\\QXuqfa`@",
        "Hed@@Dj{uZjj@FraBXUMpr",
        "KAx@@ITsUUPAlHPfES\\",
        "Hid@@LdbRQk``b@AmHPfES\\LXrF@",
        "sOp@DjWkB@@GdH",
        "sJQ@@dkU@OH@",
        "Hid@@DjU^nBBH@FtaBXUMpqahX@",
        "sGP@DiYj``mfHlZ",
        "KAx@@IToUUPAlHPfES]FLa`@",
        "qCp@AJZ`ODH",
        "Hmt@@DjU]ZxHHj@AmhPfES\\Lj",
        "sGP@DkUjPFsDVM@",
        "qC`@IVtA@",
        "Hed@@LdbJReUUP@vTHSBinFP",
        "sNp@DjuZj@[\\QXuqea`@",
        "KAx@@IUkmUPAlHPfES]FNa`@",
        "HiD@@DkVUjj`AmHPfES\\H",
        "Hed@@DkUeZjj@FraBXUMpr",
        "sNp@DkVZj@[\\QXuqea`@",
        "sJP@DiVhHKZbKFLLL@",
        "HiD@@Djuyjj`AmHPfES\\H",
        "sNp@DjUjj@[\\QXuq`a`@",
        "HeT@@DjYUXPbH`@[JDIaTwCH",
        "HiD@@DjwUjj`AmHPfES\\LXRF@",
        "sNq@@djmUPB",
        "KAx@@YIEEZjh@vDHSBincCPp@",
        "sGP@Di^V`dmfHlZ",
        "Hid@@DjYUaBHP@FtaBXUMp`",
        "sNp@DjYjj@[\\QXuqba`@",
        "sGP@Dkej`FsDVMcKC@",
        "HeT@@DjU^k``b`@[JDIaTwCH",
        "qC`@Qv`MbHmFBF@",
        "sGQ@@djmT@`",
        "qCr@XIKTAyH",
        "qC`@IVtAlQEhpPp@",
        "Hid@@LdbbQxXF@@AmHPfES\\LXjF@",
        "sGP@DkYj`FsDVMcCC@",
        "KAx@@IVsMUPAlHPfES\\",
        "qCp@AIj`ODl",
        "HiD@@DkeUjj`AmHPfES\\H",
        "HeT@@DjU[kjjjh@ZLDXSSYPaLJfxY@",
        "sJP@DkVdAmQEcFRF@",
        "HiD@@LdbJTjjh@[RDIaTwCFDa`@",
        "HiD@@DkYyjj`AmHPfES\\H",
        "sJP@DjZhAyH",
        "KAx@@IVkMUPAlHPfES]FDa`@",
        "sJX@@dkU@OI@",
        "Hed@@LdbRQUUUP@vTHSBinFXpLL@",
        "Hed@@DjuUZjj@FraBXUMpr",
        "sGP@Djfj`FsDVMcKC@",
        "sNp@DkVZj@[\\QXuqba`@",
        "sNp@DjyZj@[\\QXuqfa`@",
        "qCb@AIj`OH@",
        "sNp@DjUZj@[\\QXu`",
        "KAx@@IWOMUPAlHPfES\\",
        "Hid@@DjU^nBBH@D",
        "Hed@@DjuvZjj@FraBXUMpr",
        "sJP@DiVhHKZbKFLtL@",
        "Hmt@@DjU_Zzjjj`AhpQaLmmBDpj[aeXplL@",
        "sNp@DjuZj@[\\QXuqca`@",
        "sJP@DkfhAmQEcFJF@",
        "sNp@LdbJZj@[\\QXu`",
        "HeT@@DjU_k``b`@[JDIaTwCLXFF@",
        "KAx@@IVlmUPAlHPfES]FNa`@",
        "HeT@@LdbbRKBDQD@CYPaLJfxYcEPp@",
        "Hid@@DjUZnBBH@FtaBXUMpqcHX@",
        "qCa@CIKTA@",
        "HiD@@Dj~]jj`AmHPfES\\LXFF@",
        "sKP@Di\\Zj@[TQX`",
        "sGP@Djfj`FsDVMcEC@",
        "HiD@@DkgYjj`AmHPfES\\H",
        "sNp@DjuZj@[\\QXuqaa`@",
        "KAx@@YIMDVjh@vDHSBincDPp@",
        "sJP@DjVhHKZbKFLTL@",
        "Hid@@LdbRQk``b@AmHPfES\\LXZF@",
        "HiD@@Dj}Ujj`AmHPfES\\LXzF@",
        "HeT@@DjU_k``bP@[JDIaTwCH",
        "sNp@DkUZi@[\\QXu`",
        "HiD@@DjYfjj`AmHPfES\\H",
        "sGP@DjZj`FsDVMcAC@",
        "Hmt@@DjU_jxHHj@AmhPfES\\Lj",
        "Hid@@LdbRQk``R@AmHPfES\\H",
        "KAx@@YIDUJjh@vDHSBincDPp@",
        "qCr@XIKTAyD",
        "sOq@@drm\\@@@|`@",
        "Hed@@DjW^jjj@FraBXUMpsFBa`@",
        "HeT@@DjY]zXFB@@[JDIaTwCH",
        "Hkl@@DjU_Vk``bj`@[VDIaTwCJzX",
        "Hid@@DjY}nBHH@FtaBXUMpqcHX@",
        "sGX@@eKuT@|d@",
        "sGP@Dj^Y`FsDVM@",
        "HcL@@DjU_ZnBBJh@FqaBXUMprn`",
        "sJP@DkVdAmQEcFJF@",
        "sOq@@drm\\@@@|b@",
        "sNp@DjyZj@[\\QXuqaa`@",
        "HaD@@DjUZxHH@AyD@",
        "qC`@Qv`H",
        "Hmt@@DjU_Zzjjj`AhpQaLmmBDpj[aeXqdL@",
        "sGP@Dkej`FsDVMcMC@",
        "Hed@@DjUUjjj@FraBXUMpsFHa`@",
        "HeT@@LdbbRkBDQD@CYPaLJfxY@",
        "KAx@@IU{MUPAlHPfES]FLa`@",
        "RG@DTH",
        "sJY@DDeVhA@",
        "KAx@@YIDUJjh@vDHSBinc@Pp@",
        "sJX@@dkU@OI`",
        "sJQ@@dju@OI`",
        "HeT@@LdbbRKBDQD@CYPaLJfxYcFPp@",
        "sFp@DiTvjhAlqEcXpPp@",
        "HaD@@DjUZxHH@AyG@",
        "sNx@@eJ}UPB",
        "sNp@LddUjj@[\\QXuqca`@",
        "HaDH@@RVU[j@@@D",
        "sNp@DkgZi@[\\QXu`",
        "sGY@LDeVj`D",
        "sNp@LdbJfZBZvxbqk@",
        "sJP`@dfvhAyL",
        "sGX@AddQjhAxe`",
        "Hmt@@DjU_ZxHHj@AmhPfES\\LkFIa`@",
        "qCh@CIKTA@",
        "sNp@LdbLjj@[\\QXuq`a`@",
        "sOq@@drm\\@@@|a@",
        "KAx@@IUzmUPAlHPfES]FJa`@",
        "sNx@AddQUUPB",
        "sGP@Di]jP`mfHlZ",
        "sJP`@TeZhA@",
        "KAx@@IRjmUPHKXPaLJfx",
        "HeT@@LdbRTM\\DDT@CYPaLJfxY@",
        "HaF@@@Rfu[j@@@D",
        "Hid@@DjYUaBH`@FtaBXUMpqchX@",
        "KAx@@IUjmTpAlHPfES\\",
        "Hid@@DjU^nBBD@FtaBXUMpqcHX@",
        "sGP@DiUj``mfHl[FFF@",
        "KAx@@IUvmUPAlHPfES]FLa`@",
        "Hed@@LdbQTUUUP@vTHSBinFXqDL@",
        "sJP@DkVhA@",
        "sOx@@drm\\@@@|b@",
        "KAx@@IUkMUPAlHPfES]FDa`@",
        "HeT@@LdbRQU\\DDT@CYPaLJfxY@",
        "HiD@@Dj}Yjj`AmHPfES\\LXrF@",
        "HiD@@Dj{ujj`AmHPfES\\LXFF@",
        "KAx@@IWNmUPAlHPfES]FFa`@",
        "KAx@@IRkMUPHKXPaLJfx",
        "sJP@DjYdAmQEcFZF@",
        "sJY@LDeZhAyL",
        "HaDH@@RVU[f@@@D",
        "sJP`@deVhAyB",
        "HaD@@DjWjZjj`AlHPfES\\",
        "sGP@DkYj`FsDVMcMC@",
        "sNp@DkgZj@[\\QXuqea`@",
        "sJQ@@dlu@H",
        "HeT@@DjU]k``b`@[JDIaTwCLXrF@",
        "sJX@@dkU@OH`",
        "RFDDQFCr`",
        "sJP@DiYXIKZbKFLLL@",
        "KAx@@YIHjjjh@vDHSBincGPp@",
        "Hk\\@@DjU^ukmLHH@@@AmXPfES\\Lki`",
        "sGQ@@djmT@|b@",
        "Hid@@DjUfaBB`@FtaBXUMpqahX@",
        "sNx@@eRmUPB",
        "Hmt@@LdbRVak``ah@FvaBXUMprh",
        "qCr@XIJtA@",
        "KAx@@IWMmUPAlHPfES]FNa`@",
        "HeT@@DjYYZPbJ@@[JDIaTwCH",
        "sNp@DkfZj@[\\QXuqea`@",
        "Hid@@DjU^nBAHAEVtaBXUMp`",
        "Hmt@@DjYU^Vjjj`AhtISRmmBDpj[aeP",
        "sGP@DkejPFsDVM@",
        "sNx@@eJmUPB",
        "qCb@AIf`H",
        "HcL@@DjU_VnBBJh@FqaBXUMprnqcXX@",
        "Hid@@DjUZnBBH@FtaBXUMpqahX@",
        "sNp@LdbQZjBBvxbqkcGC@",
        "sOx@@drm\\@@@|c@",
        "sJP@H~j@^R@",
        "KAx@@YIDcFjhDElHPfES\\",
        "Hid@@DjUZnBAH@FtaBXUMp`",
        "sNp@LddUji@[\\QXu`",
        "sGP@DjfjPFsDVM@",
        "HeT@@DjYUXPbD`@[JDIaTwCH",
        "KAx@@IUoMUPAlHPfES]FDa`@",
        "sFp@DiTt@@AyaD",
        "Hed@@DjuuZjj@FraBXUMpsFIa`@",
        "HeT@@DjUghP`h`@[JDIaTwCLXfF@",
        "sOp@DjWkjj`FwDVM\\YhX@",
        "sGP@Djfj`FsDVMcIC@",
        "KAx@@IRkmUPHKXPaLJfzL]C@",
        "sNx@@djmUPB",
        "QM`AIdD",
        "sOp@DjWkB@@Gbe@",
        "sNp@DjyZj@[\\QXuqca`@",
        "QM@HuAmd`",
        "sNp@LddUjj@[\\QXuqea`@",
        "HaD@@DkeVyjj`AhrXUMuaBDpj[hpDL@",
        "qCb@AIZPH",
        "HiD@@LdbJTjjh@[RDIaTwCF@a`@",
        "Hmt@@DjU_ZxHHi@AmhPfES\\Lj",
        "HaDH@@RYWih@H@D",
        "HiD@@LdbJTjjh@[RDIaTwCFHa`@",
        "sGX@@djuT@|a@",
        "sNp@DkfZj@[\\QXuqaa`@",
        "Hid@@DjU^nBBH@GdL",
        "KAx@@IVkMUPAlHPfES]FJa`@",
        "qCr@XIKTAy@",
        "HmT@@Dj{uVjjh@[ZDIaTwCJqaXX@",
        "Hmt@@DjYWVFjjj`AhpQe\\mmBDpj[aeP",
        "Hif@@@RUe^Fh@@@P",
        "HaDH@@Rfu[j@@@GdH",
        "KAx@@IVsMUPAlHPfES]FDa`@",
        "sKP@Di\\Zj@[TQXq`a`@",
        "sJX@@eMU@OH@",
        "HeT@@DjU^k``b`@[JDIaTwCLXFF@",
        "Hmt@@LdbbRJXPbHh@FvaBXUMprh",
        "sJP@DjvhAmQEcFBF@",
        "Hmt@@LdbbRNXZjjj@FcAFUrvtHSBinFUcBpp@",
        "sJP`@dfvhAyD",
        "sGP@Di^V`dmfHl[FVF@",
        "KAx@@IVsmUPAlHPfES]FBa`@",
        "sOq@@drm\\@@@|PP",
        "sJY@BDeZhA@",
        "HeT@@LdbRbmBDED@CYPaLJfxY@",
        "Hed@@Djy[Zjj@FraBXUMpr",
        "HeT@@DjU]k``b`@[JDIaTwCLXFF@",
        "Hid@@DjUfaBB`@D",
        "qCa@CIJtA@",
        "QMPARVA@",
        "Hid@@DjUfaBB`@FtaBXUMpqcHX@",
        "sJY@BDfZhA@",
        "HeT@@DjUghP`hP@[JDIaTwCH",
        "Hed@@Dj{uZjj@FraBXUMpsFIa`@",
        "Hmt@@LdbbRUXZjjj@FcAFUrvtHSBinFUcFPp@",
        "sNp`@dfuZj@P",
        "sJQ@@dmU@OH@",
        "sJX@@dmU@H",
        "HeT@@DjU]k``b`@[JDIaTwCLXZF@",
        "HiD@@LdfbJZjh@[RDIaTwCFAa`@",
        "sOx@@drm\\@@@|a@",
        "HeT@@LdbbQgCUUU@CQhRfz[JDIaTwCH",
        "Hmt@@DjU]Zzjjj`AhpQaLmmBDpj[aeXplL@",
        "sOp@DjWkjj`FwDVM\\XHX@",
        "HcL@@LdbbRNSBDQEP@McBDpj[ae]cFpp@",
        "HiD@@Dj}Yji`AmHPfES\\H",
        "HaDH@@RYe[hB@@D",
        "Hid@@DjU^njjj@FtaBXUMpq`XX@",
        "HeT@@DkYeFVjjh@ZMaUpsYPaLJfxY@",
        "QMPARZA@",
        "sOq@@drm\\@@@|QX",
        "HaD@@DjYvxH`@A@",
        "HcL@@LdbbRNcBDQEP@McBDpj[ae]@",
        "QMhDRZA@",
        "RG@DXLHmP",
        "QM`BN`XQYd",
        "RG@DTLHmP",
        "QMHAIXFEVd",
        "QMDARVAaH",
        "RFPDXLHmP",
        "RF@Q`vRbdLEC@",
        "RF@QpvR@",
        "QO@HyjAmd`",
        "`II@B",
        "`II@CFspqJp",
        "`II@CF[@hM@prB`",
        "`H@[T[|B`XN@PdM@p|@bHrBcDk@",
		"RG@DXMj}F@",
		"QM`BN`[L~b@",
		"RG@DTMj}D@",
		"QMHAIXFt~j@",
		"QMDARVA}L@",
		"RFPDXMj}D@",
		"sKP@Di\\YZ@[TQXqaa`@",
		"RG@DXMH" };

	private static StereoMolecule[]  sKeyFragment;
	private SSSearcher			mSSSearcher;
	private StereoMolecule		mMolecule,mFragment;
	private int[]				mMoleculeIndex,mFragmentIndex;
	private byte[]				mMoleculeIDCode,mFragmentIDCode;

	public static int getNoOfKeys() {
        return cKeyIDCode.length;
	   	}

	public SSSearcherWithIndex() {
		mSSSearcher = new SSSearcher();
		init();
		}


	public SSSearcherWithIndex(int matchMode) {
		mSSSearcher = new SSSearcher(matchMode);
		init();
		}


	public StereoMolecule getKeyFragment(int no) {
		return sKeyFragment[no];
		}


	public void setFragment(StereoMolecule fragment, int[] index) {
		mFragmentIDCode = null;
		mFragment = fragment;
		if (index == null)
			mFragmentIndex = createIndex(fragment);
		else
			mFragmentIndex = index;
		}


	public void setFragment(String idcode, int[] index) {
	    setFragment(idcode.getBytes(), index);
		}


	public void setFragment(byte[] idcode, int[] index) {
		mFragmentIDCode = idcode;
		if (index == null) {
			mFragment = (new IDCodeParser(false)).getCompactMolecule(idcode);
			mFragmentIndex = createIndex(mFragment);
			}
		else {
			mFragment = null;
			mFragmentIndex = index;
			}
		}


	public void setMolecule(StereoMolecule molecule, int[] index) {
		mMoleculeIDCode = null;
		mMolecule = molecule;
		if (index == null)
			mMoleculeIndex = createIndex(molecule);
		else
			mMoleculeIndex = index;
		}


	public void setMolecule(String idcode, int[] index) {
	    setMolecule(idcode.getBytes(), index);
		}


	public void setMolecule(byte[] idcode, int[] index) {
		mMoleculeIDCode = idcode;
		if (index == null) {
			mMolecule = (new IDCodeParser(false)).getCompactMolecule(idcode);
			mMoleculeIndex = createIndex(mMolecule);
			}
		else {
			mMolecule = null;
			mMoleculeIndex = index;
			}
		}


	public int getFirstHittingIndexBlockNo() {
		for (int i=0; i<mMoleculeIndex.length; i++)
			if ((mFragmentIndex[i] & ~mMoleculeIndex[i]) != 0)
				return i;
        return -1;
        }


	/**
	 * Returns the most recently defined molecule. If the molecule was
	 * passed as idcode and was not created yet, then a molecule is
	 * constructed from the idcode and returned.
	 * @return current Molecule (null, if setMolecule() has never been called)
	 */
	public StereoMolecule getMolecule() {
        if (mMolecule == null && mMoleculeIDCode != null)
            mMolecule = (new IDCodeParser(false)).getCompactMolecule(mMoleculeIDCode);

        return mMolecule;
	    }


	public boolean isFragmentInMolecule() {
		for (int i=0; i<mMoleculeIndex.length; i++)
			if ((mFragmentIndex[i] & ~mMoleculeIndex[i]) != 0)
				return false;

		if (mMolecule == null)
			mMolecule = (new IDCodeParser(false)).getCompactMolecule(mMoleculeIDCode);
		if (mFragment == null)
			mFragment = (new IDCodeParser(false)).getCompactMolecule(mFragmentIDCode);

		mSSSearcher.setMolecule(mMolecule);
		mSSSearcher.setFragment(mFragment);
		return mSSSearcher.isFragmentInMolecule();
		}


	public int findFragmentInMolecule() {
		for (int i=0; i<mMoleculeIndex.length; i++)
			if ((mFragmentIndex[i] & ~mMoleculeIndex[i]) != 0)
				return 0;

		if (mMolecule == null)
			mMolecule = (new IDCodeParser(false)).getCompactMolecule(mMoleculeIDCode);
		if (mFragment == null)
			mFragment = (new IDCodeParser(false)).getCompactMolecule(mFragmentIDCode);

		mSSSearcher.setMolecule(mMolecule);
		mSSSearcher.setFragment(mFragment);
		return mSSSearcher.findFragmentInMolecule();
		}


	public int findFragmentInMolecule(int countMode, int matchMode) {
	    return findFragmentInMolecule(countMode, matchMode, null);
	    }

	public int findFragmentInMolecule(int countMode, int matchMode, final boolean[] atomExcluded) {
		for (int i=0; i<mMoleculeIndex.length; i++)
			if ((mFragmentIndex[i] & ~mMoleculeIndex[i]) != 0)
				return 0;

		if (mMolecule == null)
			mMolecule = (new IDCodeParser(false)).getCompactMolecule(mMoleculeIDCode);
		if (mFragment == null)
			mFragment = (new IDCodeParser(false)).getCompactMolecule(mFragmentIDCode);

		mSSSearcher.setMolecule(mMolecule);
		mSSSearcher.setFragment(mFragment);
		return mSSSearcher.findFragmentInMolecule(countMode, matchMode, atomExcluded);
		}


	public ArrayList<int[]> getMatchList() {
		return mSSSearcher.getMatchList();
		}


	public int[] createIndex(StereoMolecule mol) {
		int[] index = new int[(cKeyIDCode.length+31)/32];
		mSSSearcher.setMolecule(mol);
		for (int i=0; i<cKeyIDCode.length; i++) {
			mSSSearcher.setFragment(sKeyFragment[i]);
			if (mSSSearcher.isFragmentInMolecule(SSSearcher.cIndexMatchMode))
				index[i/32] |= (1 << (31-i%32));
			}

		return index;
		}


    public static float getSimilarityTanimoto(int[] index1, int[] index2) {
        int sharedKeys = 0;
        int allKeys = 0;
        for (int i=0; i<index1.length; i++) {
            sharedKeys += bitCount(index1[i] & index2[i]);
            allKeys += bitCount(index1[i] | index2[i]);
        	}
        return (float)sharedKeys/(float)allKeys;
        }


    public static float getSimilarityAngleCosine(int[] index1, int[] index2) {
        int sharedKeys = 0;
        int index1Keys = 0;
        int index2Keys = 0;
        for (int i=0; i<index1.length; i++) {
            sharedKeys += bitCount(index1[i] & index2[i]);
            index1Keys += bitCount(index1[i]);
            index2Keys += bitCount(index2[i]);
        	}
        return (float)sharedKeys/(float)Math.sqrt(index1Keys*index2Keys);
    	}


	public static int[] getIndexFromHexString(String hex) {
		if (hex.length() == 0 || (hex.length() & 7) != 0)
			return null;

		int[] index = new int[hex.length()/8];
		for (int i=0; i<hex.length(); i++) {
			int j = i/8;
			int code = hex.charAt(i) - '0';
			if (code > 16)
				code -= 7;
			index[j] <<= 4;
			index[j] += code;
			}

		return index;
		}


	public static int[] getIndexFromHexString(byte[] bytes) {
		if (bytes==null || bytes.length==0 || (bytes.length & 7) != 0)
			return null;

		int[] index = new int[bytes.length/8];
		for (int i=0; i<bytes.length; i++) {
			int j = i/8;
			int code = (int)bytes[i] - '0';
			if (code > 16)
				code -= 7;
			index[j] <<= 4;
			index[j] += code;
			}

		return index;
		}


	public static String getHexStringFromIndex(int[] index) {
	    if (index == null)
	        return null;

	    byte[] bytes = new byte[index.length*8];
		for (int i=0; i<index.length; i++) {
			int value = index[i];
			for (int j=7; j>=0; j--) {
				int code = value & 15;
				if (code > 9)
					code += 7;
				bytes[i*8+j] = (byte)('0'+code);
				value >>= 4;
				}
			}

		return new String(bytes);
		}

	
    public static int bitCount(int x) {
        int temp;

        temp = 0x55555555;
        x = (x & temp) + (x >>> 1 & temp);
        temp = 0x33333333;
        x = (x & temp) + (x >>> 2 & temp);
        temp = 0x07070707;
        x = (x & temp) + (x >>> 4 & temp);
        temp = 0x000F000F;
        x = (x & temp) + (x >>> 8 & temp);

        return (x & 0x1F) + (x >>> 16);
    	}


	private void init() {
		synchronized(SSSearcherWithIndex.class) {
		    if (sKeyFragment == null) {
	    		IDCodeParser theParser = new IDCodeParser(false);
	    		sKeyFragment = new StereoMolecule[cKeyIDCode.length];
	    		for (int i=0; i<cKeyIDCode.length; i++) {
	    			sKeyFragment[i] = theParser.getCompactMolecule(cKeyIDCode[i]);
	    			sKeyFragment[i].ensureHelperArrays(Molecule.cHelperNeighbours);
	    			}
	    		}
			}
	    }
	}
