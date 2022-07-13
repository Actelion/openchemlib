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
 * @author Thomas Sander
 */

package com.actelion.research.chem.coords;

import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;

public class InventorDefaultTemplateList extends ArrayList<InventorTemplate> {
	// Simple annelated rings, bi- or tri-cyclo compounds are typically generated reasonably well
	// by the CoordinateInventor without using templates. More complex highly bridged polycyclic
	// compounds gain from provided templates, which usually are optimized 2D-projections from the
	// 3-dimensional molecule.
	private static final String[] DEFAULT_TEMPLATE = {
			// fullerene
			"gkvt@@@@LddTTTrbTRTRTRRRRRRRRRRRRRrVRrIh\\IAaQxlY@gRHdJCJcRXlv_CfJx|A\\hRHejiLaQjTje^kSjtFcIhvXmVKMjt{lN{Kavy\\^wGjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjh@@vo@HBC@PhLN@bPhtFKCcpDbILaRhtzCIbsX\\nOO`JDbqDjSKdJeJmQjtz}Ahr[LVkMnpz\\nwGj{PBhBdBlBBBjBfBnBaBiBeBmBcBkBgBoB`bhbdblbbbjbfbnbabibebmbcbkbgbob`RhRdRlRbRjRfRnRaRiReRmRcRkRgRoR`rhrdrlrbrjrfrnrarirermrcrkrgror`JhJdJlJbJjJfJnJaJiJeJmJcJkJgJoJ`jhjdjljbjjjfjnjajijej` !BnkjyVwsVr|iQn|Q|goTZWPIJwbudnRkVYBez]siZymNJZUqNFBqZWxS~iCXVU]SeRjwrtSPAjkvXLpBAZauDPzq]PfMlecrMnkv|@\\SFD`m|mWiEoCXp`SIe_J[l|[XCbloTV`[Gc@FJGopyyoOlFQfUy^w\\Bgz|",
			"gcrt@@@@LdbbbbTRbRbRbRRRRRRRRRRRRVRrVQIA`HtRGAaIxZAHfShTjCIbqylQGKgqdBaXeQJeruBiPitZmFoPZLFSYbvZlVGMnsZ]vWSmr{]UUUUUUUUUUUUUUUUUUUUUUUUUUUUUT@@[G`DAA`HTFG@QHTZCEaqxBQDfPiTZ]AdqYlNWGgpEBQXbUIerEReVhuZ]^`tYMfKUfwX]NW[jkPBhBdBlBbBjBfBnBaBiBeBmBcBkBgBoB`bhbdblbbbjbfbnbabibebmbcbkbgbob`RhRdRlRbRjRfRnRaRiReRmRcRkRgRoR`rhrdrlrbrjrfrnrarirermrcrkrgror`JhJdJlJbJjJfJnJaJiJeJmJcJkJgJoJ`jhjdjljbjjjfjnjajij` !B^cR]`]Fm]QkfljE\\p\u007FUVfgOmFXsQe_gXPyXis_wF|vUUX_XbxpzU]HUFgYViwFo~@uemc@}~T\u007FIEPioYVwr]JnM~[ZEC\\g}~o_pUfdo~irsklTLiyVJshnw^iVAsZ`_~}PYkckURH{FYMImFaRaccUlCZSHMfP",

			// adamantane
			"dml@@Dje^VGiyZjjjh@vtHSBinFU@ !BPTCTy[skMzUPF`AJbBixEZHS[Il",
			"dml@@DjYVvGiyZjjjh@vtHSBinFU@ !BwLo~BJ~UquhXBinZ\\ykA@F_eMrT",
			"dml@@LdfbbQX^fUZjjj`C[PaLJfxYT !BzxIHVc{OiJVRpprePho~]}y\u007FwLl",
			"deL@@DjUYkfEijjjj@MeBDpj[ad !B\u007FaA\u007FMVr[AvkKzm_jKvVbD{sk",

			// cubane
			"dil@@LddTQRl[NX^Fjjjj@MiBDpj[a@ !BPfL@\u007Fox@M~T@\u007Fox@\u007F`C~@@",
			"daL@@DjYtKJqjynjjjj@MaBDpj[` !B`bL@_gx@@Gy~@Gx@_`@",
	};

	public InventorDefaultTemplateList() {
		SSSearcherWithIndex searcher = new SSSearcherWithIndex();
		for (String idcode:DEFAULT_TEMPLATE) {
			StereoMolecule fragment = new IDCodeParserWithoutCoordinateInvention().getCompactMolecule(idcode);
			long[] ffp = searcher.createLongIndex(fragment);
			InventorTemplate template = new InventorTemplate(fragment, ffp, false);
			template.normalizeCoordinates();
			add(template);
			}
		}
	}
