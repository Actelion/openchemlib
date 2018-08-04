package com.actelion.research.chem.coords;

import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;

public class InventorDefaultTemplateList extends ArrayList<InventorTemplate> {
	// Simple annelated rings, bi- or tri-cyclo compounds are typically generated reasonably well
	// by the CoordinateInventor without using templates. More complex highly bridged polycylcic
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
			InventorTemplate template = new InventorTemplate(fragment, ffp);
			template.normalizeCoordinates();
			add(template);
			}
		}
	}
