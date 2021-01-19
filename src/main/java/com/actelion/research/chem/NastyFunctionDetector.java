/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
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

package com.actelion.research.chem;

import java.util.ArrayList;

public class NastyFunctionDetector {
	public static final String[] cNastyFunctionsUnknown = null;

	public static final String sNastyFunction[][] = {
		{ "gGP@LdbKT`]RMdmGJFCPpWN@PBJpd@", "polar activated DB" },
		{ "gGP@DjVePNlFJXhypB@Q\\xA@HjlPd@", "twice activated DB" },
		{ "eM@HvCjFtQ[N@PBDiBHqJrlD@", "acyl-halogenide type" },
		{ "eFBJHcAbc\\@axABIVVQFIV@", "Cl,Br,I on N,O,P,S,Se,I" },
		{ "gC`HADIKRAuHubL", "allyl/benzyl chloride" },
		{ "eM@HzCBKUFTqJp", "prim. alkyl-bromide/iodide" },
		{ "gC`@H}PFrbcLfIV@", "sec./tert. alkyl-bromide/iodide" },
		{ "gC`hH`DIVtAuL`", "alkyl sulfonate/sulfate type" },
		{ "gJQ`@bdjt`P", "anhydride" },
		{ "gJXA@IczhB", "quart. ammonium" },
		{ "gChA@Icm@P", "tert. immonium" },
		{ "fHT`P", "carbenium" },
		{ "gCh`hEIWILtAuM@", "aromatic nitro" },
		{ "gCd@Adeb@p`}M@", "1,2-diamino-aryl" },
		{ "gJT@@TeXHCA@`", "1,3-diamino-aryl" },
		{ "gGT@ATeV@`LDJ", "1,4-diamino-aryl" },
		{ "gCd@ADiZDEsA@", "azo" },
		{ "gJU@h`NdiLsPdh", "azoxy" },
		{ "eMPRIncTH", "diazo" },
		{ "eMPRI^cxH", "diazo" },
		{ "gJT@@Te^lB", "1,1-dinitrile" },
		{ "eMHAIhLDhsW@H^@Pb@", "formaldehyde aduct" },
		{ "eO@HyjCYJLipB@", "oxiran/aziridine" },
		{ "eMPBchLDkR", "hydrazine" },
		{ "eM`AITLYs`D@`", "isocyanate type" },
		{ "fH@MjM~@p\\@aHZA`x^@QDYAQbU`", "unwanted atom" },
		{ "fHw`dB", "phosphonium" },
		{ "gCi@hAteIi`H", "nitrone" },
		{ "eMhHRVCZP", "nitroso" },
		{ "gCa`@lduPD", "orthoester/acid" },
		{ "eFDBcA@", "peroxo" },
		{ "gGY@HDiViPMdmEGN@PBKg@HA@", "N-acyloxy-amide" },
		{ "gC`@H{PFJVQFIV[HcDk@", "1,1-dihalo-alkene" },
		{ "gC`@DiZDEbedQbUfrHqJp", "1,2-dihalo-alkene" },
		{ "fIRPNj@", "pyrylium" },
		{ "gCaHHGAIZPMXb@", "silylenol-ether" },
		{ "gCd@ADie@y``", "dimethylene-hydrazine" },
		{ "eMPARZCJg}T@", "methanediamine" },
		{ "daFD`Bs`BLdTTIUSRpA@", "limit! methylene-thiazolidine-2,4-dione" },
		{ "gOtHLPDISOkSM@XP`", "limit! thiazol-2-ylamine" },
		{ "gGU@DPdsmRAeDx", "acyl-hydrazone" },
		{ "gCh@@eKP`lIIROFIFC`@", "imine/hydrazone of aldehyde" },
		{ "deUD@BxIrJJKPlKTmL@ZqP`", "2,3-diamino-quinone" },
		{ "difL@DBarJIPhfZif@LHh", "limit! 4-acyl-3-azoline-2-one-3-ol" },
		{ "gGT`EPTfyi`H", "limit! oxal-diamide" },
		{ "daED@DpFRYUJfjV@H", "limit! 5-methylene-imidazolidine-2,4-dione" },
		{ "gJQ@@dls@XpyDXeX", "2-halo-enone" },
		{ "gJQ@@djsBJqarHqJp", "3-halo-enone" },
		{ "gCd`i`iJyIf`H", "N-nitro" },
		{ "gChHD@aIf`LYdXN@", "thio-amide/urea" } };

private static StereoMolecule[]	sFragmentList;
	private static long[][]			sIndexList;
	private static boolean			sInitialized;

	public NastyFunctionDetector() {
		synchronized(NastyFunctionDetector.class) {
			if (!sInitialized) {
		        try {
					sFragmentList = new StereoMolecule[sNastyFunction.length];
					sIndexList = new long[sNastyFunction.length][];
					SSSearcherWithIndex sss = new SSSearcherWithIndex(SSSearcher.cMatchAtomCharge);
					for (int i=0; i<sNastyFunction.length; i++) {
						sFragmentList[i] = new IDCodeParser(false).getCompactMolecule(sNastyFunction[i][0]);
						sIndexList[i] = sss.createLongIndex(sFragmentList[i]);
						}
					sInitialized = true;
					}
				catch (Exception e) {
		            System.out.println("Unable to initialize NastyFunctionDetector");
					}
				}
			}
		}

	/**
	 * @param mol
	 * @param index FragFp descriptor
	 * @return list of names of detected nasty function or null in case of initialization error
	 */
	public String[] getNastyFunctionList(StereoMolecule mol, long[] index) {
		if (!sInitialized)
			return cNastyFunctionsUnknown;

		ArrayList<String> nastyFunctionList = new ArrayList<String>();

		SSSearcherWithIndex sss = new SSSearcherWithIndex(SSSearcher.cMatchAtomCharge);
		sss.setMolecule(mol, index);
		for (int i=0; i<sNastyFunction.length; i++) {
			sss.setFragment(sFragmentList[i], sIndexList[i]);
			if (sss.isFragmentInMolecule())
				nastyFunctionList.add(sNastyFunction[i][1]);
			}

		addPolyHaloAromaticRings(mol, nastyFunctionList);

		return nastyFunctionList.toArray(new String[0]);
		}

	/**
	 * @param mol
	 * @param index FragFp descriptor
	 * @return '; ' separated list of detected nasty function names
	 */
	public String getNastyFunctionString(StereoMolecule mol, long[] index) {
		String[] nfl = getNastyFunctionList(mol, index);
		if (nfl == null)
			return "initialization error";
		if (nfl.length == 0)
			return "";
		StringBuilder sb = new StringBuilder(nfl[0]);
		for (int i=1; i<nfl.length; i++)
			sb.append("; "+nfl[i]);
		return sb.toString();
		}

	private void addPolyHaloAromaticRings(StereoMolecule mol, ArrayList<String> nastyFunctionList) {
		RingCollection ringSet = mol.getRingSet();
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			if (ringSet.isAromatic(ring)) {
				int halogenes = 0;
				int[] ringAtom = ringSet.getRingAtoms(ring);
				for (int i=0; i<ringAtom.length; i++) {
					for (int j=0; j<mol.getConnAtoms(ringAtom[i]); j++) {
						int connAtom = mol.getConnAtom(ringAtom[i], j);
						if (!mol.isRingAtom(connAtom)
						 && (mol.getAtomicNo(connAtom) == 9
						  || mol.getAtomicNo(connAtom) == 17
						  || mol.getAtomicNo(connAtom) == 35
						  || mol.getAtomicNo(connAtom) == 53))
							halogenes++;
						}
					}

				if (halogenes > 2) {
					nastyFunctionList.add("polyhalo aromatic ring");
					}
				}
			}
		}
	}
