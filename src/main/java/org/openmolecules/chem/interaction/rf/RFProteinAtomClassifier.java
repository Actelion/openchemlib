package org.openmolecules.chem.interaction.rf;

import com.actelion.research.chem.SortedStringList;
import com.actelion.research.chem.StereoMolecule;
import org.openmolecules.chem.interaction.AtomClassifier;

public class RFProteinAtomClassifier extends AtomClassifier {
	private static final String[][] PROTEIN_TYPE = {
			{ "gGP@LdbMU@XPZpRu?RHFzj@", "C_ali_apol" },
			{ "eM@HzCB[UFhXOtdVzj{P", "C_ali_apol" },
			{ "gJP@DiZhCCBKRrW}H`[jh", "C_ali_apol" },
			{ "gC`@Dij@ppbuL?iDC]U@", "C_ali_apol" },
			{ "gJP@DjZhCCBGRSg}Hd[jh", "C_ali_apol" },
			{ "gChHA`aIZ`LHmRLYdBDS}Hb[jkl@", "C_ali_apol" },
			{ "eMB@HchLIuS^m?iHmuT", "C_ali_apol" },
			{ "gCi@LDeZ@wb{~dQMuT", "C_ali_don" },
			{ "eM`AIhLJnfw_zRS]U@", "C_ali_don" },
			{ "eM`AIhLDkso}IInj`", "C_ali_don" },
			{ "gCi@DDeZ@zR\\yb[s[?RIfzjcD`", "C_ali_don" },
			{ "eMHAIhLDmR?iIMuT", "C_ali_don" },
			{ "gCa@@duPFBDtdy?RIFzj@", "C_ali_don" },
			{ "gCh@@duTe@p`uDvfcM@a~dRMuT", "C_ali_don" },
			{ "gJ\\@AbeMU@XPURql[dozQDwUP", "C_ali_don" },
			{ "gCi@LDeZ@pb?RHfzj@", "C_ali_don" },
			{ "gJY@LDefjH`XP?iDK]U@", "C_ali_weak_pol" },
			{ "gGY@LDemjhbAaCFUF@H?RHVzj@", "C_ali_weak_pol" },
			{ "gNy@LDenZjbHFDOzQBwUP", "C_ali_weak_pol" },
			{ "eMHAIhLIuS^n?iHmuT", "C_ali_weak_pol" },
			{ "gCa@@dmPFDUXMWd{~dQMuT", "C_ali_weak_pol" },
			{ "gJX@@dkU@XZK_iDS]U@", "C_ali_weak_pol" },
			{ "eMhDRVC^n?iIMuT", "C_pi_carbonyl" },
			{ "eFHBJF}[?RR[jh", "C_pi_carbonyl" },
			{ "eM`AIhLJjfkWzRS]UQRh", "C_pi_don" },
			{ "eMPARZCBjiju~ddwUTTj@", "C_pi_don" },
			{ "eMDARZCrjibqnk~ddwUP", "C_pi_neg" },
			{ "gCh@@duPFjfek}Hd[jjJdkl@", "C_pi_phenyl" },
			{ "eM@HzCfFmZ?RP[jjJM@", "C_pi_phenyl" },
			{ "eM`AIhOJjfkWzRS]UQRh", "C_pi_phenyl" },
			{ "gCl@@lduPFnFiiZwI_tbQnj`", "C_pi_pos" },
			{ "eMhDRVCB^U~dbwUP", "N_amide_don" },
			{ "eMhDRVCrkso}IEnj`", "N_amide_don" },
			{ "gCi@DDeZ@|buLyb[s]?RHfzj@", "N_amide_don_tert" },
			{ "gKT@@Ti\\YZ@pHVOIOtRAnjbwUP", "N_arom_mix" },
			{ "gKT@Adi\\Vf@ptVyC~bPMuUfzj@", "N_arom_mix" },
			{ "diV@@@RfU|kahDB@CA`Pfx^QWzP`MuT", "N_arom_don" },
			{ "gCi@LDeZ@pfoEw}Hf[jh", "N_don_pos" },
			{ "gGX`DJdjmRAyMjYsDwbw~dQMuUuB@", "N_don_pos" },
			{ "eM`AIhLDmRwgozRC]U@", "N_don_pos" },
			{ "eM`AIhLDkso}IAnj`", "N_don_pos" },
			{ "gCl@@lduPFBVkajZVmrV|F?trAnjdwUV[jh", "N_pi_don_pos" },
			{ "eMDARVCjf}[?RQ[jkl@", "O_ali_mix" },
			{ "gCa@@duPFDfdgOzQ@wUP", "O_ali_mix" },
			{ "eMhDRVCjW}IAnj`", "O_pi_acc" },
			{ "eMDARVCB[uo}IAnjnr", "O_pi_acc" },
			{ "eMDARVCB_tTFzjVzj@", "O_pi_acc_neg" },
			{ "eFHBLFDZj?RP[jh", "O_pi_mix" },
			{ "fH`PAo[?Ranj`", "S_apol" },
			{ "fH`PAoW?Ranj`", "S_don" },
			{ "fHbHAoW?Ranj`", "Se_don" },
			{ "fI@FE~eC]U@", "Water" },
			{ "fH@MjL}@p\\@aHNG`DQ@HbXaYV_iPwUP", "Metal" },
	};

	private static volatile RFProteinAtomClassifier sDefaultInstance;
	private static volatile SortedStringList sTypeNameList;
	private static volatile StereoMolecule[] sFragment;
	private static volatile int[][] sFlaggedAtom;

	/**
	 * @return not threadsafe default instance of this class
	 */
	public static RFProteinAtomClassifier getDefaultInstance() {
		if (sDefaultInstance == null)
			sDefaultInstance = new RFProteinAtomClassifier();

		return sDefaultInstance;
	}

	public RFProteinAtomClassifier() {
		if (sTypeNameList == null) {
			synchronized (RFProteinAtomClassifier.class) {
				if (sTypeNameList == null) {
					sTypeNameList = new SortedStringList();
					for (String[] idcodeAndTypeName : PROTEIN_TYPE)
						sTypeNameList.addString(idcodeAndTypeName[1]);

					sFragment = new StereoMolecule[PROTEIN_TYPE.length];
					sFlaggedAtom = new int[PROTEIN_TYPE.length][];

					initializeStaticStuff(PROTEIN_TYPE, sFragment, sFlaggedAtom);
				}
			}
		}
		initialize(sFragment, sFlaggedAtom);
	}

	@Override
	public int getAtomTypeIndex(String typeName) {
		return (typeName == null || typeName.equals("?")) ? 0 : 1+sTypeNameList.getListIndex(typeName);
	}

	@Override
	public String getAtomTypeName(int typeIndex) {
		return typeIndex == TYPE_UNKNOWN ? "?" : sTypeNameList.getStringAt(typeIndex-1);
	}
}
