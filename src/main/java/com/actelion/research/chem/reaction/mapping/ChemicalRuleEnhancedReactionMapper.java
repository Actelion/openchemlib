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

package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.IReactionMapper;
import com.actelion.research.chem.reaction.Reaction;

public class ChemicalRuleEnhancedReactionMapper implements IReactionMapper {
	private static final int MAX_MATCH_COUNT = 512;  // Protection for combinatorial explosion, e.g. for metathesis or DielsAlder in fullerene

	// Chemical rule reactions must neither be stoichiometrically complete, nor must they be completely mapped!!!
	// If rules contains exclude atoms, these must not be mapped.
	// Of course, when a rule is applied, then only the mapped region of the rule is used as a template to
	// change bonding of the reaction the rule is applied to. Nevertheless, the rule's entire reactant is used
	// for the substructure search that identifies applicability.
	private static final ChemicalRule[] CHEMICAL_RULE = {
		new ChemicalRule("d","gGQ@@eKtrAkDH!gGQ@@djsRIKVPP#qMsT qM\\V#!R_yL@dw~l@Jp@dsNR_@", 3.5f),
		new ChemicalRule("e","daXJ@@PjyITuiX@`!dahJ@@SJYIMMfPB#IaLJfxP IaAhl[`#!ROrp?Ds|lOqNk?g?l_zLsSGp", 1.5f),
// bad	new ChemicalRule("f","gJQ`@bdjt`P!gKa`@ldfrA@#qbqh qqlP#!R_zq?dw~l_yLsXgp", 0.5f),
		new ChemicalRule("g","gBa`@lde@XS@!gCA`@mDPOdku`#qbq qJt#!R_zq?dxgFLvHB", 5.5f),
// no evidence		new ChemicalRule("h","daXL@@Klb`JSSHA@!daXL@@[dfTBijh@`#IyjFDp` IyjEL``#!R?g~HO[\\}[Lfw?K_}m?psLap", 7.5f),
// evidences covered by aldol		new ChemicalRule("aldolEnoneCyclisation","daXD@@kHh`Ttl@P!daxD@@yIeRfZZ`B#IqBbay` IqHKjXP#!R_g}wwWC]]xKCfXbBbMtshep", 6.5f),
// fishy    	new ChemicalRule("j","gBi@DDcXCnR`!gBi@HDFXCnY`#qbq qfQ#!R_vpy`W}lLvK|", 3.5f),
// no evidence		new ChemicalRule("l","daxH@@RUUjjPB!daDH@@RVU[jij@H#IGfhlR@ IGfbdK@#!R@IL@k@BS@Jp@dpBl@ILs|kp", 0.5f),
// no evidence		new ChemicalRule("m","gFQHBGAISiPH!gGQHJGAIUjhB#qNT] qNTk#!R@AL@[@@Sa_x@DsNro@", 0.5f),
// no evidence		new ChemicalRule("n","gOp`@|dRaqij`H!gOp`@tiUMij`H#qosJ` qgSqh#!RTv]`YRqg?g~HPh`}L{H|", 0.5f),
		new ChemicalRule("o","gGQHDDqIUjdB!gGQHHDqAeMjhB#qbqk qfQ]#!R_zq?dw~l_yM?kCM|?@", 0.5f),
// no evidence	new ChemicalRule("p","gGQ@@eJuRA@!gFQ@@diuPD#qqUc qrcM#!R_zp@kG~S@IM?kCLb?@", 5.5f),
		new ChemicalRule("r","gOQdEjHbAFQRBrtAaJ!gNQdEbHbABCTKPFDH#qbMwX qbM~X#!RCwpb@@M`CpL}cg@CL|jB", 0.5f),
		new ChemicalRule("s","gOp`ATigujj`H!gOp`ATiVKjj`H#qnyS` qY~eP#!R?`@_YQ|ZFBqSFHc}L{IB", 0.5f),
		new ChemicalRule("t","gOP`@dcUZdB!gNp`@tiTMjj@`#q{ir` qkLrx#!R@ANZXPAl@AL@[@@SLtj|", 0.5f),
		new ChemicalRule("u","daXB`Hrn@HrPEILt`D!daXB`DJn@HRUMjV@H#IxVaLJ` IylJhPP#!R_zL@hs`Q_zq?dw~l_yLsBkp", 0.5f),

		new ChemicalRule("Sakurai", "gOQH@wAINvZ@pdcFe@x@!gOQH@wAIgJi@pdcFe@x@#qreKx qrLkx#!R_g~HO_fQbOvw?[_|L}r\\", 4.5f),
		new ChemicalRule("Mitsunobu", "gFP`ATfRjdPp`}KEYg]d@!gFP`ATfRjd`pekL{l`#qrLk qZLn#!Rw`Bg?Hc|i}uUYcMb``", 4.5f),

		new ChemicalRule("Aldol-Addition", "gOQ@AdTAcS@^Pvb}GdThXJg@HUfI\u007FlP!gGQ@@dsuRAcJg@HUaH#qYEbp qYub#!Rw[\\\\mw?^akFC|CtwLtI\\", 1.5f),
		new ChemicalRule("Aldol-Condensation", "gOQ@AdTAcS@^Pvb}GdThXJg@HUfI\u007FlP!gFQ@@`rrpdlHHpipBEXb@#qYEbp q^aU#!Rw[\\\\mw?^akFC|CtwLtI\\", 2.5f),
		new ChemicalRule("Acetal-Aldol-Addition", "dmdB@@serQS@sJjfd@p`Xir\\@`j\\@aUJHI@!daxL@@[df[ZjT@qQdxACdxABjTQb@#qB@`OuX qBtM{#!RM?rH?C]}_`CW?Ev^@T@wwS^B_`@sHop", 1.5f),
		new ChemicalRule("Acetal-Aldol-Condensation", "dmdB@@serQS@sJjfd@p`Xir\\@`j\\@aUJHI@!gNp`CTjUiV@qQS`DJg@HUVHR@#qB@`OuX qqj{`#!RM?rH?C]}_`CW?Ev^@T@wwS^B_`@sHop", 1.5f),
		new ChemicalRule("Acetal-Aldol-Condensation-Cyclization", "dkLB@@ZURYUvUjejhHYpaHpr\\@aUJHp`!didD@@EIfU[hBA@CFAS`DJqq@#IXljNPY@@@ IXljXxT#!R_`CW@h`BuwvH_[yOb@I~@M_|bOsW_Wx@LuJb", 7.5f),
		new ChemicalRule("Enolester-Cleavage", "gOQ`@fdscT`_Qp!gOQ`@cdTASS@P#q}Kr` q}cNP#!R?g~H?[_}bOrH?WzeLyH\\", 5.3f),

		new ChemicalRule("propargylEnone", "gCa@@dmXFD@!gCa@@dkHD#qNT qLV#!RXIq`pp@sLwI|", 5.5f),
		new ChemicalRule("Arndt-Eistert", "daiDaJYBBHj^{HhAYMpAaA@!daiD`FzLBHPVsZl@p`@#IyHHZ[@ IzDGBi`#!R@W|h_U\\}X{GUJU\\}TEpsHap", 11.5f),
		new ChemicalRule("Curtius", "gO]IcVaDF[s{HhCIe@`!gN]HMWADHJfm`XP@#q~Jk` qytUX#!R?g}HoU_]U\\eWwQ@\\Lwq\\", 9.5f),
		new ChemicalRule("diazomethanHomologation", "gFU@lQioIIs\\AyH!gFU@CPdimXD#qbM^ qbqk#!Rk}rop?v~k|L@kKNB@`", 7.5f),

		// methathese
		new ChemicalRule("ene-Metathesis","daX@@LdPLSSPHEelRXwQIHXLa`ZFChXO}IL[`!daX@@LdPLSSPHEelRXwQIHXLa`ZFChXO}IL[`#qT@q@ qt@Q@#!RNIu^@O{wD^EGhkzO?aBsdcp", 3.5f),
		new ChemicalRule("yne-Metathesis","daX@@LdPLWWPAlRXwQIHXLa`ZFChXO}IL[`!daX@@LdPLWWPAlRXwQIHXLa`ZFChXO}IL[`#qT@q@ qt@Q@#!RZmoi@Fjo|SFe|IkGiUBSLop", 3.5f),
		new ChemicalRule("ene-yne-Metathesis","dcd@@LdPLPLWSSPIdulrXwKlVRFCHXFa`zFAXXMa`udqnWP!dcT@@LdbbplTsTtFPx}[MeMr{Ela`jFAhXNa`VFCXXO}[J[et#qe@N@S@ qeHP@s@#!R_c}~@Gx?QgF}bKwW@h`yoosW?Hb}usNRO@", 7.5f),
		new ChemicalRule("Alkyne-Cyclisation", "gG`@DcO|CFAFC`@!gFp@DiTt@@CFAFC`@#qi\\V qiSt#!Rb@JNyk\\Bl^{~@CORp`", 8.5f),

		// two step
		new ChemicalRule("Elimination-Claisen", "gNp`AldTQji@~a`!gOP`@teTZdCzN@#qtHUX qtSi@#!RupI~Owx@uwwW_]_|LyK|", 4.5f),
		new ChemicalRule("imineFormationAzaCope", "daZH@LAIMUjd@pRL@!daZH@HAAn]jd@p`@#IGfaLJ` IFDzfK@#!RXpAl@HYrXs}lOvL?[C|sTdH", 8.5f),
		new ChemicalRule("didehydroCopeWithAromatisation", "gNp@Di]ejDcjbrlwK`!gOp@DjWkB@@H#qrLkx q\\oQp#!R?`Bw?[\\BmpK~@K\\BL~JB", 4.5f),

		// multi step with cyclisation/condensation
		new ChemicalRule("symAldolNitrogenRing", "dovJ@GBfttf\\v\\qjViPCADGbDodnGp!doNJ@JCSmtefWTCaYjje@H#IlZXi]]yL~C IqMVCzaIim?#!R@hb}b@A~Owz}uzyl_]\\Bus}~@GxBbLfaOwzUicMbX`", 0.5f),

		// pericyclic
		new ChemicalRule("Diels-Alder","gFP@LdPLjA@!gFp@DiTujhB#qiZf qMSf#!R?`BH?X`BIo[~_sNr``", 3.5f),
		new ChemicalRule("Cope", "gGQ@DeZmRAbhcApIF@P@!gGQ@HeZmRAbhc@HIFC`@#qkNT qi\\V#!R@BM_Hu}lWrM_[COBO@", 5.5f),
		new ChemicalRule("OxyCope", "gNq@@dr}SHFD@!gNq@@djkUHD#qynZ` qykbp#!Ro`AH`c]|\\KtwoS]|LvIB", 4.5f),

		// rearrangements
		new ChemicalRule("Vinylcyclopropane", "gKP@DklVj@`!gKP@Di\\Vj@`#qaVh qTqh#!Rm?t@?h`BbOtsdop", 3.5f),
		new ChemicalRule("Furan-CH2-OH-Rearrangement", "gOp`@tiguif`H!gOp`@tfWMZZ`H#qZna@ qZtM@#!RTXC@z]BRe?s|bKx@L}KB", 6.5f),
		new ChemicalRule("rearrangement1032", "gOp`ATieMji`H!gOp`ATeekZj`H#qaSnx qa]~P#!ROh]`lkoYCONJ_quT|qJl", 5.5f),

		// 1,2-shifts
		new ChemicalRule("Pinacol-Rearrangement", "gNp`ATffjj@pPh!gOP`@tfXZhCAA`#qb^ix qb^oH#!R?m}WoRb}Og?~wu^BLsH\\", 6.5f),
		new ChemicalRule("1,3-WaterEliminationWith12Shift", "gJQ@@dmU@_SP!gKa@@`vtA}A@#qiTp qjap#!R?`ADddRm?basHdH", 6.5f),

		// oxidative rearrangements
		new ChemicalRule("Epoxydation", "gB``ADcdCB@!gC``AhtUPGtt@#qqb qtQ#!R_vsFWg}lLvK|", 6.3f),
		new ChemicalRule("oxydativePropargylAmine13Shift", "gKi@HDEZpLHOQP!gJY@BDeVXQL#qMr` qNTh#!R|Wk@H|@\\@BrStnH", 6.5f),
		new ChemicalRule("Baeyer-Villiger", "gFQ`@[dTAZ`LHP!gFQ`@jdrMPGtl@#qrak qrlK#!R?g~H?[_}AZfw?COBG@", 7.5f),
	};

	private static boolean sInitialized;

	private StereoMolecule mReactant,mProduct;
	private float mScore;
	private int mMaxRuleTries;
	private ChemicalRule mAppliedRule;
	private StringBuilder mHistory;

	public ChemicalRuleEnhancedReactionMapper() {
		if (!sInitialized)
			synchronized (this.getClass()) {
				if (!sInitialized) {
					for (ChemicalRule rule:CHEMICAL_RULE)
						rule.initialize();
					sInitialized = true;
				}
			}
		mMaxRuleTries = Integer.MAX_VALUE;
		}

	@Override
	public Reaction mapReaction(Reaction rxn, SSSearcher sss) {
		map(rxn);
		return rxn;
		}

	public void map(Reaction rxn) {
		SimilarityGraphBasedReactionMapper mapper = new SimilarityGraphBasedReactionMapper();
		mapper.mergeReactantsAndProducts(rxn);

		mReactant = mapper.getReactant();
		mProduct = mapper.getProduct();
		mReactant.ensureHelperArrays(Molecule.cHelperNeighbours);
		mProduct.ensureHelperArrays(Molecule.cHelperNeighbours);

		// TODO use indexes
		SSSearcher reactantSearcher = new SSSearcher();
		SSSearcher productSearcher = new SSSearcher();
		reactantSearcher.setMolecule(mReactant);
		productSearcher.setMolecule(mProduct);

		mScore = Integer.MIN_VALUE;
		int[] bestReactantMapNo = null;
		int[] bestProductMapNo = null;
		int bestGraphMapNoCount = 0;
		mAppliedRule = null;
		int ruleApplicationCount = 0;
		mHistory = new StringBuilder();

if (SimilarityGraphBasedReactionMapper.DEBUG)
 System.out.println("Reaction\tScore");

		StereoMolecule reactant = new StereoMolecule(); // reusable container

		for (ChemicalRule rule:CHEMICAL_RULE) {
			if (ruleApplicationCount++ == mMaxRuleTries)
				break;

			reactantSearcher.setFragment(rule.getReactant());
			reactantSearcher.setFragmentSymmetryConstraints(rule.getReactantAtomSymmetryConstraints());
			if (0 != reactantSearcher.findFragmentInMolecule(SSSearcher.cCountModeUnique, SSSearcher.cDefaultMatchMode)) {
				productSearcher.setFragment(rule.getProduct());
				if (productSearcher.isFragmentInMolecule()
				 && reactantSearcher.getMatchList().size() <= MAX_MATCH_COUNT) {
float historyScore = -10000;
					for (int[] reactantMatch:reactantSearcher.getMatchList()) {
						if (ruleApplicationCount++ >= mMaxRuleTries)
							break;

						mReactant.copyMolecule(reactant);
						rule.apply(reactant, reactantMatch);
						int[] reactantMapNo = new int[mReactant.getAtoms()];
						int[] productMapNo = new int[mProduct.getAtoms()];
//System.out.println(new MolfileCreator(reactant).getMolfile());
						mapper.map(reactant, mProduct, reactantMapNo, productMapNo);
						float score = mapper.getScore() - rule.getPanalty();

if (historyScore < score) historyScore = score;
						if (mScore < score) {
							mScore = score;
							bestReactantMapNo = reactantMapNo;
							bestProductMapNo = productMapNo;
							bestGraphMapNoCount = mapper.getGraphMapNoCount();
							mAppliedRule = rule;
							}
						}
String pairSequences = mapper.getAtomPairSequenceCount() <= 1 ? "" : " (rootPairSets:"+mapper.getAtomPairSequenceCount()+")";
mHistory.append(rule.getName()+historyScore+pairSequences+"\n");
					}
				}
			}

		// map and score the reaction without applying any rules
		int[] reactantMapNo = new int[mReactant.getAtoms()];
		int[] productMapNo = new int[mProduct.getAtoms()];
		mapper.map(mReactant, mProduct, reactantMapNo, productMapNo);
		float score = mapper.getScore();

		if (mScore <= score) {
			mAppliedRule = null;
			mScore = score;
			bestReactantMapNo = reactantMapNo;
			bestProductMapNo = productMapNo;
			bestGraphMapNoCount = mapper.getGraphMapNoCount();
			}
String pairSequences = mapper.getAtomPairSequenceCount() <= 1 ? "" : " (rootPairSets:"+mapper.getAtomPairSequenceCount()+")";
mHistory.append("no rule:"+score+pairSequences+"\n");

		if (mScore != Integer.MIN_VALUE)
			mapper.copyMapNosToReaction(rxn, bestReactantMapNo, bestProductMapNo, bestGraphMapNoCount);

if (SimilarityGraphBasedReactionMapper.DEBUG)
 System.out.println("Done; used "+ruleApplicationCount+" of "+mMaxRuleTries+" allowed rule application tries.");
		}

	/**
	 * This limits the maximum number of times a chemical rule is applied to improve the mapping.
	 * @param max
	 */
	public void setMaximumRuleTries(int max) {
		mMaxRuleTries = max;
		}

	public String getHistory() {
		return mHistory.toString();
	}

	public float getScore() {
		return mScore;
		}

	public ChemicalRule getAppliedRule() {
		return mAppliedRule;
		}
	}
