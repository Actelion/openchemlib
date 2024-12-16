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
	// If rules contain exclude atoms, these must not be mapped.
	// Of course, when a rule is applied, then only the mapped region of the rule is used as a template to
	// change bonding of the reaction the rule is applied to. Nevertheless, the rule's entire reactant is used
	// for the substructure search that identifies applicability.
	// NOTE: Rule reactions must contain one reactant and one product only (usually these consist of multiple fragments)!
	// NOTE: You may use the idcodeexplorer.jar to generate new rules.
	private static final ChemicalRule[] CHEMICAL_RULE = {
// With cleaned coordinates:
			new ChemicalRule("d","gGQ@@eKtRA@!gGQ@@djqRA@#qMsT qM\\V#!B_qL@Dw}l@Fp@Dp !B@AL@[@@S@Fp@Dp##", 3.5f),
			new ChemicalRule("e","daXJ@@PjyITuiX@`!dahJ@@SJYIMMfPB#IaLJfxP IaAhl[`#!BDpAl@AL@[@Be}aL@[@@ !B|Osp?QZR?O_}}zGze`@##", 1.5f),
// bad	new ChemicalRule("f","gJQ`@bdjt`P!gKa`@ldfrA@#qbqh qqlP#!R_zq?dw~l_yLsXgp", 0.5f),
			new ChemicalRule("g","gBa`@lde@XS@!gCA`@mDPOdku`#qbq qJt#!B[G|S_qgq !BjW}q]cga##", 5.5f),
// no evidence		new ChemicalRule("h","daXL@@Klb`JSSHA@!daXL@@[dfTBijh@`#IyjFDp` IyjEL``#!R?g~HO[\\}[Lfw?K_}m?psLap", 7.5f),
// evidences covered by aldol		new ChemicalRule("aldolEnoneCyclisation","daXD@@kHh`Ttl@P!daxD@@yIeRfZZ`B#IqBbay` IqHKjXP#!R_g}wwWC]]xKCfXbBbMtshep", 6.5f),
// fishy    	new ChemicalRule("j","gBi@DDcXCnR`!gBi@HDFXCnY`#qbq qfQ#!R_vpy`W}lLvK|", 3.5f),
// no evidence		new ChemicalRule("l","daxH@@RUUjjPB!daDH@@RVU[jij@H#IGfhlR@ IGfbdK@#!R@IL@k@BS@Jp@dpBl@ILs|kp", 0.5f),
// no evidence		new ChemicalRule("m","gFQHBGAISiPH!gGQHJGAIUjhB#qNT] qNTk#!R@AL@[@@Sa_x@DsNro@", 0.5f),
// no evidence		new ChemicalRule("n","gOp`@|dRaqij`H!gOp`@tiUMij`H#qosJ` qgSqh#!RTv]`YRqg?g~HPh`}L{H|", 0.5f),
			new ChemicalRule("o","gGQHDDqIUjdB!gGQHHDqAeMjhB#qbqk qfQ]#!B@AL@[@@S@Fp@Dp !BFQ{~_q|ZGvUSYp##", 0.5f),
// no evidence	new ChemicalRule("p","gGQ@@eJuRA@!gFQ@@diuPD#qqUc qrcM#!R_zp@kG~S@IM?kCLb?@", 5.5f),
			new ChemicalRule("r","gOQdEjHbAFQRBrtAaJ!gNQdEbHbABCTKPFDH#qbMwX qbM~X#!BNDm`ra?UjW~YBYX@ !Ba[zw?_x@?g~H?XO~##", 0.5f),
			new ChemicalRule("s","gOp`ATigujj`H!gOp`ATiVKjj`H#qnyS` qY~eP#!BTLtk^sE?BOs|]pc} !BbOvw?_y??g~H?[_}##", 0.5f),
			new ChemicalRule("t","gOP`@dcUZdB!gNp`@tiTMjj@`#q{ir` qkLrx#!Be`Bzr_wp?OC}|Osp !B?g~w@k_}m?vw@n[a##", 0.5f),
			new ChemicalRule("u","daXB`Hrn@HrPEILt`D!daXB`DJn@HRUMjV@H#IxVaLJ` IylJhPP#!B[G}l@OKyDpAl@AL@[@@ !B@Fp@DpAl@AN]?`Al@AL##", 0.5f),

			new ChemicalRule("Sakurai", "gOQH@wAINvZ@pdcFe@x@!gOQH@wAIgJi@pdcFe@x@#qreKx qrLkx#!BDpAl@IknDw|S@Fp@ !Bb@JH?_x@b@JH?Ven##", 4.5f),
			new ChemicalRule("Mitsunobu", "gFP`ATfRhdPp`}KEYg]d@!gFP`ATfRhd`pekL{l`#qrLk qZLn#!B@hc}b@C~@h`YM` !B@hc}b@C~@h`YM`##", 4.5f),

			new ChemicalRule("Aldol-Addition", "gOQ@AdTAcS@[Q^crJTLES`DJsL?vH!gGQ@@dsuRAcJg@HUaX#qYEbp qYub#!BxOyBzLKg`dG~xG~{ !Bb@K~@Hc}FBIA@@##", 1.5f),
			new ChemicalRule("Aldol-Condensation", "gOQ@AdTAcS@[Q^crJTLES`DJsL?vH!gFQ@@`rrpdlHHpipBEXf@#qYEbp q^aU#!B{ZRRqA?AQfyA@L_C !B}QFw@h`B_tnH_P##", 2.5f),
			new ChemicalRule("Acetal-Aldol-Addition", "dmdB@@serQS@sJjfd@p`Xir\\@`j\\@aUJXK@!daxL@@[df[ZjT@qQdxACdxABjTqf@#qB@`OuX qBtM{#!B_]]}rHKWw?y}uy[~GJbBu{wWqG| !BfJK~TIa]fJJghg{`pP@##", 1.5f),
			new ChemicalRule("Acetal-Aldol-Condensation", "dmdB@@serQS@sJjfd@p`Xir\\@`j\\@aUJXK@!gNp`CTjUiV@qQS`DJg@HUVXV@#qB@`OuX qqj{`#!B?[_}b@Jw_?{}m~[~[N@Bm?vwkN| !BfsK~yzPrw}m?rzQM##", 1.5f),
			new ChemicalRule("Acetal-Aldol-Condensation-Cyclization", "dkLB@@ZURYUvUjljh@paHpr\\@`!didD@@yIfUXXHL@CFNS`D@#IXljNPY@@@ IXljIyA#!BbOw~_x`Bm?vH?[_}b@JH?_y?b@Jw?Xc} !BbOvH?Oy??`BH?Xa}?`C~_p##", 7.5f),
			new ChemicalRule("Enolester-Cleavage", "gOQ`@fdscT`_Qp!gOQ`@cdTASS@P#q}Kr` q}cNP#!B@k]}mpC~@k]}mqdY !Bb@K~@Hc}BfzH@hc}##", 5.3f),

			new ChemicalRule("propargylEnone", "gCa@@dmXFD@!gCa@@dkHD#qNT qLV#!BbOw~_?y? !B@AL@[@@S##", 5.5f),
			new ChemicalRule("Arndt-Eistert", "daiDaJYBBHj^{HhAYMpAaA@!daiD`FzLBHPVsZl@p`@#IyHHZ[@ IzDGBi`#!B?`Bw?H`Bn]{~\\g?~@Ox !B@rzH?_y?b@JH?_n~bOt##", 11.5f),
			new ChemicalRule("Curtius", "gO]IcVaDF[s{HhCIe@`!gN]HMWADHJfm`XP@#q~Jk` qytUX#!B@O{|b@Jw\\o{~_?x@ !Bj]y??g?~?[^G_hc}##", 9.5f),
			new ChemicalRule("diazomethanHomologation", "gFU@lQioIIs\\AyH!gFU@CPdimXD#qbM^ qbqk#!B?X`BbFZN?`C~_p !B@AL@[@@Su[x@Dp##", 7.5f),

			// methathese
			new ChemicalRule("ene-Metathesis","daX@@LdPLSSPAlRXwQIHXLa`ZFChXO}IL[`!daX@@LdPLSSPAlRXwQIHXLa`ZFChXO}IL[`#qT@q@ qt@Q@#!BNIu^@O{wD^EGhkzO?a@ !BNIu^@O{wD^EGhkzO?a@##", 3.5f),
			new ChemicalRule("yne-Metathesis","daX@@LdPLWWPAlRXwQIHXLa`ZFChXO}IL[`!daX@@LdPLWWPAlRXwQIHXLa`ZFChXO}IL[`#qT@q@ qt@Q@#!BZmoi@Fjo|SFe|IkGiU@ !BZmoi@Fjo|SFe|IkGiU@##", 3.5f),
			new ChemicalRule("ene-yne-Metathesis","dcd@@LdPLPLWSSPAlrXwKlVRFCHXFa`zFAXXMa`udqnWP!dcT@@LdbbplTsTt@[MeMr{Ela`jFAhXNa`VFCXXO}[J[et#qe@N@S@ qeHP@s@#!B?g?~@Oy?^gG}bOvw?H`E@PJw@hc}mp !B?`BH?[_}mpJH@oy??`AuC`Jw@hc}mp##", 7.5f),
			new ChemicalRule("Alkyne-Cyclisation", "gG`@DcO|CFAFC`@!gFp@DiTt@@CFAFC`@#qi\\V qiSt#!B_y[qA`Biu^zV@@ !B?g~w@k_}m?vw@`##", 8.5f),

			// two-step
			new ChemicalRule("Elimination-Claisen", "gNp`AldTQji@~a`!gOP`@teTZdCzN@#qtHUX qtSi@#!Bm?vH?[\\B?g~H@hc} !B@AL@[@@S@Fp@DweA##", 4.5f),
			new ChemicalRule("imineFormationAzaCope", "daZH@LAIMUjd@pRL@!daZH@HAAn]jd@p`@#IGfaLJ` IFDzfK@#!BDpAl@IkqDpAl@AL@[@@ !BFaFw@h`BbOw~@H`BbOt##", 8.5f),
			new ChemicalRule("didehydroCopeWithAromatisation", "gNp@Di]ejDcjcr|wK`!gOp@DjWkB@@H#qrLkx q\\oQp#!B?g~H?K_}bGvH?H`B !Bm?w~@Hc}mpJw@ox@##", 4.5f),

			// multi-step with cyclisation/condensation
			new ChemicalRule("symAldolNitrogenRing", "dovJ@GBfttf\\v\\qjViPCADGbDodnGp!doNJ@JCSmtefWTCaYjje@H#IlZXi]]yL~C IqMVCzaIim?#!BQtl_riY?Qtl_rfuvNCQ`uZd@NCQ`uVVu}?sA]P !B?`BH@ox@bGvH@k\\Bb@JH_Xa}b@K~_rYltUr|W@##", 0.5f),

			// pericyclic
			new ChemicalRule("Diels-Alder","gFP@LdPLjA@!gFp@DiTujhB#qiZf qMSf#!B?_C}}?spIPFV@@ !B?g~w@k_}m?vw@`##", 3.5f),
			new ChemicalRule("Cope", "gGQ@DeZmRAbhcApIF@P@!gGQ@HeZmRAbhc@HIFC`@#qkNT qi\\V#!B_vp@[@@S@Fp@Dp !B_vp@[@@S@Fp@Dp##", 5.5f),
			new ChemicalRule("OxyCope", "gNq@@dr}SHFD@!gNq@@djkUHD#qynZ` qykbp#!B?g~w?Xc}mpJH@hc} !B@Fp@DpAl@AL@[@@S##", 4.5f),

			// rearrangements
			new ChemicalRule("Vinylcyclopropane", "gKP@DklVj@`!gKP@Di\\Vj@`#qaVh qTqh#!Bm?vH?PC~?K\\ !B?g|_Fa}eTv\\##", 3.5f),
			new ChemicalRule("Furan-CH2-OH-Rearrangement", "gOp`@tiguif`H!gOp`@tfWMZZ`H#qZna@ qZtM@#!BTLtk^sE?BOs|]pc} !BBOuPtdy?UGm@V]Ho##", 6.5f),
			new ChemicalRule("rearrangement1032", "gOp`ATieMji`H!gOp`ATeekZj`H#qaSnx qa]~P#!BTLtk^pc|LW?|]pc} !BBOpH?UCRSg}T^tAY##", 5.5f),

			// 1,2-shifts
			new ChemicalRule("Pinacol-Rearrangement", "gNp`ATffjj@pPh!gOP`@tfXZhCAA`#qb^ix qb^oH#!B@k^H@k_}@k_~@Hc} !BbOvH@oy??`BH?PFf##", 6.5f),
			new ChemicalRule("1,3-WaterEliminationWith12Shift", "gJQ@@dmU@_SP!gKa@@`vtA}A@#qiTp qjap#!BbOvH@ox@bOt !BJdE?[@Al@AL##", 6.5f),

			// oxidative rearrangements
			new ChemicalRule("Epoxydation", "gB``ADcdCB@!gC``AhtUPGtt@#qqb qtQ#!BjW}Y\\YX@ !B?g~w?^Va##", 6.3f),
			new ChemicalRule("oxydativePropargylAmine13Shift", "gKi@HDEZpLHOQP!gJY@BDeVXB#qMr` qNTh#!BqiXTy{U?mW| !B@Fp@DpAl@AL##", 6.5f),
			new ChemicalRule("Baeyer-Villiger", "gFQ`@[dTAZ`LHP!gFQ`@jdrMPGtl@#qrak qrlK#!B_?{}mwvHs^FVP@ !BbOvH@oy?bOuQzP##", 7.5f),

			// condensation with ring closure
			new ChemicalRule("Hantzsch Thiazol", "gOYDGaDDHRTve`H!gKXHL@aJWFe`H#qB`ip qiV`#!B_vq?Dw}lL{y?[G|S !BTqa`FbpX?`@##", 3.5f),
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
