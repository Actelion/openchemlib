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
import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.IReactionMapper;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.util.DoubleFormat;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * This class is openchemlib's most advanced reaction mapper and should be used whenever a reaction's atoms need to be mapped,
 * which means when an automatic procedure shall associate reactant atoms to their corresponding product atoms.
 * This is done by computationally mapping reactant atoms to products atoms based on similar neighbourhood environments
 * and by applying the knowledge of common chemical transformation patterns. Individual mapping solutions are scored
 * by adding penalty scores for bonds that need to be broken, changed, or added in order to create the product from the reactant
 * with the given atom mapping.<br>
 * This ChemicalRuleEnhancedReactionMapper uses a knowledge base of common organic chemistry transformations. Each of these
 * transformation rules is defined as a reactant and a product substructure including proper atom mapping.
 * These chemical rule reactions must neither be stoichiometrically complete, nor must they be completely mapped.
 * However, all atoms that exist on both sides of the reaction need to have the same mapping number. Atoms that exist on
 * one side of the reaction only, like disappearing halogenide or added water, must not and, in fact, cannot be mapped.
 * Being substructures, the rule's reactant and products may carry query features and exclude groups. In fact, query
 * features and exclude groups should be used whereever possible to make the reactant substructure specific enough to
 * match all intended reactions, but not many others. For instance, an amine for a reductive amination should be defined
 * as carrying two hydrogen atoms. Moreover, an exclude group should be used to prevent attached carbonyl- and similar groups.<br>
 * When a reaction is given to this class to be mapped, then every rule's reactant and product substructures are searches
 * in the corresponding sides of the reaction. For every match the rule is applied, mapped by the SimilarityGraphBasedReactionMapper
 * and scored. The unmodified reaction is mapped and scored as well. The best scoring mapping is finally added to the given reaction.
 */
public class ChemicalRuleEnhancedReactionMapper implements IReactionMapper {
	private static final String RULES_DATAWARRIOR_FILE = "/home/thomas/data/reactionmapping/ChemicalRuleEnhancedReactionMapper_allRules.dwar";
	private static final int MAX_MATCH_COUNT = 512;  // Protection for combinatorial explosion, e.g. for metathesis or DielsAlder in fullerene

	private static final boolean DEBUG_PRINT_REACTION_AFTER_APPLYING_RULE = false;
	private static final boolean DEBUG_PRINT_MOLFILES_AFTER_APPLYING_RULE = false;

	// Chemical rule reactions must neither be stoichiometrically complete, nor must they be completely mapped!!!
	// If rules contain exclude atoms, these must not be mapped.
	// When a rule is applied to a reaction, then the reaction is checked, whether the rule's reactant and product
	// can be found as substructures within the respective components of the reaction. If yes, then the mapped region of
	// the rule is used as a template to modify the reactant's bonds according to the rule's transformation.
	// After the change, the modified reaction is automatically mapped applying the SimilarityGraphBasedReactionMapper
	// and then scored. If a rule matches what actually happens in the reaction, then its application will modify the
	// reactant into a structure being much more similar to the product. Thus, a mapping can be found that requires
	// no or very few bond changes to create the product from the modified reactant, resulting in a high score.
	// NOTE: You may use the idcodeexplorer.jar to generate new rules. An up-to-date version should always be available
	//       from https://openmolecules.org/idcodeexplorer.jar
	// NOTE: Rule reactions must contain one reactant and one product only. If you have multiple reactants (or products),
	//       then reduce the distance between individual reactants such that the idcodeexplorer recognizes them all
	//       as one logical product indicated by one large 'A' in the background.
	// NOTE: OpenChemLib distinguishes between manual (red) and automatic (green) mapping numbers. Rules should contain
	//       exclusively manual mapping numbers (use popup menu in mapping mode to change).
	private static final String[][] CHEMICAL_RULE = {
			{"e", "daXJ@@PjyITuiX@`!dahJ@@SJYIMMfPB#IaLJfxP IaAhl[`#!BDpAl@AL@[@Be}aL@[@@ !B|Osp?QZR?O_}}zGze`@##"},
// bad		{"f", "gJQ`@bdjt`P!gKa`@ldfrA@#qbqh qqlP#!R_zq?dw~l_yLsXgp"},
			{"g", "gBa`@lde@XS@!gCA`@mDPOdku`#qbq qJt#!B[G|S_qgq !BjW}q]cga##"},
// no evidence	{"h", "daXL@@Klb`JSSHA@!daXL@@[dfTBijh@`#IyjFDp` IyjEL``#!R?g~HO[\\}[Lfw?K_}m?psLap"},
// evidences covered by aldol		{"aldolEnoneCyclisation", "daXD@@kHh`Ttl@P!daxD@@yIeRfZZ`B#IqBbay` IqHKjXP#!R_g}wwWC]]xKCfXbBbMtshep"},
// fishy    	{"j", "gBi@DDcXCnR`!gBi@HDFXCnY`#qbq qfQ#!R_vpy`W}lLvK|"},
// no evidence	{"l", "daxH@@RUUjjPB!daDH@@RVU[jij@H#IGfhlR@ IGfbdK@#!R@IL@k@BS@Jp@dpBl@ILs|kp"},
// no evidence	{"m", "gFQHBGAISiPH!gGQHJGAIUjhB#qNT] qNTk#!R@AL@[@@Sa_x@DsNro@"},
// no evidence	{"n", "gOp`@|dRaqij`H!gOp`@tiUMij`H#qosJ` qgSqh#!RTv]`YRqg?g~HPh`}L{H|"},
			{"o", "gGQHDDqIUjdB!gGQHHDqAeMjhB#qbqk qfQ]#!B@AL@[@@S@Fp@Dp !BFQ{~_q|ZGvUSYp##"},
// no evidence	{"p", "gGQ@@eJuRA@!gFQ@@diuPD#qqUc qrcM#!R_zp@kG~S@IM?kCLb?@"},
			{"r", "gOQdEjHbAFQRBrtAaJ!gNQdEbHbABCTKPFDH#qbMwX qbM~X#!BNDm`ra?UjW~YBYX@ !Ba[zw?_x@?g~H?XO~##"},
			{"s", "gOp`ATigujj`H!gOp`ATiVKjj`H#qnyS` qY~eP#!BTLtk^sE?BOs|]pc} !BbOvw?_y??g~H?[_}##"},
			{"t", "gOP`@dcUZdB!gNp`@tiTMjj@`#q{ir` qkLrx#!Be`Bzr_wp?OC}|Osp !B?g~w@k_}m?vw@n[a##"},
			{"u", "daXB`Hrn@HrPEILt`D!daXB`DJn@HRUMjV@H#IxVaLJ` IylJhPP#!B[G}l@OKyDpAl@AL@[@@ !B@Fp@DpAl@AN]?`Al@AL##"},

			{"Sakurai", "gOQH@wAINvZ@pdcFe@x@!gOQH@wAIgJi@pdcFe@x@#qreKx qrLkx#!BDpAl@IknDw|S@Fp@ !Bb@JH?_x@b@JH?Ven##"},
			{"Mitsunobu", "gFP`ATfRhdPp`}KEYg]d@!gFP`ATfRhd`pekL{l`#qrLk qZLn#!B@hc}b@C~@h`YM` !B@hc}b@C~@h`YM`##"},

			{"Aldol-Addition", "gOQ@AdTAcS@[Q^crJTLES`DJsL?vH!gGQ@@dsuRAcJg@HUaX#qYEbp qYub#!BxOyBzLKg`dG~xG~{ !Bb@K~@Hc}FBIA@@##"},
			{"Aldol-Condensation", "gOQ@AdTAcS@[Q^crJTLES`DJsL?vH!gFQ@@`rrpdlHHpipBEXf@#qYEbp q^aU#!B{ZRRqA?AQfyA@L_C !B}QFw@h`B_tnH_P##"},
			{"Acetal-Aldol-Addition", "dmdB@@serQS@sJjfd@p`Xir\\@`j\\@aUJXK@!daxL@@[df[ZjT@qQdxACdxABjTqf@#qB@`OuX qBtM{#!B_]]}rHKWw?y}uy[~GJbBu{wWqG| !BfJK~TIa]fJJghg{`pP@##"},
			{"Acetal-Aldol-Condensation", "dmdB@@serQS@sJjfd@p`Xir\\@`j\\@aUJXK@!gNp`CTjUiV@qQS`DJg@HUVXV@#qB@`OuX qqj{`#!B?[_}b@Jw_?{}m~[~[N@Bm?vwkN| !BfsK~yzPrw}m?rzQM##"},
			{"Acetal-Aldol-Condensation-Cyclization", "dkLB@@ZURYUvUjljh@paHpr\\@`!didD@@yIfUXXHL@CFNS`D@#IXljNPY@@@ IXljIyA#!BbOw~_x`Bm?vH?[_}b@JH?_y?b@Jw?Xc} !BbOvH?Oy??`BH?Xa}?`C~_p##"},
			{"Enolester-Cleavage", "gOQ`@fdscT`_Qp!gOQ`@cdTASS@P#q}Kr` q}cNP#!B@k]}mpC~@k]}mqdY !Bb@K~@Hc}BfzH@hc}##"},

			{"propargylEnone", "gCa@@dmXFD@!gCa@@dkHD#qNT qLV#!BbOw~_?y? !B@AL@[@@S##"},
			{"Arndt-Eistert", "daiDaJYBBHj^{HhAYMpAaA@!daiD`FzLBHPVsZl@p`@#IyHHZ[@ IzDGBi`#!B?`Bw?H`Bn]{~\\g?~@Ox !B@rzH?_y?b@JH?_n~bOt##"},
			{"Curtius", "gO]IcVaDF[s{HhCIe@`!gN]HMWADHJfm`XP@#q~Jk` qytUX#!B@O{|b@Jw\\o{~_?x@ !Bj]y??g?~?[^G_hc}##"},
			{"diazomethanHomologation", "gFU@lQioIIs\\AyH!gFU@CPdimXD#qbM^ qbqk#!B?X`BbFZN?`C~_p !B@AL@[@@Su[x@Dp##"},

			// methathese
			{"ene-Metathesis", "deD@@LdbEdSTu@FqHWSDda`JFChXIa`?tdKi@!deD@@Ldb`\\SUM@FqHPsDda`JF@XXIa`?tdHY@#qTEApX qQECf@#!BQzK~}ubbW`BEgcb]?a@gg[zO !BQzK~}ubbW`Ag{VVAQzJ~c?xP##"},
			{"ene-Metathesis", "deD@@LdbEdSTu@FqHWSDda`JFChXIa`?tdKi@!gC`@DiZDE@#qPDA@p qQf#!BmpK~_x`Bm?tAs[]}?`BH_[_| !B_vp@[G|S##"},
			{"yne-Metathesis", "daX@@LdPLWWPAlRXwQIHXLa`ZFChXO}IL[`!daX@@LdPLWWPAlRXwQIHXLa`ZFChXO}IL[`#qT@q@ qt@Q@#!BZmoi@Fjo|SFe|IkGiU@ !BZmoi@Fjo|SFe|IkGiU@##"},
			{"ene-yne-Metathesis", "dcd@@LdPLPLWSSPAlrXwKlVRFCHXFa`zFAXXMa`?udqnWP!dcT@@LdbbplTsTt@[MeMr{Ela`jFAhXNa`VFCXXO}[J[et#qe@N@S@ qeHP@s@#!B?g?~@Oy?^gG}bOvw?H`E@PJw@hc}mp !B?`BH?[_}mpJH@oy??`AuC`Jw@hc}mp##"},
			{"Alkyne-Cyclisation", "gG`@DcO|CFAFC`@!gFp@DiTt@@CFAFC`@#qi\\V qiSt#!B_y[qA`Biu^zV@@ !B?g~w@k_}m?vw@`##"},

			// two-step
			{"Elimination-Claisen", "gNp`AldTQji@~a`!gGQ@@djmRA@#qtHUX qtSi#!Bm?vH?[\\B?g~H@hc} !B@AL@[@@S@Fp@Dp##"},
			{"imineFormationAzaCope", "daZH@LAIMUjd@pRL@!daZH@HAAn]jd@p`@#IGfaLJ` IFDzfK@#!BDpAl@IkqDpAl@AL@[@@ !BFaFw@h`BbOw~@H`BbOt##"},
			{"didehydroCopeWithAromatisation", "gNp@Di]ejDcjcr|wK`!gOp@DjWkB@@H#qrLkx q\\oQp#!B?g~H?K_}bGvH?H`B !Bm?w~@Hc}mpJw@ox@##"},

			// multistep with cyclisation/condensation
			{"symAldolNitrogenRing", "dovJ@GBfttf\\v\\qjViPCADGbDodnGp!doNJ@JCSmtefWTCaYjje@H#IlZXi]]yL~C IqMVCzaIim?#!BQtl_riY?Qtl_rfuvNCQ`uZd@NCQ`uVVu}?sA]P !B?`BH@ox@bGvH@k\\Bb@JH_Xa}b@K~_rYltUr|W@##"},

			// pericyclic
			{"Diels-Alder", "gFP@LdPLjA@!gFp@DiTujhB#qiZf qMSf#!B?_C}}?spIPFV@@ !B?g~w@k_}m?vw@`##"},
			{"Claisen-RA", "gGQ@@eKtRA@!gGQ@@djqRA@#qMsT qM\\V#!B_qL@Dw}l@Fp@Dp !B@AL@[@@S@Fp@Dp##"},
			{"Cope", "gGQ@DeZmRAbhcApIF@P@!gGQ@HeZmRAbhc@HIFC`@#qkNT qi\\V#!B_vp@[@@S@Fp@Dp !B_vp@[@@S@Fp@Dp##"},
			{"OxyCope", "gNq@@dr}SHFD@!gNq@@djkUHD#qynZ` qykbp#!B?g~w?Xc}mpJH@hc} !B@Fp@DpAl@AL@[@@S##"},

			// rearrangements
			{"Vinylcyclopropane", "gKP@DklVj@`!gKP@Di\\Vj@`#qaVh qTqh#!Bm?vH?PC~?K\\ !B?g|_Fa}eTv\\##"},
			{"Furan-CH2-OH-Rearrangement", "gOp`@tiguif`H!gOp`@tfWMZZ`H#qZna@ qZtM@#!BTLtk^sE?BOs|]pc} !BBOuPtdy?UGm@V]Ho##"},
			{"rearrangement1032", "gOp`ATieMji`H!gOp`ATeekZj`H#qaSnx qa]~P#!BTLtk^pc|LW?|]pc} !BBOpH?UCRSg}T^tAY##"},

			// 1,2-shifts
			{"Pinacol-Rearrangement", "gNp`ATffjj@pPh!gOP`@tfXZhCAA`#qb^ix qb^oH#!B@k^H@k_}@k_~@Hc} !BbOvH@oy??`BH?PFf##"},
/*TODO check this!*/	{"1,3-WaterEliminationWith12Shift", "gJQ@@dmU@_SP!gKa@@`vtA}A@#qiTp qjap#!BbOvH@ox@bOt !BJdE?[@Al@AL##"},

			// oxidative rearrangements
			{"Epoxydation", "gB``ADcdCB@!gC``AhtUPGtt@#qqb qtQ#!BjW}Y\\YX@ !B?g~w?^Va##"},
/*TODO check this!*/	{"oxydativePropargylAmine13Shift", "gKi@HDEZpLHOQP!gJY@BDeVXB#qMr` qNTh#!BqiXTy{U?mW| !B@Fp@DpAl@AL##"},
			{"Baeyer-Villiger", "gFQ`@[dTAZ`LHP!gFQ`@jdrMPGtl@#qrak qrlK#!B_?{}mwvHs^FVP@ !BbOvH@oy?bOuQzP##"},

			// addition with ring closure
			{"Halogenation ring closure", "gGa@@dYs@XHFJIBDQbUeHPbLRl@!gFQ@@eNUPFJIBDQbUeHPbLRls@`#qbq] qfQk#!B@AOIDW}l@tD@Dp !B_qL@Dw}l_qNcDP##"},
			{"Halogenation ring closure", "gBa@@d\\`XP@!gJQ@@eOU@XpdHQFIVY`P#qbq qfQ@#!B@AOIDW}l !B_qL@Dw}l_qL##"},

			// condensation with ring closure
			{"Hantzsch Thiazol", "daZHPDp@chaIMefh@ppDzTD~hYmC^bhbcPp]dQbUg~pp!gKXHL@aJWFe`H#qNPe@ qNj`#!BvuK[KUxv_yS[k_zhvuH !BTqa`FbpX?`@##"},
			{"Oxadiazole", "gOX`BEdTASW@XQ@!gOu@HPeKNMKTA@#qrDMX qpULX#!BmpK~@K_}Mlx@?`C~ !BZ?`C}v|m_rYR[z?\\##"},
			{"Imidazole", "dmeHPNg@qJqLbTtATijZ@LLJnuDmhWtSDXUFC`?rIoTAP!gOt@ATieuej`H#qDPpM_@ q~ZM`#!BqvKGg_yOqvKGg_xphrGkLcz@_sD !BTMHkACD@BOw|B@QT##"},
			{"1,2,3-Triazole", "gB`HAbIGXFDWiM@!gF|@ADeKXmT`P#QIp Q@v#!B_vpU?g}l !BTv]`YRqg?g|XK@##"},
			{"1,2,4-Triazole", "deFD@NALbbfASUW@FD]YJZLUCAVJ}?nES@!gO|@ABeKNLuRA@#qDB@FM q@LuP#!BY?r~@F_]jDJW`j`}Iaxx[UC] !BTv]@IPqgog|hCBT_##"},

			// ester cleavage; probably not needed, because implicitly correctly handles by SimilarityGraphBasedReactionMapper
//			{"Ester cleavage", "gKa`@bdhtA@!gKa`@ldftA@#qbqh qbnH#!BDw}l_qM?i^d !B?OC}|IfVjW|##"},
	};

	private static ChemicalRule[] sChemicalRule;

	private float mScore;
	private int mMaxRuleTries;
	private ChemicalRule mAppliedRule;
	private StringBuilder mHistory;

	public ChemicalRuleEnhancedReactionMapper() {
		mMaxRuleTries = Integer.MAX_VALUE;
	}

	private static void initialize() {
		if (sChemicalRule == null)
			synchronized (ChemicalRuleEnhancedReactionMapper.class) {
				if (sChemicalRule == null) {
					ChemicalRule[] chemicalRule = new ChemicalRule[CHEMICAL_RULE.length];
					for (int i=0; i<CHEMICAL_RULE.length; i++)
						chemicalRule[i] = new ChemicalRule(CHEMICAL_RULE[i][0], CHEMICAL_RULE[i][1]);
					sChemicalRule = chemicalRule;
				}
			}
		}

	@Override
	public Reaction mapReaction(Reaction rxn, SSSearcher sss) {
		map(rxn);
		return rxn;
		}

	public void map(Reaction rxn) {
		initialize();
		SimilarityGraphBasedReactionMapper mapper = new SimilarityGraphBasedReactionMapper();
		mapper.mergeReactantsAndProducts(rxn);

		StereoMolecule reactant = mapper.getReactant();
		StereoMolecule product = mapper.getProduct();
		reactant.ensureHelperArrays(Molecule.cHelperNeighbours);
		product.ensureHelperArrays(Molecule.cHelperNeighbours);

		// TODO use indexes
		SSSearcher reactantSearcher = new SSSearcher();
		SSSearcher productSearcher = new SSSearcher();
		reactantSearcher.setMolecule(reactant);
		productSearcher.setMolecule(product);

		mScore = Integer.MIN_VALUE;
		int[] bestReactantMapNo = null;
		int[] bestProductMapNo = null;
		int bestGraphMapNoCount = 0;
		mAppliedRule = null;
		int ruleApplicationCount = 0;
		mHistory = new StringBuilder();

if (SimilarityGraphBasedReactionMapper.DEBUG)
 System.out.println("Reaction\tScore");

		StereoMolecule adaptedReactant = new StereoMolecule(); // reusable container

		for (ChemicalRule rule:sChemicalRule) {
			if (ruleApplicationCount++ == mMaxRuleTries)
				break;

			reactantSearcher.setFragment(rule.getReactant());
			reactantSearcher.setFragmentSymmetryConstraints(rule.getReactantAtomSymmetryConstraints());
			if (0 != reactantSearcher.findFragmentInMolecule(SSSearcher.cCountModeUnique, SSSearcher.cDefaultMatchMode)) {
				productSearcher.setFragment(rule.getProduct());
				if (0 != productSearcher.findFragmentInMolecule(SSSearcher.cCountModeFirstMatch, SSSearcher.cDefaultMatchMode)
				 && reactantSearcher.getMatchList().size() <= MAX_MATCH_COUNT) {
float historyScore = -10000;
					int[] productMatch = productSearcher.getMatchList().get(0);
					for (int[] reactantMatch:reactantSearcher.getMatchList()) {
						if (ruleApplicationCount++ >= mMaxRuleTries)
							break;

						// For every given rule we check, whether we can apply the rule and whether we get a mapping score better than anything else:
						// We do for all rule reactant substructure matches on the real reactant:
						// - if there is rule product substructure match on the real product then
						//   - we apply the rule to the real reactant by
						//     - changing all changing rule bonds between mapped reactants
						//     - deleting all bonds between mapped and unmapped reactant rule bonds (bonds to leaving atoms)
						//     - if the rule has leaving AND incoming (unmapped) atoms, then create a veto matrix:
						//       - for every leaving atom list all incoming atoms as forbidden mapping partners
						//         (here we assume that the substructure match on the product is the relevant one to locate the correct
						//          real product atom from the rule's unmapped atom. Of course, this is not necessarily correct.)
						//     - using the similarity based mapper we map the reaction of the modified reactant and untouched product
						//     - we score that atom mapping and finally use the best scoring mapping numbers as solution.

						reactant.copyMolecule(adaptedReactant);
						int[] originalToAdaptedAtom = rule.apply(adaptedReactant, reactantMatch);
						boolean[][] vetoMatrix = rule.createVetoMatrix(reactant.getAtoms(), reactantMatch, product.getAtoms(), productMatch);

						adaptedReactant.ensureHelperArrays(Molecule.cHelperNeighbours);
						int[] adaptedReactantMapNo = new int[adaptedReactant.getAtoms()];
						int[] productMapNo = new int[product.getAtoms()];
//System.out.println(new MolfileCreator(reactant).getMolfile());
						mapper.map(adaptedReactant, product, adaptedReactantMapNo, productMapNo, vetoMatrix);
						float score = mapper.getScore() - rule.getPanalty();

if (DEBUG_PRINT_REACTION_AFTER_APPLYING_RULE || DEBUG_PRINT_MOLFILES_AFTER_APPLYING_RULE) {
 if (DEBUG_PRINT_REACTION_AFTER_APPLYING_RULE) {
  Reaction adaptedReaction = new Reaction();
  adaptedReaction.addReactant(adaptedReactant);
  adaptedReaction.addProduct(product);
  System.out.println("Rule " + rule.getName() + " score:" + DoubleFormat.toString(score) + " applied: " + ReactionEncoder.encode(adaptedReaction, false, ReactionEncoder.INCLUDE_DEFAULT));
 }
 if (DEBUG_PRINT_MOLFILES_AFTER_APPLYING_RULE) {
  System.out.println(new MolfileCreator(adaptedReactant).getMolfile());
  System.out.println(new MolfileCreator(product).getMolfile());
 }
}

if (historyScore < score) historyScore = score;
						if (mScore < score) {
							mScore = score;
							bestReactantMapNo = mapNoToOriginal(adaptedReactantMapNo, productMapNo, originalToAdaptedAtom, reactant.getAtoms());
							bestProductMapNo = productMapNo;
							bestGraphMapNoCount = mapper.getGraphMapNoCount();
							mAppliedRule = rule;
							}
						}
String pairSequences = mapper.getAtomPairSequenceCount() <= 1 ? "" : " (rootPairSets:"+mapper.getAtomPairSequenceCount()+")";
mHistory.append(rule.getName()+ DoubleFormat.toString(historyScore)+pairSequences+"\n");
					}
				}
			}

		// map and score the reaction without applying any rules
		int[] reactantMapNo = new int[reactant.getAtoms()];
		int[] productMapNo = new int[product.getAtoms()];
		mapper.map(reactant, product, reactantMapNo, productMapNo, null);
		float score = mapper.getScore();

		if (mScore <= score) {
			mAppliedRule = null;
			mScore = score;
			bestReactantMapNo = reactantMapNo;
			bestProductMapNo = productMapNo;
			bestGraphMapNoCount = mapper.getGraphMapNoCount();
			}
String pairSequences = mapper.getAtomPairSequenceCount() <= 1 ? "" : " (rootPairSets:"+mapper.getAtomPairSequenceCount()+")";
mHistory.append("no rule:"+DoubleFormat.toString(score)+pairSequences+"\n");

		if (mScore != Integer.MIN_VALUE)
			mapper.copyMapNosToReaction(rxn, bestReactantMapNo, bestProductMapNo, bestGraphMapNoCount);

if (SimilarityGraphBasedReactionMapper.DEBUG)
 System.out.println("Done; used "+ruleApplicationCount+" of "+mMaxRuleTries+" allowed rule application tries.");
		}

	/**
	 * If the rule removed some atoms in apply(), then this method can be used to translate generated mapping number
	 * of the adapted reactant after applying a rule back to the original reactant's atoms.
	 * During rule application some atoms may be removed and others may come in. Thus, atom indexes change.
	 * If an adapted reactant atom that didn't exist in the original reactant got a mapping number, then
	 * this mapping number is removed from the product to  well.
	 * @param adaptedReactantMapNo similarity based mapping numbers applied to reactant after applying a rule
	 * @param productMapNo similarity based mapping numbers of the product of the reaction
	 * @param originalToAdaptedAtom maps atom indexes from original reactant to adapted reactant after applying a rule
	 * @return array with mapping numbers corresponding to original reactant's atom indexes
	 */
	private int[] mapNoToOriginal(int[] adaptedReactantMapNo, int[] productMapNo, int[] originalToAdaptedAtom, int origAtomCount) {
		boolean[] isUsedMapNo = new boolean[adaptedReactantMapNo.length+1];
		int[] origMapNo = new int[origAtomCount];
		for (int origAtom=0; origAtom<origAtomCount; origAtom++) {
			int adaptedAtom = (originalToAdaptedAtom == null) ? origAtom : originalToAdaptedAtom[origAtom];
			if (adaptedAtom != -1) {
				origMapNo[origAtom] = adaptedReactantMapNo[adaptedAtom];
				isUsedMapNo[origMapNo[origAtom]] = true;
			}
		}

		for (int i=0; i<productMapNo.length; i++)
			if (!isUsedMapNo[productMapNo[i]])
				productMapNo[i] = 0;

		return origMapNo;
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

	public static void main(String args[]) {
		writeRulesDWARFile();
		}

	private static void writeRulesDWARFile() {
		initialize();
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(RULES_DATAWARRIOR_FILE));
			writer.write("<datawarrior-fileinfo>\n" +
					"<version=\"3.3\">\n" +
					"</datawarrior-fileinfo>\n" +
					"<column properties>\n" +
					"<columnName=\"RxnFP\">\n" +
					"<columnProperty=\"parent\tReaction\">\n" +
					"<columnProperty=\"specialType\tRxnFP\">\n" +
					"<columnProperty=\"version\t1.0.0\">\n" +
					"<columnName=\"FragFp\">\n" +
					"<columnProperty=\"parent\tReaction\">\n" +
					"<columnProperty=\"specialType\tFragFp\">\n" +
					"<columnProperty=\"reactionPart\tproducts\">\n" +
					"<columnProperty=\"version\t1.2.1\">\n" +
					"<columnName=\"Reaction\">\n" +
					"<columnProperty=\"specialType\trxncode\">\n" +
					"</column properties>\n" +
					"RxnFP\tFragFp\tName\tReaction\tPenalty\n");
			for (int i = 0; i<CHEMICAL_RULE.length; i++)
				writer.write("\t\t" + sChemicalRule[i].getName() + "\t" + CHEMICAL_RULE[i][1] + "\t" + DoubleFormat.toString(sChemicalRule[i].getPanalty()) + "\n");
			writer.write("<datawarrior properties>\n" +
					"<columnFilter_Table=\"\">\n" +
					"<columnWidth_Table_Name=\"166\">\n" +
					"<columnWidth_Table_Penalty=\"80\">\n" +
					"<columnWidth_Table_Reaction=\"401\">\n" +
					"<detailView=\"height[Data]=0.5;height[Reaction]=0.5\">\n" +
					"<filter0=\"#string#\tName\">\n" +
					"<filter1=\"#reaction#\tReaction\">\n" +
					"<filter2=\"#double#\tPenalty\">\n" +
					"<filterAnimation2=\"state=stopped low2=80% high1=20% time=10\">\n" +
					"<headerLines_Table=\"2\">\n" +
					"<mainSplitting=\"0.71982\">\n" +
					"<mainView=\"Table\">\n" +
					"<mainViewCount=\"1\">\n" +
					"<mainViewDockInfo0=\"root\">\n" +
					"<mainViewName0=\"Table\">\n" +
					"<mainViewType0=\"tableView\">\n" +
					"<rightSplitting=\"0.66326\">\n" +
					"<rowHeight_Table=\"113\">\n" +
					"</datawarrior properties>\n");
			writer.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}
}
