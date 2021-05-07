package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;

public class ChemicalRuleEnhancedReactionMapper {
	private static final ChemicalRule[] CHEMICAL_RULE = {
		new ChemicalRule("a","gGP@DkUjPLVDXNBHp\\BQ`xLcApEFC`jLG@!gGP@DkUjPLVDXNBHp\\BQ`xLcApEFC`jLG@#qkNT qVci#!Rb@KW@gx@b@JH_SLrP`", 0),
		new ChemicalRule("b","gFP@LdPLjA@!gFp@DiTujhB#qiZf qMSf#!R?`BH?X`BIo[~_sNr``", 0),
		new ChemicalRule("c","gJP`@deVdB!gK``AddvPH#qir` qreH#!R@Jp@dpBl@ILslgp", 0),
		new ChemicalRule("d","gGQ@@eKtrAkDH!gGQ@@djsRIKVPP#qMsT qM\\V#!R_yL@dw~l@Jp@dsNR_@", 0),
		new ChemicalRule("e","daXJ@@PjyITuiX@`!dahJ@@SJYIMMfPB#IaLJfxP IaAhl[`#!ROrp?Ds|lOqNk?g?l_zLsSGp", 0),
// bad	new ChemicalRule("f","gJQ`@bdjt`P!gKa`@ldfrA@#qbqh qqlP#!R_zq?dw~l_yLsXgp", 0),
		new ChemicalRule("g","gBa`@lde@XS@!gCA`@mDPOdku`#qbq qJt#!R_zq?dxgFLvHB", 0),
		new ChemicalRule("h","daXL@@Klb`JSSHA@!daXL@@[dfTBijh@`#IyjFDp` IyjEL``#!R?g~HO[\\}[Lfw?K_}m?psLap", 0),
		new ChemicalRule("i","daXD@@kHh`Ttl@P!daxD@@yIeRfZj`B#IqBbay` IqBbnXP#!R?g}wOP`}]tKCV\\aBrCtsXep", 0),
		new ChemicalRule("j","gBi@DDcXCnR`!gBi@HDFXCnY`#qbq qfQ#!R_vpy`W}lLvK|", 0),
		new ChemicalRule("k","gB`@DcTB!gB`@DcTB#qbq qfQ#!R_vsY`W}lLvK|", 0),
		new ChemicalRule("l","daxH@@RUUjjPB!daDH@@RVU[jij@H#IGfhlR@ IGfbdK@#!R@IL@k@BS@Jp@dpBl@ILs|kp", 0),
		new ChemicalRule("m","gFQHBGAISiPH!gGQHJGAIUjhB#qNT] qNTk#!R@AL@[@@Sa_x@DsNro@", 0),
		new ChemicalRule("n","gOp`@|dRaqij`H!gOp`@tiUMij`H#qosJ` qgSqh#!RTv]`YRqg?g~HPh`}L{H|", 0),
		new ChemicalRule("o","gGQHDDqIUjdB!gGQHHDqAeMjhB#qbqk qfQ]#!R_zq?dw~l_yM?kCM|?@", 0),
// check	new ChemicalRule("p","gGQ@@eJuRA@!gFQ@@diuPD#qqUc qrcM#!R_zp@kG~S@IM?kCLb?@", 0),
		new ChemicalRule("r","gOQdEjHbAFQRBrtAaJ!gNQdEbHbABCTKPFDH#qbMwX qbM~X#!RCwpb@@M`CpL}cg@CL|jB", 0),
		new ChemicalRule("s","gOp`ATigujj`H!gOp`ATiVKjj`H#qnyS` qY~eP#!R?`@_YQ|ZFBqSFHc}L{IB", 0),
		new ChemicalRule("t","gOP`@dcUZdB!gNp`@tiTMjj@`#q{ir` qkLrx#!R@ANZXPAl@AL@[@@SLtj|", 0),
		new ChemicalRule("u","daXB`Hrn@HrPEILt`D!daXB`DJn@HRUMjV@H#IxVaLJ` IylJhPP#!R_zL@hs`Q_zq?dw~l_yLsBkp", 0),

		new ChemicalRule("aldol", "gFP`Adduf@payIzK@!gFP`ATeQfDU}K@#qisT qirc#!R@Jp@dqK~@Jp@dsNj@`", 0),
		};

	private static boolean sInitialized;

	private StereoMolecule mReactant,mProduct;
	private float mScore;
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
		mHistory = new StringBuilder();

		StereoMolecule reactant = new StereoMolecule(); // reusable container

		for (ChemicalRule rule:CHEMICAL_RULE) {
			reactantSearcher.setFragment(rule.getReactant());
			if (0 != reactantSearcher.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode)) {
				productSearcher.setFragment(rule.getProduct());
				if (productSearcher.isFragmentInMolecule()) {
mHistory.append(rule.getName()+":");
float historyScore = -10000;
					for (int[] reactantMatch:reactantSearcher.getMatchList()) {
						mReactant.copyMolecule(reactant);
						rule.apply(reactant, reactantMatch);
						int[] reactantMapNo = new int[mReactant.getAtoms()];
						int[] productMapNo = new int[mProduct.getAtoms()];
						mapper.map(reactant, mProduct, reactantMapNo, productMapNo);
						float score = mapper.calculateScore() + rule.getPanalty();
if (historyScore < score) historyScore = score;
						if (mScore < score) {
							mScore = score;
							bestReactantMapNo = reactantMapNo;
							bestProductMapNo = productMapNo;
							bestGraphMapNoCount = mapper.getGraphMapNoCount();
							mAppliedRule = rule;
							}
						}
mHistory.append(historyScore);
mHistory.append("\n");
					}
				}
			}

		// map and score the reaction without applying and rules
		int[] reactantMapNo = new int[mReactant.getAtoms()];
		int[] productMapNo = new int[mProduct.getAtoms()];
		mapper.map(mReactant, mProduct, reactantMapNo, productMapNo);
		int score = mapper.calculateScore();
mHistory.append("no rule:"+score);
mHistory.append("\n");
		if (mScore <= score) {
			mAppliedRule = null;
			mScore = score;
			bestReactantMapNo = reactantMapNo;
			bestProductMapNo = productMapNo;
			bestGraphMapNoCount = mapper.getGraphMapNoCount();
			}

		if (mScore != Integer.MIN_VALUE)
			mapper.copyMapNosToReaction(rxn, bestReactantMapNo, bestProductMapNo, bestGraphMapNoCount);
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
