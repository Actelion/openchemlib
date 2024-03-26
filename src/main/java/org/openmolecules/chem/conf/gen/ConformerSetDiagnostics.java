package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.util.DoubleFormat;
import org.openmolecules.chem.conf.so.ConformationRule;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.TreeMap;

public class ConformerSetDiagnostics {
	private TreeMap<String,ConformerDiagnostics> mDiagnosticsMap;
	private ConformerGenerator mConformerGenerator;
	private RigidFragment[] mRigidFragment;
	private RotatableBond[] mRotatableBond;
	private String mExitReason;

	protected ConformerSetDiagnostics(ConformerGenerator cg) {
		mConformerGenerator = cg;
		mRotatableBond = cg.getRotatableBonds();
		mRigidFragment = cg.getRigidFragments();
		mDiagnosticsMap = new TreeMap<>();
	}

	protected void addNew(TorsionSet ts) {
		mDiagnosticsMap.put(ts.getConformer().getName(), new ConformerDiagnostics(ts));
	}

	public Collection<ConformerDiagnostics> getDiagnostics() {
		return mDiagnosticsMap.values();
	}

	protected ConformerDiagnostics get(TorsionSet ts) {
		return mDiagnosticsMap.get(ts.getConformer().getName());
	}

	protected void setExitReason(String er) {
		mExitReason = er;
	}

//	public int getTorsion(int bond, int index) {
//		return mRotatableBond[bond].getTorsion(index);
//	}

	public int getRotatableBondCount() {
		return mRotatableBond.length;
	}

	public int getRigidFragmentCount() {
		return mRigidFragment.length;
	}

	public double getRigidFragmentLikelyhood(int conformer, int index) {
		return mRigidFragment[conformer].getConformerLikelihood(index);
	}

//	public double getTorsionLikelyhood(int bond, int index) {
//		return mRotatableBond[bond].getTorsionLikelyhood(index);
//	}

	public String getExitReason() {
		return mExitReason;
	}

	public void writeEliminationRuleFile(String path) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(path), StandardCharsets.UTF_8));
			writeDataWarriorHeader(writer, true);
			writer.write("Structure\tcoords\tconformer\telim_rules");
			writer.newLine();
			int conformer = 0;
			for (ConformerDiagnostics cd:mDiagnosticsMap.values()) {
				writer.write(cd.getIDCode());
				writer.write("\t");
				writer.write(cd.getCoords());
				writer.write("\t");
				writer.write(conformer++);
				writer.write("\t");
				for (String rule:cd.getEliminationRules()) {
					writer.write(rule);
					writer.write("<NL>");
				}
				writer.write("\n");
			}
			writeDataWarriorFooter(writer);
			writer.close();
		}
		catch (IOException ioe) {}
	}

	public void writePermutationSpace(String fileName) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), StandardCharsets.UTF_8));
			writer.write("rigid fragments");
			for (int i=1; i<=mRotatableBond.length; i++)
				writer.write("\ttorsion "+i);
			writer.write("\telimination rules\tconformers tried");
			writer.newLine();
			for (int[] rigidFragmentIndex:mConformerGenerator.getBaseConformerMap().keySet()) {
				for (int rfi:rigidFragmentIndex)
					writer.write(rfi+" ");
				writer.write("\t");
				BaseConformer bc = mConformerGenerator.getBaseConformer(rigidFragmentIndex);
				for (int rb=0; rb<mRotatableBond.length; rb++) {
					for (int ti=0; ti<mRotatableBond[rb].getTorsionCount(); ti++) {
						writer.write(bc.getTorsion(rb, ti)+" ("+DoubleFormat.toString(bc.getTorsionLikelyhood(rb, ti), 3)+")");
						writer.write(ti == mRotatableBond[rb].getTorsionCount()-1 ? "\t" : "<NL>");
					}
				}
				TorsionSetEncoder encoder = mConformerGenerator.getTorsionSetStrategy().getTorsionSetEncoder();
				for (TorsionSetEliminationRule er: bc.getEliminationRules()) {
					writer.write(encoder.createRuleString(er, bc));
					writer.write("<NL>");
				}
				writer.write("\t"+bc.getDerivedConformerCount());
				writer.newLine();
			}
			writeDataWarriorFooter(writer);
			writer.close();
		}
		catch (IOException ioe) {}
	}

	public void writeAllConformersFile(String fileName) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), StandardCharsets.UTF_8));
			writeDataWarriorHeader(writer, true);
			writer.write("Structure\tcoords\tname\tsuccess\tlikelihood\tcollision");
			for (int i=1; i<=mRigidFragment.length; i++)
				writer.write("\tfragment "+i);
			for (int i=1; i<=mRotatableBond.length; i++)
				writer.write("\ttorsion "+i);
			writer.write("\tnew elimination rules");
			writer.newLine();
			for (ConformerDiagnostics cd:mDiagnosticsMap.values()) {
				writer.write(cd.getIDCode() + "\t" + cd.getCoords() + "\t" + cd.getName() + "\t" + (cd.isSuccess() ? "yes" : "no") + "\t" + DoubleFormat.toString(cd.getLikelihood(),3) + "\t" + DoubleFormat.toString(cd.getCollisionIntensity(),3));
				int[] ci = cd.getRigidFragmentIndexes();
				for (int rf=0; rf<ci.length; rf++)
					writer.write("\t" + ci[rf] + " (" + DoubleFormat.toString(mRigidFragment[rf].getConformerLikelihood(ci[rf]),3) + ")");
				BaseConformer baseConformer = cd.getBaseConformer();
				int[] ti = cd.getTorsionIndexes();
				for (int rb=0; rb<ti.length; rb++) {
					writer.write("\t" + ti[rb] + " (" + baseConformer.getTorsion(rb, ti[rb]) + ", " + DoubleFormat.toString(baseConformer.getTorsionLikelyhood(rb, ti[rb]),3) + ")");
					int fixedTorsion = cd.getFixedTorsions()[rb];
					if (fixedTorsion != -1 && fixedTorsion != baseConformer.getTorsion(rb, ti[rb]))
						writer.write("<NL>optimized:"+fixedTorsion);
					}
				writer.write("\t");
				for (String rule:cd.getEliminationRules()) {
					writer.write(rule);
					writer.write("<NL>");
				}
				writer.newLine();
			}
			writeDataWarriorFooter(writer);
			writer.close();
		}
		catch (IOException ioe) {}
	}

	public void writeRigidFragmentFile(String path) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(path), StandardCharsets.UTF_8));
			writeDataWarriorHeader(writer, true);
			writer.write("Structure\tcoords\tfragment\tconformer\tlikelyhood\trule strain\tatom strain");
			writer.newLine();
			for (int i=0; i<mRigidFragment.length; i++) {
				RigidFragment f = mRigidFragment[i];
				for (int j=0; j<f.getConformerCount(); j++) {
					Canonizer canonizer = new Canonizer(f.getConformer(j).toMolecule());
					String idcode = canonizer.getIDCode();
					String coords = canonizer.getEncodedCoordinates();
					SelfOrganizedConformer soc = f.getConformer(j);
					writer.write(idcode + "\t" + coords + "\t" + i + "\t" + j + "\t"
							+ DoubleFormat.toString(f.getConformerLikelihood(j),3) + "\t");
					for (int r=0; r<ConformationRule.RULE_NAME.length; r++) {
						writer.write(ConformationRule.RULE_NAME[r]+":"+DoubleFormat.toString(soc.getRuleStrain(r),3));
						writer.write(r==ConformationRule.RULE_NAME.length-1 ? "\t" : "<NL>");
					}
					for (int a=0; a<soc.getSize(); a++) {
						writer.write("a"+a+":"+DoubleFormat.toString(soc.getAtomStrain(a),3));
						writer.write(a==soc.getSize()-1 ? "\t" : "<NL>");
					}
					writer.newLine();
				}
			}
			writeDataWarriorFooter(writer);
			writer.close();
		}
		catch (IOException ioe) {}
	}

	public void writeRotatableBondsFile(String path) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(path), StandardCharsets.UTF_8));
			writeDataWarriorHeader(writer, false);
			writer.write("Structure\ttorsion-ID\tfragments\ttorsions\tfrequenies\trelevance");
			writer.newLine();
			for (int i=0; i<mRotatableBond.length; i++) {
				RotatableBond rb = mRotatableBond[i];
				String torsionID = rb.getTorsionID();
				writer.write(torsionID + "\t" + torsionID + "\t");
				writer.write(rb.getFragmentNo(0)+" "+rb.getFragmentNo(1) + "\t");
				for (int j=0; j<rb.getTorsionCount(); j++) {
					writer.write(Integer.toString(rb.getDefaultTorsions()[j]));
					writer.write(j==rb.getTorsionCount()-1 ? "\t" : "<NL>");
				}
				for (int j=0; j<rb.getTorsionCount(); j++) {
					writer.write(rb.getDefaultFrequencies()[j]);
					writer.write(j==rb.getTorsionCount()-1 ? "\t" : "<NL>");
				}
				writer.write(DoubleFormat.toString(rb.getRelevance(),3));
				writer.newLine();
			}
			writeDataWarriorFooter(writer);
			writer.close();
		}
		catch (IOException ioe) {}
	}

	private void writeDataWarriorHeader(BufferedWriter writer, boolean includeCoords) throws IOException {
		writer.write("<column properties>");
		writer.newLine();
		writer.write("<columnName=\"Structure\">");
		writer.newLine();
		writer.write("<columnProperty=\"specialType\tidcode\">");
		writer.newLine();
		if (includeCoords) {
			writer.write("<columnName=\"coords\">");
			writer.newLine();
			writer.write("<columnProperty=\"specialType\tidcoordinates3D\">");
			writer.newLine();
			writer.write("<columnProperty=\"parent\tStructure\">");
			writer.newLine();
		}
		writer.write("</column properties>");
		writer.newLine();
	}

	private void writeDataWarriorFooter(BufferedWriter writer) throws IOException {
		writer.write("<datawarrior properties>");
		writer.newLine();
		writer.write("<useDefault=\"filtersOnly\">");
		writer.newLine();
		writer.write("</datawarrior properties>");
		writer.newLine();
	}
}
