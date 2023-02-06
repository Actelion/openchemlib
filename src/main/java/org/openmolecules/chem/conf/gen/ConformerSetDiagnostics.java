package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.util.DoubleFormat;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class ConformerSetDiagnostics extends ArrayList<ConformerDiagnostics> {
	private RigidFragment[] mRigidFragment;
	private RotatableBond[] mRotatableBond;
	private String mExitReason;

	protected ConformerSetDiagnostics(RotatableBond[] rb, RigidFragment[] rf) {
		mRotatableBond = rb;
		mRigidFragment = rf;
	}

	protected void addNew(TorsionSet ts) {
		add(new ConformerDiagnostics(ts));
	}

	protected ConformerDiagnostics getCurrent() {
		return get(size()-1);
	}

	protected void setExitReason(String er) {
		mExitReason = er;
	}

	public int getTorsion(int bond, int index) {
		return mRotatableBond[bond].getTorsion(index);
	}

	public int getRotatableBondCount() {
		return mRotatableBond.length;
	}

	public int getRigidFragmentCount() {
		return mRigidFragment.length;
	}

	public double getRigidFragmentLikelyhood(int conformer, int index) {
		return mRigidFragment[conformer].getConformerLikelihood(index);
	}

	public double getTorsionLikelyhood(int bond, int index) {
		return mRotatableBond[bond].getTorsionLikelyhood(index);
	}

	public String getExitReason() {
		return mExitReason;
	}

	public void writeEliminationRuleFile(String path) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(path));
			writeDataWarriorHeader(writer, true);
			writer.write("Structure\tcoords\tconformer\telim_rules");
			writer.newLine();
			int conformer = 0;
			for (ConformerDiagnostics cd:this) {
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

	public void writeAllConformersFile(String fileName) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			writeDataWarriorHeader(writer, true);
			writer.write("Structure\tcoords\tsuccess\tcollision");
			for (int i=1; i<=mRotatableBond.length; i++)
				writer.write("\ttorsion "+i);
			for (int i=1; i<=mRigidFragment.length; i++)
				writer.write("\tfragment "+i);
			writer.newLine();
			for (ConformerDiagnostics cd:this) {
				writer.write(cd.getIDCode() + "\t" + cd.getCoords() + "\t" + (cd.isSuccess() ? "yes" : "no") + "\t" + DoubleFormat.toString(cd.getCollisionIntensity()));
				int[] ti = cd.getTorsionIndexes();
				for (int rb=0; rb<ti.length; rb++)
					writer.write("\t" + ti[rb] + " (" + mRotatableBond[rb].getTorsion(ti[rb]) + ", " + DoubleFormat.toString(mRotatableBond[rb].getTorsionLikelyhood(ti[rb])) + ")");
				int[] ci = cd.getRigidFragmentIndexes();
				for (int rf=0; rf<ci.length; rf++)
					writer.write("\t" + ci[rf] + " (" + DoubleFormat.toString(mRigidFragment[rf].getConformerLikelihood(ci[rf])) + ")");
				writer.newLine();
			}
			writeDataWarriorFooter(writer);
			writer.close();
		}
		catch (IOException ioe) {}
	}

	public void writeRigidFragmentFile(String path) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(path));
			writeDataWarriorHeader(writer, true);
			writer.write("Structure\tcoords\tfragment\tconformer\tlikelyhood");
			writer.newLine();
			for (int i=0; i<mRigidFragment.length; i++) {
				RigidFragment f = mRigidFragment[i];
				for (int j=0; j<f.getConformerCount(); j++) {
					Canonizer canonizer = new Canonizer(f.getConformer(j).toMolecule());
					String idcode = canonizer.getIDCode();
					String coords = canonizer.getEncodedCoordinates();
					writer.write(idcode + "\t" + coords + "\t" + i + "\t" + j + "\t" + DoubleFormat.toString(f.getConformerLikelihood(j)));
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
			BufferedWriter writer = new BufferedWriter(new FileWriter(path));
			writeDataWarriorHeader(writer, false);
			writer.write("Structure\ttorsion-ID\tfragments\ttorsions\tlikelyhoods\trelevance");
			writer.newLine();
			for (int i=0; i<mRotatableBond.length; i++) {
				RotatableBond rb = mRotatableBond[i];
				String torsionID = rb.getTorsionID();
				writer.write(torsionID + "\t" + torsionID + "\t");
				writer.write(rb.getFragmentNo(0)+" "+rb.getFragmentNo(1) + "\t");
				for (int j=0; j<rb.getTorsionCount(); j++) {
					writer.write(Integer.toString(rb.getTorsion(j)));
					writer.write(j==rb.getTorsionCount()-1 ? "\t" : "<NL>");
				}
				for (int j=0; j<rb.getTorsionCount(); j++) {
					writer.write(DoubleFormat.toString(rb.getTorsionLikelyhood(j)));
					writer.write(j==rb.getTorsionCount()-1 ? "\t" : "<NL>");
				}
				writer.write(DoubleFormat.toString(rb.getRelevance()));
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
