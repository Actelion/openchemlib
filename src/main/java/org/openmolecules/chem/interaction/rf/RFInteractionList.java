package org.openmolecules.chem.interaction.rf;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.VDWRadii;
import org.openmolecules.chem.interaction.AtomClassifier;

import java.util.ArrayList;
import java.util.TreeSet;

public class RFInteractionList extends ArrayList<RFInteractionList.RFInteraction> {
	private final StereoMolecule mLigand, mProtein;

	public RFInteractionList(StereoMolecule ligand, StereoMolecule protein, boolean shortestInteractionOnly) {
		super();

		mLigand = ligand;
		mProtein = protein;

		mLigand.ensureHelperArrays(Molecule.cHelperNeighbours);
		mProtein.ensureHelperArrays(Molecule.cHelperNeighbours);

		ArrayList<InteractionCandidate> allCandidates = new ArrayList<>();
		TreeSet<InteractionCandidate> candidateSet = new TreeSet<>();	// candidates of one ligand atom
		for (int lAtom = 0; lAtom<mLigand.getAtoms(); lAtom++) {

			// Find close neighbors on protein side as interaction candidates
			// and sort them ascending on distance between VDW-spheres.
			double lRadius = VDWRadii.getVDWRadius(mLigand.getAtomicNo(lAtom));
			for (int pAtom = 0; pAtom<mProtein.getAtoms(); pAtom++) {
				double pRadius = VDWRadii.getVDWRadius(mProtein.getAtomicNo(pAtom));
				double distance = mLigand.getAtomCoordinates(lAtom).distance(mProtein.getAtomCoordinates(pAtom));
				double relDistance = distance - lRadius - pRadius;
				if (relDistance < 1.0)
					candidateSet.add(new InteractionCandidate(lAtom, pAtom, distance, relDistance));
			}

			// Collect closest line-of-sight neighbor as the interaction partner for this ligand atom
			if (!candidateSet.isEmpty()) {
				for (InteractionCandidate candidate : candidateSet) {
					if (candidate.isLineOfSight(candidateSet)) {
						allCandidates.add(candidate);
						if (shortestInteractionOnly)
							break;
					}
				}
				candidateSet.clear();
			}
		}

		int[] ligandAtomType = RFLigandAtomClassifier.getDefaultInstance().classifyAtoms(mLigand);
		int[] proteinAtomType = RFProteinAtomClassifier.getDefaultInstance().classifyAtoms(mProtein);

		// where a protein atom has more than one ligand contact, keep only the closest
		allCandidates.sort((o1, o2) -> (o1.pAtom != o2.pAtom) ?
			Integer.compare(o1.pAtom, o2.pAtom) : Double.compare(o1.relDistance, o2.relDistance));
		int previousProteinAtom = -1;
		for (InteractionCandidate candidate : allCandidates) {
			if (candidate.pAtom > previousProteinAtom || !shortestInteractionOnly) {
				previousProteinAtom = candidate.pAtom;
				if (ligandAtomType[candidate.lAtom] != AtomClassifier.TYPE_UNKNOWN
				 && proteinAtomType[candidate.pAtom] != AtomClassifier.TYPE_UNKNOWN)
					add(new RFInteraction(candidate.pAtom, candidate.lAtom,
							proteinAtomType[candidate.pAtom], ligandAtomType[candidate.lAtom],
							candidate.calculateAngle(), candidate.getDistance()));
			}
		}
	}

	private class InteractionCandidate implements Comparable<InteractionCandidate> {
		public int lAtom,pAtom;
		public double distance,relDistance;
		private final Coordinates p1,p2,dir;

		public InteractionCandidate(int lAtom, int pAtom, double distance, double relDistance) {
			this.lAtom = lAtom;
			this.pAtom = pAtom;
			this.distance = distance;
			this.relDistance = relDistance;
			Coordinates unitL2P = mProtein.getAtomCoordinates(pAtom).subC(mLigand.getAtomCoordinates(lAtom)).unit();
			p1 = mLigand.getAtomCoordinates(lAtom).addC(unitL2P.scaleC(VDWRadii.getVDWRadius(mLigand.getAtomicNo(lAtom))));
			p2 = mProtein.getAtomCoordinates(pAtom).subC(unitL2P.scaleC(VDWRadii.getVDWRadius(mProtein.getAtomicNo(pAtom))));
			dir = p2.subC(p1);
		}

		@Override
		public int compareTo(InteractionCandidate o) {
			return Double.compare(this.relDistance, o.relDistance);
		}

		private boolean isLineOfSight(TreeSet<InteractionCandidate> candidateSet) {
			for (InteractionCandidate candidate : candidateSet)
				if (candidate != this && isCollidingAtom(mProtein, candidate.pAtom))
					return false;

			for (int atom=0; atom<mLigand.getAtoms(); atom++)
				if (isCollidingAtom(mLigand, atom))
					return false;

			return true;
		}

		private boolean isCollidingAtom(StereoMolecule mol, int atom) {
			if ((mol == mProtein && atom == pAtom) || (mol == mLigand && atom == lAtom))
				return false;

			Coordinates pa = mol.getAtomCoordinates(atom);
			double radiusSquare = VDWRadii.getVDWRadius(mol.getAtomicNo(atom));
			radiusSquare = radiusSquare * radiusSquare;

/*if (mLigand.getAtomicNo(lAtom) == 6 && mLigand.getConnAtoms(lAtom) == 4
 && mProtein.getAtomicNo(pAtom) == 6 && mProtein.getAtomPi(pAtom) == 1
 && mol == mProtein && mol.getAtomicNo(atom) == 8 && mol.getAtomPi(atom) == 1) {
	System.out.println("RFInteractionList.isCollidingAtom(" + (mol == mLigand ? "ligand" : "protein") + "," + mol.getAtomLabel(atom) + atom + ") lAtomDistance:"+DoubleFormat.toString(mLigand.getAtomCoordinates(lAtom).distance(mol.getAtomCoordinates(atom)))+" LAtom:"+mLigand.getAtomLabel(lAtom)+lAtom+" PAtom:"+mProtein.getAtomLabel(pAtom) + pAtom+" -------------- below");
	System.out.println("aX:"+ DoubleFormat.toString(pa.x)+" aY:"+ DoubleFormat.toString(pa.y)+" aZ:"+ DoubleFormat.toString(pa.z));
	Coordinates pl = mLigand.getAtomCoordinates(lAtom);
	System.out.println("lX:"+ DoubleFormat.toString(pl.x)+" lY:"+ DoubleFormat.toString(pl.y)+" lZ:"+ DoubleFormat.toString(pl.z));
	Coordinates pp = mProtein.getAtomCoordinates(pAtom);
	System.out.println("pX:"+ DoubleFormat.toString(pp.x)+" pY:"+ DoubleFormat.toString(pp.y)+" pZ:"+ DoubleFormat.toString(pp.z));
	System.out.println("p1:"+ DoubleFormat.toString(p1.x)+" p1:"+ DoubleFormat.toString(p1.y)+" p1:"+ DoubleFormat.toString(p1.z));
	System.out.println("p2:"+ DoubleFormat.toString(p2.x)+" p2:"+ DoubleFormat.toString(p2.y)+" p2:"+ DoubleFormat.toString(p2.z));
	String collides = (relDistance <= 0) ? (pa.distanceSquared(p1)<radiusSquare && pa.distanceSquared(p2)<radiusSquare ? "yes" : "no") : "next";
	System.out.println("#1 collides:"+collides+" relDistance:"+DoubleFormat.toString(relDistance)+" pa.distance(p1):"+DoubleFormat.toString(pa.distance(p1))+" pa.distance(p2):"+DoubleFormat.toString(pa.distance(p2))+" radius:"+DoubleFormat.toString(Math.sqrt(radiusSquare)));

	double distanceSquared = pa.subC(p1).cross(dir).getLength() / relDistance;
	distanceSquared = distanceSquared * distanceSquared;
	collides = (distanceSquared>radiusSquare) ? "no" : "next";
	System.out.println("#2 collides:"+collides+" distanceToP1P2:"+DoubleFormat.toString(Math.sqrt(distanceSquared))+" pa.distance(p1):"+DoubleFormat.toString(pa.distance(p1))+" pa.distance(p2):"+DoubleFormat.toString(pa.distance(p2))+" radius:"+DoubleFormat.toString(Math.sqrt(radiusSquare)));

	double max = relDistance * relDistance + distanceSquared;
	collides = (pa.distanceSquared(p1) > max) ? (pa.distanceSquared(p2)<radiusSquare ? "yes" : "no") : "next";
	System.out.println("#3 collides:"+collides+" max:"+DoubleFormat.toString(Math.sqrt(max))+" pa.distance(p1)):"+DoubleFormat.toString(pa.distance(p1)));

	collides = (pa.distanceSquared(p2) > max) ? (pa.distanceSquared(p1)<radiusSquare ? "yes" : "no") : "next yes";
	System.out.println("#4 collides:"+collides+" max:"+DoubleFormat.toString(Math.sqrt(max))+" pa.distance(p2)):"+DoubleFormat.toString(pa.distance(p2)));
}*/

			if (relDistance<=0)
				return pa.distanceSquared(p1)<radiusSquare && pa.distanceSquared(p2)<radiusSquare;

			// distance from straight line (p1->p2) to atom
			double distanceSquared = pa.subC(p1).cross(dir).getLength() / relDistance;
			distanceSquared = distanceSquared * distanceSquared;
			if (distanceSquared>radiusSquare)
				return false;

			double max = relDistance * relDistance + distanceSquared;
			if (pa.distanceSquared(p1) > max)
				return pa.distanceSquared(p2)<radiusSquare;
			if (pa.distanceSquared(p2) > max)
				return pa.distanceSquared(p1)<radiusSquare;
			return true;
		}

		public double calculateAngle() {
			return 0.0;
		}

		public double getDistance() {
			// TODO possibly calculate distance rectangular to pi-systems
			return distance;
		}
	}

	public static class RFInteraction {
		public int pAtom,lAtom,pType,lType;
		public double angle,relDistance;

		RFInteraction(int pAtom, int lAtom, int pType, int lType, double angle, double relDistance) {
			this.pAtom = pAtom;
			this.lAtom = lAtom;
			this.pType = pType;
			this.lType = lType;
			this.angle = angle;
			this.relDistance = relDistance;
		}
	}
}