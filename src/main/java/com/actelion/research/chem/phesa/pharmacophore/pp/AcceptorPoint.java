package com.actelion.research.chem.phesa.pharmacophore.pp;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.InteractionSimilarityTable;
import com.actelion.research.chem.phesaflex.MathHelper;

public class AcceptorPoint implements IPharmacophorePoint {
	private int acceptorAtom;
	private List<Integer> neighbours;
	private Coordinates directionality;
	private int interactionClass;
	private Coordinates center;
	private int acceptorID; //necessary to assign different directionalities to two acceptor points in sp2 oxygen
	
	public AcceptorPoint(StereoMolecule mol, int a, List<Integer> neighbours, int interactionClass) {
		this(mol, a, neighbours, interactionClass, 0);
	}
	
	public AcceptorPoint(AcceptorPoint aP) {
		acceptorAtom = aP.acceptorAtom;
		neighbours = new ArrayList<Integer>();
		for(int neighbour : aP.neighbours) {
			neighbours.add(neighbour);
		}
		
		directionality = new Coordinates(aP.directionality);
		interactionClass = aP.interactionClass;
		center = new Coordinates(aP.center);
		acceptorID = aP.acceptorID;
	}
	
	public AcceptorPoint(StereoMolecule mol, int a, List<Integer> neighbours, int interactionClass, int acceptorID) {
		acceptorAtom = a;
		this.neighbours = neighbours;
		this.interactionClass = interactionClass;
		this.acceptorID = acceptorID;
		updateCoordinates(mol.getAtomCoordinates());
	}
	
	private AcceptorPoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public static AcceptorPoint fromString(String ppString, StereoMolecule mol) {
		return new AcceptorPoint(ppString,mol);
	}


	
	@Override
	public void updateCoordinates(Coordinates[] coords) {
		center = new Coordinates(coords[acceptorAtom].x ,coords[acceptorAtom].y, coords[acceptorAtom].z);
		if(neighbours.size()==1) {
			int aa1 = neighbours.get(0);
			directionality = center.subC(coords[aa1]);
		}
			
			
		else if(neighbours.size()==2 && acceptorID!=0) {
			int aa1 = neighbours.get(0);
			Coordinates v1 = center.subC(coords[aa1]);
			int aa2 = neighbours.get(1);
			Coordinates v2 = coords[aa2].subC(center);
			Coordinates rotAxis = v1.cross(v2).unit();
			double theta = acceptorID == 1 ? 45.0/180.0*Math.PI :  -45.0/180.0*Math.PI;
			directionality = v1.rotate(rotAxis, theta);
		}
		
		
		else if(neighbours.size()==3) {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			int aa3 = neighbours.get(2);
			Coordinates v1 = center.subC(coords[aa1]).unit();
			Coordinates v2 = center.subC(coords[aa2]).unit();
			Coordinates v3 = center.subC(coords[aa3]).unit();
			directionality = v3.add(v2).add(v1);
		}
		
		
		else {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			Coordinates v1 = center.subC(coords[aa1]).unit();
			Coordinates v2 = center.subC(coords[aa2]).unit();
			directionality = v1.addC(v2);
		}
		directionality.unit();
		// TODO Auto-generated method stub
		
	}
	


	@Override
	public Coordinates getCenter() {
		return center;
	}

	@Override
	public Coordinates getDirectionality() {
		return directionality;
	}

	@Override
	public String encode() {
		StringBuilder molVolString = new StringBuilder();
		molVolString.append("a");
		molVolString.append(" ");
		molVolString.append(Integer.toString(acceptorAtom));
		molVolString.append(" ");
		molVolString.append(Integer.toString(interactionClass));
		molVolString.append(" ");
		molVolString.append(Integer.toString(acceptorID));
		molVolString.append(" ");
		//molVolString.append(Integer.toString(neighbours.size()));
		//molVolString.append(" ");
		for(Integer neighbour : neighbours) {
			molVolString.append(neighbour);
			molVolString.append(" ");
		}
		return molVolString.toString().trim();
	}
	

	private void decode(String ppString, StereoMolecule mol) {
		String[] strings = ppString.split(" ");
		acceptorAtom = Integer.decode(strings[1]);
		interactionClass = Integer.decode(strings[2]);
		acceptorID = Integer.decode(strings[3]);
		neighbours = new ArrayList<Integer>();
		for(int i=4;i<strings.length;i++) {
			neighbours.add(Integer.decode(strings[i]));
		}
		updateCoordinates(mol.getAtomCoordinates());
	}
	
	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		if(pp instanceof AcceptorPoint) {
			return 1.0;
		}
		return 0.0;
	}

	public int getInteractionClass() {
		return interactionClass;
	}

	@Override
	public int getCenterID() {
		return acceptorAtom;
	}
	
	@Override
	public void setCenterID(int centerID) {
		acceptorAtom = centerID;
	}
	
	@Override
	public void setDirectionality(Coordinates directionality) {
		this.directionality = directionality;
		
	}
	
	public int getAcceptorID() {
		return acceptorID;
	}
	
	@Override
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim) {
		
		 /*
		  * di and dj are the the directionality vectors of the reference and fit pharmacophores respectively
		  * the similarity of the directionalities is calculated by the crossproduct of di and dj
		  * Ssim,vec = dix*djx + diy*djy + diz*djz
		  * the derivative of the similarity with respect to the atomic coordinates of the reference molecule xj:
		  *  dSsim,vec/dxj = dix*djx/dxj + diy*djy/dyj + diz*djz/dzj 
		  */
       
		
		/*
		 * AA1-----A ------>d
		 * 
		 * d = A-AA1
		 * 
		 * 
		 * 
		 * 
		 */
		
		
		
		if(neighbours.size()==1) {
			int aa1 = neighbours.get(0);
			grad[3*acceptorAtom] = sim*di.x/3.0;
			grad[3*acceptorAtom+1] = sim*di.y/3.0;
			grad[3*acceptorAtom+2] = sim*di.z/3.0;
			grad[3*aa1] = sim*-di.x/3.0;
			grad[3*aa1+1] = sim*-di.y/3.0;
			grad[3*aa1+2] = sim*-di.z/3.0;
			
			}
			
			
		/*
		 *                        d: directionality vector of lone pair
		 C1                       v2: C1-C2  
		  \              /        v1: A-C2 
		   \ v2       d /         d is constructed by rotating v1 by the mount theta around the axis u that
		    \          /          is perpendicular to v2 and v1 -> u = v1xv2 
		     C2=======A           d depends on C1,C2,A -> d(C1,C2,A)
		    /     v1              we therefore need to calculate the derivatives dd/dC1, dd/dC2, dd/dA 
		   /                      d = R(u,theta)*v1    theta is constant
		  /  *                    dd/dC1 = dR/du*du/dC1*v1+R*d(O4-O2)/dC1
		 C3                       dd/dC1 = dR/du*du/dC2*v1  (since v1= A-C2 is not a function of C1
		                          u = v1xv2 
		                          du/dC1 = dv2/dC1 x v1

		*/
		else if(neighbours.size()==2 && acceptorID!=0) { //sp2 oxygen 
			Coordinates centerCoords = new Coordinates(v[3*acceptorAtom],v[3*acceptorAtom+1],v[3*acceptorAtom+2]);
			int c2 = neighbours.get(0);
			Coordinates c2Coords = new Coordinates(v[3*c2],v[3*c2+1],v[3*c2+2]);
			Coordinates v1 = centerCoords.subC(c2Coords);
			int c1 = neighbours.get(1);
			Coordinates c1Coords = new Coordinates(v[3*c1],v[3*c1+1],v[3*c1+2]);
			Coordinates v2 = c1Coords.subC(centerCoords);
			Coordinates u =  v1.cross(v2).unit();
			double theta = acceptorID == 1 ? 45.0/180.0*Math.PI :  -45.0/180.0*Math.PI;
			double[][] r = new double[3][3];
			Coordinates[][] drdu = new Coordinates[3][3];
			double[][] drdc1 = new double[3][3];
			MathHelper.getRotMatrix(u, theta, r);
			MathHelper.getRotMatrixDerivative(u, theta, drdu);
			Coordinates dv2dc1 = new Coordinates(1,1,1);
			Coordinates dudc1 = dv2dc1.cross(v1);
			drdc1[0][0] = drdu[0][0].dot(dudc1);
			drdc1[0][1] = drdu[0][1].dot(dudc1);
			drdc1[0][2] = drdu[0][2].dot(dudc1);
			drdc1[1][0] = drdu[1][0].dot(dudc1);
			drdc1[1][1] = drdu[1][1].dot(dudc1);
			drdc1[1][2] = drdu[1][2].dot(dudc1);
			drdc1[2][0] = drdu[2][0].dot(dudc1);
			drdc1[2][1] = drdu[2][2].dot(dudc1);
			drdc1[2][2] = drdu[2][2].dot(dudc1);
			Coordinates dddc1 = new Coordinates(drdc1[0][0]*v1.x+drdc1[0][1]*v1.y + drdc1[0][2]*v1.z, // dd/dc1
					drdc1[1][0]*v1.x+drdc1[1][1]*v1.y + drdc1[2][1]*v1.z,
					drdc1[2][0]*v1.x+drdc1[2][1]*v1.y + drdc1[2][2]*v1.z);
			grad[3*c1] += sim*dddc1.x/3.0;
			grad[3*c1+1] += sim*dddc1.y/3.0;
			grad[3*c1+2] += sim*dddc1.z/3.0;
			Coordinates dv1dc2 = new Coordinates(-1,-1,-1);
			Coordinates dv2dc2 = dv1dc2;
			Coordinates dudc2 = dv1dc2.cross(v2).add(dv2dc2.cross(v1));
			double [][]drdc2 = drdc1;
			drdc2[0][0] = drdu[0][0].dot(dudc2);
			drdc2[0][1] = drdu[0][1].dot(dudc2);
			drdc2[0][2] = drdu[0][2].dot(dudc2);
			drdc2[1][0] = drdu[1][0].dot(dudc2);
			drdc2[1][1] = drdu[1][1].dot(dudc2);
			drdc2[1][2] = drdu[1][2].dot(dudc2);
			drdc2[2][0] = drdu[2][0].dot(dudc2);
			drdc2[2][1] = drdu[2][2].dot(dudc2);
			drdc2[2][2] = drdu[2][2].dot(dudc2);
			Coordinates dddc2 = new Coordinates(drdc1[0][0]*v1.x+drdc1[0][1]*v1.y + drdc1[0][2]*v1.z + 
					r[0][0]*dv1dc2.x + r[0][1]*dv1dc2.y + r[0][2]*dv1dc2.z, 
					drdc1[1][0]*v1.x+drdc1[1][1]*v1.y + drdc1[1][2]*v1.z + 
					r[1][0]*dv1dc2.x + r[1][1]*dv1dc2.y + r[1][2]*dv1dc2.z, 
					drdc1[2][0]*v1.x+drdc1[2][1]*v1.y + drdc1[2][2]*v1.z + 
					r[2][0]*dv1dc2.x + r[2][1]*dv1dc2.y + r[2][2]*dv1dc2.z);
			grad[3*c2] += sim*dddc2.x/3.0;  //directionality is scaled by 1/3
			grad[3*c2+1] += sim*dddc2.y/3.0;
			grad[3*c2+2] += sim*dddc2.z/3.0;
			Coordinates dv1da = dv2dc1;
			Coordinates duda = dv1da.cross(v2);
			double [][]drda = drdc1;
			drda[0][0] = drdu[0][0].dot(duda);
			drda[0][1] = drdu[0][1].dot(duda);
			drda[0][2] = drdu[0][2].dot(duda);
			drda[1][0] = drdu[1][0].dot(duda);
			drda[1][1] = drdu[1][1].dot(duda);
			drda[1][2] = drdu[1][2].dot(duda);
			drda[2][0] = drdu[2][0].dot(duda);
			drda[2][1] = drdu[2][2].dot(duda);
			drda[2][2] = drdu[2][2].dot(duda);
			Coordinates ddda = new Coordinates(drda[0][0]*v1.x+drda[0][1]*v1.y + drda[0][2]*v1.z + 
					r[0][0]*dv1da.x + r[0][1]*dv1da.y + r[0][2]*dv1da.z, 
					drda[1][0]*v1.x+drda[1][1]*v1.y + drda[1][2]*v1.z + 
					r[1][0]*dv1da.x + r[1][1]*dv1da.y + r[1][2]*dv1da.z, 
					drda[2][0]*v1.x+drda[2][1]*v1.y + drda[2][2]*v1.z + 
					r[2][0]*dv1da.x + r[2][1]*dv1da.y + r[2][2]*dv1da.z);
			
			grad[3*acceptorAtom] += sim*ddda.x/3.0;
			grad[3*acceptorAtom+1] += sim*ddda.y/3.0;
			grad[3*acceptorAtom+2] += sim*ddda.z/3.0;
			
		}
		
		else if(neighbours.size()==3) {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			int aa3 = neighbours.get(2);
			grad[3*acceptorAtom] += sim*3*di.x/3.0;
			grad[3*acceptorAtom+1] += sim*3*di.y/3.0;
			grad[3*acceptorAtom+2] += sim*3*di.z/3.0;
			grad[3*aa1] += sim*-di.x/3.0;
			grad[3*aa1+1] += sim*-di.y/3.0;
			grad[3*aa1+2] += sim*-di.z/3.0;
			grad[3*aa2] += sim*-di.x/3.0;
			grad[3*aa2+1] += sim*-di.y/3.0;
			grad[3*aa2+2] += sim*-di.z/3.0;
			grad[3*aa3] += sim*-di.x/3.0;
			grad[3*aa3+1] += sim*-di.y/3.0;
			grad[3*aa3+2] += sim*-di.z/3.0;

		}
		
		else { //two neighbours, sp3 acceptor
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			grad[3*acceptorAtom] += sim*2*di.x/3.0;
			grad[3*acceptorAtom+1] += sim*2*di.y/3.0;
			grad[3*acceptorAtom+2] += sim*2*di.z/3.0;
			grad[3*aa1] += sim*-di.x/3.0;
			grad[3*aa1+1] += sim*-di.y/3.0;
			grad[3*aa1+2] += sim*-di.z/3.0;
			grad[3*aa2] += sim*-di.x/3.0;
			grad[3*aa2+1] += sim*-di.y/3.0;
			grad[3*aa2+2] += sim*-di.z/3.0;
			
			
		}
	}

	@Override
	public void updateAtomIndices(int[] map) {
		acceptorAtom = map[acceptorAtom];
		for(int i=0;i<neighbours.size();i++) {
			int neighbour = map[neighbours.get(i)];
			neighbours.set(i, neighbour);
		}
	}

	@Override
	public int[] getAtomIndices() {
		int [] a = {acceptorAtom};
		return a;
	}

	@Override
	public IPharmacophorePoint copyPharmacophorePoint() {
		return new AcceptorPoint(this);
	}

	@Override
	public int getFunctionalityIndex() {
		return IPharmacophorePoint.Functionality.ACCEPTOR.getIndex();
	}

	@Override
	public Coordinates getRotatedDirectionality(double[][] rotMatrix, double scaleFactor) {
		Coordinates directMod = new Coordinates();
		directMod = directionality.rotateC(rotMatrix);
		directMod.scale(scaleFactor); // scale by the invers
		return directMod;
	}


	
		
	
}

