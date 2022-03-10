package com.actelion.research.chem.optimization;

import java.util.Random;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.DockingUtils;

public class MCHelper {
	
	private BondRotationHelper torsionHelper;
	private int[] mcsRotBondIndeces;
	private Random random;
	private static double TEMPERATURE = 0.0043; //move that reduces Tanimoto by 0.01 has 10% chance of acceptance
	
	public MCHelper(BondRotationHelper torsionHelper, int[] mcsRotBondIndeces, Random random) {
		this.torsionHelper = torsionHelper;
		this.mcsRotBondIndeces = mcsRotBondIndeces;
		this.random = random;
	}
	
	
	private static double MOVE_AMPLITUDE = 2.0;
	
	public void randomPerturbation(Conformer conf, double[] state) {
		if(mcsRotBondIndeces==null) {
			int num = (int) (3*random.nextDouble());
			if(num==0) { //translation
				Coordinates shift = DockingUtils.randomVectorInSphere(random).scale(MOVE_AMPLITUDE);
				state[0]+=shift.x;
				state[1]+=shift.y;
				state[2]+=shift.z;
				/*
				for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) { 
					Coordinates c = ligConf.getCoordinates(a);
					c.add(shift);
				}
				*/
			}
			else if(num==1) {
				double r = getGyrationRadius(conf);
				Coordinates rot = DockingUtils.randomVectorInSphere(random).scale(MOVE_AMPLITUDE/r);
				double angle = rot.dist();
				Quaternion q = new Quaternion(1.0,0.0,0.0,0.0);
				if(angle>0.0001) {
					Coordinates axis = rot.scale(1.0/angle);
					q = new Quaternion(axis,angle);
				}
				ExponentialMap em = new ExponentialMap(state[3],state[4],state[5]);
				Quaternion qOrig = em.toQuaternion();
				q.multiply(qOrig);
				ExponentialMap emNew = new ExponentialMap(q);
				state[3] = emNew.getP().x;
				state[4] = emNew.getP().y;
				state[5] = emNew.getP().z;
	
				
			}
			else  {
				torsionPerturbation(conf,state);
				
			}
		}
		//mcs docking
		else {
			if(torsionHelper.getRotatableBonds().length==0)
				return;
			else if(mcsRotBondIndeces.length==0)
				return;
			double rnd = random.nextDouble();
			rnd*=180.0;
			double rotateBy = random.nextBoolean() ? rnd : -rnd;
			rotateBy = rotateBy*Math.PI/180.0;
			int bond = mcsRotBondIndeces[random.nextInt(mcsRotBondIndeces.length)];
			double previousAngle = state[6+bond];
			double newAngle = previousAngle+rotateBy; 
			if(newAngle>Math.PI) {
				newAngle -= 2*Math.PI;
			}
			state[6+bond] = newAngle;
			
			
		}
		
		
	}
	
	public void torsionPerturbation(Conformer ligConf, double[] state) {
		if(torsionHelper.getRotatableBonds().length==0)
			return;
		double rnd = random.nextDouble();
		rnd*=180.0;
		double rotateBy = random.nextBoolean() ? rnd : -rnd;
		rotateBy = rotateBy*Math.PI/180.0;
		int bond = random.nextInt(torsionHelper.getRotatableBonds().length);
		double previousAngle = state[6+bond];
		double newAngle = previousAngle+rotateBy; 
		if(newAngle>Math.PI) {
			newAngle -= 2*Math.PI;
		}
		state[6+bond] = newAngle;
	}
	
	public double getGyrationRadius(Conformer conf) {
		Coordinates com = DockingUtils.getCOM(conf);
		double r = 0.0;
		int counter = 0;
		for(int a=0;a<conf.getMolecule().getAtoms();a++) {
			Coordinates c = conf.getCoordinates(a);
			r+= c.distanceSquared(com);
			counter++;
		}
		r/=counter;
		return Math.sqrt(r);
	}

	public void setTorsionHelper(BondRotationHelper torsionHelper) {
		this.torsionHelper = torsionHelper;
	}

	public void setMcsRotBondIndeces(int[] mcsRotBondIndeces) {
		this.mcsRotBondIndeces = mcsRotBondIndeces;
	}
	
	public boolean accept(double oldScore, double newScore) {
		boolean accept = false;
		if(newScore>oldScore)
			accept = true;
		else {
			double delta = -(newScore-oldScore);
			double p = Math.exp(-delta/TEMPERATURE);
			double rnd = random.nextDouble();
			if(rnd<p)
				accept = true;
		}
		return accept;
	}

}
