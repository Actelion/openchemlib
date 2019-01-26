package org.openmolecules.chem.conf.render;

import com.actelion.research.chem.Coordinates;

import java.util.ArrayList;
import java.util.Arrays;

public class CachingMoleculeBuilder implements MoleculeBuilder {
	private MoleculeBuilder mRenderer;
	private ArrayList<CachedSphere> mSphereList;
	private ArrayList<CachedCylinder> mCylinderList;

	public CachingMoleculeBuilder(MoleculeBuilder renderer) {
		mRenderer = renderer;
		mSphereList = new ArrayList<CachedSphere>();
		mCylinderList = new  ArrayList<CachedCylinder>();
		}

	@Override
	public void init() {
		mSphereList.clear();
		mCylinderList.clear();
		}

		@Override
	public void addSphere(int atom, int bond, Coordinates c, double radius, int argb) {
		mSphereList.add(new CachedSphere(atom, bond, c, radius, argb));
		}

	@Override
	public void addCylinder(int bond, double radius, double length, Coordinates c, double rotationY, double rotationZ, int argb) {
		mCylinderList.add(new CachedCylinder(bond, radius, length, c, rotationY, rotationZ, argb));
		}

	@Override
	public void done() {
		CachedCylinder[] cylinders = mCylinderList.toArray(new CachedCylinder[0]);
		Arrays.sort(cylinders);
		for (CachedCylinder c:cylinders)
			mRenderer.addCylinder(c.bond, c.radius, c.length, c.c, c.rotationY, c.rotationZ, c.argb);

		CachedSphere[] spheres = mSphereList.toArray(new CachedSphere[0]);
		Arrays.sort(spheres);
		for (CachedSphere s:spheres)
			mRenderer.addSphere(s.atom, s.bond, s.c, s.radius, s.argb);
		}

	private class CachedSphere implements Comparable<CachedSphere> {
		int atom,bond,argb;
		Coordinates c;
		double radius;

		public CachedSphere(int atom, int bond, Coordinates c, double radius, int argb) {
			this.atom = atom;
			this.bond = bond;
			this.c = new Coordinates(c);
			this.radius = radius;
			this.argb = argb;
			}

		@Override
		public int compareTo(CachedSphere o) {
			return (argb == o.argb) ? 0 : (argb < o.argb) ? -1 : 1;
			}
		}

	private class CachedCylinder implements Comparable<CachedCylinder> {
		int bond,argb;
		Coordinates c;
		double radius,length,rotationY,rotationZ;

		public CachedCylinder(int bond, double radius, double length, Coordinates c, double rotationY, double rotationZ, int argb) {
			this.bond = bond;
			this.radius = radius;
			this.length = length;
			this.c = new Coordinates(c);
			this.rotationY = rotationY;
			this.rotationZ = rotationZ;
			this.argb = argb;
			}

		@Override
		public int compareTo(CachedCylinder o) {
			return (argb == o.argb) ? 0 : (argb < o.argb) ? -1 : 1;
			}
		}
	}
