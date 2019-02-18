/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of ActelionMMFF94.
 * 
 * ActelionMMFF94 is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * ActelionMMFF94 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with ActelionMMFF94.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Paolo Tosco,Daniel Bergmann
 */

package mmff;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.AbstractForceField;

/**
 * The MMFF ForceField class is the top level class used to perform
 * energy calculations/minimisation on a molecule. It accepts an
 * ExtendedMolecule and the string name of the parameter tables to use.
 * Options can also be passed to control which energy terms should be
 * included in the force field equation (default: all 7 terms are
 * included), which non-bonded cutoff should be used (default: 100.0
 * angstrom), which dielectric constant should be used (default: 1.0),
 * and whether a constant or a distance-dependent dielectric model should
 * be used (default: linear).
 *
 * Tables: The Tables object contains a collection of parameter tables.
 * These are added statically to the ForceField class to allow all
 * ForceField instances to select which single instance of parameter
 * Tables to use. For example, there would be one tables instance for
 * MMFF94 and another tables instance for MMFF94s. Then when constructing
 * a ForceField object either MMFF94 or MMFF94s can be selected by passing
 * the string identifier used when loading that table with "loadTable()".
 * Options:
 *  - "nonbonded cutoff": A double, limits the number of Electrostatic and
 *      VanDerWaals energy terms to those involving atoms not further
 *      apart than "nonbonded cutoff" (default: 100.0 angstrom).
 *  - "dielectric constant": A double, the dielectric constant value
 *      (default: 1.0)
 *  - "dielectric model": A string for the dielectric model. "distance"
 *      selects a distance-dependent dielectric model, everything else
*       selects the constant dielectric model (default: "constant")
 *  - "angle bend": A boolean, default True, for whether to include angle
 *      bending energy terms.
 *  - "bond stretch": A boolean, default True, for whether to include bond
 *      stretching energy terms.
 *  - "electrostatic": A boolean, default True, for whether to include the
 *      nonbonded electrostatic energy terms.
 *  - "out of plane": A boolean, default True, for whether to include out
        of plane energy terms.
 *  - "stretch bend": A boolean, default True, for whether to include
        stretch bending energy terms.
 *  - "torsion angle": A boolean, default True, for whether to include
 *      torsional angle energy terms.
 *  - "van der waals": A boolean, default True, for whether to include the
 *      nonbonded van der Waals energy terms.
 */
public final class ForceFieldMMFF94 extends AbstractForceField {
	public static final String MMFF94 = "MMFF94";
	public static final String MMFF94S = "MMFF94s";
    public static final String MMFF94SPLUS = "MMFF94s+";


    private final MMFFMolecule mMMFFMol;
    public static Map<String, Tables> mTables = new HashMap<String, Tables>();
    private List<EnergyTerm> mEnergies = new ArrayList<EnergyTerm>();
    private RingBoolean[] mRingArom;
    private int[] mAtomTypes;
    private int[] mHydrogenMap;
    
    /**
     * Forcefield constructor.
     *  @param m The molecule to construct the forcefield on.
     *  @param tablename The string name for the Tables to be used. There
     *      must be a table with this name that has been loaded with
     *      "loadTable()".
     *  @param options A Map containing the ForceField options and values.
     *      See class description of a list of options.
     */
    
    public ForceFieldMMFF94(StereoMolecule m, String tablename,
            Map<String, Object> options) {
    	super(m);
    	mMMFFMol = new mmff.MMFFMolecule(m);
        Tables table = mTables.get(tablename);

        double nonBondedThresh = options.containsKey("nonbonded cutoff")
            ? ((Double)options.get("nonbonded cutoff")).doubleValue()
            : 100.0;

        double dielConst = options.containsKey("dielectric constant")
            ? ((Double)options.get("dielectric constant")).doubleValue() : 1.0;

        boolean dielModel = options.containsKey("dielectric model")
            ? ((String)options.get("dielectric model")).equals("distance")
            : false;



        Separation sep = new Separation(mMMFFMol);

        if (!options.containsKey("angle bend")
                || (Boolean)options.get("angle bend"))
            mEnergies.addAll(AngleBend.findIn(table, mMMFFMol));

        if (!options.containsKey("bond stretch")
                || (Boolean)options.get("bond stretch"))
        	mEnergies.addAll(BondStretch.findIn(table, mMMFFMol));

        if (!options.containsKey("electrostatic")
                || (Boolean)options.get("electrostatic"))
        	mEnergies.addAll(Electrostatic.findIn(table, mMMFFMol, sep,
                        nonBondedThresh, dielModel, dielConst));

        if (!options.containsKey("out of plane")
                || (Boolean)options.get("out of plane"))
        	mEnergies.addAll(OutOfPlane.findIn(table, mMMFFMol));

        if (!options.containsKey("stretch bend")
                || (Boolean)options.get("stretch bend"))
        	mEnergies.addAll(StretchBend.findIn(table,mMMFFMol));

        if (!options.containsKey("torsion angle")
                || (Boolean)options.get("torsion angle"))
        	mEnergies.addAll(TorsionAngle.findIn(table, mMMFFMol));

        if (!options.containsKey("van der waals")
                || (Boolean)options.get("van der waals"))
        	mEnergies.addAll(VanDerWaals.findIn(table, mMMFFMol, sep, nonBondedThresh));
    }

    /**
     * Forcefield constructor. Overloaded to pass the default (empty)
     * options to a ForceField.
     *  @param mol The molecule to construct the forcefield on.
     *  @param tablename The string name for the Tables to be used. There
     *      must be a table with this name that has been loaded with
     *      "loadTable()".
     */
    public ForceFieldMMFF94(StereoMolecule mol, String tablename) {
        this(mol, tablename, new HashMap<String, Object>());
    }

    /**
     * Returns the total number of atoms in this force field.
     *  @return Total number of atoms.
     */
    public int size() {
        return mMMFFMol.getAllAtoms();
    }


    /**
     * Minimise the current molecule using default parameter values for
     * the number of iterations, energy tolerance and gradient tolerance.
     *  @return Return code, 0 on success.
     */
    public int minimise() {
        return minimise(200, 1e-4, 1e-6);
    }
    
    

    @Override
    public double updateGradient() {
        mGrad = new double[mDim];

        for (EnergyTerm engy : mEnergies) 
            engy.getGradient(mPos, mGrad);
        double maxGrad = -1e8;
        double gradScale = 0.1;
        for (int i=0; i<mDim; i++) {
            mGrad[i] *= gradScale;
            if (mGrad[i] > maxGrad)
                maxGrad = mGrad[i];
        }

        if (maxGrad > 10.0) {
            while (maxGrad*gradScale > 10.0)
                gradScale *= 0.5;

            for (int i=0; i<mDim;i++)
                mGrad[i] *= gradScale;
        }
        return gradScale;
    }

    /**
     * Gets the total energy of the molecule as the sum of the energy terms.
     *  @return The total force field energy.
     */
    @Override
    public double[] getCurrentPositionsMapped() {
    	int[] atomMap = mMMFFMol.getHydrogenMap();
    	double[] pos = Arrays.copyOf(mPos, mPos.length);
    	for(int i=0;i<atomMap.length;i++) {
    		pos[3*i] = mPos[3*atomMap[i]];
    		pos[3*i+1] = mPos[3*atomMap[i]+1];
    		pos[3*i+2] = mPos[3*atomMap[i]+2];
    		
    	}
    	return pos;
    	
    }
    
    @Override
    public double[] getCurrentPositions() {
    	return Arrays.copyOf(mPos, mPos.length);
    }
    
    
    
    public double getTotalEnergy(double[] pos) {
        double total = 0.0;
        for (EnergyTerm term : mEnergies)
            total += term.getEnergy(pos);

        return total;
    }

    /**
     * Gets the total energy of the molecule as the sum of the energy
     * terms. This function passes the force fields `pos` array to
     * getTotalEnergy().
     *  @return The total force field energy.
     */
    public double getTotalEnergy() {
        return getTotalEnergy(mPos);
    }

    public static void initialize(String tableSet) {
        ForceFieldMMFF94.loadTable(tableSet, Tables.newMMFF94(tableSet));
    }

    /**
     * Loads and registers a tables object with the ForceField class so it
     * can be used by new ForceField instances.
     *  @param name The string name used to identifiy the tables object.
     *  @param table The tables object.
     */
    public static synchronized void loadTable(String name, Tables table) {
        if (!mTables.containsKey(name))
            mTables.put(name, table);
    }

    /**
     * Returns a table given a table name.
     *  @param name The string name of a table.
     *  @return The tables object.
     */
    public static Tables table(String name) {
        return mTables.get(name);
    }


	
	public MMFFMolecule getMMFFMolecule() {
		return mMMFFMol;
	}
}
