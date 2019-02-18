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

import mmff.table.*;

/**
 * A collection of tables with the official parameters of the
 * MMFF94/MMFF94s force field variants.
 */
public final class Tables {
	private static final String RESOURCE_PATH = "/resources/forcefield/mmff94/";

	public final Angle                  angle;    // Angle bending.
    public final Atom                   atom;     // Atom.
    public final Bndk                   bndk;     // Bond stretching.
    public final Bond                   bond;     // Bond stretching.
    public final CovRad                 covrad;   // Bond stretching.
    public final Dfsb                   dfsb;     // Bond stretching.
    public final Def                    def;      // Equivalent types table.
    public final HerschbachLaurie       hblaurie; // Stretch bending.
    public final mmff.table.OutOfPlane  oop;      // Out of plane.
    public final mmff.table.Charge      chge;     // Charge.
    public final Stbn                   stbn;     // Stretch bending.
    public final Torsion                torsion;  // Torsional angles.
    public final mmff.table.VanDerWaals vdws;     // Van der Waals.

    /**
     * Construct a new Tables object. Takes several string paths to CSV
     * files to be used for the different parameter tables.
     *  @param csv_angle Path to the angles parameter table. The angle
     *      table holds ideal angle and force constant values.
     *  @param csv_atom Path to the atoms parameter table. This table
     *      specifies the chemical, topological and geometrical properties
     *      associated with each of the MMFF atom types.
     *  @param csv_bci Path to the bond charge increments table. This
     *      table is used for the MMFF partial charge computation.
     *  @param csv_bndk Path to the bond length and default force
     *      constants table. Used as fallback values if bond parameters
     *      are not found in csv_bond.
     *  @param csv_bond Path to the bond parameters table. Holds the ideal
     *      length and force constants for bonds between MMFF atom typed
     *      atoms.
     *  @param csv_covrad Path to the CovRad CSV table.
     *  @param csv_dfsb Path to the Dfsb CSV table. Holds the default
     *      stretch bend parameters.
     *  @param csv_def Path to the default types table. Holds fallback
     *      MMFF types for atoms. Used to help find matches for energy
     *      terms in their parameter tables.
     *  @param csv_hblaurie Path to the HerschbachLaurie table. Parameters
     *      table for Badger's rule.
     *  @param csv_oop Path to the out of plane parameters table. Holds
     *      the force constant parameter.
     *  @param csv_pbci Path to the partial bond charge increments and
     *      formal charge adjustments table.
     *  @param csv_stbn Path to the stretch bend parameters table. Holds
     *      the force constants for both i-j-k and k-j-i atom
     *      configurations.
     *  @param csv_torsion Path to the torsional parameters table. Holds
     *      the three force constants for given MMFF atom types.
     *  @param csv_vdws Path to the van der Waals parameters table.
     */
    public Tables(
                String csv_angle,
                String csv_atom,
                String csv_bci,
                String csv_bndk,
                String csv_bond,
                String csv_covrad,
                String csv_dfsb,
                String csv_def,
                String csv_hblaurie,
                String csv_oop,
                String csv_pbci,
                String csv_stbn,
                String csv_torsion,
                String csv_vdws) {

        angle    = new                  Angle(this, csv_angle);
        chge     = new      mmff.table.Charge(this, csv_pbci, csv_bci);
        atom     = new                   Atom(this, csv_atom);
        bndk     = new                   Bndk(this, csv_bndk);
        bond     = new                   Bond(this, csv_bond);
        covrad   = new                 CovRad(this, csv_covrad);
        dfsb     = new                   Dfsb(this, csv_dfsb);
        def      = new                    Def(this, csv_def);
        hblaurie = new       HerschbachLaurie(this, csv_hblaurie);
        oop      = new  mmff.table.OutOfPlane(this, csv_oop);
        stbn     = new                   Stbn(this, csv_stbn);
        torsion  = new                Torsion(this, csv_torsion);
        vdws     = new mmff.table.VanDerWaals(this, csv_vdws);
    }

    /**
     * Returns a new MMFF94 table. The paths provided here are correct
     * when running the project using ant/eclipse or from the same
     * directory as build.xml and the README.
     */
    public static Tables newMMFF94(String tableSet) {
        return new mmff.Tables(
        		RESOURCE_PATH + "angle.csv",
        		RESOURCE_PATH + "atom.csv",
        		RESOURCE_PATH + "bci.csv",
        		RESOURCE_PATH + "bndk.csv",
        		RESOURCE_PATH + "bond.csv",
        		RESOURCE_PATH + "covrad.csv",
        		RESOURCE_PATH + "dfsb.csv",
        		RESOURCE_PATH + "def.csv",
        		RESOURCE_PATH + "herschbachlaurie.csv",
        		RESOURCE_PATH + (tableSet.equals(ForceFieldMMFF94.MMFF94S) || tableSet.equals(ForceFieldMMFF94.MMFF94SPLUS)?"94s/outofplane.csv":"outofplane.csv"),
        		RESOURCE_PATH + "pbci.csv",
        		RESOURCE_PATH + "stbn.csv",
        		RESOURCE_PATH + (tableSet.equals(ForceFieldMMFF94.MMFF94S)?"94s/torsion.csv":tableSet.equals(ForceFieldMMFF94.MMFF94SPLUS)?"94s/torsionPlus.csv":"torsion.csv"),
        		RESOURCE_PATH + "vanderwaals.csv");
    }
}
