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
*/

package com.actelion.research.chem.forcefield.mmff;

import com.actelion.research.chem.forcefield.mmff.table.*;

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
    public final com.actelion.research.chem.forcefield.mmff.table.OutOfPlane  oop;      // Out of plane.
    public final com.actelion.research.chem.forcefield.mmff.table.Charge      chge;     // Charge.
    public final Stbn                   stbn;     // Stretch bending.
    public final Torsion                torsion;  // Torsional angles.
    public final com.actelion.research.chem.forcefield.mmff.table.VanDerWaals vdws;     // Van der Waals.

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
        chge     = new      com.actelion.research.chem.forcefield.mmff.table.Charge(this, csv_pbci, csv_bci);
        atom     = new                   Atom(this, csv_atom);
        bndk     = new                   Bndk(this, csv_bndk);
        bond     = new                   Bond(this, csv_bond);
        covrad   = new                 CovRad(this, csv_covrad);
        dfsb     = new                   Dfsb(this, csv_dfsb);
        def      = new                    Def(this, csv_def);
        hblaurie = new       HerschbachLaurie(this, csv_hblaurie);
        oop      = new  com.actelion.research.chem.forcefield.mmff.table.OutOfPlane(this, csv_oop);
        stbn     = new                   Stbn(this, csv_stbn);
        torsion  = new                Torsion(this, csv_torsion);
        vdws     = new com.actelion.research.chem.forcefield.mmff.table.VanDerWaals(this, csv_vdws);
    }

    /**
     * Returns a new MMFF94 table. The paths provided here are correct
     * when running the project using ant/eclipse or from the same
     * directory as build.xml and the README.
     */
    public static Tables newMMFF94(String tableSet) {
        return new com.actelion.research.chem.forcefield.mmff.Tables(
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
