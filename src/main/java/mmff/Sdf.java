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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;

import mmff.Vector3;

/**
 * The Sdf class provides helper functions for dealing with the SDF file
 * format.
 */
public class Sdf {
    public interface OnMolecule {
        boolean check(String sdfpath, String refpath);
    }

    /**
     * Tests an entire pair of folders, one containing SDF files and the
     * other containing reference files. Each folder should contain either
     * all SDF files or all reference files. There should be a
     * corresponding reference file for each SDF file, and they should
     * have the same name (excluding extension).
     *  @param sdfbase The base folder which contains all of the SDF
     *      files.
     *  @param refbase The base folder which contains all of the reference
     *      files.
     */
    public static void testFolder(String sdfbase, String refbase, String ext,
            OnMolecule cb) {
        int failed = 0;
        int total = 0;
        File dir = new File(sdfbase);
        File[] ls = dir.listFiles();

        for (File fp : ls) {
            String sdfpath = fp.getAbsolutePath();

            if (!sdfpath.endsWith(".sdf"))
                continue;

            String refpath = new File(refbase, fp.getName()).getPath()
                .replace(".sdf", ext);
            failed += cb.check(sdfpath, refpath) ? 0 : 1;
            total++;
        }

        System.out.println("Failed molecules: "+failed+" / "+total+"    ("
                +((failed*100.0)/total)+"%)");
    }

    /**
     * Parses a SDF file to extract the positions of each atom. Positions
     * are extracted as an (X, Y, Z) Vector3 object.
     *  @param sdfpath The string path to the SDF file to be loaded.
     *  @return Returns an array of Vector3's which contain the X, Y, Z
     *      coordinates of each atom.
     */
    public static Vector3[] getPositions(String sdfpath) {
        BufferedReader br = null;

        try {
            String line;
            br = new BufferedReader(new FileReader(sdfpath));

            // We don't need these lines.
            br.readLine();
            br.readLine();
            br.readLine();

            int natoms = Integer.parseInt(br.readLine().substring(0, 3).trim());
            Vector3[] res = new Vector3[natoms];

            for (int i=0; i<natoms && (line = br.readLine()) != null; i++) {
                res[i] = new Vector3(
                    Double.parseDouble(line.substring( 0,10).trim()),
                    Double.parseDouble(line.substring(10,20).trim()),
                    Double.parseDouble(line.substring(20,30).trim())
                );
            }
            return res;
        } catch (FileNotFoundException e) {
            System.out.println("Could not find file: '"+sdfpath+"'");
        } catch (IOException e) {
            System.out.println("IO Exception!");
        } finally {
            try {
                if (br != null)
                    br.close();
            } catch (IOException e) {
                System.out.println("Couldn't close buffered reader!");
            }
        }

        return new Vector3[1];
    }
}
