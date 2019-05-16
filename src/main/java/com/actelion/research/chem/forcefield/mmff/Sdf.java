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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import com.actelion.research.chem.forcefield.mmff.Vector3;

import java.io.File;

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
