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

package com.actelion.research.chem.properties.fractaldimension;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.CommandLineParser;
import com.actelion.research.util.IO;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * FractalDimensionMoleculeMain
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 27.08.18.
 */
public class FractalDimensionMoleculeMain {

    public static final String PREFIX_OUTPUT = "fractalDimensionMoleculeResult";

    public static final int TOTAL_CAPACITY = (int)(30*Math.pow(10,6));

    private static final String USAGE =

        "FractalDimensionMoleculeMain.\n" +
            "Modest von Korff, Thomas Sander.\n" +
            "2018.\n" +
            "Calculates the fractal dimension for molecules. The " +
                "application is memory demanding. 1 GB RAM should be " +
                "assigned. Input is a file with a list of SMILES strings. Output is a tab separated file with a header " +
                "line. One line is produced for each SMILES string in the result file. \n" +
            "-h this help text.\n" +
            "-i String path to SMILES input file. Text file with with SMILES, separated by newline.\n" +
            "-w String path to output directory.\n" +
            "[-c] integer maximum capacity times 1000 for fragments, default is " + TOTAL_CAPACITY + " .\n" +
            "[-D] boolean elusive output" + ".\n" +
            "\n";


    public static void main(String[] args) throws Throwable {

        CommandLineParser clp = new CommandLineParser();

        clp.parse(args);

        if (clp.getNumArguments()==0 || clp.help()) {
            System.out.println(USAGE);
            System.exit(0);
        }



        clp.checkCommandWithValue("-i");

        File fiMolecule = new File(clp.get("-i"));

        if (!fiMolecule.isFile()) {
            System.err.println("Not a file " + fiMolecule.getAbsolutePath() + ".");
            System.out.println(USAGE);
            System.exit(0);
        }

        clp.checkCommandWithValue("-w");

        File workdir = new File(clp.get("-w"));
        if (!workdir.isDirectory()) {
            System.err.println("Not a directory '" + workdir.getAbsolutePath() + "'.");
            System.out.println(USAGE);
            System.exit(0);
        }


        File fiOutTbl = new File(workdir, PREFIX_OUTPUT + ".txt");

        fiOutTbl.createNewFile();

        if(!fiOutTbl.canWrite()){
            throw new IOException("No permission to write output file " + fiOutTbl.getAbsolutePath() + ".");
        }

        List<String> liSMILES = IO.readLines2List(fiMolecule);

        List<InputObjectFracDimCalc> liMolecule = new ArrayList<InputObjectFracDimCalc>();

        SmilesParser smilesParser = new SmilesParser();


        System.out.println("Parsing SMILES");

        for (int i = 0; i < liSMILES.size(); i++) {

            String smile = liSMILES.get(i);

            try {
                StereoMolecule mol = new StereoMolecule();

                smilesParser.parse(mol, smile);

                mol.ensureHelperArrays(Molecule.cHelperRings);

                InputObjectFracDimCalc inputObjectFracDimCalc = new InputObjectFracDimCalc(mol, i, smile);

                liMolecule.add(inputObjectFracDimCalc);

            } catch (Exception e) {
                System.err.println("SMILES parsing error in line " + i + " for " + smile);
                e.printStackTrace();
            }
        }

        System.out.println("Parsed " + liSMILES.size() + " SMILES. Succeeded for " + liMolecule.size() + ".");

        int totalCapacity = TOTAL_CAPACITY;

        if (clp.contains("-c")) {

            String s = clp.get("-c");

            if(s==null){
                System.err.println("No capacity value given for -c.");
            } else {
                totalCapacity = Integer.parseInt(s) * 1000;
            }
        }

        System.out.println("Total capacity " + totalCapacity);

        boolean elusive = false;
        if(clp.contains("-D")){
            elusive = true;
            System.out.println("Elusive output.");
        }

        FractalDimensionMolecule fractalDimensionMolecule = new FractalDimensionMolecule(totalCapacity, elusive);

        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fiOutTbl), StandardCharsets.UTF_8));

        ResultFracDimCalcHeaderTags resultFracDimCalcHeaderTags = new ResultFracDimCalcHeaderTags();

        bw.write(resultFracDimCalcHeaderTags.toStringHeader());

        for (InputObjectFracDimCalc inputObjectFracDimCalc : liMolecule) {

            ResultFracDimCalc result = fractalDimensionMolecule.process(inputObjectFracDimCalc);

            bw.write("\n");
            bw.write(result.toString());
            bw.flush();
        }

        bw.close();

        fractalDimensionMolecule.finalizeThreads();

        System.out.println("Wrote output into " + fiOutTbl.getAbsolutePath() + ".");

        System.out.println("Finished");

    }




}