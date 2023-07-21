package com.actelion.research.chem.mcs;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.properties.complexity.FragmentDefinedByBondsIdCode;
import com.actelion.research.chem.properties.complexity.IBitArray;
import com.actelion.research.util.Pipeline;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

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
public class RunBondVector2IdCode implements Runnable {

    private static final long SLEEP = 10;

    private int id;

    private Pipeline<IBitArray> pipeInputFragIndexListsFromEFG;

    private Pipeline<FragmentDefinedByBondsIdCode> pipeOutputFragmentDefinedByBondsIdCode;

    private AtomicBoolean endOfRunReached;

    private BondVector2IdCode bondVector2IdCode;

    private AtomicInteger processedFragments;

    public RunBondVector2IdCode(int id, Pipeline<IBitArray> pipeInputFragIndexListsFromEFG, Pipeline<FragmentDefinedByBondsIdCode> pipeOutputFragmentDefinedByBondsIdCode) {

        this.id = id;

        this.pipeInputFragIndexListsFromEFG = pipeInputFragIndexListsFromEFG;

        this.pipeOutputFragmentDefinedByBondsIdCode = pipeOutputFragmentDefinedByBondsIdCode;

        endOfRunReached = new AtomicBoolean(false);
    }

    public void init(StereoMolecule mol){
        StereoMolecule molCopy = new StereoMolecule(mol);
        molCopy.ensureHelperArrays(Molecule.cHelperCIP);
        bondVector2IdCode = new BondVector2IdCode(molCopy);
        processedFragments = new AtomicInteger();
    }


    public void run() {

        try {

            while(!pipeInputFragIndexListsFromEFG.wereAllDataFetched()){
                IBitArray livIndexBond = pipeInputFragIndexListsFromEFG.pollData();
                if(livIndexBond == null){
                    try {Thread.sleep(SLEEP);} catch (InterruptedException e) {e.printStackTrace();}
                    continue;
                }

                try {
                    FragmentDefinedByBondsIdCode livIdCode = new FragmentDefinedByBondsIdCode(livIndexBond);
                    String idcodeFrag = bondVector2IdCode.getFragmentIdCode(livIndexBond);
                    livIdCode.setIdCode(idcodeFrag);
                    pipeOutputFragmentDefinedByBondsIdCode.addData(livIdCode);

                } catch (Throwable e) {
                    e.printStackTrace();
                } finally {
                    processedFragments.incrementAndGet();
                }

            } // End while

        } catch (Exception e){
          e.printStackTrace();
        } finally {
            // System.out.println("RunBondVector2IdCode finally reached.");
            endOfRunReached.set(true);
        }

    } // End run()

    public boolean isEndOfRunReached() {
        return endOfRunReached.get();
    }

}




