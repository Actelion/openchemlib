package com.actelion.research.chem.descriptor.flexophore.generator;

/*
 
 Copyright (c) 2024 Alipheron AG. All rights reserved.
 
 This file is part of the Alipheron AG software suite.
 
 Licensed under the Alipheron AG Software License Agreement (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at the company's official website or upon request.
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License. 
 
 Created by Modest von Korff 
 21/02/2024
 
 */

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.util.TimeDelta;
import org.openmolecules.chem.conf.gen.ConformerGenerator;
import org.openmolecules.chem.conf.gen.RigidFragmentCache;

import java.util.Date;

/**
 * Tries different initializations for the conformation generator.
 */
public class ConformerGeneratorStageTries {
    public static final long TIMEOUT_CONFORMER_CALCULATION_MS = TimeDelta.MS_MINUTE * 6;
    public static final long SEED = 123456789;
    private static final int MAX_TORSION_SETS = 100000;
    private static final int MAX_TRIES_CONFORMERS = 10;

    private static final int MAX_INITIALIZATION_STAGE = 5;
    private ConformerGenerator conformerGenerator;
    private int initializationStage;

    private long seed;
    private long t0ConformerCalcStarted;
    private Exception recentException = null;

    private int ccGeneratedConformers;

    private Molecule3D molInPlace;

    public void resetInitializationStage() {
        initializationStage = 0;
    }

    public ConformerGeneratorStageTries() {
        seed = SEED;
        conformerGenerator = new ConformerGenerator(seed, false);
        RigidFragmentCache.getDefaultInstance().loadDefaultCache();
        initializationStage = 0;
        t0ConformerCalcStarted = 0;
    }

    private void initializeHelper() {
        t0ConformerCalcStarted = System.currentTimeMillis();
        ccGeneratedConformers = 0;
    }

    public void setThreadMaster(ThreadMaster threadMaster) {
        conformerGenerator.setThreadMaster(threadMaster);
    }

    public boolean incrementInitializationStage() {
        initializationStage++;
        if (initializationStage > MAX_INITIALIZATION_STAGE) {
            throw new RuntimeException("Maximum initialization stage exceeded!");
        }
        return initializeConformers();
    }

    public boolean canIncrementInitializationStage() {
        if (initializationStage < MAX_INITIALIZATION_STAGE) {
            return true;
        }
        return false;
    }

    public int getPotentialConformerCount() {
        return conformerGenerator.getPotentialConformerCount();
    }

    public boolean setMolecule(Molecule3D molInPlace) {
        this.molInPlace = molInPlace;
        resetInitializationStage();
        boolean successfulInitialization = initializeConformers();
        return successfulInitialization;
    }

    /**
     * setMolecule(Molecule3D molInPlace) first!
     * @return
     */
    public boolean initializeConformers(){
        initializeHelper();
        boolean successfulInitialization = false;
        Exception exception = null;

        while (!successfulInitialization) {

            try {
                if (initializationStage == 0) { // default
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_TORSION_SETS, false);
                } else if (initializationStage == 1) {
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_TORSION_SETS, true);
                }  else if (initializationStage == 2) {
                    RigidFragmentCache.getDefaultInstance().clear();
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_TORSION_SETS, true);
                }  else if (initializationStage == 3) {
                    RigidFragmentCache.getDefaultInstance().clear();
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_RANDOM, MAX_TORSION_SETS, true);
                } else if (initializationStage == 4) {
                    RigidFragmentCache.getDefaultInstance().clear();
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_SYSTEMATIC, MAX_TORSION_SETS, true);
                } else if (initializationStage == 5) {
                    RigidFragmentCache.getDefaultInstance().clear();
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_ADAPTIVE_RANDOM, MAX_TORSION_SETS, true);
                } else if (initializationStage > MAX_INITIALIZATION_STAGE) {
                    break;
                }
            } catch (Exception e) {
                exception = e;
                successfulInitialization=false;
            }

            if (!successfulInitialization) {
                System.out.println("Initialization failed for stage " + initializationStage);

                if(canIncrementInitializationStage())
                    if(!incrementInitializationStage()){
                        successfulInitialization = false;
                        break;
                    }
                else {
                    successfulInitialization = false;
                    break;
                }
            }

//            if(initializationStage>0)
//                System.out.println("CreatorMolDistHistViz initialization stage " + initializationStage);

        }

        if(!successfulInitialization){
            System.err.println("CreatorMolDistHistViz initializeConformers(...) failed for " + molInPlace.getIDCode());
            recentException = exception;
        }

        return successfulInitialization;
    }

    public Exception getRecentException() {
        return recentException;
    }

    public boolean generateConformerAndSetCoordinates(int nAtoms, Molecule3D molInPlace){

        boolean nextConformerAvailable = false;

        Conformer conformer = conformerGenerator.getNextConformer();

        int ccTries=0;
        while (conformer==null && ccGeneratedConformers==0) {
            if(ccTries>MAX_TRIES_CONFORMERS){
                if(canIncrementInitializationStage()){
                    // System.out.println("incrementInitializationStage");
                    incrementInitializationStage();
                    conformer = conformerGenerator.getNextConformer();
                    ccTries=0;
                } else {
                    // System.out.println("Cannot incrementInitializationStage, break!");
                    break;
                }
            } else {
                // System.out.println("Inject new seed, ccTries " + ccTries);
                injectNewSeed();
                conformer = conformerGenerator.getNextConformer();
                ccTries++;
            }

        }

        if(conformer != null){
            ccGeneratedConformers++;
            // System.out.println(ConformerUtils.toStringDistances(ConformerUtils.getDistanceMatrix(conformer)));

            for (int i = 0; i < nAtoms; i++) {
                double x = conformer.getX(i);
                double y = conformer.getY(i);
                double z = conformer.getZ(i);

                molInPlace.setAtomX(i,x);
                molInPlace.setAtomY(i,y);
                molInPlace.setAtomZ(i,z);
            }
            nextConformerAvailable = true;
        } else {
            // System.out.println("Conformer null, ccGeneratedConformers " + ccGeneratedConformers);
        }

        return nextConformerAvailable;
    }
    public void injectNewSeed(){
        seed = new Date().getTime();
        conformerGenerator = new ConformerGenerator(seed, false);
        initializeConformers();
    }
}
