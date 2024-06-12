package com.actelion.research.chem.descriptor.flexophore.generator;

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
    public static final long TIMEOUT_CONFORMER_CALCULATION_MS = TimeDelta.MS_MINUTE * 3;
    public static final long SEED = 123456789;
    private static final int MAX_TORSION_SETS = 100000;
    private static final int MAX_TRIES_CONFORMERS = 10;

    private static final int MAX_INITIALIZATION_STAGE = 3;
    private ConformerGenerator conformerGenerator;
    private int initializationStage;

    private long seed;

    private Exception recentException = null;

    private int ccGeneratedConformers;

    private Molecule3D molInPlace;

    public void resetInitializationStage() {
        initializationStage = 0;
    }

    public ConformerGeneratorStageTries() {
        seed = SEED;
        conformerGenerator = new ConformerGenerator(seed, false);
        conformerGenerator.setTimeOut(TIMEOUT_CONFORMER_CALCULATION_MS);
        RigidFragmentCache.getDefaultInstance().loadDefaultCache();
        initializationStage = 0;
    }

    private void initializeHelper() {
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
        return initializeConformers();
    }

    /**
     * setMolecule(Molecule3D molInPlace) first!
     * @return
     */
    private boolean initializeConformers(){
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
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_LIKELY_SYSTEMATIC, MAX_TORSION_SETS, true);
                } else if (initializationStage == 3) {
                    conformerGenerator = new ConformerGenerator();
                    successfulInitialization = conformerGenerator.initializeConformers(molInPlace, ConformerGenerator.STRATEGY_ADAPTIVE_RANDOM, MAX_TORSION_SETS, true);
                } else if (initializationStage > MAX_INITIALIZATION_STAGE) {
                    break;
                }
            } catch (Exception e) {
                exception = e;
            }

            if (!successfulInitialization) {
                System.out.println("Initialization failed for stage " + initializationStage);

                if(canIncrementInitializationStage()) {
                    if(!incrementInitializationStage())
                        break;
                    }
                else {
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
