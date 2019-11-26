package com.actelion.research.chem.phesa;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.IntStream;
import com.actelion.research.chem.StereoMolecule;



/**
 * @author joel
 * Algorithm to create to align several active compounds in abscence of a bioactive conformation
 */
public class PheSABindingHypothesis {
	private static final int SOLUTIONS = 5; //report 5 best hypotheses
	private double[][] similarityMatrix;
	private List<StereoMolecule> actives;
	private int[] nConfs;
	private int nProcesses;
	
	public PheSABindingHypothesis(List<StereoMolecule> actives) {
		this(actives,1);
	}
	
	public PheSABindingHypothesis(List<StereoMolecule> actives, int nProcesses) {
		this.actives = actives;
		this.nProcesses = nProcesses;
	}
	
	public StereoMolecule[][] generate() {
		List<PheSAMolecule> phesaModels = new ArrayList<PheSAMolecule>();
		DescriptorHandlerShape dhs = new DescriptorHandlerShape();
		actives.forEach(e -> phesaModels.add(dhs.createDescriptor(e)));
		nConfs = new int[actives.size()];
		IntStream.range(0,actives.size()).forEach(e -> nConfs[e]=phesaModels.get(e).getVolumes().size());
		int totNrStrucs = Arrays.stream(nConfs).sum();
		similarityMatrix = new double[totNrStrucs][totNrStrucs];
		StereoMolecule[][] alignments = new StereoMolecule[SOLUTIONS][actives.size()];
		if(nProcesses==1)
			calculateSimilarityMatrix(phesaModels);
		else 
			calculateSimilarityMatrixParallel(phesaModels);
		int[][] bestEnsembles = getBestEnsembles();
		int counter = 0;
		for(int[] bestEnsemble : bestEnsembles) {
			int referenceConfIndex = bestEnsemble[bestEnsemble.length-1]; //index in similarity matrix
			int[] molIndexconfIndex = confIndexToMolIndex(referenceConfIndex);
			int molIndex = molIndexconfIndex[0];
			int confIndex = molIndexconfIndex[1]; 
			PheSAMolecule refModel = phesaModels.get(molIndex);
			MolecularVolume refVol = refModel.getVolumes().get(confIndex);
			StereoMolecule refMol = refModel.getConformer(refVol);
			StereoMolecule[] alignment = alignments[counter];
			alignment[0] = refMol;
			int i=0;
			for(int index : bestEnsemble) {
				if(index==referenceConfIndex)
					continue;
				else {
					i++;
					molIndexconfIndex = confIndexToMolIndex(index);
					molIndex = molIndexconfIndex[0];
					confIndex = molIndexconfIndex[1]; 
					PheSAMolecule fitModel = phesaModels.get(molIndex);
					MolecularVolume fitVol = fitModel.getVolumes().get(confIndex);
					StereoMolecule fitMol = fitModel.getConformer(fitVol);
					PheSAAlignmentOptimizer.align(fitModel, refVol, fitVol, fitMol);
					alignment[i] = fitMol;
				
				}
			}
			counter++;
		}

		return alignments;
	}
	
	private int[][] getBestEnsembles() {
		double[] bestScores = new double[SOLUTIONS];
		IntStream.range(0, bestScores.length).forEach(i -> bestScores[i]=-Double.MAX_VALUE);
		int[][] bestEnsembles = new int[SOLUTIONS][nConfs.length+1]; //the last one is the reference
		for(int refMol=0;refMol<actives.size();refMol++) {
			for(int refConf=0;refConf<nConfs[refMol];refConf++) {
				int indexRefConf=0;
				for(int m=0;m<refMol;m++) 
					indexRefConf+=nConfs[m];
				indexRefConf+=refConf;
				double score = 0.0;
				int[] ensemble = new int[nConfs.length+1];
				ensemble[refMol] = indexRefConf;
				ensemble[ensemble.length-1] = indexRefConf;
				for(int fitMol=0;fitMol<actives.size();fitMol++) {
					if(fitMol==refMol)
						continue;
					int indexStart=0;
					for(int n=0;n<fitMol;n++) 
						indexStart+=nConfs[n];
					int indexEnd = indexStart+nConfs[fitMol]; 
					int bestMatchingConfIndex = findIndexMaximumValueSubArray(similarityMatrix[indexRefConf],indexStart,indexEnd);
						score+=similarityMatrix[indexRefConf][bestMatchingConfIndex];

					ensemble[fitMol] = bestMatchingConfIndex;	
				}
				addToSolutions(bestEnsembles, ensemble, bestScores,score);
			}
		}
		return bestEnsembles;
	}
	
	private void addToSolutions(int[][] bestEnsembles, int[] ensemble, double[] bestScores, double score) {
		int n = bestScores.length;
		int rank = -1;
		if(score<bestScores[n-1])
			return;
		else {
			for(int i=0;i<bestScores.length;i++) {
				if(score>bestScores[i]) {
					rank = i;
					break;
				}
			}
			for(int j=bestScores.length-1;j>rank;j--) {
				bestScores[j] = bestScores[j-1];
				bestEnsembles[j] = bestEnsembles[j-1];
			}
			bestScores[rank] = score;
			bestEnsembles[rank] = ensemble;
				
		}
	}
					
	
	private void calculateSimilarityMatrix(List<PheSAMolecule> phesaModels) {
		for(int i=0;i<actives.size();i++) {
			for(int j=0;j<nConfs[i];j++) {
				MolecularVolume refVol = phesaModels.get(i).getVolumes().get(j);
				int index1=0;
				for(int m=0;m<i;m++) 
					index1+=nConfs[m];
				index1+=j;
				for(int k=i+1;k<actives.size();k++) {

					if(i==k) //conf from same molecule
						continue;
					for(int l=0;l<nConfs[k];l++) {
						MolecularVolume fitVol = new MolecularVolume(phesaModels.get(k).getVolumes().get(l));
						int index2=0;
						for(int n=0;n<k;n++) 
							index2+=nConfs[n];
						index2+=l;
						StereoMolecule fitMol = phesaModels.get(k).getConformer(fitVol); 
						similarityMatrix[index1][index2] = PheSAAlignmentOptimizer.align(phesaModels.get(k), refVol, fitVol, fitMol);
						similarityMatrix[index2][index1] = similarityMatrix[index1][index2];
					}
				}
			}
		}
	}
	
	private void calculateSimilarityMatrixParallel(List<PheSAMolecule> phesaModels) {
		ExecutorService executorService = Executors.newFixedThreadPool(nProcesses);
		for(int i=0;i<actives.size();i++) {
			for(int j=0;j<nConfs[i];j++) {
				MolecularVolume refVol = phesaModels.get(i).getVolumes().get(j);
				int index1=0;
				for(int m=0;m<i;m++) 
					index1+=nConfs[m];
				index1+=j;
				final int molIndex = i;
				final int confIndex = index1;
				executorService.execute(() -> compareOneConfWithAll(refVol,molIndex,confIndex,phesaModels));
			}
		}
		executorService.shutdown();

		while(!executorService.isTerminated()){
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
					
	
	
	private void compareOneConfWithAll(MolecularVolume refVol,int molIndex,int confIndex, List<PheSAMolecule> phesaModels) {
		for(int k=molIndex+1;k<actives.size();k++) {

			if(molIndex==k) //conf from same molecule
				continue;
			for(int l=0;l<nConfs[k];l++) {
				MolecularVolume fitVol = new MolecularVolume(phesaModels.get(k).getVolumes().get(l));
				int index2=0;
				for(int n=0;n<k;n++) 
					index2+=nConfs[n];
				index2+=l;
				StereoMolecule fitMol = phesaModels.get(k).getConformer(fitVol); 
				similarityMatrix[confIndex][index2] = PheSAAlignmentOptimizer.align(phesaModels.get(k), refVol, fitVol, fitMol);
				similarityMatrix[index2][confIndex] = similarityMatrix[confIndex][index2];
			}
		}
		
	}
	
	private  int findIndexMaximumValueSubArray(double[] arr, int startIncl, int endExcl) {
		double maxValue = -Double.MAX_VALUE;
		int maxIndex=-1;
		for(int i=startIncl;i<endExcl;i++) {
			if(arr[i]>maxValue) {
				maxValue = arr[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	
	/* returns the index of the molecule and the index of the conformer with 
	   *respect to this very molecule
	   */
	
	private int[] confIndexToMolIndex(int confIndex) { 												   
		int sum = 0;
		int indexMol = 0;
		int indexConf = 0;
		for(int i=0;i<nConfs.length;i++) {
			sum+=nConfs[i];
			if(confIndex<sum) {
				indexMol = i;
				indexConf = nConfs[i]-(sum-confIndex);
				return new int[] {indexMol,indexConf};
			}
		}
		return new int[] {indexMol,indexConf};
			
	}
	
	

}
