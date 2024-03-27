package com.actelion.research.chem.chemicalspaces;

import com.actelion.research.chem.*;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonCreator;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.descriptor.DescriptorHandlerLongFFP512;
import com.actelion.research.chem.io.DWARFileCreator;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.Reactor;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.stream.Collectors;


public class ChemicalSpaceCreator {

	private File outdirectory;
	private Set<String> bbs;
	private Map<String,Map<String,List<String>>> bbData;
	private List<Reaction> reactions;
	private final List<Reaction> functionalizations;

	public ChemicalSpaceCreator(Set<String> bbs, List<Reaction> reactions, File outdirectory) {
		this.bbs = bbs;
		this.reactions = reactions;
		this.outdirectory = outdirectory;
		this.functionalizations = new ArrayList<>(); //reactions that only consist of modifications ("functionalizations") of one building block are treated specially
		for(Reaction rxn : reactions ) {
			if(rxn.getReactants()==1)
				functionalizations.add(rxn);
		}
		reactions.removeAll(functionalizations);
	}
	
	public void setOutdirectory(File outdirectory) {
		this.outdirectory = outdirectory;
	}

	public void setBBs(Set<String> bbs) {
		this.bbs = bbs;
	}

	public void setBBData(Map<String, Map<String, List<String>>> bbData) {
		this.bbData = bbData;
	}

	public void setReactions(List<Reaction> reactions) {
		this.reactions = reactions;
	}
	
	public void create() {
		Map<String,List<Reaction>> allSynthonTransformations = new HashMap<String,List<Reaction>>();
		generateSynthonTransformations(reactions,allSynthonTransformations);
		ConcurrentMap<String,String> processedToOrigIDCode = new ConcurrentHashMap<String,String>();
		ConcurrentMap<String,List<Map<String,String>>> reactionsWithSynthons = new ConcurrentHashMap<String,List<Map<String,String>>>(); 
		processBuildingBlocks(this.bbs,processedToOrigIDCode,functionalizations);
		ConcurrentMap<String, long[]> fps = new ConcurrentHashMap<>();
		calcFragFPs(processedToOrigIDCode.keySet(), fps);
		generateSynthons(reactions, processedToOrigIDCode, reactionsWithSynthons, fps,allSynthonTransformations);
		generateCombinatoriaLibraries(reactionsWithSynthons, allSynthonTransformations);
	}
	
	private static void calcFragFPs(Collection<String> idcodes, ConcurrentMap<String,long[]> fps) {
		idcodes.parallelStream().forEach(idc -> {
			IDCodeParser parser = new IDCodeParser();
			DescriptorHandlerLongFFP512 dhf = new DescriptorHandlerLongFFP512();
			StereoMolecule mol = new StereoMolecule();
			try {
				parser.parse(mol, idc);
				long[] desc = dhf.createDescriptor(mol);
				fps.put(idc, desc);
			}
			catch(Exception e) {}
		});
	}
	
	private static void processBuildingBlocks(Collection<String> bbs, ConcurrentMap<String,String> processedToOrigIDCode, List<Reaction> functionalizations) {
		bbs.parallelStream().forEach( idcode -> {
			StereoMolecule mol = new IDCodeParser().getCompactMolecule(idcode);
			if (mol != null) {
				mol.ensureHelperArrays(Molecule.cHelperCIP);
				mol.stripSmallFragments();
				mol.normalizeAmbiguousBonds();
				if(processedToOrigIDCode.get(mol.getIDCode())!=null)
					return;
				processedToOrigIDCode.put(mol.getIDCode(), idcode);
				
				for(Reaction functionalization : functionalizations) {
					SSSearcher searcher = new SSSearcher();
					searcher.setFragment(functionalization.getReactant(0));
					searcher.setMolecule(mol);
					
					if (searcher.isFragmentInMolecule()) {
						StereoMolecule product = getProduct(functionalization,Arrays.asList(mol));
						if(product!=null)
							processedToOrigIDCode.put(product.getIDCode(), idcode);
					}
				}
			}
		});
	}
	
	private static void generateSynthons(List<Reaction> reactions,
		ConcurrentMap<String,String> processedBBToBB,
		ConcurrentMap<String,List<Map<String,String>>> reactionsWithSynthons,ConcurrentMap<String,long[]> fps,
		Map<String,List<Reaction>> allSynthonTransformations) {

		for(Reaction rxn : reactions) {
			processReaction(rxn, processedBBToBB, reactionsWithSynthons, fps, allSynthonTransformations);
        }
  	}
	
	private static void generateSynthonTransformations(List<Reaction> reactions, Map<String,List<Reaction>> allSynthonTransformations) {
		for(Reaction rxn : reactions) {
			List<Reaction> synthonTransformations = new ArrayList<>();
	     	try {	
	     		Reaction[] transformations = SynthonCreator.create(rxn);
	     		synthonTransformations = Arrays.asList(transformations);
	     	}
	     	catch(Exception e) {
	     		e.printStackTrace();
	     	}
	     	allSynthonTransformations.put(rxn.getName(), synthonTransformations);
		}
	}
	
	private static void processReaction(Reaction rxn,  ConcurrentMap<String,String> processedToOrigBB, 
			ConcurrentMap<String,List<Map<String,String>>> reactionsWithSynthons, ConcurrentMap<String,long[]> fps,
			Map<String,List<Reaction>> synthonTransformations) {

		List<List<String>> reactants = new ArrayList<>();
     	//System.out.println(rxn);
     	getReactants(rxn, processedToOrigBB, reactants,fps);

     	reactionsWithSynthons.putIfAbsent(rxn.getName(), new ArrayList<>());
		
     	//System.out.println("bbs");

     	for(int i=0;i<reactants.size();i++) {
     		List<String> rList = reactants.get(i);
     		ConcurrentMap<String,String> synthonList = new ConcurrentHashMap<String,String>();
     		reactionsWithSynthons.get(rxn.getName()).add(synthonList);
     		Reaction synthonTransformation = synthonTransformations.get(rxn.getName()).get(i);

     		//for(String idcode : rList) {

     		rList.parallelStream().forEach(idcode -> {
     			IDCodeParser parser = new IDCodeParser();
     			StereoMolecule bb = new StereoMolecule();
     			parser.parse(bb, idcode);
     			bb.ensureHelperArrays(Molecule.cHelperCIP);
     			String synthonIDCode = transformToSynthon(synthonTransformation,bb);
     			String origIDCode = processedToOrigBB.get(idcode);
     			if(synthonIDCode!=null) {
     				synthonList.put(synthonIDCode, origIDCode);
     			}
     		});
     	}
   	}

	private static void getReactants(Reaction rxn,  ConcurrentMap<String,String> processedToOrigIDCode,
			List<List<String>> reactants, Map<String,long[]> fps) {
		
		SSSearcherWithIndex[] searchers = new SSSearcherWithIndex[rxn.getReactants()];
		for(int i=0;i<rxn.getReactants();i++) {
			DescriptorHandlerLongFFP512 dhf = new DescriptorHandlerLongFFP512();
			StereoMolecule r = rxn.getReactant(i);
			long[] reactantFFP = dhf.createDescriptor(r);
			SSSearcherWithIndex searcher = new SSSearcherWithIndex();
			searcher.setFragment(r, reactantFFP);
			searchers[i] = searcher;
		}
		for(int i=0;i<rxn.getReactants();i++) {
			List<String> reactantList = Collections.synchronizedList(new ArrayList<String>());
			final int ii = i;
			
			processedToOrigIDCode.keySet().stream().forEach(processedIDcode -> { 
				StereoMolecule mol = new IDCodeParser().getCompactMolecule(processedIDcode);
				if (mol != null) {
					long[] fp = fps.get(processedIDcode);
					if(fp==null)
						return;
					boolean matchesReaction = matchesReactionRole(rxn, searchers,ii,mol,fp);
					if(matchesReaction)
						reactantList.add(processedIDcode);
						
				
				}	
			});
			reactants.add(reactantList);
		}
	}
	
	private static String transformToSynthon(Reaction synthonTransformation, StereoMolecule bb) {
		String synthonIDCode = null;
		Reactor reactor = new Reactor(synthonTransformation, Reactor.MODE_RETAIN_COORDINATES
					+Reactor.MODE_FULLY_MAP_REACTIONS+Reactor.MODE_REMOVE_DUPLICATE_PRODUCTS+Reactor.MODE_ALLOW_CHARGE_CORRECTIONS, Integer.MAX_VALUE);
		bb.ensureHelperArrays(Molecule.cHelperCIP);
		reactor.setReactant(0, bb);

		String[][] productCodes = reactor.getProductIDCodes();
		int productCount = getNormalizedProductsCount(reactor.getProducts());
		
		if(productCount==0 || productCount>1)
			return synthonIDCode;
		else {
     		synthonIDCode = productCodes[0][0]; 
			}
		return synthonIDCode;
	}
	
	private static StereoMolecule getProduct(Reaction rxn, List<StereoMolecule> reactants) {
		Reactor reactor = new Reactor(rxn, Reactor.MODE_RETAIN_COORDINATES
					+Reactor.MODE_FULLY_MAP_REACTIONS+Reactor.MODE_REMOVE_DUPLICATE_PRODUCTS+Reactor.MODE_ALLOW_CHARGE_CORRECTIONS, Integer.MAX_VALUE);
		for(int i=0;i<reactants.size();i++)
			reactor.setReactant(i, reactants.get(i));
		StereoMolecule[][] products = reactor.getProducts();
		return products.length == 1 && products[0].length == 1 ? products[0][0] : null;
	}
	
	private static int getNormalizedProductsCount(StereoMolecule[][] products) {
		Set<Long> normalizedProducts = new HashSet<>();
		for(int i=0;i<products.length;i++) {
			for(int j=0;j<products[i].length;j++) {
				normalizedProducts.add(CanonizerUtil.getTautomerHash(products[i][j],true));
			}
		}
		return normalizedProducts.size();
	}
	
	private void generateCombinatoriaLibraries(ConcurrentMap<String,List<Map<String,String>>> reactionsWithSynthons,
			Map<String,List<Reaction>> synthonTransformations) {

		// iterate over reactions twice (inner loop, outer loop)
		// check if synthons from the first reaction match the generic substructure of the second
		Map<String,List<String>> productsWithSynthons = new HashMap<String,List<String>>();
		Map<String,List<String>> productsWithBBs = new HashMap<String,List<String>>();
		Map<String,List<String>> productsWithReactions = new HashMap<String,List<String>>();
		Set<Reaction> functionalizations = new HashSet<Reaction>();
		for(Reaction rxn : reactions) {
			if(rxn.getReactants()==1)
				functionalizations.add(rxn);
		}
		List<Long> sizes = new ArrayList<>();
		reactions.stream().forEach(reaction -> {
			if(reaction.getReactants()<2)
				return;
			IDCodeParser parser = new IDCodeParser();
			CombinatorialLibrary combiLibrary = new CombinatorialLibrary();
			combiLibrary.reaction = reaction;
			String libraryReaction = reaction.getName();
			combiLibrary.bbSynthons = reactionsWithSynthons.get(libraryReaction);
			List<Map<String,List<Map<String,String>>>> precursorLibs = new ArrayList<>();
			combiLibrary.precursorLibs = precursorLibs;
			if(reaction.getReactants()>2) { //no addditional steps for 3-cmpd reactions
				combiLibrary.cleanup();
				sizes.add(combiLibrary.getSize());
				try {
					writeCombinatorialLibrary(combiLibrary);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				return;
			}
			SSSearcherWithIndex[] searchers = new SSSearcherWithIndex[reaction.getReactants()];
			for(int i=0;i<searchers.length;i++) {
				SSSearcherWithIndex searcher = new SSSearcherWithIndex();
				searcher.setFragment(reaction.getReactant(i), (long[]) null);
				searchers[i] = searcher;
			}

			for(int libReactantID=0;libReactantID<reactionsWithSynthons.get(libraryReaction).size();libReactantID++) {
				Map<String,List<Map<String,String>>> precursorLib = new HashMap<String,List<Map<String,String>>>();
				precursorLibs.add(precursorLib);
				Reaction synthonTransformation = synthonTransformations.get(libraryReaction).get(libReactantID);
				StereoMolecule genericProduct = synthonTransformation.getProduct(0);
				int offset = 0; //offset is required to correctly mutate connector atoms
				for(int a=0;a<genericProduct.getAtoms();a++) {
					if(genericProduct.getAtomicNo(a)>=92)
						offset++;
				}
				
				for(Reaction reactionPrec : reactions) { // precursor reaction
					if(reactionPrec.getReactants()!=2) 
						continue;
					String precursorReaction = reactionPrec.getName();

					for(int reactantID=0;reactantID<reactionsWithSynthons.get(precursorReaction).size();reactantID++) {
						List<Map<String,String>> libSynthons = new ArrayList<>();
						String precursorName = null;
						precursorName = precursorReaction + "_" + reactantID;
						precursorLib.put(precursorName, libSynthons);
						Map<Integer,StereoMolecule> dummyReactants = new HashMap<Integer,StereoMolecule>();
						for(int l=0;l<reactionPrec.getReactants();l++) {
							if(l==reactantID)
								continue;
							else {
								dummyReactants.put(l,reactionPrec.getReactant(l));
							}
						}
						for(int i=0;i<reactionsWithSynthons.get(precursorReaction).size();i++) {
							Map<String,String> precursorSynthons = new HashMap<String,String>();
							libSynthons.add(precursorSynthons);
						}		
						Map<String,String> twoConnectorSynthons = libSynthons.get(reactantID); 
					
						Map<String,String> synthons = reactionsWithSynthons.get(precursorReaction).get(reactantID);
						for(String synthon : synthons.keySet()) {
							StereoMolecule mol = new StereoMolecule();
							String bb = synthons.get(synthon);
							parser.parse(mol, synthon);
//if (reactionPrec.getReactants() <= reactantID)  // to prevent OutOfBoundsExceptions; TLS 25-Mar-2024 in dummyReaction
//	continue;
							StereoMolecule reactedBB = dummyReaction(bb,reactionPrec,dummyReactants,reactantID); 
							if(reactedBB==null)
								continue;
							// we create candidate products by a "dummy reaction" and see if it still matches the substructure query of rxn1
							if (matchesReactionRole(reaction, searchers, libReactantID, reactedBB,null)) { // synthon also matches library rxn
								mutateConnectorAtoms(mol,offset);
								String transformedSynthon = transformToSynthon(synthonTransformation,mol);
								if(transformedSynthon!=null)
									twoConnectorSynthons.put(transformedSynthon,bb);
							}
							//now check if a functionalization/deprotection step can lead to additional matches to the library rxn
							else {
								for(Reaction functionalization : functionalizations) {
									if(functionalizations.contains(reactionPrec))
										continue; //don't couple two functionalization reactions
									SSSearcher ssearcher = new SSSearcher();
									ssearcher.setFragment(functionalization.getReactant(0));
									ssearcher.setMolecule(reactedBB);
									if (ssearcher.isFragmentInMolecule()) {// can be functionalized
										StereoMolecule product = getProduct(functionalization,Arrays.asList(reactedBB));
										if(product==null)
											continue;
										else {
											if(matchesReactionRole(reaction, searchers,libReactantID, product,null)) { //functionalized BB matches rxn1
												//now functionalize synthon BB and create second linker
												StereoMolecule prod = getProduct(functionalization,Arrays.asList(mol));
												if(prod==null)
													continue;
												mutateConnectorAtoms(prod,offset);
						
												String transformedSynthon = transformToSynthon(synthonTransformation,prod);
												if(!precursorLib.containsKey(precursorName)) {
													List<Map<String,String>> libSynthons2 = new ArrayList<>();
													precursorLib.put(precursorName, libSynthons2);
													for(int i=0;i<reactionsWithSynthons.get(precursorReaction).size();i++) {
														Map<String,String> precursorSynthons2 = new HashMap<String,String>();
														libSynthons2.add(precursorSynthons2);
													}
												}
												precursorLib.get(precursorName).get(reactantID).put(transformedSynthon, bb);
											}
										}
									}
								}
							}
						}
						
						// if the BB that is used in the precursor step also matches any role of the second reaction, it is excluded (otherwise multiple products may be formed)
						for(int j=0;j<reactionsWithSynthons.get(precursorReaction).size();j++) { //get synthons from rxn1 that are compatible with the two-connector synthon
							if(j==reactantID) {
								continue;
							}
							else {
								Map<String,String> compatibleSynthons = reactionsWithSynthons.get(precursorReaction).get(j);
								Map<String,String> mutatedSynthons = new HashMap<String,String>();
								for(String s: compatibleSynthons.keySet()) {
									StereoMolecule mutatedSynthon = new StereoMolecule();
									boolean compatible = true;
									parser.parse(mutatedSynthon, s);
									for(SSSearcherWithIndex sssearcher : searchers) {
										sssearcher.setMolecule(mutatedSynthon,(long[]) null);
										if(sssearcher.isFragmentInMolecule())
											compatible = false;
									}
									if(!compatible) //matches also role from second reaction --> don't add to library
										continue;
	
									mutateConnectorAtoms(mutatedSynthon,offset);
									
									mutatedSynthons.put(mutatedSynthon.getIDCode(), compatibleSynthons.get(s));
								}
								libSynthons.get(j).putAll(mutatedSynthons);
							}
						}
					}
				}
			}
			combiLibrary.cleanup();
			sizes.add(combiLibrary.getSize());
			System.out.println(reaction.getName());
			combiLibrary.generateRandomProducts(1000,productsWithSynthons,productsWithBBs,productsWithReactions);
			try {
				writeCombinatorialLibrary(combiLibrary);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		});
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("space.dwar"), StandardCharsets.UTF_8));
			DWARFileCreator creator = new DWARFileCreator(writer);
			List<Integer> synthonColumns = new ArrayList<>();
			for(int i=0;i<4;i++) {
				synthonColumns.add(creator.addStructureColumn("Synthon" + (i + 1), "IDcode"));
			} 
			List<Integer> bbColumns = new ArrayList<>();
			for(int i=0;i<4;i++) {
				bbColumns.add(creator.addStructureColumn("BB" + (i + 1), "IDcode"));
			}
			List<Integer> reactionColumns = new ArrayList<>();
			for(int i=0;i<3;i++) {
				reactionColumns.add(creator.addAlphanumericalColumn("Reaction" + (i + 1)));
			}
			
			int structureColumn = creator.addStructureColumn("Product", "IDcode");
			creator.writeHeader(-1);
			productsWithSynthons.entrySet().stream().forEach(e -> {
				creator.setRowStructure(e.getKey(), structureColumn);
				List<String> synthons = e.getValue();
				for(int i=0;i<synthons.size();i++) {
					creator.setRowStructure(synthons.get(i), synthonColumns.get(i));			
				}
				List<String> bbs = productsWithBBs.get(e.getKey());
				for(int i=0;i<bbs.size();i++) {
					creator.setRowStructure(bbs.get(i), bbColumns.get(i));			
				}
				List<String> reactions = productsWithReactions.get(e.getKey());
				for(int i=0;i<reactions.size();i++) {
					creator.setRowValue(reactions.get(i), reactionColumns.get(i));			
				}
				try {
					creator.writeCurrentRow();
				} catch (IOException e3) {
					// TODO Auto-generated catch block
					e3.printStackTrace();
				}
			}
			);
			creator.writeEnd();
		}
		catch(Exception e) {
			System.out.println("########");
			e.printStackTrace();
		}

		long size = sizes.stream().reduce(0L,(e1,e2) -> e1+e2);
		System.out.println("totSize");
		System.out.println(size);
	}
	
	private static StereoMolecule dummyReaction(String bbIDCode, Reaction rxn, Map<Integer,StereoMolecule> reactants, int reactantID){
		StereoMolecule product = null;
		SSSearcher searcher = new SSSearcher();
		searcher.setFragment(rxn.getReactant(reactantID));
		StereoMolecule mol = new IDCodeParser().getCompactMolecule(bbIDCode);
		if (mol != null) {
			searcher.setMolecule(mol);
			if (searcher.isFragmentInMolecule()) {
				Reactor reactor = new Reactor(rxn, Reactor.MODE_RETAIN_COORDINATES
     					+Reactor.MODE_FULLY_MAP_REACTIONS+Reactor.MODE_REMOVE_DUPLICATE_PRODUCTS+Reactor.MODE_ALLOW_CHARGE_CORRECTIONS, Integer.MAX_VALUE);
				reactor.setReactant(reactantID, mol);
				for(int i : reactants.keySet()) {
					reactor.setReactant(i, reactants.get(i));
				}
				StereoMolecule[][] products = reactor.getProducts();
				if(products.length == 1 && products[0].length == 1)
					product = products[0][0];
			}
		}
		return product;
	}
	
	private static void mutateConnectorAtoms(StereoMolecule mol, int offset) {
		for(int a=0;a<mol.getAtoms();a++) {
			int atomicNo = mol.getAtomicNo(a);
			if(atomicNo>=92) {
				atomicNo += offset;
				mol.setAtomicNo(a, atomicNo);
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}

	private  void writeCombinatorialLibrary(CombinatorialLibrary combiLib) throws IOException {
		File htmcdir = new File(this.outdirectory + "/CombinatorialLibraries/");
		htmcdir.mkdir();
		File htmcSubDir = new File(htmcdir.getAbsolutePath() + "/" +  combiLib.reaction.getName());
		htmcSubDir.mkdir();
		File dirA = new File(htmcSubDir.getAbsolutePath() + "/A");
		dirA.mkdir();
		File dirB = new File(htmcSubDir.getAbsolutePath() + "/B");
		dirB.mkdir();
		File dirC = null;
		File dirD = null;
		if(combiLib.bbSynthons.size()>2) {
			dirC = new File(htmcSubDir.getAbsolutePath() + "/C");
			dirC.mkdir();
		}
		if(combiLib.bbSynthons.size()>3) {
			dirD = new File(htmcSubDir.getAbsolutePath() + "/D");
			dirD.mkdir();
		}

		if(combiLib.precursorLibs.size()==2) {
			for(int i=0;i<combiLib.precursorLibs.size();i++) {
				File dir = null;
				if(i==0)
					dir = new File(dirA.getAbsolutePath() + "/virtual_bbs");
				else 
					dir = new File(dirB.getAbsolutePath() + "/virtual_bbs");
				dir.mkdir();
				
				for(String flowReaction : combiLib.precursorLibs.get(i).keySet()) {
					String flowDir = dir.getAbsoluteFile() + "/" + flowReaction;
					List<Map<String,String>> flowSynthons = combiLib.precursorLibs.get(i).get(flowReaction);
					int counter=1;
					for(Map<String,String> synthons : flowSynthons) {
						if(synthons.isEmpty())
							continue;
						BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
								new FileOutputStream(flowDir + "_" + counter + ".dwar"), StandardCharsets.UTF_8));
						DWARFileCreator creator = new DWARFileCreator(writer);
						int synthonColumn = creator.addStructureColumn("Synthon", "IDcode");
						int structureColumn = creator.addStructureColumn("Building Block", "IDcode");
						List<String> fields = new ArrayList<>();
						if (bbData != null) {
							bbData.values().stream().forEach(e -> {
								for (String key : e.keySet()) {
									if (!fields.contains(key))
										fields.add(key);
								}
							});
						}
						List<Integer> fieldIndexes = new ArrayList<>();
						for(String field : fields) {
							int columnIndex = creator.addAlphanumericalColumn(field);
							fieldIndexes.add(columnIndex);
						}

						creator.writeHeader(-1);
						for(String s : synthons.keySet()) {
							String origIDCode = synthons.get(s);
							if (bbData != null) {
								Map<String,List<String>> propertyMap = bbData.get(origIDCode);
								for(int j=0;j<fields.size();j++) {
									String field = fields.get(j);
									StringBuilder sb = new StringBuilder();
									List<String> values = propertyMap.get(field);
									if(values!=null) {
										values.forEach(e -> sb.append(e).append(";"));
									}
									creator.setRowValue(sb.toString(), fieldIndexes.get(j));
								}
							}
							creator.setRowStructure(s, synthonColumn);
							creator.setRowStructure(origIDCode,structureColumn);
							creator.writeCurrentRow();
						}
						creator.writeEnd();
						counter++;
					}
				}
			}
		}

		int i=0;
		for(Map<String,String> synthons : combiLib.bbSynthons) {
			if(synthons.isEmpty())
				continue;

			File dir = null;
			switch(i) {
				case 0:
					dir = dirA;
					break;
				case 1: 
					dir = dirB;
					break;
				case 2: 
					dir = dirC;
					break;
				case 3: 
					dir = dirD;
					break;
				default: 
					break;
			}

			String outfile = dir + "/" + combiLib.reaction.getName() + ".dwar";
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile), StandardCharsets.UTF_8));
			DWARFileCreator creator = new DWARFileCreator(writer);
			int synthonColumn = creator.addStructureColumn("Synthon", "IDcode");
			int structureColumn = creator.addStructureColumn("Building Block", "IDcode");
			List<String> fields = new ArrayList<>();
			if (bbData != null) {
				bbData.values().stream().forEach(e -> {
					for (String key : e.keySet()) {
						if (!fields.contains(key))
							fields.add(key);
					}
				});
			}
			List<Integer> fieldIndeces = new ArrayList<>();
			for(String field : fields) {
				int columnIndex = creator.addAlphanumericalColumn(field);
				fieldIndeces.add(columnIndex);
			}
			creator.writeHeader(-1);
			for(String s : synthons.keySet()) {
				String origIDCode = synthons.get(s);
				if (bbData != null) {
					Map<String,List<String>> propertyMap = bbData.get(origIDCode);
					for(int j=0;j<fields.size();j++) {
						String field = fields.get(j);
						StringBuilder sb = new StringBuilder("");
						List<String> values = propertyMap.get(field);
						if(values!=null) {
							values.forEach(e -> sb.append(e + ";"));
						}
						creator.setRowValue(sb.toString(), fieldIndeces.get(j));
					}
				}
				creator.setRowStructure(s, synthonColumn);
				creator.setRowStructure(synthons.get(s),structureColumn);
				creator.writeCurrentRow();
			}
			creator.writeEnd();
			i++;
		}
	}
	
	private static class CombinatorialLibrary {
		public Reaction reaction;
		public List<Map<String,String>> bbSynthons;
		public List<Map<String,List<Map<String,String>>>> precursorLibs;
		
		List<String> toRemove = new ArrayList<>();
		
		private void cleanup() { //remove precursor reactions with missing synthon sets
			for(Map<String,List<Map<String,String>>> precursorLib : precursorLibs) {
				for(String reaction : precursorLib.keySet()) {
					for(Map<String,String> synthons : precursorLib.get(reaction)) {
						if(synthons==null)
							toRemove.add(reaction);
						else if(synthons.isEmpty())
							toRemove.add(reaction);
					}
				}
				for(String r : toRemove) {
					precursorLib.remove(r);
				}
				toRemove = new ArrayList<>();
			}
		}
		
		public long getSize() {
			long sizeOneStep = 1;
			// first calculate size of the single-step space
			for(Map<String,String> synthons : bbSynthons) {
				sizeOneStep*=synthons.size();
			}
			if(precursorLibs.size()==2) {
				long reactantsA = bbSynthons.get(0).size();
				Map<String,List<Map<String,String>>> precursorLibsA = precursorLibs.get(0);
				long reactantsB = bbSynthons.get(1).size();
				Map<String,List<Map<String,String>>> precursorLibsB = precursorLibs.get(1);
				for(String precReaction : precursorLibsA.keySet()) {
					List<Map<String,String>> precLibA = precursorLibsA.get(precReaction);
					int precSize = 1;
					for(Map<String,String> precLib : precLibA) {
						precSize*=precLib.size();
					}
					reactantsA += precSize;
				}
				for(String precReaction : precursorLibsB.keySet()) {
					List<Map<String,String>> precLibB = precursorLibsB.get(precReaction);
					int precSize = 1;
					for(Map<String,String> precLib : precLibB) {
						precSize*=precLib.size();
					}
					reactantsB += precSize;
				}
				return reactantsA*reactantsB;
			}
			else 
				return sizeOneStep;
		}

		
		public void generateRandomProducts(int nProducts, Map<String,List<String>> productsWithSynthons, Map<String,List<String>> productsWithBBs, Map<String,
				List<String>> productsWithReactions ) {
			Random rnd = new Random();
			IDCodeParser parser = new IDCodeParser();
			long max = productsWithSynthons.keySet().size() +  Math.min(getSize(), nProducts);
			while(productsWithSynthons.keySet().size()<max) {
				List<StereoMolecule> synthons = new ArrayList<>();
				List<String> bbs = new ArrayList<>();
				List<String> reactions = new ArrayList<>();
				List<StereoMolecule> preCoupledSynthons = new ArrayList<>();
				String precursorName = "";
				reactions.add(reaction.getName());
				for(int i=0;i<reaction.getReactants();i++) {
					int r = rnd.nextInt(2);
					Map<String,List<Map<String,String>>> libs = precursorLibs.get(i);
					if(libs.isEmpty() || r==0) { //pick commercial bb
						if(bbSynthons.get(i).isEmpty())
							return ;
						int r2 = rnd.nextInt(bbSynthons.get(i).size());
						StereoMolecule mol = new StereoMolecule();
						Map<String,String> s = bbSynthons.get(i);
						List<String> keys = new ArrayList<>(s.keySet());
						parser.parse(mol, keys.get(r2));
						synthons.add(mol);
						bbs.add(s.get(keys.get(r2)));
						preCoupledSynthons.add(mol);
					}
					else { //create random virtual bb
						//pick random precursor reaction
						List<StereoMolecule> precursors = new ArrayList<>();
						int reactionNr = rnd.nextInt(libs.keySet().size());
						List<String> keys = new ArrayList<>(libs.keySet());
						List<Map<String,String>> lib = libs.get(keys.get(reactionNr));
						precursorName = keys.get(reactionNr);
						reactions.add(precursorName);
						for(Map<String,String> precursorBBs : lib) {
							List<String> keySet2 = new ArrayList<>(precursorBBs.keySet());
							int r3 = rnd.nextInt(keySet2.size());
							StereoMolecule mol = new StereoMolecule();
							parser.parse(mol, keySet2.get(r3));
							synthons.add(mol);
							bbs.add(precursorBBs.get(keySet2.get(r3)));
							precursors.add(mol);
						}
						try {
							preCoupledSynthons.add(SynthonReactor.react(precursors));
						}
						catch(Exception e) {
							e.printStackTrace();
						}
					}
				}
				
				try {
					StereoMolecule product = SynthonReactor.react(preCoupledSynthons);
					if(product == null || productsWithSynthons.containsKey(product.getIDCode()))
						continue;
					List<String> bbIDCodes = synthons.stream().map(e -> e.getIDCode()).collect(Collectors.toList());
					productsWithSynthons.put(product.getIDCode(),bbIDCodes);
					productsWithBBs.put(product.getIDCode(),bbs);
					productsWithReactions.put(product.getIDCode(),reactions);
				}
				catch(Exception e) {
					e.printStackTrace();
				}
			}
		}
	}
	

	public static boolean matchesReactionRole(Reaction rxn, SSSearcherWithIndex[] searchers, int component, StereoMolecule reactant,
			long[] index) {
		boolean isMatch = true;
		SSSearcherWithIndex searcher = searchers[component];
		searcher.setMolecule(reactant,index);
		if (searcher.isFragmentInMolecule()) {
			//check if reactant also matches other roles in the reaction, if yes, exclude it (to prevent self-polymerization)
			for(int j=0;j<rxn.getReactants();j++) {
				if (component != j
				 && !rxn.getReactant(component).getIDCode().equals(rxn.getReactant(j).getIDCode())) {
					SSSearcherWithIndex searcher2 = searchers[j];
					searcher2.setMolecule(reactant,index);
					if(searcher2.isFragmentInMolecule()) {
						isMatch=false;
					}
				}
			}
		}
		else {
			isMatch = false;
		}
		return isMatch;
	}
}
