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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.util.IntArrayComparator;
import com.actelion.research.util.SortedList;

import java.io.*;
import java.net.URI;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.AbstractMap.SimpleEntry;
import java.util.zip.GZIPInputStream;

/**
 * PDBFileParser
 * Created by korffmo1 on 20.03.18.
 */
public class PDBFileParser {

    // 05-FEB-18
    public static final String DATE_FORMAT = "dd-MMM-yy";

    public static final String TAG_HEADER = "HEADER";

    public static final String TAG_OBSOLTE = "OBSLTE";

    public static final String TAG_TITLE = "TITLE";

    public static final String TAG_SPLIT = "SPLIT";

    public static final String TAG_CAVEAT = "CAVEAT";

    // Mandatory
    public static final String TAG_COMPND = "COMPND";
    // Mandatory
    public static final String TAG_SOURCE = "SOURCE";
    // Mandatory
    public static final String TAG_KEYWDS = "KEYWDS";

    // Mandatory
    public static final String TAG_EXPDTA = "EXPDTA";

    // Optional Mandatory for NMR ensemble entries.
    public static final String TAG_NUMMDL = "NUMMDL";

    // Optional Mandatory for NMR minimized average Structures or when the entire polymer chain contains C alpha or P atoms only.
    public static final String TAG_MDLTYP = "MDLTYP";

    // Mandatory
    public static final String TAG_AUTHOR = "AUTHOR";
    // Mandatory
    public static final String TAG_REVDAT = "REVDAT";
    // Optional Mandatory for a replacement entry.
    public static final String TAG_SPRSDE = "SPRSDE";
    // Optional Mandatory for a publication describes the experiment.
    public static final String TAG_JRNL = "JRNL";

    // Remark has format 'other'
    public static final String TAG_REMARK = "REMARK";

    // Optional Mandatory for a re-refined structure
    public static final String TAG_REMARK0 = "REMARK 0";

    // Optional
    public static final String TAG_REMARK1 = "REMARK 1";
    // Mandatory
    public static final String TAG_REMARK2 = "REMARK 2";
    // Mandatory
    public static final String TAG_REMARK3 = "REMARK 3";
    // Optional
    public static final String TAG_REMARK_N = "REMARK N";

    // Optional Mandatory for all polymers.
    public static final String TAG_DBREF = "DBREF";
    // Optional Mandatory when certain sequence database accession and/or sequence numbering does not fit preceding DBREF format.
    public static final String TAG_DBREF1_DBREF2 = "DBREF1/DBREF2";
    // Optional Mandatory if sequence conflict exists.
    public static final String TAG_SEQADV = "SEQADV";
    // Mandatory if ATOM records exist.
    public static final String TAG_SEQRES = "SEQRES";
    // Optional Mandatory if modified group exists in the coordinates.
    public static final String TAG_MODRES = "MODRES";
    // Optional Mandatory if a non-standard group other than water appears in the coordinates.
    public static final String TAG_HET = "HET";
    // Optional Mandatory if
    public static final String TAG_HETNAM = "HETNAM";
    // Optional
    public static final String TAG_HETSYN = "HETSYN";
    // Optional
    public static final String TAG_FORMUL = "FORMUL";
    // Optional
    public static final String TAG_HELIX = "HELIX";

    // Optional
    public static final String TAG_SHEET = "SHEET";
    // Optional Mandatory if a disulfide bond is present.
    public static final String TAG_SSBOND = "SSBOND";
    // Optional Mandatory if non-standard residues appear in a polymer
    public static final String TAG_LINK = "LINK";

    // Optional
    public static final String TAG_CISPEP = "CISPEP";
    // Optional
    public static final String TAG_SITE = "SITE";
    // Mandatory
    public static final String TAG_CRYST1 = "CRYST1";
    // Mandatory
    public static final String TAG_ORIGX1 = "ORIGX1";
    // Mandatory
    public static final String TAG_ORIGX2 = "ORIGX2";
    // Mandatory
    public static final String TAG_ORIGX3 = "ORIGX3";
    // Mandatory
    public static final String TAG_SCALE1 = "SCALE1";
    // Mandatory
    public static final String TAG_SCALE2 = "SCALE2";
    // Mandatory
    public static final String TAG_SCALE3 = "SCALE3";
    //
    public static final String TAG_MTRIX1 = "MTRIX1";

    public static final String TAG_MTRIX2 = "MTRIX2";
    // Optional Mandatory if the complete asymmetric unit must be generated from the given coordinates using non-crystallographic symmetry.
    public static final String TAG_MTRIX3 = "MTRIX3";

    // Optional Mandatory if more than one model is present in the entry.
    public static final String TAG_MODEL = "MODEL";

    // Optional Mandatory if standard residues exist.
    public static final String TAG_ATOM = "ATOM";
    // Optional
    public static final String TAG_ANISOU = "ANISOU";
    // Optional Mandatory if ATOM records exist.
    public static final String TAG_TER = "TER";
    // Optional Mandatory if non-standard group exists.
    public static final String TAG_HETATM = "HETATM";
    // Optional Mandatory if MODEL appears.
    public static final String TAG_ENDMDL = "ENDMDL";
    // Optional Mandatory if non-standard group appears and if LINK or SSBOND records exist.
    public static final String TAG_CONECT = "CONECT";

    // Mandatory
    public static final String TAG_MASTER = "MASTER";
    // Mandatory
    public static final String TAG_END = "END";



    private final DateFormat dfDateDeposition;

	private final HetNameParser hetNameParser;

    private final HetSynonymParser hetSynonymParser;
    private final FormulaParser formulaParser;

    private final SiteParser siteParser;

    private final ModelParser modelParser;

	public PDBCoordEntryFile getFromPDB(String pdbID) throws Exception {
		URLConnection con = new URI("https://files.rcsb.org/download/"+pdbID+".pdb.gz").toURL().openConnection();
		return new PDBFileParser().parse(new BufferedReader(new InputStreamReader(new GZIPInputStream(con.getInputStream()))));
	}

    public PDBFileParser() {
        dfDateDeposition = new SimpleDateFormat(DATE_FORMAT);
	    RemarkParser remarkParser = new RemarkParser();
        hetNameParser = new HetNameParser();
        hetSynonymParser = new HetSynonymParser();
        formulaParser = new FormulaParser();
        siteParser = new SiteParser();
        modelParser = new ModelParser();
    }

    public PDBCoordEntryFile parse(File fiPDB) throws IOException, ParseException {
		InputStream stream = fiPDB.getName().toLowerCase().endsWith(".pdb.gz") ?
				new GZIPInputStream(new FileInputStream(fiPDB))
			  : new FileInputStream(fiPDB);
		return parse(new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8)));
    }

    public PDBCoordEntryFile parse(BufferedReader br) throws IOException, ParseException {
        PDBCoordEntryFile pdbCoordEntryFile = new PDBCoordEntryFile();
        
		ArrayList<String> liRaw = new ArrayList<>();

		String sCurrentLine;
		while ((sCurrentLine = br.readLine()) != null) {
			liRaw.add(sCurrentLine);
		}

        int indexLine = 0;
        
        String lHeader = liRaw.get(indexLine);

        if(lHeader.startsWith(TAG_HEADER)) {
        	try {
        		indexLine = parseHeader(lHeader, pdbCoordEntryFile);}
        	catch(Exception e) {
        		indexLine++;
        	}
    	}

	    while(indexLine<liRaw.size()
		   && !liRaw.get(indexLine).startsWith(TAG_ATOM)
		   && !liRaw.get(indexLine).startsWith(TAG_HETATM)) {
	        // Not mandatory
	    	if(liRaw.get(indexLine).startsWith(TAG_OBSOLTE)) {
	            SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_OBSOLTE);
	            pdbCoordEntryFile.setObsolete(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_TITLE)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_TITLE);
	            pdbCoordEntryFile.setTitle(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }

	        // Not mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_SPLIT)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_SPLIT);
	            pdbCoordEntryFile.setSplit(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        // Not mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_CAVEAT)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_CAVEAT);
	            pdbCoordEntryFile.setSplit(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        // Mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_COMPND)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_COMPND);
	            pdbCoordEntryFile.setCompound(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	        //} else {
	        //    throw new RuntimeException("Missing " + TAG_COMPND);
	        //}
	
	        // Mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_SOURCE)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_SOURCE);
	            pdbCoordEntryFile.setSource(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	        //} else {
	        //   throw new RuntimeException("Missing " + TAG_SOURCE);
	        
	
	        // Mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_KEYWDS)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_KEYWDS);
	            pdbCoordEntryFile.setKeywords(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        } //else {
	            //throw new RuntimeException("Missing " + TAG_KEYWDS);
	        //}
	
	        // Mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_EXPDTA)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_EXPDTA);
	            pdbCoordEntryFile.setExpdata(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        } //else {
	        //    throw new RuntimeException("Missing " + TAG_EXPDTA);
	        //}
	
	        // Optional Mandatory for NMR ensemble entries.
	        if(liRaw.get(indexLine).startsWith(TAG_NUMMDL)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_NUMMDL);
	            pdbCoordEntryFile.setNummdl(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        // Optional Mandatory for NMR minimized average Structures or when the entire polymer chain contains C alpha or P atoms only.
	        if(liRaw.get(indexLine).startsWith(TAG_MDLTYP)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_MDLTYP);
	            pdbCoordEntryFile.setMdltyp(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        // Mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_AUTHOR)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_AUTHOR);
	            pdbCoordEntryFile.setAuthor(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        // Mandatory
	        if(liRaw.get(indexLine).startsWith(TAG_REVDAT)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_REVDAT);
	            pdbCoordEntryFile.setRevdat(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        // Optional Mandatory for a replacement entry.
	        if(liRaw.get(indexLine).startsWith(TAG_SPRSDE)) {
	        	SimpleEntry<String, Integer> siIndex = parseOneTimeMultipleLines(liRaw, indexLine, TAG_SPRSDE);
	            pdbCoordEntryFile.setSprsde(siIndex.getKey());
	            indexLine = siIndex.getValue();
	        }
	
	        // Optional Mandatory for a publication describes the experiment.
	        if(liRaw.get(indexLine).startsWith(TAG_JRNL)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_JRNL);
	            pdbCoordEntryFile.setJrnl(liIndex);
	            indexLine += liIndex.size();
	        }

			while (liRaw.get(indexLine).startsWith(TAG_REMARK)) {
				String line = liRaw.get(indexLine).trim();
				if (line.contains("RESOLUTION.") && line.endsWith("ANGSTROMS."))
					pdbCoordEntryFile.setResolution(line.substring(12+line.indexOf("RESOLUTION."), line.length()-11).trim());
				indexLine++;
			}
//	        remarkParser.parse(liRaw, indexLine);
//	        HashMap<Integer, String> hmNo_Remark = remarkParser.getHmNo_Remark();
//	        indexLine = remarkParser.getIndexLine();
//	        pdbCoordEntryFile.setRemarks(hmNo_Remark);

	        // Mandatory for all polymers.
	        if(liRaw.get(indexLine).startsWith(TAG_DBREF)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_DBREF);
	            pdbCoordEntryFile.setDBRef(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        // Mandatory for all polymers.
	        if(liRaw.get(indexLine).startsWith(TAG_DBREF1_DBREF2)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_DBREF1_DBREF2);
	            pdbCoordEntryFile.setDBRef1DBRef2(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        // Mandatory if sequence conflict exists.
	        if(liRaw.get(indexLine).startsWith(TAG_SEQADV)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_SEQADV);
	            pdbCoordEntryFile.setSEQADV(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        // Mandatory if ATOM records exist.
	        // Primary sequence of backbone residues.
	        if(liRaw.get(indexLine).startsWith(TAG_SEQRES)) {
	            ListInteger<String> liIndexChains = parseMultipleTimesMultipleLinesSEQRES(liRaw, indexLine, TAG_SEQRES);
	            pdbCoordEntryFile.setSEQRES(liIndexChains.getLi());
	            indexLine = liIndexChains.getId();
	        }
	
	        // Mandatory if modified group exists in the coordinates.
	        // Identification of modifications to standard residues.
	        if(liRaw.get(indexLine).startsWith(TAG_MODRES)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_MODRES);
	            pdbCoordEntryFile.setModRes(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        // Mandatory if a non-standard group other than water appears in the coordinates.
	        // Identification of non-standard groups heterogens).
	
	        String patternHET = TAG_HET + " ";
	
	        if(liRaw.get(indexLine).startsWith(patternHET)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, patternHET);
	            pdbCoordEntryFile.setHet(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        //
	        // Compound name of the heterogens.
	        if(liRaw.get(indexLine).startsWith(TAG_HETNAM)) {
	            hetNameParser.parse(liRaw, indexLine);
	            HashMap<String, String> hmId_Name = hetNameParser.getHMId_Name();
	            pdbCoordEntryFile.setHmId_Name(hmId_Name);
	            indexLine = hetNameParser.getIndexLine();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_HETSYN)) {
	            hetSynonymParser.parse(liRaw, indexLine);
	            HashMap<String, String> hmId_Synonyms = hetSynonymParser.getHMId_Synonyms();
	            pdbCoordEntryFile.setHmId_Synonyms(hmId_Synonyms);
	            indexLine = hetSynonymParser.getIndexLine();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_FORMUL)) {
	            formulaParser.parse(liRaw, indexLine);
	            HashMap<String, String> hmId_Formula = formulaParser.getHMId_Formula();
	            pdbCoordEntryFile.setHmId_Formula(hmId_Formula);
	            indexLine = formulaParser.getIndexLine();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_HELIX)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_HELIX);
	            pdbCoordEntryFile.setHelix(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_SHEET)) {
	            List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_SHEET);
	            pdbCoordEntryFile.setSheet(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_SSBOND)) {
				List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_SSBOND);
	            pdbCoordEntryFile.setSSBond(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_LINK)) {
				List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_LINK);
	            pdbCoordEntryFile.setLink(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_CISPEP)) {
				List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_CISPEP);
	            pdbCoordEntryFile.setCisPep(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_SITE)) {
	            siteParser.parse(liRaw, indexLine);
	            HashMap<String, String> hmId_Site = siteParser.getHMId_Site();
	            pdbCoordEntryFile.setHmId_Site(hmId_Site);
	            indexLine = siteParser.getIndexLine();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_CRYST1)) {
	            pdbCoordEntryFile.setCryst1(liRaw.get(indexLine).substring(6));
	            indexLine++;
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_ORIGX1)) {
	            pdbCoordEntryFile.setOrigX1(liRaw.get(indexLine).substring(10));
	            indexLine++;
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_ORIGX2)) {
	            pdbCoordEntryFile.setOrigX2(liRaw.get(indexLine).substring(10));
	            indexLine++;
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_ORIGX3)) {
	            pdbCoordEntryFile.setOrigX3(liRaw.get(indexLine).substring(10));
	            indexLine++;
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_SCALE1)) {
	            pdbCoordEntryFile.setScale1(liRaw.get(indexLine).substring(10));
	            indexLine++;
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_SCALE2)) {
	            pdbCoordEntryFile.setScale2(liRaw.get(indexLine).substring(10));
	            indexLine++;
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_SCALE3)) {
	            pdbCoordEntryFile.setScale3(liRaw.get(indexLine).substring(10));
	            indexLine++;
	        }

			if(liRaw.get(indexLine).startsWith(TAG_MODEL)) {
				pdbCoordEntryFile.setModel(liRaw.get(indexLine).substring(10));
				indexLine++;
			}

			if(liRaw.get(indexLine).startsWith(TAG_MTRIX1)) {
				List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_MTRIX1);
	            pdbCoordEntryFile.setMtrix1(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_MTRIX2)) {
				List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_MTRIX2);
	            pdbCoordEntryFile.setMtrix2(liIndex);
	            indexLine += liIndex.size();
	        }
	
	        if(liRaw.get(indexLine).startsWith(TAG_MTRIX3)) {
				List<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_MTRIX3);
	            pdbCoordEntryFile.setMtrix3(liIndex);
	            indexLine += liIndex.size();
	        }
	    }
	    
	    //indexLine--;
        TreeSet<AtomRecord> hetAtomRecords = new TreeSet<>();
		TreeSet<AtomRecord> protAtomRecords = new TreeSet<AtomRecord>();
        modelParser.parse(liRaw, indexLine,protAtomRecords,hetAtomRecords);

		ArrayList<AtomRecord> protAtomList = new ArrayList<>();
		for (AtomRecord ar : protAtomRecords)
			protAtomList.add(ar);

		ArrayList<AtomRecord> hetAtomList = new ArrayList<>();
		for (AtomRecord ar : hetAtomRecords)
			hetAtomList.add(ar);

		pdbCoordEntryFile.setProtAtomRecords(protAtomList);
        pdbCoordEntryFile.setHetAtomRecords(hetAtomList);
        //List<ModelModel> liModelModel = modelParser.getLiModelModel();
        //pdbCoordEntryFile.setLiModelModel(liModelModel);

        indexLine = modelParser.getIndexLine();

        //
        // Parsing atom connections
        //
		SortedList<int[]> bonds = new SortedList<>(new IntArrayComparator());
		indexLine = parseCONECTLines(liRaw, indexLine, bonds);
		pdbCoordEntryFile.setLiConnect(bonds);

 /* replaced by something more efficient, because the original was limited to atom indexes <= 9999; TLS 6Nov2024
        if(liRaw.get(indexLine).startsWith(TAG_CONECT)) {
            ListInteger<String> liIndex = parseMultipleTimesOneLine(liRaw, indexLine, TAG_CONECT);
            for(String bondInfo:liIndex.getLi()) {
            	bondInfo = bondInfo.trim();
            	String[] strArr = bondInfo.split("\\s+");
            	try {
            	int firstAtom = Integer.parseInt(strArr[0]);
            	IntStream.range(1,strArr.length).forEach(e -> {
            		int[] bond = new int[2];
            		bond[0] = firstAtom;
            		bond[1] = Integer.parseInt(strArr[e]);
            		bonds.add(bond);
            	});
            	}
            	catch(Exception e) {
            		continue;
            	}

            }
            indexLine = liIndex.getId();
        }
        pdbCoordEntryFile.setLiConnect(bonds);

  */

        if(liRaw.get(indexLine).startsWith(TAG_MASTER)) {
            pdbCoordEntryFile.setMaster(liRaw.get(indexLine).substring(10).trim());
            indexLine++;
        } 
        if(liRaw.get(indexLine).startsWith(TAG_END)) {
            pdbCoordEntryFile.setEnd(true);
        } else {
            pdbCoordEntryFile.setEnd(false);
        }
        
        return pdbCoordEntryFile;

    }

    private int parseHeader(String lHeader, PDBCoordEntryFile pdbCoordEntryFile) throws ParseException {

         int length = lHeader.length();

         pdbCoordEntryFile.setClassification(lHeader.substring (10, Math.min(length,50)).trim() );

         Date date = dfDateDeposition.parse(lHeader.substring (50, Math.min(length,59)).trim() );

         pdbCoordEntryFile.setDateDeposition(date);

         pdbCoordEntryFile.setID(lHeader.substring (62, Math.min(length,66)).trim() );

         return 1;

    }

    /**
     * One time, multiple lines: There are records that conceptually exist only once in an entry, but the
     information content may exceed the number of columns available. These records are therefore
     continued on subsequent lines.
     * @param liRaw
     * @param indexLine
     * @return
     * @throws ParseException
     */

    private AbstractMap.SimpleEntry<String, Integer> parseOneTimeMultipleLines(List<String> liRaw, int indexLine, String tag) throws ParseException {

        String l0 = liRaw.get(indexLine);
        if(!l0.startsWith(tag)) {
            throw new RuntimeException("Error in parsing " + tag);
        }
        String titleSub0 = l0.substring(tag.length()).trim();

        StringBuilder sb = new StringBuilder(titleSub0);

        indexLine++;

        int start = indexLine;

        for (int i = start; i < liRaw.size(); i++) {

            String l = liRaw.get(i);
            if(l.startsWith(tag)) {

                String [] arr = l.split("[ ]+");
                sb.append(" ");

				int first = (arr.length >= 2 && arr[1].equals(Integer.toString(i-start+2)) ? 2 : 1);
                for (int j = first; j < arr.length; j++) {
                    sb.append(arr[j]);

                    if(j < arr.length-1) {
                        sb.append(" ");
                    }
                }

                indexLine++;
            } else {
                break;
            }
        }

        AbstractMap.SimpleEntry<String, Integer> siTextIndex = new AbstractMap.SimpleEntry<>(sb.toString(), indexLine);

        return siTextIndex;
    }

    private List<String> parseMultipleTimesOneLine(List<String> liRaw, int indexLine, String tag) throws ParseException {

        String l0 = liRaw.get(indexLine);
        if(!l0.startsWith(tag)) {
            throw new RuntimeException("Error in parsing " + tag);
        }
        String titleSub0 = l0.substring(tag.length()).trim();

        List<String> liTxt = new ArrayList<>();

        liTxt.add(titleSub0);

        for (int i = indexLine+1; i < liRaw.size(); i++) {

            String l = liRaw.get(i);
            if(l.startsWith(tag)) {
                StringBuilder sb = new StringBuilder();

                String [] arr = l.split("[ ]+");
                sb.append(" ");
                for (int j = 1; j < arr.length; j++) {
                    sb.append(arr[j]);

                    if(j < arr.length-1) {
                        sb.append(" ");
                    }
                }

                liTxt.add(sb.toString());
            } else {
                break;
            }
        }

        return liTxt;
    }

	private int parseCONECTLines(List<String> liRaw, int lineIndex, SortedList<int[]> bondList) throws ParseException {
		while (liRaw.get(lineIndex).startsWith(TAG_CONECT)) {
			String line = liRaw.get(lineIndex++);
			if (line.length() >= 16) {
				int atom1 = Integer.parseInt(line.substring(6,11).trim());
				int index = 16;
				while (line.length() >= index) {
					String s = line.substring(index-5, index).trim();
					if (s.isEmpty())
						break;
					int atom2 = Integer.parseInt(s);
					int[] atoms = new int[2];
					if (atom1 < atom2) {
						atoms[0] = atom1;
						atoms[1] = atom2;
					}
					else {
						atoms[0] = atom2;
						atoms[1] = atom1;
					}
					bondList.add(atoms);
					index += 5;
				}
			}
		}

		return lineIndex;
	}

	private static ListInteger<String> parseMultipleTimesMultipleLinesSEQRES(List<String> liRaw, int indexLine, String tag) throws ParseException {

        String l0 = liRaw.get(indexLine);
        if(!l0.startsWith(tag)) {
            throw new RuntimeException("Error in parsing " + tag);
        }

        StringBuilder sb = new StringBuilder();

        int start = indexLine;

        String chainId = l0.substring(11, 12);

        int numResidues = Integer.parseInt(l0.substring(13,17).trim());

        List<String> liChain = new ArrayList<>();

        for (int i = start; i < liRaw.size(); i++) {

            String l = liRaw.get(i);

            if(!l.startsWith(tag)) {
                break;
            }

            String chainIdLine = l.substring(11, 12);

            int numResiduesLine = Integer.parseInt(l.substring(13,17).trim());


            if(!chainId.equals(chainIdLine)) {

                String chain = sb.toString();

                liChain.add(chain);

                chainId = chainIdLine;

                numResidues = numResiduesLine;

                sb = new StringBuilder();
            }

            if(sb.length()>0){
                sb.append(" ");
            }
            /*
            if(numResidues!=numResiduesLine) {
                throw new RuntimeException("Number of residues differs!");
            }
			*/
            String chainLine = l.substring(19).trim();
            sb.append(chainLine);

            indexLine++;
        }

        if(sb.length()>0){
            String chain = sb.toString();

            liChain.add(chain);
        }

        ListInteger<String> listIndexLineChain = new ListInteger<String>(liChain, indexLine);

        return listIndexLineChain;

    }


}
