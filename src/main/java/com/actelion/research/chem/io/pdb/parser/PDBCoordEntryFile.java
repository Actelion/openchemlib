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

import com.actelion.research.chem.Molecule3D;
import com.actelion.research.util.SortedList;

import java.util.*;

/**
 * PDBCoordEntryFile
 * Created by korffmo1 on 20.03.18.
 */
public class PDBCoordEntryFile {
	
	/**
	 *  description of pdb-file format: http://www.wwpdb.org/documentation/file-format
	 */

    //
    // Header information
    //
	
	/**
	*  Classification may be based on function, metabolic role, molecule type, cellular location, etc.
	*  This record can describe dual functions of a molecules, and when applicable, separated by a
	*  comma “
	*/
	
    private String classification;

    private String pdbID;

    private Date dateDeposition;

    private String obsolete;

    private String title;

    private String split;

    private String caveat;

    private String compound;

    private String source;

    private String keywords;

    private String expdata;
    private String nummdl;
    private String mdltyp;
    private String author;

    private List<String> liRevdat;

    private String sprsde;
    private List<String> liJRNL;

    private HashMap<Integer, String> hmNo_Remark;

    private List<String> liDBRef;

    private List<String> liDBRef1DBRef2;

    private List<String> liSEQADV;

    private List<String> liSEQRES;

    private List<String> liModRes;

    private List<String> liHet;

    private HashMap<String, String> hmId_Name;

    private HashMap<String, String> hmId_Synonyms;

    private HashMap<String, String> hmId_Formula;

    private List<String> liHelix;
    private List<String> liSheet;
    private List<String> liSSBond;
    private List<String> liLink;
    private List<String> liCisPep;

    private HashMap<String, String> hmId_Site;

    private String cryst1;
    private String origX1;
    private String origX2;
    private String origX3;
    private String scale1;
    private String scale2;
    private String scale3;

    private List<String> liMtrix1;
    private List<String> liMtrix2;
    private List<String> liMtrix3;
    
    private List<AtomRecord> protAtomRecords;
    private List<AtomRecord> hetAtomRecords;


    private SortedList<int[]> liConnect;

    private String master;

    private boolean end;

    public String getClassification() {
        return classification;
    }

    public void setClassification(String classification) {
        this.classification = classification;
    }

    public String getID() {
        return pdbID;
    }

    public void setID(String pdbID) {
        this.pdbID = pdbID;
    }

    public Date getDateDeposition() {
        return dateDeposition;
    }

    public void setDateDeposition(Date dateDeposition) {
        this.dateDeposition = dateDeposition;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }
    
    public List<AtomRecord> getProtAtomRecords() {
    	return protAtomRecords;
    }
    
    public void setProtAtomRecords(List<AtomRecord> protAtomRecords) {
    	this.protAtomRecords = protAtomRecords;
    }
    
    public List<AtomRecord> getHetAtomRecords() {
    	return hetAtomRecords;
    }
    
    public void setHetAtomRecords(List<AtomRecord> hetAtomRecords) {
    	this.hetAtomRecords = hetAtomRecords;
    }
  

    public String getObsolete() {
        return obsolete;
    }

    public void setObsolete(String obsolete) {
        this.obsolete = obsolete;
    }

    public String getSplit() {
        return split;
    }

    public void setSplit(String split) {
        this.split = split;
    }

    public String getCaveat() {
        return caveat;
    }

    public void setCaveat(String caveat) {
        this.caveat = caveat;
    }

    public String getCompound() {
        return compound;
    }

    public void setCompound(String compound) {
        this.compound = compound;
    }

    public String getSource() {
        return source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public String getKeywords() {
        return keywords;
    }

    public void setKeywords(String keywords) {
        this.keywords = keywords;
    }

    public String getExpdata() {
        return expdata;
    }

    public void setExpdata(String expdata) {
        this.expdata = expdata;
    }

    public String getNummdl() {
        return nummdl;
    }

    public void setNummdl(String nummdl) {
        this.nummdl = nummdl;
    }

    public String getMdltyp() {
        return mdltyp;
    }

    public void setMdltyp(String mdltyp) {
        this.mdltyp = mdltyp;
    }

    public String getAuthor() {
        return author;
    }

    public void setAuthor(String author) {
        this.author = author;
    }

    public List<String> getRevdat() {
        return liRevdat;
    }

    public void setRevdat(List<String> liRevdat) {
        this.liRevdat = liRevdat;
    }

    public String getSprsde() {
        return sprsde;
    }

    public void setSprsde(String sprsde) {
        this.sprsde = sprsde;
    }

    public List<String> getJrnl() {
        return liJRNL;
    }

    public void setJrnl(List<String> jrnl) {
        this.liJRNL = jrnl;
    }

    public String getRemark0() {
        return hmNo_Remark.get(0);
    }


    public String getRemark1() {
        return hmNo_Remark.get(1);
    }


    public String getRemark2() {
        return hmNo_Remark.get(2);
    }


    public String getRemark3() {
        return hmNo_Remark.get(3);
    }

    public void setRemarks(HashMap<Integer, String> hmNo_Remark) {
        this.hmNo_Remark = hmNo_Remark;
    }

    public String getRemark(int n) {
        return hmNo_Remark.get(n);
    }


    public List<Integer> getRemarks(){
        List<Integer> liRemarkNo = new ArrayList<>(hmNo_Remark.keySet());

        Collections.sort(liRemarkNo);

        return liRemarkNo;
    }

    public List<String> getDBRef() {
        return liDBRef;
    }

    public void setDBRef(List<String> liDBRef) {
        this.liDBRef = liDBRef;
    }

    public List<String> getDBRef1DBRef2() {
        return liDBRef1DBRef2;
    }

    public void setDBRef1DBRef2(List<String> liDBRef1DBRef2) {
        this.liDBRef1DBRef2 = liDBRef1DBRef2;
    }

    public List<String> getSEQADV() {
        return liSEQADV;
    }

    public void setSEQADV(List<String> liSEQADV) {
        this.liSEQADV = liSEQADV;
    }

    public List<String> getSEQRES() {
        return liSEQRES;
    }

    public void setSEQRES(List<String> liSEQRES) {
        this.liSEQRES = liSEQRES;
    }

    public List<String> getModRes() {
        return liModRes;
    }

    public void setModRes(List<String> liModRes) {
        this.liModRes = liModRes;
    }

    public List<String> getHet() {
        return liHet;
    }

    public void setHet(List<String> liHet) {
        this.liHet = liHet;
    }

    public HashMap<String, String> getHmId_Name() {
        return hmId_Name;
    }

    public void setHmId_Name(HashMap<String, String> hmId_Name) {
        this.hmId_Name = hmId_Name;
    }

    public List<String> getNameIDs(){
        List<String> li = new ArrayList<>(hmId_Name.keySet());
        Collections.sort(li);
        return li;
    }

    public String getName(String nameId){
        return hmId_Name.get(nameId);
    }

    public void setHmId_Synonyms(HashMap<String, String> hm) {
        this.hmId_Synonyms = hm;
    }

    public List<String> getSynonymIDs(){
        List<String> li = new ArrayList<>(hmId_Synonyms.keySet());
        Collections.sort(li);
        return li;
    }

    public String getSynonyms(String synonymId){
        return hmId_Synonyms.get(synonymId);
    }

    public void setHmId_Formula(HashMap<String, String> hm) {
        this.hmId_Formula = hm;
    }

    public List<String> getFormulaIDs(){
        List<String> li = new ArrayList<>(hmId_Formula.keySet());
        Collections.sort(li);
        return li;
    }

    public String getFormula(String formulaId){
        return hmId_Formula.get(formulaId);
    }

    public List<String> getHelix() {
        return liHelix;
    }

    public void setHelix(List<String> liHelix) {
        this.liHelix = liHelix;
    }

    public List<String> getSheet() {
        return liSheet;
    }

    public void setSheet(List<String> liSheet) {
        this.liSheet = liSheet;
    }

    public List<String> getSSBond() {
        return liSSBond;
    }

    public void setSSBond(List<String> liSSBond) {
        this.liSSBond = liSSBond;
    }

    public List<String> getLink() {
        return liLink;
    }

    public void setLink(List<String> liLink) {
        this.liLink = liLink;
    }

    public List<String> getCisPep() {
        return liCisPep;
    }

    public void setCisPep(List<String> liCisPep) {
        this.liCisPep = liCisPep;
    }

    public void setHmId_Site(HashMap<String, String> hm) {
        this.hmId_Site = hm;
    }

    public List<String> getSiteIDs(){
        List<String> li = new ArrayList<>(hmId_Site.keySet());
        Collections.sort(li);
        return li;
    }

    public String getSite(String id){
        return hmId_Site.get(id);
    }

    public String getCryst1() {
        return cryst1;
    }

    public void setCryst1(String cryst1) {
        this.cryst1 = cryst1;
    }

    public String getOrigX1() {
        return origX1;
    }

    public void setOrigX1(String origX1) {
        this.origX1 = origX1;
    }

    public String getOrigX2() {
        return origX2;
    }

    public void setOrigX2(String origX2) {
        this.origX2 = origX2;
    }

    public String getOrigX3() {
        return origX3;
    }

    public void setOrigX3(String origX3) {
        this.origX3 = origX3;
    }

    public String getScale1() {
        return scale1;
    }

    public void setScale1(String scale1) {
        this.scale1 = scale1;
    }

    public String getScale2() {
        return scale2;
    }

    public void setScale2(String scale2) {
        this.scale2 = scale2;
    }

    public String getScale3() {
        return scale3;
    }

    public void setScale3(String scale3) {
        this.scale3 = scale3;
    }

    public List<String> getMtrix1() {
        return liMtrix1;
    }

    public void setMtrix1(List<String> liMtrix1) {
        this.liMtrix1 = liMtrix1;
    }

    public List<String> getMtrix2() {
        return liMtrix2;
    }

    public void setMtrix2(List<String> liMtrix2) {
        this.liMtrix2 = liMtrix2;
    }

    public List<String> getMtrix3() {
        return liMtrix3;
    }

    public void setMtrix3(List<String> liMtrix3) {
        this.liMtrix3 = liMtrix3;
    }

    public SortedList<int[]> getLiConnect() {
        return liConnect;
    }

    public void setLiConnect(SortedList<int[]> liConnect) {
        this.liConnect = liConnect;
    }

    public String getMaster() {
        return master;
    }

    public void setMaster(String master) {
        this.master = master;
    }

    public boolean isEnd() {
        return end;
    }

    public void setEnd(boolean end) {
        this.end = end;
    }
    
    public Map<String,List<Molecule3D>> extractMols() {
        return extractMols(false);
    }

    public Map<String,List<Molecule3D>> extractMols(boolean detachCovalentLigands) {
        StructureAssembler assembler = new StructureAssembler(liConnect,protAtomRecords,hetAtomRecords);
        assembler.setDetachCovalentLigands(detachCovalentLigands);
        return assembler.assemble();
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("PDBCoordEntryFile{");
        sb.append("classification='").append(classification).append('\'');
        sb.append(", pdbID'").append(pdbID).append('\'');
        sb.append(", dateDeposition=").append(dateDeposition);
        sb.append(", obsolete='").append(obsolete).append('\'');
        sb.append(", title='").append(title).append('\'');
        sb.append(", split='").append(split).append('\'');
        sb.append(", caveat='").append(caveat).append('\'');
        sb.append(", compound='").append(compound).append('\'');
        sb.append(", source='").append(source).append('\'');
        sb.append(", keywords='").append(keywords).append('\'');
        sb.append(", expdata='").append(expdata).append('\'');
        sb.append(", nummdl='").append(nummdl).append('\'');
        sb.append(", mdltyp='").append(mdltyp).append('\'');
        sb.append(", author='").append(author).append('\'');
        sb.append(", revdat='").append(liRevdat).append('\'');
        sb.append(", sprsde='").append(sprsde).append('\'');
        sb.append(", jrnl='").append(liJRNL).append('\'');

        List<Integer> liKeyRemark = getRemarks();

        for (int key : liKeyRemark) {

            String s = getRemark(key);

            sb.append(", remark");
            sb.append(key);
            sb.append("='").append(s).append('\'');;
        }

        sb.append('}');
        return sb.toString();
    }
}
