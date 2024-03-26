/*
 * Created on Dec 20, 2004
 *
 */
package com.actelion.research.chem.io;

import com.actelion.research.chem.Molecule3D;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;


/**
 * A parser is used to load and save molecules from the filesystem.
 * The save procedure is optional
 * 
 * @author freyssj
 */
@SuppressWarnings("resource")
public abstract class AbstractParser {

	protected static final String NEWLINE = System.getProperty("line.separator");

	protected List<String> errors = new ArrayList<String>();
	protected boolean optimize3D = true;
	
	/**
	 * 
	 * @param fileName
	 * @throws Exception
	 */
	public final List<Molecule3D> loadGroup(String fileName) throws Exception {
		Reader r = null;
		try {
			
			if(fileName.toUpperCase().endsWith(".GZ")) {
				GZIPInputStream is = new GZIPInputStream(new FileInputStream(fileName));
				r = new InputStreamReader(is, StandardCharsets.UTF_8);
			} else if(fileName.toUpperCase().endsWith(".ZIP")) {
				ZipInputStream is = new ZipInputStream(new FileInputStream(fileName));
				r = new InputStreamReader(is, StandardCharsets.UTF_8);
			} else {
				r = new BufferedReader(new FileReader(fileName)); 				
			}			
			List<Molecule3D> mols = loadGroup(fileName, r, -1, -1);
			int count = 0;
			for (Molecule3D mol: mols) {
				count++;
				if(mol!=null && (mol.getName()==null || mol.getName().length()==0)) mol.setName(new File(fileName).getName() + (mols.size()>1? "#"+ count:""));
			}			
			return mols;
		} finally {
			if(r!=null) try {r.close();} catch(Exception e){}
		}
	}
	
	public List<Molecule3D> loadGroup(String fileName, Reader in) throws Exception{
		return loadGroup(fileName, in, -1, -1);
	}
	public abstract List<Molecule3D> loadGroup(String fileName, Reader in, int from, int to) throws Exception;

	/**
	 *
	 * @param file
	 * @throws Exception
	 */
	public final Molecule3D load(File file) throws Exception {
		return load(file.getPath());
	}

	/**
	 * 
	 * @param fileName
	 * @throws Exception
	 */
	public final Molecule3D load(String fileName) throws Exception {
		List<Molecule3D> list = loadGroup(fileName);
		return list.size()>0? list.get(0): null;
	}

	/**
	 *
	 * @param fileName
	 * @param in
	 * @throws Exception
	 */
	public final Molecule3D load(String fileName, Reader in) throws Exception {
		List<Molecule3D> list = loadGroup(fileName, in, -1, -1);
		return list.size()>0? list.get(0): null;
	}
	
	

	/**
	 * 
	 * @param mol
	 * @param fileName
	 * @throws Exception
	 */ 
	public final void save(List<Molecule3D> mol, String fileName) throws Exception {
		Writer w = null;
		try {
			if(fileName.toUpperCase().endsWith(".GZ")) {
				GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(fileName));
				w = new OutputStreamWriter(os);				
			} else if(fileName.toUpperCase().endsWith(".ZIP")) {
				ZipOutputStream os = new ZipOutputStream(new FileOutputStream(fileName));
				w = new OutputStreamWriter(os);
			} else {
				w = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fileName), StandardCharsets.UTF_8));
			}
			
			save(mol, w);
			w.close();
		} finally {
			if(w!=null) try {w.close();} catch(Exception e){}
		}
	}
	
		
	/**
	 * 
	 * @param mol
	 * @param fileName
	 * @throws Exception
	 */
	public final void save(Molecule3D mol, String fileName) throws Exception {
		Writer w = null;
		try {
			if(fileName.toUpperCase().endsWith(".GZ")) {
				GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(fileName));
				w = new OutputStreamWriter(os);				
			} else if(fileName.toUpperCase().endsWith(".ZIP")) {
				ZipOutputStream os = new ZipOutputStream(new FileOutputStream(fileName));
				w = new OutputStreamWriter(os);
			} else {
				w = new FileWriter(fileName); 
			}
			
			save(mol, w);
			w.close();
		} finally {
			if(w!=null) try {w.close();} catch(Exception e){}
		}
	}
	
	/**
	 * 
	 * @param mol
	 * @param writer
	 * @throws Exception
	 */
	public void save(Molecule3D mol, Writer writer) throws Exception {
		throw new IllegalAccessError("Not implemented");
	}

	
	/**
	 * If not subclassed, save the files separately
	 * @param mols
	 * @param writer
	 * @throws Exception
	 */
	public void save(List<Molecule3D> mols, Writer writer) throws Exception {
		if(mols==null || mols.size()==0) return;
		if(mols.size()==1) save(mols.get(0), writer);
		else {
			throw new IllegalAccessError("Cannot save more than one file in  this format");
		}
		
	}

	public static final void convertDataToPrimitiveTypes(List<Molecule3D> res) {
		//Convert the field values to integer or double if possible
		if(res.size()==0) return;
		Set<String> fieldNames = new HashSet<String>();
		for (Molecule3D m : res) {
			fieldNames.addAll(m.getAuxiliaryInfos().keySet());
		}
		loop: for(String fieldName: fieldNames) {
			int type = 0;
			String val = "";
			for (Molecule3D m : res) {
				Object o = m.getAuxiliaryInfos().get(fieldName);
				if(o==null) continue;
				if(!(o instanceof String)) {
					continue loop; 
				}
				val = (String) o; 
				if(type==0) try {Integer.parseInt(val);} catch (Exception e) {type++;}
				if(type==1) try {Double.parseDouble(val);} catch (Exception e) {type++;}
				if(type==2) break; //String
			}
			if(type==0) {
				for (Molecule3D m : res) {
					val = (String)m.getAuxiliaryInfos().get(fieldName);
					if(val!=null) m.getAuxiliaryInfos().put(fieldName, Integer.parseInt((String)m.getAuxiliaryInfos().get(fieldName)));
				}
			} else if(type==1) {
				for (Molecule3D m : res) {
					val = (String)m.getAuxiliaryInfos().get(fieldName);
					if(val!=null) m.getAuxiliaryInfos().put(fieldName, Double.parseDouble(val));
				}
			}
		}				
	}
	

	protected static void writeR(Writer writer, String data, int len) throws IOException {
		if(data==null) data = "";
		int l = Math.max(0, len - data.length());
		for (int i = 0; i < l; i++)
			writer.write(' ');
		writer.write(data);
	}

	protected static void writeL(Writer writer, String data, int len) throws IOException {
		if(data==null) data="";
		writer.write(data);
		int l = Math.max(0, len - data.length());
		for (int i = 0; i < l; i++)
			writer.write(' ');
	}

	public List<String> getErrors() {
		return errors;
	}
	
	protected static boolean is3D(Molecule3D m) {
		for(int a=0; a<m.getAllAtoms();a++) {
			if(Math.abs(m.getAtomZ(a))>0.1) return true;
		}
		return false;
	}

	public void setOptimize3D(boolean optimize3d) {
		optimize3D = optimize3d;
	}
	
	public boolean isOptimize3D() {
		return optimize3D;
	}

	
}