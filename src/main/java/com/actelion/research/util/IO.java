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

package com.actelion.research.util;


import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import com.actelion.research.io.StringReadChannel;
import com.actelion.research.util.datamodel.DoubleArray;
import com.actelion.research.util.datamodel.IntArray;

/**
 * IO
 * 2003 MvK: Start implementation
 */
public class IO {

	// private static final NumberFormat NF = new DecimalFormat("000");

	private static final int LIMIT_FILES_DIR = 10000;

	public static final String SEP = System.getProperty("file.separator");



	public static boolean canWriteAndDeleteInPath(File dir) throws IOException {
		boolean w = false;

		File f = File.createTempFile("test", ".txt", dir);

		if(f.canWrite()){
			w=true;
		}

		w = f.delete();

		return w;
	}


	/**
	 * Do not forget to close BufferedReader
	 * 
	 * @param sAbsolutePathIn
	 *            path
	 * @return BufferedReader
	 */
	public static BufferedReader getBufferedReader(String sAbsolutePathIn) throws FileNotFoundException {

		BufferedReader bufferedReader = null;

		if (sAbsolutePathIn.length() > 0) {
			FileInputStream fis = new FileInputStream(sAbsolutePathIn);
			InputStreamReader isr = new InputStreamReader(fis, StandardCharsets.UTF_8);
			bufferedReader = new BufferedReader(isr);
		}

		return bufferedReader;
	}


	/**
	 * A number is added to the base name of the file
	 * 
	 * @param totalpath
	 * @return
	 */
	public static String getUniqueFileName(String totalpath) {
		
		File fi = new File(totalpath);

		return getUniqueFileName(fi, fi.getParentFile(), Formatter.dfI3).getAbsolutePath();

	}
	
	public static File getUniqueFileName(File fiIn) {
		return getUniqueFileName(fiIn, fiIn.getParentFile());
	}
	
	public static File getUniqueFileName(String sFilename, String sDirDestination) {
		return getUniqueFileName(new File(sFilename), new File(sDirDestination), Formatter.dfI3);
	}
	
	public static File getUniqueFileName(File file, File dirDestination) {
		return getUniqueFileName(file, dirDestination, Formatter.dfI3);
	}
	
	/**
	 * If the file does not exists the input file is returned.
	 * @param file
	 * @param dirDestination
	 * @param df
	 * @return
	 */
	public static File getUniqueFileName(File file, File dirDestination, DecimalFormat df) {

		String sFileName = file.getName();
		
		File fiNew = new File(dirDestination, sFileName);
		
		if(!fiNew.exists()){
			return fiNew;
		}
		
		String sBaseName = getBaseName(file);
		
		String sExtension = getExtension(file);
		
		String sBaseNameNumber = sBaseName;

		int cc = 1;

		while (fiNew.isFile()) {

			sBaseNameNumber = next(sBaseNameNumber, df);
			
			sFileName = sBaseNameNumber + sExtension;
			
			fiNew = new File(dirDestination, sFileName);

			cc++;

			if (cc > LIMIT_FILES_DIR) {
				sFileName = null;
				IOException ex = new IOException("To many files in this dir.");
				ex.printStackTrace();
				break;
			}
		}
		
		return fiNew;
	}

	public static File getUniqueUserDir() throws IOException {
    	String sWorkdir = System.getProperty("user.home");
    	File dirParent = new File(sWorkdir);
    	return getUniqueDateDir(dirParent);
	}
	
	public static File getUniqueDateDir(File dirParent) throws IOException {
		return getUniqueDateDir(dirParent, ""); 
	}
	
	public static File getUniqueDateDir(String appendix) throws IOException {
    	String sWorkdir = System.getProperty("user.home");
    	File dirParent = new File(sWorkdir);
    	return getUniqueDateDir(dirParent, appendix);
	}
	public static File getUniqueDateDir(File dirParent, String appendix) throws IOException {
		
		Date d = new Date();
		DateFormat df = new SimpleDateFormat("yyMMdd");

		String sDate = df.format(d);
		
		String name = sDate + appendix;
		
		File dir = new File(dirParent, name);
		
		int index = 1;
		while(dir.exists()){
			name = sDate + "_" + index + "_" + appendix;
			index++;
			dir = new File(dirParent, name);
		}
		
		if(!dir.mkdirs()){
			String e = "Not possible to make dir " + dir.getAbsolutePath() + ".";
			throw new IOException(e);
		}
		
		return dir;
	}
	
	public static File getUniqueDir(File dirParent, String suffix) throws IOException {
		
		String name = suffix;
		
		File dir = new File(dirParent, name);
		
		int index = 1;
		while(dir.exists()){
			name = suffix+index;
			index++;
			dir = new File(dirParent, name);
		}
		
		if(!dir.mkdirs()){
			String e = "Not possible to make dir " + dir.getAbsolutePath() + ".";
			throw new IOException(e);
		}
		
		return dir;
	}

	/**
	 * Has to be the total path of the file or there will be errors.
	 * @param totalpath
	 * @return
	 */
	public static String getNextFileName(String totalpath) {

		File fi = new File(totalpath);

		return getUniqueFileName(fi, fi.getParentFile(), Formatter.dfI3).getAbsolutePath();
	}
	
	public static String next(String txt) {
		return next(txt, Formatter.dfI3);
	}
	
	/**
	 * 
	 * @param txt String with an integer num at the end. (blabla567) 
	 * @return Adds one to the last number (blabla568). 
	 * If no number in text 001 is added. 
	 */
	public static String next(String txt, DecimalFormat dfExtern) {
		String s = "";

		int ind = -1;
		for (int i = txt.length() - 1; i >= 0; i--) {
			if (!Character.isDigit(txt.charAt(i))) {
				ind = i;
				break;
			}
		}

		s = txt.substring(0, ind + 1);

		String sN = txt.substring(ind + 1, txt.length());
		
		NumberFormat nfName = dfExtern;
		
		int num = 0;
		
		if (sN.length() > 0) {
			num = Integer.parseInt(sN);
			
			String patFormat = "";
			for (int i = 0; i < sN.length(); i++) {
				patFormat += "0";
			}
			
			nfName = new DecimalFormat(patFormat);
		}
		
		num++;

		
		s += nfName.format(num);
		

		return s;
	}

	public static BufferedWriter getBuffWrite(String sAbsolutePathOut,
			boolean bAppend) throws IOException {

		BufferedWriter bufferedWriter = null;
		if (sAbsolutePathOut.length() > 0) {
			FileOutputStream fos = new FileOutputStream(sAbsolutePathOut, bAppend);
			OutputStreamWriter isw = new OutputStreamWriter(fos);
			bufferedWriter = new BufferedWriter(isw);
		} else {
			OutputStreamWriter isw = new OutputStreamWriter(System.out);
			bufferedWriter = new BufferedWriter(isw);
		}
		return bufferedWriter;
	}

	/**
	 * 
	 * @param str
	 * @return base name without extension.
	 */
	public static String getBaseName(String str) {

		File fi = new File(str);

		String name = fi.getName();

		int iIndexStart = name.lastIndexOf('\\');
		if(iIndexStart==-1){
			iIndexStart = name.lastIndexOf(File.separator);
		}
		if(iIndexStart==-1){
			iIndexStart = 0;
		}

		int iIndexEnd = name.lastIndexOf('.');

		String sBaseName = "";
				
		if(iIndexEnd == -1)
			iIndexEnd = name.length();
				
		sBaseName = name.substring(iIndexStart, iIndexEnd);
		
		return sBaseName;
	}

	public static String getBaseName(File file) {
		return getBaseName(file.getAbsolutePath());
	}

	/**
	 * 
	 * @param file
	 * @return the part after the the last '.' inclusive the '.'. 
	 * Returns String with length null when no extension found.
	 */
	public static String getExtension(File file) {
		
		String sName = file.getName();
		
		int iIndexStartName = sName.lastIndexOf(SEP) + 1;
		
		int iIndexStartExtension = sName.lastIndexOf('.');
		
		if(iIndexStartExtension==-1){
			return "";
		}
		
		String sExtension = sName.substring(iIndexStartExtension);
		
		return sExtension;
	}

	public static void mkdirs(String path) throws IOException {
		
		mkdirs(new File(path));

	}
	
	public static void mkdirs(File dir) throws IOException {
		
		if (!dir.exists()) {
			if(!dir.mkdirs()){
				throw new IOException("Not possible to make dir " + dir.getAbsolutePath() + ".");
			}
		}
		
	}
	
	public static void readBetweenTags(String sAbsolutePathIn,
		String sTagStartRegEx, 
		String sTagEndRegEx, 
		Vector<String> vecStringContent) {

		try {
			StringReadChannel channel = new StringReadChannel(Channels.newChannel(new FileInputStream(new File(sAbsolutePathIn))));
			
			boolean bTagStartFound = false;
			
			boolean bTagEndFound = false;
			
			while (channel.hasMoreLines()) {
				String sLine = channel.readLine();

				if (sLine != null) {
					if (sLine.matches(sTagStartRegEx)) {
						bTagStartFound = true;
						break;
					}
				}
			}

			if (bTagStartFound == true) {
				while (channel.hasMoreLines()) {
					String sLine = channel.readLine();
					if (sLine != null) {
						if (sLine.matches(sTagEndRegEx)) {
							bTagEndFound = true;
							break;
						} else {
							vecStringContent.addElement(sLine);
						}
					}
				}
			} else {
				System.err.println("Tag: " + sTagStartRegEx + " not found");
				(new RuntimeException()).printStackTrace();
			}

			channel.close();
			
			if (!bTagEndFound) {
				System.err.println("Tag: " + sTagEndRegEx + " not found");
				(new RuntimeException()).printStackTrace();
			}
			
		} catch (IOException e) {
			System.err.println(e);
		} catch (Exception e) {
			System.err.println(e);
		}

	}

	/**
	 * Reads all lines after a given tag and stores the lines as Strings in a
	 * vector object.
	 * 
	 * @param sAbsolutePathIn
	 *            path of the input file
	 * @param sTagRegEx
	 *            the tag as regular expression
	 * @param vecStringContent
	 *            contains the result of the read in.
	 */
	public static void readFromTag(String sAbsolutePathIn, String sTagRegEx, Vector<String> vecStringContent) {
		
		try {
			StringReadChannel channel = new StringReadChannel(Channels.newChannel(new FileInputStream(new File(sAbsolutePathIn))));
			boolean bTagFound = false;
			String sLine = "";
			
			while (channel.hasMoreLines()) {
				sLine = channel.readLine();

				if (sLine != null) {
					if (sLine.matches(sTagRegEx)) {
						bTagFound = true;
						break;
					}
				}
			}
			if (bTagFound == true) {
				while (channel.hasMoreLines()) {
					sLine = channel.readLine();

					if (sLine != null) {
						if (sLine.length() < 1) {
							break;
						} else {
							vecStringContent.addElement(sLine);
						}
					}
				}
			} else {
				System.err.println("Tag: " + sTagRegEx + " not found");
				(new RuntimeException()).printStackTrace();
			}

			channel.close();
		} catch (IOException e) {
			System.err.println(e);
		} catch (Exception e) {
			System.err.println(e);
		}
	}
	
    public static void skipUntilLineMatchesRegEx(InputStream in, String regex)throws Exception {
    	
    	int limit = 10000;
    	int cc=0;
    	String line = readLine(in);
    	while(!line.matches(regex)){
    		line = readLine(in);
    		if(cc==limit) {
    			throw new Exception("Limit of " + limit + " lines exceeded.");
    		}
    		cc++;
    	}
    }

	public static String read(InputStream is) throws IOException {
		BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
		String line=null;
		StringBuilder sb = new StringBuilder();
		while ((line = reader.readLine())!=null) {
			if(sb.length()>0)
				sb.append("\n");
			sb.append(line);
		}
		reader.close();
		return sb.toString();
	}
	
	public static List<Integer> readListInteger(File fiTxt) throws IOException {
		Scanner scanner = new Scanner(fiTxt);
	    List<Integer> li = new ArrayList<Integer>();
	    while (scanner.hasNextInt()) {
	        li.add(scanner.nextInt());
	    }
	    scanner.close();
	    return li;
	}

	/**
	 * Reads a file that contains one int per line.
	 * @param fiIntLineWise
	 * @return
	 * @throws IOException
	 */
	public static int [] readLines2IntArray(File fiIntLineWise) throws IOException{

		IntArray intArray = new IntArray();

		BufferedReader br = new BufferedReader(new FileReader(fiIntLineWise));

		String line = "";

		while((line=br.readLine()) != null){
			int val = Integer.parseInt(line);

			intArray.add(val);
		}

		br.close();

		return intArray.get();
	}

	public static double [] readLines2DoubleArray(File fiIntLineWise) throws IOException{

		DoubleArray doubleArray = new DoubleArray();

		BufferedReader br = new BufferedReader(new FileReader(fiIntLineWise));

		String line = "";

		while((line=br.readLine()) != null){
			double v = Double.parseDouble(line);

			doubleArray.add(v);
		}

		br.close();

		return doubleArray.get();
	}

	public static String readLine(InputStream is) throws IOException {
		StringBuilder sb = new StringBuilder();
		int c = -1;
		while((c = is.read()) != '\n'){
			
			if(c==-1)
				break;
			
			sb.append((char)c);
		}
		
    	String str = StringFunctions.removeCharacter(sb, '\r');
    	
    	return str;

	}


	public static String read(File file) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line=null;
		StringBuilder sb = new StringBuilder();
		while ((line = reader.readLine())!=null) {
			if(sb.length()>0)
				sb.append("\n");
			sb.append(line);
		}
		reader.close();
		return sb.toString();
	}




    public static String readLine(Reader is) throws IOException {
		StringBuffer b = new StringBuffer();
		int c=0;
		c=is.read();
		while((c != -1)&&(c != '\n')) {
			b.append((char)c);
			c=is.read();
		}
		return b.toString();
	}
    

	
	public static String readLine(FileChannel fc) throws IOException {
		ByteBuffer buf = ByteBuffer.allocate(1000);
		
		StringBuilder sb = new StringBuilder();
		
		boolean bEOF=false;
		boolean bEOLine=false;
		
		long posFCPos = fc.position();
		
		
		while (!bEOF){
			buf.position(0);
			int size = fc.read(buf);
			if(size==-1){
				break;
			}
			
			for (int i = 0; i < size; i++) {
			
				int c = buf.get(i);
				posFCPos++;
				
				if(c==-1) {
					bEOF=true;
					
				} else {
					if(c=='\n') {
						fc.position(posFCPos);
						bEOLine=true;
						break;
					} else {
						sb.append((char)c);
					}
				}
			}
			if(bEOLine)
				break;
		}
		
		return sb.toString();
	}


	
	public static String readFirstLine(File file) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = reader.readLine();
		reader.close();
		return line;
	}
	public static List<String> readLines2List(File file) throws IOException {
		List<String> li = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(new FileReader(file));

		String line = null;

		while ((line = reader.readLine())!=null) {
			li.add(line);
		}

		reader.close();

		return li;
	}

	/**
	 * The stream is not closed.
	 * @param is
	 * @return
	 * @throws IOException
	 */
	public static List<String> readLines2List(InputStream is) throws IOException {
		List<String> li = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));

		String line = null;

		while ((line = reader.readLine())!=null) {
			li.add(line);
		}

		return li;
	}




	public static List<String> readLines2List(List<File> liFile) throws IOException {
		
		List<String> li = new ArrayList<String>();
		
		for (File fi : liFile) {
			li.addAll(readLines2List(fi));
		}
		
		return li;
	}
	


	public static void write(String sAbsolutePathOut, String sContent) {
		write(sAbsolutePathOut, sContent, true);
	}


	public static void write(String sAbsolutePathOut, String sContent, boolean bAppend) {
		write(new File(sAbsolutePathOut), sContent, bAppend);
	}
	
	public static void write(File file, String sContent) {
		write(file, sContent, false);
	}

	public static void write(File file, String sContent, boolean bAppend) {
		try {
			BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file, bAppend));
			bufferedWriter.write(sContent);
			bufferedWriter.close();
		} catch (IOException ex) {
			throw new RuntimeException(ex);
		}
	}
	
	public void write2Channel(FileChannel fc, String str) throws IOException {
    	ByteBuffer buf = ByteBuffer.allocate(str.length());
        byte[] bytes = str.getBytes();
        buf.put(bytes);
        buf.flip();
        fc.write(buf);
    }

	/**
	 * Writes each string in a separate line
	 * @param file
	 * @param li
	 * @throws IOException
	 */
	public static void write(File file, List<String> li) throws IOException {
		
		FileWriter fw = new FileWriter(file);
		for (int i = 0; i < li.size(); i++) {
			fw.write(li.get(i));
			if(i<li.size()-1)
				fw.write("\n");	
		}
		fw.close();
	}
	
	public static void writeIntegerList(File file, List<Integer> li) throws IOException {
		
		FileWriter fw = new FileWriter(file);
		for (int i = 0; i < li.size(); i++) {
			fw.write(li.get(i).toString());
			if(i<li.size()-1)
				fw.write("\n");	
		}
		fw.close();
	}

	public static void write(File file, int [] arr) throws IOException {

		BufferedWriter fw = new BufferedWriter(new FileWriter(file));

		for (int i = 0; i < arr.length; i++) {
			fw.write(Integer.toString(arr[i]));
			if(i<arr.length-1)
				fw.write("\n");
		}
		fw.close();
	}

	public static void write(File file, double [] arr) throws IOException {

		BufferedWriter fw = new BufferedWriter(new FileWriter(file));

		for (int i = 0; i < arr.length; i++) {
			fw.write(Double.toString(arr[i]));
			if(i<arr.length-1)
				fw.write("\n");
		}
		fw.close();
	}

	public static void saveObject(String file, Object o) throws IOException {
		FileOutputStream out = new FileOutputStream(file);
		ObjectOutputStream s = new ObjectOutputStream(out);
		s.writeObject(o);
		out.close();
	}

	public static Object loadObject(String file) throws IOException, ClassNotFoundException {
		FileInputStream in = new FileInputStream(file);
		ObjectInputStream s = new ObjectInputStream(in);
		Object o = s.readObject();
		in.close();
		return o;
	}

	public static String readerToString(Reader reader) throws IOException {
		return  readerToString(reader, Integer.MAX_VALUE);
	}

	public static String readerToString(Reader reader, int maxSize) throws IOException {
		if(maxSize<=0) maxSize = Integer.MAX_VALUE;
		char[] buf = new char[512];
		int c;
		StringBuilder sb = new StringBuilder();
		while(sb.length()<maxSize && ( c = reader.read(buf, 0, Math.min(buf.length, maxSize-sb.length()))) > 0) {
			sb.append(buf, 0, c);
		}
		return sb.toString();
	}

	public static String streamToString(InputStream is) throws IOException {
		byte[] buf = new byte[512];
		int c;
		StringBuilder sb = new StringBuilder();
		while(( c = is.read(buf)) > 0) {
			sb.append(new String(buf, 0, c));
		}
		return sb.toString();
	}

	public static byte[] getBytes(File f) throws IOException {
		FileInputStream is = new FileInputStream(f);
		byte[] res = new byte[(int) f.length()];
		is.read(res, 0, res.length);
		is.close();
		return res;
	}

	public static void bytesToFile(byte[] bytes, File f) throws IOException {
		FileOutputStream os = new FileOutputStream(f);
		os.write(bytes);
		os.close();
	}

	public static void stringToFile(String s, File f) throws IOException {
		try (FileWriter os = new FileWriter(f)) {
			os.write(s);
		}
	}

	public static String fileToString(File f) throws IOException {
		return fileToString(f, Integer.MAX_VALUE);
	}

	public static String fileToString(File f, int maxSize) throws IOException {
		try (FileReader reader = new FileReader(f)) {
			return readerToString(reader, maxSize);
		}
	}

	/**
	 * Redirects the streams, without closing them
	 * @param is
	 * @param os
	 * @throws IOException
	 */
	public static void redirect(InputStream is, OutputStream os) throws IOException {
		byte[] buf = new byte[512];
		int c;
		while((c=is.read(buf))>0) {
			os.write(buf, 0, c);
		}
		is.close();
	}

	/**
	 * Redirects the reader to the writer, without closing them
	 * @param is
	 * @param os
	 * @throws IOException
	 */
	public static void redirect(Reader is, Writer os) throws IOException {
		char[] buf = new char[512];
		int c;
		while((c=is.read(buf))>0) {
			os.write(buf, 0, c);
		}
		is.close();
	}

	public static void copy(File src, File dest) throws IOException {
		InputStream is = null;
		OutputStream os = null;
		try {
			is = new BufferedInputStream(new FileInputStream(src));
			os = new BufferedOutputStream(new FileOutputStream(dest));
			redirect(is, os);

		} finally {
			try {if(is!=null) is.close();}catch(Exception e){}
			try {if(os!=null) os.close();}catch(Exception e){}
		}
	}

	/**
	 * Debug instruction to view the relevant stacktrace
	 * @throws IOException
	 */
	public static void dumpStack() {
		Exception e = new Exception("StackTrace");
		Set<Throwable> dejaVu = Collections.newSetFromMap(new IdentityHashMap<Throwable, Boolean>());
		dejaVu.add(e);

		synchronized (System.err) {
			StackTraceElement[] trace = e.getStackTrace();
			System.err.println("Exception "+e.getMessage());
			for (StackTraceElement traceElement : trace) {
				if(traceElement.getClassName().contains("com.actelion.research.util.IOUtils")) continue;
				if(traceElement.getClassName().contains("com.actelion")) {
					System.err.println("\tat " + traceElement);
				}
			}
			if(e.getCause()!=null) {
				System.err.println("Caused by "+e.getCause());
				for (StackTraceElement traceElement : trace) {
					if(traceElement.getClassName().contains("com.actelion.research.util.IOUtils.dumpStack")) continue;
					if(traceElement.getClassName().contains("com.actelion")) {
						System.err.println("\tat " + traceElement);
					}
				}
			}
		}

	}

}