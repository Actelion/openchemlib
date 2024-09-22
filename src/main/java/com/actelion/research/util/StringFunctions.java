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

import com.actelion.research.util.datamodel.DoubleArray;

import java.awt.Point;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class StringFunctions {

    // http://stackoverflow.com/questions/4731055/whitespace-matching-regex-java
	public static final String PAT_WHITESPACE = "[\\s\\u0085\\p{Z}]";

    // http://stackoverflow.com/questions/1805518/replacing-all-non-alphanumeric-characters-with-empty-strings
	public static final String PAT_NOT_ALPHANUMERIC = "[^\\p{IsAlphabetic}^\\p{IsDigit}]";


	public static final String [] REGEX_META_CHARACTERS = {"*","%","@","&","+", "(", ")"};

	public static final String SEP = "; ";



	public static double [] parse2Double(String s, String sepRegEx){
		String [] a = s.split(sepRegEx);
		double [] arrDouble = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			arrDouble[i]=Double.parseDouble(a[i]);
		}
		return arrDouble;
	}


	public static String getAppendedSorted(String s1, String s2) {

		StringBuilder sb = new StringBuilder();

		if(s1.compareTo(s2)>0){
			sb.append(s1);
			sb.append(s2);
		} else {
			sb.append(s2);
			sb.append(s1);
		}

		return sb.toString();
	}

	public static String min(String s1, String s2) {

		StringBuilder sb = new StringBuilder();

		if(s1.compareTo(s2)<0){
			sb.append(s1);
		} else {
			sb.append(s2);
		}

		return sb.toString();
	}

	public static String max(String s1, String s2) {

		StringBuilder sb = new StringBuilder();

		if(s1.compareTo(s2)>0){
			sb.append(s1);
		} else {
			sb.append(s2);
		}

		return sb.toString();
	}

	public static int countIntegerInText(String txt) {

		int nInt = 0;

		boolean started=false;

		for(int i=0; i<txt.length(); i++) {

			char c = txt.charAt(i);

			if(Character.isDigit(c)){
				if(!started){
					nInt++;
					started=true;
				}
			} else {
				started=false;
			}
		}

		return nInt;
	}

	public static int countWordInText(String txt) {

		int nWords = 0;

		boolean started=false;

		for(int i=0; i<txt.length(); i++) {

			char c = txt.charAt(i);

			if(Character.isLetter(c)){
				if(!started){
					nWords++;
					started=true;
				}
			} else {
				started=false;
			}
		}

		return nWords;
	}


	public static boolean equal(byte [] b1, byte [] b2) {

		if(b1 == null && b2==null) {
			return true;
		}

		if(b1 != null && b2==null) {
			return false;
		}

		if(b1 == null && b2!=null) {
			return false;
		}

		String s1 = new String(b1);
		String s2 = new String(b2);

		if(!s1.equals(s2)) {
			return false;
		}

		return true;
	}


	public static String encodeHTML(String txt) {
		
	    StringBuilder sb = new StringBuilder();
	    
	    for(int i=0; i<txt.length(); i++) {
	    	
	        char c = txt.charAt(i);
	        
	        if(c > 127 || c=='"' || c=='<' || c=='>') {
	           sb.append("&#"+(int)c+";");
	        } else {
	            sb.append(c);
	        }
	    }
	    
	    return sb.toString();
	}


	public static Comparator<String> getComparatorLength(){

		return new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {

				int cmp = 0;

				if(o1.length() > o2.length()){
					cmp = 1;
				}else if(o1.length() < o2.length()){
					cmp = -1;
				}

				return cmp;
			}
		};
	}



	public static DecimalFormat getDecimalFormat(int precision){
		
		StringBuilder sbFormPat = new StringBuilder("0");
		
		if(precision>0)
			sbFormPat.append(".");
			
		for (int i = 0; i < precision; i++) {
			sbFormPat.append("0");
		}
		
		return new DecimalFormat(sbFormPat.toString());
	}
	
	public static List<String> getAllSubStrings(String str, int minsize){
		
		TreeSet<String> ts = new TreeSet<String>();
		
		for (int i = 0; i < str.length(); i++) {
			for (int j = i+minsize; j < str.length()+1; j++) {
				String pat = str.substring(i, j);
				ts.add(pat);
			}
		}
		return new ArrayList<String>(ts);
	}
	
	/**
	 * 
	 * @param min minimum length
	 * @param max maximum length
	 * @return
	 */
	public static String getRandom(int min, int max){
		Random rnd = new Random();
		
		int s = rnd.nextInt(max-min)+min;
		
		StringBuilder sb = new StringBuilder();
		
		int startASCII = 97;
		
		int endASCII = 122+1;
		
		for (int i = 0; i < s; i++) {
			int c = rnd.nextInt(endASCII-startASCII)+startASCII;
			sb.append(Character.toString((char)c));
		}
		
		return sb.toString();
	}
	
	/**
	 * 
	 * @param str has to be of this form [1,2,3][2,3,4]. The seperator has to be given.
	 * @param seperator
	 * @return
	 */
	public static int [][] getMatrixFromString(String str, String seperator){
		
		List<List<Integer>> lili = new ArrayList<List<Integer>>();
		
		int start = -1;
		
		int end = -1;
		
		String strTmp = str;
		
		int pos = 0;
		while(strTmp.length()>0) {
			start = str.indexOf("[", pos);
			
			if(start==-1){
				break;
			}
			
			end = StringFunctions.nextClosingBracket(str, start);
			
			pos = end;
			
			String sVector = str.substring(start+1, end);
			
			StringTokenizer st = new StringTokenizer(sVector, seperator);
			
			List<Integer> li = new ArrayList<Integer>();
			
			while(st.hasMoreTokens()){
				String sVal = st.nextToken();
				int v = Integer.parseInt(sVal);
				li.add(v);
			}
			lili.add(li);
		}
		
		int [] [] arr = new int [lili.size()][];
		
		for (int i=0; i < lili.size(); i++) {
			
			List<Integer> li = lili.get(i);
			
			int [] a = new int [li.size()];

			for (int j = 0; j < a.length; j++) {
				a[j] = li.get(j);
			}
			
			arr[i]=a;
		}
		
		return arr;
	}
	

	/**
	 * Finds the maximum common String in all Strings. Position independent.
	 * @param li
	 * @return
	 */
	public static String getMaximumOverlap(List<String> li, int minsize){
		String maxo = "";
		
		List<String> liPat = getAllSubStrings(li.get(0), minsize);
		for (int i = 0; i < liPat.size(); i++) {
			String pat = liPat.get(i);
			Pattern p = Pattern.compile(pat);
			boolean match=true;
			for (int j = 1; j < li.size(); j++) {
				String line = li.get(j);
				Matcher m = p.matcher(line);
				if(!m.find()){
					match=false;
					break;
				}
			}
			
			if(match){
				if(pat.length() > maxo.length()){
					maxo = pat;
				}
			}
		}
		
		return maxo;
	}

	
	public static String removeCharacter(String str, char c) {
		StringBuilder sb = new StringBuilder();
		
		for (int i = 0; i < str.length(); i++) {
			if(str.charAt(i) != c)
				sb.append(str.charAt(i));
		}
		
		return sb.toString();
	}
	
	public static String removeCharacter(StringBuilder str, char c) {
		StringBuilder sb = new StringBuilder();
		
		for (int i = 0; i < str.length(); i++) {
			if(str.charAt(i) != c)
				sb.append(str.charAt(i));
		}
		
		return sb.toString();
	}

	public static int countOccurence(String str, char c){
		
		int cc = 0;
	    for (int i=0; i < str.length(); i++) {
	        if (str.charAt(i) == c){
	             cc++;
	        }
	    }
	    return cc;
	}
	
	/**
	 * Not allowed are: \ / : * ? < > |
	 * @param str input string
	 * @return string with -X- instead of the not allowed characters.
	 * 10.09.2003 MK
	 */
	public static String convertToValidFileNameCharacters(String str) {
		String strFormated = "";
		for (int jj = 0; jj < str.length(); jj++) {
			char cChar = str.charAt(jj);
			String sInsert = "";

			if (cChar == '\\' || cChar == '/' || cChar == ':' || cChar == '*'
					|| cChar == '?' || cChar == '<' || cChar == '>') {
				sInsert = "-X-";

			} else {
				sInsert += cChar;
			}

			strFormated += sInsert;
		}
		return strFormated;
	}

	public static String toStringFileNameCompatible(double d) {
		String str = Double.toString(d);
		str = str.replace('.', '-');
		return str;
	}

	//
	// 05.08.2003 MvK
	/**
	 * This function was implemented because in AxoSOMSampleView was a new line
	 * in SMILES molConvert from ChemAxon that is not detected by replaceAll("\n", "");
	 * @param str input String
	 * @return a String only with printable ASCII characters. No extended ASCII
	 * characters.
	 */
	public static String formatToPrintableCharactersOnly(String str) {
		String strFormated = "";
		for (int i = 0; i < str.length(); i++) {
			int iChar = str.charAt(i);
			if (iChar > 31 && iChar < 127)
				strFormated += str.charAt(i);
		}
		return strFormated;
	}
	
	public static String formatToCharactersAndDigits(String str) {
		
		String regex = "[0-9a-zA-Z ]";
		Pattern pa = Pattern.compile(regex);
    	Matcher ma = pa.matcher(str);
    	StringBuilder sb = new StringBuilder();
    	int pos = 0;
    	while(ma.find(pos)){
    		MatchResult mr = ma.toMatchResult();
    		int start = mr.start();
    		int end = mr.end();
    		pos = end;
    		sb.append(str.substring(start,end));
    	}
		return sb.toString();
	}
	
	public static String format2DefinedLengthTrailing(String s, int length){
		StringBuilder sb = new StringBuilder(s);
		while(sb.length() < length)
			sb.append(' ');
		return sb.toString();
	}
	
	public static String format2DefinedLengthLeading(String s, int length){
		StringBuilder sb = new StringBuilder(s);
		while(sb.length() < length)
			sb.insert(0, ' ');
		return sb.toString();
	}
	
	/**
	 * Keeps the minus. Every other ASCII character, not a letter nor a digit is replaced with '_'. 
	 * @param str
	 * @return
	 */
	public static String format(String str) {
		return format(str, '_');
	}
	
	
	/**
	 * Keeps the minus. Every other ASCII character, not a letter nor a digit is replaced with <code>replacement</code>. 
	 * @param str
	 * @return
	 */
	public static String format(String str, char replacement) {
		StringBuilder sb = new StringBuilder();
		
		int cLastLetter = 0;
		for (int i = 0; i < str.length(); i++) {
			char ch = str.charAt(i);
			if (Character.isLetterOrDigit(ch)) {
				sb.append(ch);
				cLastLetter = ch;
			} else if (ch == '-') {
				sb.append(ch);
				cLastLetter = ch;
			} else if(cLastLetter != replacement) {
				sb.append(replacement);
				cLastLetter = replacement;
			}
		}
		return sb.toString();
	}
	
	/**
	 *
	 * @param sLine input string
	 * @param sStart start tag
	 * @param sEnd end tag
	 * @return string between the two tags, if one the tags is not founds a string
	 * with the length 0 is returned.
	 */
	public static String getString(String sLine, String sStart, String sEnd) {
		String str = "";
		int iStart = sLine.indexOf(sStart) + sStart.length();
		int iEnd = sLine.indexOf(sEnd);
		if (iStart > -1 && iEnd > -1)
			str = sLine.substring(iStart, iEnd);
		return str;
	}

	/**
	 * 
	 * @param str
	 * @param regex
	 * @return null if substring not found.
	 */
	public static String getStringFromRegEx(String str, String regex) {

		Pattern pa = Pattern.compile(regex);
    	Matcher ma = pa.matcher(str);
		int start = -1;
		int end = -1;
    	if(ma.find()){
    		MatchResult mr = ma.toMatchResult();
    		start = mr.start();
    		end = mr.end();
    	} else {
    		return null;
    	}

		return str.substring(start,end);
	}

	public static Point getStartEnd(String str, String regex) {
		Pattern pa = Pattern.compile(regex);
		Matcher ma = pa.matcher(str);
		int start = -1;
		int end = -1;
		if(ma.find()){
			MatchResult mr = ma.toMatchResult();
			start = mr.start();
			end = mr.end();
		} else {
			return null;
		}
		return new Point(start,end);
	}

	public static List<String> splitIncludeSEP(String str, String regex) {
		List<String> li=new ArrayList<>();
		Pattern pa = Pattern.compile(regex);
		Matcher ma = pa.matcher(str);
		int start = 0;
		int end = -1;
		while(ma.find()){
			MatchResult mr = ma.toMatchResult();
			end = mr.start();
			String s = str.substring(start, end);
			start = end;
			li.add(s);
		}
		return li;
	}

	
	public static boolean isRegexInString(String str, String regex) {
		Pattern pa = Pattern.compile(regex);
    	Matcher ma = pa.matcher(str);
    	if(ma.find()){
    		return true;
    	} else {
    		return false;
    	}
	}
	
	/**
	 * 
	 * @param str
	 * @param regex
	 * @return expression which was matched by regex.
	 */
	public static String extract(String str, String regex) {
		String substring = "";
		Pattern pa = Pattern.compile(regex);
    	Matcher ma = pa.matcher(str);
    	if(ma.find()) {
			MatchResult mr = ma.toMatchResult();
			substring = mr.group();
    	}
    	return substring;
	}
	
	/**
	 * 
	 * @param str
	 * @param regex
	 * @return the combined not matching parts of the string.
	 */
	public static String extractInverse(String str, String regex) {

		String substring = "";
		Pattern pa = Pattern.compile(regex);
    	Matcher ma = pa.matcher(str);
    	if(ma.find()) {
			MatchResult mr = ma.toMatchResult();
			int start = mr.start();
			int end = mr.end();
			String rest1 = str.substring(0, start);
			String rest2 = str.substring(end);
			substring = rest1 + rest2;
		}
    	return substring;
	}

	/**
	 *
	 * @param sLine input string
	 * @param sStart start tag
	 * @param sEnd end tag
	 * @param iFromIndex start index
	 * @return string between the two tags, if one the tags is not founds a string
	 * with the length 0 is returned.
	 * @return
	 */
	public static String getString(String sLine, String sStart, String sEnd,
			int iFromIndex) {
		String str = "";

		int iStart = sLine.indexOf(sStart, iFromIndex);
		int iEnd = sLine.indexOf(sEnd, iStart);

		if (iStart > -1 && iEnd > -1) {
			iStart += sStart.length();
			str = sLine.substring(iStart, iEnd);
		}

		return str;
	}

	/**
	 * Removes all non characters and digits.
	 * @param txt
	 * @return
	 */
	public static List<String> getWordsFormatted(String txt) {
		String txtFormatted = StringFunctions.formatToCharactersAndDigits(txt);
		
		List<String> li = new ArrayList<String>();
		
		StringTokenizer st = new StringTokenizer(txtFormatted);
		
		while(st.hasMoreTokens()) {
			String s = st.nextToken().trim();
			li.add(s);
		}
		
		return li;
	}
	
	/**
	 * 
	 * @param txt
	 * @return formatted, unique and lower case
	 */
	public static List<String> getWordsFormattedUniqueLowerCase(String txt) {
		String txtFormatted = StringFunctions.formatToCharactersAndDigits(txt);
		
		HashSet<String> hs = new HashSet<String>();
		
		StringTokenizer st = new StringTokenizer(txtFormatted);
		
		while(st.hasMoreTokens()) {
			String s = st.nextToken().trim().toLowerCase();
			hs.add(s);
		}
		
		return new ArrayList<String>(hs);
	}
	
	/**
	 * Get a list from quoted and comma or otherwise separated phrases.
	 * @param txt
	 * @return
	 */
	public static List<String> getTokenizedQuoted(String txt) {
		
		List<String> li = new ArrayList<String>();
		
		int startPos = 0;
		
		while (startPos < txt.length()) {
			int indexBegin = txt.indexOf('"', startPos);
			
			if(indexBegin==-1){
				break;
			}
			
			int indexEnd = txt.indexOf('"', indexBegin+1);
			
			if(indexEnd==-1){
				break;
			}
			
			String sub = txt.substring(indexBegin+1, indexEnd);
			
			li.add(sub);
			
			startPos = indexEnd+1;
		}
		
		if(startPos==0){
			li.add(txt);
		}
		
		return li;
	}
	
	
	/**
	 * Returns the tokenized and trimmed values.
	 * @param txt
	 * @param separator
	 * @return
	 */
	public static List<String> getTokenized(String txt, String separator) {
		
		List<String> li = new ArrayList<String>();
		
		StringTokenizer st = new StringTokenizer(txt, separator);
		
		while(st.hasMoreTokens()){
			String tok = st.nextToken().trim();
			
			li.add(tok);
		}
		
		
		return li;

	}
	
	public static List<String> getTokenizedBySeperatorRegex(String txt, String regex) {
				
		List<Point> liSep = match(txt, regex);
		
		List<String> li = new ArrayList<String>();
		
		if(liSep.size()==0){
			li.add(txt);
			return li;
		}
		
		int start = 0;
		
		for (int i = 0; i < liSep.size(); i++) {
			Point p = liSep.get(i);
			
			String tok = txt.substring(start, p.x);
			
			li.add(tok.trim());
			
			start = p.y;
			
		}
				
		Point p = liSep.get(liSep.size() - 1);
		
		start = p.y;
		
		String tok = txt.substring(start);
		
		li.add(tok.trim());
				
		return li;

	}

	/**
	 * Generates a list with overlap
	 * @param liWords
	 * @param lenSubText so many words are in each entry.
	 * @param lenOverlap 
	 * @return
	 */
	public static List<String> getSplittedOverlappingText(List<String> liWords, int lenSubText, int lenOverlap) {
		List<String> li = new ArrayList<String>();

		int steplen = lenSubText-lenOverlap;
		for (int i = 0; i < liWords.size(); i+=steplen) {
			
			int end = Math.min(i+lenSubText, liWords.size());
			
			StringBuilder sb = new StringBuilder();
			for (int j = i; j < end; j++) {
				sb.append(liWords.get(j));
				if(j<end-1) {
					sb.append(" ");
				}
			}
			
			li.add(sb.toString());
		}
		
		return li;
	}

	/**
	 * https://stackoverflow.com/questions/4385623/bytes-of-a-string-in-java
	 * sizeof(string) =
	 * 8 + // object header used by the VM
	 * 8 + // 64-bit reference to char array (value)
	 * 8 + string.length() * 2 + // character array itself (object header + 16-bit chars)
	 * 4 + // offset integer
	 * 4 + // count integer
	 * 4 + // cached hash code
	 * @param s
	 * @return
	 */
	public static int sizeOf(String s) {
		return 36 + s.length() * 2;
	}

	public static int sizeOf(List<String> l) {
		int n = 0;

		for (String s : l) {
			n+=sizeOf(s);
		}

		return n;
	}

	
	public static String toString(double [] arr, NumberFormat nf){
		return toString(arr, nf, 0);
	}

	public static String toString(double [] arr, NumberFormat nf, int width){

		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {


			String str = nf.format(arr[i]);

			while (str.length()<width){
				str = " " + str;
			}

			sb.append(str);
			if(i < arr.length-1){
				sb.append(ConstantsDWAR.SEP_VALUE);
			}
		}
		return sb.toString();
	}

	public static String toString(float [] arr, NumberFormat nf){

		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			sb.append(nf.format(arr[i]));
			if(i < arr.length-1){
				sb.append(ConstantsDWAR.SEP_VALUE);
			}
		}
		return sb.toString();
	}

	public static String toString(byte [] arr){

		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			int v = arr[i];
			sb.append(v);
			if(i < arr.length-1){
				sb.append(ConstantsDWAR.SEP_VALUE);
			}
		}
		return sb.toString();
	}
	public static String toStringShort(byte [] arr){

		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			int v = arr[i];
			sb.append(v);
			if(i < arr.length-1){
				sb.append(" ");
			}
		}
		return sb.toString();
	}

	public static String toString(boolean [] arr){

		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			if(arr[i]) {
				sb.append(1);
			} else {
				sb.append(0);
			}
		}

		return sb.toString();
	}
	public static String toString(boolean [] [] arr){

		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			for (int j = 0; j < arr[i].length; j++) {
				if (arr[i][j]) {
					sb.append(1);
				} else {
					sb.append(0);
				}
			}
			sb.append("\n");
		}

		return sb.toString();
	}

    public static String toString(List<Double> li, NumberFormat nf){

        StringBuilder sb = new  StringBuilder();
        for (int i = 0; i < li.size(); i++) {
            sb.append(nf.format(li.get(i)));
            if(i < li.size()-1){
                sb.append(ConstantsDWAR.SEP_VALUE);
            }
        }
        return sb.toString();
    }

    public static String toStringInteger(List<Integer> li, String sep){

        StringBuilder sb = new  StringBuilder();
        for (int i = 0; i < li.size(); i++) {
            sb.append(li.get(i));
            if(i < li.size()-1){
                sb.append(sep);
            }
        }
        return sb.toString();
    }

//    public static String toStringInt(List<Integer> li){
//
//        StringBuilder sb = new  StringBuilder();
//        for (int i = 0; i < li.size(); i++) {
//            sb.append(li.get(i));
//            if(i < li.size()-1){
//                sb.append(ConstantsDWAR.SEP_VALUE);
//            }
//        }
//        return sb.toString();
//    }

	public static String toString(double [] arr){

        NumberFormat nf = new DecimalFormat("0.000");

		return toString(arr, nf);
	}


	public static String toString(int [][] arr, String seperator){
		
		StringBuilder sb = new  StringBuilder();
		
		for (int i = 0; i < arr.length; i++) {
			
			sb.append("[");
			
			for (int j = 0; j < arr[i].length; j++) {
				sb.append(arr[i][j]);
				
				if(j<arr[i].length-1){
					sb.append(seperator);
				}
			}
			
			sb.append("]");
			
		}
		
		return sb.toString();
	}

	public static String toString(String [] arr, String seperator){

		StringBuilder sb = new  StringBuilder();

		for (int i = 0; i < arr.length; i++) {

			if(sb.length()>0)
				sb.append(seperator);

			sb.append(arr[i]);

		}

		return sb.toString();
	}

	public static String toString(int [] arr, String seperator){
		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			if(sb.length()>0)
				sb.append(seperator);
			sb.append(Integer.toString(arr[i]));
		}
		return sb.toString();
	}

	/**
	 * Elements are separated by tabs and rows are separated by newline.
	 * @param arr
	 * @return
	 */
	public static String toStringTabNL(String [][] arr){
		StringBuilder sb = new  StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			for (int j = 0; j < arr[i].length; j++) {
				sb.append(arr[i][j]);
				if(j <arr [i].length-1){
					sb.append("\t");
				}
			}
			if(i < arr.length-1){
				sb.append("\n");
			}
		}
		return sb.toString();
	}


	public static String toString(Exception ex) {

		StringWriter sw = new StringWriter();
		ex.printStackTrace(new PrintWriter(sw));
		String exceptionAsString = sw.toString();

		return exceptionAsString;
	}


	public static String toString(List<String> li) {
		return toString(li, " ");
	}

	public static String toString(Collection<String> li, String sep) {
		StringBuilder sb = new StringBuilder();
		boolean started=false;
		for (String s : li) {
			if(started){
				sb.append(sep);
			}
			sb.append(s);
			started=true;
		}

		return sb.toString();
	}
	
	public static String toStringLong(List<Long> li, String sep) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < li.size(); i++) {
			sb.append(li.get(i));
			if(i < li.size()-1){
				sb.append(sep);
			}
		}
		return sb.toString();
	}
	public static String toStringInt(List<Integer> li, String sep) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < li.size(); i++) {
			sb.append(li.get(i));
			if(i < li.size()-1){
				sb.append(sep);
			}
		}
		return sb.toString();
	}
	public static String toStringInt(List<Integer> li) {
		return toStringInt(li, SEP);
	}

	public static String toSortedString(List<String> li) {
		StringBuilder sb = new StringBuilder();
		Collections.sort(li);
		for (String s : li) {
			sb.append(s);
		}
		return sb.toString();
	}

	public static String toStringBinary(int v) {
		String str = "";

		int len = Integer.SIZE;

		for (int ii = 0; ii < len; ii++) {
			if ((v & 1) == 1) {
				str = "1 " + str;
			} else {
				str = "0 " + str;
			}
			v = v >> 1;
		}

		return str.trim();
	}

	public static String toStringBinary(long v) {
		String str = "";

		int len = Long.SIZE;

		for (int ii = 0; ii < len; ii++) {
			if ((v & 1) == 1) {
				str = "1 " + str;
			} else {
				str = "0 " + str;
			}
			v = v >> 1;
		}

		return str.trim();
	}

	public static String toStringHex(String s) {
		String sHex = "";
		for (int i = 0; i < s.length(); i++) {
			String h = Integer.toHexString(s.charAt(i));
			while (h.length() < 2) {
				h = '0' + h;
			}
			sHex += h;
		}
		return sHex;
	}

	public static String hex2String(String hex) {
		String s = "";
		for (int i = 0; i < hex.length(); i += 2) {
			String h = hex.substring(i, i + 2);
			char c = (char) Integer.parseInt(h, 16);
			String sC = Character.toString(c);
			s += sC;
		}
		return s;
	}

	/**
	 * finds the next balanced closing bracket "]" to the first open bracket
	 * "[" in the string.
	 * @param txt String
	 * @param iIndexStart start index
	 * @return index of the next corresponding bracket
	 */
	public static int nextClosingBracket(String txt, int iIndexStart) {
		// sLine = "[[(Com)][(.*I)]][ny)it)]";
		return nextClosing(txt, iIndexStart, '[', ']');
	}

	/**
	 * Escapes the meta characters in a regular expression pattern with \\. 
	 * @param pattern
	 * @return
	 */
	public static String escapeDanglingMetaCharacters(String pattern){
		
		for (int i = 0; i < StringFunctions.REGEX_META_CHARACTERS.length; i++) {
			pattern = pattern.replace(StringFunctions.REGEX_META_CHARACTERS[i], "\\"+StringFunctions.REGEX_META_CHARACTERS[i]);
		}
		
		return pattern;
	}
	
	/**
	 * 
	 * @param str
	 * @param regex
	 * @return list with points, x start, y end of matching string (offset after the last character matched).
	 */
    public static final List<Point> match(String str, String regex) {
    	
    	Pattern pa = Pattern.compile(regex);
    	
    	Matcher ma = pa.matcher(str);
    	
    	List<Point> li = new ArrayList<Point>();
    	
    	while(ma.find()){
    		MatchResult mr = ma.toMatchResult();
    		int start = mr.start();
    		int end = mr.end();
    		
    		li.add(new Point(start, end));
    		
    	}
    	
    	return li;
    }
    
    public static final Point matchFirst(String str, String regex) {
    	
    	Pattern pa = Pattern.compile(regex);
    	
    	Matcher ma = pa.matcher(str);
    
    	Point p = null;
    	
    	if(ma.find()){
    		MatchResult mr = ma.toMatchResult();
    		int start = mr.start();
    		int end = mr.end();
    		
    		p = new Point(start, end);
    		
    	}
    	
    	return p;
    }

	
	/**
	 * finds the next corresponding closing bracket char to the first open char
	 * @param txt
	 * @param indexStart
	 * @param cOpen
	 * @param cClose
	 * @return
	 */	
	public static int nextClosing(String txt, int indexStart, char cOpen, char cClose) {
		// sLine = "[[(Com)][(.*I)]][ny)it)]";

		int indexBracket = -1;
		
		int indexOpenBracket = txt.indexOf(cOpen, indexStart);
		
		int numberOpenBrackets = 1;
		
		for (int i = indexOpenBracket + 1; i < txt.length(); i++) {
			
			if (txt.charAt(i) == cOpen) {
				
				numberOpenBrackets++;
				
			} else if (txt.charAt(i) ==  cClose) {
				
				numberOpenBrackets--;
			}
			
			if (numberOpenBrackets == 0) {
				
				indexBracket = i;
				break;
			}
		}
		
		return indexBracket;
	}
	
	public static boolean isAllLetter(String s) {
		
	    for(char c : s.toCharArray()) {
	    	
	       if(!Character.isLetter(c)) {
	           return false;
	        }
	    }
	    
	    return true;
	}
	
	public static boolean isAllUpperCase(String s) {
		
	    for(char c : s.toCharArray()) {
	    	
	       if(Character.isLetter(c) && Character.isLowerCase(c)) {
	           return false;
	        }
	    }
	    
	    return true;
	}

	public static boolean isAllLowerCase(String s) {

	    for(char c : s.toCharArray()) {

	       if(Character.isLetter(c) && Character.isUpperCase(c)) {
	           return false;
	        }
	    }

	    return true;
	}

	public static boolean containsUpperCase(String s) {

		for(char c : s.toCharArray()) {

			if(Character.isLetter(c) && Character.isUpperCase(c)) {
				return true;
			}
		}

		return false;
	}

	public static boolean containsLowerCase(String s) {

		for(char c : s.toCharArray()) {

			if(Character.isLetter(c) && Character.isLowerCase(c)) {
				return true;
			}
		}

		return false;
	}


	public static boolean isAlphaNumeric(char char1) {
		return (char1 >= 'a' && char1 <= 'z') || (char1 >= 'A' && char1 <= 'Z') || (char1 >= '0' && char1 <= '9');
	}

	public static boolean isAlphaNumeric(String s) {
		boolean a = true;
		for (int i = 0; i < s.length(); i++) {
			Character c = s.charAt(i);
			if(!isAlphaNumeric(c)){
				a = false;
				break;
			}
		}
		return a;
	}

	/**
	 * 
	 * @param s
	 * @return true only if the first letter is capitalized and all other words are lower case letters.
	 */
	public static boolean isCapitalizedWord(String s) {
		
		if(s.length() < 1){
			return false;
		}
		
		if(!Character.isUpperCase(s.charAt(0))){
			return false;
		}
		
		// Only one letter
		if(s.length() < 2){
			return true;
		}
		
		String sub = s.substring(1);
		
	    for(char c : sub.toCharArray()) {
	    	
	       if(!Character.isLetter(c) || Character.isUpperCase(c)) {
	           return false;
	        }
	    }
	    
	    return true;
	}
	
	public static boolean isUpperAndLowerCase(String s) {
		
		boolean upper=false;
		boolean lower=false;
		
	    for(char c : s.toCharArray()) {
	       if(Character.isUpperCase(c)) {
	    	   upper = true;
	        } else if(Character.isLowerCase(c)) {
	    	   lower = true;
	        }
	       if(upper && lower){
	    	   break;
	       }
	    }
	    
	    if(upper && lower){
	    	return true;
	    }
	    return false;
	}
	
	
	
	/**
	 * 
	 * @param name
	 * @return false if for each opening parenthesis none closing one is present.
	 */
	public static boolean isMissingParenthesis(String name) {
		
		int ccOpen=0;
		int ccClose=0;
		for (int i = 0; i < name.length(); i++) {
			if(name.charAt(i) == '(') {
				ccOpen++;
			} else if(name.charAt(i) == ')') {
				ccClose++;
			}
		}
		
		if(ccOpen==ccClose) {
			return false;
		}
		
		return true;
	}

	public static String toStringStackTrace(Exception ex){
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		ex.printStackTrace(pw);
		return sw.toString();

	}
	public static void main(String[] args) {

		String sLine = "Pos3ition: 8 15 StartName:XXXEn7890dName0";

		System.out.println(countWordInText(sLine));

		/*
		 String sLine = "Position: 8 15 StartName:XXXEndName";

		 sLine.replaceAll("\\s", "");

		 System.out.println(sLine);

		 String sStart = "StartName:";
		 String sEnd = "EndName";

		 String str = getString(sLine, sStart, sEnd);
		 System.out.println(str);
		 */
	}


}
