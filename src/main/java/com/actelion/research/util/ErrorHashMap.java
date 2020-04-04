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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

/**
 * 
 * ErrorHashMap
 * 2005 MvK: Start implementation

 */
public class ErrorHashMap implements Serializable {

	private static final long serialVersionUID = 20122011;
	
	private HashMap<String,ExceptFreq> hmMessage_ExceptionFreq;
	
	public ErrorHashMap(){
		hmMessage_ExceptionFreq = new HashMap<String, ExceptFreq>();
		
	}
	
	public void add(ErrorHashMap ehm){
		List<ExceptFreq> liExceptFreq = new ArrayList<ExceptFreq>(ehm.hmMessage_ExceptionFreq.values());
		for (ExceptFreq exceptFreq : liExceptFreq) {
			add(exceptFreq.get());
		}
	}

	public void addAll(List<Exception> li){
		for (Exception exception : li) {
			add(exception);
		}
	}

	public void add(Exception ex){
		
		StackTraceElement [] arrSTr = ex.getStackTrace();
		
		if(arrSTr.length==0){
			String key = ex.getMessage();
			ExceptFreq ef = hmMessage_ExceptionFreq.get(key);
			if(ef==null){
				ef = new ExceptFreq(ex);
				hmMessage_ExceptionFreq.put(key, ef);
			} else {
				ef.increase();
			}
			
		} else {
			int ln = arrSTr[0].getLineNumber();
			String method = arrSTr[0].getMethodName();
			
			String key = method + ln; 
			
			ExceptFreq ef = hmMessage_ExceptionFreq.get(key);
			
			if(ef==null){
				ef = new ExceptFreq(ex);
				hmMessage_ExceptionFreq.put(key, ef);
			} else {
				ef.increase();
			}
		}
	}
	
	public void clear(){
		hmMessage_ExceptionFreq.clear();
	}
	
	public boolean hasErrors(){
		if(hmMessage_ExceptionFreq.size()>0)
			return true;
		else return false;
	}
	
	public static void main(String [] args){
		
		ErrorHashMap errorHashMap = new ErrorHashMap();
		
		for (int i = 0; i < 10; i++) {
			
			try {
				Integer.parseInt("abc");
			} catch (NumberFormatException e) {
				errorHashMap.add(e);
			}
		}
		
		System.err.println(errorHashMap.toString());
		
	}
	
	public int getSize(){
		return hmMessage_ExceptionFreq.size();
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		List<ExceptFreq> li = new ArrayList<ExceptFreq>(hmMessage_ExceptionFreq.values());
		for (ExceptFreq exfreq : li) {
			sb.append(exfreq.toString());
			sb.append("\n");
		}
		return sb.toString();
	}
	public String toStringShortMessage(){
		StringBuilder sb = new StringBuilder();
		List<ExceptFreq> li = new ArrayList<ExceptFreq>(hmMessage_ExceptionFreq.values());
		for (ExceptFreq exfreq : li) {
			sb.append(exfreq.toStringShortMessage());
			sb.append("\n");
		}
		return sb.toString();
	}
	
	
	
}
class ExceptFreq implements Serializable {
	
	private static final long serialVersionUID = 20122011;

	private Exception exception;
	
	private int frequency;
	
	private List<Date> liDate;
	
	public ExceptFreq(Exception ex){
		exception = ex;
		liDate = new ArrayList<Date>();
		increase();
	}
	
	public void increase(){
		frequency++;
		liDate.add(new Date());
	}
	
	public Exception get(){
		return exception;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		
		sb.append("The exception:\n");
		sb.append(exception.getMessage()+"\n");
		StackTraceElement [] arrSTr = exception.getStackTrace();
		for (int i = 0; i < arrSTr.length; i++) {
			sb.append(arrSTr[i].toString());
			sb.append("\n");
		}
		
		sb.append("\nWas called " + frequency + " times.");

		return sb.toString();
	}
	
	public String toStringShortMessage(){
		StringBuilder sb = new StringBuilder();
		sb.append("The exception:\n");
		sb.append(exception.getMessage()+"\n");
		StackTraceElement [] arrSTr = exception.getStackTrace();
		sb.append(arrSTr[0].toString());
		sb.append("\n");
		return sb.toString();
	}
	
	public String toStringDates(){
		StringBuilder sb = new StringBuilder();
		
		sb.append("The exception:\n");
		sb.append(exception.getMessage()+"\n");
		StackTraceElement [] arrSTr = exception.getStackTrace();
		for (int i = 0; i < arrSTr.length; i++) {
			sb.append(arrSTr[i].toString());
			sb.append("\n");
		}
		
		sb.append("\nWas called " + frequency + " times: ");
		if(liDate.size()<1000){
			for (Date d : liDate) {
				sb.append(d);	
				sb.append(", ");	
			}
		}

		sb.deleteCharAt(sb.length()-1);
		sb.deleteCharAt(sb.length()-1);
		
		return sb.toString();
	}
	
	
}