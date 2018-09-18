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

package com.actelion.research.io;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.channels.ReadableByteChannel;
import java.util.Arrays;

import com.actelion.research.util.ConstantsDWAR;
import com.actelion.research.util.Pipeline;

/**
 *
 * StringReadChannel
 * 2007 MvK: Start implementation
 * 25.06.2009 MvK: implementation changed
 * 12.02.2014 MvK: added charset encoding to handle Umlaute.
 * 24.04.2014 MvK: Pipeline replaced simple LinkedList because of needed concurrent access.
 * 29.01.2015 MvK: Increased capacity CAPACITY_LINE_BUFFER to 10,000,000 because of overflow when reading PubMed records.
 * 03.06.2015 MvK: Increased capacity CAPACITY_LINE_BUFFER to 50,000,000 because of overflow when reading g2dDiseasePublicationSlope.dwar
 */
public class StringReadChannel {

	private static final int CAPACITY_LINE_BUFFER = 50000000;
	
	private static final int CAPACITY_READ_BUFFER = 100000;
	
	private static final int CAPACITY_PIPE = 1000;
	
	
	
	private ReadableByteChannel byteChannel;
		
	private ByteBuffer buffer;
	
	private ByteBuffer byteBufferLine;
	
	private Pipeline<String> pipeline;
		
	public StringReadChannel(ReadableByteChannel ch)throws IOException{
		init(ch);
	}
	
	private void init(ReadableByteChannel bc) throws IOException{
		
		byteChannel = bc;
		

		// We will fill up this ByteBuffer instance later.
		byteBufferLine = ByteBuffer.allocate(CAPACITY_LINE_BUFFER);
		
		buffer = ByteBuffer.allocate(CAPACITY_READ_BUFFER);
		
		pipeline = new Pipeline<String>();
		
		readLine2List();
	}
	
	public boolean hasMoreLines() throws IOException {
		
		return !pipeline.wereAllDataFetched();
		
	}
	
	/**
	 * 
	 * @return null if EOF reached.
	 * @throws IOException
	 */
	public String readLine() throws IOException {

		int maxCycles = 100;
		
		String str = null;
		
		
		if(!pipeline.isEmpty()){
			
			str = pipeline.pollData();
			
		} else if(!pipeline.isAllDataIn()){
			
			int ccCycle=0;
			
			while(pipeline.isEmpty()){
			
				try {Thread.sleep(100);} catch (InterruptedException e) {e.printStackTrace();}
			
				ccCycle++;
				
				if(ccCycle>maxCycles){
					
					RuntimeException ex = new RuntimeException("Arbitrary break. Max number of cycles (" + maxCycles + ") exceeded.");
					
					// ex.printStackTrace();
					
					byteChannel.close();
					
					throw ex;
				}
				
				str = pipeline.pollData();
			}
		}
		
		
		if(!pipeline.isAllDataIn()){
			if(pipeline.sizePipe() < CAPACITY_PIPE){
				readLine2List();	
			}
			
			
		}
		
		
		return str;
	}
	
	public void finalize()throws IOException{
		close();
	}
	
	private int readLine2List() throws IOException {
				
		boolean lineFinished=false;
		
		int sizeLine = 0;
		
		while(!lineFinished) {
			
			((Buffer)buffer).clear();
			
			int size = byteChannel.read(buffer);
					
			if(size==-1){ // end of stream
				
				if(byteBufferLine.position() > 0){
				
					writeBuffer2Pipe();
					lineFinished = true;
					
				}
				
				pipeline.setAllDataIn(true);
							
				return -1;
			}
			
			//
			// Copying from one buffer to the other
			//
			for (int i = 0; i < size; i++) {
				
				byte c = buffer.get(i);
							
				if(c=='\n') {
										
					writeBuffer2Pipe();
						
					lineFinished = true;
					
				} else {
					if(c!='\r') {
					
						byteBufferLine.put(c);
						
						sizeLine++;
						
					}
				}
			}
			
		}
		
		return sizeLine;
	}
	
	private void writeBuffer2Pipe() throws UnsupportedEncodingException{
		
		byte [] contentsOnly = Arrays.copyOf(byteBufferLine.array(), byteBufferLine.position());
		
		String str = new String(contentsOnly, ConstantsDWAR.CHARSET_ENCODING);
		
		StringBuilder sb = new StringBuilder(str);
		
		pipeline.addData(sb.toString());
		
		((Buffer)byteBufferLine).clear();
	}
	
	
	
//	private int readLine2List() throws IOException {
//		
//		ByteBuffer buffer = ByteBuffer.allocate(CAPACITY);
//		
//		int size = byteChannel.read(buffer);
//				
//		if(VERBOSE){
//			System.out.println("StringReadChannel readLine2List() size " + size);
//		}
//		
//		
//		if(size==-1){
//			
//			if(byteBufferLine.position() > 0){
//			
//				byte [] contentsOnly = Arrays.copyOf(byteBufferLine.array(), byteBufferLine.position());
//							
//				String str = new String(contentsOnly, ConstantsDWAR.CHARSET_ENCODING);
//				
//				StringBuilder sb = new StringBuilder(str);
//				
//				liLine.add(sb);
//			
//			}
//			
//			bEOF = true;
//			return -1;
//		}
//		
//		for (int i = 0; i < size; i++) {
//			byte c = buffer.get(i);
//			
//			if(c==-1) {
//				bEOF=true;
//			} else {
//				if(c=='\n') {
//										
//					byte [] contentsOnly = Arrays.copyOf(byteBufferLine.array(), byteBufferLine.position());
//					
//					byteBufferLine.clear();
//					
//					String str = new String(contentsOnly, ConstantsDWAR.CHARSET_ENCODING);
//					
//					StringBuilder sb = new StringBuilder(str);
//					
//					liLine.add(sb);
//										
//				} else {
//					if(c!='\r') {
//					
//						byteBufferLine.put(c);
//											
//					}
//				}
//			}
//		}
//		
//		return size;
//	}
	
    public static void skipUntilLineMatchesRegEx(StringReadChannel src, String regex) throws NoSuchFieldException, IOException {
    	int limit = 10000;
    	
    	skipUntilLineMatchesRegEx(src, regex, limit);
    }
    
    public static String skipUntilLineMatchesRegEx(StringReadChannel src, String regex, int limit) throws NoSuchFieldException, IOException {
    	    	    		
    	String line = src.readLine();
    	
    	boolean match = false;
    	
    	if(line.matches(regex)){
    		match = true;
    	}
    	
    	int cc=0;
    	
		while(!match){
    		
			if(!src.hasMoreLines()){
				break;
			}
			
			line = src.readLine();
			
	    	if(line.matches(regex)){
	    		match = true;
	    	} 
			
	    	cc++;
	    	
	    	if(cc > limit){
	    		break;
	    	}
	    	
    		
    	}	
		
		String lineMatch = null;
		
		if(match){
			
			lineMatch = line;
			
		} else {
			throw new NoSuchFieldException("Regex " + regex + " was not found.");
		}
    	
		return lineMatch;
		
    	
    }

    public void close() throws IOException {
    	byteChannel.close();
    }

}
