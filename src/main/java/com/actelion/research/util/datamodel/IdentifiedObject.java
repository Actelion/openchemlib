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

package com.actelion.research.util.datamodel;

import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;


/**
 * 
 * 
 * IdentifiedObject
 * Nov 2, 2011 MvK: Start implementation
 */
public class IdentifiedObject<T> implements IIdentifiedObject<T>, Comparable<IdentifiedObject<T>> {
	
	
	private T data;
	
	private long id;
	
	public IdentifiedObject() {
		
	}
	
	public IdentifiedObject(T data, long id) {
		
		this.data = data;
		
		this.id = id;
		
	}

	
	
	public T getData() {
		return data;
	}

	public void setData(T data) {
		this.data = data;
	}

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
				
		sb.append(id);
				
		return sb.toString();
	}

	
	@Override
	public int compareTo(IdentifiedObject<T> io) {
		
		if(id > io.id){
			return 1;
		} else if(id < io.id){
			return -1;
		}
		
		return 0;
	}

	public static<T> HashMap<Long, T> getHashMap(List<IdentifiedObject<T>> li){
		
		HashMap<Long, T> hmId_Descriptor = new HashMap<Long, T>(li.size());

		for (IdentifiedObject<T> idObj : li) {
			hmId_Descriptor.put(idObj.getId(), idObj.getData());
		}
		
		return hmId_Descriptor;
	}

	// Implementing Fisher–Yates shuffle
	public static void shuffleArray(IdentifiedObject[] ar) {
		Random rnd = ThreadLocalRandom.current();
		for (int i = ar.length - 1; i > 0; i--) {
			int index = rnd.nextInt(i + 1);
			// Simple swap
			IdentifiedObject a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}
	}

	public static Comparator<IdentifiedObject> getComparatorId() {

		return new Comparator<IdentifiedObject>() {

			public int compare(IdentifiedObject o1, IdentifiedObject o2) {

				if(o1.id > o2.id) {
					return 1;
				} else if(o2.id > o1.id) {
					return -1;
				}

				return o2.compareTo(o1);
			}
		};
	}


}
