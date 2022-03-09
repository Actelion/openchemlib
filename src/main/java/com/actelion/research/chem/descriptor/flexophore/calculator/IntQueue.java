/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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
 */

package com.actelion.research.chem.descriptor.flexophore.calculator;

import com.actelion.research.util.ArrayUtils;

/**
 * Implementation for a FIFO of int. This class is particularly useful
 * for the implementation of BFS algorithm 
 * @author freyssj
 */
public final class IntQueue {
	
	private int[] s;
	private int begin = 0;
	private int end = 0;
	
	/**
	 * Creeates a new Queue with an initial size of 100
	 *
	 */
	public IntQueue() {
		 this(100);
	}
	
	public IntQueue(int size) {
		 s = new int[size];
	}
	
	public final int get(int ind) {
		return s[ind];
	}

	public final int getBegin() {
		return begin;
	}

	public final int getEnd() {
		return end;
	}

	public final boolean isEmpty() {
		return end <= begin;
	}
	
	public final void clear() {
		begin = end = 0;
	}
	
	public final int getSize() {
		return end - begin;
	}
	
	/**
	 * Pop the first element of the Queue
	 * @return
	 */
	public final int pop() {
		return s[begin++];
	}

	/**
	 * Peek the first element of the Queue
	 * (do not increment the queue index)
	 * @return
	 */
	public final int peek() {
		return s[begin++];
	}
	
	/**
	 * Push an element on top of the queue
	 * @param i
	 */
	public void push(int i) {
		if(end>=s.length) {
			ArrayUtils.resize(s, s.length*2);
			s = (int[]) ArrayUtils.resize(s, s.length*2);
		} 
		s[end++] = i;
	}
	
	public int indexOf(int ind) {
		for(int i=0; i<end; i++) {
			if(s[i]==ind) return i;
		}
		return -1;
	}
		
}
