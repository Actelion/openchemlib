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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;

/**
 * 
 * 
 * Pipeline
 * Enables concurrent access to a queue.
 * @author Modest von Korff
 * Mar 27, 2012 MvK: Start implementation
 * Oct 9 2012 MvK: bug fix, added reset()
 */
public class Pipeline<T> implements IPipeline<T>, Supplier<T> {

	private AtomicBoolean allDataIn;
	
	private ConcurrentLinkedQueue<T> queue;
	
	private AtomicLong added;
	
	private AtomicLong polled;
	
	public Pipeline() {
		allDataIn = new AtomicBoolean(false);
		queue = new ConcurrentLinkedQueue<T>();
		added = new AtomicLong();
		polled = new AtomicLong();
	}
	
	/**
	 * The 'all data in' flag is set true.
	 * @param li
	 */
	public Pipeline(List<T> li) {
		this();
		queue.addAll(li);
		setAllDataIn(true);
	}

	/**
	 * Sets all to 0 and allDataIn to false.
	 */
	public void reset(){
		allDataIn.set(false);
		added.set(0);
		polled.set(0);
		queue.clear();
	}
	
	public boolean isAllDataIn() {
		return allDataIn.get();
	}

	/**
	 * has to be set true or <code>wereAllDataFetched()</code> will never become true. 
	 */
	public void setAllDataIn(boolean allDataIn) {
		this.allDataIn.set(allDataIn);
	}

	public void setAllDataIn() {
		this.allDataIn.set(true);
	}

	public void addData(T t) {
		queue.add(t);
		added.incrementAndGet();
	}
	
	public void addData(List<T> li) {
		queue.addAll(li);
		added.addAndGet(li.size());
	}

	/**
	 * 
	 * @return null if nothing is in the queue.
	 */
	public T pollData() {
		T t = queue.poll();
		if(t!=null)
			polled.incrementAndGet();
		return t;
	}

	@Override
	public T get() {

		if(wereAllDataFetched())
			return null;

		T row = null;

		while(row == null){
			row = pollData();
			if(row==null){
				try {Thread.sleep(10);} catch (InterruptedException e) {e.printStackTrace();}
				continue;
			}
		}

		return row;
	}

	public int sizePipe(){
		return queue.size();
	}

	public boolean isEmpty(){
		return queue.isEmpty();
	}
	
	public long getAdded() {
		return added.get();
	}

	public long getPolled() {
		return polled.get();
	}

	/**
	 * all data in flag has to be set.
	 * @return all data
	 */
	public List<T> pollAllWithWait(){
		List<T> li = new ArrayList<>();
		while(!wereAllDataFetched()){
			T row = pollData();
			if(row==null){
				try {Thread.sleep(10);} catch (InterruptedException e) {e.printStackTrace();}
				continue;
			}
			li.add(row);
		}
		return li;
	}

	public List<T> pollBatchWithWait(int sizeBatch){
		List<T> li = new ArrayList<>();
		while(!wereAllDataFetched()){
			T row = pollData();
			if(row==null){
				try {Thread.sleep(10);} catch (InterruptedException e) {e.printStackTrace();}
				continue;
			}
			li.add(row);
			if(li.size()==sizeBatch){
				break;
			}
		}
		return li;
	}

	public List<T> pollAll(){
		if(!isAllDataIn()){
			throw new RuntimeException("all_data_in flag not set.");
		}
		List<T> li = new ArrayList<T>();
		while(!isEmpty()){
			li.add(pollData());
		}
		return li;
	}
	
	/**
	 * Returns true if all data in was set and the queue is empty.
	 */
	public boolean wereAllDataFetched() {
		if(!isAllDataIn()){
			return false;
		}
		return queue.isEmpty();
	}

	public void clear(){
		queue.clear();
	}


}
