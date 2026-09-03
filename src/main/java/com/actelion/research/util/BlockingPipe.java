/*
* Copyright (c) 2026
* Modest von Korff
* Zwischeb den Wegen 9
* 79540 Lörrach
* Germany
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
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;

/**
 * BlockingPipe
 * Enables concurrent access to a queue.
 * @author Modest von Korff
 * Replaces class Pipeline
 * Oct 9 2012 MvK: bug fix, added reset()
 */
public class BlockingPipe<T> implements IPipeline<T>, Supplier<T> {

	public static final int CAPACITY = 10;
	public static final long TIMEOUT = 10;
	public static final TimeUnit MS = TimeUnit.MILLISECONDS;

	private volatile AtomicBoolean allDataIn;

	private ArrayBlockingQueue<T> queue;

	private AtomicLong added;

	private AtomicLong polled;

	private TimeUnit unit;
	private long timeout;

	public BlockingPipe() {
		this(CAPACITY);
	}

	public BlockingPipe(int capacity) {
		allDataIn = new AtomicBoolean(false);
		queue = new ArrayBlockingQueue<T>(capacity);
		added = new AtomicLong();
		polled = new AtomicLong();
		setTimeOut(TIMEOUT, MS);
	}

	/**
	 * The 'all data in' flag is set true.
	 * @param li
	 */
	public BlockingPipe(List<T> li) throws InterruptedException {
		this(li.size());
		put(li);
		setAllDataIn(true);
	}

	public void setTimeOut(long timeout, TimeUnit unit){
		this.timeout = timeout;
		this.unit = unit;
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

	/**
	 * Waits
	 * @param t
	 * @throws InterruptedException
	 */
	public void put(T t) throws InterruptedException {
		if(isAllDataIn())
			throw new RuntimeException("All data in flag set!");
		queue.put(t);
		added.incrementAndGet();
	}

	private void put(List<T> li) throws InterruptedException {
		for (T t : li) {
			put(t);
		}
	}

	/**
	 *
	 * @return null if nothing is in the queue.
	 */
	public T poll() throws InterruptedException {
		if(wereAllDataFetched()){
			throw new RuntimeException("All data already fetched!");
		}
		T t = queue.poll(timeout, unit);
		if(t!=null)
			polled.incrementAndGet();
		return t;
	}

	@Override
	public T get() {
        T t = null;
        try {
			while (!wereAllDataFetched()) {
				t = poll();
				if(t!=null)
					break;
			}
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return t;
	}

	public int sizePipe(){
		return (int)(getAdded()-getPolled());
		// size queze is O(n)
		// return queue.size();
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
	public List<T> pollAll(){
		List<T> li = new ArrayList<>();
		while(!wereAllDataFetched()){
            T row = null;
            try {
                row = poll();
				if(row==null)
					continue;
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
			li.add(row);
		}
		return li;
	}

	public List<T> pollBatch(int sizeBatch){
		List<T> li = new ArrayList<>();
		while(!wereAllDataFetched()){
			T row = null;
			try {
				row = poll();
				if(row==null)
					continue;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			li.add(row);
			if(li.size()==sizeBatch){
				break;
			}
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
