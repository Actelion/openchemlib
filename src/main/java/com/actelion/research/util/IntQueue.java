package com.actelion.research.util;

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
