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

package com.actelion.research.calc;

import java.util.ArrayList;

public class DataProcessor {
	private ArrayList<ProgressListener> mListener;
	private ThreadMaster	mThreadMaster;
    private boolean         mVerbose;
	
	public DataProcessor() {
		mListener = new ArrayList<ProgressListener>();
        mVerbose = true;
		}

	public void addProgressListener(ProgressListener l) {
		mListener.add(l);
		}

	public void removeProgressListener(ProgressListener l) {
		mListener.remove(l);
		}

	public void setThreadMaster(ThreadMaster t) {
		mThreadMaster = t;
		}

	public void startProgress(String message, int min, int max) {
	    if (mListener.size() == 0 && mVerbose)
			System.out.println(message);
		for (int i=0; i<mListener.size(); i++)
			mListener.get(i).startProgress(message, min, max);
		}

	public void updateProgress(int value) {
		for (int i=0; i<mListener.size(); i++)
			mListener.get(i).updateProgress(value);
		}

	public void stopProgress(String message) {
	    if (mListener.size() == 0 && mVerbose)
			System.out.println(message);
		for (int i=0; i<mListener.size(); i++)
			mListener.get(i).stopProgress();
		}

	public boolean threadMustDie() {
		return (mThreadMaster == null) ? false : mThreadMaster.threadMustDie();
		}
	
	public void setVerbose(boolean verbose) {
        mVerbose = verbose;
	    }
	}
