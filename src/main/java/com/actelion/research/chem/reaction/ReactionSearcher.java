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

package com.actelion.research.chem.reaction;


// TODO Implement Index!
public class ReactionSearcher {
	public static final String cIndexVersion = "1.0.0";

	private String      mQueryString,mReactionString;
    private Reaction    mQuery,mReaction;
    private int[]       mQueryIndex,mReactionIndex;
    

	public ReactionSearcher() {
		}


	public void setQuery(Reaction query, int[] index) {
        mQueryString = null;
        mQuery = query;
		if (index == null)
            mQueryIndex = createIndex(mQuery);
		else
            mQueryIndex = index;
		}


	public void setQuery(String rxncode, int[] index) {
        mQueryString = rxncode;
		if (index == null) {
            mQuery = ReactionEncoder.decode(rxncode, false);
			mQueryIndex = createIndex(mQuery);
			}
		else {
			mQuery = null;
            mQueryIndex = index;
			}
		}


    public void setReaction(Reaction reaction, int[] index) {
        mReactionString = null;
        mReaction = reaction;
        if (index == null)
            mReactionIndex = createIndex(mReaction);
        else
            mReactionIndex = index;
        }


    public void setReaction(String rxncode, int[] index) {
        mReactionString = rxncode;
        if (index == null) {
            mReaction = ReactionEncoder.decode(rxncode, false);
            mReactionIndex = createIndex(mReaction);
            }
        else {
            mReaction = null;
            mReactionIndex = index;
            }
        }


	public boolean isHit() {
        return false;
		}


	public int[] createIndex(Reaction rxn) {
		return null;
		}


    public static double getSimilarity(int[] index1, int[] index2) {
        return 0.0;
        }


	public static int[] getIndexFromHexString(String hex) {
		if (hex.length() == 0 || (hex.length() & 7) != 0)
			return null;

		int[] index = new int[hex.length()/8];
		for (int i=0; i<hex.length(); i++) {
			int j = i/8;
			int code = hex.charAt(i) - '0';
			if (code > 16)
				code -= 7;
			index[j] <<= 4;
			index[j] += code;
			}

		return index;
		}


	public static String getHexStringFromIndex(int[] index) {
	    if (index == null)
	        return null;

	    byte[] bytes = new byte[index.length*8];
		for (int i=0; i<index.length; i++) {
			int value = index[i];
			for (int j=7; j>=0; j--) {
				int code = value & 15;
				if (code > 9)
					code += 7;
				bytes[i*8+j] = (byte)('0'+code);
				value >>= 4;
				}
			}

		return new String(bytes);
		}

	
    public static int bitCount(int x) {
        int temp;

        temp = 0x55555555;
        x = (x & temp) + (x >>> 1 & temp);
        temp = 0x33333333;
        x = (x & temp) + (x >>> 2 & temp);
        temp = 0x07070707;
        x = (x & temp) + (x >>> 4 & temp);
        temp = 0x000F000F;
        x = (x & temp) + (x >>> 8 & temp);

        return (x & 0x1F) + (x >>> 16);
    	}
	}
