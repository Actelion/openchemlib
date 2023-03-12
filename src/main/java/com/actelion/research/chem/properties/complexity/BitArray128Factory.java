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

package com.actelion.research.chem.properties.complexity;

import com.actelion.research.util.BurtleHasher;
import com.actelion.research.util.BurtleHasherABC;

public class BitArray128Factory implements IBitArrayFactory<BitArray128> {
	
	private BurtleHasherABC burtleHasherABC;
	
	/**
	 * 
	 */
	public BitArray128Factory() {
		burtleHasherABC = new BurtleHasherABC(0, 0, 0);


	}
	
	/* (non-Javadoc)
	 * @see com.actelion.research.chem.properties.complexity.IBitArrayCreator#getNew()
	 */
	@Override
	public BitArray128 getNew(int index) {
		return new BitArray128(index);
	};
	
	public void calculateHash(BitArray128 f){
		
		burtleHasherABC.a = f.l1;
		burtleHasherABC.b = f.l2;
		burtleHasherABC.c = 0;
		
		BurtleHasher.mix64(burtleHasherABC);

		int h = (int)burtleHasherABC.c;

		if(BitArray128.HASH_NOT_SET == h){
			h = BitArray128.HASH_SIDE_STEP;
			System.out.println("BitArray128Factory calculateHash(...) hash calculated from long-->int with " + BitArray128.HASH_NOT_SET + ", set to " + BitArray128.HASH_SIDE_STEP+ ".");
		}

		f.setHash(h);
		
	}


}
