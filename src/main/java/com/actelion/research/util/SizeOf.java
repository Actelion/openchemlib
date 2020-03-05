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

import java.util.Hashtable;

/**
 * 13.02.2006 MvK
 * @author vonkorm
 *
 */
public class SizeOf
{
    private static final Runtime RUNTIME = Runtime.getRuntime ();

    public static void main (String [] args) throws Exception
    {
    	
    	int numMols = 500000;
    	int par = 3;
    	
    	Hashtable ht = new Hashtable(numMols);
    	System.out.println("Used memory " + SizeOf.usedMemory());        

    	for (int i = 0; i < numMols; i++) {
    		float[] arr = new float[par];
    		ht.put(new Integer(i), arr);
		}
     	System.out.println("Used memory " + SizeOf.usedMemory());        
    	
    	System.out.println("Finished");
    	System.exit(0);
    	
    }


    private static void runGC() throws Exception {

        long usedMem1 = usedMemory();

        long usedMem2 = Long.MAX_VALUE;

        for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++ i)  {

            RUNTIME.runFinalization();

            RUNTIME.gc();

            Thread.currentThread().yield ();
            
            usedMem2 = usedMem1;

            usedMem1 = usedMemory ();
        }
    }

    public static long usedMemory ()
    {
        return RUNTIME.totalMemory () - RUNTIME.freeMemory ();
    }

    public static long usedMemoryMB () {
        return (long)((RUNTIME.totalMemory () - RUNTIME.freeMemory ())/1000000.0 + 0.5);
    }




} // End of class
