package com.actelion.research.util;

import java.util.Hashtable;

/**
 * 13.02.2006 MvK
 * http://www.javaworld.com/javaworld/javatips/jw-javatip130.html
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
