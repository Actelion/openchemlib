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

import com.actelion.research.util.Formatter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ContainerFragBondsSolutions {
	
	public static boolean ELUSIVE = false; 

	/** Strychnine number of solutions 20.06.2017
	1	31
            2	51
            3	100
            4	219
            5	505
            6	1172
            7	2709
            8	6167
            9	13666
            10	29323
            11	60560
            12	119880
            13	226229
            14	404703
            15	682196
            16	1076336
            17	1577022
            18	2126335
            19	2612329
            20	2886786
            21	2822260
     **/

    private static int [] ARR_CAPACITY = {
	        10, // 1 bond
            100,
            200,
            400,
            800, // 5
            1600,
            3200,
            6400,
            18 * 1000,
            36 * 1000, // 10
            70 * 1000,
            140 * 1000,
            280 * 1000,
            560 * 1000,
            1000 * 1000, // 15
            2 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000, // 20
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000, // 30
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000, // 40
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000,
            4 * 1000 * 1000 // 50

    };

	protected static double FACTOR_CAPACITY = 1.7;

	protected static int START_CAPACITY = 100;

	protected static int DEFAULT_CAPACITY = 4 * 1000 * 1000;

	protected static int MAX_NUM_BONDS = 100;

	
	private ContainerBitArray containerListFragmentDefinedByBonds;
	
	// For each number of bonds one hash map in the list.
	private List<HashMap<IBitArray, IBitArray>> liHMFragmentDefinedByBonds;
	
	private int bondsMolecule;

	private int maximumNumberBondsInFragment;

	// Equivalent to the maximum number of bonds.
	private int bits;

	
	/**
	 * Fragments are represented as bit arrays. Each bit represents a bond. The index of the bit equals the index of the
	 * bond in {@link com.actelion.research.chem.Molecule}
	 * @param bits the maximum number of bonds in the Molecule that can be stored.
	 * @param totalCapacity this is the capacity for all records. Memory is acquired until the maximum capacity is
	 * reached.
	 */
	public ContainerFragBondsSolutions(int bits, int totalCapacity) {
		this.bits = bits;
		int [] arrHashMapCapacity = getHashMapCapacity(totalCapacity);
		init(arrHashMapCapacity, totalCapacity);
	}
	
	private int [] getHashMapCapacity(int totalMaximumCapacity) {

        int totalCapacity = 0;
		int [] arrCapacity = new int [MAX_NUM_BONDS+1];
        for (int i = 0; i < arrCapacity.length; i++) {
            arrCapacity[i]=1;
        }

        for (int i = 0; i < ARR_CAPACITY.length; i++) {
            arrCapacity[i]= ARR_CAPACITY[i];
            totalCapacity += arrCapacity[i];
			maximumNumberBondsInFragment = i;
            if(totalCapacity>=totalMaximumCapacity) {
                break;
            }
        }

        if(totalCapacity < totalMaximumCapacity) {

            for (int i = ARR_CAPACITY.length; i < MAX_NUM_BONDS + 1; i++) {

                arrCapacity[i] = DEFAULT_CAPACITY;

                totalCapacity += arrCapacity[i];

				maximumNumberBondsInFragment = i;

                if (totalCapacity >= totalMaximumCapacity) {
                    break;
                }
            }
        }

		if(ELUSIVE)
			System.out.println("ContainerFragBondsSolutions maximum number of bonds in a single fragment " + maximumNumberBondsInFragment + ".");

        return arrCapacity;
	}

	private void init(int [] arrCapacity, int totalCapacity){
		
		
		int cumulatedInitHashMapCapacity = 0;
		
		if(ELUSIVE) {
			System.out.println("ContainerFragBondsSolutions Capacity");
			System.out.println("Bonds\tCapacity");
		}
		
		for (int i = 0; i < arrCapacity.length; i++) {
			
			cumulatedInitHashMapCapacity += arrCapacity[i];
			
			if(ELUSIVE) {
				int bonds = i+1;
				System.out.println(bonds + "\t" + arrCapacity[i]);
			}
		}
		
		if(ELUSIVE)
			System.out.println("ContainerFragBondsSolutions initialized cumulated hash map capacity " + Formatter.group(cumulatedInitHashMapCapacity) + ".");
		
		liHMFragmentDefinedByBonds = new ArrayList<HashMap<IBitArray,IBitArray>>();
		
		liHMFragmentDefinedByBonds.add(new HashMap<IBitArray, IBitArray>());

		for (int i = 0; i < arrCapacity.length; i++) {
			liHMFragmentDefinedByBonds.add(new HashMap<IBitArray, IBitArray>(arrCapacity[i]));
		}
		
		containerListFragmentDefinedByBonds = new ContainerBitArray(bits, totalCapacity);
		
		if(ELUSIVE)
			System.out.println("ContainerFragBondsSolutions constructor finished.");
	}

	/**
	 * The fragment is added if it is not already in the hash map.
	 * @param f
	 * @return
	 */
	public boolean addFacultative(IBitArray f){
		
		boolean added = false;
		
		calculateHash(f);

		// Number of bits set is equivalent to the number of bonds in the fragment.
		int bits = getBitsSet(f);
		
		HashMap<IBitArray, IBitArray> hm = liHMFragmentDefinedByBonds.get(bits);
		
		if(hm.containsKey(f)){
			containerListFragmentDefinedByBonds.receycle(f);
		} else {
			hm.put(f, f);
			added = true;
		}
		
		return added;
	}

	public int getBitsSet(IBitArray f){
		int bits = 0;
		for (int i = 0; i < bondsMolecule; i++) {
			if(f.isBitSet(i)){
				bits++;
			}
		}
		return bits;
	}
	
	public IBitArray getWithCopy(IBitArray orign){
		IBitArray f = containerListFragmentDefinedByBonds.get();
		f.copyIntoThis(orign);
		return f;
	}
	
	public List<IBitArray> getList(int bonds){
		HashMap<IBitArray, IBitArray> hm = liHMFragmentDefinedByBonds.get(bonds);
		return new ArrayList<>(hm.values());
	}
	

	/**
	 * @return the bondsMolecule
	 */
	public int getBondsMolecule() {
		return bondsMolecule;
	}

	/**
	 * @param bondsMolecule the bondsMolecule to set
	 */
	public void setBondsMolecule(int bondsMolecule) {
		this.bondsMolecule = bondsMolecule;
	}

	public void calculateHash(IBitArray f){
		containerListFragmentDefinedByBonds.calculateHash(f);
	}
	
	public IBitArray get() {
		return containerListFragmentDefinedByBonds.get();
	}
	
	public int getSizeBinaryArray(){
		return containerListFragmentDefinedByBonds.getSizeBinaryArray();
	}
	

	public int getTotalSizeResults (){
		int size = 0;
		
		for (HashMap<IBitArray, IBitArray> hm : liHMFragmentDefinedByBonds) {
			size += hm.size();
		}
				
		return size;
	}
	
	public void reset(){
		for (HashMap<IBitArray, IBitArray> hm : liHMFragmentDefinedByBonds) {
			hm.clear();
		}
		containerListFragmentDefinedByBonds.reset();
	}

	/**
	 * Clears the hash map with the records that contain the given number of bits set.
	 * @param bits
	 */
	public void reset(int bits){

		HashMap<IBitArray, IBitArray> hm = liHMFragmentDefinedByBonds.get(bits);

		List<IBitArray> li = new ArrayList<>(hm.keySet());

		for (IBitArray bitArray : li) {
			containerListFragmentDefinedByBonds.receycle(bitArray);
		}

		hm.clear();
	}

	/**
	 * @return the maximumNumberBondsInFragment
	 */
	public int getMaximumCapacityBondsInFragment() {
		return maximumNumberBondsInFragment;
	}
	
	public int getCapacity(){
		return containerListFragmentDefinedByBonds.getCapacity();
	}
	
	public int getAvailable(){
		return containerListFragmentDefinedByBonds.getAvailable();
	}

	public int getMaximumNumberBondsInMolecule(){
		return bits;
	}

}
