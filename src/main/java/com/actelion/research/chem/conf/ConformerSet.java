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
 * @author Thomas Sander
 */

package com.actelion.research.chem.conf;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParserWithoutCoordinateInvention;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ArrayUtils;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * A ConformerSet is a set of any or just of distinct Conformers, depending on whether a TorsionDescriptorHelper
 * is provided. If one is given, then it is used to calculate a torsion descriptor for every added Copnformer
 * and to check, whether it is equal (actually highly similar) to any other conformer in the set and added only,
 * if it is substantially different. If no TorsionDescriptorHelper is given, then any Conformers are added.
 */
public class ConformerSet extends ArrayList<Conformer> {
	public static final String PLUS_ENANTIOMERS = " plusEnantiomers";

	private TorsionDescriptorHelper mDescriptorHelper;

	/**
	 * Instantiates an empty ConformerSet for new Conformers of the same molecule to be added.
	 * @param descriptorHelper null or helper used to generate torsion descriptors and to check for novelty
	 */
	public ConformerSet(TorsionDescriptorHelper descriptorHelper) {
		super();
		mDescriptorHelper = descriptorHelper;
	}

	/**
	 * Instantiates a ConformerSet which is populated by one or more Conformers encoded as SPACE-separated String
	 * of one idcode, multiple idcoordinates, and optionally a 'plusEnantiomers', which would generate an enantiomeric
	 * set of all the conformers given before that.
	 * @param s String starting with an idcode followed multiple times by ' '+id_coordinates and possible by ' plusEnantiomers'
	 * @param descriptorHelper null or helper used to generate torsion descriptors and to check for novelty
	 */
	public ConformerSet(String s, TorsionDescriptorHelper descriptorHelper) {
		super();
		mDescriptorHelper = descriptorHelper;
		if (s != null) {
			boolean plusEnantiomers = s.endsWith(PLUS_ENANTIOMERS);
			if (plusEnantiomers)
				s = s.substring(0, s.length()-PLUS_ENANTIOMERS.length());

			String[] confString = s.split(" ");
			if (confString.length >= 2) {
				IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
				StereoMolecule mol = parser.getCompactMolecule(confString[0], confString[1]);
				add(new Conformer(mol));

				for(int i=2; i<confString.length;i++) {
					Conformer conformer = new Conformer(mol, false);
					try {
						parser.parseCoordinates(confString[i].getBytes(StandardCharsets.UTF_8), 0, mol, conformer.getCoordinates());
						add(conformer);
					} catch (Exception e) {}
				}

				if (plusEnantiomers) {
					int originalCount = size();
					for (int i=0; i<originalCount; i++)
						add(get(i).createEnantiomer());
				}
			}
		}
	}

	/**
	 * Instantiates a ConformerSet which is populated by one or more Conformers encoded as idcode plus
	 * SPACE-separated byte[] with multiple conformer idcoordinates optionally followed by 'plusEnantiomers',
	 * which would generate an enantiomeric set of all the conformers given before that.
	 * @param idcode Canonizer encoded molecule
	 * @param coords multiple times Canonizer encoded 3D-coordinates, possible followed by ' plusEnantiomers'
	 * @param descriptorHelper null or helper used to generate torsion descriptors and to check for novelty
	 */
	public ConformerSet(byte[] idcode, byte[] coords, TorsionDescriptorHelper descriptorHelper) {
		super();
		if (idcode != null && coords != null) {
			IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
			StereoMolecule mol = parser.getCompactMolecule(idcode, coords);
			add(new Conformer(mol));

			int index = ArrayUtils.indexOf(coords, (byte)32,0);
			while (index != -1) {
				Conformer conformer = new Conformer(mol, false);
				try {
					parser.parseCoordinates(coords, index+1, mol, conformer.getCoordinates());
					add(conformer);
				} catch (Exception e) {}
				index = ArrayUtils.indexOf(coords, (byte)32,index+1);

				boolean isPlusEnantioners = true;
				for (int i=0; i<PLUS_ENANTIOMERS.length(); i++) {
					if (coords.length <= index+i || coords[index+i] != (byte)PLUS_ENANTIOMERS.charAt(i)) {
						isPlusEnantioners = false;
						break;
					}
				}
				if (isPlusEnantioners) {
					int originalCount = size();
					for (int i=0; i<originalCount; i++)
						add(get(i).createEnantiomer());
					break;
				}
			}
		}
	}

	/**
	 * @param size
	 * @return shallow subset
	 */
	public ConformerSet getSubset(int size) {
		ConformerSet subset = new ConformerSet(mDescriptorHelper);
		subset.addAll(this.subList(0, size));
		return subset;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (!isEmpty()) {
			Iterator<Conformer> iterator = iterator();
			StereoMolecule mol = iterator.next().toMolecule();
			Canonizer can = new Canonizer(mol, Canonizer.COORDS_ARE_3D);
			sb.append(can.getIDCode());
			sb.append(" ");
			sb.append(can.getEncodedCoordinates(true));
			boolean isPlusEnantioners = false;
			while (iterator.hasNext()) {
				Conformer conformer = iterator.next();
				if (conformer.isEnantiomer()) {
					isPlusEnantioners = true;
				}
				else {
					can.invalidateCoordinates();
					conformer.toMolecule(mol);
					sb.append(" ");
					sb.append(can.getEncodedCoordinates(true));
				}
			}
			if (isPlusEnantioners)
				sb.append(PLUS_ENANTIOMERS);
		}
		return sb.toString();
	}

	@Override
	public boolean add(Conformer c) {
		if (mDescriptorHelper != null) {
			if (!c.hasTorsionDescriptor())
				c.calculateDescriptor(mDescriptorHelper);

			for (Conformer existing : this)
				if (c.equals(existing))
					return false;
		}

		return super.add(c);
	}
}
