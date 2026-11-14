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
import java.util.TreeSet;


public class ConformerSet extends TreeSet<Conformer> {
	public static final String PLUS_ENANTIOMERS = " plusEnantiomers";

	public ConformerSet() {
		this(null);
	}

	public ConformerSet(String s) {
		super();
		if (s != null) {
			boolean plusEnantiomers = s.endsWith(PLUS_ENANTIOMERS);
			if (plusEnantiomers)
				s = s.substring(0, s.length()-PLUS_ENANTIOMERS.length());

			String[] confString = s.split(" ");
			if (confString.length >= 2) {
				IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
				StereoMolecule mol = parser.getCompactMolecule(confString[0], confString[1]);
				TorsionDescriptorHelper tdh = new TorsionDescriptorHelper(mol, true);
				Conformer firstConformer = new Conformer(mol);
				firstConformer.calculateDescriptor(tdh);
				add(firstConformer);

				for(int i=2; i<confString.length;i++) {
					Conformer conformer = new Conformer(firstConformer);
					try {
						parser.parseCoordinates(confString[i].getBytes(StandardCharsets.UTF_8), 0, mol, conformer.getCoordinates());
						conformer.calculateDescriptor(tdh);
						add(conformer);
					} catch (Exception e) {}
				}

				if (plusEnantiomers) {
					ArrayList<Conformer> enantiomerList = new ArrayList<>();
					for (Conformer c : this) {
						Conformer enantiomer = c.createEnantiomer();
						enantiomer.calculateDescriptor(tdh);
						enantiomerList.add(enantiomer);
					}
					addAll(enantiomerList);
				}
			}
		}
	}

	public ConformerSet(byte[] idcode, byte[] coords) {
		super();
		if (idcode != null && coords != null) {
			IDCodeParserWithoutCoordinateInvention parser = new IDCodeParserWithoutCoordinateInvention();
			StereoMolecule mol = parser.getCompactMolecule(idcode, coords);
			TorsionDescriptorHelper tdh = new TorsionDescriptorHelper(mol, true);
			Conformer firstConformer = new Conformer(mol);
			firstConformer.calculateDescriptor(tdh);
			add(firstConformer);

			int index = ArrayUtils.indexOf(coords, (byte)32,0);
			while (index != -1) {
				Conformer conformer = new Conformer(firstConformer);
				try {
					parser.parseCoordinates(coords, index+1, mol, conformer.getCoordinates());
					conformer.calculateDescriptor(tdh);
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
					ArrayList<Conformer> enantiomerList = new ArrayList<>();
					for (Conformer c : this) {
						Conformer enantiomer = c.createEnantiomer();
						enantiomer.calculateDescriptor(tdh);
						enantiomerList.add(enantiomer);
					}
					addAll(enantiomerList);
					break;
				}
			}
		}
	}

	public ConformerSet getSubset(int size) {
		ConformerSet treeSet = new ConformerSet();
		this.stream().limit(size).forEach(c->treeSet.add(c));
		return treeSet;

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
}
