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

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;

public class MCSReactionMapper implements IReactionMapper
{
    private static int COLORTABLE[] = {
        Molecule.cAtomColorNone,
        Molecule.cAtomColorBlue,
        Molecule.cAtomColorRed,
        Molecule.cAtomColorGreen,
        Molecule.cAtomColorMagenta,
        Molecule.cAtomColorOrange,
        Molecule.cAtomColorDarkGreen,
    };
    static final int PRODUCTFLAG = 200;
    static final int REACTANTFLAG_ATOMNUMBER = 105;
    static final int PRODUCTFLAG_ATOMNUMBER = 106;
    static final int REACTANTFLAG = 100;


    private int mapIndex = 0;
    private static int colorIndex = 1;

    public static int getNextColor()
    {
        colorIndex = ++colorIndex % COLORTABLE.length;
        if (colorIndex == 0)
            ++colorIndex;
        return COLORTABLE[colorIndex];
    }

    public Reaction matchReaction(Reaction reaction)
    {
        return mapReaction(reaction,null);
    }

    @Override
    public Reaction mapReaction(Reaction reaction, SSSearcher sss)
    {
        mapIndex = 0;
        try {
            CommonSubGraphHelper.Result mcsResult = CommonSubGraphHelper.getMCS(reaction,null,sss);
            while (mcsResult != null) {
                int color = getNextColor();
                StereoMolecule fragment = mcsResult.getMolecule();
                fragment.setFragment(true);

                // Do we need this???
//                CoordinateInventor coordinateInventor = new CoordinateInventor();
//                coordinateInventor.invent(fragment);

                StereoMolecule target = reaction.getReactant(mcsResult.getReactant());
                // matchlist returns the indices on the
                // A matchlist contains the matched atom index for each atom on the target
                // (Atom-) Array[0..n] (of fragment) contains indexes of the matched target atoms
                // Fragment atom(f)[x] => atom(t)[matchList[x]] in target
                int[] matchList = findMatches(fragment, target, sss);
                if (matchList != null && matchList.length > 0) {
                    applyMaps(fragment, target, matchList);

                    // We mark temporarily the atom with an invalid atom number, so it gets excluded from the SSS in MCS
                    // TLS 11Feb2021: this is an awful hack that destroys atom lists in the SSS matching. We need to do that better in the new due ReactionMapper!!!
                    // TLS 11Feb2021: Introduced restoration of atom mass to be able to use atom mass matching by the SSSearcher
                    for (int i = 0; i < matchList.length; i++) {
                        int atom = matchList[i];
                        int t = REACTANTFLAG + target.getAtomicNo(atom);
                        target.setAtomList(atom, new int[]{t});
                        int atomMass = target.getAtomMass(atom);    // retain original atom mass
                        target.setAtomicNo(atom, REACTANTFLAG_ATOMNUMBER);
                        target.setAtomMass(atom, atomMass);
                    }

                    target = reaction.getProduct(mcsResult.getProduct());
                    fragment.setFragment(true);
                    matchList = findMatches(fragment, target, sss);
                    applySourceMapOnTarget(fragment, target, matchList);

                    for (int i = 0; i < matchList.length; i++) {
                        int atom = matchList[i];
                        int t = PRODUCTFLAG + target.getAtomicNo(atom);
                        target.setAtomList(atom, new int[]{t});
                        int atomMass = target.getAtomMass(atom);    // retain original atom mass
                        target.setAtomicNo(atom, PRODUCTFLAG_ATOMNUMBER);
                        target.setAtomMass(atom, atomMass);
                    }

                    System.out.println(mcsResult);
                    mcsResult = CommonSubGraphHelper.getMCS(reaction, null, sss);
                } else {
                    // If we did not have any matches, we could not mark any atoms with dummy atom no,
                    // so any further attempt would lead to the same result, hence endless loop!

                    break;
                }
            }

            // Reset the temporarily marked atoms
            for (int i = 0; i < reaction.getMolecules(); i++) {
                StereoMolecule mol = reaction.getMolecule(i);
                for (int j = 0; j < mol.getAllAtoms(); j++) {
                    if (mol.getAtomicNo(j) == REACTANTFLAG_ATOMNUMBER) {
                        if(mol.getAtomList(j) != null) {
                            int atomMass = mol.getAtomMass(j);    // retain original atom mass
                            mol.setAtomicNo(j, mol.getAtomList(j)[0] - REACTANTFLAG);
                            mol.setAtomMass(j, atomMass);
                            mol.setAtomList(j, null);
                        }
                    } else if (mol.getAtomicNo(j) == PRODUCTFLAG_ATOMNUMBER) {
                        if(mol.getAtomList(j) != null) {
                            int atomMass = mol.getAtomMass(j);    // retain original atom mass
                            mol.setAtomicNo(j, mol.getAtomList(j)[0] - PRODUCTFLAG);
                            mol.setAtomMass(j, atomMass);
                            mol.setAtomList(j, null);
                        }
                    }
                }
            }
            reaction.setFragment(false);
            return reaction;
        } catch (Exception e1) {
            e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return null;
        }
    }


    public void resetFragments(Reaction rxn)
    {
        for (int i = 0; i < rxn.getMolecules(); i++) {
            rxn.getMolecule(i).setFragment(false);
        }
    }

    public void removeMapping(Reaction reaction)
    {
        for (int i = 0; i < reaction.getMolecules(); i++) {
            StereoMolecule mol = reaction.getMolecule(i);
            for (int j = 0; j < mol.getAllAtoms(); j++) {
                mol.setAtomMapNo(j, 0, false);
            }
        }
    }

    private int[] findMatches(StereoMolecule fragment, StereoMolecule molecule, SSSearcher searcher)
    {
        searcher.setMol(fragment, molecule); // fragment, molecule
        int a = searcher.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode);

        ArrayList<int[]> matchList = searcher.getMatchList(); // found match list of target
        if (matchList != null && matchList.size() != 0) {
            return matchList.get(0);
        }
        return null;
    }

    private int[] highlightQuery(StereoMolecule fragment, StereoMolecule molecule, int color,SSSearcher searcher)
    {
        searcher.setMol(fragment, molecule); // fragment, molecule
        int a = searcher.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode);

        ArrayList<int[]> matchList = searcher.getMatchList(); // found match list of target
        if (matchList != null && matchList.size() != 0) {
            for (int[] matching : matchList) {
                for (int k = 0; k < matching.length; k++) {
                    molecule.setAtomColor(matching[k], color);
                }
            }
            return matchList.get(0);
        }
        return null;
    }

    private void applyMaps(StereoMolecule fragment, StereoMolecule molecule, int matchList[])
    {
        for (int atom = 0; atom < matchList.length; atom++) {
            fragment.setAtomMapNo(atom, ++mapIndex, true);
            molecule.setAtomMapNo(matchList[atom], mapIndex, true);
        }
    }

    private void applySourceMapOnTarget(StereoMolecule source, StereoMolecule target, int matchList[])
    {
        for (int atom = 0; atom < matchList.length; atom++) {
            int map = source.getAtomMapNo(atom);
            target.setAtomMapNo(matchList[atom], map, true);
        }
    }


    private boolean validateMapping(Reaction reaction)
    {
        int reactants = reaction.getReactants();
        int products = reaction.getProducts();

        int[][] rMaps = new int[reactants][];
        int[][] pMaps = new int[products][];

        for (int r = 0; r < reactants; r++) {
            StereoMolecule src = reaction.getReactant(r);
            rMaps[r] = getMappingNos(src);
        }
        for (int p = 0; p < products; p++) {
            StereoMolecule target = reaction.getProduct(p);
            pMaps[p] = getMappingNos(target);
        }


        if (!matchMappings(rMaps, pMaps))
            return false;
        return true;
    }

    private boolean matchMappings(int[][] rMaps, int[][] pMaps)
    {
        for (int i = 0; i < rMaps.length; i++) {
            for (int sourceMapNo : rMaps[i]) {
                if (sourceMapNo == 0)
                    continue;
                boolean found = false;
                for (int j = 0; j < pMaps.length; j++) {
                    for (int targetMapNo : pMaps[j]) {
                        if (sourceMapNo == targetMapNo)
                            found = true;
                    }
                }
                if (!found)
                    return false;
            }
        }
        return true;
    }

    private int[] getMappingNos(StereoMolecule molecule)
    {
        int ret[] = new int[molecule.getAllAtoms()];
        int atoms = molecule.getAllAtoms();
        for (int atom = 0; atom < atoms; atom++) {
            ret[atoms] = molecule.getAtomMapNo(atom);
        }
        return ret;
    }


    private void remove(StereoMolecule fragment, StereoMolecule reactant, StereoMolecule product)
    {
        SSSearcher searcher = new SSSearcher();
        StereoMolecule candidates[] = {
            reactant,
            product
        };
        for (StereoMolecule candidate : candidates) {
            searcher.setMol(fragment, candidate);
            searcher.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode);
            ArrayList<int[]> matchList = searcher.getMatchList();
            if (matchList != null && matchList.size() != 0) {
                int count = 0;
                for (int[] matching : matchList) {
                    for (int k = 0; k < matching.length; k++) {
                        count++;
                    }
                }
                int remove[] = new int[count];
                int idx = 0;
                for (int[] matching : matchList) {
                    for (int k = 0; k < matching.length; k++) {
                        remove[idx++] = matching[k];
                    }
                }
                Arrays.sort(remove);
                for (int i = remove.length - 1; i >= 0; i--) {
                    candidate.deleteAtom(remove[i]);
                }
            }
        }
    }
}
