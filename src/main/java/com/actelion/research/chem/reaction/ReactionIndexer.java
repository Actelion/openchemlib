/*
 * Project: DD_core
 * @(#)ReactionIndexer.java
 *
 * Copyright (c) 1997- 2014
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.RXNFileParser;

public class ReactionIndexer
{
    private static boolean debug = false;
/*
    public static final int BOND_BREAK_1 = 0x0001;
    public static final int BOND_BREAK_2 = 0x0002;
    public static final int BOND_BREAK_3 = 0x0004;
    public static final int BOND_BREAK_M = 0x0008;

    public static final int BOND_CREATE_1 = 0x0010;
    public static final int BOND_CREATE_2 = 0x0020;
    public static final int BOND_CREATE_3 = 0x0040;
    public static final int BOND_CHANGE_13 = 0x0080;

    public static final int BOND_CHANGE_12 = 0x0100;
    public static final int BOND_CHANGE_21 = 0x0200;
    public static final int BOND_CHANGE_23 = 0x0400;
    public static final int BOND_CHANGE_32 = 0x0800;
    public static final int BOND_CHANGE_31 = 0x1000;
    public static final int BOND_CHANGE_A2 = 0x2000;
    public static final int BOND_CHANGE_A1 = 0x4000;
    public static final int BOND_CHANGE_O = 0x8000;

    public static final int ATOM_C_C = 0x00010000;
    public static final int ATOM_C_O = 0x00020000;
    public static final int ATOM_C_N = 0x00040000;
    public static final int ATOM_C_X = 0x00080000;
    public static final int ATOM_C_M = 0x00100000;
*/

    private static final byte CR_ONLY = 0;
    // Carbon Create Keys
    private static final byte CR_CC1 = 1;
    private static final byte CR_CC2 = 2;
    private static final byte CR_CC3 = 3;
    private static final byte CR_CCA = 4;
    private static final byte CR_CN1 = 5;
    private static final byte CR_CN2 = 6;
    private static final byte CR_CN3 = 7;
    private static final byte CR_CO1 = 8;
    private static final byte CR_CO2 = 9;
    private static final byte CR_CCl = 10;
    private static final byte CR_CX = 11;

    // Nitrogen Create Keys
    private static final byte CR_NN1 = 12;
    private static final byte CR_NN2 = 13;
    private static final byte CR_NO1 = 14;
    private static final byte CR_NO2 = 15;

    // Change Carbon Keys
    private static final byte CH_CCSD = 16;
    private static final byte CH_CCDS = 17;
    private static final byte CH_CCST = 18;
    private static final byte CH_CCTS = 19;
    private static final byte CH_CCDT = 20;
    private static final byte CH_CCTD = 21;
    private static final byte CH_CCSA = 22;
    private static final byte CH_CCAS = 23;
    private static final byte CH_CCDA = 24;
    private static final byte CH_CCAD = 25;
    private static final byte CH_CCTA = 26;
    private static final byte CH_CCAT = 27;

    private static final byte CH_CNSD = 28;
    private static final byte CH_CNDS = 29;
    private static final byte CH_CNST = 30;
    private static final byte CH_CNTS = 31;
    private static final byte CH_CNDT = 32;
    private static final byte CH_CNTD = 33;
    private static final byte CH_CNSA = 34;
    private static final byte CH_CNAS = 35;
    private static final byte CH_CNDA = 36;
    private static final byte CH_CNAD = 37;

    private static final byte CH_COSD = 38;
    private static final byte CH_CODS = 39;


    private static final byte DEL_ONLY = 40;
    private static final byte DEL_CCS = 41;
    private static final byte DEL_CCD = 42;
    private static final byte DEL_CCT = 43;
    private static final byte DEL_CCA = 44;
    private static final byte DEL_CNS = 45;
    private static final byte DEL_CND = 46;
    private static final byte DEL_CNT = 47;
    private static final byte DEL_COS = 48;
    private static final byte DEL_COT = 49;
    private static final byte DEL_CF = 50;
    private static final byte DEL_CCl = 51;
    private static final byte DEL_CBR = 52;
    private static final byte DEL_CX = 53;

    // Nitrogen Create Keys
    private static final byte DEL_NN1 = 54;
    private static final byte DEL_NN2 = 55;
    private static final byte DEL_NO1 = 56;
    private static final byte DEL_NO2 = 57;

    private static final byte CH_ONLY = 58;

    public static final byte NUMKEYS = 1 + CH_ONLY;

    public static final int NITROGEN = 7;
    public static final int OXYGEN = 8;
    public static final int FLUORINE = 9;
    public static final int CHLORINE = 17;
    public static final int BROMINE = 35;
    public static final int IODINE = 53;
    public static final int ASTATINE = 85;

    private static final int CARBON = 6;

    private static final String[] KEY_STRING = {
        "CR_ONLY",
        "CR_CC1",
        "CR_CC2",
        "CR_CC3",
        "CR_CCA",
        "CR_CN1",
        "CR_CN2",
        "CR_CN3",
        "CR_CO1",
        "CR_CO2",
        "CR_CCl",
        "CR_CX",
        "CR_NN1",
        "CR_NN2",
        "CR_NO1",
        "CR_NO2",
        "CH_CCSD",
        "CH_CCDS",
        "CH_CCST",
        "CH_CCTS",
        "CH_CCDT",
        "CH_CCTD",
        "CH_CCSA",
        "CH_CCAS",
        "CH_CCDA",
        "CH_CCAD",
        "CH_CCTA",
        "CH_CCAT",
        "CH_CNSD",
        "CH_CNDS",
        "CH_CNST",
        "CH_CNTS",
        "CH_CNDT",
        "CH_CNTD",
        "CH_CNSA",
        "CH_CNAS",
        "CH_CNDA",
        "CH_CNAD",
        "CH_COSD",
        "CH_CODS",
        "DEL_ONLY",
        "DEL_CC1",
        "DEL_CC2",
        "DEL_CC3",
        "DEL_CCA",
        "DEL_CN1",
        "DEL_CN2",
        "DEL_CN3",
        "DEL_CO1",
        "DEL_CO2",
        "DEL_CF ",
        "DEL_CCl",
        "DEL_CBR",
        "DEL_CX ",
        "DEL_NN1",
        "DEL_NN2",
        "DEL_NO1",
        "DEL_NO2",
        "CH_ONLY",
    };

    private byte[] rxnKeys = null;

    public ReactionIndexer(boolean debug)
    {
        this.debug = debug;
        init();
    }

    public ReactionIndexer()
    {
        this(false);
    }

    private void init()
    {
        rxnKeys = new byte[NUMKEYS];
    }

    public String getKeysString(Reaction rxn)
    {
        init();
        generateKeys(rxn);
        return getKeysString();

    }

    public static String getKeyName(int key)
    {
        return KEY_STRING[key];
    }
    public byte[] getKeys(Reaction rxn)
    {
        init();
        generateKeys(rxn);
        byte[] ret = new byte[rxnKeys.length];
        System.arraycopy(rxnKeys, 0, ret, 0, rxnKeys.length);
        return ret;
    }

    public boolean hasKeys()
    {
        for (int k : rxnKeys) {
            if (k != 0) {
                return true;
            }
        }
        return false;
    }

    private void addKey(int index)
    {
        if (index >= 0 && rxnKeys[index] < Byte.MAX_VALUE) {
            rxnKeys[index]++;
        }
    }

    private void initKeys()
    {
        init();
    }


    public String getFoundKeys()
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < NUMKEYS; i++) {
            if (rxnKeys[i] > 0) {
                sb.append(String.format("%s (%d)\n",KEY_STRING[i],(int)rxnKeys[i]));
            }
        }
        return sb.toString();

    }
    public static byte[] getKeysFromString(String s)
    {
        byte[] ret = new byte[NUMKEYS];
        if (s != null) {
            if (s.length() != NUMKEYS)
                throw new RuntimeException("Invalid KeyString Length");
            for (int i = 0; i < NUMKEYS; i++) {
                ret[i] = (byte)(s.charAt(i)-'0');
            }
        }

        return ret;
    }


    public static String getKeysString(byte[] keys)
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < NUMKEYS; i++) {
            sb.append(String.format("%1d", (int) keys[i]));
        }
        return sb.toString();
    }

    private String getKeysString()
    {
        return getKeysString(rxnKeys);
    }

    private int addCreateKey(int bondtype, int atomSymbol1, int atomSymbol2)
    {
        int key = CR_ONLY;
        int a1 = atomSymbol1;
        int a2 = atomSymbol2;
        if (atomSymbol2 < atomSymbol1) {
            a1 = atomSymbol2;
            a2 = atomSymbol1;
        }
        switch (a1) {
            case CARBON:
                key = getCarbonCreateKey(bondtype, a2);
                break;
            case NITROGEN:
                key = getNitroCreateKey(bondtype, a2);
                break;

        }
        return key;
    }

    private int getCarbonCreateKey(int bondtype, int atomNo)
    {
        int key = CR_ONLY;
        switch (atomNo) {
            case CARBON:
                key = getCreateCarbonCarbonKey(bondtype);
                break;

            case NITROGEN:
                key = getCreateCarbonNitroKey(bondtype);
                break;

            case OXYGEN:
                key = getCreateCarbonOxygenKey(bondtype);
                break;

            case CHLORINE:
                key = CR_CCl;
                break;

            case FLUORINE:
            case BROMINE:
            case IODINE:
            case ASTATINE:
                key = CR_CX;
                break;

        }
        return key;
    }

    private int getNitroCreateKey(int bondtype, int atomNo)
    {
        int key = CR_ONLY;
        switch (atomNo) {
            case NITROGEN:
                if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                    key = CR_NN1;
                } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                    key = CR_NN2;
                }
                break;

            case OXYGEN:
                if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                    key = CR_NO1;
                } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                    key = CR_NO2;
                }
                break;
        }
        return key;
    }

    private int getCreateCarbonCarbonKey(int bondtype)
    {
        if (bondtype == Molecule.cBondTypeDelocalized) {
            return CR_CCA;
        } else if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            return CR_CC1;
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            return CR_CC2;
        } else if ((bondtype & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            return CR_CC3;
        }
        return CR_ONLY;
    }

    private int getCreateCarbonNitroKey(int bondtype)
    {

        if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            return CR_CN1;
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            return CR_CN2;
        } else if ((bondtype & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            return CR_CN3;
        }
        return CR_ONLY;
    }

    private int getCreateCarbonOxygenKey(int bondtype)
    {
        if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            return CR_CO1;
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            return CR_CO2;
        }
        return CR_ONLY;
    }

    ///////// Change////////
    private void addChangeKey(int atomNo1, int atomNo2, int bt1, int bt2)
    {
        addKey(CH_ONLY);
        int a1 = atomNo1;
        int a2 = atomNo2;
        if (atomNo2 < atomNo1) {
            a1 = atomNo2;
            a2 = atomNo1;
        }
        switch (a1) {
            case CARBON:
                addCarbonChangeKey(a2, bt1, bt2);
                break;
        }
        //return key;
    }

    private void addCarbonChangeKey(int atomNo, int bt1, int bt2)
    {
        int key = CH_ONLY;
        switch (atomNo) {
            case CARBON:
                addChangeCarbonCarbonKey(bt1, bt2);
                break;

            case NITROGEN:
                addChangeCarbonNitroKey(bt1, bt2);
                break;

            case OXYGEN:
                addChangeCarbonOxygenKey(bt1, bt2);
                break;

        }
//        return key;
    }


    private void addChangeCarbonCarbonKey(int bt1, int bt2)
    {
//        int key = CH_ONLY;
        if (bt1 == Molecule.cBondTypeDelocalized) {
            if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CCAS);
            } else if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_CCAD);
            } else if ((bt2 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
                addKey(CH_CCAT);
            }
        } else if ((bt1 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            if (bt2 == Molecule.cBondTypeDelocalized) {
                addKey(CH_CCSA);
            } else if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_CCSD);
            } else if ((bt2 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
                addKey(CH_CCST);
            }
        } else if ((bt1 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            if (bt2 == Molecule.cBondTypeDelocalized) {
                addKey(CH_CCDA);
            } else if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CCDS);
            } else if ((bt2 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
                addKey(CH_CCDT);
            }
        } else if ((bt1 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            if (bt2 == Molecule.cBondTypeDelocalized) {
                addKey(CH_CCTA);
            } else if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CCTS);
            } else if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_CCTD);
            }
        }

//        return key;
    }

    private void addChangeCarbonNitroKey(int bt1, int bt2)
    {
//        int addKey(CH_ONLY;
        if (bt1 == Molecule.cBondTypeDelocalized) {
            if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CNAS);
            } else if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_CNAD);
            }
        } else if ((bt1 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            if (bt2 == Molecule.cBondTypeDelocalized) {
                addKey(CH_CNSA);
            } else if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_CNSD);
            } else if ((bt2 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
                addKey(CH_CNST);
            }
        } else if ((bt1 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            if (bt2 == Molecule.cBondTypeDelocalized) {
                addKey(CH_CNDA);
            } else if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CNDS);
            } else if ((bt2 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
                addKey(CH_CNDT);
            }
        } else if ((bt1 & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CNTS);
            } else if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_CNTD);
            }
        }

//        return key;
    }

    private void addChangeCarbonOxygenKey(int bt1, int bt2)
    {
//        int addKey(CH_ONLY;
        if ((bt1 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            if ((bt2 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                addKey(CH_COSD);
            }
        } else if ((bt1 & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            if ((bt2 & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                addKey(CH_CODS);
            }
        }
        //return key;
    }


    /*
        private int getChangeKey(int bt1, int bt2)
        {
            int addKey(-1;
            if (bt1 == Molecule.cBondTypeDelocalized) {
                if ((bt2 & Molecule.cBondTypeSingle) ==Molecule.cBondTypeSingle) {
                    addKey(CH_AS;
                }else if ((bt2 & Molecule.cBondTypeDouble) ==Molecule.cBondTypeDouble) {
                    addKey(CH_AD;
                }else if ((bt2 & Molecule.cBondTypeTriple) ==Molecule.cBondTypeTriple) {
                    addKey(CH_AT;
                }
            } else if ((bt1 & Molecule.cBondTypeSingle) ==Molecule.cBondTypeSingle) {
                if (bt2 == Molecule.cBondTypeDelocalized) {
                    addKey(CH_SA;
                } else if ((bt2 & Molecule.cBondTypeDouble) ==Molecule.cBondTypeDouble) {
                    addKey(CH_SD;
                }else if ((bt2 & Molecule.cBondTypeTriple) ==Molecule.cBondTypeTriple) {
                    addKey(CH_ST;
                }
            } else if ((bt1 & Molecule.cBondTypeDouble) ==Molecule.cBondTypeDouble) {
                if (bt2 == Molecule.cBondTypeDelocalized) {
                    addKey(CH_DA;
                } else if ((bt2 & Molecule.cBondTypeSingle) ==Molecule.cBondTypeSingle) {
                    addKey(CH_DS;
                }else if ((bt2 & Molecule.cBondTypeTriple) ==Molecule.cBondTypeTriple) {
                    addKey(CH_DT;
                }
            } else if ((bt1 & Molecule.cBondTypeTriple) ==Molecule.cBondTypeTriple) {
                if (bt2 == Molecule.cBondTypeDelocalized) {
                    addKey(CH_TA;
                } else if ((bt2 & Molecule.cBondTypeSingle) ==Molecule.cBondTypeSingle) {
                    addKey(CH_TS;
                }else if ((bt2 & Molecule.cBondTypeDouble) ==Molecule.cBondTypeDouble) {
                    addKey(CH_TD;
                }
            }

            return key;
        }
    */
    private void addDeleteKey(int bondtype, int atomNo1, int atomNo2)
    {
        addKey(DEL_ONLY);
        int a1 = atomNo1;
        int a2 = atomNo2;
        if (atomNo2 < atomNo1) {
            a1 = atomNo2;
            a2 = atomNo1;
        }
        switch (a1) {
            case CARBON:
                addCarbonDeleteKey(bondtype, a2);
                break;
            case NITROGEN:
                addNitroDeleteKey(bondtype, a2);
                break;

        }
    }

    private void addCarbonDeleteKey(int bondtype, int atomNo)
    {
//          int addKey(DEL_ONLY;
        switch (atomNo) {
            case CARBON:
                addDeleteCarbonCarbonKey(bondtype);
                break;

            case NITROGEN:
                addDeleteCarbonNitroKey(bondtype);
                break;

            case OXYGEN:
                addDeleteCarbonOxygenKey(bondtype);
                break;
            case FLUORINE:
            case CHLORINE:
            case 35:
            case IODINE:
                addHalogenDeleteKey(atomNo);
                break;

        }
//        addKey();
    }

    private void addDeleteCarbonCarbonKey(int bondtype)
    {
        if (bondtype == Molecule.cBondTypeDelocalized) {
            addKey(DEL_CCA);
        } else if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            addKey(DEL_CCS);
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            addKey(DEL_CCD);
        } else if ((bondtype & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            addKey(DEL_CCT);
        }
        //addKey(DEL_ONLY;
    }

    private void addDeleteCarbonNitroKey(int bondtype)
    {

        if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            addKey(DEL_CNS);
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            addKey(DEL_CND);
        } else if ((bondtype & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            addKey(DEL_CNT);
        }
        //return DEL_ONLY;
    }

    private void addDeleteCarbonOxygenKey(int bondtype)
    {
        if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            addKey(DEL_COS);
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            addKey(DEL_COT);
        }
        //return DEL_ONLY;
    }

    private void addNitroDeleteKey(int bondtype, int atomNo)
    {
//        int addKey(DEL_ONLY;
        switch (atomNo) {
            case NITROGEN:
                if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                    addKey(DEL_NN1);
                } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                    addKey(DEL_NN2);
                }
                break;

            case 8:
                if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
                    addKey(DEL_NO1);
                } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
                    addKey(DEL_NO2);
                }
                break;
        }
//        return key;
    }

    private void addHalogenDeleteKey(int atomNo)
    {
        addKey(DEL_CX);
        switch (atomNo) {
            case FLUORINE:
                addKey(DEL_CF);
                break;
            case CHLORINE:
                addKey(DEL_CCl);
                break;
            case BROMINE:
                addKey(DEL_CBR);
                break;
            case IODINE:
                break;
        }
        //return key;
    }


/*
    C1             C5                                  C1        C5
      \           /                                      \      /
       C2-Cl    N4                 -->                    C2 - N4
      /           \                                      /      \
     C3            C6                                   C3       C6

    C1                                                 C1
      \                                                  \
       C2 = O3                     -->                   C2 - 03
      /                                                  /
     C3                                                 C3

    C1         C5                                       C1            C5
      \       /                                           \          /
       C2 - N4                     -->                    C2       N4
      /      \                                           /           \
     C3       C6                                        C3           C6


    C1         C         C5                            C1         C5
      \       /        /                                 \       /
       C2 - N        N4            -->                    C2 - N4
      /      \         \                                 /       \
     C3       C         C6                             C3         C6


*/
//  C1 -> Map No = 1, pAtomNo = x, cConntectedAtoms = 1, pNumConnectedAtoms=1
//       cBondAtom = (C2), bondType = 1, cOtherMapNo = 2, pBondType=1 = > NO CHANGE
//  C2 -> Map No = 2, pAtomNo = x, cConntextedAtoms = 3, pNumConnectedAtoms = 3
//        C1 == OK
//        C3 == OK
//       cBondAtom (Cl),cBondType=1,cOtherMapNo=0
//

    /*
                            // For all atoms connected to Component Atoms (cAtomNo)
                            for (int cbi = 0; cbi < cNumConnectedAtoms; cbi++) {
                                // get the other atom at this bond
                                int cBondAtom = source.getConnAtom(cAtomNo, cbi);
                                // Get the bond index from these two atoms
                                int sourceBondIndex = source.getConnBond(cAtomNo, cbi);
                                if (cVisitedBonds[sourceBondIndex])
                                    continue;
                                // Find the mapping number of the bond's other end
                                int cOtherMapNo = source.getAtomMapNo(cBondAtom);
                                // get the bond type of this atom
                                int cBondType = getBondType(source, sourceBondIndex);
                                ////sb.append(String.format("Source Bond From %d(%d) - %d(%d) \n",cAtomNo,cMapNo,cBondAtom,cOtherMapNo));

                                cVisitedBonds[sourceBondIndex] = true;
                                // the other end is mapped
                                if (cOtherMapNo != 0) {
                                    // find the cont type of the two mapped atoms on the preoduct
                                    int pBondType = getMappedBondType(target, cMapNo, cOtherMapNo);
                                    //sb.append(String.format("Proiduct Bond Type (%d) - (%d) = %d \n",cMapNo,cOtherMapNo,pBondType));
                                    // If there's no bond, then the bond has been brokwn
                                    if (pBondType == -1) {
                                        if (!reverse) {
                                            System.out.println(String.format("###BOND_BREAK(1) at %d-%d \n", cOtherMapNo, cMapNo));
                                            int key = addDeleteKey(cBondType, source.getAtomicNo(cAtomNo), source.getAtomicNo(cBondAtom));
                                            if (debug)
                                                System.out.println("Delete Key is " + key);
                                            addKey(key);
                                            addKey(DEL_ONLY);
                                        } else {
                                            System.out.println(String.format("###BOND_CREATE(1) at %d-%d \n", cOtherMapNo, cMapNo));
                                            int key = addCreateKey(cBondType, source.getAtomicNo(cAtomNo), source.getAtomicNo(cBondAtom));
                                            addKey(key);
                                            addKey(CR_ONLY);
                                        }
                                        if (debug)
                                            sb.append(String.format("Bond %d %s-%s\n", sourceBondIndex, source.getAtomLabel(cAtomNo), source.getAtomLabel(cBondAtom)));
                                        // BOND_BREAK;
                                        // If the bond type is the same, nothing happend
                                    } else if (cBondType == pBondType) {
                                        ;
                                        // the bond type is different, so it's a change
                                    } else {
                                        if (!reverse) {
                                            System.out.println(String.format("###BOND_CHANGE(1) Type=%d Map %d-%d \n", pBondType, cOtherMapNo, cMapNo));
                                            int key = getChangeKey(source.getAtomicNo(cAtomNo), source.getAtomicNo(cBondAtom), cBondType, pBondType);
                                            addKey(key);
                                            addKey(CH_ONLY);
                                        }
                                        // BOND_CHANGE
                                    }

                                } else {
                                    System.out.println(String.format("not found other map no %d %d", listIndex, cMapNo));
                                    //if (!cVisitedBonds[sourceBondIndex]) {
                                    //    System.out.println(String.format("## UNDEFINED",listIndex,cMapNo));
                                    //}
                                    // Does this bond break
                                }
                            }
                            if (!reverse) {
                                if (cNumConnectedAtoms < pNumConnectedAtoms) {
                                    System.out.println(String.format("searching backward? %d/%d (%d)", cNumConnectedAtoms, pNumConnectedAtoms, cMapNo));
                                    int bondIndex = findDeletedBond(target, pAtomNo, source, cAtomNo);
                                    if (bondIndex != -1 && !cVisitedBonds[bondIndex]) {
                                        int ba = target.getConnAtom(pAtomNo, bondIndex);
                                        int type = target.getAtomicNo(ba);
                                        int bt = getBondType(target, bondIndex);
                                        System.out.println(String.format("### BOND_CREATE(2) at Target Atom %d(%d)[%d] from Source Atom %d(%d)", pAtomNo, cMapNo, type, cAtomNo, cMapNo));
                                        // BOND_CREATE
                                        int key = addCreateKey(bt, target.getAtomicNo(pAtomNo), target.getAtomicNo(ba));
                                        addKey(key);
                                        addKey(CR_ONLY);
                                    }
                                } else if (cNumConnectedAtoms > pNumConnectedAtoms) {
                                    System.out.println(String.format("searching forward? %d/%d (%d)", cNumConnectedAtoms, pNumConnectedAtoms, cMapNo));
                                    int bondIndex = findDeletedBond(source, cAtomNo, target, pAtomNo);
                                    if (bondIndex != -1 && !cVisitedBonds[bondIndex]) {
                                        int ba = source.getConnAtom(cAtomNo, bondIndex);
                                        int type = source.getAtomicNo(ba);
                                        int bt = getBondType(target, bondIndex);
                                        System.out.println(String.format("### BOND_DELETE(2) at Target Atom %d(%d)[%d] from Source Atom %d(%d)", pAtomNo, cMapNo, type, cAtomNo, cMapNo));
                                        // BOND_CREATE
                                        int key = addDeleteKey(bt, source.getAtomicNo(cAtomNo), source.getAtomicNo(ba));
                                        addKey(key);
                                        addKey(DEL_ONLY);
                                    }
                                }
                            }
                        } else {
                        }
                    }
                }
            } else {

            }
            sb.append(getKeysString());
            return sb.toString();
        }
    */
    private int findDeletedBond(ExtendedMolecule source, int sourceAtom, ExtendedMolecule target, int targetAtom)
    {
        int sourceNumConnectedAtoms = source.getConnAtoms(sourceAtom);
        int originAtomType = source.getAtomicNo(sourceAtom);
        int targetNumConnectedAtoms = target.getConnAtoms(targetAtom);
        boolean visitedSourceBonds[] = new boolean[sourceNumConnectedAtoms];
        if (sourceNumConnectedAtoms > targetNumConnectedAtoms) {
            if (debug) {
                System.out.println(String.format("findDeletedBond() walk source atom %d", sourceAtom));
            }
// Walk all bonds at this atom
            boolean visitedTargetBonds[] = new boolean[targetNumConnectedAtoms];
            for (int i = 0; i < sourceNumConnectedAtoms; i++) {
                if (!visitedSourceBonds[i]) {
                    // This is the atom of the bond
                    int sourceBondAtom = source.getConnAtom(sourceAtom, i);
// And this is the bond
                    int sourceBond = source.getConnBond(sourceAtom, i);
// The bond type
                    int sourceBondType = getBondType(source, sourceBond);
// The atom type
                    int sourceAtomType = source.getAtomicNo(sourceBondAtom);
                    if (debug) {
                        System.out.println(String.format("findDeletedBond() walk source bond atom %d bond %d type %d atomtype %d", sourceBondAtom, i, sourceBondType, sourceAtomType));
                    }
                    for (int j = 0; j < targetNumConnectedAtoms; j++) {
                        if (!visitedTargetBonds[j]) {
                            int bondAtom = target.getConnAtom(targetAtom, j);
                            int bond = target.getConnBond(targetAtom, j);
                            int targetBondType = getBondType(target, bond);
                            int targetAtomType = target.getAtomicNo(bondAtom);
                            if (debug) {
                                System.out.println(String.format("findDeletedBond() walk checking bond %d->%d (%d/%d) %d/%d", sourceAtomType, targetAtomType, i, j, sourceBondType, targetBondType));
                            }
                            if (sourceAtomType == targetAtomType && sourceBondType == targetBondType) {
                                visitedSourceBonds[i] = true;
                                visitedTargetBonds[j] = true;
                                System.out.println("Found Bond");
                                break;
                            }
                        }
                    }
                    if (!visitedSourceBonds[i]) {
                        if (debug) {
                            System.out.println(String.format("findDeletedBond() walk source bond not found atom %d bond %d type %d atomtype %d->%d",

                                sourceBondAtom, i, sourceBondType, originAtomType, sourceAtomType));
                        }
                    }
                }
            }
            for (int i = 0; i < sourceNumConnectedAtoms; i++) {
                if (!visitedSourceBonds[i]) {
                    if (debug) {
                        System.out.println(String.format("Did not find bond: %d(%d)", sourceAtom, i));
                    }
// this bond has not been seen, so it's breacking
                    return i;
                }
            }
        }
        return -1;
    }

    private boolean findBond(ExtendedMolecule mol, int atomIndex, int atomType, int bondType)
    {
        int pNumConnectedAtoms = mol.getConnAtoms(atomIndex);
        return false;
    }

    /*
        ExtendedMolecule findMappedAtomMol(int mapno)
        {
            ExtendedMolecule ret = null;
            int prod = rxn.getProducts();
            for (int i = 0; i < prod; i++) {
                ExtendedMolecule m = rxn.getMolecule(i);
                for (int a = 0; a < m.getAllAtoms(); a++) {
                    if (m.getAtomMapNo(a) == mapno) {
                        ret = m;
                        break;
                    }
                }
            }
            return ret;
        }
    */
    int getMappedAtom(ExtendedMolecule m, int map1)
    {
        int a1 = -1;
        for (int a = 0; a < m.getAllAtoms(); a++) {
            if (m.getAtomMapNo(a) == map1) {
                a1 = a;
                break;
            }
        }
        return a1;
    }

    int getMappedBondType(ExtendedMolecule m, int map1, int map2)
    {
        int ret = -1; // Bond Type
        int bond = getMappedBond(m, map1, map2);
        if (bond != -1) {
            return getBondType(m, bond);
        }
        return ret;
    }

    private int getBondType(ExtendedMolecule m, int bond)
    {
        int ret = m.getBondType(bond);
        if (m.isAromaticBond(bond) || m.isDelocalizedBond(bond)) {
            ret = Molecule.cBondTypeDelocalized;
        }
        return getSimpleBondType(ret);
    }

    private int getSimpleBondType(int bondtype)
    {
        if (bondtype == Molecule.cBondTypeDelocalized) {
            return Molecule.cBondTypeDelocalized;
        } else if ((bondtype & Molecule.cBondTypeSingle) == Molecule.cBondTypeSingle) {
            return Molecule.cBondTypeSingle;
        } else if ((bondtype & Molecule.cBondTypeDouble) == Molecule.cBondTypeDouble) {
            return Molecule.cBondTypeDouble;
        } else if ((bondtype & Molecule.cBondTypeTriple) == Molecule.cBondTypeTriple) {
            return bondtype & Molecule.cBondTypeTriple;
        }
        return bondtype;
    }

    int getMappedBond(ExtendedMolecule m, int map1, int map2)
    {
        int ret = -1; // Bond
        int a1 = getMappedAtom(m, map1);
        if (a1 != -1) {
            int bcnt = m.getConnAtoms(a1);
            for (int b = 0; b < bcnt; b++) {
                int a2 = m.getConnAtom(a1, b);
//                int a2 = m.getBondAtom(a1,b);
                int mt = m.getAtomMapNo(a2);
                if (debug) {
                    System.out.println(String.format("Checking Product bond %d-%d : Product Atom %d has map no %d", a1, a2, a2, mt));
                }
                if (mt == map2) {
                    ret = m.getConnBond(a1, b);
                    break;
                }
            }
        }
        return ret;
    }

//    public int getKey(ExtendedMolecule s, int a1, int a2, int bond)
//    {
//        return 0;
//    }


    private int getSameBond(ExtendedMolecule source, int cBondAtom, int cBondIndex, ExtendedMolecule target, int targetAtom)
    {
        int cAtomSymbol = source.getAtomicNo(cBondAtom);
        int map = source.getAtomMapNo(cBondAtom);
//        int cBondAtom = source.getConnAtom(sourceAtom,cbi);
//		int cBondIndex = source.getConnBond(sourceAtom,cbi);
        int cBondType = getBondType(source, cBondIndex);
        int targetBondCount = target.getConnAtoms(targetAtom);
        for (int tbi = 0; tbi < targetBondCount; tbi++) {
            int tBondAtom = target.getConnAtom(targetAtom, tbi);
            int tBondIndex = target.getConnBond(targetAtom, tbi);
            int tBondType = getBondType(target, tBondIndex);
            int tAtomSymbol = target.getAtomicNo(tBondAtom);
            int tmap = target.getAtomMapNo(tBondAtom);
            if (cAtomSymbol == tAtomSymbol && cBondType == tBondType && map == tmap) {
                return tBondIndex;
            }
        }
        return -1;
    }

    private int getBond(ExtendedMolecule mol, int atom1, int atom2)
    {
        int ret = -1;
        int bondCount = mol.getConnAtoms(atom1);
        for (int tbi = 0; tbi < bondCount; tbi++) {
            int tBondAtom = mol.getConnAtom(atom1, tbi);
            if (tBondAtom == atom2) {
                ret = mol.getConnBond(atom1, tbi);
                break;
            }
        }
        return ret;
    }

    private int findMapAt(ExtendedMolecule mol, int atom, int map)
    {
        int ret = -1;
        int bondCount = mol.getConnAtoms(atom);
        int m = mol.getAtomMapNo(atom);
        if (m == map) {
            ret = atom;
        }
        if (ret == -1) {
            for (int tbi = 0; tbi < bondCount; tbi++) {
                int tBondAtom = mol.getConnAtom(atom, tbi);
                m = mol.getAtomMapNo(tBondAtom);
                if (m == map) {
                    ret = tBondAtom;
                    break;
                }
            }
        }
        return ret;
    }

    private void generateKeys(Reaction rxn)
    {

        // (1)
        // FIND ATOM WITH MAP NO IN TARGET
        // FOR ALL ATOMS ON TARGET
        // IF MAP(ATOM)== MAP
        // TATOM = ATOM
        // END FOR

        //(2)
        // for all bonds at targetatom
        // if map(atom)== map
        // return atom
        // end for

        StringBuilder sb = new StringBuilder();
        ExtendedMolecule m = new ExtendedMolecule();
        initKeys();
        if (debug) {
            System.out.println(" " + m.toString());
        }
        int comp = rxn.getReactants();
        int prod = rxn.getProducts();
        int flags = Molecule.cHelperBitNeighbours | Molecule.cHelperBitRings;

        if (comp != 0 && prod != 0) {
            if (debug) {
                System.out.println("RXNSearcher.search()");
            }
            ExtendedMolecule component = new StereoMolecule(rxn.getMolecule(0));
            for (int c = 1; c < comp; c++) {
                ExtendedMolecule mol = rxn.getMolecule(c);
                component.addMolecule(mol);
            }
            component.ensureHelperArrays(flags);
            //sb.append(String.format("Component has %d atoms\n",component.getAllAtoms()));

            ExtendedMolecule product = new StereoMolecule(rxn.getProduct(0));
            for (int c = 1; c < prod; c++) {
                product.addMolecule(rxn.getProduct(c));
            }
            product.ensureHelperArrays(flags);

            ExtendedMolecule source = component;
            ExtendedMolecule target = product;

            // All Mapped atom indizes
            int mappedSourceAtoms[] = new int[source.getAllAtoms()];
            int mappedTargetAtoms[] = new int[target.getAllAtoms()];
            boolean cVisitedBonds[] = new boolean[source.getAllBonds()];

            //sb.append(String.format("source has %d atoms\n",source.getAllAtoms()));
            int noMappedSourceAtoms = 0;
            for (int a = 0; a < source.getAllAtoms(); a++) {
                int map = source.getAtomMapNo(a);
                if (map > 0) {
                    mappedSourceAtoms[noMappedSourceAtoms++] = a;
                }

            }


            int noMappedTargetAtoms = 0;
            for (int a = 0; a < target.getAllAtoms(); a++) {
                int map = target.getAtomMapNo(a);
                if (map > 0) {
                    mappedTargetAtoms[noMappedTargetAtoms++] = a;
                }

            }

            // setup SOURCEBOND and TARGETBOND ARRAY based on MAPPED ATOMS
            boolean sourceBonds[] = new boolean[source.getAllBonds()];
            boolean targetBonds[] = new boolean[target.getAllBonds()];


//            for (int k = 0; k < idx; k++) {
////                System.out.printf("Map [%d] -> [%d]\n",k,source.getAtomMapNo(k));
//                int s = source.getConnAtoms(k);
//                for (int b = 0; b < s; b++) {
//                    int a2 = source.getConnAtom(k, b);
//                    System.out.printf("Map #[%d] -> M[%d] = #[%d] -> M[%d]\n",
//                        k,source.getAtomMapNo(k),a2,
//                        source.getAtomMapNo(a2));
//                }
//            }
//            if (true)
//                return "";

            // for all mapped atoms
            for (int k = 0; k < noMappedSourceAtoms; k++) {
                // sourceAtom
                int sourceAtom = mappedSourceAtoms[k];
                int sourceMapNo = source.getAtomMapNo(sourceAtom);
                // targetAtom = find atom with same map no as sourceAtom (1)
                int targetAtom = getMappedAtom(target, sourceMapNo);
                debug("Working on atom %d (%d) target atom is %d\n", sourceMapNo, sourceAtom, targetAtom);

                // if not found
                // break;
                if (targetAtom == -1) {
                    break;
                }

//                debug("Map first (%d) (%d-%d)\n",sourceMapNo, sourceAtom,targetAtom);

                // boundCounts at source and target atom
                int sourceBondCount = source.getConnAtoms(sourceAtom);
                int targetBondCount = target.getConnAtoms(targetAtom);
                // for all bonds at source Atom
                for (int cbi = 0; cbi < sourceBondCount; cbi++) {
                    int sourceBondConnectedAtom = source.getConnAtom(sourceAtom, cbi);
                    int sourceBondIndex = source.getConnBond(sourceAtom, cbi);
                    int sourceBondType = getBondType(source, sourceBondIndex);
                    // if sourceBond is not visited
                    if (!sourceBonds[sourceBondIndex]) {
                        // amap = map(bondAtom)
                        int sourceBondConnectedAtomMap = source.getAtomMapNo(sourceBondConnectedAtom);
                        debug("\tChecking Source Bond on atom %d (%d)\n",
                            sourceBondConnectedAtomMap, sourceBondConnectedAtom);

                        // if bondAtom is mapped
                        if (sourceBondConnectedAtomMap != 0) {
                            // find atom with amap as at targetAtom bonds
                            int targetBondAtom = findMapAt(target, targetAtom, sourceBondConnectedAtomMap);
//                            debug("Mapping source bond %d (%d-%d) (%d-%d)\n",
//                                cbi,
//                                sourceMapNo,
//                                sourceBondConnectedAtomMap,
//                                sourceAtom,
//                                sourceBondConnectedAtom
//                            );
                            // if same map is found
                            if (targetBondAtom != -1) {
                                int targetBondIndex = getBond(target, targetAtom, targetBondAtom);
                                // TARGETBOND is visited
                                targetBonds[targetBondIndex] = true;
                                // if bondtype differs
                                int targetBondType = getBondType(target, targetBondIndex);
                                // The bond types don't match
                                if (targetBondType != sourceBondType) {
                                    // its a change
                                    // sourceBond[i] is visited
                                    if (debug) {
                                        System.out.println(String.format("Map: %d-%d bond CHANGE %d/%d", sourceMapNo, sourceBondConnectedAtomMap, sourceBondType, targetBondType));
                                    }
                                    sourceBonds[sourceBondIndex] = true;
                                    addChangeKey(source.getAtomicNo(sourceAtom), source.getAtomicNo(sourceBondConnectedAtom), sourceBondType, targetBondType);
                                    //addKey(key);
                                    //addKey(CH_ONLY);
                                    // else
                                } else {
                                    // sourceBond[i] is visited
                                    sourceBonds[sourceBondIndex] = true;
                                    // no change
                                    // end if
                                }
                                // else (map not found)
                            } else {
                                // its a break;
                                // sourceBond[i] is visited
                                if (debug) {
                                    System.out.println(String.format("Map: %d-%d bond BREAK %d", sourceMapNo, sourceBondConnectedAtomMap, sourceBondType));
                                }
                                sourceBonds[sourceBondIndex] = true;
                                addDeleteKey(sourceBondType, source.getAtomicNo(sourceAtom), source.getAtomicNo(sourceBondConnectedAtom));
//                                addKey(key);
//                                addKey(DEL_ONLY);
                                // ??? don't mark bond as visited we want to find the CREATE
                                // end if
                            }
                        } else {
                            ;
                            // Not a mapped atom, we will deal with this in the next section
                            // System.out.println(String.format("Map: %d Not mapped atom", sourceMapNo));
                            // not mapped atom
                            // ignore
                        }
                    } // end if
                } // end for

                // Deal with the non-mapped atoms...
                // for all bonds at source Atom
                for (int cbi = 0; cbi < sourceBondCount; cbi++) {
                    int cBondAtom = source.getConnAtom(sourceAtom, cbi);
                    int cBondIndex = source.getConnBond(sourceAtom, cbi);
                    int cBondType = getBondType(source, cBondIndex);
//                    int cAtomSymbol = source.getAtomicNo(cBondAtom);
                    // if sourceBond is not visited
                    // System.out.println(String.format("Map %d running at atom %d map %d",sourceMapNo,cBondAtom,source.getAtomMapNo(cBondAtom)));
                    if (!sourceBonds[cBondIndex]) {
                        int cbMap = source.getAtomMapNo(cBondAtom);
                        // if bondAtom is not  mapped
                        if (cbMap == 0) {
                            int tBondIndex = getSameBond(source, cBondAtom, cBondIndex, target, targetAtom);
                            if (tBondIndex != -1) {
                                // sourcebond is visited
                                // TARGETBOND IS VISITED
                                sourceBonds[cBondIndex] = true;
                                targetBonds[tBondIndex] = true;
                                // System.out.println(String.format("No Map/No change -> found bond index %d ", sourceMapNo));
                                // no change
                                // else
                            } else {
                                // BREAK_BOND
                                // sourcebond is visited
//                                if (debug)
//                                    System.out.println(String.format("Map: %d (%d-%d) Atom Type: %d  bond at atom " +
//                                        "BREAK",
//                                        sourceMapNo,
//                                        sourceAtom,cBondAtom,
//                                        source.getAtomicNo(cBondAtom)));
                                sourceBonds[cBondIndex] = true;
                                addDeleteKey(cBondType, source.getAtomicNo(sourceAtom), source.getAtomicNo(cBondAtom));
                                if (debug) {
                                    System.out.println(String.format("Map: %d (%d-%d) Atom Type: %d  bond at atom " +
                                            "DELETE",
                                        sourceMapNo,
                                        sourceAtom, cBondAtom,
                                        source.getAtomicNo(cBondAtom)));
                                }

//                                addKey(key);
                                //addKey(DEL_ONLY);
                                // end if
                            }
                        } else {
                            // This should not happen, since we dealt with this already before
                            if (debug) {
                                System.out.println(String.format("Map: %d Found SHOULD NOT", cbMap));
                            }
                        } // end if
                    } else {
                        // The bond has already been
                        // else
                        //???
                    }
                    // end if
                    // end for
                }
                // for all target bounds which have not yet been visited
                // they might be created, however, we'd need to check whether there's a mapping
//                if (sourceBondCount < targetBondCount) {
                // for all bonds at target atom
                for (int tbi = 0; tbi < targetBondCount; tbi++) {
                    int tBondIndex = target.getConnBond(targetAtom, tbi);
                    int tBondType = getBondType(target, tBondIndex);
                    int tAtom = target.getConnAtom(targetAtom, tbi);
                    int tmap = target.getAtomMapNo(tAtom);
                    int an = target.getAtomicNo(tAtom);
                    // if targetbond is not visited
                    if (!targetBonds[tBondIndex]) {
                        //those are CREATE_BONDS
                        if (debug) {
                            System.out.println(String.format("Map: %d-%d (%d) bond at atom CREATE_2", sourceMapNo, tmap, an));
                        }
                        targetBonds[tBondIndex] = true;
                        int key = addCreateKey(tBondType, target.getAtomicNo(targetAtom), target.getAtomicNo(tAtom));
                        addKey(key);
                        addKey(CR_ONLY);

                    } //endif
                } // end for
//                } //endif
            } //end for all mapped atoms
        }
    }

    public void debug(String format, Object... args)
    {
        if (debug)
            System.out.printf(format, args);

    }

    public static void main(String args[])
    {
        ReactionIndexer rs = new ReactionIndexer();
        try {

            for (int i = 1; i <= 11; i++) {
                RXNFileParser p = new RXNFileParser();
                Reaction r = new Reaction();
                p.parse(r, new java.io.File("RXN" + i + ".rxn"));
                if (debug) {
                    System.out.println("RXN: " + r.getMolecules());
                }
                rs.getKeys(r);
            }
        } catch (Throwable e) {
            System.err.println("Error parsing RDFILE " + e);
            e.printStackTrace();
        }
    }


}


/*
class RDFileReader
{
    BufferedReader rd = null;
    StringBuilder rxn = null;
    StringBuilder data = null;

    public RDFileReader(InputStream is) throws IOException
    {
        rd = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
        String line;
        line = rd.readLine();
        if (line == null || !line.startsWith("$RDFILE "))
            throw new IOException("Invalid File Header");
        line = rd.readLine();
        if (line == null)
            throw new IOException("Invalid File Header");
        line = rd.readLine();
        if (line == null || !line.startsWith("$RFMT"))
            throw new IOException("File is empty");
    }

    boolean hasNext()
    {
        boolean ret = false;
        StringBuilder r = new StringBuilder();
        StringBuilder d = new StringBuilder();
        // The default buffer is the reaction buffer
        StringBuilder sb = r;
        String line;
        try {
            boolean eof = true;
            while ((line = rd.readLine()) != null) {
                eof = false;
                if (line.startsWith("$RFMT")) {
//                        System.out.println("End of This RXN..." +line);
                    break;
                }
                // Switch the buffer to the data buffer
                if (line.startsWith("$DTYPE")) {
//                        System.out.println("Switching To Data Buffer " + line);
                    sb = data;
                }
                sb.append(line);
                sb.append("\n");
            }
            if (!eof) {
                data = d;
                rxn = r;
                ret = true;
            }
        } catch (IOException e) {
            ret = false;
        }
        return ret;
    }

    public Reaction getReaction()
    {
        Reaction r = new Reaction();
        RXNFileParser p = new RXNFileParser();
        try {
            r = p.getReaction(rxn.toString());
        } catch (Exception e) {
            System.err.println("Error parsing reaction...");
            r = null;
        }
        return r;
    }

}
*/
