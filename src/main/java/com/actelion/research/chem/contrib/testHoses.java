package com.actelion.research.chem.contrib;

import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.contrib.HoseCodeCreator;

public class testHoses {

    public static void main(String[] args) {
        String[] hoses=HoseCodeCreator.getHoseCodesFromDiaID(
                "fdyA`@@LTdlmNs}Bd{sMUSUUU`a@bID@_iAHNET",
                5,
                0
        );



        for (String hose : hoses) {
            System.out.println(hose);
        }


    }

}

