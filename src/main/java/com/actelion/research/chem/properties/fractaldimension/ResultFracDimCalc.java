package com.actelion.research.chem.properties.fractaldimension;

import com.actelion.research.util.Formatter;

/**
 * ResultFractalDimensionCalculation
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 28.08.18.
 */
public class ResultFracDimCalc extends InputObjectFracDimCalc {

    public static final String TAG_SMILES = "SMILES";

    public static final String TAG_ID = "MoleculeId";

    public static final String TAG_SUM_UNIQUE_FRAGMENTS_CALC = "SumUniqueFragmentsCalculated";

    public static final String TAG_BONDS_AT_MAX_FRAGS_CALC = "BondNumberAtMaxNumFragCalculated";

    public static final String TAG_MAX_NUM_UNIQUE_FRAGMENTS_CALC = "MaxNumUniqueFragmentsCalculated";

    public static final String TAG_FRACTAL_DIM = "FractalDimension";
    public static final String TAG_MESSAGE = "Message";

    public static final String [] ARR_TAGS = {
            TAG_SMILES,
            TAG_ID,
            TAG_SUM_UNIQUE_FRAGMENTS_CALC,
            TAG_BONDS_AT_MAX_FRAGS_CALC,
            TAG_MAX_NUM_UNIQUE_FRAGMENTS_CALC,
            TAG_FRACTAL_DIM,
            TAG_MESSAGE};


    public static final String SEP = "\t";

    int idMolecule;

    double fractalDimension;

    int bondsAtMaxFrag;

    int maxNumUniqueFrags;

    int sumUniqueFrags;

    String message;

    public ResultFracDimCalc(InputObjectFracDimCalc inputObjectFracDimCalc) {
        super(inputObjectFracDimCalc);

        idMolecule = -1;

        fractalDimension = Double.NaN;

        bondsAtMaxFrag = -1;

        maxNumUniqueFrags = -1;

        sumUniqueFrags = -1;

        message = "";
    }

    public double getFractalDimension() {
        return fractalDimension;
    }

    public int getBondsAtMaxFrag() {
        return bondsAtMaxFrag;
    }

    public int getMaxNumUniqueFrags() {
        return maxNumUniqueFrags;
    }

    public int getSumUniqueFrags() {
        return sumUniqueFrags;
    }

    public String getMessage() {
        return message;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append(getSmiles());
        sb.append(SEP);
        sb.append(getId());
        sb.append(SEP);
        sb.append(getSumUniqueFrags());
        sb.append(SEP);
        sb.append(getBondsAtMaxFrag());
        sb.append(SEP);
        sb.append(getMaxNumUniqueFrags());
        sb.append(SEP);
        sb.append(Formatter.format3(getFractalDimension()));
        sb.append(SEP);
        sb.append(getMessage());

        return sb.toString();
    }

    public static String toStringHeader() {
        StringBuilder sb = new StringBuilder();

        sb.append(TAG_SMILES);
        sb.append(SEP);
        sb.append(TAG_ID);
        sb.append(SEP);
        sb.append(TAG_SUM_UNIQUE_FRAGMENTS_CALC);
        sb.append(SEP);
        sb.append(TAG_BONDS_AT_MAX_FRAGS_CALC);
        sb.append(SEP);
        sb.append(TAG_MAX_NUM_UNIQUE_FRAGMENTS_CALC);
        sb.append(SEP);
        sb.append(TAG_FRACTAL_DIM);
        sb.append(SEP);
        sb.append(TAG_MESSAGE);

        return sb.toString();
    }


}
