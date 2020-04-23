package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.calc.VectorSimilarity;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorConstants;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.chem.descriptor.DescriptorInfo;
import com.actelion.research.chem.descriptor.flexophore.generator.CreatorMolDistHistViz;
import com.actelion.research.chem.descriptor.flexophore.generator.VectophoreGenerator;
import com.actelion.research.util.EncoderFloatingPointNumbers;

/**
 * Vectorized Flexophore descriptor
 */
public class DescriptorHandlerVectophore implements DescriptorHandler<float[], StereoMolecule> {


    public static final String VERSION = DescriptorConstants.DESCRIPTOR_Vectophore.version;

    private static final float[] FAILED_OBJECT = new float[0];

//    Min and max values in the flaot descriptor vector. Analyzed  approx 10'000 cmpds from the DUDE dataset.
//    private static final double MIN_DESCRIPTOR_VALUE = 0;
//    private static final double MAX_DESCRIPTOR_VALUE = 800 * 1000;

    private static final int ENCODING_PRECISION_BITS = 24;

    private static DescriptorHandlerVectophore INSTANCE;


    private VectophoreGenerator vectophoreGenerator;

    public DescriptorHandlerVectophore() {
        vectophoreGenerator = new VectophoreGenerator();
    }

    @Override
    public float getSimilarity(float[] d1, float[] d2) {
        return (float) VectorSimilarity.getTanimotoSimilarity(d1, d2);
    }

    @Override
    public DescriptorInfo getInfo() {
        return DescriptorConstants.DESCRIPTOR_Vectophore;
    }

    @Override
    public String getVersion() {
        return VERSION;
    }

    @Override
    public String encode(float[] f) {
        return EncoderFloatingPointNumbers.encode(f, ENCODING_PRECISION_BITS);
    }

    @Override
    public float[] decode(String s) {
        double [] d = EncoderFloatingPointNumbers.decode(s);

        float [] f = new float[d.length];

        for (int i = 0; i < d.length; i++) {
            f[i]=(float)d[i];
        }

        return f;
    }

    @Override
    public float[] decode(byte[] bytes) {
        return decode(new String(bytes));
    }

    @Override
    public float[] createDescriptor(StereoMolecule mol) {
        float [] arrVectophore = FAILED_OBJECT;

        try {
            CreatorMolDistHistViz creatorMolDistHistViz = new CreatorMolDistHistViz();

            MolDistHistViz mdhv = creatorMolDistHistViz.createMultipleConformations(mol);

            MolDistHist mdh = mdhv.getMolDistHist();

            arrVectophore = vectophoreGenerator.create(mdh);

        } catch (Exception e) {
            e.printStackTrace();
        }

        return arrVectophore;
    }

    @Override
    public boolean calculationFailed(float[] o) {
        return o==null || o.length == 0;
    }

    public static DescriptorHandlerVectophore getDefaultInstance() {

        if (INSTANCE == null) {
            synchronized (DescriptorHandlerVectophore.class) {
                INSTANCE = new DescriptorHandlerVectophore();
            }
        }

        return INSTANCE;

    }
        @Override
    public DescriptorHandler<float[], StereoMolecule> getThreadSafeCopy() {
        return new DescriptorHandlerVectophore();
    }
}
