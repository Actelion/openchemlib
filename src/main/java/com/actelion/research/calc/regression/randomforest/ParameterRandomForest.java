package com.actelion.research.calc.regression.randomforest;

import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

/**
 * ParameterRandomForest
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 14.01.19.
 */
public class ParameterRandomForest extends ParameterRegressionMethod {


    public static final int NUM_TREES = 100;
    public static final double MTRY = 0.333;
    public static final int NODE_SIZE = 5;
    public static final int MAX_NODES = 50;

    public static final String TAG_TREES="NumTrees";
    public static final String TAG_MTRY="MTry";
    public static final String TAG_NODE_SIZE="NodeSize";
    public static final String TAG_MAX_NODES="MaxNodes";



    private int nTrees;

    /**
     Number of variables available for splitting at each tree node. In the random forests literature, this is referred
     to as the mtry parameter. default for regression models, it is the number of predictor variables divided by 3 (rounded down).

     There is extensive discussion in the literature about the influence of mtry. Cutler et al. (2007) reported that
     different values of mtry did not affect the correct classification rates of their model and that other performance
     metrics (sensitivity, specificity, kappa, and ROC AUC) were stable under different values of mtry. On the other
     hand, Strobl et al. (2008) reported that mtry had a strong influence on predictor variable importance estimates.
     Due to the conflicting evidence reported in the literature, we suggest you start with the default value (i.e.
     leave this parameter blank) but review the literature carefully and form your own opinion about what value might
     be suitable for your specific model.
     http://code.env.duke.edu/projects/mget/export/HEAD/MGET/Trunk/PythonPackage/dist/TracOnlineDocumentation/Documentation/ArcGISReference/RandomForestModel.FitToArcGISTable.html
     */

    // mtry the number of input variables to be used to determine the decision (split) at a node of the tree.
    // p/3 seems to give generally good performance, where p is the number of variables.
    // public int mTry;


    // The number of variables is not fix. So we take the fraction.
    private double fractionMTry;

    // The number of train instances in a node below which the tree will not split,
    // setting nodeSize = 5 generally gives good results.
    // https://stats.stackexchange.com/questions/158583/what-does-node-size-refer-to-in-the-random-forest
    private int nodeSize;

    // the maximum number of leaf nodes in the tree.
    private int maxNodes;


    public ParameterRandomForest() {
        super(ConstantsRegressionMethods.MODEL_RND_FOREST);
        initialize();
    }

    private void initialize(){
        setFractionMTry(MTRY);

        setMaxNodes(MAX_NODES);

        setNodeSize(NODE_SIZE);

        setNumTrees(NUM_TREES);
    }


    public int getNumberOfTrees() {
        return nTrees;
    }

//    public int getNumberInputVariables() {
//        return mTry;
//    }



    public double getFractionMTry() {
        return fractionMTry;
    }

    public int getNodeSize() {
        return nodeSize;
    }

    public int getMaxNodes() {
        return maxNodes;
    }


    public void setFractionMTry(double fractionMTry) {

        this.fractionMTry = fractionMTry;

        properties.put(TAG_MTRY, Double.toString(fractionMTry));

    }

    public void setNumTrees(int nTrees) {

        this.nTrees = nTrees;

        properties.put(TAG_TREES, Integer.toString(nTrees));

    }

//    public void setMTry(int mtry) {
//        this.mTry = mtry;
//    }

    public void setNodeSize(int nodeSize) {
        this.nodeSize = nodeSize;

        properties.put(TAG_NODE_SIZE, Integer.toString(nodeSize));

    }

    public void setMaxNodes(int maxNodes) {
        this.maxNodes = maxNodes;

        properties.put(TAG_MAX_NODES, Integer.toString(maxNodes));

    }

    @Override
    protected void decodeProperties2Parameter() {

        nTrees = Integer.parseInt(properties.getProperty(TAG_TREES));

        fractionMTry = Double.parseDouble(properties.getProperty(TAG_MTRY));

        nodeSize = Integer.parseInt(properties.getProperty(TAG_NODE_SIZE));

        maxNodes = Integer.parseInt(properties.getProperty(TAG_MAX_NODES));

    }

    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterRandomForest p = (ParameterRandomForest)o;

        if(nTrees > p.nTrees) {
            cmp=1;
        } else if(nTrees < p.nTrees) {
            cmp=-1;
        }

        if(cmp==0){
            if(fractionMTry > p.fractionMTry) {
                cmp=1;
            } else if(fractionMTry < p.fractionMTry) {
                cmp=-1;
            }
        }

        if(cmp==0){
            if(nodeSize > p.nodeSize) {
                cmp=1;
            } else if(nodeSize < p.nodeSize) {
                cmp=-1;
            }
        }

        if(cmp==0){
            if(maxNodes > p.maxNodes) {
                cmp=1;
            } else if(maxNodes < p.maxNodes) {
                cmp=-1;
            }
        }


        return cmp;
    }

    @Override
    public String toString() {

        DecimalFormat df = new DecimalFormat("0.0###");

        final StringBuilder sb = new StringBuilder("ParameterRNDForest{");
        sb.append("name=").append(getName());
        sb.append(" trees=").append(nTrees);
        sb.append(" frac variables=").append(fractionMTry);
        sb.append(" nodeSize=").append(nodeSize);
        sb.append(" maxNodes=").append(maxNodes);
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        li.add(TAG_TREES);
        li.add(TAG_MTRY);
        li.add(TAG_NODE_SIZE);
        li.add(TAG_MAX_NODES);

        return li;
    }

    public static void main(String[] args) throws IOException {

        File dir = new File("/home/korffmo1/tmp/tmp00");

        File fiProp = new File(dir, "randomForest.properties");

        ParameterRandomForest parameter = new ParameterRandomForest();

        parameter.nTrees = 1000;
        parameter.fractionMTry = 0.321;
        parameter.nodeSize = 7;
        parameter.maxNodes = 42;

        parameter.write(fiProp);

        ParameterRandomForest parameterIn = new ParameterRandomForest();

        parameterIn.read(fiProp);

        System.out.println(parameterIn.toString());

    }


}
