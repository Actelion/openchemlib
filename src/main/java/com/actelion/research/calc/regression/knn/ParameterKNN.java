package com.actelion.research.calc.regression.knn;

import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * ParameterKNN
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 06.12.18.
 */
public class ParameterKNN extends ParameterRegressionMethod {

    public static final String TAG_NEIGHBOURS="Neighbours";

    private int neighbours;

    public ParameterKNN() {
        super(ConstantsRegressionMethods.MODEL_KNN);
        setNeighbours(KNNRegression.NEIGHBOURS);
    }

    public ParameterKNN(int neighbours) {
        super(ConstantsRegressionMethods.MODEL_KNN);
        setNeighbours(neighbours);
    }

    public int getNeighbours() {
        return neighbours;
    }

    public void setNeighbours(int neighbours) {
        this.neighbours = neighbours;
        properties.setProperty(TAG_NEIGHBOURS, Integer.toString(neighbours));
    }

    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterKNN parameterKNN = (ParameterKNN)o;

        if(neighbours>parameterKNN.neighbours) {
            cmp=1;
        }else if(neighbours<parameterKNN.neighbours) {
            cmp=-1;
        }

        return cmp;
    }

    @Override
    protected void decodeProperties2Parameter() {
        neighbours = Integer.parseInt(properties.getProperty(TAG_NEIGHBOURS));
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ParameterKNN{");
        sb.append("name=").append(getName());
        sb.append("neighbours=").append(neighbours);
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        li.add(TAG_NEIGHBOURS);

        return li;
    }


    public static void main(String[] args) throws IOException {

        File dir = new File("/home/korffmo1/tmp/tmp00");

        File fiProp = new File(dir, "knn.properties");

        ParameterKNN parameterKNN = new ParameterKNN();

        parameterKNN.setNeighbours(3);

        parameterKNN.write(fiProp);


        ParameterKNN parameterKNNIn = new ParameterKNN();

        parameterKNNIn.read(fiProp);

        System.out.println(parameterKNNIn.toString());
    }

}
