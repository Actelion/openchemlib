package com.actelion.research.calc.regression.linear.pls;

import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Properties;

/**
 * ParameterPLS
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 06.12.18.
 */
public class ParameterPLS extends ParameterRegressionMethod {

    public static final int FACTORS= 3;

    public static final String TAG_FACTORS="Factors";

    public static final boolean CENTER_DATA = true;

    protected int factors;

    private boolean centerData;

    public ParameterPLS() {
        super(ConstantsRegressionMethods.MODEL_PLS);

        setFactors(FACTORS);

        centerData = CENTER_DATA;
    }

    public ParameterPLS(int factors) {
        super(ConstantsRegressionMethods.MODEL_PLS);

        setFactors(factors);

        centerData = CENTER_DATA;
    }

    public ParameterPLS(ParameterPLS p) {
        super(p);
        copy(p);
    }

    public void copy(ParameterPLS orig){
        super.copy(orig);
        setFactors(orig.getFactors());
        setCenterData(orig.isCenterData());
    }

    @Override
    public boolean equals(Object obj) {

        if(!(obj instanceof ParameterPLS)){
            return false;
        }

        boolean eq = true;

        ParameterPLS p = (ParameterPLS)obj;

        if(!getName().equals(p.getName())){
            eq = false;
        } if(getFactors() != p.getFactors()){
            eq = false;
        }

        return eq;
    }

    /**
     * For inheritance of BxCox power transformation.
     * @param name
     * @param factors
     */
    protected ParameterPLS(String name, int factors) {
        super(name);

        setFactors(factors);
        centerData = CENTER_DATA;
    }

    public int getFactors() {
        return factors;
    }

    public void setFactors(int factors) {
        this.factors = factors;
        properties.put(TAG_FACTORS, Integer.toString(factors));
    }

    public boolean isCenterData() {
        return centerData;
    }

    public void setCenterData(boolean centerData) {
        this.centerData = centerData;
    }

    @Override
    protected void decodeProperties2Parameter() {
        factors = Integer.parseInt(properties.getProperty(TAG_FACTORS));
    }

    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterPLS parameterPLS = (ParameterPLS)o;

        if(factors>parameterPLS.factors) {
            cmp=1;
        }else if(factors<parameterPLS.factors) {
            cmp=-1;
        }

        return cmp;
    }

    @Override
    public void read(File fiProperties)throws IOException {

        Properties properties = new Properties();

        properties.load(new FileReader(fiProperties));

        String name = properties.getProperty(TAG_NAME);

        if(!getName().equals(name)) {
            throw  new RuntimeException("Name confusion!");
        }

        String sFactors = properties.getProperty(TAG_FACTORS);

        factors = Integer.parseInt(sFactors);

        centerData = CENTER_DATA;

    }

    @Override
    public void write(File fiProperties) throws IOException {

        Properties properties = new Properties();

        properties.put(TAG_NAME, getName());

        properties.put(TAG_FACTORS, Integer.toString(factors));

        properties.store(new FileWriter(fiProperties), "");

    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ParameterPLS{");
        sb.append("name=").append(getName());
        sb.append(" factors=").append(factors);
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        li.add(TAG_FACTORS);

        return li;
    }

    public static void main(String[] args) throws IOException {

        File dir = new File("/home/korffmo1/tmp/tmp00");

        File fiProp = new File(dir, "pls.properties");

        ParameterPLS parameter = new ParameterPLS();

        parameter.factors = 3;

        parameter.write(fiProp);

        ParameterPLS parameterIn = new ParameterPLS();

        parameterIn.read(fiProp);

        System.out.println(parameterIn.toString());

    }

}
