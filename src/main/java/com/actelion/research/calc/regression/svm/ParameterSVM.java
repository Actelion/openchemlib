package com.actelion.research.calc.regression.svm;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.regression.ConstantsRegressionMethods;
import com.actelion.research.calc.regression.ParameterRegressionMethod;
import org.machinelearning.svm.libsvm.svm_parameter;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;

/**
 * ParameterSVM
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 06.12.18.
 */
public class ParameterSVM extends ParameterRegressionMethod {

    public static final String TAG_SVM_TYPE="SVMRegressionType";
    public static final String TAG_KERNEL="Kernel";
    public static final String TAG_DEGREE="Degree";
    public static final String TAG_GAMMA="Gamma";
    public static final String TAG_COEF0="Coef0";
    public static final String TAG_CACHE_SIZE="CacheSize";
    public static final String TAG_EPSILON="Epsilon";
    public static final String TAG_C="C";
    public static final String TAG_NR_WEIGHT="nrWeight";
    public static final String TAG_NU="Nu";
    public static final String TAG_P="p";
    public static final String TAG_SHRINKING="Shrinking";
    public static final String TAG_PROBABILITY="Probability";


    private int svmType;
    private int kernelType;
    private int degree;	// for poly
    private double gamma;	// for poly/rbf/sigmoid
    private double coef0;	// for poly/sigmoid

    // these are for training only
    private double cache_size; // in MB
    private double eps;	// stopping criteria
    private double C;	// for C_SVC, EPSILON_SVR and NU_SVR
    private int nr_weight;		// for C_SVC
    private int[] weight_label;	// for C_SVC
    private double[] weight;		// for C_SVC
    private double nu;	// for NU_SVC, ONE_CLASS, and NU_SVR
    private double p;	// for EPSILON_SVR
    private int shrinking;	// use the shrinking heuristics
    private int probability; // do probability estimates





    public ParameterSVM() {
        this(SVMParameterHelper.regressionEpsilonSVR());
        initializeProperties();
    }

    public ParameterSVM(double nu) {
        this(SVMParameterHelper.standard());
        setNu(nu);
    }

    public ParameterSVM(ParameterSVM parameterSVM) {
        super(ConstantsRegressionMethods.MODEL_SVM);
        this.svmType = parameterSVM.svmType;
        this.kernelType = parameterSVM.kernelType;
        this.degree = parameterSVM.degree;
        this.gamma = parameterSVM.gamma;
        this.coef0 = parameterSVM.coef0;
        this.cache_size = parameterSVM.cache_size;
        this.eps = parameterSVM.eps;
        this.C = parameterSVM.C;
        this.nr_weight = parameterSVM.nr_weight;
        this.weight_label = parameterSVM.weight_label;
        this.weight = parameterSVM.weight;
        this.nu = parameterSVM.nu;
        this.p = parameterSVM.p;
        this.shrinking = parameterSVM.shrinking;
        this.probability = parameterSVM.probability;
        initializeProperties();
    }

    public ParameterSVM(svm_parameter svmParameter) {
        super(ConstantsRegressionMethods.MODEL_SVM);
        this.svmType = svmParameter.svm_type;
        this.kernelType = svmParameter.kernel_type;
        this.degree = svmParameter.degree;
        this.gamma = svmParameter.gamma;
        this.coef0 = svmParameter.coef0;
        this.cache_size = svmParameter.cache_size;
        this.eps = svmParameter.eps;
        this.C = svmParameter.C;
        this.nr_weight = svmParameter.nr_weight;
        this.weight_label = svmParameter.weight_label;
        this.weight = svmParameter.weight;
        this.nu = svmParameter.nu;
        this.p = svmParameter.p;
        this.shrinking = svmParameter.shrinking;
        this.probability = svmParameter.probability;
        initializeProperties();
    }

    private void initializeProperties(){

        if(svmType != svm_parameter.EPSILON_SVR) {
            System.out.println("Error: only Epsilon support vector regression possible!");
            System.out.println("Error: opsilon support vector regression possible!");
            System.out.println("Error: opsilon support vector regression possible!");
            throw new RuntimeException("Only Epsilon support vector regression possible!");
        }

        setSVMRegressionType(svmType);
        setKernelType(kernelType);

        setEpsilon(eps);

        setC(C);

        setGamma(gamma);
    }

    public svm_parameter getSvmParameter(){

        svm_parameter svmParameter = new svm_parameter();

        svmParameter.svm_type = this.svmType;
        svmParameter.kernel_type = this.kernelType;
        svmParameter.degree = this.degree;
        svmParameter.gamma = this.gamma;
        svmParameter.coef0 = this.coef0;
        svmParameter.cache_size = this.cache_size;
        svmParameter.eps = this.eps;
        svmParameter.C = this.C;
        svmParameter.nr_weight = this.nr_weight;
        svmParameter.weight_label = this.weight_label;
        svmParameter.weight = this.weight;
        svmParameter.nu = this.nu;
        svmParameter.p = this.p;
        svmParameter.shrinking = this.shrinking;
        svmParameter.probability = this.probability;

        return svmParameter;
    }

    protected ParameterSVM(String nameRegressionMethod) {
        super(nameRegressionMethod);
    }

    public double getGamma(){
        return gamma;
    }

    public double getNu(){
        return nu;
    }

    public void setGamma(double gamma){
        this.gamma = gamma;
        properties.put(TAG_GAMMA, Double.toString(gamma));
    }

    public void setNu(double nu){
        this.nu = nu;
        properties.put(TAG_NU, Double.toString(nu));
    }

    public void setSVMRegressionType(int svmType) {
        if(svmType != svm_parameter.EPSILON_SVR) {
            throw new RuntimeException("Only Epsilon support vector regression possible!");
        }
        this.svmType = svmType;
        properties.put(TAG_SVM_TYPE, Integer.toString(svmType));
    }

    public void setEpsilon(double eps) {
        this.eps = eps;
        properties.put(TAG_EPSILON, Double.toString(eps));
    }

    public void setKernelType(int kernelType) {
        this.kernelType = kernelType;
        properties.put(TAG_KERNEL, Integer.toString(kernelType));
    }

    public void setC(double C) {
        this.C = C;
        properties.put(TAG_C, Double.toString(C));
    }

    @Override
    public int compareTo(ParameterRegressionMethod o) {

        int cmp = 0;

        ParameterSVM parameterSVM = (ParameterSVM)o;

        if(svmType > parameterSVM.svmType) {
            cmp=1;
        } else if(svmType < parameterSVM.svmType) {
            cmp=-1;
        }

        if(cmp==0) {

            if (kernelType > parameterSVM.kernelType) {
                cmp = 1;
            } else if (kernelType < parameterSVM.kernelType) {
                cmp = -1;
            }
        }

        if(cmp==0) {

            if (degree > parameterSVM.degree) {
                cmp = 1;
            } else if (degree < parameterSVM.degree) {
                cmp = -1;
            }
        }

        if(cmp==0) {

            if (gamma > parameterSVM.gamma) {
                cmp = 1;
            } else if (gamma < parameterSVM.gamma) {
                cmp = -1;
            }
        }

        if(cmp==0) {

            if (nu > parameterSVM.nu) {
                cmp = 1;
            } else if (nu < parameterSVM.nu) {
                cmp = -1;
            }
        }

        if(cmp==0) {

            if (C > parameterSVM.C) {
                cmp = 1;
            } else if (C < parameterSVM.C) {
                cmp = -1;
            }
        }

        if(cmp==0) {

            if (eps > parameterSVM.eps) {
                cmp = 1;
            } else if (eps < parameterSVM.eps) {
                cmp = -1;
            }
        }

        return cmp;
    }


    @Override
    public boolean equals(Object obj) {
        boolean eq = super.equals(obj);

        if(!eq){
            return false;
        }

        if(!(obj instanceof ParameterSVM)){
            return false;
        }

        ParameterSVM parameter = (ParameterSVM)obj;

        if(nr_weight!=parameter.nr_weight){
            return false;
        } else if(Math.abs(nr_weight-parameter.nr_weight) > Matrix.TINY08){
            return false;
        } else if(Math.abs(gamma -parameter.gamma) > Matrix.TINY08){
            return false;
        } else if(Math.abs(eps -parameter.eps) > Matrix.TINY08){
            return false;
        } else if(Math.abs(C -parameter.C) > Matrix.TINY08){
            return false;
        } else if(Math.abs(nu -parameter.nu) > Matrix.TINY08){
            return false;
        } else if(Math.abs(cache_size -parameter.cache_size) > Matrix.TINY08){
            return false;
        } else if(Math.abs(coef0 -parameter.coef0) > Matrix.TINY08){
            return false;
        } else if(Math.abs(this.p-parameter.p) > Matrix.TINY08){
            return false;
        } else if(svmType != parameter.svmType){
            return false;
        }else if(degree != parameter.degree){
            return false;
        }else if(kernelType != parameter.kernelType){
            return false;
        }else if(probability != parameter.probability){
            return false;
        }else if(shrinking != parameter.shrinking){
            return false;
        } else if(!Arrays.equals(weight, parameter.weight)){
            return false;
        } else if(!Arrays.equals(weight_label, parameter.weight_label)){
            return false;
        }

        return eq;
    }

    @Override
    protected void decodeProperties2Parameter() {
        svmType = Integer.parseInt(properties.getProperty(TAG_SVM_TYPE));

        kernelType = Integer.parseInt(properties.getProperty(TAG_KERNEL));

        C = Double.parseDouble(properties.getProperty(TAG_C));

        gamma = Double.parseDouble(properties.getProperty(TAG_GAMMA));

        eps = Double.parseDouble(properties.getProperty(TAG_EPSILON));
    }

    @Override
    public String toString() {

        DecimalFormat df = new DecimalFormat("0.0###");

        final StringBuilder sb = new StringBuilder("ParameterSVM{");
        sb.append("name=").append(getName());
        sb.append(" svm type=").append(SVMParameterHelper.getSVMType(svmType));
        sb.append(" kernel_type=").append(SVMParameterHelper.getKernelType(kernelType));
        sb.append(" degree=").append(degree);
        sb.append(" gamma=").append(df.format(gamma));
        sb.append(" coef0=").append(df.format(coef0));
        sb.append(" cache_size=").append(df.format(cache_size));
        sb.append(" eps=").append(df.format(eps));
        sb.append(" C=").append(df.format(C));
        sb.append(" nr_weight=").append(nr_weight);
        sb.append(" nu=").append(df.format(nu));
        sb.append(" p=").append(df.format(p));
        sb.append(" shrinking=").append(shrinking);
        sb.append(" probability=").append(probability);
        sb.append('}');
        return sb.toString();
    }

    public static List<String> getHeader(){

        List<String> li = ParameterRegressionMethod.getHeader();

        li.add(TAG_SVM_TYPE);
        li.add(TAG_KERNEL);
        li.add(TAG_DEGREE);
        li.add(TAG_GAMMA);
        li.add(TAG_COEF0);
        li.add(TAG_CACHE_SIZE);
        li.add(TAG_EPSILON);
        li.add(TAG_C);
        li.add(TAG_NR_WEIGHT);
        li.add(TAG_NU);
        li.add(TAG_P);
        li.add(TAG_SHRINKING);
        li.add(TAG_PROBABILITY);

        return li;
    }
}
