package com.actelion.research.calc.statistics;

import java.io.Serializable;

/**
 * ModelStatisticsOverview
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 29.11.18.
 */
public class ModelStatisticsOverview  implements Serializable {


    public double min;
    public double avr;
    public double max;

    public double sdv;

    public ModelStatisticsOverview(double min, double avr, double max, double sdv) {
        this.min = min;
        this.avr = avr;
        this.max = max;
        this.sdv = sdv;
    }

    public ModelStatisticsOverview() {
    }

    public boolean areAllFinite(){
        if(!Double.isFinite(min)) {
            return false;
        } else if(!Double.isFinite(avr)) {
            return false;
        } else if(!Double.isFinite(max)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ModelStatisticsOverview{");
        sb.append("min=").append(min);
        sb.append(", avr=").append(avr);
        sb.append(", max=").append(max);
        sb.append(", sdv=").append(sdv);
        sb.append('}');
        return sb.toString();
    }


}
