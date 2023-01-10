package com.actelion.research.calc.statistics;

import com.actelion.research.util.Formatter;

/**
 * ModelStatisticsOverviewMedian
 * <p>Modest v. Korff</p>
 * <p>
 * Created by korffmo1 on 27.03.19.
 */
public class ModelStatisticsOverviewMedian {


    public double percentile05;
    public double percentile25;
    public double median;
    public double percentile75;
    public double percentile95;

    public ModelStatisticsOverviewMedian(double percentile05, double percentile25, double median, double percentile75, double percentile95) {
        this.percentile05 = percentile05;
        this.percentile25 = percentile25;
        this.median = median;
        this.percentile75 = percentile75;
        this.percentile95 = percentile95;
    }

    public ModelStatisticsOverviewMedian() {
    }

    public boolean areAllFinite(){
        if(!Double.isFinite(percentile05)) {
            return false;
        } else if(!Double.isFinite(percentile25)) {
            return false;
        } else if(!Double.isFinite(median)) {
            return false;
        }else if(!Double.isFinite(percentile75)) {
            return false;
        }else if(!Double.isFinite(percentile95)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ModelStatisticsOverviewMedian{");
        sb.append("percentile05=").append(StatisticsOverview.DF4.format(percentile05));
        sb.append(", percentile25=").append(StatisticsOverview.DF4.format(percentile25));
        sb.append(", median=").append(StatisticsOverview.DF4.format(median));
        sb.append(", percentile75=").append(StatisticsOverview.DF4.format(percentile75));
        sb.append(", percentile95=").append(StatisticsOverview.DF4.format(percentile95));
        sb.append('}');
        return sb.toString();
    }
}
