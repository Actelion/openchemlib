package com.actelion.research.share.gui;

/**
 * User: rufenec
 * Creation Date: 8/24/2016
 */
public abstract class DrawConfig
{
    public final long createColor(double r, double g, double b, double alpha)
    {
        long col = (long)(r * 255) << 24 | (long)(g * 255) << 16 | (long)(b * 255) << 8 | (long)(alpha*255);
        return col;
    }

    public abstract long getHighLightColor();
    public abstract long getMapToolColor();
    public abstract long getSelectionColor();
    public abstract long getForegroundColor();
    public abstract long getBackgroundColor();
}
