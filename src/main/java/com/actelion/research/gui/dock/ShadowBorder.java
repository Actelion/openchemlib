package com.actelion.research.gui.dock;

import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.util.ColorHelper;

import javax.swing.border.AbstractBorder;
import java.awt.*;

public class ShadowBorder extends AbstractBorder {
    private static final long serialVersionUID = 0x20070807;

    private static int sBorderWidth,sThinSize;
    private Insets mInsets;
    private Color mColor;

    /**
     * Creates a border with a grey 1 pixel line on left and top
     * and a 3 pixel (scaled to Hi-dpi factor) wide shadow on right and bottom.
     * Insets are as large as needed for the drawn border.
     */
    public ShadowBorder() {
    	mInsets = new Insets(thinSize(), thinSize(), shadowSize(), shadowSize());
    	}

	public static int thinSize() {
		if (sThinSize == 0)
			sThinSize = HiDPIHelper.scale(1);

		return sThinSize;
		}

	public static int shadowSize() {
    	if (sBorderWidth == 0)
    		sBorderWidth = HiDPIHelper.scale(3);

    	return sBorderWidth;
		}

	public void setColor(Color c) {
    	mColor = c;
		}

    /**
     * Creates a border with a grey 1 pixel line on left and top
     * and a 3 pixel wide shadow on right and bottom.
     * Insets can be specified to be larger than the drawn border,
     * which allows to adjust the spacing outside the drawn border.
     */
    public ShadowBorder(int top, int left, int bottom, int right) {
    	mInsets = new Insets(HiDPIHelper.scale(top), HiDPIHelper.scale(left), HiDPIHelper.scale(bottom), HiDPIHelper.scale(right));
    	}

    public Insets getBorderInsets(Component c) {
        return mInsets;
        }

    @Override
    public void paintBorder(Component c, Graphics g, int x, int y, int w, int h) {
    	super.paintBorder(c, g, x, y, w, h);

    	Color[] shadow = new Color[shadowSize()];
        shadow[0] = mColor != null ? mColor : c.getBackground().darker();
        Color parentBG = c.getParent().getBackground();
        for (int i=1; i<shadow.length; i++)
	        shadow[i] = ColorHelper.intermediateColor(shadow[0], parentBG, 0.3f+0.7f*(float)Math.pow((float)i/shadow.length, 1.6));

        g.translate(x, y);

        int t = mInsets.top;
        int l = mInsets.left;
        h -= mInsets.bottom;
        w -= mInsets.right;
        int b = h;
        int r = w;
        h -= mInsets.top;
        w -= mInsets.left;

		g.setColor(parentBG);
		g.fillRect(r+1, t-1, shadow.length-1, shadow.length);
		g.fillRect(l-1, b+1, shadow.length, shadow.length-1);
//		g.fillRect(r+shadow.length-1, t+shadow.length-2, 1, 1);
//		g.fillRect(l+shadow.length-2, b+shadow.length-1, 1, 1);
		g.fillRect(r+shadow.length-1, b+shadow.length-1, 1, 1);

		g.setColor(shadow[0]);
		g.fillRect(l-1, t-1, w+1, 1);
		g.fillRect(l-1, t-1, 1, h+1);
		g.fillRect(r, t, 1, h);
		g.fillRect(l, b, w, 1);

		for (int i=1; i<shadow.length; i++) {
			g.setColor(shadow[i]);
			g.fillRect(r+i, t+i, 1, h);
			g.fillRect(l+i, b+i, w, 1);
			g.fillRect(r+i-1, t+i-2, 1, 1);
			g.fillRect(l+i-2, b+i-1, 1, 1);
			g.fillRect(r+i-1, b+i-1, 1, 1);
			}

        g.translate(-x, -y);
        }
    }
