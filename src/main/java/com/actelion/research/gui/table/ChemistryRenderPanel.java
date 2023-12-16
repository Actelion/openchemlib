package com.actelion.research.gui.table;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.Depictor2D;
import com.actelion.research.chem.ExtendedDepictor;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.LookAndFeelHelper;
import com.actelion.research.gui.generic.GenericDepictor;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import com.actelion.research.gui.swing.SwingDrawContext;
import com.actelion.research.util.ColorHelper;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

public class ChemistryRenderPanel extends JPanel {
    static final long serialVersionUID = 0x20070312;

    private Object  mChemistry;
    private int mDisplayMode,mOverruleForegroundARGB;
    private double mTextSizeFactor;
    private boolean mAlternateBackground;

    public ChemistryRenderPanel() {
        mOverruleForegroundARGB = 0;
        mDisplayMode = 0;
        mTextSizeFactor = 1.0;
        }

    public void setChemistry(Object chemistry) {
        mChemistry = chemistry;
        repaint();
        }

    public void update(Graphics g) {
        paint(g);
        }

    public void setSelected(boolean isSelected) {
        if (isSelected) {
            setForeground(UIManager.getColor("Table.selectionForeground"));
            setBackground(UIManager.getColor("Table.selectionBackground"));
            }
        else {
            setForeground(UIManager.getColor("Table.foreground"));
            Color bg = UIManager.getColor("Table.background");
            setBackground(!mAlternateBackground ? bg : ColorHelper.darker(bg, 0.94f));
            }
        }

    public void setFocus(boolean hasFocus) {
        if (hasFocus)
            setBorder( UIManager.getBorder("Table.focusCellHighlightBorder") );
        else
            setBorder(new EmptyBorder(1, 1, 1, 1));
        }

    @Deprecated
    public void setOverruleForeground(Color fg) {
        mOverruleForegroundARGB = fg == null ? 0 : fg.getRGB();
    	}

    public void setOverruleForeground(int argb) {
        mOverruleForegroundARGB = argb;
    }

    public void setAlternateBackground(boolean b) {
        mAlternateBackground = b;
        }

    public void setDisplayMode(int mode) {
        mDisplayMode = mode;
        }

    public void setTextSizeFactor(double factor) {
        mTextSizeFactor = factor;
        }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        ((Graphics2D)g).setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

        Rectangle r = new Rectangle(new java.awt.Point(0,0), getSize());

        // Substance Graphite LaF does not consider the defined background
        if (LookAndFeelHelper.isNewSubstance() || LookAndFeelHelper.isRadiance()) {
            g.setColor(getBackground());
            ((Graphics2D) g).fill(r);
            }

        r.grow(-2, -2);

        Insets insets = getInsets();
        r.x += insets.left;
        r.y += insets.top;
        r.width -= insets.left + insets.right;
        r.height -= insets.top + insets.bottom;

        if (mChemistry != null && r.width > 0 && r.height > 0) {
            if (mChemistry instanceof StereoMolecule) {
                GenericDepictor d = new GenericDepictor((StereoMolecule)mChemistry, mDisplayMode | Depictor2D.cDModeSuppressChiralText);
                d.setFactorTextSize(mTextSizeFactor);
                d.setForegroundColor(getForeground().getRGB(), getBackground().getRGB());
                if (mOverruleForegroundARGB != 0)
                	d.setOverruleColor(mOverruleForegroundARGB, getBackground().getRGB());
                int avbl = HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen);
                SwingDrawContext context = new SwingDrawContext((Graphics2D)g);
                d.validateView(context, new GenericRectangle(r.x, r.y, r.width, r.height), AbstractDepictor.cModeInflateToMaxAVBL | avbl);
                d.paint(context);
                }
            if (mChemistry instanceof Reaction) {
            	Reaction rxn = (Reaction)mChemistry;
                ExtendedDepictor d = new ExtendedDepictor(rxn, rxn.getDrawingObjects(), rxn.isReactionLayoutRequired());
                d.setFactorTextSize(mTextSizeFactor);
                d.setDisplayMode(mDisplayMode | Depictor2D.cDModeSuppressChiralText);
                d.setForegroundColor(getForeground().getRGB(), getBackground().getRGB());
                if (mOverruleForegroundARGB != 0)
                	d.setOverruleColor(mOverruleForegroundARGB, getBackground().getRGB());
                int avbl = HiDPIHelper.scale(AbstractDepictor.cOptAvBondLen);
                GenericDrawContext context = new SwingDrawContext((Graphics2D)g);
                d.validateView(context, new GenericRectangle(r.x, r.y, r.width, r.height), AbstractDepictor.cModeInflateToMaxAVBL | avbl);
                d.paint(context);
                }
            }
        }
    }
