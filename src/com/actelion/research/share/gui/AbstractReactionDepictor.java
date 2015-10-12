/*
 * Project: DD_jfx
 * @(#)AbstractReactionDepictor.java
 *
 * Copyright (c) 1997- 2015
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * Author: Christian Rufener
 */



package com.actelion.research.share.gui;


//import java.awt.*;
import java.awt.geom.*;

import com.actelion.research.chem.*;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.share.gui.editor.geom.IDrawContext;

public abstract class AbstractReactionDepictor
{
    private Reaction rxn_ = null;
    public static final int ARROWWIDTH = 40;
    public static final int ARROWHEIGHT = 20;
    public static final int PLUSSIZE = 20;
    private static final String PLUSSTRING = "+";
    private String PLUSFONT = "Helvetica";
    private int PLUSFONTSIZE = 12;
    private static Rectangle2D dim_ = new Rectangle2D.Double(0, 0, 600, 400);
    private DrawingObjectList drwobj_ = null;
    private int displaymode = 0;
    private java.awt.Color foreGround = null;

/*
    public AbstractReactionDepictor()
    {
    }
*/


    public AbstractReactionDepictor(Reaction r, DrawingObjectList drwobj)
    {
        rxn_ = r;
        drwobj_ = drwobj;
    }

    public AbstractReactionDepictor(Reaction r)
    {
        this(r, null);
    }

    public void setReaction(Reaction r)
    {
        rxn_ = r;
    }

    public void setDisplayMode(int mode)
    {
        displaymode = mode;
    }

    public void setDrawingObjects(DrawingObjectList drwobj)
    {
        drwobj_ = drwobj;
    }

//    private boolean useExtendedDepictor = true;

    public DepictorTransformation paint(IDrawContext g)
    {
        if (rxn_ != null)
            drawReaction(g, rxn_, drwobj_);
        return null;
    }
/*
    private ReactionArrow getReactionArrow()
    {
        java.awt.geom.Rectangle2D ar = ChemistryHelper.getArrowBoundingRect(rxn_);
        if (ar != null) {
            ReactionArrow arrow = new ReactionArrow();
            arrow.setCoordinates(ar.getX(),ar.getY(),ar.getX()+ar.getWidth(),ar.getY()+ ar.getHeight());
            return arrow;
        }
        return null;
    }


    private Dimension getSize()
    {
        return new Dimension(dim_.width,dim_.height);
    }
*/

    private void drawReaction(IDrawContext ctx, Reaction rxn, DrawingObjectList dobs)
    {
        java.awt.geom.Rectangle2D rxnBounds = ChemistryHelper.getBoundingRect(rxn, true);
        int numProd = rxn.getProducts();
        int numReact = rxn.getReactants();
        int mols = rxn.getMolecules();
        if (rxnBounds != null) {

            double offsetx = -rxnBounds.getX();
            double offsety = -rxnBounds.getY();
            double scaler = rxnBounds.getHeight() / rxnBounds.getWidth();
            Rectangle2D d = dim_;
            double insetx = dim_.getMinX() + 4;
            double insety = dim_.getMinX() + 4;
            double insetwidth = dim_.getWidth() - 8;
            double insetheight = dim_.getHeight() - 8;

            double scalec = (double) d.getHeight() / (double) d.getWidth();
            double scale = 1;

            if (scalec > scaler) {                   // scale on x-dimension
                scale = insetwidth / rxnBounds.getWidth();
                double t = (insetheight / scale - rxnBounds.getHeight()) / 2;
                offsety += t;
            } else {                                // scale on y-dimension
                scale = insetheight / rxnBounds.getHeight();
                double t = (insetwidth / scale - rxnBounds.getWidth()) / 2;
                offsetx += t;
            }
            System.out.println("Transforming reaction " + offsetx + " " + offsety + " " + scale);
            ChemistryHelper.transformReaction(rxn, offsetx, offsety, scale);
            java.awt.geom.Rectangle2D rn = ChemistryHelper.getBoundingRect(rxn, false);
            System.out.println("Transformed reaction " + rn);

            for (int i = 0; i < mols; i++) {
                AbstractDepictor depict = createDepictor(rxn.getMolecule(i), displaymode);
//                if (foreGround != null)
//                    depict.setOverruleColor(foreGround,background);
                java.awt.geom.Rectangle2D o = ChemistryHelper.getBoundingRect(rxn.getMolecule(i));
                if (o != null) {
                    DepictorTransformation tm = depict.validateView(ctx.getNative(),
                            new java.awt.geom.Rectangle2D.Float((float) o.getX() + (float) insetx, (float) o.getY() + (float) insety, (float) o.getWidth(), (float) o.getHeight()),
                            AbstractDepictor.cModeInflateToMaxAVBL);
                    depict.paint(ctx.getNative());
                }
            }
            java.awt.geom.Rectangle2D ar = ChemistryHelper.getArrowBoundingRect(rxn);
            if (ar != null) {
                Arrow a = null;
                if (ar.getWidth() == 0) {
                    a = new Arrow((int) insetx + (int) ar.getX() - ARROWWIDTH, (int) insety + (int) ar.getY(), ARROWWIDTH, ARROWHEIGHT);
                } else if (ar.getWidth() < ARROWWIDTH) {
                    a = new Arrow((int) insetx + (int) ar.getX(), (int) insety + (int) ar.getY(), (int) ar.getWidth(), ARROWHEIGHT);
                } else {
                    a = new Arrow((int) insetx + (int) ar.getX() + (int) (ar.getWidth() - ARROWWIDTH) / 2, (int) insety + (int) ar.getY(), ARROWWIDTH, ARROWHEIGHT);
                }
                a.paint(ctx);
            }
//            FontMetrics fm = g.getFontMetrics(PLUSFONT);
//            java.awt.font.LineMetrics lm = fm.getLineMetrics(PLUSSTRING, g);

            for (int i = 1; i < numReact; i++) {
                java.awt.geom.Rectangle2D r1 = ChemistryHelper.getBoundingRect(rxn.getMolecule(i - 1));
                java.awt.geom.Rectangle2D r2 = ChemistryHelper.getBoundingRect(rxn.getMolecule(i));
                java.awt.geom.Rectangle2D rp = null;
                if (r1.intersects(r2.getX(), r2.getY(), r2.getWidth(), r2.getHeight())) {
                    rp = r1.createIntersection(r2);
                } else
                    rp = ChemistryHelper.getDiffRect(r1, r2);
                int x = (int) rp.getCenterX();
                int y = (int) rp.getCenterY();
                ctx.setFont(PLUSFONT,PLUSFONTSIZE,true);
                ctx.drawText(PLUSSTRING, (int) insetx + x, (int) insety + y,true,true);
            }

            for (int i = 1; i < numProd; i++) {
                java.awt.geom.Rectangle2D r1 = ChemistryHelper.getBoundingRect(rxn.getMolecule(numReact + i - 1));
                java.awt.geom.Rectangle2D r2 = ChemistryHelper.getBoundingRect(rxn.getMolecule(numReact + i));
                java.awt.geom.Rectangle2D rp = null;
                if (r1.intersects(r2.getX(), r2.getY(), r2.getWidth(), r2.getHeight())) {
                    rp = r1.createIntersection(r2);
                } else
                    rp = ChemistryHelper.getDiffRect(r1, r2);
                int x = (int) rp.getCenterX();
                int y = (int) rp.getCenterY();
                ctx.setFont(PLUSFONT,PLUSFONTSIZE,true);
                ctx.drawText(PLUSSTRING, (int) insetx + x, (int) insety + y,true, true);
            }

    /*
            if (dobs != null) {
                int size = dobs.size();
                for (int i = 0; i < size; i ++) {
                    AbstractDrawingObject o = (AbstractDrawingObject)dobs.get(i);
                    java.awt.geom.Rectangle2D r1 = o.getBoundingRect();
                    System.out.println("Drawing Object " + o + "\nBounds: " + r1 + "\n Offsets (" + offsetx +","+offsety + ")");
                    DepictorTransformation tm = new DepictorTransformation();
         //           tm.move(offsetx*scale,offsety*scale);
          //          tm.setScaling(scale);
                    o.draw(g,tm);
                }
            }
     */
        }
    }

    public abstract AbstractDepictor createDepictor(StereoMolecule molecule, int displaymode);

    public boolean updateCoords(IDrawContext g, double x1, double y1, double width, double height, int mode)
    {
        ExtendedDepictor dep = new ExtendedDepictor(rxn_, null,false, false);
        dep.updateCoords(g.getNative(), new java.awt.geom.Rectangle2D.Float((float)x1, (float)y1, (float)width, (float)height), AbstractDepictor.cModeInflateToMaxAVBL);

        dim_ = new Rectangle2D.Double((int) x1, (int) y1, (int) width, (int) height);
        return true;
    }

    private DrawingObjectList getPlusOjects(Reaction rxn)
    {
        DrawingObjectList res = new DrawingObjectList();
        if (rxn != null) {
            int mols = rxn.getMolecules();
            int rn = rxn.getReactants();
            int pn = rxn.getProducts();

            for (int i = 1; i < mols; i++) {
                if (i != rn) {
                    ExtendedMolecule m = rxn.getMolecule(i - 1);
                    java.awt.geom.Rectangle2D rc1 = ChemistryHelper.getBoundingRect(m);
                    m = rxn.getMolecule(i);
                    java.awt.geom.Rectangle2D rc2 = ChemistryHelper.getBoundingRect(m);
                    java.awt.geom.Rectangle2D.Double d = ChemistryHelper.getDiffRect(rc1, rc2);
                    if (d != null) {
                        TextDrawingObject o = new TextDrawingObject();
                        o.setCoordinates((float)d.getCenterX(), (float)d.getCenterY());
                        o.setValues("+", 10, TextDrawingObject.DEFAULT_STYLE);
                        res.add(o);
                    }
                }
            }
        }
        return res;
    }

    public void setColor(java.awt.Color c)
    {
        foreGround = c;
    }
}
