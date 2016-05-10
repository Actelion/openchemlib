/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.chem;

import com.actelion.research.chem.reaction.Reaction;

import java.awt.*;
import java.awt.geom.Rectangle2D;


public class ChemistryHelper
{
    public static final int REACTION_TYPE_NOMOLS = 0;
    public static final int REACTION_TYPE_NOPRODUCTS = 1;
    public static final int REACTION_TYPE_REACTANTS = 2;
    public static final int REACTION_TYPE_NORMAL = 3;

    private ChemistryHelper()
    {
    }

    public static int getReactionType(Reaction r)
    {
        int mols = r.getMolecules();
        int nr = r.getReactants();
        int np = r.getProducts();
        if (mols == 0)
            return REACTION_TYPE_NOMOLS;
        else if (nr == 0)
            return REACTION_TYPE_REACTANTS;
        else if (np == 0)
            return REACTION_TYPE_NOPRODUCTS;
        else
            return REACTION_TYPE_NORMAL;
    }

    public static Rectangle2D.Double getBoundingRect(Reaction r, boolean includearrows)
    {
        Rectangle2D.Double res = null;
        if (r == null)
            return null;
        int nr = r.getReactants();
        int np = r.getProducts();
        int m = r.getMolecules();
        double width = 0;
        double avgwidth = 0;

        for (int i = 0; i < m; i++) {
            ExtendedMolecule mol = r.getMolecule(i);
            int nb = mol.getAllBonds();
            Rectangle2D.Double d = getBoundingRect(mol);
//			System.out.println("getBoundingRect " + i + "=" + d);
            if (res == null) {
                if (d != null && nb > 0) {
                    res = d;
                    width += d.getWidth();
                } else {
//					System.out.println("Empty bounding rect for molecule " + i);
                }
            } else {
                if (d != null && nb > 0) {
                    res = (Rectangle2D.Double) res.createUnion(d);
                    width += d.getWidth();
                } else {
//					System.out.println("Empty bounding rect for molecule " + i);
                }
            }
        }
        if (includearrows && nr == 0 && m > 0) {
            avgwidth = width / m;
            res = new Rectangle2D.Double(res.getX() - avgwidth, res.getY(), res.getWidth() + avgwidth, res.getHeight());
        } else if (includearrows && np == 0 && m > 0) {
            avgwidth = width / m;
            res = new Rectangle2D.Double(res.getX(), res.getY(), res.getWidth() + avgwidth, res.getHeight());
        }
//		System.out.println("Reaction av bnd length: " + width + " getBoundingBox " + res);

        return res;
    }

    public static Rectangle2D.Double getBoundingRect(ExtendedMolecule m)
    {
        double xmax = Double.MIN_VALUE;
        double ymax = Double.MIN_VALUE;
        double xmin = Double.MAX_VALUE;
        double ymin = Double.MAX_VALUE;

        if (m == null)
            return null;

        int na = m.getAllAtoms();
        double bl = m.getAverageBondLength();
        for (int i = 0; i < na; i++) {
            xmax = Math.max(xmax, m.getAtomX(i));
            xmin = Math.min(xmin, m.getAtomX(i));
            ymax = Math.max(ymax, m.getAtomY(i));
            ymin = Math.min(ymin, m.getAtomY(i));
        }

        return (na > 0) ? new Rectangle2D.Double(xmin, ymin, Math.max(xmax - xmin, bl), Math.max(ymax - ymin, bl)) : null;
    }


    public static double getAverageBondLength(Reaction r)
    {
        int rn = r.getMolecules();
        double avg = 0;
        double t = 0;
        for (int i = 0; i < rn; i++) {
            t += r.getMolecule(i).getAverageBondLength();
        }
        if (rn > 0)
            return t / rn;
        else
            return 0;
    }

    public static Rectangle2D.Double getReactantsBoundingRect(Reaction r)
    {
        if (r == null)
            throw new NullPointerException("Cannot pass null reaction");
        int rn = r.getReactants();
        Rectangle2D.Double rc = null;
        for (int i = 0; i < rn; i++) {
            Rectangle2D.Double rt = getBoundingRect(r.getReactant(i));
            if (rc != null && rt != null) {
                Rectangle2D tmp = rc.createUnion(rt);
                rc = new Rectangle2D.Double(tmp.getX(), tmp.getY(), tmp.getWidth(), tmp.getHeight());
            } else
                rc = rt;
        }
        return rc;
    }


    public static Rectangle2D.Double getArrowBoundingRect(Reaction r)
    {
        Rectangle2D rr = getReactantsBoundingRect(r);
        Rectangle2D rp = getProductsBoundingRect(r);
        if (rr != null && rp != null) {
            Rectangle2D union = rr.createUnion(rp);
            double y = union.getHeight() / 2 + union.getY();
            double rx = rr.getMaxX();
            double px = rp.getMinX();
            return new Rectangle2D.Double(rx < px ? rx : px, y, Math.abs(px - rx), 0);
        } else if (rr != null) {
            double y = rr.getHeight() / 2 + rr.getY();
            double rx = rr.getMaxX();
            double width = rr.getWidth();
            return new Rectangle2D.Double(rx, y, width, 0);

        } else if (rp != null) {
            double y = rp.getHeight() / 2 + rp.getY();
            double width = 0;
            double rx = rp.getMinX();
            return new Rectangle2D.Double(rx, y, width, 0);
        } else {
            return new Rectangle2D.Double(0, 0, 0, 0);
        }
    }

    public static Rectangle2D.Double getDiffRect(Rectangle2D rr, Rectangle2D rp)
    {
        if (rr != null && rp != null) {
            Rectangle2D union = rr.createUnion(rp);
            double y = union.getMinY();
            double rx = rr.getMaxX();
            double px = rp.getMinX();
            double ry = union.getMinY();
            double py = union.getMaxY();
            Rectangle2D.Double res = new Rectangle2D.Double(
                rx < px ? rx : px,
                ry < py ? ry : py,
                Math.abs(px - rx),
                Math.abs(py - ry));
//            System.out.println("getDiffRect (1)" + rr);
//            System.out.println("getDiffRect (2)" + rp);
//            System.out.println("getDiffRect (out)" + res);
            return res;
        } else {
            return null;
        }
    }

    private static Rectangle2D.Double scaleTo(Rectangle2D.Double src, double width, double height)
    {
        double x = (src.getWidth() - width) / 2 + src.getX();
        double y = (src.getHeight() - height) / 2 + src.getY();
        return new Rectangle2D.Double(x, y, width, height);
    }

    public static Rectangle2D getProductsBoundingRect(Reaction r)
    {
        if (r == null)
            throw new NullPointerException("Cannot pass null reaction");
        int rn = r.getProducts();
        Rectangle2D rc = null;
        for (int i = 0; i < rn; i++) {
            Rectangle2D rt = getBoundingRect(r.getProduct(i));
            if (rc != null && rt != null)
                rc = rc.createUnion(rt);
            else
                rc = rt;
        }
        return rc;
    }


    public static void transformReaction(Reaction r, double offsetx, double offsety, double scale)
    {
        int mols = r.getMolecules();
        for (int i = 0; i < mols; i++) {
            transformMolecule(r.getMolecule(i), offsetx, offsety, scale);
        }
    }

    public static void transformMolecules(Molecule[] mols, double offsetx, double offsety, double scale)
    {
        if (mols != null) {
            for (Molecule m : mols) {
                transformMolecule(m,offsetx,offsety,scale);
            }
        }
    }
    public static void transformMolecule(Molecule m, double offsetx, double offsety, double scale)
    {
        int atoms = m.getAllAtoms();
        for (int i = 0; i < atoms; i++) {
            m.setAtomX(i, (m.getAtomX(i) + offsetx) * scale);
            m.setAtomY(i, (m.getAtomY(i) + offsety) * scale);
        }
    }

    public static void scaleInto(Reaction reaction, double x, double y, double w, double h)
    {
        Rectangle2D bounds = getBoundingRect(reaction, false);
        if (bounds == null)
            return;

        double scaleHorizontal = w / bounds.getWidth();
        double scaleVertical = h / bounds.getHeight();
        double scale;

        if (scaleHorizontal < scaleVertical) {   // scale on x-dimension
            scale = w / (float) bounds.getWidth();
        } else {                                // scale on y-dimension
            scale = h / (float) bounds.getHeight();
        }
        transformReaction(reaction, 0, 0, scale);
        bounds = getBoundingRect(reaction, false);
        double dx = x - bounds.getMinX() + (w - bounds.getWidth()) / 2;
        double dy = y - bounds.getMinY() + (h - bounds.getHeight()) / 2;
        transformReaction(reaction, dx, dy, 1);

    }

    public static void scaleInto(Reaction reaction, double x, double y, double width, double height, double arrowSize)
    {

        if (width > arrowSize * 2) {
            ExtendedMolecule[] reactants = getReactants(reaction);
            ExtendedMolecule[] products = getProducts(reaction);

            Rectangle2D rb = getBoundingRect(reactants);
            Rectangle2D pb = getBoundingRect(products);

            // left and right space
            double w = (width - arrowSize) / 2;
            double h = height;
            double scaleHorizontal = w / Math.max(rb.getWidth(), pb.getWidth());
            double scaleVertical = h / Math.max(rb.getHeight(), pb.getHeight());
            //        System.out.printf("Scaling %f vs %f\n", scaleHorizontal, scaleVertical);

            double scale;
            if (scaleHorizontal < scaleVertical) {   // scale on x-dimension
                scale = w / Math.max((float) rb.getWidth(), (float) pb.getWidth()) ;
                //            System.out.printf("Scaling horiz %f\n", scale);
            } else {                                // scale on y-dimension
                scale = h / Math.max((float) rb.getHeight(), (float) pb.getHeight()) ;
                //            System.out.printf("Scaling vert %f\n", scale);
            }

            double avbl = getAverageBondLength(reaction);
//            double sc = AbstractDepictor.cOptAvBondLen/avbl;
            scale = Math.min(AbstractDepictor.cOptAvBondLen/avbl,scale);
//            if (reactants.length <= 1 && products.length <= 1)
//                scale /= 4;
            {

                transformMolecules(reactants, 0, 0, scale);
                rb = getBoundingRect(reactants);
                double dx = x - rb.getMinX() + (w - rb.getWidth()) / 2;
                double dy = y - rb.getMinY() + (h - rb.getHeight()) / 2;
                transformMolecules(reactants, dx, dy, 1);
            }

            {
                transformMolecules(products, 0, 0, scale);
                pb = getBoundingRect(products);
                double dx = x + w + arrowSize - pb.getMinX() + (w - pb.getWidth()) / 2;
                double dy = y - pb.getMinY() + (h - pb.getHeight()) / 2;
                transformMolecules(products, dx, dy, 1);
            }

        }
    }


    public static ExtendedMolecule[] getReactants(Reaction r)
    {
        ExtendedMolecule[] ret = new ExtendedMolecule[r.getReactants()];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = r.getReactant(i);
        }
        return ret;
    }

    public static ExtendedMolecule[] getProducts(Reaction r)
    {
        ExtendedMolecule[] ret = new ExtendedMolecule[r.getProducts()];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = r.getProduct(i);
        }
        return ret;
    }

    public static Rectangle2D getBoundingRect(ExtendedMolecule[] mols)
    {
        if (mols == null || mols.length == 0) {
            return new Rectangle2D.Double(0,0,0,0);
        }

        Rectangle2D r = getBoundingRect(mols[0]);
        for (int i = 1; i < mols.length; i++) {
            Rectangle2D t = getBoundingRect(mols[i]);
            if (t != null) {
                if (r == null)
                    r = t;
                else
                    r = r.createUnion(getBoundingRect(mols[i]));
            }
        }
        if (r == null)
            r = new Rectangle2D.Double(0,0,0,0);
        return r;
    }

}
