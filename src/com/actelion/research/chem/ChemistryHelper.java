/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Christian Rufener
 */

package com.actelion.research.chem;

import com.actelion.research.chem.reaction.Reaction;

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
        double avgblen = getAverageBondLength(r);
        Rectangle2D rr = getReactantsBoundingRect(r);
        Rectangle2D rp = getProductsBoundingRect(r);
        if (rr != null && rp != null) {
            Rectangle2D union = rr.createUnion(rp);
            double y = union.getHeight() / 2 + union.getY();
            double rx = rr.getMaxX();
            double px = rp.getMinX();
//            return scaleTo(new Rectangle2D.Double(rx < px ? rx : px,y,Math.abs(px - rx),0),avgblen * 3,0);
            return new Rectangle2D.Double(rx < px ? rx : px, y, Math.abs(px - rx), 0);
        } else if (rr != null) {
            double y = rr.getHeight() / 2 + rr.getY();
            double rx = rr.getMaxX();
            double width = rr.getWidth();
            return new Rectangle2D.Double(rx, y, width, 0);

        } else if (rp != null) {
            double y = rp.getHeight() / 2 + rp.getY();
//            double width = rp.getWidth();
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

//        System.out.printf("Clean RXN Into %f %f %f %f\n",x,y,width,height);
        if (width > arrowSize * 2) {
            ExtendedMolecule[] reactants = getReactants(reaction);
            ExtendedMolecule[] products = getProducts(reaction);

            Rectangle2D rb = getBoundingRect(reactants);
            Rectangle2D pb = getBoundingRect(products);

//            System.out.printf("Reactant rectangle %s\n", rb);
//            System.out.printf("Product rectangle %s\n", pb);
            // left and right space
            double w = (width - arrowSize) / 2;
            double h = height;
            double scaleHorizontal = w / Math.max(rb.getWidth(), pb.getWidth());
            double scaleVertical = h / Math.max(rb.getHeight(), pb.getHeight());
            //        System.out.printf("Scaling %f vs %f\n", scaleHorizontal, scaleVertical);

            double scale;
            if (scaleHorizontal < scaleVertical) {   // scale on x-dimension
                scale = w / Math.max((float) rb.getWidth(), (float) pb.getWidth()) / 2;
                //            System.out.printf("Scaling horiz %f\n", scale);
            } else {                                // scale on y-dimension
                scale = h / Math.max((float) rb.getHeight(), (float) pb.getHeight()) / 2;
                //            System.out.printf("Scaling vert %f\n", scale);
            }


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

            //        for (int i = 0; i < reactants.length; i++) {
            //            ExtendedMolecule m = reactants[i];
            //            transformMolecule(m, 0, 0, scale);
            //            Rectangle2D bounds = getBoundingRect(m);
            //            double dx = x - bounds.getMinX() + (w - bounds.getWidth()) / 2;
            //            double dy = y - bounds.getMinY() + (h - bounds.getHeight()) / 2;
            //            transformMolecule(m, dx, dy, 1);
            //        }

            //        for (int i = 0; i < products.length; i++) {
            //            ExtendedMolecule m = products[i];
            //            transformMolecule(m, 0, 0, scale);
            //            Rectangle2D bounds = getBoundingRect(m);
            //            double dx = x + w + arrowSize - bounds.getMinX() + (w - bounds.getWidth()) / 2;
            //            double dy = y - bounds.getMinY() + (h - bounds.getHeight()) / 2;
            //            transformMolecule(m, dx, dy, 1);
            //        }
        }
    }

//    private static void foo(ExtendedMolecule[] mols, double x, double y, double width, double height, double arrowSize)
//    {
//        {
//            double w = (width - arrowSize) / 2;
//            double h = height;
//            Rectangle2D rb = getBoundingRect(mols);
//            if (rb != null) {
//                double scaleHorizontal = width / rb.getWidth();
//                double scaleVertical = width / rb.getHeight();
//                double scale;
//                if (scaleHorizontal < scaleVertical) {   // scale on x-dimension
//                    scale = w / (float) rb.getWidth();
//                } else {                                // scale on y-dimension
//                    scale = h / (float) rb.getHeight();
//                }
//                for (int i = 0; i < mols.length; i++) {
//                    ExtendedMolecule m = mols[i];
//                    transformMolecule(m, 0, 0, scale / 3);
//                    Rectangle2D bounds = getBoundingRect(m);
//                    double dx = x - bounds.getMinX() + (w - bounds.getWidth()) / 2;
//                    double dy = y - bounds.getMinY() + (h - bounds.getHeight()) / 2;
//                    transformMolecule(m, dx, dy, 1);
//                }
//            }
//        }
//    }

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
