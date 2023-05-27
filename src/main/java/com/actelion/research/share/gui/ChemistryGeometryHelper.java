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

package com.actelion.research.share.gui;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.generic.GenericRectangle;

import java.awt.*;


public class ChemistryGeometryHelper
{
    public static final int REACTION_TYPE_NOMOLS = 0;
    public static final int REACTION_TYPE_NOPRODUCTS = 1;
    public static final int REACTION_TYPE_REACTANTS = 2;
    public static final int REACTION_TYPE_NORMAL = 3;

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

    public static GenericRectangle getBoundingRect(Reaction r, boolean includearrows)
    {
        GenericRectangle res = null;
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
            GenericRectangle d = getBoundingRect(mol);
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
                    res = res.union(d);
                    width += d.getWidth();
                } else {
//					System.out.println("Empty bounding rect for molecule " + i);
                }
            }
        }
        if (includearrows && nr == 0 && m > 0) {
            avgwidth = width / m;
            res = new GenericRectangle(res.getX() - avgwidth, res.getY(), res.getWidth() + avgwidth, res.getHeight());
        } else if (includearrows && np == 0 && m > 0) {
            avgwidth = width / m;
            res = new GenericRectangle(res.getX(), res.getY(), res.getWidth() + avgwidth, res.getHeight());
        }
//		System.out.println("Reaction av bnd length: " + width + " getBoundingBox " + res);

        return res;
    }

    public static GenericRectangle getBoundingRect(ExtendedMolecule m)
    {
        double xmax = Double.MIN_VALUE;
        double ymax = Double.MIN_VALUE;
        double xmin = Double.MAX_VALUE;
        double ymin = Double.MAX_VALUE;

        if (m == null)
            return null;

        int na = m.getAllAtoms();
        double bl = 0;
//        if (na > 1)
            bl = m.getAverageBondLength();
        for (int i = 0; i < na; i++) {
            xmax = Math.max(xmax, m.getAtomX(i));
            xmin = Math.min(xmin, m.getAtomX(i));
            ymax = Math.max(ymax, m.getAtomY(i));
            ymin = Math.min(ymin, m.getAtomY(i));
        }

        return (na > 0) ? new GenericRectangle(xmin, ymin, Math.max(xmax - xmin, bl), Math.max(ymax - ymin, bl)) : null;
    }


    public static double getAverageBondLength(Reaction r)
    {
        int rn = r.getMolecules();
        double avg = 0;
        double t = 0;
        for (int i = 0; i < rn; i++) {
            StereoMolecule molecule = r.getMolecule(i);
            if (molecule.getAllAtoms() > 1)
                t += molecule.getAverageBondLength();
        }
        if (rn > 0)
            return t / rn;
        else
            return 0;
    }



    public static  void setAverageBondLength(Reaction rxn, double bndlen)
    {
        double dx = 0;
        int len = rxn.getMolecules();
        for (int fragment=0; fragment<len; fragment++) {
            ExtendedMolecule m = rxn.getMolecule(fragment);
            dx = m.getAverageBondLength();
            double scale = bndlen / dx;
            ChemistryGeometryHelper.transformMolecule(m,0,0,scale);

        }
    }

    public static GenericRectangle getReactantsBoundingRect(Reaction r)
    {
        if (r == null)
            throw new NullPointerException("Cannot pass null reaction");
        int rn = r.getReactants();
        GenericRectangle rc = null;
        for (int i = 0; i < rn; i++) {
            GenericRectangle rt = getBoundingRect(r.getReactant(i));
            if (rc != null && rt != null) {
                GenericRectangle tmp = rc.union(rt);
                rc = new GenericRectangle(tmp.getX(), tmp.getY(), tmp.getWidth(), tmp.getHeight());
            } else
                rc = rt;
        }
        return rc;
    }


    public static GenericRectangle getArrowBoundingRect(Reaction r)
    {
        GenericRectangle rr = getReactantsBoundingRect(r);
        GenericRectangle rp = getProductsBoundingRect(r);
        if (rr != null && rp != null) {
            GenericRectangle union = rr.union(rp);
            double y = union.getHeight() / 2 + union.getY();
            double rx = rr.x+rr.width;
            double px = rp.x;
            return new GenericRectangle(rx < px ? rx : px, y, Math.abs(px - rx), 0);
        } else if (rr != null) {
            double y = rr.getHeight() / 2 + rr.getY();
            double rx = rr.x+rr.width;
            double width = rr.getWidth();
            return new GenericRectangle(rx, y, width, 0);

        } else if (rp != null) {
            double y = rp.getHeight() / 2 + rp.getY();
            double width = 0;
            double rx = rp.x;
            return new GenericRectangle(rx, y, width, 0);
        } else {
            return new GenericRectangle(0, 0, 0, 0);
        }
    }

    public static GenericRectangle getDiffRect(GenericRectangle rr, GenericRectangle rp)
    {
        if (rr != null && rp != null) {
            GenericRectangle union = rr.union(rp);
            double y = union.y;
            double rx = rr.x+rr.width;
            double px = rp.x;
            double ry = union.y;
            double py = union.y+union.height;
            GenericRectangle res = new GenericRectangle(
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

    private static GenericRectangle scaleTo(GenericRectangle src, double width, double height)
    {
        double x = (src.getWidth() - width) / 2 + src.getX();
        double y = (src.getHeight() - height) / 2 + src.getY();
        return new GenericRectangle(x, y, width, height);
    }

    public static GenericRectangle getProductsBoundingRect(Reaction r)
    {
        if (r == null)
            throw new NullPointerException("Cannot pass null reaction");
        int rn = r.getProducts();
        GenericRectangle rc = null;
        for (int i = 0; i < rn; i++) {
            GenericRectangle rt = getBoundingRect(r.getProduct(i));
            if (rc != null && rt != null)
                rc = rc.union(rt);
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

/*
    public static void scaleInto(Reaction reaction, double x, double y, double w, double h)
    {
        GenericRectangle bounds = getBoundingRect(reaction, false);
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
*/
    public static void scaleIntoF(Reaction reaction, double x, double y, double width, double height, double arrowSize)
    {
        GenericRectangle rr = getBoundingRect(reaction,true);
        if (rr != null) {
            double cx = -rr.x;
            double cy = -rr.y;
            ChemistryGeometryHelper.transformReaction(reaction, cx, cy, 1);
            rr = ChemistryGeometryHelper.getBoundingRect(reaction, true);

            double sumWidth = rr.getWidth(), sumHeight = rr.getHeight();

            double scH = width / sumWidth;
            double scV = height / sumHeight;

            double scale = scH;
            if (scH > scV)
                scale = scV;

//            System.out.printf("Scaleinto scale %s\n",scale);
            ChemistryGeometryHelper.transformReaction(reaction, 0, 0, scH);
        }
    }

    public static void scaleInto(Reaction reaction, double x, double y, double width, double height, double arrowSize)
    {

        if (width > arrowSize * 2) {


            ExtendedMolecule[] reactants = getReactants(reaction);
            ExtendedMolecule[] products = getProducts(reaction);

            double sumWidth = 0, sumHeight = 0;
            int numMols = reaction.getMolecules();
            for (int i = 0; i < numMols; i++) {
                ExtendedMolecule m = reaction.getMolecule(i);
                if (m.getAllAtoms() > 1)
                {
                    GenericRectangle r = ChemistryGeometryHelper.getBoundingRect(m);
                    if (r != null) {
                        //                    System.out.printf("MoleculeID %s bounds: %s\n",System.identityHashCode(m),r);
                        sumHeight += r.getHeight();
                        sumWidth += r.getWidth();
                    }
                }
            }

            if (sumHeight == 0 || sumWidth == 0)
                return;

            double scH = width/sumWidth;
            double scV = height/sumHeight;

            double scale = scH;
            if  (scH > scV)
                scale = scV;

            double avbl = getAverageBondLength(reaction);
         //   setAverageBondLength(reaction,avbl);

            double w = (width - arrowSize) / 2;
            double h = height;
//            System.out.printf("Scale = %s vs %s\n",AbstractDepictor.cOptAvBondLen/avbl,scale);
            scale = Math.min(AbstractDepictor.cOptAvBondLen/avbl,scale);
//            if (reactants.length <= 1 && products.length <= 1)
//                scale /= 4;
            {

                GenericRectangle rb = getBoundingRect(reactants);
//                System.out.printf("W/h before scaling prods %s\n",rb);
                transformMolecules(reactants, 0, 0, scale);
                rb = getBoundingRect(reactants);
                double dx = x - rb.x + (w - rb.getWidth()) / 2;
                double dy = y - rb.y + (h - rb.getHeight()) / 2;
//                System.out.printf("W/h after scaling reactants %s,%s $s\n",dx,dy,rb);
                transformMolecules(reactants, dx, dy, 1);
            }

            {
                GenericRectangle pb = getBoundingRect(products);
//                System.out.printf("W/h before scaling prods %s\n",pb);
                transformMolecules(products, 0, 0, scale);
                pb = getBoundingRect(products);
                double dx = x + w + arrowSize - pb.x + (w - pb.getWidth()) / 2;
                double dy = y - pb.y + (h - pb.getHeight()) / 2;
//                System.out.printf("W/h after scaling prods %s,%s %s\n",dx,dy,pb);
                transformMolecules(products, dx, dy, 1);
            }

        }
    }

    public static void scaleIntoOld(Reaction reaction, double x, double y, double width, double height, double arrowSize)
    {

        if (width > arrowSize * 2) {
            ExtendedMolecule[] reactants = getReactants(reaction);
            ExtendedMolecule[] products = getProducts(reaction);

            GenericRectangle rb = getBoundingRect(reactants);
            GenericRectangle pb = getBoundingRect(products);

            System.out.printf("Reactants bounds %s %s %s\n",width,rb.getWidth(),reactants.length);
            System.out.printf("Product bounds %s %s\n",height,pb.getHeight());

            // left and right space
            double w = (width - arrowSize) / 2;
            double h = height;
            double scaleHorizontal = w / Math.max(rb.getWidth(), pb.getWidth());
            double scaleVertical = h / Math.max(rb.getHeight(), pb.getHeight());
                    System.out.printf("Scaling %f vs %f\n", scaleHorizontal, scaleVertical);

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
//            scale *=  10;
//            System.out.printf("Scale = %s\n",scale);
//            if (reactants.length <= 1 && products.length <= 1)
//                scale /= 4;
            {

                transformMolecules(reactants, 0, 0, scale);
                rb = getBoundingRect(reactants);
                double dx = x - rb.x + (w - rb.getWidth()) / 2;
                double dy = y - rb.y + (h - rb.getHeight()) / 2;
                transformMolecules(reactants, dx, dy, 1);
            }

            {
                transformMolecules(products, 0, 0, scale);
                pb = getBoundingRect(products);
                double dx = x + w + arrowSize - pb.x + (w - pb.getWidth()) / 2;
                double dy = y - pb.y + (h - pb.getHeight()) / 2;
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

    public static GenericRectangle getBoundingRect(ExtendedMolecule[] mols)
    {
        if (mols == null || mols.length == 0) {
            return new GenericRectangle(0,0,0,0);
        }

        GenericRectangle r = getBoundingRect(mols[0]);
        for (int i = 1; i < mols.length; i++) {
            GenericRectangle t = getBoundingRect(mols[i]);
            if (t != null) {
                if (r == null)
                    r = t;
                else
                    r = r.union(getBoundingRect(mols[i]));
            }
        }
        if (r == null)
            r = new GenericRectangle(0,0,0,0);
        return r;
    }



    public static void arrangeReaction(Reaction rxn,Dimension size)
    {
        double cx = 0,cy=0;
        if (rxn != null) {
            int len = rxn.getMolecules();
            if (len > 0) {
                double avBndLen = getAverageBondLength(rxn);
                setAverageBondLength(rxn,avBndLen);

                for (int fragment=0; fragment<len; fragment++) {
                    ExtendedMolecule m = rxn.getMolecule(fragment);
                    GenericRectangle rc = ChemistryGeometryHelper.getBoundingRect(m);
                    if (rc != null) {
                        cx += rc.width;
                        cy += rc.height;
                    }
                }
                GenericRectangle rxnBounding = ChemistryGeometryHelper.getArrowBoundingRect(rxn);

                int plus = Math.max(0,len-2);
                double plusWidth = size.width / (len + plus + 1) / 2;
                double dx = cx / len;
                cx = cx + plus * dx / 2 + rxnBounding.getWidth();
                cx += (cx/len);
                double avY = cy / len;
                java.awt.geom.Rectangle2D.Double rb = new java.awt.geom.Rectangle2D.Double(0,0,(double)size.width,(double)size.height);
                double scalex = rb.width / cx;
                double scaley = rb.height / avY;
                double scale = Math.min(scalex,scaley);
                double offsetx = 0;
                for (int fragment=0; fragment<len; fragment++) {
                    ExtendedMolecule m = rxn.getMolecule(fragment);
                    // make the same scale
                    ChemistryGeometryHelper.transformMolecule(m,0,0,scale);

                    GenericRectangle rectBefore = ChemistryGeometryHelper.getBoundingRect(m);
                    // calculate the offset
                    if (rectBefore != null && rb != null) {
                        double offsety = (rb.height - rectBefore.height) / 2;
                        double moveX =  -rectBefore.x + offsetx;
                        double moveY =  -rectBefore.y + offsety;

                        ChemistryGeometryHelper.transformMolecule(m,moveX,moveY,1);
                        GenericRectangle db = ChemistryGeometryHelper.getBoundingRect(m);

                        if (fragment == rxn.getReactants()-1) {
                            offsetx += (db.width + plusWidth * 2);
                        } else {
                            offsetx += (db.width + plusWidth);
                        }
                    }
                }
            }
        }
    }

}
