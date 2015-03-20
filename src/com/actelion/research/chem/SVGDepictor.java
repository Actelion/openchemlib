/*
* Copyright (c) 1997 - 2015
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

import java.awt.*;
import java.awt.font.GlyphVector;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class SVGDepictor extends AbstractDepictor
{
    public static final int DEFAULT_ELEM_WIDTH = 8;

    private static final String FONTNAME = "Helvetica";
    private static int instanceCnt = 0;

    private float lineWidth = 1;
    private int textSize = 10;
    private int width = 400;
    private int height = 400;

    private String currentColor = "black";
    private final java.util.List<String> bonds = new ArrayList<String>();
    private final java.util.List<String> atoms = new ArrayList<String>();

    private String id;
    private StringBuffer buffer = new StringBuffer();

    private Font currentFont = new Font(FONTNAME, Font.PLAIN, 12);
    private Graphics2D graphics;

    public SVGDepictor(StereoMolecule mol, String id)
    {
        this(mol, 0, id);
    }

    public SVGDepictor(StereoMolecule mol, int displayMode, String id)
    {
        super(mol, displayMode);
        this.id = id;
        instanceCnt++;
    }

    public static final String makeColor(int r, int g, int b)
    {
        return "rgb(" + r + "," + g + "," + b + ")";
    }

    public String getId()
    {
        return id != null ? id : ("mol" + instanceCnt);
    }

    private void write(String s)
    {
        buffer.append("\t");
        buffer.append(s);
        buffer.append("\n");
    }

    @Override
    protected void drawBlackLine(DepictorLine theLine)
    {
        int x1 = (int) theLine.x1;
        int x2 = (int) theLine.x2;
        int y1 = (int) theLine.y1;
        int y2 = (int) theLine.y2;
        String s = "<line " +
                "x1=\"" + x1 + "\" " +
                "y1=\"" + y1 + "\" " +
                "x2=\"" + x2 + "\" " +
                "y2=\"" + y2 + "\" " +
                "style=\"stroke:" + currentColor + ";" +
                "stroke-width:" + (int) (lineWidth) + "\"/>";
        write(s);
    }

    @Override
    protected void drawDottedLine(DepictorLine theLine)
    {
        int x1 = (int) theLine.x1;
        int x2 = (int) theLine.x2;
        int y1 = (int) theLine.y1;
        int y2 = (int) theLine.y2;
        String s = "<line stroke-dasharray=\"3, 3\" " +
                "x1=\"" + x1 + "\" " +
                "y1=\"" + y1 + "\" " +
                "x2=\"" + x2 + "\" " +
                "y2=\"" + y2 + "\" " +
                "stroke=\"" + currentColor + "\" " +
                "stroke-width:" + (int) (lineWidth) + "\"/>";

        write(s);
    }

    @Override
    protected void drawPolygon(float[] x, float[] y, int count)
    {
        StringBuilder s = new StringBuilder("<polygon points=\"");
        for (int i = 0; i < count; i++) {
            s.append((int) x[i]);
            s.append(",");
            s.append((int) y[i]);
            s.append(" ");
        }
        s.append("\" " +
                "style=\"fill:" + currentColor + ";" +
                "stroke:" + currentColor + ";" +
                "stroke-width:1\"/>");
        write(s.toString());
    }

    @Override
    protected void drawString(String theString, float x, float y)
    {

        float strWidth = getStringWidth(theString);
        String s = "<text " +
                "x=\"" + (int) (x - strWidth / 2.0) + "\" " +
                "y=\"" + (int) (y + textSize / 3) + "\" " +
                "font-family=\" " + currentFont.getName() + "\" " +
                "font-size=\"" + currentFont.getSize() + "\" " +
                "fill=\"" + currentColor + "\">" + theString +
                "</text>";
        write(s);
    }

    @Override
    protected void fillCircle(float x, float y, float r)
    {
        String s = "<circle " +
                "cx=\"" + (int) x + "\" " +
                "cy=\"" + (int) y + "\" " +
                "r=\"" + (int) r + "\" " +
                "fill=\"" + currentColor + "\" />";
        write(s);
    }

    @Override
    protected float getLineWidth()
    {
        return lineWidth;
    }

    @Override
    protected float getStringWidth(String theString)
    {
        GlyphVector mCurrentGlyphVector = currentFont.createGlyphVector(graphics.getFontRenderContext(), theString);
        return (float) mCurrentGlyphVector.getLogicalBounds().getWidth();
    }

    @Override
    protected int getTextSize()
    {
        return textSize;
    }

    @Override
    protected void setTextSize(int theSize)
    {
        if (textSize != theSize) {
            textSize = theSize;
            currentFont = new Font(FONTNAME, Font.PLAIN, theSize);
        }
    }

    @Override
    protected void setLineWidth(float width)
    {
        lineWidth = width;
    }

    @Override
    protected void setColor(Color theColor)
    {
        currentColor = makeColor(theColor.getRed(), theColor.getGreen(), theColor.getBlue());
    }

    @Override
    protected void onDrawBond(int atom1, int atom2, float x1, float y1, float x2, float y2)
    {
        String s = "<line " +
                "id=\"" + getId() + ":Bond:" + atom1 + "-" + atom2 + "\" " +
                "class=\"event\" " +
                "x1=\"" + (int) (x1) + "\" " +
                "y1=\"" + (int) (y1) + "\" " +
                "x2=\"" + (int) (x2) + "\" " +
                "y2=\"" + (int) (y2) + "\" " +
                "stroke-width=\"" + DEFAULT_ELEM_WIDTH + "\" " +
                "stroke-opacity=\"0\"" +
                "/>";
        bonds.add(s);
    }

    @Override
    protected void onDrawAtom(int atom, String symbol, float x, float y)
    {
        int r = DEFAULT_ELEM_WIDTH;
        String s = "<circle " +
                "id=\"" + getId() + ":Atom:" + atom + "\" " +
                "class=\"event\" " +
                "cx=\"" + (int) (x) + "\" " +
                "cy=\"" + (int) (y) + "\" " +
                "r=\"" + r + "\" " +
                "fill-opacity=\"0\"/>";
        atoms.add(s);
    }


    @Override
    public String toString()
    {
        String header =
                "<svg " +
                        "id=\"" + getId() + "\" " +
                        "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" " +
                        "width=\"" + width + "px\" " +
                        "height=\"" + height + "px\" " +
                        "viewBox=\"0 0 " + width + " " + height + "\">\n";

        String footer = "</svg>";

        String style = "<style>" +
                " #" + getId() +
                " {pointer-events:none; } " +
                " #" + getId() + " .event " +
                " { pointer-events:all;} " +
                " </style>\n";
        String rect =
                "<rect " +
                        "x=\"" + 0 + "\" " +
                        "y=\"" + 0 + "\" " +
                        "width='" + width + "' " +
                        "height='" + height + "' style='fill-opacity:0;stroke:red;stroke-width:3'/>\n";
        header += "\t";
        header += style;
//        header += rect;


        // Append the (invisible) bond lines
        for (String b : bonds) {
            write(b);
        }
        // Append the (invisible) atom circles
        for (String a : atoms) {
            write(a);
        }

        return header + buffer.toString() + footer;
    }

    @Override
    protected DepictorTransformation simpleValidateView(Rectangle2D.Float viewRect, int mode)
    {

        width = (int) viewRect.getWidth();
        height = (int) viewRect.getHeight();
        BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        graphics = img.createGraphics();

        return super.simpleValidateView(viewRect, mode);
    }


}
