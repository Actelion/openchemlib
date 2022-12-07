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
* 3. Neither the name of the copyright holder nor the
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

import com.actelion.research.gui.generic.GenericPolygon;
import com.actelion.research.gui.generic.GenericRectangle;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

public class SVGDepictor extends AbstractDepictor {
    public static final int DEFAULT_ELEM_WIDTH = 8;

    private static final String FONTNAME = "Helvetica";
    private static int instanceCnt = 0;

    private double lineWidth = 1;
    private int textSize = 10;
    private int width = 400;
    private int height = 400;
    private boolean legacyMode = true;

    private String currentColor = "black";
    private final java.util.List<String> bonds = new ArrayList<String>();
    private final java.util.List<String> atoms = new ArrayList<String>();

    private String id;
    private StringBuilder buffer = new StringBuilder();

    private Font currentFont = new Font(FONTNAME, Font.PLAIN, 12);
    private Graphics2D graphics;

    public SVGDepictor(StereoMolecule mol, String id)
    {
        this(mol, 0, id);
    }

    public SVGDepictor(StereoMolecule mol, int displayMode, String id) {
        super(mol, displayMode);
        this.id = id;
        instanceCnt++;
    }

    /**
     * For legacy reasons the default is to include invisible event atoms and bonds.
     * Call this method after contruction to skip these non-visible contributions to the SVG.
     * @param b
     */
    public void setLegacyMode(boolean b) {
        legacyMode = b;
    }

    public static final String makeColor(int r, int g, int b) {
        return "rgb(" + r + "," + g + "," + b + ")";
    }

    public String getId() {
        return id != null ? id : ("mol" + instanceCnt);
    }

    private void write(String s) {
        buffer.append("\t");
        buffer.append(s);
        buffer.append("\n");
    }

    @Override
    protected void drawBlackLine(DepictorLine theLine) {
        int x1 = (int) theLine.x1;
        int x2 = (int) theLine.x2;
        int y1 = (int) theLine.y1;
        int y2 = (int) theLine.y2;
        String s = "<line " +
                "x1=\"" + x1 + "\" " +
                "y1=\"" + y1 + "\" " +
                "x2=\"" + x2 + "\" " +
                "y2=\"" + y2 + "\" " +
                "style=\"stroke:" + currentColor + "; stroke-width:" + lineWidth + "\"/>";
        write(s);
    }

    @Override
    protected void drawDottedLine(DepictorLine theLine) {
        int x1 = (int) theLine.x1;
        int x2 = (int) theLine.x2;
        int y1 = (int) theLine.y1;
        int y2 = (int) theLine.y2;
        int d = (int)(3*lineWidth);
        String s = "<line stroke-dasharray=\""+d+","+d+"\" " +
                "x1=\"" + x1 + "\" " +
                "y1=\"" + y1 + "\" " +
                "x2=\"" + x2 + "\" " +
                "y2=\"" + y2 + "\" " +
                "style=\"stroke:" + currentColor + "; stroke-width:" + lineWidth + "\"/>";

        write(s);
    }

    @Override
    protected void drawPolygon(GenericPolygon p) {
        StringBuilder s = new StringBuilder("<polygon points=\"");
        for (int i=0; i<p.getSize(); i++) {
            s.append(Math.round(p.getX(i)));
            s.append(",");
            s.append(Math.round(p.getY(i)));
            s.append(" ");
        }
        s.append("\" " +
                "style=\"fill:" + currentColor + "; stroke:" + currentColor + "; stroke-width:"+lineWidth+"\"/>");
        write(s.toString());
    }

    @Override
    protected void drawString(String theString, double x, double y) {
        double strWidth = getStringWidth(theString);
        String s = "<text " +
                "x=\"" + (int) (x - strWidth / 2.0) + "\" " +
                "y=\"" + (int) (y + textSize / 3) + "\" " +
                "stroke=\"none\" " +
//                "font-family=\" " + currentFont.getName() + "\" " +
                "font-size=\"" + currentFont.getSize() + "\" " +
                "fill=\"" + currentColor + "\">" + theString +
                "</text>";
        write(s);
    }

    @Override
    protected void fillCircle(double x, double y, double d) {
        String s = "<circle " +
                "cx=\"" + (int) (x+d/2) + "\" " +
                "cy=\"" + (int) (y+d/2) + "\" " +
                "r=\"" + (int) (d/2) + "\" " +
                "fill=\"" + currentColor + "\" />";
        write(s);
    }

    @Override
    protected double getLineWidth()
    {
        return lineWidth;
    }

    @Override
    protected double getStringWidth(String theString) {
        float ret =  (float)currentFont.getStringBounds(theString,graphics.getFontRenderContext()).getWidth();
        return ret;
    }

    @Override
    protected int getTextSize() {
        return textSize;
    }

    @Override
    protected void setTextSize(int theSize) {
        if (textSize != theSize) {
            textSize = theSize;
            currentFont = new Font(FONTNAME, Font.PLAIN, theSize);
        }
    }

    @Override
    protected void setLineWidth(double width) {
        lineWidth = Math.round(100 * Math.max(width, 1.0)) / 100;
    }

    @Override
    protected void setRGB(int rgb) {
        currentColor = makeColor((rgb & 0x00FF0000) >> 16, (rgb & 0x0000FF00) >> 8, rgb & 0x000000FF);
    }

    @Override
    protected void onDrawBond(int bond, double x1, double y1, double x2, double y2) {
        String s = "<line " +
                "id=\"" + getId() + ":Bond:" + bond + "\" " +
                "class=\"event\" " +	// class to respond to the mouse event
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
    protected void onDrawAtom(int atom, String symbol, double x, double y) {
        int r = DEFAULT_ELEM_WIDTH;
        String s = "<circle " +
                "id=\"" + getId() + ":Atom:" + atom + "\" " +
                "class=\"event\" " + // class to respond to the mouse event
                "cx=\"" + (int) (x) + "\" " +
                "cy=\"" + (int) (y) + "\" " +
                "r=\"" + r + "\" " +
                "fill-opacity=\"0\"/>";
        atoms.add(s);
    }


    @Override
    public String toString() {
        String header = "<svg " +
                        "id=\"" + getId() + "\" " +
                        "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" " +
                        "width=\"" + width + "px\" " +
                        "height=\"" + height + "px\" " +
                        "viewBox=\"0 0 " + width + " " + height + "\">\n";

        String footer = "</svg>";

        String style = legacyMode ?
                "<style>" +
                " #" + getId() +
                " {pointer-events:none; } " +	// Disable Mouse events on the root element so they get passed to the childs
                " #" + getId() + " .event " +
                " { pointer-events:all;} " +	// Enable Mouse events for elements possessing the class "event"
                " </style>\n"
              : "<g style=\"font-size:"+getTextSize()+"px; fill-opacity:1; stroke-opacity:1; fill:black; stroke:black;"
              + " font-weight:normal; text-rendering:optimizeLegibility; font-family:sans-serif;"
              + " stroke-linejoin:round; stroke-linecap:round; stroke-dashoffset:0;\">";

//        String rect = "<rect " +
//                      "x=\"" + 0 + "\" " +
//                      "y=\"" + 0 + "\" " +
//                      "width='" + width + "' " +
//                      "height='" + height + "' style='fill-opacity:0;stroke:red;stroke-width:3'/>\n";
        header += "\t";
        header += style;
//        header += rect;

        if (legacyMode) {
            // Append the (invisible) bond lines
            for (String b : bonds)
                write(b);

            // Append the (invisible) atom circles
            for (String a : atoms)
                write(a);
        }

        if (!legacyMode)
            write("</g>");

        return header + buffer.toString() + footer;
    }

    @Override
    public DepictorTransformation simpleValidateView(GenericRectangle viewRect, int mode) {
        width = (int) viewRect.getWidth();
        height = (int) viewRect.getHeight();
        BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        graphics = img.createGraphics();

        return super.simpleValidateView(viewRect, mode);
    }
}
