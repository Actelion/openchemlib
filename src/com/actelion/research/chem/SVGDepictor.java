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
import java.awt.geom.Rectangle2D;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class SVGDepictor extends AbstractDepictor
{
    public static final int DEFAULT_ELEM_WIDTH = 8;
    private static int instanceCnt  = 0;

    private float lineWidth = 1;
    private int textSize = 10;
    private String currentColor = "black";
    private int width = 400;
    private int height = 400;

    private final java.util.List<String> bonds = new ArrayList<String>();
    private final java.util.List<String> atoms = new ArrayList<String>();
    private String id;

    private StringBuffer buffer = new StringBuffer();

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

    private String getId()
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
        String s ="<line " +
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
                "x=\"" + (int) (x - strWidth / 2.0 - 2) + "\" " +
                "y=\"" + (int) (y + textSize / 3 ) + "\" " +
                "font-family=\"Helvetica\" font-size=\"" + textSize + "\" " +
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
        return theString.length() * 5;
    }

    @Override
    protected int getTextSize()
    {
        return textSize;
    }

    @Override
    protected void setTextSize(int theSize)
    {
        textSize = theSize;
//        System.out.printf("Test size set = %d\n", textSize);
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

        return super.simpleValidateView(viewRect, mode);
    }

    public static void main(String... args) throws IOException
    {
        test();
    }

    private static final String SCRIPT =
            "var id = '%s';\n" +
                    "var svgElement = document.getElementById(id);\n" +
                    "\n" +
                    "svgElement.addEventListener('mouseover', function (e) {\n" +
                    "    var id = e.target.id;\n" +
                    "    var split = id.split(':');\n" +
                    "    var element = document.getElementById(id);\n" +
                    "    if (split[1] === 'Bond') {\n" +
                    "        element.setAttribute('stroke-opacity', '0.4');\n" +
                    "        element.setAttribute('stroke', 'red');\n" +
                    "        console.log('Over a bond : ' + split[2]);\n" +
                    "    } else if (split[1] === 'Atom') {\n" +
                    "        element.setAttribute('fill-opacity', '0.4');\n" +
                    "        element.setAttribute('fill', 'navy');\n" +
                    "        console.log('Over an atom : ' + split[2]);\n" +
                    "    }\n" +
                    "});\n" +
                    "\n" +
                    "svgElement.addEventListener('mouseout', function (e) {\n" +
                    "    var id = e.target.id;\n" +
                    "    var split = id.split(':');\n" +
                    "    var element = document.getElementById(id);\n" +
                    "    if (split[1] === 'Bond') {\n" +
                    "        element.setAttribute('stroke-opacity', '0');\n" +
                    "        console.log('Out of bond : ' + split[2]);\n" +
                    "    } else if (split[1] === 'Atom') {\n" +
                    "        element.setAttribute('fill-opacity', '0');\n" +
                    "        console.log('Out of atom : ' + split[2]);\n" +
                    "    }\n" +
                    "});\n";

    private static void test() throws IOException
    {
        String idcode = "fby``@M@QdTRbRvRRbAdijB@jjh`y@@ !BCCW~VpqJCCW~I@pu\u007FelLR`qJCCTLR`puCDhLMPqJ\u007FelLMP@r\u007Fox";
        IDCodeParser p = new IDCodeParser();
        StereoMolecule mol = new StereoMolecule();
        p.parse(mol, idcode);
        CoordinateInventor c = new CoordinateInventor();
        c.invent(mol);
        SVGDepictor d = new SVGDepictor(mol,null);
        d.validateView(null, new Rectangle2D.Float(0, 0, 400, 400), AbstractDepictor.cModeInflateToMaxAVBL);
        d.paint(null);
        String svg = d.toString();
        String script = String.format(SCRIPT, d.getId());
        BufferedWriter os = new BufferedWriter(new FileWriter("svgtest.html"));
        os.write("<!DOCTYPE html " +
                "PUBLIC " +
                "\"-//W3C//DTD XHTML 1.0 Strict//EN\"\n" +
                "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" +
                "\n" +
                "<html>\n" +
                "<head>\n" +
                "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />\n" +
                "<title>Test</title>\n" +
                "</head>\n" +
                "<body>\n\n" +
                "<h1>Molecule SVG Test</h1>\n");
        os.write(svg);
        os.write("\n<script>" + script + "</script>\n");
        os.write("\n</body>\n</html>");
        os.close();


        if (!svg.equals(RESULT)) {
            System.err.println("Test failed!");
        } else
            System.out.println("Test passed!");
        //System.out.println(svg);
    }

    private static final String RESULT = "<svg id=\"mol1\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"400px\" height=\"400px\" viewBox=\"0 0 400 400\">\n" +
            "\t<style> #mol1 {pointer-events:none; }  #mol1 .event  { pointer-events:all;}  </style>\n" +
            "\t<text x=\"162\" y=\"294\" font-family=\"Helvetica\" font-size=\"12\" fill=\"rgb(255,0,0)\">this enantiomer</text>\n" +
            "\t<text x=\"299\" y=\"168\" font-family=\"Helvetica\" font-size=\"14\" fill=\"rgb(244,12,12)\">O</text>\n" +
            "\t<text x=\"278\" y=\"204\" font-family=\"Helvetica\" font-size=\"14\" fill=\"rgb(244,12,12)\">O</text>\n" +
            "\t<text x=\"195\" y=\"132\" font-family=\"Helvetica\" font-size=\"14\" fill=\"rgb(46,76,236)\">N</text>\n" +
            "\t<text x=\"195\" y=\"118\" font-family=\"Helvetica\" font-size=\"14\" fill=\"rgb(46,76,236)\">H</text>\n" +
            "\t<text x=\"130\" y=\"148\" font-family=\"Helvetica\" font-size=\"9\" fill=\"rgb(255,0,0)\">S</text>\n" +
            "\t<text x=\"128\" y=\"136\" font-family=\"Helvetica\" font-size=\"9\" fill=\"rgb(160,0,0)\">abs</text>\n" +
            "\t<line x1=\"283\" y1=\"178\" x2=\"301\" y2=\"167\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"281\" y1=\"174\" x2=\"299\" y2=\"164\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"262\" y1=\"164\" x2=\"262\" y2=\"140\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"258\" y1=\"161\" x2=\"258\" y2=\"142\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"241\" y1=\"176\" x2=\"220\" y2=\"164\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"241\" y1=\"171\" x2=\"224\" y2=\"162\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"241\" y1=\"128\" x2=\"220\" y2=\"140\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"241\" y1=\"132\" x2=\"224\" y2=\"141\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"283\" y1=\"176\" x2=\"283\" y2=\"193\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"283\" y1=\"176\" x2=\"262\" y2=\"164\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"279\" y1=\"202\" x2=\"262\" y2=\"212\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"262\" y1=\"164\" x2=\"241\" y2=\"176\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"262\" y1=\"212\" x2=\"262\" y2=\"236\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"262\" y1=\"140\" x2=\"241\" y2=\"128\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"220\" y1=\"140\" x2=\"203\" y2=\"130\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"196\" y1=\"130\" x2=\"179\" y2=\"140\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"179\" y1=\"140\" x2=\"158\" y2=\"128\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"158\" y1=\"128\" x2=\"137\" y2=\"140\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"137\" y1=\"140\" x2=\"116\" y2=\"128\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<polygon points=\"137,140 136,152 138,152 \" style=\"fill:rgb(160,0,0);stroke:rgb(160,0,0);stroke-width:1\"/>\n" +
            "\t<polygon points=\"138,152 136,152 134,164 140,164 \" style=\"fill:rgb(160,0,0);stroke:rgb(160,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"116\" y1=\"128\" x2=\"96\" y2=\"140\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"213\" y1=\"260\" x2=\"192\" y2=\"272\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line x1=\"220\" y1=\"164\" x2=\"220\" y2=\"140\" style=\"stroke:rgb(0,0,0);stroke-width:1\"/>\n" +
            "\t<line id=\"mol1:Bond:0-1\" class=\"event\" x1=\"303\" y1=\"164\" x2=\"283\" y2=\"176\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:3-5\" class=\"event\" x1=\"262\" y1=\"164\" x2=\"262\" y2=\"140\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:6-9\" class=\"event\" x1=\"241\" y1=\"176\" x2=\"220\" y2=\"164\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:8-10\" class=\"event\" x1=\"241\" y1=\"128\" x2=\"220\" y2=\"140\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:1-2\" class=\"event\" x1=\"283\" y1=\"176\" x2=\"283\" y2=\"200\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:1-3\" class=\"event\" x1=\"283\" y1=\"176\" x2=\"262\" y2=\"164\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:2-4\" class=\"event\" x1=\"283\" y1=\"200\" x2=\"262\" y2=\"212\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:3-6\" class=\"event\" x1=\"262\" y1=\"164\" x2=\"241\" y2=\"176\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:4-7\" class=\"event\" x1=\"262\" y1=\"212\" x2=\"262\" y2=\"236\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:5-8\" class=\"event\" x1=\"262\" y1=\"140\" x2=\"241\" y2=\"128\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:10-11\" class=\"event\" x1=\"220\" y1=\"140\" x2=\"200\" y2=\"128\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:11-12\" class=\"event\" x1=\"200\" y1=\"128\" x2=\"179\" y2=\"140\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:12-13\" class=\"event\" x1=\"179\" y1=\"140\" x2=\"158\" y2=\"128\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:13-14\" class=\"event\" x1=\"158\" y1=\"128\" x2=\"137\" y2=\"140\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:14-15\" class=\"event\" x1=\"137\" y1=\"140\" x2=\"116\" y2=\"128\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:14-16\" class=\"event\" x1=\"137\" y1=\"140\" x2=\"137\" y2=\"164\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:15-17\" class=\"event\" x1=\"116\" y1=\"128\" x2=\"96\" y2=\"140\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:18-19\" class=\"event\" x1=\"213\" y1=\"260\" x2=\"192\" y2=\"272\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<line id=\"mol1:Bond:9-10\" class=\"event\" x1=\"220\" y1=\"164\" x2=\"220\" y2=\"140\" stroke-width=\"8\" stroke-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:0\" class=\"event\" cx=\"303\" cy=\"164\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:1\" class=\"event\" cx=\"283\" cy=\"176\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:2\" class=\"event\" cx=\"283\" cy=\"200\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:3\" class=\"event\" cx=\"262\" cy=\"164\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:4\" class=\"event\" cx=\"262\" cy=\"212\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:5\" class=\"event\" cx=\"262\" cy=\"140\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:6\" class=\"event\" cx=\"241\" cy=\"176\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:7\" class=\"event\" cx=\"262\" cy=\"236\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:8\" class=\"event\" cx=\"241\" cy=\"128\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:9\" class=\"event\" cx=\"220\" cy=\"164\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:10\" class=\"event\" cx=\"220\" cy=\"140\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:11\" class=\"event\" cx=\"200\" cy=\"128\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:12\" class=\"event\" cx=\"179\" cy=\"140\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:13\" class=\"event\" cx=\"158\" cy=\"128\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:14\" class=\"event\" cx=\"137\" cy=\"140\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:15\" class=\"event\" cx=\"116\" cy=\"128\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:16\" class=\"event\" cx=\"137\" cy=\"164\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:17\" class=\"event\" cx=\"96\" cy=\"140\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:18\" class=\"event\" cx=\"213\" cy=\"260\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "\t<circle id=\"mol1:Atom:19\" class=\"event\" cx=\"192\" cy=\"272\" r=\"8\" fill-opacity=\"0\"/>\n" +
            "</svg>";


}
