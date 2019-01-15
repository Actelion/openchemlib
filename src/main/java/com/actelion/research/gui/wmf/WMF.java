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
package com.actelion.research.gui.wmf;

import java.awt.*;
import java.io.*;
import java.util.Locale;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;

@SuppressWarnings("deprecation")
public class WMF extends MetaFile implements WMFConstants
{

    private int maxObjectSize;
    private ByteArrayOutputStream wmf;
    private Vector<Boolean> handles;
    private static boolean loaded = false;
    private static String[][] fontnames = {
        {
            "helvetica",
            "Arial"
        },
        {
            "timesroman",
            "Times New Roman"
        },
        {
            "courier",
            "Courier New"
        },
        {
            "zapfdingsbat",
            "Windings"
        },
        {
            "dialog",
            "Arial"
        },
        {
            "dialoginput",
            "Courier New"
        },
        {
            "serif",
            "Times New Roman"
        },
        {
            "sansserif",
            "Arial"
        },
        {
            "monospaced",
            "Courier New"
        }
    };

//    private byte[] word;

    public static final int MFCOMMENT = 15;

    public WMF()
    {
        maxObjectSize = 0;
        wmf = new ByteArrayOutputStream();
        handles = new Vector<Boolean>();
        setupFontNames();
    }

    protected int addHandle()
    {
        for (int i = 0; i < handles.size(); i++)
            if (!handles.elementAt(i)) {
                handles.setElementAt(Boolean.TRUE, i);

                return i;
            }

        handles.addElement(Boolean.TRUE);

        return handles.size() - 1;
    }

    @Override
    public void arc(int left, int top, int right, int bottom, int xstart, int ystart, int xend, int yend)
    {
        writeRecordHeader(META_ARC, 8);
        writeWord(yend);
        writeWord(xend);
        writeWord(ystart);
        writeWord(xstart);
        writeWord(bottom);
        writeWord(right);
        writeWord(top);
        writeWord(left);
    }


    @Override
    public int createBrushIndirect(int style, Color color, int hatch)
    {
        writeRecordHeader(META_CREATEBRUSHINDIRECT, 4);
        writeWord(style);
        writeColor(color);
        writeWord(hatch);

        return addHandle();
    }

    @Override
    public int createFont(Font font, int i, boolean flag, boolean flag1)
    {
        char c = '\u0190';
        if (font.isBold()) {
            c = '\u02BC';
        }

        return createFont(
            -font.getSize(), 0, i, 0, c, font.isItalic(), flag, flag1, (byte) 0, (byte) 0, (byte) 0,
            (byte) 0, (byte) 0, translateFontName(font.getName()));
    }

    @Override
    public int createFont(
        int height, int with, int esc, int orient, int weight, boolean italic, boolean underline, boolean strikeOut, byte charSet,
        byte outPrecision, byte clipPrecision, byte quality, byte pitchAndFamily, String s)
    {
        writeRecordHeader(META_CREATEFONTINDIRECT, 9 + ((s.length() + 2) / 2));
        writeWord(height);
        writeWord(with);
        writeWord(esc);
        writeWord(orient);
        writeWord(weight);

        int value = 0;

        if (italic) {
            value = 1;
        }

        if (underline) {
            value += 256;
        }

        writeWord(value);
        value = (charSet << 8) & 0xff00;

        if (strikeOut) {
            value++;
        }

        writeWord(value);
        writeWord(outPrecision | ((clipPrecision << 8) & 0xff00));
        writeWord(quality | ((pitchAndFamily << 8) & 0xff00));

        byte[] buffer = new byte[s.length() + 2];
        s.getBytes(0, s.length(), buffer, 0);
        buffer[buffer.length - 2] = 0;
        buffer[buffer.length - 1] = 0;

        for (int i = 0; i < (buffer.length / 2); i++) {
            writeWord(buffer[i * 2] | ((buffer[(i * 2) + 1] << 8) & 0xff00));
        }
        return addHandle();
    }

    @Override
    public int createPatternBrush(int[] ai, int i, int j)
    {
        int k = (((i * 3) + 3) / 4) * 4;
        writeRecordHeader(META_DIBCREATEPATTERNBRUSH, 22 + ((k / 2) * j));
        writeWord(3);
        writeWord(0);
        writeBitmap(ai, i, j);

        return addHandle();
    }

    @Override
    public int createPenIndirect(int style, int width, Color color)
    {
        writeRecordHeader(META_CREATEPENINDIRECT, 5);
        writeWord(style);
        writeInteger(width);
        writeColor(color);

        return addHandle();
    }

    @Override
    public void deleteObject(int i)
    {
        if ((i < handles.size()) && handles.elementAt(i)) {
            writeRecordHeader(META_DELETEOBJECT, 1);
            writeWord(i);
            handles.setElementAt(Boolean.FALSE, i);
        } else {
            throw new ArrayIndexOutOfBoundsException();
        }
    }

//    @Override
//    private void deleteObjects()
//    {
//        for (int i = 0; i < handles.size(); i++)
//            if (handles.elementAt(i)) {
//                deleteObject(i);
//            }
//    }

    @Override
    public void ellipse(int left, int top, int right, int bottom)
    {
        writeRecordHeader(META_ELLIPSE, 4);
        writeWord(bottom);
        writeWord(right);
        writeWord(top);
        writeWord(left);
    }

    @Override
    public void escape(int function, byte[] data)
    {
        writeRecordHeader(META_ESCAPE, 2 + ((data.length + 1) / 2));
        writeWord(function);
        writeWord(data.length);

        byte[] copy = new byte[((data.length + 1) / 2) * 2];
        System.arraycopy(data, 0, copy, 0, data.length);
        for (int j = 0; j < copy.length; j += 2) {
            // Fixed the following line by and-ing the first byte with 0xFF
            int t = (copy[j] & 0x00ff) | ((copy[j + 1] << 8) & 0xff00);
            writeWord(t);
        }

    }

    private int getBodySize()
    {
        return wmf.size() / 2;
    }

    private String getFontProperty(Properties properties, String s) throws Exception
    {
        String s1 = properties.getProperty(s + ".0");

        if (s1 == null) {
            s1 = properties.getProperty(s + ".plain.0");
        }

        if (s1 == null) {
            throw new Exception(s + " not found");
        } else {
            StringTokenizer stringtokenizer = new StringTokenizer(s1, ",");

            return stringtokenizer.nextToken().trim();
        }
    }

    private int highWord(int i)
    {
        return (i & 0xffff0000) >> 16;
    }

    public void intersectClipRect(int i, int j, int k, int l)
    {
        writeRecordHeader(META_INTERSECTCLIPRECT, 4);
        writeWord(l);
        writeWord(k);
        writeWord(j);
        writeWord(i);
    }

    public void lineTo(int x, int y)
    {
        writeRecordHeader(META_LINETO, 2);
        writeWord(y);
        writeWord(x);
    }

    private void setupFontNames()
    {
        if (!loaded) {
            try {
                String language = Locale.getDefault().getLanguage();
                String country = Locale.getDefault().getCountry();
                String propsLib = "lib" + File.separatorChar + "font.properties";
                File file = new File(System.getProperty("java.home"), propsLib + '.' + language + '_' + country);

                if (!file.exists()) {
                    file = new File(System.getProperty("java.home"), propsLib + '.' + language);
                }

                if (!file.exists()) {
                    file = new File(System.getProperty("java.home"), propsLib);
                }

                if (file.exists()) {
                    Properties properties = new Properties();
                    properties.load(new FileInputStream(file));

                    String[][] as = new String[5][2];
                    int i = 0;
                    as[i][0] = "dialog";
                    as[i][1] = getFontProperty(properties, as[i][0]);
                    i++;
                    as[i][0] = "dialoginput";
                    as[i][1] = getFontProperty(properties, as[i][0]);
                    i++;
                    as[i][0] = "serif";
                    as[i][1] = getFontProperty(properties, as[i][0]);
                    i++;
                    as[i][0] = "sansserif";
                    as[i][1] = getFontProperty(properties, as[i][0]);
                    i++;
                    as[i][0] = "monospaced";
                    as[i][1] = getFontProperty(properties, as[i][0]);
                    fontnames = as;
                    loaded = true;
                }
            } catch (Exception _ex) {
            }
        }
    }

    private int lowWord(int i)
    {
        return i & 0xffff;
    }

    private void maxObjectSize(int i)
    {
        if (i > maxObjectSize) {
            maxObjectSize = i;
        }
    }

    protected void writeRecordHeader(int record, int size)
    {
        int totalSize = size + 3;
        writeInteger(totalSize);
        writeWord(record);
        maxObjectSize(totalSize);
    }

    public void moveTo(int x, int y)
    {
        writeRecordHeader(META_MOVETO, 2);
        writeWord(y);
        writeWord(x);
    }

    private void outputInteger(OutputStream outputstream, int i) throws IOException
    {
        outputWord(outputstream, lowWord(i));
        outputWord(outputstream, highWord(i));
    }

    private void outputWord(OutputStream outputstream, int i) throws IOException
    {
        byte[] word = new byte[2];
        word[0] = (byte) (i & 0x00ff);
        word[1] = (byte) ((i & 0xff00) >> 8);
        outputstream.write(word);
    }


    public void pie(int left, int top, int right, int bottom, int xR1, int yR1, int xR2, int yR2)
    {
        writeRecordHeader(META_PIE, 8);
        writeWord(yR2);
        writeWord(xR2);
        writeWord(yR1);
        writeWord(xR1);
        writeWord(bottom);
        writeWord(right);
        writeWord(top);
        writeWord(left);
    }

    public void polygon(int[] ptx, int[] pty, int count)
    {
        writeRecordHeader(META_POLYGON, 1 + (2 * count));
        writeWord(count);

        for (int j = 0; j < count; j++) {
            writeWord(ptx[j]);
            writeWord(pty[j]);
        }
    }

    @Override
    public void polyline(int[] ptx, int[] pty, int count)
    {
        writeRecordHeader(META_POLYLINE, 1 + (2 * count));
        writeWord(count);

        for (int j = 0; j < count; j++) {
            writeWord(ptx[j]);
            writeWord(pty[j]);
        }
    }

    @Override
    public void polypolygon(Polygon[] apolygon)
    {
        int i = 0;

        for (int j = 0; j < apolygon.length; j++)
            i += apolygon[j].npoints;

        writeRecordHeader(META_POLYPOLYGON, 1 + apolygon.length + (2 * i));
        writeWord(apolygon.length);

        for (int k = 0; k < apolygon.length; k++)
            writeWord(apolygon[k].npoints);

        for (int l = 0; l < apolygon.length; l++) {
            for (int i1 = 0; i1 < apolygon[l].npoints; i1++) {
                writeWord(apolygon[l].xpoints[i1]);
                writeWord(apolygon[l].ypoints[i1]);
            }
        }
    }

    @Override
    public void rectangle(int left, int top, int right, int bottom)
    {
        writeRecordHeader(META_RECTANGLE, 4);
        writeWord(bottom);
        writeWord(right);
        writeWord(top);
        writeWord(left);
    }

    @Override
    public void roundRect(int left, int top, int right, int bottom, int width, int height)
    {
        writeRecordHeader(META_ROUNDRECT, 6);
        writeWord(height);
        writeWord(width);
        writeWord(bottom);
        writeWord(right);
        writeWord(top);
        writeWord(left);
    }


    @Override
    public void selectObject(int handle)
    {
        if ((handle < handles.size()) && ((Boolean) handles.elementAt(handle)).booleanValue()) {
            writeRecordHeader(META_SELECTOBJECT, 1);
            writeWord(handle);
        } else {
            throw new ArrayIndexOutOfBoundsException();
        }
    }

    @Override
    public void setBKColor(Color color)
    {
        writeRecordHeader(META_SETBKCOLOR, 2);
        writeColor(color);
    }

    @Override
    public void setBKMode(int mode)
    {
        writeRecordHeader(META_SETBKMODE, 1);
        writeWord(mode);
    }

    @Override
    public void setClipRgn()
    {
        writeRecordHeader(META_SELECTCLIPREGION, 1);
        writeWord(0);
    }

    @Override
    public void setMapMode(int mode)
    {
        writeRecordHeader(META_SETMAPMODE, 1);
        writeWord(mode);
    }

    @Override
    public void setPixel(int x, int y, Color color)
    {
        writeRecordHeader(META_SETPIXEL, 4);
        writeColor(color);
        writeWord(y);
        writeWord(x);
    }

    @Override
    public void setPolyFillMode(int mode)
    {
        writeRecordHeader(META_SETPOLYFILLMODE, 1);
        writeWord(mode);
    }

    @Override
    public void setROP2(int mode)
    {
        writeRecordHeader(META_SETROP2, 1);
        writeWord(mode);
    }

    @Override
    public void setStretchBltMode(int mode)
    {
        writeRecordHeader(META_SETSTRETCHBLTMODE, 1);
        writeWord(mode);
    }

    @Override
    public void setTextAlign(int i)
    {
        writeRecordHeader(META_SETTEXTALIGN, 1);
        writeWord(i);
    }

    @Override
    public void setTextCharacterExtra(int i)
    {
        writeRecordHeader(META_SETTEXTCHAREXTRA, 1);
        writeWord(i);
    }

    @Override
    public void setTextColor(Color color)
    {
        writeRecordHeader(META_SETTEXTCOLOR, 2);
        writeColor(color);
    }

//    public void setTextJustification(int i, int j)
//    {
//        writeRecordHeader(META_SETTEXTJUSTIFICATION, 2);
//        writeWord(j);
//        writeWord(i);
//    }
//
//    public void setTranslateFontNames(String[][] as)
//    {
//        fontnames = as;
//    }

    @Override
    public void setViewportExt(int i, int j)
    {
//        System.out.println("setViewportExt " + i + " " + j);
        writeRecordHeader(META_SETVIEWPORTEXT, 2);
        writeWord(j);
        writeWord(i);
    }

    @Override
    public void setWindowExt(int cx, int cy)
    {
//        System.out.println("setWindowExt " + cx + " " + cy);
        writeRecordHeader(META_SETWINDOWEXT, 2);
        writeWord(cy);
        writeWord(cx);
    }

    @Override
    public void setWindowOrg(int x, int y)
    {
        writeRecordHeader(META_SETWINDOWORG, 2);
        writeWord(y);
        writeWord(x);
    }

    @Override
    public void stretchBlt(
        int xOrigDest, int yOrigDest, int widthDest, int heightDest, int xOrigSrc, int yOrigSrc, int widthSrc, int heightSrc, int rasterOp, int[] pixelData, int imageWidth, int imageHeight)
    {
        int l2 = (((imageWidth * 3) + 3) / 4) * 4;
        writeRecordHeader(META_DIBSTRETCHBLT, 30 + ((l2 / 2) * imageHeight));
        writeInteger(rasterOp);
        writeWord(heightSrc);
        writeWord(widthSrc);
        writeWord(yOrigSrc);
        writeWord(xOrigSrc);
        writeWord(heightDest);
        writeWord(widthDest);
        writeWord(yOrigDest);
        writeWord(xOrigDest);
        writeBitmap(pixelData, imageWidth, imageHeight);
    }

    @Override
    public void textOut(int x, int y, String s)
    {
        writeRecordHeader(META_TEXTOUT, 3 + ((s.length() + 1) / 2));
        writeWord(s.length());

        byte[] buffer = new byte[s.length() + 1];
        s.getBytes(0, s.length(), buffer, 0);
        buffer[buffer.length - 1] = 0;

        for (int k = 0; k < (buffer.length / 2); k++) {
            int l = 0;
            l = buffer[k * 2] | ((buffer[(k * 2) + 1] << 8) & 0xff00);
            writeWord(l);
        }

        writeWord(y);
        writeWord(x);
    }

    @Override
    public String translateFontName(String s)
    {
        String s1 = s;
        String s2 = s.toLowerCase();

        for (int i = 0; i < fontnames.length; i++) {
            if (!s2.equals(fontnames[i][0])) {
                continue;
            }

            s1 = fontnames[i][1];

            break;
        }

        return s1;
    }

    protected void writeBitmap(int[] ai, int imageWidth, int imageHeight)
    {
        int alignedWidth = (((imageWidth * 3) + 3) / 4) * 4;
        byte[] buffer = new byte[alignedWidth * imageHeight];

        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++) {
                int index = (x * 3) + (y * alignedWidth);
                int temp = ai[x + (y * imageWidth)];
                buffer[index + 2] = (byte) ((temp >> 16) & 0xff);
                buffer[index + 1] = (byte) ((temp >> 8) & 0xff);
                buffer[index] = (byte) (temp & 0xff);
            }
        }

        writeInteger(40);
        writeInteger(imageWidth);
        writeInteger(imageHeight);
        writeWord(1);
        writeWord(24);
        writeInteger(0);
        writeInteger(0);
        writeInteger(0);
        writeInteger(0);
        writeInteger(0);
        writeInteger(0);

        for (int y = imageHeight - 1; y >= 0; y--) {
            int index = y * alignedWidth;
            for (int x = 0; x < alignedWidth; x += 2) {
                int k2 = x + index;
                writeWord(((buffer[k2 + 1] << 8) & 0xff00) | (buffer[k2] & 0xff));
            }
        }
    }

    private void writeBody(OutputStream outputstream) throws IOException
    {
        wmf.writeTo(outputstream);
    }

    protected void writeColor(Color color)
    {
        writeInteger(
            (color.getRed() & 0xff) | ((color.getGreen() << 8) & 0xff00)
                | ((color.getBlue() << 16) & 0xff0000));
    }

    private void writeHeader(OutputStream outputstream) throws IOException
    {
        writeRecordHeader(0, 0);
        outputWord(outputstream, 1);
        outputWord(outputstream, 9);
        outputWord(outputstream, 768);
        outputInteger(outputstream, getBodySize() + 9);
        outputWord(outputstream, handles.size());
        outputInteger(outputstream, maxObjectSize);
        outputWord(outputstream, 0);

    }

    protected void writeInteger(int i)
    {
        writeWord(lowWord(i));
        writeWord(highWord(i));
    }

    public void writeWMF(OutputStream outputstream) throws IOException
    {
        writeHeader(outputstream);
        writeBody(outputstream);
    }

    protected void writeWord(int i)
    {
        try {
            outputWord(wmf, i);
        } catch (IOException ioexception) {
            System.out.println(ioexception);
        }
    }

}

