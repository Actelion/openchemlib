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
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.font.TextLayout;
import java.awt.geom.*;
import java.awt.image.*;
import java.awt.image.renderable.RenderableImage;
import java.text.AttributedCharacterIterator;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;

public class WMFGraphics2D extends Graphics2D
{
    private BufferedImage img;
    private WMFGraphics wmfg;
    private Graphics2D g2D;
    private AffineTransform trans;
    private Stroke stroke;
    private Paint paint;
    private Color color;
    private Shape deviceclip;
    private boolean gdifontdrawing;
    private boolean gdipendrawing;
    private boolean gdipenwidthdrawing;
    private double flatness;

    public WMFGraphics2D(WMF wmf, int width, int height, Color fg, Color bg)
    {
        gdifontdrawing = true;
        gdipendrawing = true;
        gdipenwidthdrawing = true;
        flatness = 0.1D;
        if (fg == null)
            fg = Color.black;
        if (bg == null)
            bg = Color.white;
        wmfg = new WMFGraphics(wmf, width, height, fg, bg);
        img = new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB);
        g2D = (Graphics2D) img.getGraphics();
        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2D.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

        trans = new AffineTransform();
        stroke = g2D.getStroke();
        paint = fg;
        color = fg;
        deviceclip = null;
    }


    private WMFGraphics2D(WMFGraphics2D src)
    {
        wmfg = (WMFGraphics) src.wmfg.create();
        trans = (AffineTransform) src.trans.clone();
        stroke = src.stroke;
        paint = src.paint;
        color = src.color;
        deviceclip = src.deviceclip;
        img = new BufferedImage(src.img.getWidth(), src.img.getHeight(), BufferedImage.TYPE_INT_ARGB);
        flatness = src.flatness;
        gdifontdrawing = src.gdifontdrawing;
        gdipendrawing = src.gdipendrawing;
        gdipenwidthdrawing = src.gdipenwidthdrawing;
        g2D = (Graphics2D) src.img.getGraphics();
        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2D.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);

    }

    //    public void addRenderingHints(Map map)
    public void addRenderingHints(Map<?, ?> map)
    {
        g2D.addRenderingHints(map);
    }

    public void clearRect(int x, int y, int width, int height)
    {
        Paint paint1 = paint;
        setColor(getBackground());
        fillRect(x, y, width, height);
        setPaint(paint1);
    }

    public void clip(Shape shape)
    {
        if (deviceclip != null) {
            Area area = new Area(getClip());
            if (shape != null)
                area.intersect(new Area(shape));
            shape = area;
        }
        setClip(shape);
    }

    public void clipRect(int x, int y, int w, int h)
    {
        clip(new Rectangle(x, y, w, h));
    }

    public void copyArea(int x, int y, int width, int height, int dx, int dy)
    {
        g2D.copyArea(x, y, width, height, dx, dy);
    }

    public Graphics create()
    {
        WMFGraphics2D wmfgraphics2d = new WMFGraphics2D(this);
        return wmfgraphics2d;
    }

    public void dispose()
    {
        wmfg.dispose();
        g2D.dispose();
        img.flush();
    }

    private void draw(Shape shape, int i)
    {
        Shape area = null;
        shape = trans.createTransformedShape(shape);
        if (deviceclip != null) {
            area = deviceclip;
        }
        AffineTransform affinetransform = new AffineTransform();
        PathIterator pathiterator = shape.getPathIterator(affinetransform, flatness);
        float af[] = new float[6];
        Vector vector = new Vector();
        Polygon polygon = null;
        float af1[] = new float[6];
        for (; !pathiterator.isDone(); pathiterator.next()) {
            int segment = pathiterator.currentSegment(af);
            switch (segment) {

                case PathIterator.SEG_MOVETO: // '\0'
                    if (polygon != null)
                        vector.add(polygon);
                    polygon = new Polygon();
                    polygon.addPoint((int) ((double) af[0] + 0.5D), (int) ((double) af[1] + 0.5D));
                    System.arraycopy(af, 0, af1, 0, af.length);
                    break;

                case PathIterator.SEG_CUBICTO: // '\003'
                    if (polygon != null)
                        polygon.addPoint((int) ((double) af[2] + 0.5D), (int) ((double) af[3] + 0.5D));
                    break;

                case PathIterator.SEG_QUADTO: // '\002'
                    if (polygon != null)
                        polygon.addPoint((int) ((double) af[4] + 0.5D), (int) ((double) af[5] + 0.5D));
                    break;

                case PathIterator.SEG_LINETO: // '\001'
                    if (polygon != null)
                        polygon.addPoint((int) ((double) af[0] + 0.5D), (int) ((double) af[1] + 0.5D));
                    break;

                case PathIterator.SEG_CLOSE: // '\004'
                    if (polygon != null) {
                        polygon.addPoint((int) af1[0], (int) af1[1]);
                        vector.add(polygon);
                        polygon = null;
                    }
                    break;
                default:
                    break;
            }
        }

        if (polygon != null)
            vector.add(polygon);
        Polygon apolygon[] = new Polygon[vector.size()];
        int k = 0;
        for (Iterator iterator = vector.iterator(); iterator.hasNext(); ) {
            apolygon[k] = (Polygon) iterator.next();
            k++;
        }

        int penStyle = wmfg.getPenStyle();
        int brushFillStyle = wmfg.getBrushFillStyle();
        wmfg.setClip(area);
        wmfg.setColor(getColor());
        wmfg.setPenStyle(WMF.PS_SOLID);
        if (gdipenwidthdrawing)
            wmfg.setPenWidth(i);
        else
            wmfg.setPenWidth(0);
        wmfg.setGDIHollowBrush();
        for (int l = 0; l < apolygon.length; l++)
            wmfg.drawPolyline(apolygon[l].xpoints, apolygon[l].ypoints, apolygon[l].npoints);

        wmfg.setBrushFillStyle(brushFillStyle);
        wmfg.setPenStyle(penStyle);
    }


    public void fill(Shape shape)
    {
        Shape save = deviceclip;
        wmfg.setClip(save);
        Shape theShape = trans.createTransformedShape(shape);
        Rectangle2D area = theShape.getBounds2D();
        if (deviceclip != null) {
            area = area.createIntersection(deviceclip.getBounds2D());
        }
        if (!(getPaint() instanceof Color)) {
            Rectangle bounds = area.getBounds();
            PaintContext context = paint.createContext(img.getColorModel(), bounds, shape.getBounds(), trans, getRenderingHints());
            Raster raster = context.getRaster(bounds.x, bounds.y, bounds.width, bounds.height);
            BufferedImage bufferedimage = new BufferedImage(context.getColorModel(), raster.createCompatibleWritableRaster(), false, null);
            bufferedimage.setData(raster);
            BufferedImage bufferedimage1 = new BufferedImage(bounds.width, bounds.height, 2);
            Graphics2D graphics2d = (Graphics2D) bufferedimage1.getGraphics();
            graphics2d.addRenderingHints(getRenderingHints());
            graphics2d.setTransform(AffineTransform.getTranslateInstance(-bounds.x, -bounds.y));
            graphics2d.setClip(area);
            graphics2d.drawImage(bufferedimage, bounds.x, bounds.y, null);
            wmfg.drawImage(bufferedimage1, bounds.x, bounds.y, bounds.width, bounds.height, null);
            return;
        }
        AffineTransform affinetransform = new AffineTransform();
        PathIterator pathiterator = theShape.getPathIterator(affinetransform, flatness);
        if (pathiterator.getWindingRule() == 0)
            wmfg.setPolyFillMode(WMF.ALTERNATE);
        else
            wmfg.setPolyFillMode(WMF.WINDING);
        int i = wmfg.getPenStyle();
        wmfg.setPenStyle(WMF.PS_NULL);
        wmfg.setColor(getColor());
        float af[] = new float[6];
        Vector vector = new Vector();
        Polygon polygon = null;
        for (; !pathiterator.isDone(); pathiterator.next()) {
            int j = pathiterator.currentSegment(af);
            switch (j) {
                case PathIterator.SEG_MOVETO:
                    if (polygon != null)
                        vector.add(polygon);
                    polygon = new Polygon();
                    polygon.addPoint((int) ((double) af[0] + 0.5D), (int) ((double) af[1] + 0.5D));
                    break;

                case PathIterator.SEG_CUBICTO:
                    if (polygon != null)
                        polygon.addPoint((int) ((double) af[2] + 0.5D), (int) ((double) af[3] + 0.5D));
                    break;

                case PathIterator.SEG_QUADTO:
                    if (polygon != null)
                        polygon.addPoint((int) ((double) af[4] + 0.5D), (int) ((double) af[5] + 0.5D));
                    break;

                case PathIterator.SEG_LINETO:
                    if (polygon != null)
                        polygon.addPoint((int) ((double) af[0] + 0.5D), (int) ((double) af[1] + 0.5D));
                    break;

                case PathIterator.SEG_CLOSE:
                    if (polygon != null) {
                        vector.add(polygon);
                        polygon = null;
                    }
                    break;

                default:
                    break;

            }
        }

        if (polygon != null)
            vector.add(polygon);
        Polygon apolygon[] = new Polygon[vector.size()];
        int k = 0;
        for (Iterator iterator = vector.iterator(); iterator.hasNext(); ) {
            apolygon[k] = (Polygon) iterator.next();
            k++;
        }

        boolean flag1 = false;
        if (shape instanceof Arc2D)
            flag1 = true;
        if (!flag1) {
            wmfg.GDIPolyPolygon(apolygon);
        } else {
            for (int l = 0; l < apolygon.length; l++)
                wmfg.fillPolygon(apolygon[l]);

        }
        wmfg.setPenStyle(i);
    }

    public void draw(Shape shape)
    {
        if (gdipendrawing && (getPaint() instanceof Color)
            && (stroke instanceof BasicStroke) &&
            (((BasicStroke) stroke).getDashArray() == null || ((BasicStroke) stroke).getDashArray().length == 0)) {
            draw(shape, (int) ((BasicStroke) stroke).getLineWidth());
        }
    }

    public void drawArc(int x, int y, int width, int height, int startAngle, int arcAngle)
    {
        draw(new java.awt.geom.Arc2D.Float(x, y, width, height, startAngle, arcAngle, Arc2D.OPEN));
    }

    public void drawGlyphVector(GlyphVector glyphvector, float f, float f1)
    {
        fill(glyphvector.getOutline(f, f1));
    }

    public boolean drawImage(Image image, int dx1, int dy1, int dx2, int dy2, int sx1, int sy1,
                             int sx2, int sy2, Color color1, ImageObserver imageobserver)
    {
        int ai[] = {
            dx1,
            dy1,
            dx2,
            dy2,
            sx1,
            sy1,
            sx2,
            sy2
        };
        Rectangle rectangle = (new java.awt.geom.Line2D.Float(dx1, dy1, dx2, dy2)).getBounds();
        Shape shape = trans.createTransformedShape(rectangle);
        Rectangle rectangle1 = shape.getBounds();
        Image image1 = transformImage(image, ai, rectangle1, imageobserver, color1);
        wmfg.setClip(null);
        wmfg.drawImage(image1, rectangle1.x, rectangle1.y, rectangle1.width, rectangle1.height, imageobserver);
        return true;
    }

    public boolean drawImage(Image image, int dx1, int dy1, int dx2, int dy2, int sx1, int sy1,
                             int sx2, int sy2, ImageObserver imageobserver)
    {
        return drawImage(image, dx1, dy1, dx2, dy2, sx1, sy1, sx2, sy2, null, imageobserver);
    }

    public boolean drawImage(Image image, int x, int y, int width, int height, Color color1, ImageObserver imageobserver)
    {
        Rectangle rc = new Rectangle(x, y, width, height);
        Shape shape = trans.createTransformedShape(rc);
        Rectangle bound = shape.getBounds();
        Image image1 = transformImage(image, rc, bound, imageobserver, color1);
        wmfg.setClip(null);
        wmfg.drawImage(image1, bound.x, bound.y, bound.width, bound.height, imageobserver);
        return true;
    }

    public boolean drawImage(Image image, int x, int y, int width, int height, ImageObserver imageobserver)
    {
        return drawImage(image, x, y, width, height, null, imageobserver);
    }

    public boolean drawImage(Image image, int x, int y, Color color1, ImageObserver imageobserver)
    {
        return drawImage(image, x, y, image.getWidth(imageobserver), image.getHeight(imageobserver), color1, imageobserver);
    }

    public boolean drawImage(Image image, int x, int y, ImageObserver imageobserver)
    {
        return drawImage(image, x, y, image.getWidth(imageobserver), image.getHeight(imageobserver), imageobserver);
    }

    public boolean drawImage(Image image, AffineTransform affinetransform, ImageObserver imageobserver)
    {
        AffineTransform affinetransform1 = (AffineTransform) trans.clone();
        trans.concatenate(affinetransform);
        drawImage(image, 0, 0, imageobserver);
        trans = affinetransform1;
        return true;
    }

    public void drawImage(BufferedImage bufferedimage, BufferedImageOp bufferedimageop, int i, int j)
    {
        BufferedImage bufferedimage1 = bufferedimageop.filter(bufferedimage, null);
        drawImage(((Image) (bufferedimage1)), new AffineTransform(1.0F, 0.0F, 0.0F, 1.0F, i, j), null);
    }

    public void drawLine(int x1, int y1, int x2, int y2)
    {
        draw(new GeneralPath(new java.awt.geom.Line2D.Float(x1, y1, x2, y2)));
    }

    public void drawOval(int x, int y, int width, int height)
    {
        draw(new java.awt.geom.Ellipse2D.Float(x, y, width, height));
    }

    public void drawPolygon(int xPoints[], int yPoints[], int i)
    {
        draw(new Polygon(xPoints, yPoints, i));
    }

    public void drawPolyline(int xPoints[], int yPoints[], int i)
    {
        if (i > 0) {
            GeneralPath generalpath = new GeneralPath();
            generalpath.moveTo(xPoints[0], yPoints[0]);
            for (int j = 1; j < i; j++)
                generalpath.lineTo(xPoints[j], yPoints[j]);

            draw(generalpath);
        }
    }

    public void drawRect(int x, int y, int width, int height)
    {
        Rectangle rectangle = new Rectangle(x, y, width, height);
        draw(rectangle);
    }

    public void drawRenderableImage(RenderableImage img, AffineTransform transform)
    {
        drawRenderedImage(img.createDefaultRendering(), transform);
    }

    public void drawRenderedImage(RenderedImage img, AffineTransform transform)
    {
        BufferedImage bufferedimage =
            new BufferedImage(img.getColorModel(), img.getData().createCompatibleWritableRaster(), false, null);
        bufferedimage.setData(img.getData());
        drawImage(bufferedimage, transform, null);
    }

    public void drawRoundRect(int x, int y, int width, int height, int arcWidth, int arcHeight)
    {
        draw(new java.awt.geom.RoundRectangle2D.Float(x, y, width, height, arcWidth, arcHeight));
    }

    public void drawString(String s, float x, float y)
    {
        if (isGDIFontDrawing() && (getPaint() instanceof Color)) {
            boolean noClipping = deviceclip == null;
            if (!noClipping) {
                GlyphVector glyphvector = getFont().createGlyphVector(getFontRenderContext(), s);
                java.awt.geom.Rectangle2D rectangle2d = glyphvector.getOutline(x, y).getBounds2D();
                noClipping = deviceclip.contains(trans.createTransformedShape(rectangle2d).getBounds2D());
            }
            if (noClipping) {
                boolean isIdent = trans.isIdentity() && getFont().getTransform().isIdentity();
                double scaleX = 1.0D;
                double scale = 0.0D;
                if (!isIdent) {
                    AffineTransform affinetransform = getFont().getTransform();
                    affinetransform.preConcatenate(trans);
                    scaleX = affinetransform.getScaleX();
                    isIdent = affinetransform.getShearX() == -affinetransform.getShearY() && affinetransform.getScaleX() == affinetransform.getScaleY();
                    if (isIdent && (affinetransform.getShearX() != 0.0D || affinetransform.getScaleX() < 0.0D)) {
                        scaleX = Math.sqrt(affinetransform.getScaleX() * affinetransform.getScaleX() + affinetransform.getShearX() * affinetransform.getShearX());
                        scale = Math.acos(affinetransform.getScaleX() / scaleX);
                        if (affinetransform.getShearX() > 0.0D)
                            scale = -scale;
                    }
                }
                if (isIdent) {
//                    boolean flag2 = true;
//                    if (flag2)
                    {
                        wmfg.setColor(getColor());
                        float f2 = (float) ((double) getFont().getSize2D() * scaleX);
                        Font font = getFont();
                        wmfg.setFont(getFont().deriveFont(f2));
                        if (scale != 0.0D)
                            wmfg.setFontEscapement((int) ((-scale * 1800D) / Math.PI));
                        java.awt.geom.Point2D.Double double1 = new java.awt.geom.Point2D.Double(
                            (double) x + getFont().getTransform().getTranslateX(), (double) y + getFont().getTransform().getTranslateY());
                        trans.transform(double1, double1);
                        wmfg.drawString(s, (int) double1.getX(), (int) double1.getY());
                        if (scale != 0.0D)
                            wmfg.setFontEscapement(0);
                        wmfg.setFont(font);
                        return;
                    }
                }
            }
        }
        drawGlyphVector(getFont().createGlyphVector(getFontRenderContext(), s), x, y);
    }

    public void drawString(String s, int x, int y)
    {
        drawString(s, (float) x, (float) y);
    }

    public void drawString(AttributedCharacterIterator attributedcharacteriterator, float f, float f1)
    {
        TextLayout textlayout = new TextLayout(attributedcharacteriterator, getFontRenderContext());
        Paint paint1 = getPaint();
        setColor(getColor());
        fill(textlayout.getOutline(AffineTransform.getTranslateInstance(f, f1)));
        setPaint(paint1);
    }

    public void drawString(AttributedCharacterIterator it, int i, int j)
    {
        drawString(it, i, j);
    }

//    public void fill(Shape shape)
//    {
//        fill(shape, false);
//    }

    public void fillArc(int x, int y, int width, int height, int startAngle, int arcAngle)
    {
        fill(new java.awt.geom.Arc2D.Float(x, y, width, height, startAngle, arcAngle, Arc2D.PIE));
    }

    public void fillOval(int x, int y, int width, int height)
    {
        fill(new java.awt.geom.Ellipse2D.Float(x, y, width, height));
    }

    public void fillPolygon(int xPoints[], int yPoints[], int i)
    {
        fill(new Polygon(xPoints, yPoints, i));
    }

    public void fillRect(int x, int y, int width, int height)
    {
        fill(new Rectangle(x, y, width, height));
    }

    public void fillRoundRect(int x, int y, int width, int height, int arcWidth, int arcHeight)
    {
        fill(new java.awt.geom.RoundRectangle2D.Float(x, y, width, height, arcWidth, arcHeight));
    }

    public Color getBackground()
    {
        return wmfg.getBackground();
    }

    public Shape getClip()
    {
        try {
            return trans.createInverse().createTransformedShape(deviceclip);
        } catch (Exception _ex) {
            return null;
        }
    }

    public Rectangle getClipBounds()
    {
        if (deviceclip != null)
            return getClip().getBounds();
        else
            return null;
    }

    public Color getColor()
    {
        return color;
    }

    public Composite getComposite()
    {
        return g2D.getComposite();
    }

    public GraphicsConfiguration getDeviceConfiguration()
    {
        return g2D.getDeviceConfiguration();
    }

    public double getFlatness()
    {
        return flatness;
    }

    public Font getFont()
    {
        return wmfg.getFont();
    }

    public FontMetrics getFontMetrics(Font font)
    {
        return wmfg.getFontMetrics(font);
    }

    public FontRenderContext getFontRenderContext()
    {
        g2D.setTransform(trans);
        return g2D.getFontRenderContext();
    }

    public Paint getPaint()
    {
        return paint;
    }

    public Object getRenderingHint(java.awt.RenderingHints.Key key)
    {
        return g2D.getRenderingHint(key);
    }

    public RenderingHints getRenderingHints()
    {
        return g2D.getRenderingHints();
    }

    public Stroke getStroke()
    {
        return stroke;
    }

    public AffineTransform getTransform()
    {
        return (AffineTransform) trans.clone();
    }

    public boolean hit(Rectangle rectangle, Shape shape, boolean flag)
    {
        g2D.setTransform(trans);
        g2D.setStroke(getStroke());
        g2D.setClip(getClip());
        return g2D.hit(rectangle, shape, flag);
    }

    public boolean isGDIFontDrawing()
    {
        return gdifontdrawing;
    }

    public boolean isGDIPenDrawing()
    {
        return gdipendrawing;
    }

    public boolean isGDIPenWidthDrawing()
    {
        return gdipenwidthdrawing;
    }

    public void rotate(double theta)
    {
        trans.rotate(theta);
    }

    public void rotate(double theta, double x, double y)
    {
        trans.rotate(theta, x, y);
    }

    public void scale(double sx, double sy)
    {
        trans.scale(sx, sy);
    }

    public void setBackground(Color color)
    {
        wmfg.setBackground(color);
    }

    public void setClip(int x, int y, int width, int height)
    {
        setClip(((Shape) (new Rectangle(x, y, width, height))));
    }

    public void setClip(Shape shape)
    {
//        System.out.println("SetClip " + (shape != null ? shape.getBounds() : null));
        deviceclip = trans.createTransformedShape(shape);
    }

    public void setColor(Color color1)
    {
        setPaint(color1);
    }

    public void setComposite(Composite composite)
    {
        g2D.setComposite(composite);
    }

    public void setFlatness(double d)
    {
        flatness = d;
    }

    public void setFont(Font font)
    {
        wmfg.setFont(font);
    }

    public void setGDIFontDrawing(boolean flag)
    {
        gdifontdrawing = flag;
    }

    public void setGDIPenDrawing(boolean flag)
    {
        gdipendrawing = flag;
    }

    public void setGDIPenWidthDrawing(boolean flag)
    {
        gdipenwidthdrawing = flag;
    }

    public void setPaint(Paint p)
    {
        if (p != null) {
            paint = p;
            if (p instanceof Color)
                color = (Color) p;
        }
    }

    public void setPaintMode()
    {
        wmfg.setPaintMode();
    }

    public void setRenderingHint(java.awt.RenderingHints.Key key, Object obj)
    {
        g2D.setRenderingHint(key, obj);
    }

    public void setRenderingHints(Map<?, ?> map)
    {
        g2D.setRenderingHints(map);
    }

    public void setStroke(Stroke stroke1)
    {
        stroke = stroke1;
    }

    public void setTransform(AffineTransform affinetransform)
    {
        trans = (AffineTransform) affinetransform.clone();
    }

    public void setXORMode(Color color1)
    {
        wmfg.setXORMode(color1);
    }

    public void shear(double shx, double shy)
    {
        trans.shear(shx, shy);
    }

    public void transform(AffineTransform affinetransform)
    {
        trans.concatenate(affinetransform);
    }

    private Image transformImage(Image image,
                                 Rectangle rectangle, Rectangle rectangle1,
                                 ImageObserver imageobserver, Color bgcolor)
    {
        if (trans.isIdentity() && getClip() == null && bgcolor == null)
            return image;
        if (deviceclip != null) {
            Area area = new Area(deviceclip);
            area.intersect(new Area(trans.createTransformedShape(rectangle)));
            rectangle1.setBounds(area.getBounds());
        }
        BufferedImage bufferedimage = new BufferedImage(rectangle1.width, rectangle1.height, 2);
        Graphics2D graphics2d = (Graphics2D) bufferedimage.getGraphics();
        graphics2d.addRenderingHints(getRenderingHints());
        graphics2d.translate(-rectangle1.x, -rectangle1.y);
        graphics2d.transform((AffineTransform) getTransform().clone());
        graphics2d.setClip(getClip());
        if (bgcolor != null)
            graphics2d.drawImage(image, rectangle.x, rectangle.y, rectangle.width, rectangle.height, bgcolor, imageobserver);
        else
            graphics2d.drawImage(image, rectangle.x, rectangle.y, rectangle.width, rectangle.height, imageobserver);
        graphics2d.dispose();
        return bufferedimage;
    }

    private Image transformImage(Image image, int ai[], Rectangle rectangle, ImageObserver imageobserver, Color bgcolor)
    {
        BufferedImage bufferedimage = new BufferedImage(rectangle.width, rectangle.height, 2);
        Graphics2D graphics2d = (Graphics2D) bufferedimage.getGraphics();
        graphics2d.addRenderingHints(getRenderingHints());
        graphics2d.translate(-rectangle.x, -rectangle.y);
        graphics2d.transform((AffineTransform) getTransform().clone());
        graphics2d.setClip(getClip());
        if (bgcolor != null)
            graphics2d.drawImage(image, ai[0], ai[1], ai[2], ai[3], ai[4], ai[5], ai[6], ai[7], bgcolor, imageobserver);
        else
            graphics2d.drawImage(image, ai[0], ai[1], ai[2], ai[3], ai[4], ai[5], ai[6], ai[7], imageobserver);
        return bufferedimage;
    }

    public void translate(double tx, double ty)
    {
        trans.translate(tx, ty);
    }

    public void translate(int x, int y)
    {
        trans.translate(x, y);
    }

}
