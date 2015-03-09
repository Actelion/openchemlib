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
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;

public class TextDrawingObject extends AbstractDrawingObject {
	public static final String TYPE_STRING = "text";
	public static final float DEFAULT_SIZE = 9.0f;
	public static final int DEFAULT_STYLE = Font.PLAIN;

	private static final float LINE_SPACING = 1.4f;

	float				mSize,mZoomReferenceSize;
	String				mText;
	int					mStyle;
	boolean				mHilite;

	public TextDrawingObject() {
        this("",new Point2D.Float(),DEFAULT_SIZE, DEFAULT_STYLE);
		}

    public TextDrawingObject(String text, float x, float y)
    {
        this(text,x,y,DEFAULT_SIZE, DEFAULT_STYLE);
    }
    
    public TextDrawingObject(String text, float x, float y, float size, int style)
    {
        this(text,new Point2D.Float(x,y),size,style);
    }
    
	public TextDrawingObject(String descriptorDetail) {
		this();

		int index1 = 0;
		while (index1 != -1) {
            // text="this is a test"
            // 012345678901234567890
			int index2 = descriptorDetail.indexOf("=\"", index1);
			if (index2 == -1)
				break;	// should never happen
			String key = descriptorDetail.substring(index1+1, index2);
			index1 = descriptorDetail.indexOf("\"", index2+2);
			String value = (index1 == -1) ? descriptorDetail.substring(index2+1)
								: descriptorDetail.substring(index2+1, index1);
			if (key.equals("text"))
                setText(value);
//				mText = decodeText(value);
			else if (key.equals("x"))
                setX(value);
//				try { mPoint[0].x = Float.parseDouble(value); } catch (NumberFormatException nfe) {}
			else if (key.equals("y"))
                setY(value);
//				try { mPoint[0].y = Float.parseDouble(value); } catch (NumberFormatException nfe) {}
			else if (key.equals("size"))
                setSize(value);
//				try { mSize = Integer.parseInt(value); } catch (NumberFormatException nfe) {}
			else if (key.equals("style"))
                setStyle(value);
//				try { mStyle = Integer.parseInt(value); } catch (NumberFormatException nfe) {}
			}
		}
    protected void setText(String value)
    {
				mText = decodeText(value);
    }
    
    protected void setX(String value)
    {
        try { mPoint[0].x = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
    }

    protected void setY(String value) 
    {
            try { mPoint[0].y = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
    }
    
    protected void setSize(String value)
    {
        try { mSize = Float.parseFloat(value); } catch (NumberFormatException nfe) {}
    }
    
    protected void setStyle(String value)
    {
        try { mStyle = Integer.parseInt(value); } catch (NumberFormatException nfe) {}
        
    }



    private TextDrawingObject(String text, Point2D.Float pt,  float size, int style) {
		mText = text;
		mSize = size;
		mStyle = style;
		mPoint = new Point2D.Float[1];
		mPoint[0] = pt;        
    	}

    public String getTypeString() {
        return TYPE_STRING;   
    	}

    public String getDescriptorDetail() {
		StringBuilder detail = new StringBuilder();
		detail.append(" text=\""+encodeText(mText) + "\"");
		detail.append(" x=\""+mPoint[0].x + "\"");
		detail.append(" y=\""+mPoint[0].y + "\"");
		if (mSize != DEFAULT_SIZE)
			detail.append(String.format(" size=\"%.4f\"", new Float(mSize)));
		if (mStyle != DEFAULT_STYLE)
			detail.append(" style=\""+mStyle+ "\"");

		return detail.toString();
		}

	public AbstractDrawingObject clone() {
		TextDrawingObject duplicate = new TextDrawingObject();
		duplicate.setValues(mText, mSize, mStyle);
		duplicate.setCoordinates(mPoint[0].x, mPoint[0].y);
		duplicate.mIsSelected = this.mIsSelected;
		return duplicate;
		}

	public void setCoordinates(float x, float y) {
		mPoint[0].x = x;
		mPoint[0].y = y;
		}

	public void scale(float f) {
		super.scale(f);
		mSize *= f;
		}

	public void zoomAndRotateInit(float x, float y) {
		super.zoomAndRotateInit(x, y);
		mZoomReferenceSize = mSize;
		}

	public void zoomAndRotate(float zoom,float angle) {
		super.zoomAndRotate(zoom, angle);
		mSize = mZoomReferenceSize * zoom;
		}

	public void draw(Graphics g, DepictorTransformation t) {
		float size = (t == null) ? mSize : t.getScaling()*mSize;
		g.setFont(new Font("Helvetica", mStyle, (int)size));
		g.setColor(mIsSelected ? Color.red : Color.black);

		ArrayList<String> textList = getTextLineList();
		Rectangle2D.Float bounds = calculateBoundingRect(/* g.getFontMetrics() */ g.getFont(), textList);
		if (t != null)
			t.applyTo(bounds);

		for (int i=0; i<textList.size(); i++)
			g.drawString(textList.get(i),
						 (int)bounds.x,
						 (int)(bounds.y+1+size*5/6+size*LINE_SPACING*i));
		}

	public void draw2D(Graphics2D g, DepictorTransformation t) {
		draw(g, t);
		}

	private Rectangle2D.Float calculateBoundingRect(/* FontMetrics fm */Font font, ArrayList<String> textList) {
		FontRenderContext frc = new FontRenderContext(font.getTransform(), true, true);
				
		float maxWidth = 0;
		for (String text:textList) {
			if (text.length() != 0) {
				TextLayout textLayout = new TextLayout(text, font, frc);
				float width = (float)textLayout.getBounds().getWidth();
	//			float width = fm.stringWidth(text);
				if (maxWidth < width)
					maxWidth = width;
				}
			}
		float height = mSize * LINE_SPACING * (textList.size()-1) + mSize;
		return new Rectangle2D.Float(mPoint[0].x, mPoint[0].y-mSize/2, maxWidth, height);
		}

	public Rectangle2D.Float getBoundingRect() {
//		FontMetrics fm = Toolkit.getDefaultToolkit().getFontMetrics(new Font("Helvetica", mStyle, (int)mSize));
		ArrayList<String> textList = getTextLineList();
		return calculateBoundingRect(/* fm */new Font("Helvetica", mStyle, (int)mSize), textList);
		}

	private ArrayList<String> getTextLineList() {
		ArrayList<String> textList = new ArrayList<String>();
		int lineBreak = mText.indexOf('\n');
		if (lineBreak == -1) {
			textList.add(mText);
			}
		else {
			int textStart = 0;
			while (lineBreak != -1) {
				textList.add(mText.substring(textStart, lineBreak));
				textStart = lineBreak+1;
				lineBreak =  mText.indexOf('\n', textStart);
				}
			textList.add(mText.substring(textStart));
			}
		return textList;
		}

	public void hilite(Graphics g) {
		Rectangle2D.Float bounds = getBoundingRect();
		g.setColor(SELECTION_COLOR);
		g.fillRect((int)bounds.x, (int)bounds.y,
				   (int)bounds.width, (int)bounds.height);
		}

	public boolean checkHiliting(float x, float y) {
		Rectangle2D.Float bounds = getBoundingRect();
		mHilite = bounds.contains(x, y);
		return mHilite;
		}

	public boolean contains(float x, float y) {
		Rectangle2D.Float bounds = getBoundingRect();
		return bounds.contains(x, y);
		}

	public void clearHiliting() {
		mHilite = false;
		}

	public void setValues(String text, float size, int style) {
		mText = text;
		mSize = size;
		mStyle = style;
		}

	public String getText() {
		return mText;
		}

	public float getSize() {
		return mSize;
		}

	public int getStyle() {
		return mStyle;
		}

	private String encodeText(String source) {
		StringBuffer text = new StringBuffer();
		for (int i=0; i<source.length(); i++) {
			switch (source.charAt(i)) {
			case '&':
				text.append("&&");
				break;
			case '\t':
				text.append("&09");
				break;
			case '\n':
				text.append("&0A");
				break;
			case ' ':
				text.append("&20");
				break;
			default:
				text.append(source.charAt(i));
				break;
				}
			}
		return text.toString();
		}

	private String decodeText(String source) {
		int index = source.indexOf('&');
		if (index == -1)
			return source;

		int startIndex = 0;
		StringBuffer text = new StringBuffer();
		while (index != -1) {
			text.append(source.substring(startIndex, index));
			if (source.charAt(index+1) == '&') {
				text.append('&');
				startIndex = index+2;
				}
			else {
				int h = source.charAt(index+1);
				h = h - ((h < 'A') ? '0' : (h < 'a') ? 'A' : 'a');
				int l = source.charAt(index+2);
				l = l - ((l < 'A') ? '0' : (l < 'a') ? 'A' : 'a');
				text.append((char)(16*h+l));
				startIndex = index+3;
				}
			index = source.indexOf('&', startIndex);
			}
		text.append(source.substring(startIndex));
		return text.toString();
		}
	}
