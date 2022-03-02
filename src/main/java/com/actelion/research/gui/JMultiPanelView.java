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
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import com.actelion.research.util.DoubleFormat;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.util.ArrayList;

public class JMultiPanelView extends JPanel implements MultiPanelDragListener {
    static final long serialVersionUID = 0x20070314;

    public static final String VIEW_HEIGHT = "height";
	private ArrayList<MultiPanelObject>	mPanelList;
	private TableLayout	mLayout;
	private int			mDragTitle,mDragStart;
	private int[]		mDragStartPanelHeight;

	public JMultiPanelView() {
		double[][] size = { { TableLayout.FILL }, {} };
		mLayout = new TableLayout(size);
		setLayout(mLayout);
        setOpaque(false);
		mPanelList = new ArrayList<>();
		}

	public void add(JComponent c, String title) {
		add(c, title, mPanelList.size());
		}

	public void add(JComponent c, String title, int position) {
		double height = suggestHeight();
		JMultiPanelTitle titlePanel = new JMultiPanelTitle(this, title);
		if (position == 0) {
			titlePanel.setDragEnabled(false);
			for (int i=0; i<mPanelList.size(); i++)
				mPanelList.get(i).titlePanel.setDragEnabled(true);
			}

		mLayout.insertRow(2*position, JMultiPanelTitle.HEIGHT);
		mLayout.insertRow(2*position+1, height);

		super.add(titlePanel, "0, "+(2*position));
		super.add(c, "0, "+(2*position+1));

		mPanelList.add(position, new MultiPanelObject(c, titlePanel, title));
		validateHeights();
		}

	public int getViewCount() {
		return mPanelList.size();
		}

	public JComponent getView(int index) {
		return mPanelList.get(index).panel;
		}

	public void setTitle(int position, String title) {
		mPanelList.get(position).titlePanel.setTitle(title);
		}

	public void remove(JComponent c) {
		super.remove(c);
		for (int i=0; i<mPanelList.size(); i++) {
			if (mPanelList.get(i).panel == c) {
			    super.remove(mPanelList.get(i).titlePanel);
				mPanelList.remove(i);
				mLayout.deleteRow(2*i);
				mLayout.deleteRow(2*i);
				validateHeights();
				return;
				}
			}
		}

	public String getProperties() {
		StringBuffer properties = new StringBuffer();
		for (int i=0; i<mPanelList.size(); i++) {
			MultiPanelObject o = mPanelList.get(i);
			if (mLayout.getRow(2*i+1) > 0.0) {
				if (properties.length() > 0)
					properties.append(';');
				properties.append(VIEW_HEIGHT+"["+o.title+"]="+DoubleFormat.toString(mLayout.getRow(2*i+1)));
				}
			}
		return properties.toString();
		}

	public void setProperties(String properties) {
		for (int i=0; i<mPanelList.size(); i++)
			mLayout.setRow(2*i+1, 0.0);

		if (properties.startsWith(VIEW_HEIGHT+"[")) {
			boolean[] isAssigned = new boolean[mPanelList.size()];
			int titleStart = VIEW_HEIGHT.length()+1;
			while (titleStart < properties.length()) {
				int titleEnd = properties.indexOf("]=", titleStart);
				if (titleEnd == -1)
					break;

				String title = properties.substring(titleStart, titleEnd);
				int numberEnd = properties.indexOf(";", titleEnd);
				if (numberEnd == -1)
					numberEnd = properties.length();

				for (int i=0; i<mPanelList.size(); i++) {
					if (!isAssigned[i] && mPanelList.get(i).titleMatches(title)) {
						isAssigned[i] = true;
						try {
							double height = Double.parseDouble(properties.substring(titleEnd+2, numberEnd));
							mLayout.setRow(2*i+1, Math.min(0.999999, height));
							}
						catch (NumberFormatException nfe) {}
						}
					}
				titleStart = numberEnd + VIEW_HEIGHT.length()+2;
				}
			}

		validateHeights();
		}

	public void dragStarted(int titleY, JMultiPanelTitle t) {
		mDragTitle = 0;	// don't act on first title dragging
		while (mDragTitle<mPanelList.size() && mPanelList.get(mDragTitle).titlePanel != t)
			mDragTitle++;

		if (mDragTitle != 0) {
			mDragStartPanelHeight = getPanelHeights();
			mDragStart = getAbsoluteY(titleY, mDragStartPanelHeight);
			}
		}

	public void dragContinued(int titleY) {
		if (mDragTitle != 0) {
			int[] panelHeight = getPanelHeights();
			int dy = getAbsoluteY(titleY, panelHeight) - mDragStart;
			int[] newPanelHeight = new int[mPanelList.size()];
			for (int i=0; i<mPanelList.size(); i++)
				newPanelHeight[i] = mDragStartPanelHeight[i];

			if (dy < 0) {
				int heightUsed = 0;
				for (int panel=mDragTitle-1; panel>=0; panel--) {
					if (-dy <= mDragStartPanelHeight[panel] - 6) {
						newPanelHeight[panel] += dy;
						heightUsed -= dy;
						break;
						}
					else if (-dy <= mDragStartPanelHeight[panel]) {
						newPanelHeight[panel] = 0;
						heightUsed += mDragStartPanelHeight[panel];
						break;
						}
					dy += mDragStartPanelHeight[panel];
					heightUsed += mDragStartPanelHeight[panel];
					newPanelHeight[panel] = 0;
					}
				newPanelHeight[mDragTitle] += heightUsed;
				}
			else if (dy > 0) {
				int heightUsed = 0;
				for (int panel=mDragTitle; panel<mPanelList.size(); panel++) {
					if (dy <= mDragStartPanelHeight[panel] - 6) {
						newPanelHeight[panel] -= dy;
						heightUsed += dy;
						break;
						}
					else if (dy <= mDragStartPanelHeight[panel]) {
						newPanelHeight[panel] = 0;
						heightUsed += mDragStartPanelHeight[panel];
						break;
						}
					dy -= mDragStartPanelHeight[panel];
					heightUsed += mDragStartPanelHeight[panel];
					newPanelHeight[panel] = 0;
					}
				newPanelHeight[mDragTitle-1] += heightUsed;
				}

			int totalHeight = 0;
			for (int panel=0; panel<mPanelList.size(); panel++)
				totalHeight += panelHeight[panel];
			for (int panel=0; panel<mPanelList.size(); panel++) {
				if (newPanelHeight[panel] != panelHeight[panel])
					mLayout.setRow(panel*2+1, Math.min(0.999999, (double)newPanelHeight[panel]/(double)totalHeight));
				revalidate();
				}
			}
		}

	private int[] getPanelHeights() {
		int[] size = new int[mPanelList.size()];
		for (int i=0; i<mPanelList.size(); i++)
			size[i] = mPanelList.get(i).panel.getHeight();
		return size;
		}

	private int getAbsoluteY(int y, int[] panelHeight) {
		int absoluteY = y;
		for (int i=0; i<mDragTitle; i++)
			absoluteY += JMultiPanelTitle.HEIGHT + panelHeight[i];
		return absoluteY;
		}

	private double suggestHeight() {
		int count = 0;
		double totalHeight = 0.0;
		for (int i=0; i<mPanelList.size(); i++) {
			if (mLayout.getRow(2*i+1) != 0.0) {
				totalHeight += mLayout.getRow(2*i+1);
				count++;
				}
			}
		return (count == 0) ? 1.0 : totalHeight / count;
		}

	private void validateHeights() {
		double totalHeight = 0.0;
		for (int i=0; i<mPanelList.size(); i++)
			totalHeight += mLayout.getRow(2*i+1);

		if (totalHeight != 0.0) {
			double factor = 1.0 / totalHeight;
			for (int i=0; i<mPanelList.size(); i++)
				mLayout.setRow(2*i+1, Math.min(0.999999, factor * mLayout.getRow(2*i+1)));
			}
		revalidate();
		}
	}

class MultiPanelObject {
	public JComponent		panel;
	public JMultiPanelTitle	titlePanel;
	public String			title;

	public MultiPanelObject(JComponent c, JMultiPanelTitle titlePanel, String title) {
		this.panel = c;
		this.titlePanel = titlePanel;
		this.title = title;
		}

	/**
	 * The given name matches the title of this panel object if<br>
	 *     - name equals title
	 *     - if title starts with name followed by ' ('+someString+')'
	 * @param name
	 * @return
	 */
	public boolean titleMatches(String name) {
		if (title.equals(name))
			return true;
		if (!title.startsWith(name))
			return false;
		int nameLength = name.length();
		return (title.length() > nameLength+3 && title.substring(nameLength, nameLength+2).equals(" (") && title.endsWith(")"));
		}
	}
