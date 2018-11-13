package com.actelion.research.gui;

import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.LayoutManager;

/**
 *  A vertical flow layout is similar to a flow layuot but it layouts the
 *  components vertically instead of horizontally.
 *
 *@author     vassilidzuba
 *@created    24 mars 2001
 */
public class VerticalFlowLayout implements LayoutManager, java.io.Serializable {
    static final long serialVersionUID = 0x20110413;

	int _halign;
	int _valign;
	int _hgap;
	int _vgap;
	boolean _hfill;

	/**
	 *  Description of the Field
	 */
	public final static int TOP = 0;
	/**
	 *  Description of the Field
	 */
	public final static int CENTER = 1;
	/**
	 *  Description of the Field
	 */
	public final static int BOTTOM = 2;
	/**
	 *  Description of the Field
	 */
	public final static int LEFT = 3;
	/**
	 *  Description of the Field
	 */
	public final static int RIGHT = 4;


	/**
	 *  Constructor for the VerticalFlowLayout object
	 */
	public VerticalFlowLayout() {
		this(LEFT, TOP, 5, 5, true);
	}


	/**
	 *  Constructor for the VerticalFlowLayout object
	 *
	 *@param  halign  Description of Parameter
	 *@param  valign  Description of Parameter
	 */
	public VerticalFlowLayout(int halign, int valign) {
		this(halign, valign, 5, 5, true);
	}


	/**
	 *  Constructor for the VerticalFlowLayout object
	 *
	 *@param  halign  Description of Parameter
	 *@param  valign  Description of Parameter
	 *@param  hgap    Description of Parameter
	 *@param  vgap    Description of Parameter
	 */
	public VerticalFlowLayout(int halign, int valign, int hgap, int vgap, boolean hfill) {
		_hgap = hgap;
		_vgap = vgap;
		_hfill = hfill;
		setAlignment(halign, valign);
	}


	/**
	 *  Sets the Alignment attribute of the VerticalFlowLayout object
	 *
	 *@param  halign  The new Alignment value
	 *@param  valign  The new Alignment value
	 */
	public void setAlignment(int halign, int valign) {
		_halign = halign;
		_valign = valign;
	}


	/**
	 *  Sets the Hgap attribute of the VerticalFlowLayout object
	 *
	 *@param  hgap  The new Hgap value
	 */
	public void setHgap(int hgap) {
		_hgap = hgap;
	}


	/**
	 *  Sets the Vgap attribute of the VerticalFlowLayout object
	 *
	 *@param  vgap  The new Vgap value
	 */
	public void setVgap(int vgap) {
		_vgap = vgap;
	}


	/**
	 *  Gets the Halignment attribute of the VerticalFlowLayout object
	 *
	 *@return    The Halignment value
	 */
	public int getHalignment() {
		return _halign;
	}


	/**
	 *  Gets the Valignment attribute of the VerticalFlowLayout object
	 *
	 *@return    The Valignment value
	 */
	public int getValignment() {
		return _valign;
	}


	/**
	 *  Gets the Hgap attribute of the VerticalFlowLayout object
	 *
	 *@return    The Hgap value
	 */
	public int getHgap() {
		return _hgap;
	}


	/**
	 *  Gets the Vgap attribute of the VerticalFlowLayout object
	 *
	 *@return    The Vgap value
	 */
	public int getVgap() {
		return _vgap;
	}


	/**
	 *  Sets the Hfill attribute of the VerticalFlowLayout object
	 *
	 * @param hfill
	 */
	public void setHorizontalFill(boolean hfill) {
		_hfill = hfill;
	}

	/**
	 *  Gets the Hfill attribute of the VerticalFlowLayout object
	 *
	 *@return    The Hfill value
	 */
	public boolean getHorizontalFill() {
		return _hfill;
	}


	/**
	 *  Adds a feature to the LayoutComponent attribute of the VerticalFlowLayout
	 *  object
	 *
	 *@param  name  The feature to be added to the LayoutComponent attribute
	 *@param  comp  The feature to be added to the LayoutComponent attribute
	 */
	public void addLayoutComponent(String name, Component comp) {
	}


	/**
	 *  Description of the Method
	 *
	 *@param  comp  Description of Parameter
	 */
	public void removeLayoutComponent(Component comp) {
	}


	/**
	 *  Description of the Method
	 *
	 *@param  target  Description of Parameter
	 *@return         Description of the Returned Value
	 */
	public Dimension preferredLayoutSize(Container target) {
		synchronized (target.getTreeLock()) {
			Dimension dim = new Dimension(0, 0);
			int nmembers = target.getComponentCount();
			boolean firstVisibleComponent = true;

			for (int ii = 0; ii < nmembers; ii++) {
				Component m = target.getComponent(ii);
				if (m.isVisible()) {
					Dimension d = m.getPreferredSize();
					dim.width = Math.max(dim.width, d.width);
					if (firstVisibleComponent) {
						firstVisibleComponent = false;
					}
					else {
						dim.height += _vgap;
					}
					dim.height += d.height;
				}
			}
			Insets insets = target.getInsets();
			dim.width += insets.left + insets.right + _hgap * 2;
			dim.height += insets.top + insets.bottom + _vgap * 2;
			return dim;
		}
	}


	/**
	 *  Description of the Method
	 *
	 *@param  target  Description of Parameter
	 *@return         Description of the Returned Value
	 */
	public Dimension minimumLayoutSize(Container target) {
		synchronized (target.getTreeLock()) {
			Dimension dim = new Dimension(0, 0);
			int nmembers = target.getComponentCount();
			boolean firstVisibleComponent = true;

			for (int ii = 0; ii < nmembers; ii++) {
				Component m = target.getComponent(ii);
				if (m.isVisible()) {
					Dimension d = m.getPreferredSize();
					dim.width = Math.max(dim.width, d.width);
					if (firstVisibleComponent) {
						firstVisibleComponent = false;
					}
					else {
						dim.height += _vgap;
					}
					dim.height += d.height;
				}
			}
			Insets insets = target.getInsets();
			dim.width += insets.left + insets.right + _hgap * 2;
			dim.height += insets.top + insets.bottom + _vgap * 2;
			return dim;
		}
	}


	/**
	 *  Description of the Method
	 *
	 *@param  target  Description of Parameter
	 */
	public void layoutContainer(Container target) {
		synchronized (target.getTreeLock()) {
			Insets insets = target.getInsets();
			int maxheight = target.getHeight() - (insets.top + insets.bottom) - 2 * _vgap;
		    int maxwidth = target.getWidth() - (insets.left + insets.right) - 2 * _hgap;
			int nmembers = target.getComponentCount();

			Dimension preferredSize = preferredLayoutSize(target);
			Dimension targetSize = target.getSize();

			int y = (_valign == TOP)    ? insets.top
				  : (_valign == CENTER) ? (targetSize.height - preferredSize.height) / 2
				  :						  targetSize.height - preferredSize.height - insets.bottom;

			for (int i = 0; i < nmembers; i++) {
				Component m = target.getComponent(i);
				if (m.isVisible()) {
					Dimension d = m.getPreferredSize();
			        if (_hfill ) {
			        	m.setSize(maxwidth, d.height);
			            d.width = maxwidth;
			        }
			        else {
			        	m.setSize(d.width, d.height);
			        }

					if ((y + d.height) <= maxheight) {
						if (y > 0) {
							y += _vgap;
						}

						int x = (_halign == LEFT) ? insets.left
							: (_halign == CENTER) ? (targetSize.width - d.width) / 2
							:						targetSize.width - d.width - insets.right;

						m.setLocation(x + _hgap, y + _vgap);

						y += d.getHeight();

					}
					else {
						break;
					}
				}
			}
		}
	}


	/**
	 *  Description of the Method
	 *
	 *@return    Description of the Returned Value
	 */
	public String toString() {
		String halign = "";
		switch (_halign) {
			case TOP:
				halign = "top";
				break;
			case CENTER:
				halign = "center";
				break;
			case BOTTOM:
				halign = "bottom";
				break;
		}
		String valign = "";
		switch (_valign) {
			case TOP:
				valign = "top";
				break;
			case CENTER:
				valign = "center";
				break;
			case BOTTOM:
				valign = "bottom";
				break;
		}
		return getClass().getName() + "[hgap=" + _hgap + ",vgap=" + _vgap + ",halign=" + halign + ",valign=" + valign + "]";
	}
}
