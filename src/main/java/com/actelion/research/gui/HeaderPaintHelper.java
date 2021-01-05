/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
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

package com.actelion.research.gui;

import javax.swing.*;
import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Paint;
import java.security.AccessControlException;

public class HeaderPaintHelper {
	private static Color sColor1,sColor2;

	public static Paint getHeaderPaint(boolean isSelected, int headerHeight) {
		// for the development we use a yellow paint
		boolean isDevelopment = false;
		try {
			isDevelopment = (System.getProperty("development") != null);
			}
		catch (AccessControlException ace) {}

		// the lighter color on the top
        Color color1 = LookAndFeelHelper.isDarkLookAndFeel() ?
		        (!isSelected ? new Color(0x606060) : isDevelopment ? new Color(0xC0C000) : sColor1 != null ? sColor1 : new Color(0x3838C0))
			  :	(!isSelected ? new Color(0xF0F0F0) : isDevelopment ? new Color(0xFFFFCD) : sColor1 != null ? sColor1 : new Color(0xAEDBFF));

		// the darker color on the bottom
		Color color2 = LookAndFeelHelper.isDarkLookAndFeel() ?
				(!isSelected ? new Color(0x404040) : isDevelopment ? new Color(0x404000) : sColor2 != null ? sColor2 : new Color(0x252560))
			  : (!isSelected ? new Color(0xD0D0D0) : isDevelopment ? new Color(0xAD9C00) : sColor2 != null ? sColor2 : new Color(0x0060FF));

		return new GradientPaint(0, -1, color1, 0, headerHeight, color2);
		}

    /**
     * Determines and returns the header's background color. Tries to lookup a special color from the L&F.
     * In case it is absent, it uses the standard internal frame background.
     *
     * @return the color of the header's background
     */
    private static Color getHeaderBackground(boolean selected) {
        return UIManager.getColor(selected ? "InternalFrame.activeTitleBackground" : "InternalFrame.inactiveTitleBackground");
        }

    public static void setSpotColors(Color c1, Color c2) {
    	sColor1 = c1;
    	sColor2 = c2;
		}

    /**
     * Determines and returns the header's text foreground color. Tries to lookup a special color from the
     * L&amp;F. In case it is absent, it uses the standard internal frame forground.
     *
     * @param selected true to lookup the active color, false for the inactive
     * @return the color of the foreground text
     */
    private static Color getTextForeground(boolean selected) {
        return UIManager.getColor(selected ? "InternalFrame.activeTitleForeground" : "InternalFrame.inactiveTitleForeground");
        }
	}
