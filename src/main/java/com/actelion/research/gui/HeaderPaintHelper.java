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

import javax.swing.*;
import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Paint;

public class HeaderPaintHelper {
	private static int[] sRGB = null;

	public static Paint getHeaderPaint(boolean isSelected, int headerHeight) {
		// for the development we use a yellow paint
		// the lighter color on the top
        int rgb1 = LookAndFeelHelper.isDarkLookAndFeel() ?
		        (!isSelected ? 0x606060 : sRGB != null ? sRGB[0] : 0x3838C0)
			  :	(!isSelected ? 0xF0F0F0 : sRGB != null ? sRGB[0] : 0xAEDBFF);

		// the darker color on the bottom
		int rgb2 = LookAndFeelHelper.isDarkLookAndFeel() ?
				(!isSelected ? 0x404040 : sRGB != null ? sRGB[1] : 0x252560)
			  : (!isSelected ? 0xD0D0D0 : sRGB != null ? sRGB[1] : 0x0060FF);

		return new GradientPaint(0, -1, new Color(rgb1), 0, headerHeight, new Color(rgb2));
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

	/**
	 * @return rgb values to paint header gradient; lighter color first
	 */
    public static int[] getThemeColors() {
    	return sRGB;
		}

	/**
	 * @param rgb values to paint header gradient; lighter color first
	 */
	public static void setThemeColors(int[] rgb) {
		sRGB = rgb;
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
