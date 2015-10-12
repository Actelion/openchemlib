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

package com.actelion.research.util;

import java.awt.Color;

public class ColorHelper {
	/**
	 * Creates an intermediate color between color c1 and color c2
	 * in the RGB color space.
	 * @param c1
	 * @param c2
	 * @param ratio 0.0 -> returns c1; 1.0 -> returns c2
	 * @return color in between c1 and c2
	 */
	public static Color intermediateColor(Color c1, Color c2, float ratio) {
		return new Color((int)(c1.getRed()+ratio*(c2.getRed()-c1.getRed())),
						 (int)(c1.getGreen()+ratio*(c2.getGreen()-c1.getGreen())),
						 (int)(c1.getBlue()+ratio*(c2.getBlue()-c1.getBlue())));
		}

	/**
	 * This is a color's perceived brightness by the human eye
	 * @param c
	 * @return brightness from 0.0 to 1.0
	 */
	public static float perceivedBrightness(Color c) {
		return (c == null) ? 1.0f : ((float)(299 * c.getRed() + 587 * c.getGreen() + 114 * c.getBlue())) / 255000f;
		}
	}
