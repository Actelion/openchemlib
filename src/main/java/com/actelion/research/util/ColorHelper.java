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

package com.actelion.research.util;

import java.awt.Color;

public class ColorHelper {
	private static final float[] PERCEIVED_BRIGHTNESS = { 0.299f, 0.587f, 0.114f };
	// the desired perceived brightness difference between foreground color and background
	private static final float MIN_CONTRAST_TO_BACKGROUND = 0.30f;

	/**
	 * Creates an intermediate color between color c1 and color c2 in the RGB color space.
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
	 * Creates an intermediate color between rgb1 and rgb2 in the RGB color space.
	 * @param rgb1
	 * @param rgb2
	 * @param ratio 0.0 -> returns c1; 1.0 -> returns c2
	 * @return color in between c1 and c2
	 */
	public static int intermediateColor(int rgb1, int rgb2, float ratio) {
		int r1 = (rgb1 & 0x00FF0000) >> 16;
		int g1 = (rgb1 & 0x0000FF00) >> 8;
		int b1 = rgb1 & 0x000000FF;
		int r2 = (rgb2 & 0x00FF0000) >> 16;
		int g2 = (rgb2 & 0x0000FF00) >> 8;
		int b2 = rgb2 & 0x000000FF;
		return (r1+(Math.round(ratio*(r2-r1))) << 16)
			 + (g1+(Math.round(ratio*(g2-g1))) << 8)
			 +  b1+(Math.round(ratio*(b2-b1)));
		}

	/**
	 * Creates a new <code>Color</code> that is a brighter version of this
	 * <code>Color</code>. This is a copy of Color.brighter(), but lets you choose
	 * the factor.
	 * @param c the color to be brightened
	 * @param factor value < 1.0; Color.brighter() uses 0.7
	 * @return     a new <code>Color</code> object that is
	 *                 a brighter version of this <code>Color</code>
	 *                 with the same {@code alpha} value.
	 */
	public static Color brighter(Color c, float factor) {
		int r = c.getRed();
		int g = c.getGreen();
		int b = c.getBlue();
		int alpha = c.getAlpha();

        /* From 2D group:
         * 1. black.brighter() should return grey
         * 2. applying brighter to blue will always return blue, brighter
         * 3. non pure color (non zero rgb) will eventually return white
         */
		int i = (int)(1.0/(1.0-factor));
		if ( r == 0 && g == 0 && b == 0) {
			return new Color(i, i, i, alpha);
		}
		if ( r > 0 && r < i ) r = i;
		if ( g > 0 && g < i ) g = i;
		if ( b > 0 && b < i ) b = i;

		return new Color(Math.min((int)(r/factor), 255),
				Math.min((int)(g/factor), 255),
				Math.min((int)(b/factor), 255),
				alpha);
		}

	/**
	 * Creates a new <code>Color</code> that is a brighter version of this
	 * <code>Color</code>. This is a copy of Color.brighter(), but lets you choose
	 * the factor.
	 * @param argb the color to be brightened
	 * @param factor value < 1.0; Color.brighter() uses 0.7
	 * @return     a new <code>Color</code> object that is
	 *                 a brighter version of this <code>Color</code>
	 *                 with the same {@code alpha} value.
	 */
	public static int brighter(int argb, float factor) {
		int alpha = argb & 0xFF000000;
		int r = (argb & 0x00FF0000) >> 16;
		int g = (argb & 0x0000FF00) >> 8;
		int b = argb & 0x000000FF;

		/* From 2D group:
		 * 1. black.brighter() should return grey
		 * 2. applying brighter to blue will always return blue, brighter
		 * 3. non pure color (non zero rgb) will eventually return white
		 */
		int i = (int)(1.0/(1.0-factor));
		if ( r == 0 && g == 0 && b == 0)
			return alpha | (i << 16) | (i << 8) | i;

		if ( r > 0 && r < i ) r = i;
		if ( g > 0 && g < i ) g = i;
		if ( b > 0 && b < i ) b = i;

		return alpha
			 | (Math.min((int)(r/factor), 255) << 16)
			 | (Math.min((int)(g/factor), 255) << 8)
			 | Math.min((int)(b/factor), 255);
		}

	/**
	 * Creates a new <code>Color</code> that is a darker version of this
	 * <code>Color</code>. This is a copy of Color.darker(), but lets you choose
	 * the factor.
	 * @param c the color to be darkened
	 * @param factor value < 1.0; Color.darker() uses 0.7
	 * @return  a new <code>Color</code> object that is
	 *                    a darker version of this <code>Color</code>
	 *                    with the same {@code alpha} value.
	 */
	public static Color darker(Color c, float factor) {
		return new Color(Math.max((int)(c.getRed()  *factor), 0),
						 Math.max((int)(c.getGreen()*factor), 0),
						 Math.max((int)(c.getBlue() *factor), 0),
						 c.getAlpha());
		}

	/**
	 * Creates a new <code>Color</code> that is a darker version of this
	 * <code>Color</code>. This is a copy of Color.darker(), but lets you choose
	 * the factor.
	 * @param argb the color to be darkened
	 * @param factor value < 1.0; Color.darker() uses 0.7
	 * @return  a new <code>Color</code> object that is
	 *                    a darker version of this <code>Color</code>
	 *                    with the same {@code alpha} value.
	 */
	public static int darker(int argb, float factor) {
		return (argb & 0xFF000000)
			 | (Math.round(factor * ((argb & 0x00FF0000) >> 16)) << 16)
			 | (Math.round(factor * ((argb & 0x0000FF00) >> 8)) << 8)
			 | Math.round(factor * (argb & 0x000000FF));
		}

	/**
	 * This is a color's perceived brightness by the human eye
	 * @param c
	 * @return brightness from 0.0 to 1.0
	 */
	public static float perceivedBrightness(Color c) {
		return (c == null) ? 1.0f : (PERCEIVED_BRIGHTNESS[0]*c.getRed()
									+PERCEIVED_BRIGHTNESS[1]*c.getGreen()
									+PERCEIVED_BRIGHTNESS[2]*c.getBlue()) / 255f;
		}

	/**
	 * This is a color's perceived brightness by the human eye
	 * @param cc color components r,g,b,a
	 * @return brightness from 0.0 to 1.0
	 */
	public static float perceivedBrightness(float[] cc) {
		return PERCEIVED_BRIGHTNESS[0]*cc[0]
			 + PERCEIVED_BRIGHTNESS[1]*cc[1]
			 + PERCEIVED_BRIGHTNESS[2]*cc[2];
	}

	/**
	 * This is a color's perceived brightness by the human eye
	 * @param argb
	 * @return brightness from 0.0 to 1.0
	 */
	public static float perceivedBrightness(int argb) {
		return (PERCEIVED_BRIGHTNESS[0] * ((argb & 0x00FF0000) >> 16)
				+PERCEIVED_BRIGHTNESS[1]*((argb & 0x0000FF00) >> 8)
				+PERCEIVED_BRIGHTNESS[2]*(argb & 0x000000FF)) / 255f;
		}

	/**
	 * Creates a new color with the hue taken from color <code>c</code>, but adjusted
	 * in brightness to match the desired perceived brightness.
	 * @param c
	 * @param perceivedBrightness
	 * @return
	 */
	public static Color createColor(Color c, float perceivedBrightness) {
		float[] cc = c.getRGBComponents(null);
		createColor(cc, perceivedBrightness);
		return new Color(cc[0], cc[1], cc[2], cc[3]);
		}

	/**
	 * Creates a new color with the hue taken from color <code>c</code>, but adjusted
	 * in brightness to match the desired perceived brightness.
	 * @param argb
	 * @param perceivedBrightness
	 * @return
	 */
	public static int createColor(int argb, float perceivedBrightness) {
		float[] cc = new float[4];
		float f = 1f / 255f;
		cc[0] = f * ((argb & 0x00FF0000) >> 16);
		cc[1] = f * ((argb & 0x0000FF00) >> 8);
		cc[2] = f *  (argb & 0x000000FF);
		createColor(cc, perceivedBrightness);
		return (argb & 0xFF000000)
			 | (Math.round(cc[0] * 255) << 16)
			 | (Math.round(cc[1] * 255) << 8)
			 |  Math.round(cc[2] * 255);
		}

	/**
	 * Creates a new color with the hue taken from color <code>cc</code>, but adjusted
	 * in brightness to match the desired perceived brightness.
	 * @param cc color components r,g,b,a
	 * @param perceivedBrightness
	 */
	private static void createColor(float[] cc, float perceivedBrightness) {
		float pb = perceivedBrightness(cc);
		if (pb == 0f) {
			cc[0] = 0f;
			cc[1] = 0f;
			cc[2] = 0f;
			return;
			}

		float f = perceivedBrightness / pb;
		float surplusBrightness = 0f;
		float sum = 0f;
		for (int i=0; i<3; i++) {
			cc[i] *= f;
			if (cc[i] < 1f) {
				sum += PERCEIVED_BRIGHTNESS[i];
				}
			else {
				surplusBrightness += (cc[i] - 1f) * PERCEIVED_BRIGHTNESS[i];
				cc[i] = 1f;
				}
			}
		if (surplusBrightness != 0) {
			float remainingBrightness = 0f;
			for (int i=0; i<3; i++) {
				if (cc[i] < 1f) {
					cc[i] += surplusBrightness / sum;
					if (cc[i] > 1f) {
						remainingBrightness += (cc[i] - 1f) * PERCEIVED_BRIGHTNESS[i];
						cc[i] = 1f;
						}
					}
				}
			if (remainingBrightness != 0f) {
				for (int i=0; i<3; i++) {
					if (cc[i] < 1f) {
						cc[i] += remainingBrightness / PERCEIVED_BRIGHTNESS[i];
						if (cc[i] > 1f) {
							cc[i] = 1f;
							}
						}
					}
				}
			}
		}

	/**
	 * Based on the differences of hue and perceived brightness of foreground color
	 * <code>fg</code> and background color <code>bg</code>, this method checks and
	 * possibly adjusts the given foreground color <code>fg</code> such that its hue
	 * stays unchanged, but its brightness is adapted to make it better perceivable
	 * on the background.
	 * @param fg foreground color
	 * @param bg background color
//	 * @param minContrast the minimum desired contrast (0 to 0.5)
	 * @return unchanged or adjusted fg
	 */
	public static Color getContrastColor(Color fg, Color bg) {
		float bgb = ColorHelper.perceivedBrightness(bg);
		float fgb = ColorHelper.perceivedBrightness(fg);

		float contrast = Math.abs(bgb - fgb);
		if (contrast > MIN_CONTRAST_TO_BACKGROUND)
			return fg;

		float[] hsbBG = Color.RGBtoHSB(bg.getRed(), bg.getGreen(), bg.getBlue(), null);
		float[] hsbFG = Color.RGBtoHSB(fg.getRed(), fg.getGreen(), fg.getBlue(), null);

		double hueDif = Math.abs(hsbFG[0] - hsbBG[0]);
		if (hueDif > 0.5)
			hueDif = 1.0 - hueDif;

		float saturationFactor = 1-Math.max(hsbFG[1], hsbBG[1]);
		float brightnessFactor = Math.abs(fgb + bgb - 1);
		float hueDifferenceFactor = (float)Math.cos(Math.PI*hueDif*3);

		float neededContrast = MIN_CONTRAST_TO_BACKGROUND * Math.max(saturationFactor, Math.max(brightnessFactor, hueDifferenceFactor));

		if (contrast > neededContrast)
			return fg;

		boolean darken = (fgb > bgb) ? (fgb + neededContrast > 1f) : (fgb - neededContrast > 0f);

		return createColor(fg, darken ? bgb - neededContrast : bgb + neededContrast);
		}

	/**
	 * Based on the differences of hue and perceived brightness of foreground color
	 * <code>fg</code> and background color <code>bg</code>, this method checks and
	 * possibly adjusts the given foreground color <code>fg</code> such that its hue
	 * stays unchanged, but its brightness is adapted to make it better perceivable
	 * on the background.
	 * @param fgRGB foreground color
	 * @param bgRGB background color
	//	 * @param minContrast the minimum desired contrast (0 to 0.5)
	 * @return unchanged or adjusted fg
	 */
	public static int getContrastColor(int fgRGB, int bgRGB) {
		float bgb = ColorHelper.perceivedBrightness(bgRGB);
		float fgb = ColorHelper.perceivedBrightness(fgRGB);

		float contrast = Math.abs(bgb - fgb);
		if (contrast > MIN_CONTRAST_TO_BACKGROUND)
			return fgRGB;

		float[] hsbBG = Color.RGBtoHSB((bgRGB & 0x00FF0000) >> 16, (bgRGB & 0x0000FF00) >> 8, bgRGB & 0x000000FF, null);
		float[] hsbFG = Color.RGBtoHSB((fgRGB & 0x00FF0000) >> 16, (fgRGB & 0x0000FF00) >> 8, fgRGB & 0x000000FF, null);

		double hueDif = Math.abs(hsbFG[0] - hsbBG[0]);
		if (hueDif > 0.5)
			hueDif = 1.0 - hueDif;

		float saturationFactor = 1-Math.max(hsbFG[1], hsbBG[1]);
		float brightnessFactor = Math.abs(fgb + bgb - 1);
		float hueDifferenceFactor = (float)Math.cos(Math.PI*hueDif*3);

		float neededContrast = MIN_CONTRAST_TO_BACKGROUND * Math.max(saturationFactor, Math.max(brightnessFactor, hueDifferenceFactor));

		if (contrast > neededContrast)
			return fgRGB;

		boolean darken = (fgb > bgb) ? (fgb + neededContrast > 1f) : (fgb - neededContrast > 0f);

		return createColor(fgRGB, darken ? bgb - neededContrast : bgb + neededContrast);
		}
	}
