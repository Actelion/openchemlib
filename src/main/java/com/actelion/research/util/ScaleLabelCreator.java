/*
 * Copyright (c) 1997 - 2022
 * Idorsia Pharmaceuticals Ltd.
 * Hegenheimermattweg 91
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
 * 3. Neither the name of the copyright holder nor the
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

package com.actelion.research.util;

import java.util.ArrayList;

public class ScaleLabelCreator {

	/**
	 * Creates a list of scale labels with their relative positions for a linear
	 * numerical scale with given numerical value range.
	 * @param rangeLow low limit of numerical range
	 * @param rangeHigh high limit of numerical range
	 * @return list of labels with relative positions
	 */
	public static ArrayList<ScaleLabel> createLinearLabelList(double rangeLow, double rangeHigh) {
		if (rangeHigh <= rangeLow)
			return null;

		double range = rangeHigh - rangeLow;

		int exponent = 0;
		while (range >= 50.0) {
			rangeLow /= 10;
			range /= 10;
			exponent++;
			}
		while (range < 5.0) {
			rangeLow *= 10;
			range *= 10.0;
			exponent--;
			}

		int gridSpacing = (int)(range / 10);
	    if (gridSpacing < 1)
			gridSpacing = 1;
		else if (gridSpacing < 2)
			gridSpacing = 2;
		else
			gridSpacing = 5;

	    int value = (rangeLow < 0) ?
		      (int)(rangeLow - 0.0000001 - (rangeLow % gridSpacing))
			: (int)(rangeLow + 0.0000001 + gridSpacing - (rangeLow % gridSpacing));

		ArrayList<ScaleLabel> labelList = new ArrayList<ScaleLabel>();
		while (value < (rangeLow + range)) {
			double position = (value-rangeLow) / range;

			labelList.add(new ScaleLabel(DoubleFormat.toShortString(value, exponent), position, value*Math.pow(10, exponent)));

			value += gridSpacing;
			}

		return labelList;
		}

	/**
	 * Creates a list of scale labels with their relative positions for a
	 * logarithmical numerical scale with given numerical value range.
	 * @param rangeLow log10 of low limit of numerical range
	 * @param rangeHigh log10 of high limit of numerical range
	 * @return list of labels with relative positions
	 */
	public static ArrayList<ScaleLabel> createLogarithmicLabelList(double rangeLow, double rangeHigh) {
		if (rangeHigh <= rangeLow)
			return null;

		double range = rangeHigh - rangeLow;

        int intMin = (int)Math.floor(rangeLow);
        int intMax = (int)Math.floor(rangeHigh);
        
		ArrayList<ScaleLabel> labelList = new ArrayList<ScaleLabel>();
        if (range > 5.4) {
            int step = 1 + (int)range/10;
            for (int i=intMin; i<=intMax; i+=step)
                addLogarithmicScaleLabel(labelList, i, rangeLow, range);
            }
        else if (range > 3.6) {
            for (int i=intMin; i<=intMax; i++) {
                addLogarithmicScaleLabel(labelList, i, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.47712125472f, rangeLow, range);
                }
            }
        else if (range > 1.8) {
            for (int i=intMin; i<=intMax; i++) {
                addLogarithmicScaleLabel(labelList, i, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.301029996f, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.698970004f, rangeLow, range);
                }
            }
        else if (range > 1.0) {
            for (int i=intMin; i<=intMax; i++) {
                addLogarithmicScaleLabel(labelList, i, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.176091259f, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.301029996f, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.477121255f, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.698970004f, rangeLow, range);
                addLogarithmicScaleLabel(labelList, i + 0.84509804f, rangeLow, range);
                }
            }
		else if (range > 0.6) {
			for (int i=intMin; i<=intMax; i++) {
				addLogarithmicScaleLabel(labelList, i, rangeLow, range);						// 1.0
				addLogarithmicScaleLabel(labelList, i + 0.113943352f, rangeLow, range);	// 1.3
				addLogarithmicScaleLabel(labelList, i + 0.204119983f, rangeLow, range);	// 1.6
				addLogarithmicScaleLabel(labelList, i + 0.301029996f, rangeLow, range);	// 2.0
				addLogarithmicScaleLabel(labelList, i + 0.397940009f, rangeLow, range);	// 2.5
				addLogarithmicScaleLabel(labelList, i + 0.505149978f, rangeLow, range);	// 3.2
				addLogarithmicScaleLabel(labelList, i + 0.602059991f, rangeLow, range);	// 4.0
				addLogarithmicScaleLabel(labelList, i + 0.698970004f, rangeLow, range);	// 5.0
				addLogarithmicScaleLabel(labelList, i + 0.806179974f, rangeLow, range);	// 6.4
				addLogarithmicScaleLabel(labelList, i + 0.903089987f, rangeLow, range);	// 8.0
				}
			}
        else {
			double start = Math.pow(10, rangeLow);
			double length = Math.pow(10, rangeLow+range) - start;

            int exponent = 0;
            while (length >= 50.0) {
                start /= 10;
                length /= 10;
                exponent++;
                }
            while (length < 5.0) {
                start *= 10;
                length *= 10.0;
                exponent--;
                }

            int gridSpacing = (int)(length / 10);
            if (gridSpacing < 1)
                gridSpacing = 1;
            else if (gridSpacing < 2)
                gridSpacing = 2;
            else
                gridSpacing = 5;

            int value = (start < 0) ?
                  (int)(start - 0.0000001 - (start % gridSpacing))
                : (int)(start + 0.0000001 + gridSpacing - (start % gridSpacing));
            while (value < (start + length)) {
				double log = Math.log10(value) + exponent;
				double position = (log-rangeLow) / range;
                labelList.add(new ScaleLabel(DoubleFormat.toShortString(value, exponent), position, log));
                value += gridSpacing;
                }
            }

        return labelList;
		}

    private static void addLogarithmicScaleLabel(ArrayList<ScaleLabel> labelList, double value, double rangeLow, double range) {
        if (value >= rangeLow && value <= rangeLow+range) {
			double position = (value-rangeLow) / range;
            labelList.add(new ScaleLabel(DoubleFormat.toString(Math.pow(10, value), 3, true), position, value));
            }
        }
	}
