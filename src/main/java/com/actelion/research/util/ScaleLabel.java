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

package com.actelion.research.util;

/**
 * Contains shown label, the relative position on the scale and the numerical value
 */
public class ScaleLabel {
	public String label;
	public double position;
	public double value;

	/**
	 * Create an ScaleLabel object consisting of visible label,
	 * its relative position on the scale and the numerical value behind.
	 * @param label the label to be shown
	 * @param position relative label position, >= 0.0 and <= 1.0
	 * @param value numerical value, which in case of a logarithmic scale is log10 of the shown value
	 */
	public ScaleLabel(String label, double position, double value) {
		this.label = label;
		this.position = position;
		this.value = value;
		}
	}
