/*
 * @(#)FittingFunction.java
 *
 * contains the fitting function 'funcs' passed to NumericalRecipes.mrqmin()
 *
 * Copyright 2003 Actelion Ltd., Inc. All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 */

package com.actelion.research.util;

public class FittingFunction {
	public double chisq;
	public double alamda;
	public double fittingFunction(double x, double[] a, double[] dyda, int ma) {
		// handles fitting of a straight line
		// for a more useful equation overwrite this function
		dyda[0] = x;
		dyda[1] = 1;
		return a[0]*x+a[1];
		}
	}
