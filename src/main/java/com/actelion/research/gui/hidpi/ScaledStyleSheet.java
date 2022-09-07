package com.actelion.research.gui.hidpi;

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
		import javax.swing.text.AttributeSet;
import javax.swing.text.html.CSS;
import javax.swing.text.html.StyleSheet;
import java.awt.*;


public class ScaledStyleSheet extends StyleSheet{

	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	public Font getFont(AttributeSet a) {
		final Font font = super.getFont(a);
		final float fontScaleFactor = getFontScaleFactor(a);
		return super.getFont(font.getFamily(), font.getStyle(), Math.round(font.getSize2D() * fontScaleFactor));
	}

	private float getFontScaleFactor(AttributeSet a) {
		final Object attribute = a.getAttribute(CSS.Attribute.FONT_SIZE);
		if(attribute == null)
			return ScaledEditorKit.getFontScaleFactor();
		final String fontSize = attribute.toString();
		final int fsLength = fontSize.length();
		if(fsLength <= 1
				|| Character.isDigit(fontSize.charAt(fsLength-1))
				|| fontSize.endsWith("pt"))
			return ScaledEditorKit.getFontScaleFactor();
		if(fontSize.endsWith("px"))
			return 1/1.3f;
//		if(fontSize.endsWith("%") || fontSize.endsWith("em") || fontSize.endsWith("ex")
//				|| fontSize.endsWith("er"))
//			return getFontScaleFactor(a);   outcommented, because of stack overflow; TLS 2022Aug12
		return ScaledEditorKit.getFontScaleFactor();
	}
}
