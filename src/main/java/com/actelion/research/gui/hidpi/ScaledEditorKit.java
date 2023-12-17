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


		import java.awt.Color;
		import java.awt.Font;
		import java.io.StringReader;

		import javax.swing.text.Document;
		import javax.swing.text.Element;
		import javax.swing.text.View;
		import javax.swing.text.ViewFactory;
		import javax.swing.text.html.HTMLDocument;
		import javax.swing.text.html.HTMLEditorKit;
		import javax.swing.text.html.ImageView;
		import javax.swing.text.html.StyleSheet;

@SuppressWarnings("serial")
public class ScaledEditorKit extends HTMLEditorKit {
	public static float getFontScaleFactor() {
		return HiDPIHelper.getUIScaleFactor();
	}

	/** Shared base style for all documents created by us use. */
	private static StyleSheet defaultStyles;

	public ScaledEditorKit() {
	};

	/**
	 * Overriden to return our own slimmed down style sheet.
	 */
	public StyleSheet getStyleSheet() {
		if (defaultStyles == null) {
			defaultStyles = new StyleSheet();
			StringReader r = new StringReader(ScaledHTML.styleChanges);
			try {
				defaultStyles.loadRules(r, null);
			}
			catch (Throwable e) {
				// don't want to die in static initialization...
				// just display things wrong.
			}
			r.close();
			defaultStyles.addStyleSheet(super.getStyleSheet());
		}
		return defaultStyles;
	}

	@Override
	public Document createDefaultDocument() {
		StyleSheet styles = getStyleSheet();
		StyleSheet ss = new ScaledStyleSheet();
		ss.addStyleSheet(styles);
		HTMLDocument doc = new HTMLDocument(ss);
		doc.setParser(getParser());
		doc.setAsynchronousLoadPriority(4);
		doc.setTokenThreshold(100);
		return doc;
	}

	/**
	 * Sets the async policy to flush everything in one chunk, and
	 * to not display unknown tags.
	 */
	Document createDefaultDocument(Font defaultFont, Color foreground) {
		StyleSheet styles = getStyleSheet();
		StyleSheet ss = new ScaledStyleSheet();
		ss.addStyleSheet(styles);
		HTMLDocument doc = new HTMLDocument(ss);
		doc.setPreservesUnknownTags(false);
		doc.getStyleSheet().addRule(displayPropertiesToCSS(defaultFont, foreground));
		doc.setParser(getParser());
		doc.setAsynchronousLoadPriority(Integer.MAX_VALUE);
		doc.setPreservesUnknownTags(false);
		return doc;
	}

	private String displayPropertiesToCSS(Font font, Color fg) {
		StringBuffer rule = new StringBuffer("body {");
		if (font != null) {
			rule.append(" font-family: ");
			rule.append(font.getFamily());
			rule.append(" ; ");
			rule.append(" font-size: ");
			final int fontSize = Math.round(font.getSize() / getFontScaleFactor());
			rule.append(fontSize);
			rule.append("pt ;");
			if (font.isBold()) {
				rule.append(" font-weight: bold ; ");
			}
			if (font.isItalic()) {
				rule.append(" font-style: italic ; ");
			}
		}
		if (fg != null) {
			rule.append(" color: ").append(colorToString(fg)).append(" ; ");
		}
		rule.append(" }");
		return rule.toString();
	}

	private String colorToString(final Color col) {
		if (col == null) {
			return null;
		}
		return String.format("#%02x%02x%02x", col.getRed(), col.getGreen(), col.getBlue());
	}

	/**
	 * Returns the ViewFactory that is used to make sure the Views don't
	 * load in the background.
	 */
	public ViewFactory getViewFactory() {
		if (basicHTMLViewFactory == null) {
			basicHTMLViewFactory = new BasicHTMLViewFactory();
		}
		return basicHTMLViewFactory;
	}

	/**
	 * BasicHTMLViewFactory extends HTMLFactory to force images to be loaded
	 * synchronously.
	 */
	static class BasicHTMLViewFactory extends HTMLEditorKit.HTMLFactory {
		public View create(Element elem) {
			View view = super.create(elem);

			if (view instanceof ImageView) {
				((ImageView)view).setLoadsSynchronously(true);
			}
			return view;
		}
	}


	static public ScaledEditorKit create() {
		if (kit == null) {
			kit = new ScaledEditorKit();
		}
		return kit;
	}

	/**
	 * The source of the html renderers
	 */
	private static ScaledEditorKit kit;
	/**
	 * Creates the Views that visually represent the model.
	 */
	static ViewFactory basicHTMLViewFactory;
}
