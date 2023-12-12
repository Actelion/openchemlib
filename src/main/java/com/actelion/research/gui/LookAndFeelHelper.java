package com.actelion.research.gui;

import javax.swing.*;

public class LookAndFeelHelper {
	public static boolean isQuaQua() {
		return UIManager.getLookAndFeel().getName().contains("Quaqua");
		}

	public static boolean isAqua() {
		return UIManager.getLookAndFeel().getName().equals("Mac OS X");
	}

	/**
	 * This is Substance version 3 or 4
	 * @return
	 */
	public static boolean isOldSubstance() {
		return UIManager.getLookAndFeel().getName().equals("Substance");
		}

	/**
	 * This is Substance from version 5
	 * @return
	 */
	public static boolean isNewSubstance() {
		return UIManager.getLookAndFeel().getName().startsWith("Substance ");
		}

	public static boolean isSubstance() {
		return UIManager.getLookAndFeel().getName().startsWith("Substance");
	}

	public static boolean isRadiance() {
		return UIManager.getLookAndFeel().getName().startsWith("Radiance");
	}

	public static boolean isDarkLookAndFeel() {
		String lafName = UIManager.getLookAndFeel().getName();
		return lafName.startsWith("Substance Graphite")
			|| lafName.startsWith("Radiance Graphite")
			|| lafName.startsWith("Radiance Night");
		}
	}
