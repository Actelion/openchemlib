package com.actelion.research.gui;

import javax.swing.*;

public class LookAndFeelHelper {
	public static boolean isQuaQua() {
		return UIManager.getLookAndFeel().getName().contains("Quaqua");
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

	public static boolean isDarkLookAndFeel() {
		return UIManager.getLookAndFeel().getName().startsWith("Substance Graphite");
		}
	}
