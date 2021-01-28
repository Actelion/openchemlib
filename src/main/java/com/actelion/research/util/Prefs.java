package com.actelion.research.util;

import java.util.prefs.Preferences;

public class Prefs {
	private static final String PREFERENCES_ROOT = "org.openmolecules.openchemlib";

	public static Preferences get() {
		return java.util.prefs.Preferences.userRoot().node(PREFERENCES_ROOT);
		}

	public static String getString(String key, String defaultValue) {
		return get().get(key, defaultValue);
		}

	public static void setString(String key, String value) {
		get().put(key, value);
		}
	}
