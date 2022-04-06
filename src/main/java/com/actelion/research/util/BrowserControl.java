package com.actelion.research.util;

import java.awt.*;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

public class BrowserControl {
	/**
	 * Display a file in the system browser.  If you want to display a
	 * file, you must include the absolute path name.
	 *
	 * @param url the file's url (the url must start with either "http://"
	 * or "file://").
	 */
	public static void displayURL(String url) {
		new Thread(() -> {
			String cmd = null;
			try {
				if (Platform.isWindows()) {
					// cmd = 'rundll32 url.dll,FileProtocolHandler http://...'
					cmd = WIN_PATH + " " + WIN_FLAG + " " + url;
					Runtime.getRuntime().exec(cmd);
				} else if (Platform.isMacintosh()) {
					cmd = OSX_PATH + " " + url;
					Runtime.getRuntime().exec(cmd);
					}
				else if (Desktop.isDesktopSupported()) {
					Desktop desktop = Desktop.getDesktop();
					if (desktop.isSupported(Desktop.Action.BROWSE)) {
						try {
							desktop.browse(new URI(url));
							}
						catch(URISyntaxException use) {
							use.printStackTrace();
							}
						}
					}
/*				else {
					// Under Unix, Netscape has to be running for the "-remote"
					// command to work.  So, we try sending the command and
					// check for an exit value.  If the exit command is 0,
					// it worked, otherwise we need to start the browser.
					// cmd = 'netscape -remote openURL(http://www.javaworld.com)'
					cmd = UNIX_PATH + " " + UNIX_PARAM_START + url + UNIX_PARAM_END;
					Process p = Runtime.getRuntime().exec(cmd);
					try {
						// wait for exit code -- if it's 0, command worked,
						// otherwise we need to start the browser up.
						int exitCode = p.waitFor();
						if (exitCode != 0) {
							// Command failed, start up the browser
							// cmd = 'firefox http://www.javaworld.com'
							cmd = UNIX_PATH + " " + url;
							p = Runtime.getRuntime().exec(cmd);
						}
					} catch (InterruptedException x) {
						System.err.println("Error bringing up browser, cmd='" + cmd
								+ "'");
						System.err.println("Caught: " + x);
					}
				}*/
			} catch (IOException x) {
				// couldn't exec browser
				System.err.println("Could not invoke browser, command=" + cmd);
				System.err.println("Caught: " + x);
			}
		} ).start();
	}

	/**
	 * Simple example.
	 */
	public static void main(String[] args) {
		displayURL("http://www.javaworld.com");
	}

	// The default system browser under windows.
	private static final String WIN_PATH = "rundll32";

	// The flag to display a url.
	private static final String WIN_FLAG = "url.dll,FileProtocolHandler";

	// The default browser under unix.
	private static final String UNIX_PATH = "firefox";

	// The flag to display a url.
//	private static final String UNIX_PARAM_START = "-remote \"openURL(";
//  private static final String UNIX_PARAM_END = ")\"";

	// The open command MacOSX.
	private static final String OSX_PATH = "open";
}