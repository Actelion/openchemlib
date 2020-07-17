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


import java.awt.*;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * Requires JRE 1.6
 */
// Think about renaming this class
public class Platform
{
    private static boolean isWin = (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0);
    private static boolean isMac = (System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0);
    private static boolean isLinux = (System.getProperty("os.name").toLowerCase().indexOf("nux") >= 0);
    private static boolean isUnix = (System.getProperty("os.name").toLowerCase().indexOf("nix") >= 0);

    private static Properties PATH_PROPERTIES = null;

    private static String WINPATH = "\\\\idorsia.com\\app\\DDSC\\config";
    private static String WINPROPERTIES = WINPATH + "\\ApplicationDatabase.db";

    private static final String[] WINDOWS_APP_DIR = { "C:\\Program Files\\", "C:\\Program Files (x86)\\" };

	private static final String[][] WINDOWS_APPL_NAME = { { "datawarrior" , "DataWarrior\\DataWarrior.exe" },
														  { "orbit" , "Orbit\\Orbit.exe" },
														  { "pymol" , "PyMOL\\pymol.exe" } };

    private static final String[][] MACINTOSH_APP_NAME = { { "orbit" , "Orbit Image Analysis" },
    													   { "datawarrior", "DataWarrior" },
    													   { "actelion3d", "Actelion3D" },
                                                           { "pymol", "PyMol" },
														   { "spirit", "Spirit Database" } };

    public static boolean isWindows()
    {
        return isWin;
    }

    public static boolean isMacintosh()
    {
        return isMac;
    }

    public static boolean isLinux()
    {
        return isLinux;
    }

    public static boolean isUnix()
    {
        return isUnix;
    }

	public static boolean is64BitJRE() { return System.getProperty("os.arch").indexOf("64")!=-1; }

    /**
     * Start an executable with parameters. In non-windows environment this is searching the PATH
     * In the windows environment @Actelion this first checks the standard locations where applications
     * are installed usually \\actelch02\pgm. This is achieved by consulting the ApplicationDatabase.db file
     * located under \\actelch02\pgm\ActelionResearch. The file format is a standard Java properties file
     * where the key contains the name of the application and the value the respective absolute path of the executable
     * e.g.
     * datawarrior=//actelch02/pgm/Datawarrior/DataWarrior.exe
     * @param program
     * @param args
     * @throws IOException
     */
    public static Process execute(String program, String... args) throws IOException
    {
        String executable = findExecutable(program);
        List<String> arguments = new ArrayList<String>();
        arguments.add(executable);
        if (args != null && args.length > 0) {
            for (String a : args)
                arguments.add(a);
        }
        final Process command = systemExec(arguments.toArray(new String[0]));

        return command;
    }

	/**
	 * Start an executable with parameters. In non-windows environment this is searching the PATH
	 * In the windows environment @Actelion this first checks the standard locations where applications
	 * are installed usually \\actelch02\pgm. This is achieved by consulting the ApplicationDatabase.db file
	 * located under \\actelch02\pgm\ActelionResearch. The file format is a standard Java properties file
	 * where the key contains the name of the application and the value the respective absolute path of the executable
	 * e.g.
	 * datawarrior=//actelch02/pgm/Datawarrior/DataWarrior.exe
	 * @param programAndArgs
	 * @throws IOException
	 */
	public static void execute(String[] programAndArgs) throws IOException
	{
		programAndArgs[0] = findExecutable(programAndArgs[0]);
        systemExec(programAndArgs);
	}

    /**
     * Wrapper for Runtime.exec() that allows different handling by platform
     * This method fixes forking issues in Windows when starting Mercury (by default new processes do not seem to be forked correctly)
     * @param programAndArgs
     * @return
     * @throws IOException
     */
	public static Process systemExec(String[] programAndArgs) throws IOException {
	    if(isWindows()) {
	        // prepend command args with: cmd /c start ""
	        String[] modifiedProgramAndArgs = new String[programAndArgs.length+4];
	        modifiedProgramAndArgs[0] = "cmd";
	        modifiedProgramAndArgs[1] = "/c";
	        modifiedProgramAndArgs[2] = "start";
            // window title
	        modifiedProgramAndArgs[3] = "\"\"";
	        System.arraycopy(programAndArgs, 0, modifiedProgramAndArgs, 4, programAndArgs.length);
	        programAndArgs = modifiedProgramAndArgs;
        }
        return Runtime.getRuntime().exec(programAndArgs);
    }

	/**
     * Given a filename, open this file with the default application
     * @param doc File name to open
     * @throws IOException if the file is unavailable
     * @throws UnsupportedOperationException, if Desktop integration of the platform is unavailable
     */
    public static void openDocument(String doc) throws IOException, UnsupportedOperationException
    {
        if (!Desktop.isDesktopSupported()) {
            if (isLinux())
                throw new UnsupportedOperationException("Please check your OS installation: Is the libgnome2 library missing?");
            else
                throw new UnsupportedOperationException("Desktop integration on this OS/Machine not supported");
        } else {
            Desktop.getDesktop().open(new File(doc));
        }
    }

    /**
     * In a platform specific way tries to find the path to the application.
     * @param name
     * @return valid full path to application or null
     */
    private static String findExecutable(String name)
    {
        String res = name;
        if (isWindows()) {
            for (String[] appKeyAndName: WINDOWS_APPL_NAME) {
	            if (appKeyAndName[0].equals(name)) {
		            for (String appDir: WINDOWS_APP_DIR) {
			            String path = appDir.concat(appKeyAndName[1]);
			            if (new File(path).exists())
				            return path;
		            }
	            }
            }

            // If we don't have it predefined, just try C:\Program Files\name\name.exe
	        for (String appDir: WINDOWS_APP_DIR) {
		        String path = appDir.concat(name).concat("\\").concat(name).concat(".exe");
		        if (new File(path).exists())
			        return path;
	        }

	        try {
                //If no local exe is found, load a properties file
                if (PATH_PROPERTIES == null) {
                    PATH_PROPERTIES = new Properties();
                    PATH_PROPERTIES.load(new FileInputStream(WINPROPERTIES));
                }
                String value = PATH_PROPERTIES.getProperty(name);
                if (value != null)
                    res = value;
            } catch (Throwable e) {
                System.err.println("Error reading Application Database file: " + e);
            }
        }
        if (isMacintosh()) {
        	for (String[] appKeyAndName: MACINTOSH_APP_NAME) {
        		if (appKeyAndName[0].equals(name)) {
        		    // we assume that the name of the launcher is equal to parameter name
                    String path = "/Applications/"+appKeyAndName[1]+".app/Contents/MacOS/"+name;
                    if (new File(path).exists())
                        return path;

                    // if the JRE7+ way doesn's work, check if we still have a JRE6 based app
        			path = "/Applications/"+appKeyAndName[1]+".app/Contents/MacOS/JavaApplicationStub";
        			return new File(path).exists() ? path : res;
        		}
        	}
        }
        if (isLinux()) {
        	String path = "/opt/actelion/"+name+"/"+name;
        	if (new File(path).exists())
        		return path;
			path = "/opt/idorsia/"+name+"/"+name;
			if (new File(path).exists())
				return path;
        	path = "/opt/"+name+"/"+name;
        	if (new File(path).exists())
        		return path;
        	return res;
        }
        return res;
    }


/*
    public static void main(String args[]) throws Exception
    {
        if (isWindows()) {
            openDocument("C:\\temp\\supplier.txt");
            openDocument("C:\\temp\\out.elb");
            openDocument("C:\\temp\\foo.sdf");
            execute("Notepad");
            execute("cmd", "/c", "start");
            execute("nemo");
            execute("datawarrior","c:\\temp\\foo.sdf");
        }
        if (isWindows()) {
        	
        }
    }
     */
}
