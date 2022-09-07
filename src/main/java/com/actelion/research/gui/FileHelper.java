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

package com.actelion.research.gui;

import com.actelion.research.chem.io.CompoundFileHelper;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicBoolean;

public class FileHelper extends CompoundFileHelper {
	private static final long TIME_OUT = 5000;

	private Component mParent;

	public FileHelper(Component parent) {
		mParent = parent;
		}

	public String selectOption(String message, String title, String[] option) {
		if (SwingUtilities.isEventDispatchThread())
	        return (String)JOptionPane.showInputDialog(mParent, message, title,
            										   JOptionPane.QUESTION_MESSAGE, null, option, option[0]);
		else
			return null;
		}

	public void showMessage(String message) {
		JOptionPane.showMessageDialog(mParent, message);
		}

	/**
	 * For compatibility reasons...
	 * @param parent
	 * @param title
	 * @param filetypes
	 * @return
	 */
	public static File getFile(Component parent, String title, int filetypes) {
		return new FileHelper(parent).selectFileToOpen(title, filetypes);
		}

	public static ArrayList<File> getCompatibleFileList(File directory, int filetypes) {
		ArrayList<File> list = new ArrayList<>();
		if (fileExists(directory)) {
			javax.swing.filechooser.FileFilter ff = FileHelper.createFileFilter(filetypes, false);
			for (File file:directory.listFiles((File pathname) -> ff.accept(pathname) ))
				list.add(file);
			}

		return list;
		}

	public File selectFileToOpen(String title, int filetypes) {
		return selectFileToOpen(title, filetypes, null);
		}

	/**
	 * Shows a file-open-dialog, lets the user choose and returns the selected file.
	 * @param title of the dialog shown
	 * @param filetypes one or more file types defined in CompoundFileHelper
	 * @param initialFileName null or a suggested file name with or without complete path
	 * @return null or selected file
	 */
	public File selectFileToOpen(String title, int filetypes, String initialFileName) {
		JFileChooser fileChooser = new JFileChooser();

		// file chooser height does not automatically grow with UI scale factor
		if (HiDPIHelper.getUIScaleFactor() > 1)
			fileChooser.setPreferredSize(new Dimension(fileChooser.getPreferredSize().width,
					HiDPIHelper.scale(fileChooser.getPreferredSize().height)));

		fileChooser.setDialogTitle(title);
		fileChooser.setCurrentDirectory(getCurrentDirectory());
		if (filetypes == cFileTypeDirectory)
			fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		else if (filetypes != cFileTypeUnknown)
			fileChooser.setFileFilter(createFileFilter(filetypes, false));
		if (initialFileName != null) {
			int index = initialFileName.lastIndexOf(File.separatorChar);
			if (index == -1) {
				fileChooser.setSelectedFile(new File(FileHelper.getCurrentDirectory(), initialFileName));
				}
			else {
				String directory = initialFileName.substring(0, index);
				if (new File(directory).exists())
					fileChooser.setSelectedFile(new File(initialFileName));
				else
					fileChooser.setSelectedFile(new File(FileHelper.getCurrentDirectory(), initialFileName.substring(index+1)));
				}
			}
		int option = fileChooser.showOpenDialog(mParent);
		setCurrentDirectory(fileChooser.getCurrentDirectory());
		if (option != JFileChooser.APPROVE_OPTION)
			return null;
		File selectedFile = fileChooser.getSelectedFile();
		if (selectedFile.exists())
			return selectedFile;
		if (selectedFile.getName().contains(".") || filetypes == cFileTypeDirectory)
			return selectedFile;
		ArrayList<String> list = getExtensionList(filetypes);
		for (String extension:list) {
			File file = new File(selectedFile.getPath()+extension);
			if (file.exists())
				return file;
			}
		return selectedFile;
		}

	/**
	 * Shows a file-save-dialog, lets the user choose and return the file's path and name.
	 * @param title of the dialog shown
	 * @param filetype one of the file types defined in CompoundFileHelper
	 * @param newFileName null or a suggested file name with or without extension
	 * @return null or complete file path and name
	 */
	public String selectFileToSave(String title, int filetype, String newFileName) {
		JFileChooserOverwrite fileChooser = new JFileChooserOverwrite();

		// file chooser height does not automatically grow with UI scale factor
		if (HiDPIHelper.getUIScaleFactor() > 1)
			fileChooser.setPreferredSize(new Dimension(fileChooser.getPreferredSize().width,
					HiDPIHelper.scale(fileChooser.getPreferredSize().height)));

		fileChooser.setCurrentDirectory(getCurrentDirectory());
		fileChooser.setDialogTitle(title);
		fileChooser.setFileFilter(createFileFilter(filetype, true));
		fileChooser.setExtensions(FileHelper.getExtensions(filetype));
		if (newFileName == null) {
			fileChooser.setSelectedFile(new File(FileHelper.getCurrentDirectory(), "Untitled"));
			}
		else {
			int index = newFileName.lastIndexOf(File.separatorChar);
			if (index == -1) {
				fileChooser.setSelectedFile(new File(FileHelper.getCurrentDirectory(), newFileName));
				}
			else {
				String directory = newFileName.substring(0, index);
				if (new File(directory).exists())
					fileChooser.setSelectedFile(new File(newFileName));
				else
					fileChooser.setSelectedFile(new File(FileHelper.getCurrentDirectory(), newFileName.substring(index+1)));
				}
			}
		int option = fileChooser.showSaveDialog(mParent);
		setCurrentDirectory(fileChooser.getCurrentDirectory());
		return (option != JFileChooser.APPROVE_OPTION) ? null : fileChooser.getFile().getPath();
		}

	/**
	 * java.io.File.exists() and java.nio.file.Files.exists() may cause minutes of delay,
	 * if a file is/was on a network share which is currently unmounted. This version returns quickly.
	 */
	public static boolean fileExists(final File file) {
		return fileExists(file, TIME_OUT);
		}

		/**
		 * java.io.File.exists() and java.nio.file.Files.exists() may cause minutes of delay,
		 * if a file is/was on a network share which is currently unmounted. This version returns quickly.
		  */
	public static boolean fileExists(final File file, final long timeOutMillis) {
		final AtomicBoolean exists = new AtomicBoolean(false);
		final AtomicBoolean done = new AtomicBoolean(false);

		new Thread(() -> {
			exists.set(file.exists());
			done.set(true);
			} ).start();

		long start = System.currentTimeMillis();
		do {
			try { Thread.sleep(1); } catch (InterruptedException ie) {}

			} while (!done.get() && System.currentTimeMillis() < start + timeOutMillis);

		return exists.get();
		}
	}
