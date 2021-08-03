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
*	list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*	this list of conditions and the following disclaimer in the documentation
*	and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*	names of its contributors may be used to endorse or promote products
*	derived from this software without specific prior written permission.
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

package com.actelion.research.chem.io;

import javax.swing.filechooser.FileFilter;
import java.io.File;
import java.util.Iterator;
import java.util.TreeMap;

/**
 * A convenience implementation of FileFilter that filters out
 * all files except for those type extensions that it knows about.
 *
 * Extensions are of the type ".foo", which is typically found on
 * Windows and Unix boxes, but not on Macinthosh. Case is ignored.
 *
 * Example - create a new filter that filerts out all files
 * but gif and jpg image files:
 *
 *	 JFileChooser chooser = new JFileChooser();
 *	 CompoundFileFilter filter = new CompoundFileFilter(
 *				   new String{"gif", "jpg"}, "JPEG & GIF Images")
 *	 chooser.addChoosableFileFilter(filter);
 *	 chooser.showOpenDialog(this);
 */
public class CompoundFileFilter extends FileFilter {
	private TreeMap<String,CompoundFileFilter> mFilterMap;
	private String mDescription = null;
	private String mFullDescription = null;
	private boolean mUseExtensionsInDescription = true;

	/**
	 * Creates a file filter. If no filters are added, then all
	 * files are accepted.
	 */
	public CompoundFileFilter() {
		mFilterMap = new TreeMap();
		}

	/**
	 * Creates a file filter that accepts files with the given extension.
	 * Example: new CompoundFileFilter("jpg");
	 */
	public CompoundFileFilter(String extension) {
		this(extension,null);
		}

	/**
	 * Creates a file filter that accepts the given file type.
	 * Example: new CompoundFileFilter("jpg", "JPEG Image Images");
	 *
	 * Note that the "." before the extension is not needed. If
	 * provided, it will be ignored.
	 */
	public CompoundFileFilter(String extension, String description) {
		this();
		if (extension != null)
			addExtension(extension);
		if (description!=null)
			setDescription(description);
		}

	/**
	 * Creates a file filter from the given string array.
	 * Example: new CompoundFileFilter(String {"gif", "jpg"});
	 *
	 * Note that the "." before the extension is not needed and will be ignored.
	 */
	public CompoundFileFilter(String[] filters) {
		this(filters, null);
		}

	/**
	 * Creates a file filter from the given string array and description.
	 * Example: new CompoundFileFilter(String {"gif", "jpg"}, "Gif and JPG Images");
	 *
	 * Note that the "." before the extension is not needed and will be ignored.
	 */
	public CompoundFileFilter(String[] filters, String description) {
		this();
		for (int i=0; i < filters.length; i++)
			// add filters one by one
			addExtension(filters[i]);

		if(description!=null)
		    setDescription(description);
	}

	/**
	 * Return true if this file should be shown in the directory pane, false if it shouldn't.
	 * Files that begin with "." are ignored.
	 */
	@Override
	public boolean accept(File f) {
		if (f != null) {
			if (f.isDirectory())
				return true;

			String extension = getExtension(f);
			if (extension != null && mFilterMap.get(getExtension(f)) != null)
				return true;
			}

		return false;
		}

	/**
	 * Return the extension portion of the file's name .
	 *
	 * @see #getExtension
	 * @see FileFilter#accept
	 */
	public String getExtension(File f) {
		return CompoundFileHelper.getExtension(f);
		}

	/**
	 * Adds a filetype "dot" extension to filter against.
	 *
	 * For example: the following code will create a filter that filters
	 * out all files except those that end in ".jpg" and ".tif":
	 *
	 *   CompoundFileFilter filter = new CompoundFileFilter();
	 *   filter.addExtension("jpg");
	 *   filter.addExtension("tif");
	 *
	 * Note that the "." before the extension is not needed and will be ignored.
	 */
	public void addExtension(String extension) {
		if (mFilterMap == null)
			mFilterMap = new TreeMap<>();

		mFilterMap.put(extension.toLowerCase(), this);
		mFullDescription = null;
		}


	/**
	 * Returns the human readable description of this filter. For
	 * example: "JPEG and GIF Image Files (*.jpg, *.gif)"
	 */
	public String getDescription() {
		if (mFullDescription == null) {
			if (mDescription == null || isExtensionListInDescription()) {
				mFullDescription = (mDescription == null) ? "(" : mDescription + " (";

				// build the description from the extension list
				Iterator<String> iterator = mFilterMap.keySet().stream().sorted().iterator();
				mFullDescription += iterator.next();
				while (iterator.hasNext())
					mFullDescription += ", " + iterator.next();

				mFullDescription += ")";
				}
			else {
				mFullDescription = mDescription;
				}
			}
		return mFullDescription;
		}

	/**
	 * Sets the human readable description of this filter. For
	 * example: filter.setDescription("Gif and JPG Images");
	 */
	public void setDescription(String description) {
		mDescription = description;
		mFullDescription = null;
		}

	public void addDescription(String description) {
	if (mDescription == null)
		mDescription = description;
	else
		mDescription = mDescription.concat(", "+description);
		}

	/**
	 * Determines whether the extension list (.jpg, .gif, etc) should
	 * show up in the human readable description.
	 *
	 * Only relevent if a description was provided in the constructor
	 * or using setDescription();
	 */
	public void setExtensionListInDescription(boolean b) {
		mUseExtensionsInDescription = b;
		mFullDescription = null;
		}

	/**
	 * Returns whether the extension list (.jpg, .gif, etc) should
	 * show up in the human readable description.
	 *
	 * Only relevent if a description was provided in the constructor
	 * or using setDescription();
	 */
	public boolean isExtensionListInDescription() {
		return mUseExtensionsInDescription;
		}
	}
