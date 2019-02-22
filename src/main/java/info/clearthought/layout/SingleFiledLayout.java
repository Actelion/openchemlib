/*
 * ====================================================================
 *
 * The Clearthought Software License, Version 2.0
 *
 * Copyright (c) 2001 Daniel Barbalace.  All rights reserved.
 *
 * Project maintained at https://tablelayout.dev.java.net
 *
 * I. Terms for redistribution of original source and binaries
 *
 * Redistribution and use of unmodified source and/or binaries are
 * permitted provided that the following condition is met:
 *
 * 1. Redistributions of original source code must retain the above
 *    copyright notice and license.  You are not required to redistribute
 *    the original source; you may choose to redistribute only the
 *    binaries.
 *
 * Basically, if you distribute unmodified source, you meet
 * automatically comply with the license with no additional effort on
 * your part.
 *
 * II. Terms for distribution of derived works via subclassing and/or
 *     composition.
 *
 * You may generate derived works by means of subclassing and/or
 * composition (e.g., the Adaptor Pattern), provided that the following
 * conditions are met:
 *
 * 1. Redistributions of original source code must retain the above
 *    copyright notice and license.  You are not required to redistribute
 *    the original source; you may choose to redistribute only the
 *    binaries.
 *
 * 2. The original software is not altered.
 *
 * 3. Derived works are not contained in the info.clearthought
 *    namespace/package or any subpackage of info.clearthought.
 *
 * 4. Derived works do not use the class or interface names from the
 *    info.clearthought... packages
 *
 * For example, you may define a class with the following full name:
 *    org.nameOfMyOrganization.layouts.RowMajorTableLayout
 *
 * However, you may not define a class with the either of the
 * following full names:
 *    info.clearthought.layout.RowMajorTableLayout
 *    org.nameOfMyOrganization.layouts.TableLayout
 *
 * III. Terms for redistribution of source modified via patch files.
 *
 * You may generate derived works by means of patch files provided
 * that the following conditions are met:
 *
 * 1. Redistributions of original source code must retain the above
 *    copyright notice and license.  You are not required to
 *    redistribute the original source; you may choose to redistribute
 *    only the binaries resulting from the patch files.
 *
 * 2. The original source files are not altered.  All alteration is
 *    done in patch files.
 *
 * 3. Derived works are not contained in the info.clearthought
 *    namespace/package or any subpackage of info.clearthought.  This
 *    means that your patch files must change the namespace/package
 *    for the derived work.  See section II for examples.
 *
 * 4. Derived works do not use the class or interface names from the
 *    info.clearthought... packages.  This means your patch files
 *    must change the names of the interfaces and classes they alter.
 *    See section II for examples.
 *
 * 5. Derived works must include the following disclaimer.
 *    "This work is derived from Clearthought's TableLayout,
 *     https://tablelayout.dev.java.net, by means of patch files
 *     rather than subclassing or composition.  Therefore, this work
 *     might not contain the latest fixes and features of TableLayout."
 *
 * IV. Terms for repackaging, transcoding, and compiling of binaries.
 *
 * You may do any of the following with the binaries of the
 * original software.
 *
 * 1. You may move binaries (.class files) from the original .jar file
 *    to your own .jar file.
 *
 * 2. You may move binaries from the original .jar file to other
 *    resource containing files, including but not limited to .zip,
 *    .gz, .tar, .dll, .exe files.
 *
 * 3. You may backend compile the binaries to any platform, including
 *    but not limited to Win32, Win64, MAC OS, Linux, Palm OS, any
 *    handheld or embedded platform.
 *
 * 4. You may transcribe the binaries to other virtual machine byte
 *    code protocols, including but not limited to .NET.
 *
 * V. License Disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR, AFFILATED BUSINESSES,
 * OR ANYONE ELSE BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 */

package info.clearthought.layout;



import java.awt.*;



/**
 * SingleFiledLayout lays out components singled filed.  This layout manager is
 * like FlowLayout except that all components are placed in a single row or
 * column.
 *
 * @version 1.1 April 4, 2002
 * @author  Daniel E. Barbalace
 */

public class SingleFiledLayout implements
    java.awt.LayoutManager,
    java.io.Serializable
{



/** Align components in a column */
public static final int COLUMN = 0;

/** Align components in a row */
public static final int ROW = 1;

/** Left justify components */
public static final int LEFT = 0;

/** Top justify components */
public static final int TOP = 0;

/** Center components */
public static final int CENTER = 1;

/** Full justify components */
public static final int FULL = 2;

/** Bottom justify components */
public static final int BOTTOM = 3;

/** Right justify components */
public static final int RIGHT = 4;

/** Default gap -- derived classes may override */
public static int DEFAULT_GAP = 5;



/** ROW or COLUMN -- should the components be aligned in a row or column */
protected int orientation;

/** LEFT, TOP, CENTER, FULL, BOTTOM, RIGHT -- how should components of different
    sizes be aligned */
protected int justification;

/** Space between components in pixels */
protected int gap;



/**
 * Constructs an instance of SingleFiledLayout that will align components in a
 * column using the default gap and LEFT justification.
 */

public SingleFiledLayout ()
{
    this (COLUMN, LEFT, DEFAULT_GAP);
}


/**
 * Constructs an instance of SingleFiledLayout using the default gap and LEFT
 * or TOP justification.
 *
 * @param orientation    ROW or COLUMN -- should the components be aligned in
 *                       a row or column
 */

public SingleFiledLayout (int orientation)
{
    this (orientation, LEFT, DEFAULT_GAP);
}



/**
 * Constructs an instance of SingleFiledLayout.
 *
 * @param orientation      ROW or COLUMN -- should the components be aligned in
 *                         a row or column
 * @param justification    LEFT, TOP, CENTER, FULL, BOTTOM, RIGHT -- how should
 *                         components of different sizes be aligned
 * @param gap              space between components in pixels
 */

public SingleFiledLayout (int orientation, int justification, int gap)
{
    // Validate parameters
    if (orientation != ROW)
        orientation = COLUMN;

    if ((justification != CENTER) && (justification != FULL) &&
        (justification != RIGHT))
        justification = LEFT;

    if (gap < 0)
        gap = 0;

    // Copy parameters
    this.orientation = orientation;
    this.justification = justification;
    this.gap = gap;
}



//******************************************************************************
//** java.awt.event.LayoutManager methods                                    ***
//******************************************************************************



/**
 * To lay out the specified container using this layout.  This method
 * repositions the components in the specified target container.
 *
 * <p>User code should not have to call this method directly.</p>
 *
 * @param container    container being served by this layout manager
 */

public void layoutContainer (Container container)
{
    // Use preferred size to get maximum width or height
    Dimension size = container.getSize();

    // Get the container's insets
    Insets inset = container.getInsets();

    // Start at top left of container
    int x = inset.left;
    int y = inset.top;

    // Get the components inside the container
    Component component[] = container.getComponents();

    // Components arranged in a column
    if (orientation == COLUMN)
        // Add each component
        for (int counter = 0; counter < component.length; counter++)
        {
            // Use preferred size of component
            Dimension d = component[counter].getPreferredSize();

            // Align component
            switch (justification)
            {
                case LEFT :
                    x = inset.left;
                break;

                case CENTER :
                    x = ((size.width - d.width) >> 1) +
                        inset.left - inset.right;
                break;

                case FULL :
                    x = inset.left;
                    d.width = size.width - inset.left - inset.right;
                break;

                case RIGHT :
                    x = size.width - d.width - inset.right;
                break;
            }

            // Set size and location
            component[counter].setBounds (x, y, d.width, d.height);

            // Increment y
            y += d.height + gap;
        }
    // Components arranged in a row
    else
        // Add each component
        for (int counter = 0; counter < component.length; counter++)
        {
            // Use preferred size of component
            Dimension d = component[counter].getPreferredSize();

            // Align component
            switch (justification)
            {
                case TOP :
                    y = inset.top;
                break;

                case CENTER :
                    y = ((size.height - d.height) >> 1) +
                        inset.top - inset.bottom;
                break;

                case FULL :
                    y = inset.top;
                    d.height = size.height - inset.top - inset.bottom;
                break;

                case BOTTOM :
                    y = size.height - d.height - inset.bottom;
                break;
            }

            // Set size and location
            component[counter].setBounds (x, y, d.width, d.height);

            // Increment x
            x += d.width + gap;
        }
}



/**
 * Determines the preferred size of the container argument using this layout.
 * The preferred size is the smallest size that, if used for the container's
 * size, will ensure that no component is truncated when the component is it's
 * preferred size.
 *
 * @param container    container being served by this layout manager
 *
 * @return a dimension indicating the container's preferred size
 */

public Dimension preferredLayoutSize (Container container)
{
    int totalWidth = 0;  // Width of all components
    int totalHeight = 0; // Height of all components

    // Get the components inside the container
    Component component[] = container.getComponents();

    // Components arranged in a column
    if (orientation == COLUMN)
    {
        // Add each component
        for (int counter = 0; counter < component.length; counter++)
        {
            Dimension d = component[counter].getPreferredSize();

            if (totalWidth < d.width)
                totalWidth = d.width;

            totalHeight += d.height + gap;
        }

        // Subtract extra gap
        totalHeight -= gap;
    }
    // Components arranged in a row
    else
    {
        // Add each component
        for (int counter = 0; counter < component.length; counter++)
        {
            Dimension d = component[counter].getPreferredSize();

            totalWidth += d.width + gap;

            if (totalHeight < d.height)
                totalHeight = d.height;
        }

        // Subtract extra gap
        totalWidth -= gap;
    }

    // Add insets to preferred size
    Insets inset = container.getInsets();
    totalWidth += inset.left + inset.right;
    totalHeight += inset.top + inset.bottom;

    return new Dimension(totalWidth, totalHeight);
}



/**
 * Determines the minimum size of the container argument using this layout.
 * The minimum size is the smallest size that, if used for the container's
 * size, will ensure that no component is truncated.  The minimum size is the
 * preferred size.
 *
 * @param container    container being served by this layout manager
 *
 * @return a dimension indicating the container's minimum size
 */

public Dimension minimumLayoutSize (Container container)
{
    int totalWidth = 0;  // Width of all components
    int totalHeight = 0; // Height of all components

    // Get the components inside the container
    Component component[] = container.getComponents();

    // Components arranged in a column
    if (orientation == COLUMN)
    {
        // Add each component
        for (int counter = 0; counter < component.length; counter++)
        {
            Dimension d = component[counter].getMinimumSize();

            if (totalWidth < d.width)
                totalWidth = d.width;

            totalHeight += d.height + gap;
        }

        // Subtract extra gap
        totalHeight -= gap;
    }
    // Components arranged in a row
    else
    {
        // Add each component
        for (int counter = 0; counter < component.length; counter++)
        {
            Dimension d = component[counter].getMinimumSize();

            totalWidth += d.width + gap;

            if (totalHeight < d.height)
                totalHeight = d.height;
        }

        // Subtract extra gap
        totalWidth =- gap;
    }

    // Add insets to preferred size
    Insets inset = container.getInsets();
    totalWidth += inset.left + inset.right;
    totalHeight += inset.top + inset.bottom;

    return new Dimension(totalWidth, totalHeight);
}



/**
 * Adds the specified component with the specified name to the layout.
 *
 * @param name         dummy parameter
 * @param component    component to add
 */

public void addLayoutComponent (String name, Component component)
{
}



/**
 * Removes the specified component with the specified name from the layout.
 *
 * @param component    component being removed
 */

public void removeLayoutComponent (Component component)
{
}



}
