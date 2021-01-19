/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.gui;

import javax.swing.*;
import java.awt.*;

public class JComboBoxWithColor extends JComboBox {
    static final long serialVersionUID = 0x20070316;

	public JComboBoxWithColor() {
        setRenderer(new DefaultListCellRenderer() {
            static final long serialVersionUID = 0x20070316;
            private Color mForeground = UIManager.getColor("ComboBox.foreground");

            @Override
            public void paintComponent(Graphics g) {
                setForeground(mForeground);
                super.paintComponent(g);
                }
            
            public Component getListCellRendererComponent(JList list, Object value,
                                                          int index, boolean isSelected,
                                                          boolean hasFocus) {
                if (value instanceof ComboBoxColorItem) {
                    Component c = super.getListCellRendererComponent(
                            list, ((ComboBoxColorItem)value).getText(),
                            index, isSelected, hasFocus);
                    c.setForeground(((ComboBoxColorItem)value).getColor());
                    mForeground = ((ComboBoxColorItem)value).getColor();
                    return c;
                    }
                return super.getListCellRendererComponent(
                        list, value, index, isSelected, hasFocus);
                }
            } );
		}

    public void addItem(String itemText) {
        addItem(itemText, UIManager.getColor("ComboBox.foreground"));
		}

    public void addItem(String itemText, Color color) {
        addItem(new ComboBoxColorItem(itemText, color));
		}

    public String getSelectedItemText() {
        if (getSelectedItem() instanceof ComboBoxColorItem)
            return ((ComboBoxColorItem)getSelectedItem()).getText();

        return getSelectedItem().toString();
        }
    }