package com.actelion.research.gui.table;

import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.gui.LookAndFeelHelper;

import javax.swing.*;
import javax.swing.table.TableCellRenderer;
import java.awt.*;

public class ChemistryCellRenderer implements ListCellRenderer,TableCellRenderer {
    private ChemistryRenderPanel    mRenderPanel;
    private boolean					mIsEnabled,mAlternateBackground;

	public ChemistryCellRenderer() {
		this(null);
		}

	public ChemistryCellRenderer(Dimension preferredSize) {
		mRenderPanel = new ChemistryRenderPanel();
		if (preferredSize != null)
			mRenderPanel.setPreferredSize(preferredSize);
		}

	public void setAlternateRowBackground(boolean b) {
		mAlternateBackground = b;
		}

	public void setDisplayMode(int mode) {
		mRenderPanel.setDisplayMode(mode);
		}

	public void setTextSizeFactor(double factor) {
		mRenderPanel.setTextSizeFactor(factor);
		}

	public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean hasFocus) {
    	mIsEnabled = list.isEnabled();
    	return getCellRendererComponent(value, isSelected, hasFocus, index);
    	}

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,int row, int col) {
    	mIsEnabled = table.isEnabled();
    	return getCellRendererComponent(value, isSelected, hasFocus, row);
        }

    private Component getCellRendererComponent(Object value, boolean isSelected, boolean hasFocus, int row) {
	    if (LookAndFeelHelper.isAqua()
	     || LookAndFeelHelper.isQuaQua())
		    mRenderPanel.setOpaque(true);
	    else
		    mRenderPanel.setOpaque(false);

        if (value == null) {
            mRenderPanel.setChemistry(null);
            }
        else if (value instanceof String) {
            String s = (String)value;
            if (s.length() == 0) {
                mRenderPanel.setChemistry(null);
                }
            else {
            	// If we have a PRODUCT_IDENTIFIER we have a reaction,
            	// unless we have an idcode+SPACE+coords with coords starting with PRODUCT_IDENTIFIER.
            	int productIndex = s.indexOf(ReactionEncoder.PRODUCT_IDENTIFIER);
            	if (productIndex > 0 && s.charAt(productIndex-1) == ' ')
            		productIndex = -1;

            	if (productIndex != -1) {
            		mRenderPanel.setChemistry(ReactionEncoder.decode((String)value, true));
                	}
	            else {
	                int index = s.indexOf('\n');
	                if (index == -1) {
	                    index = s.indexOf(' ');
	                    if (index == -1)
	                        mRenderPanel.setChemistry(new IDCodeParser(true).getCompactMolecule(s));
	                    else
	                        mRenderPanel.setChemistry(new IDCodeParser(true).getCompactMolecule(
	                                                     s.substring(0, index),
	                                                     s.substring(index+1)));
	                    }
	                else {
	                	StereoMolecule mol = new StereoMolecule();
	                    new IDCodeParser(true).parse(mol, s.substring(0, index));
	                    do {
	                        s = s.substring(index+1);
	                        index = s.indexOf('\n');
	                        mol.addMolecule(new IDCodeParser(true).getCompactMolecule(index == -1 ? s : s.substring(0, index)));
	                        } while (index != -1);
	                    new CoordinateInventor().invent(mol);
	                    mRenderPanel.setChemistry(mol);
	                    }
	                }
            	}
            }
        else {
            mRenderPanel.setChemistry(value);
            }

        mRenderPanel.setAlternateBackground(mAlternateBackground && (row & 1) == 1);
        mRenderPanel.setSelected(isSelected);
        mRenderPanel.setFocus(hasFocus);
       	mRenderPanel.setOverruleForeground(mIsEnabled ? 0 : Color.GRAY.getRGB());
        return mRenderPanel;
        }
	}
