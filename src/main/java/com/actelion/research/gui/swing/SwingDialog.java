package com.actelion.research.gui.swing;

import com.actelion.research.gui.generic.*;
import com.actelion.research.gui.hidpi.HiDPIHelper;
import info.clearthought.layout.TableLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class SwingDialog extends JDialog implements ActionListener,GenericDialog {
	private Component mParent;
	private JPanel  mContent;
	private GenericEventListener<GenericActionEvent> mConsumer;

	public SwingDialog(Frame parent, String title) {
		super(parent, title, true);
		mParent = parent;
		}

	public SwingDialog(Dialog parent, String title) {
		super(parent, title, true);
		mParent = parent;
		}

	@Override
	public void setEventConsumer(GenericEventListener<GenericActionEvent> consumer) {
		mConsumer = consumer;
	}

	@Override
	public void setLayout(int[] hLayout, int[] vLayout) {
		double[][] size = new double[2][];
		size[0] = new double[hLayout.length];
		size[1] = new double[vLayout.length];
		for (int i=0; i<hLayout.length; i++)
			size[0][i] = (hLayout[i] > 0) ? HiDPIHelper.scale(hLayout[i]) : hLayout[i];
		for (int i=0; i<vLayout.length; i++)
			size[1][i] = (vLayout[i] > 0) ? HiDPIHelper.scale(vLayout[i]) : vLayout[i];

		mContent = new JPanel();
		mContent.setLayout(new TableLayout(size));
		}

	@Override
	public void add(GenericComponent c, int x, int y) {
		mContent.add(((SwingComponent)c).getComponent(), x+","+y);
	}

	@Override
	public void add(GenericComponent c, int x1, int y1, int x2, int y2) {
		mContent.add(((SwingComponent)c).getComponent(), x1+","+y1+","+x2+","+y2);
	}

	@Override
	public void showDialog() {
		JPanel buttonpanel = new JPanel();
		int gap = HiDPIHelper.scale(8);
		buttonpanel.setBorder(BorderFactory.createEmptyBorder(gap*3/2, gap, gap, gap));
		buttonpanel.setLayout(new BorderLayout());
		JPanel ibp = new JPanel();
		ibp.setLayout(new GridLayout(1, 2, gap, 0));
		JButton bcancel = new JButton("Cancel");
		bcancel.addActionListener(this);
		ibp.add(bcancel);
		JButton bok = new JButton("OK");
		bok.addActionListener(this);
		ibp.add(bok);
		buttonpanel.add(ibp, BorderLayout.EAST);

		getContentPane().add(mContent, BorderLayout.CENTER);
		getContentPane().add(buttonpanel, BorderLayout.SOUTH);
		getRootPane().setDefaultButton(bok);

		pack();
		setLocationRelativeTo(mParent);
		setVisible(true);
		}

	@Override
	public void disposeDialog() {
		dispose();
		}

	@Override
	public void actionPerformed(ActionEvent e) {
		if ("OK".equals(e.getActionCommand()))
			mConsumer.eventHappened(new GenericActionEvent(this, GenericActionEvent.WHAT_OK, 0));
		else if ("Cancel".equals(e.getActionCommand()))
			mConsumer.eventHappened(new GenericActionEvent(this, GenericActionEvent.WHAT_CANCEL, 0));
		}

	@Override
	public void showMessage(String message) {
		JOptionPane.showMessageDialog(mParent, message);
		}

	@Override
	public GenericLabel createLabel(String text) {
		return new SwingLabel(text);
		}

	@Override
	public GenericTextField createTextField(int width, int height) {
		return new SwingTextField(width, height);
	}

	@Override
	public GenericCheckBox createCheckBox(String text) {
		return new SwingCheckBox(text);
	}

	@Override
	public GenericComboBox createComboBox() {
		return new SwingComboBox();
	}
	}
