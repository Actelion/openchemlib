package com.actelion.research.gui;

import com.actelion.research.chem.DrawingObjectList;
import com.actelion.research.chem.ExtendedDepictor;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.gui.editor.SwingEditorDialog;
import com.actelion.research.gui.hidpi.HiDPIHelper;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;

public class JEditableChemistryView extends JChemistryView {
	private static final String EDIT_MESSAGE = "<double click or drag & drop>";
	private StereoMolecule[] mMolecules;
	private Reaction mReaction;

	/**
	 * Creates a new JEditableChemistryView for showing & editing a reaction or molecule(s).
	 * This default implementation will support copy/paste and drag&drop.
	 * @param chemistryType one of the ExtendedDepictor.TYPE_... options
	 */
	public JEditableChemistryView(int chemistryType) {
		super(chemistryType);
		if (chemistryType == ExtendedDepictor.TYPE_REACTION) {
			mReaction = new Reaction();
			super.setContent(mReaction);
			}
		else {
			mMolecules = new StereoMolecule[1];
			mMolecules[0] = new StereoMolecule();
			super.setContent(mMolecules);
			}
		setEditable(true);
		}

	/**
	 * If chemistryType is
	 * @return all molecules
	 */
	public StereoMolecule[] getStructures() {
		return mMolecules;
	}

	public Reaction getReaction() {
		return mReaction;
		}

	@Override
	public void setContent(StereoMolecule[] mol) {
		setContent(mol, null);
		}

	@Override
	public void setContent(Reaction rxn) {
		setContent(rxn, null);
		}

	@Override
	public void setContent(StereoMolecule[] mol, DrawingObjectList drawingObjectList) {
		mMolecules = mol;
		super.setContent(mol, drawingObjectList);
		}

	@Override
	public void setContent(Reaction rxn, DrawingObjectList drawingObjectList) {
		if (rxn == null)
			mReaction.clear();
		else
			mReaction = rxn;
		super.setContent(mReaction, drawingObjectList);
		}

	@Override
	public void mouseClicked(MouseEvent e) {
		if (e.getClickCount() == 2 && isEnabled() && isEditable()) {
			if (getChemistryType() == ExtendedDepictor.TYPE_REACTION)
				editReaction();
			else
				editMolecules();
			return;
			}

		super.mouseClicked(e);
		}

	private boolean isEmpty() {
		if (getChemistryType() == ExtendedDepictor.TYPE_REACTION)
			return mReaction == null || mReaction.isEmpty();

		if (mMolecules != null)
			for (StereoMolecule mol:mMolecules)
				if (mol.getAllAtoms() != 0)
					return false;

		return true;
		}

	@Override
	public void paintComponent(Graphics g) {
		Dimension theSize = getSize();
		Insets insets = getInsets();
		theSize.width -= insets.left + insets.right;
		theSize.height -= insets.top + insets.bottom;

		super.paintComponent(g);

		if (isEnabled() && isEditable() && isEmpty()) {
			g.setFont(g.getFont().deriveFont(Font.PLAIN, HiDPIHelper.scale(10)));
			FontMetrics metrics = g.getFontMetrics();
			Rectangle2D bounds = metrics.getStringBounds(EDIT_MESSAGE, g);
			g.drawString(EDIT_MESSAGE, insets.left+(int)(theSize.width-bounds.getWidth())/2,
					insets.top+(theSize.height-metrics.getHeight())/2+metrics.getAscent());
			}
		}

	private void editReaction() {
		SwingEditorDialog dialog = createDrawDialog("Edit Reaction", new Reaction(mReaction));
		dialog.setVisible(true);
		if (!dialog.isCancelled()) {
			Reaction newRxn = dialog.getReactionAndDrawings();

			// rescue catalysts, because the editor doesn't handle them
			if (newRxn != null) {
				for (int i=0; i<newRxn.getMolecules(); i++)
					newRxn.getMolecule(i).removeAtomSelection();
				for (int i = 0; i < mReaction.getCatalysts(); i++)
					newRxn.addCatalyst(mReaction.getCatalyst(i));
				}

			setContent(newRxn);
			informListeners();
			}
		}

	private void editMolecules() {
		StereoMolecule[] mol = new StereoMolecule[mMolecules.length];
		for (int i=0; i<mMolecules.length; i++)
			mol[i] = new StereoMolecule(mMolecules[i]);
		SwingEditorDialog dialog = createDrawDialog("Edit Molecules", mol);
		dialog.setVisible(true);
		if (!dialog.isCancelled()) {
			StereoMolecule[] newMols = dialog.getDrawArea().getFragments();
			setContent(newMols);
			informListeners();
			}
		}

	protected SwingEditorDialog createDrawDialog(String title, Reaction reaction) {
		Component c = this;
		while (!(c instanceof Frame || c instanceof Dialog))
			c = c.getParent();

		SwingEditorDialog d = (c instanceof Frame) ?
			  new SwingEditorDialog((Frame)c, reaction, Dialog.ModalityType.DOCUMENT_MODAL)
			: new SwingEditorDialog((Dialog)c, reaction, Dialog.ModalityType.DOCUMENT_MODAL);

		if (title != null)
			d.setTitle(title);

		return d;
		}

	protected SwingEditorDialog createDrawDialog(String title, StereoMolecule[] mol) {
		Component c = this;
		while (!(c instanceof Frame || c instanceof Dialog))
			c = c.getParent();

		SwingEditorDialog d = (c instanceof Frame) ?
			  new SwingEditorDialog((Frame)c, mol, Dialog.ModalityType.DOCUMENT_MODAL)
			: new SwingEditorDialog((Dialog)c, mol, Dialog.ModalityType.DOCUMENT_MODAL);

		if (title != null)
			d.setTitle(title);

		return d;
	}
}
