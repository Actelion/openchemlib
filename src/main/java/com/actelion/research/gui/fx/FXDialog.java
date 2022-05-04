package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.*;
import javafx.scene.control.Alert;
import javafx.scene.control.ButtonType;
import javafx.scene.control.Dialog;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.RowConstraints;
import javafx.stage.Modality;
import javafx.stage.StageStyle;
import javafx.stage.Window;

public class FXDialog extends Dialog<String> implements GenericDialog {
	private GridPane mContent;

	public FXDialog(Window parent, String title) {
		initOwner(parent);
		initStyle(StageStyle.UNDECORATED);
		initModality(Modality.WINDOW_MODAL);
		setTitle(title);

		getDialogPane().getButtonTypes().addAll(ButtonType.CANCEL, ButtonType.OK);
	}

	@Override
	public void setEventConsumer(GenericEventListener<GenericActionEvent> consumer) {
		setResultConverter(dialogButton -> {
			if (dialogButton == ButtonType.OK)
				consumer.eventHappened(new GenericActionEvent(this, GenericActionEvent.WHAT_OK,0));
			else if (dialogButton == ButtonType.CANCEL)
				consumer.eventHappened(new GenericActionEvent(this, GenericActionEvent.WHAT_CANCEL,0));
			return null;
		} );
	}

	@Override
	public void setLayout(int[] hLayout, int[] vLayout) {
		mContent = new GridPane();
		for (int hl:hLayout)
			mContent.getColumnConstraints().add(hl > 0 ? new ColumnConstraints(hl) : new ColumnConstraints());
		for (int vl:vLayout)
			mContent.getRowConstraints().add(vl > 0 ? new RowConstraints(vl) : new RowConstraints());

		getDialogPane().setContent(mContent);
	}

	@Override
	public void add(GenericComponent c, int x, int y) {
		mContent.add(((FXComponent)c).getNode(), x, y);
	}

	@Override
	public void add(GenericComponent c, int x1, int y1, int x2, int y2) {
		mContent.add(((FXComponent)c).getNode(), x1, y1, x2-x1+1, y2-y1+1);
	}

	@Override
	public GenericCheckBox createCheckBox(String text) {
		return new FXCheckBox(text);
	}

	@Override
	public GenericComboBox createComboBox() {
		return new FXComboBox();
	}

	@Override
	public GenericLabel createLabel(String text) {
		return new FXLabel(text);
	}

	@Override
	public GenericTextField createTextField(int width, int height) {
		return new FXTextField(width, height);
	}

	@Override
	public void showDialog() {
		showAndWait();
	}

	@Override
	public void disposeDialog() {

	}

	@Override
	public void showMessage(String message) {
		new Alert(Alert.AlertType.INFORMATION, message).showAndWait();
	}
}
