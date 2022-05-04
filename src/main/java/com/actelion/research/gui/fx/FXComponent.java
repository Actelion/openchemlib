package com.actelion.research.gui.fx;

import com.actelion.research.gui.generic.GenericActionEvent;
import com.actelion.research.gui.generic.GenericComponent;
import com.actelion.research.gui.generic.GenericEventListener;
import javafx.scene.Node;

import java.util.ArrayList;

public class FXComponent implements GenericComponent {
	private Node mNode;
	private ArrayList<GenericEventListener<GenericActionEvent>> mConsumerList;

	public FXComponent(Node c) {
		mNode = c;
		mConsumerList = new ArrayList<>();
	}

	@Override
	public void addEventConsumer(GenericEventListener<GenericActionEvent> consumer) {
		mConsumerList.add(consumer);
	}

	@Override
	public void removeEventConsumer(GenericEventListener<GenericActionEvent> consumer) {
		mConsumerList.remove(consumer);
	}

	@Override
	public void setEnabled(boolean b) {
		if (mNode != null)
			mNode.disableProperty().set(!b);
	}

	public Node getNode() {
		return mNode;
	}

	@Override
	public void fireEvent(GenericActionEvent event) {
		for (GenericEventListener<GenericActionEvent> consumer:mConsumerList)
			consumer.eventHappened(event);
	}
}
