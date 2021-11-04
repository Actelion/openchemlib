package com.actelion.research.gui.generic;

import com.actelion.research.gui.editor.DialogEventConsumer;

public interface GenericComponent {
	void setEnabled(boolean b);
	DialogEventConsumer getEventConsumer();
	void setEventConsumer(DialogEventConsumer consumer);
}
