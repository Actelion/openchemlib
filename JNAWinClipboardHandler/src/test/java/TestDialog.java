
import com.actelion.research.chem.ExtendedDepictor;
import com.actelion.research.gui.JEditableChemistryView;
import com.actelion.research.gui.JEditableStructureView;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.clipboard.ClipboardHandler2;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.util.Collections;


public class TestDialog {

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.Y_AXIS));
        f.add(p);

        JPanel nativePanel = new JPanel();
        nativePanel.setLayout(new BoxLayout(nativePanel, BoxLayout.X_AXIS));
        JPanel jnaPanel = new JPanel();
        jnaPanel.setLayout(new BoxLayout(jnaPanel, BoxLayout.X_AXIS));


        JEditableChemistryView natMol = new JEditableChemistryView(ExtendedDepictor.TYPE_MOLECULES);
        JEditableChemistryView natRxn = new JEditableChemistryView(ExtendedDepictor.TYPE_REACTION);
        JEditableStructureView natStruct = new JEditableStructureView();
        natStruct.setClipboardHandler(new ClipboardHandler());

        nativePanel.add(natMol);
        nativePanel.add(natStruct);
        nativePanel.add(natRxn);

        TestJEditableChemistryView jnaMol = new TestJEditableChemistryView(ExtendedDepictor.TYPE_MOLECULES);
        TestJEditableChemistryView jnaRxn = new TestJEditableChemistryView(ExtendedDepictor.TYPE_REACTION);
        JEditableStructureView jnaStruct = new JEditableStructureView();
        jnaStruct.setClipboardHandler(new ClipboardHandler2());

        jnaPanel.add(jnaMol);
        jnaPanel.add(jnaStruct);
        jnaPanel.add(jnaRxn);

        nativePanel.setBorder(new TitledBorder("Native"));
        jnaPanel.setBorder(new TitledBorder("JNA"));

        p.add(nativePanel);
        p.add(jnaPanel);

        JButton toggleButton = new JButton("Toggle Native/JNA");
        toggleButton.addActionListener(a -> {
            Collections.reverse(ClipboardHandler2.getWindowsNativeCliphandler());
        });

        p.add(toggleButton);

        f.setSize(new Dimension(800, 800));
        f.setVisible(true);

    }
}
