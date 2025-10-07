
import com.actelion.research.chem.ExtendedDepictor;
import com.actelion.research.gui.JEditableChemistryView;
import com.actelion.research.gui.JEditableStructureView;
import com.actelion.research.gui.clipboard.ClipboardHandler;
import com.actelion.research.gui.clipboard.ClipboardHandler2;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;


public class TestDialog {

    public static void main(String[] args) {
        JFrame f = new JFrame();
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.Y_AXIS));
        f.add(p);

        JPanel nativePanel = new JPanel();
        nativePanel.setLayout(new BoxLayout(nativePanel, BoxLayout.X_AXIS));
        JPanel jnaPanel = new JPanel();
        jnaPanel.setLayout(new BoxLayout(jnaPanel, BoxLayout.X_AXIS));


        JPanel jnaPanel2 = new JPanel();
        jnaPanel2.setLayout(new BoxLayout(jnaPanel2, BoxLayout.X_AXIS));

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

        TestJEditableChemistryView jnaMol2 = new TestJEditableChemistryView(ExtendedDepictor.TYPE_MOLECULES);
        TestJEditableChemistryView jnaRxn2 = new TestJEditableChemistryView(ExtendedDepictor.TYPE_REACTION);
        JEditableStructureView jnaStruct2 = new JEditableStructureView();
        ClipboardHandler2 ch2 = new ClipboardHandler2();
        ch2.setJnaOverNative(false);
        jnaMol2.setJna(false);
        jnaRxn2.setJna(false);
        jnaStruct2.setClipboardHandler(ch2);

        jnaPanel2.add(jnaMol2);
        jnaPanel2.add(jnaStruct2);
        jnaPanel2.add(jnaRxn2);

        p.add(nativePanel);
        nativePanel.setBorder(new TitledBorder("Native"));
        jnaPanel.setBorder(new TitledBorder("JNA"));
        jnaPanel2.setBorder(new TitledBorder("Native2"));

        p.add(nativePanel);
        p.add(jnaPanel);
        p.add(jnaPanel2);

        f.setSize(new Dimension(800, 800));
        f.setVisible(true);

    }
}
