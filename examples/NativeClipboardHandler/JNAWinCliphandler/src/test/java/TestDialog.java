
import com.actelion.research.chem.ExtendedDepictor;
import com.actelion.research.gui.JEditableChemistryView;
import com.actelion.research.gui.JEditableStructureView;
import com.actelion.research.gui.clipboard.ClipboardHandler;

import javax.swing.*;
import java.awt.*;
import java.util.Collections;

/*
 *  Manual Test Dialog to test c&p between OCL editors and ChemDraw
 */
public class TestDialog {

    public static void main(String[] args) {
        JFrame f = new JFrame("ClipboardHandler Test");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel panely = new JPanel();
        panely.setLayout(new BoxLayout(panely, BoxLayout.Y_AXIS));
        f.add(panely);

        JPanel panelx = new JPanel();
        panelx.setLayout(new BoxLayout(panelx, BoxLayout.X_AXIS));

        JEditableChemistryView natMol = new JEditableChemistryView(ExtendedDepictor.TYPE_MOLECULES);
        JEditableChemistryView natRxn = new JEditableChemistryView(ExtendedDepictor.TYPE_REACTION);
        JEditableStructureView natStruct = new JEditableStructureView();
        natStruct.setClipboardHandler(new ClipboardHandler());

        panelx.add(natMol);
        panelx.add(natStruct);
        panelx.add(natRxn);

        panely.add(panelx);


        JButton toggleButton = new JButton("Toggle Native/JNA");
        toggleButton.addActionListener(a -> {
            Collections.reverse(ClipboardHandler.getNativeCliphandler());
        });

        panely.add(toggleButton);

        f.setSize(new Dimension(800, 800));
        f.setVisible(true);

    }
}
