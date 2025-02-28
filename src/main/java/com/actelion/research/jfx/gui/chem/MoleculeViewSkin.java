/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/
package com.actelion.research.jfx.gui.chem;

import com.actelion.research.chem.AbstractDepictor;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.gui.fx.FXDrawContext;
import com.actelion.research.gui.generic.GenericDepictor;
import com.actelion.research.gui.generic.GenericDrawContext;
import com.actelion.research.gui.generic.GenericRectangle;
import com.actelion.research.jfx.gui.misc.ClipboardHelper;
import com.actelion.research.jfx.gui.misc.Selector;
import com.actelion.research.util.ColorHelper;
import javafx.beans.InvalidationListener;
import javafx.beans.property.ObjectProperty;
import javafx.beans.property.SimpleObjectProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.ObservableList;
import javafx.event.ActionEvent;
import javafx.event.Event;
import javafx.event.EventHandler;
import javafx.geometry.Bounds;
import javafx.geometry.Point2D;
import javafx.scene.Node;
import javafx.scene.Parent;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Skin;
import javafx.scene.input.*;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;

/**
 * Project:
 * User: rufenec
 * Date: 10/12/11
 * Time: 4:10 PM
 */


public class MoleculeViewSkin //extends SkinBase<MoleculeView,MoleculeViewBehavior> //, BehaviorBase<MoleculeView>>
    implements Skin<MoleculeView>, IMoleculeView
{
    public static final int DEFAULT_IMAGE_WIDTH = 200;
    public static final int DEFAULT_IMAGE_HEIGHT = 200;
    public static final Color DEFAULT_BG = new Color(1, 1, 1, 1);

    private Canvas canvas = new Canvas();
    private MoleculeView control = null;
    private Canvas dragCanvas = null;
    private GenericDepictor dragDepictor;
    private Pane glassPane;
    private Color backgroundColor = DEFAULT_BG;
    private Color foregroundColor = null;   // determine automatically
    private Color overruleForeground,overruleBackground;

    private boolean draggingStartedHere = false;

    private Pane rootNode = new Pane()
    {
        protected void layoutChildren()
        {
            super.layoutChildren();
            draw();
        }
    };
    private Color borderColor = Color.GRAY;


    public MoleculeViewSkin(MoleculeView control)
    {
        this.control = control;
        setValue(new StereoMolecule());
        initControl();
    }


    private void initControl()
    {
        valueProperty().addListener(new ChangeListener<StereoMolecule>()
        {
            @Override
            public void changed(ObservableValue<? extends StereoMolecule> ov, StereoMolecule t, StereoMolecule t1)
            {
                //draw();
//                System.err.println("Request Layout");
                rootNode.requestLayout();
            }
        });

        rootNode.setOnDragEntered(new EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent dragEvent)
            {
                // noop
            }
        });
        rootNode.setOnDragDetected(new EventHandler<MouseEvent>()
        {
            @Override
            public void handle(MouseEvent mouseEvent)
            {
                draggingStartedHere = true;
                setupGlassPane();
                if (glassPane != null) {    // if() added to prevent NullPointerExceptions TLS 4-Aug-2017
                    setupDragImage(mouseEvent);
                    setupDragDropContent(mouseEvent);
                }
            }
        });

        rootNode.setOnDragDropped(new EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent event)
            {
                if (event.getSource() != event.getTarget()) {
                    Dragboard db = event.getDragboard();
                    boolean success = control.putContent(db);
                    /* let the source know whether the string was successfully
                     * transferred and used */
                    event.setDropCompleted(success);
                    event.consume();
                }
            }
        });

        rootNode.setOnDragDone(new EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent dragEvent)
            {
                draggingStartedHere = false;
            }
        });

        control.setOnDragDone(new EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent dragEvent)
            {
                glassPane.setMouseTransparent(true);

                releaseGlassPane();
            }
        });

        rootNode.setOnDragOver(new EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent event)
            {
                /* data is dragged over the target */
                Dragboard db = event.getDragboard();
                if (event.getGestureSource() != this && ClipboardHelper.getAcceptedFormats(db).size() > 0) {
                    /* allow for both copying and moving, whatever user chooses */
                    event.acceptTransferModes(TransferMode.COPY_OR_MOVE);
                }
                event.consume();
            }
        });
// Test only
//        rootNode.setStyle("-fx-background-color: rgba(255,0,0,.5);");
    }

    @Override
    public MoleculeView getSkinnable()
    {
        return control;
    }

    @Override
    public Node getNode()
    {
        return rootNode;
    }

    public void dispose()
    {
    }

    public void draw()
    {
        double w = control.getWidth();
        double h = control.getHeight();
        GraphicsContext ctx = (canvas.getGraphicsContext2D());
        rootNode.getChildren().setAll(canvas);
        canvas.setWidth(w);
        canvas.setHeight(h);
        drawBackground(ctx, w, h);
        drawMolecule(ctx, w, h);
        drawBorder(ctx, w, h);
    }

    public void setMolecule(StereoMolecule mol)
    {
        setValue(mol);
        draw();
    }

    private void drawMolecule(GraphicsContext ctx, double w, double h)
    {
        double d = getSkinnable().sizeContent() ? Math.min(w, h) / AbstractDepictor.cOptAvBondLen * 2 : 0;

        AbstractDepictor depictor = getSkinnable().createDepictor(getValue());
        if (overruleForeground != null) {
            java.awt.Color fg = new java.awt.Color((float)overruleForeground.getRed(), (float)overruleForeground.getGreen(), (float)overruleForeground.getBlue());
            java.awt.Color bg = (overruleBackground == null) ? null :
                    new java.awt.Color((float)overruleBackground.getRed(), (float)overruleBackground.getGreen(), (float)overruleBackground.getBlue());
            depictor.setOverruleColor(fg, bg);
        }
        else {
            java.awt.Color bg = new java.awt.Color((float)backgroundColor.getRed(), (float)backgroundColor.getGreen(), (float)backgroundColor.getBlue());
            java.awt.Color fg = (ColorHelper.perceivedBrightness(bg)>0.5) ? java.awt.Color.BLACK : java.awt.Color.WHITE;
            depictor.setForegroundColor(fg, bg);
        }
        GenericDrawContext context = new FXDrawContext(ctx);
	    depictor.validateView(context, new GenericRectangle(0, 0, (float) w, (float) h), AbstractDepictor.cModeInflateToMaxAVBL + (int) (d));
        depictor.paint(context);
    }


    private void drawBackground(GraphicsContext ctx, double w, double h)
    {
        ctx.save();
        ctx.clearRect(0, 0, w, h);
        if (!control.isTransparent()) {
            ctx.setFill(getBackgroundColor());
            ctx.fillRect(0, 0, w, h);
        } else {
//            System.out.println("Controle is Transparent");
        }
        ctx.restore();
    }

    private void drawBorder(GraphicsContext ctx, double w, double h)
    {
        if (borderColor != null) {
            ctx.save();
            ctx.setStroke(borderColor);
            ctx.strokeRect(0, 0, w, h);
            ctx.restore();
        }
    }

    public Color getBorderColor()
    {
        return borderColor;
    }

    public void setBorderColor(Color borderColor)
    {
        this.borderColor = borderColor;
    }


    public ObjectProperty<StereoMolecule> valueProperty()
    {
        return value;
    }

    public final void setValue(StereoMolecule value)
    {
        valueProperty().set(value);
    }

    public final StereoMolecule getValue()
    {
        return valueProperty().get();
    }

    void valueInvalidated()
    {
        control.fireEvent(new ActionEvent());
    }

    private ObjectProperty<StereoMolecule> value = new SimpleObjectProperty<StereoMolecule>(this, "value")
    {
        StereoMolecule oldValue;
        private ObservableValue<? extends StereoMolecule> observable = null;
        private InvalidationListener listener = null;

        @Override
        protected void invalidated()
        {
            super.invalidated();
            StereoMolecule newValue = get();

/*
            if ((oldValue == null && newValue != null) ||
                    oldValue != null && !oldValue.equals(newValue))
*/
            {
                valueInvalidated();
            }

            oldValue = newValue;
        }

        @Override
        public void set(StereoMolecule molecule)
        {
            super.set(molecule);
        }

        @Override
        public void bind(ObservableValue<? extends StereoMolecule> observable)
        {
            //System.out.printf("bind()\n");
            super.bind(observable);
        }

        @Override
        public void unbind()
        {
            //System.out.printf("Super.unbind()\n");
            super.unbind();
        }
    };

    protected void registerDragEvent()
    {
        if (dragCanvas == null)
            dragCanvas = new Canvas();
        dragCanvas.setHeight(DEFAULT_IMAGE_WIDTH);
        dragCanvas.setWidth(DEFAULT_IMAGE_HEIGHT);

        //Add the drag over listener to this component's PARENT, so that the drag over events will be processed even
        //after the cursor leaves the bounds of this component.
        glassPane.setOnDragOver(new javafx.event.EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent e)
            {
                if (draggingStartedHere) {
                    //                Point2D localPoint = getControl().getScene().getRoot().sceneToLocal(new Point2D(e.getSceneX(), e.getSceneY()));
                    Point2D localPoint = glassPane.sceneToLocal(new Point2D(e.getSceneX(), e.getSceneY()));
                    dragCanvas.relocate(
                        (int) (localPoint.getX() - dragCanvas.getBoundsInLocal().getWidth() / 2),
                        (int) (localPoint.getY() - dragCanvas.getBoundsInLocal().getHeight() / 2));
                    GraphicsContext ctx = (dragCanvas.getGraphicsContext2D());
                    ctx.clearRect(0, 0, dragCanvas.getWidth(), dragCanvas.getHeight());
                    dragDepictor.paint(new FXDrawContext(ctx));
                    Node n = findDragParentNode(e);
                    if (n != null) {
                        Event.fireEvent(n, e);
                    }
                }
            }
        });

        glassPane.setOnDragDropped(new EventHandler<DragEvent>()
        {
            @Override
            public void handle(DragEvent e)
            {
                Node n = findDragParentNode(e);
                if (n != null) {
                    Event.fireEvent(n, e);
                }
                //To change body of implemented methods use File | Settings | File Templates.
            }
        });

    }

    private Node findDragParentNode(DragEvent e)
    {
        return findLeaf(control.getScene().getRoot(), new Point2D(e.getSceneX(), e.getSceneY()), new Selector<Node>()
        {

            @Override
            public boolean match(Node node)
            {
                return node.getOnDragOver() != null && node.isVisible() && !(node instanceof DragGlassPane);
            }
        });
    }

    private Node findLeaf(Parent p, Point2D pt, Selector<Node> sel)
    {
        ObservableList<Node> childs = p.getChildrenUnmodifiable();
        for (Node n : childs) {
            Bounds b = n.localToScene(n.getBoundsInLocal());
            if (n != glassPane && b.contains(pt)) {
                if (sel.match(n))
                    return n;
                if (n instanceof Parent) {
                    p = (Parent) n;
                    n = findLeaf(p, pt, sel);
                    if (n != null)
                        return n;
                }
            }
        }
        return null;
    }

    private void setupDragDropContent(MouseEvent mouseEvent)
    {
        Dragboard db = control.startDragAndDrop(TransferMode.ANY);
        ClipboardContent content = control.getContent(mouseEvent);
        db.setContent(content);
        mouseEvent.consume();
    }

    private void setupDragImage(MouseEvent mouseEvent)
    {
        StereoMolecule mol = control.getValue();
        if (mol != null) {
            glassPane.setMouseTransparent(false);
            if (!glassPane.getChildren().contains(dragCanvas)) {
                glassPane.getChildren().add(dragCanvas);
            }
            dragDepictor = new GenericDepictor(mol);
            dragDepictor.validateView(new FXDrawContext(canvas.getGraphicsContext2D()), new GenericRectangle(0, 0, DEFAULT_IMAGE_WIDTH, DEFAULT_IMAGE_HEIGHT), JFXCanvasDepictor.cModeInflateToMaxAVBL);
            dragCanvas.setOpacity(0.7);
            dragCanvas.toFront();
            dragCanvas.setMouseTransparent(true);
            dragCanvas.setVisible(true);
            dragCanvas.relocate(
                (int) (mouseEvent.getSceneX() - dragCanvas.getBoundsInLocal().getWidth() / 2),
                (int) (mouseEvent.getSceneY() - dragCanvas.getBoundsInLocal().getHeight() / 2));
        }
    }

    private void setupGlassPane()
    {
        if (glassPane == null) {
            if (rootNode.getScene() != null) {
                final Parent parent = rootNode.getScene().getRoot();
                glassPane = new DragGlassPane();
                glassPane.setOpacity(0.5);
                glassPane.setMouseTransparent(true);
                if (parent instanceof Pane) {
                    Pane p = (Pane) parent;
                    p.getChildren().add(glassPane);
                }
                registerDragEvent();
            }
        }
    }

    private void releaseGlassPane()
    {
        dragCanvas.setVisible(false);
        dragCanvas = null;

        if (glassPane != null) {
            if (rootNode.getScene() != null) {
                final Parent parent = rootNode.getScene().getRoot();
                if (parent instanceof Pane) {
                    Pane p = (Pane) parent;
                    p.getChildren().remove(glassPane);

                    glassPane.setOnDragOver(null);
                    glassPane.setOnDragDropped(null);
                    glassPane = null;
                }
            }
        }
    }

    private StringProperty idcode = new SimpleStringProperty()
    {
        @Override
        public String get()
        {
            StereoMolecule mol = value.get();
            if (mol != null && mol.getAllAtoms() > 0)
                return mol.getIDCode();
            return null;
        }
    };

    public StringProperty idcodeProperty()
    {
        return idcode;
    }

    public Color getBackgroundColor()
    {
        return backgroundColor;
    }

    public void setBackgroundColor(Color backgroundColor)
    {
        this.backgroundColor = backgroundColor;
    }

    /**
     * Colors to be passed to Depictor. If the foreground is set, then all bonds and C- and H-atoms
     * are drawn in the foreground color. All other atoms are drawn in the atomicNo specific color,
     * which is adapted to have a minimum contrast on the given background.
     * @param foregroundColor null to use black or white, whatever contrasts better with background
     * @param backgroundColor
     */
    public void setColors(Color foregroundColor, Color backgroundColor)
    {
        this.foregroundColor = foregroundColor;
        this.backgroundColor = backgroundColor;
    }

    /**
     * Colors to be passed to Depictor. If the overrule foreground is set, then the entire molecule
     * is drawn in this color. The overrule background color is used to construct a proper bond
     * hilite color, if bond hiliting is used.
     * @param overruleForeground null to use foreground defined by setColors() or automatic foreground
     * @param overruleBackground
     */
    public void setOverruleColors(Color overruleForeground, Color overruleBackground)
    {
        this.overruleForeground = overruleForeground;
        this.overruleBackground = overruleBackground;
    }
}


