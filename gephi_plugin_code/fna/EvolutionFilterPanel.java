package org.msk.supernodehierarchy;

import org.gephi.filters.spi.Filter;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.awt.event.WindowEvent;
import java.awt.Dimension;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;


import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JRadioButton;
import javax.swing.ButtonGroup;
import javax.swing.ScrollPaneConstants;
import javax.swing.JFrame;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.Icon;

import org.gephi.filters.spi.FilterProperty;
import java.util.Hashtable;
import java.util.Random;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class EvolutionFilterPanel extends javax.swing.JPanel implements ChangeListener, ItemListener, ActionListener{
 
    private EvolutionFilter filter;
    private javax.swing.JCheckBox[] checkbox;
    private ColorIcon[] coloricons;
    private javax.swing.JSlider slider0;
    private javax.swing.JSlider slider1;
    private javax.swing.JSlider slider2;
    private javax.swing.JSlider slider00;
    private JScrollPane scrollPane;
    private JFrame legend;

    public EvolutionFilterPanel(EvolutionFilter filter) {
        this.filter = filter;
 
        // Boolean state = checkbox.isSelected();
        legend = new JFrame();        
        JPanel checkboxes_pane = new JPanel();
        checkboxes_pane.setLayout(new BoxLayout(checkboxes_pane, javax.swing.BoxLayout.Y_AXIS));
        // checkboxes_pane.setPreferredSize(new Dimension(450, 6500));
        scrollPane = new JScrollPane(checkboxes_pane);
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        scrollPane.setPreferredSize(new Dimension(450, 4000));
        this.setMinimumSize(new Dimension(200,700));
        this.setPreferredSize(new Dimension(450, 2100));
        scrollPane.getVerticalScrollBar().setUnitIncrement(16);

        JPanel aspect_pane = new JPanel();
        aspect_pane.setLayout(new BoxLayout(aspect_pane, javax.swing.BoxLayout.X_AXIS));
        JCheckBox option1 = new JCheckBox("Component");
        JCheckBox option2 = new JCheckBox("Function");
        JCheckBox option3 = new JCheckBox("Process");
        option1.addItemListener(this);
        option2.addItemListener(this);
        option3.addItemListener(this);
 
        aspect_pane.add(option1);
        aspect_pane.add(option2);
        aspect_pane.add(option3);

        JButton figure_generate = new JButton("Generate legend");
        figure_generate.setActionCommand("generate legend");
        figure_generate.addActionListener(this);
        figure_generate.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);

        // this.add(figure_generate);
        // this.add(aspect_pane);
        // this.add(scrollPane);
        // add(new JLabel("Component, Function, Process (100 each)"));

        // Hashtable labelTable = new Hashtable();
        // labelTable.put( new Integer( 0 ), new JLabel("-2") );
        // labelTable.put( new Integer( 100 ), new JLabel("0") );
        // labelTable.put( new Integer( 200 ), new JLabel("2") );
        setLayout(new BoxLayout(this, javax.swing.BoxLayout.Y_AXIS));

        slider00 = new javax.swing.JSlider(0, 200);
        // slider0.setLabelTable(labelTable);
        add(new JLabel("Calculated weights vs. Constant")); // vs. Correlation
        add(slider00);
 
        slider0 = new javax.swing.JSlider(0, 200);
        //slider0.setLabelTable(labelTable);
        add(new JLabel("Level"));
        add(slider0);

        slider1 = new javax.swing.JSlider(0, 200);
        //slider1.setLabelTable(labelTable);
        add(new JLabel("Bandwidth"));
        add(slider1);

        slider2 = new javax.swing.JSlider(0, 200);
        //slider2.setLabelTable(labelTable);
        add(new JLabel("Ordinary edge weight"));
        add(slider2);

        JButton grabber = new JButton("Grab annotations");
        grabber.setActionCommand("grab");
        grabber.addActionListener(this);
        grabber.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
        JPanel frow = new JPanel();
        FlowLayout f = new FlowLayout(FlowLayout.LEFT);
        f.setVgap(0);
        frow.setLayout(f);
        frow.add(new JLabel(""));
        frow.add(grabber);
        checkboxes_pane.add(frow);

        checkbox = new JCheckBox[300];
        coloricons = new ColorIcon[300];
        for(int i=0; i<300; i++) {
            String prepend = "";
            if(i == 0) {
                prepend = "Component";
            }
            if(i == 100) {
                prepend = "Function";
            }
            if(i == 200) {
                prepend = "Process";
            }
            if(i == 0 || i == 100 || i == 200) {
                JPanel header_row = new JPanel();
                f = new FlowLayout(FlowLayout.LEFT);
                f.setVgap(0);
                header_row.setLayout(f);
                header_row.add(new JLabel(""));

                JLabel header = new JLabel(prepend);
                header_row.add(header);

                checkboxes_pane.add(header_row);
            }

            checkbox[i] = new JCheckBox(" annotation", false);
            checkbox[i].addItemListener(this);

            JPanel row = new JPanel();
            FlowLayout flayout = new FlowLayout(FlowLayout.LEFT);
            flayout.setVgap(0);
            row.setLayout(flayout);
            // ColorIcon icon = new ColorIcon();
            coloricons[i] = new ColorIcon();
            Random rand = new Random();
            float r = rand.nextFloat();
            float g = rand.nextFloat();
            float b = rand.nextFloat();

            Color color_i = new Color(r, g, b);
            coloricons[i].setColor(color_i);
            filter.setColorOf(i+1, color_i);

            JButton colorpick = new JButton(coloricons[i]);
            colorpick.setPreferredSize(new Dimension(10,10));
            colorpick.setEnabled(true);
            colorpick.setActionCommand("color item "+Integer.toString(i+1));
            colorpick.addActionListener(this);

            row.add(colorpick);
            row.add(new JLabel(Integer.toString( ((i%100)+1) ) ));
            row.add(checkbox[i]);
            checkboxes_pane.add(row);
        }
 
        // checkbox.addChangeListener(this);
        //slider00.addChangeListener(this);
        slider0.addChangeListener(this);
        slider1.addChangeListener(this);
        slider2.addChangeListener(this);
        revalidate();
        repaint();
    }

    public void actionPerformed(ActionEvent e){ 
        if(e.getActionCommand().equals("grab")) {
            grabAnnotations();
            return;
        }
        if(e.getActionCommand().equals("generate legend")) {
            generateLegend();
            return;
        }
        String line = e.getActionCommand();
        String pattern = "color item (\\d+)";
        Pattern r = Pattern.compile(pattern);
        Matcher m = r.matcher(line);

        if(m.find()) {
            Integer item_index = Integer.parseInt(m.group(1));
            Color chosen = JColorChooser.showDialog(this, "Choose color for annotation " + m.group(1), new Color(0.5f,0.5f,0.5f));
            // set color item_index to be chosen
            filter.setColorOf(item_index, chosen);
            JButton pressed = (JButton)e.getSource();
            ColorIcon icon = (ColorIcon)pressed.getIcon();
            icon.setColor(chosen);
            pressed.requestFocus();
        }
    }

    public void stateChanged(ChangeEvent e) {
        double t;
        FilterProperty p0 = filter.getProperties()[0];
        t = slider0.getValue()/200.0;
        p0.setValue(filter.getMin() + t*(filter.getMax()-filter.getMin()));

        FilterProperty p1 = filter.getProperties()[1];
        t = slider1.getValue()/200.0;
        p1.setValue(filter.getBandwidthMin() + t*(filter.getBandwidthMax()-filter.getBandwidthMin()));

        FilterProperty p2 = filter.getProperties()[2];
        t = slider2.getValue()/200.0;
        p2.setValue(filter.getOrdinaryWeightMin() + t*(filter.getOrdinaryWeightMax()-filter.getOrdinaryWeightMin()));

        FilterProperty p3 = filter.getProperties()[3];
        t = slider00.getValue()/200.0;
        p3.setValue(filter.getUsingCorrelationMin() + t*(filter.getUsingCorrelationMax()-filter.getUsingCorrelationMin()));
    }

    public void itemStateChanged(ItemEvent e) {
        if(e.getStateChange() == ItemEvent.SELECTED) {
            System.out.println("selected item");
            JCheckBox c = (JCheckBox)e.getItem();
            String l = c.getLabel();
            filter.setNewStatus(l, true);
        }
        if(e.getStateChange() == ItemEvent.DESELECTED) {
            System.out.println("deselected item");
            JCheckBox c = (JCheckBox)e.getItem();
            String l = c.getLabel();
            filter.setNewStatus(l, false);
        }
    }

    public void grabAnnotations() {
        String[] annotations = filter.getAnnotations();
        String a = "";
        int i=0;
        try{
            for(i=0; i<300; i++) {
                a = annotations[i];
                JCheckBox cb = checkbox[i];
                cb.setText(a);
            }
        } catch(Exception e){
            System.err.println(Integer.toString(i));
            System.err.println(a);
            System.err.println(e);
        }
        //...
    }

    public void generateLegend() {
        legend.dispatchEvent(new WindowEvent(legend, WindowEvent.WINDOW_CLOSING));
        legend = new JFrame("Legend");
        legend.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        panel.setOpaque(true);

        for(int i=0; i<300; i++) {
            if(checkbox[i].isSelected()) {
                JPanel row = new JPanel();
                FlowLayout flayout = new FlowLayout(FlowLayout.LEFT);
                flayout.setVgap(0);
                row.setLayout(flayout);

                JButton colorpick = new JButton(new ColorIcon(30, 30, coloricons[i].getColor()));
                colorpick.setPreferredSize(new Dimension(30,30));
                colorpick.setEnabled(true);
                colorpick.setActionCommand("color item "+Integer.toString(i+1));
                colorpick.addActionListener(this);

                row.add(colorpick);
                
                JLabel label = new JLabel(checkbox[i].getText());
                label.setFont(new Font("Arial", Font.BOLD, 22));
                row.add(label);
                panel.add(row);
            }
        }

        legend.getContentPane().add(BorderLayout.CENTER, panel);
        legend.pack();
        legend.setLocationByPlatform(true);
        legend.setVisible(true);
        legend.setResizable(true);
        panel.requestFocus();        
    }
}


