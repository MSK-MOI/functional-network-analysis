package grapeplots;

import org.openide.util.Lookup;

import org.gephi.filters.spi.FilterProperty;
import org.gephi.graph.api.DirectedGraph;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.Graph;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.Table;
import org.gephi.graph.api.Column;

import org.gephi.graph.api.Node;
import org.gephi.filters.spi.EdgeFilter;

import org.openide.util.Exceptions;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.awt.Color;

import javax.swing.JOptionPane;

public class GRAPE implements EdgeFilter {
    private boolean c_filter=false;
    private boolean f_filter=false;
    private boolean p_filter=false;
    private GraphModel graphModel;
    private Map<String, Color> map;
    private Color gray;
    private Color[] allColors;
    private String[] annotations;
    private Boolean[] selected;

    public void setColorOf(Integer item_index, Color chosen){
        allColors[item_index-1] = chosen;
        updateColors();
    }

    public Color colorAverage(ArrayList<Color> al) {
        float[] comps = new float[3];
        for(int j=0; j<al.size(); j++) {
            float[] comps_j = al.get(j).getColorComponents(null);
            comps[0] += comps_j[0];
            comps[1] += comps_j[1];
            comps[2] += comps_j[2];
        }
        for(int i=0; i<3; i++) {
            comps[i] = comps[i] / al.size();
        }
        return new Color(comps[0], comps[1], comps[2]);
    }

    public void updateColors(){
        Table nt = graphModel.getNodeTable();

        Iterator<Node> ni = graphModel.getGraph().getNodes().iterator();
        while(ni.hasNext()) {
            Node node = ni.next();
            ArrayList<Color> al = new ArrayList<Color>();

            for(int i=0; i<300; i++) {
                if(i>=0 && i <= 99 && c_filter == false) {
                    continue;
                }
                if(i>=100 && i <= 199 && f_filter == false) {
                    continue;
                }
                if(i>=200 && i <= 299 && p_filter == false) {
                    continue;
                }
                if(!selected[i]) {
                    continue;
                }
                String node_attribute = "v_goa" + Integer.toString(i+1) + "::" + annotations[i].toLowerCase();
                String present = (String)node.getAttribute(nt.getColumn(node_attribute));
                
                if(present.equals("TRUE")) {
                    al.add(allColors[i]);
                }
            }
            if(al.size() > 0) {
                node.setColor(colorAverage(al));
            } else {
                node.setColor(gray);
            }
        }
    }

    public void clearColors(){
        for(int i=0; i<300; i++) {
            allColors[i] = new Color(0.9f, 0.9f, 0.9f);
        }
        updateColors();
    }

    public void setNewStatus(String annotation, Boolean status) {
        if(annotation == "Component") {
            c_filter=status;
        } else if(annotation == "Function") {
            f_filter=status;
        } else if(annotation == "Process") {
            p_filter=status;
        } else {
            int item_index = -1;
            for(int i=0; i<300; i++) {
                if(annotations[i].equals(annotation)) {
                    item_index = i;
                    break;
                }
            }
            selected[item_index] = status;
        }
        updateColors();
    }

    public void initializeColorsAndAnnotations(Graph graph) {
        map = new HashMap<String, Color>();

        Table nt = graphModel.getNodeTable();
        Iterator<Column> it = nt.iterator();

        while(it.hasNext()) {
            Column c = it.next();
            String line = c.getId();
            String pattern = "v_goa(\\d+)::(.*)";
            Pattern r = Pattern.compile(pattern);
            Matcher m = r.matcher(line);

            if(m.find()) {
                Integer item_index = Integer.parseInt(m.group(1));
                String annotation = m.group(2);
                annotations[item_index-1] = annotation;
            }
        }
    }

    public String[] getAnnotations() {
        return annotations;
    }

    @Override
    public boolean init(Graph graph) {
        this.graphModel = graph.getModel();

        gray = new Color(0.9f, 0.9f, 0.9f);
        allColors = new Color[300];
        for(int i=0; i<300; i++) {
            allColors[i] = gray;
        }

        annotations = new String[300];
        selected = new Boolean[300];
        for(int i=0; i<300; i++){
            annotations[i] = "no annotation found";
            selected[i] = false;
        }

        initializeColorsAndAnnotations(graph);
        return true;
    }

    @Override
    public boolean evaluate(Graph graph, Edge edge) {
        return true;
    }

    @Override
    public void finish() {
    }

    @Override
    public String getName() {
        return "GRAPE Plot";
    }

    @Override
    public FilterProperty[] getProperties() {
        FilterProperty[] filterProperties = new FilterProperty[0];
        return filterProperties;
    }
}

