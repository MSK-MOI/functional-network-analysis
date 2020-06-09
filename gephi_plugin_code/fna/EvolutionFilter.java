package org.msk.supernodehierarchy;

// import org.gephi.visualization.api.VisualizationController;
// import org.gephi.visualization.*;
// import org.gephi.visualization.api.*;
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

public class EvolutionFilter implements EdgeFilter {
    private boolean initialized=false;
    private double level=0.0;
    private double bandwidth=0.1;
    private double ordinaryweight=0.1;
    private double usingcorrelation=0.0;
    private boolean c_filter=false;
    private boolean f_filter=false;
    private boolean p_filter=false;
    private double mean_vw;
    private double std_vw;
    private double mean_c;
    private double std_c;
    private GraphModel graphModel;
    private Map<String, Color> map;
    private Color gray;
    private Color[] allColors;
    private String[] annotations;
    private Boolean[] selected;

    public void setColorOf(Integer item_index, Color chosen){
        allColors[item_index-1] = chosen;
        if(initialized) {
            updateColors();
        }
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

        // String allcolumns = "";
        // Iterator<Column> ii = nt.iterator();
        // while(ii.hasNext()){
        //     allcolumns+=ii.next().getId() + " \n";
        // }
        // JOptionPane.showMessageDialog(null, allcolumns);

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

    public double getMin() {
        return -3;
    }

    public double getMax() {
        return 7;
    }

    public double getBandwidthMin() {
        return 0.0001;
    }

    public double getBandwidthMax() {
        return 3.5;
    }

    public double getOrdinaryWeightMin() {
        return 0;
    }

    public double getOrdinaryWeightMax() {
        return 1;
    }

    public double getUsingCorrelationMin() {
        return 0;
    }

    public double getUsingCorrelationMax() {
        return 1;
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

        // map.put("binaryattribute_...", new Color(0.1f, 0.1f, 0.1f));
    }

    public String[] getAnnotations() {
        return annotations;
    }

    public void setParameters(Graph graph) {
        // if(initialized) {
        //     return;
        // }

        Table et = graphModel.getEdgeTable();
        Column typecol;
        try {
            if(et.hasColumn("e_type")) {
                typecol = et.getColumn("e_type");
            } else {
                throw new Exception("The type column is not present in graph model. There are "+Integer.toString(et.countColumns())+" columns.");
            }
            Iterator<Edge> ei = graph.getEdges().iterator();
            ArrayList<Double> al = new ArrayList();
            ArrayList<Double> al2 = new ArrayList();
            while(ei.hasNext()) {
                Edge edge = ei.next();
                String type = (String)edge.getAttribute(typecol);
                if(type.equals(new String("not original/virtual")) || type.equals(new String("original/virtual"))) {
                    al.add((Double)edge.getAttribute("e_average_distance"));
                }
                // if(type.equals(new String("original/virtual")) || type.equals(new String("original/not virtual"))) {
                //     al2.add((Double)edge.getAttribute("e_correlations"));
                // }
            }
            Double[] fa = new Double[1];
            fa[0] = 0.0;
            BasicStatistics stats = new BasicStatistics(al.toArray(fa));
            mean_vw = (Float)(float)(double)stats.getMean();
            std_vw = (Float)(float)(double)stats.getStdDev();
            // BasicStatistics stats2 = new BasicStatistics(al2.toArray(fa));
            // mean_c = (Float)(float)(double)stats.getMean();
            // std_c = (Float)(float)(double)stats.getStdDev();

        } catch(Exception e) {
            System.out.println(e);
        }
        initialized=true;
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

        setParameters(graph);
        initializeColorsAndAnnotations(graph);


        // VisualizationController controller = Lookup.getDefault().lookup(VisualizationController.class);
        // SelectionManager sm = controller.getSelectionManager();
        

        return graph instanceof DirectedGraph; //? Why not just true?
    }

    @Override
    public boolean evaluate(Graph graph, Edge edge) {
        if(edge.getAttribute("e_type").equals("not original/virtual") || edge.getAttribute("e_type").equals("original/virtual")) {
            return true;
        }

        // if(edge.getAttribute("e_type").equals("original/not virtual") || edge.getAttribute("e_type").equals("original/virtual")) {
        //     return true;
        // }
        return false;
    }

    @Override
    public void finish() {
    }

    @Override
    public String getName() {
        return "Feature Network Reduction";
    }

    @Override
    public FilterProperty[] getProperties() {
        FilterProperty[] filterProperties = new FilterProperty[0];
        try {
            filterProperties = new FilterProperty[]{
            FilterProperty.createProperty(this, Double.class, "level"),
            FilterProperty.createProperty(this, Double.class, "bandwidth"),
            FilterProperty.createProperty(this, Double.class, "ordinaryweight"),
            FilterProperty.createProperty(this, Double.class, "usingcorrelation"),
            };
            // FilterProperty.createProperty(this, Boolean.class, "usingcorrelation")
        }
        catch(NoSuchMethodException e) {
            Exceptions.printStackTrace(e);
        }
        return filterProperties;
    }

    public Double getLevel() {
        return level;
    }

    public Double getBandwidth() {
        return bandwidth;
    }

    public Double getOrdinaryweight() {
        return ordinaryweight;
    }

    public Double getUsingcorrelation() {
        return usingcorrelation;
    }

    public void setLevel(Double level) {
        this.level = level;
        updateWeights();
    }

    public void setBandwidth(Double bandwidth) {
        this.bandwidth = bandwidth;
        updateWeights();
    }

    public void setOrdinaryweight(Double ordinaryweight) {
        this.ordinaryweight = ordinaryweight;
        updateWeights();
    }

    public void setUsingcorrelation(Double usingcorrelation) {
        this.usingcorrelation = usingcorrelation;
        updateWeights();
    }

    public void updateWeights() {
        Table et = graphModel.getEdgeTable();
        Graph g = graphModel.getGraph();
        Edge[] edges = g.getEdges().toArray();

        Column typecol;
        Column distcol;
        // Column corcol;
        try{
            if(et.hasColumn("e_type") && et.hasColumn("e_average_distance")) {  //&& et.hasColumn("e_correlations")
                typecol = et.getColumn("e_type");
                distcol = et.getColumn("e_average_distance");
                // corcol = et.getColumn("e_correlations");
            } else {
                String runner = "";
                Iterator<Column> ci = et.iterator();
                while(ci.hasNext()) {
                    Column c = ci.next();
                    runner+="c.getId() returns "+c.getId()+"\n";
                    runner+="c.getTitle() returns "+c.getTitle()+"\n";
                }
                throw new Exception("The type or average_distance column is not present in graph model. There are "+Integer.toString(et.countColumns())+" columns."+runner);
            }

             // if(usingcorrelation < 0.5) {
                double running_sum = 0.0;
                int virtual_count = 0;
                for(Edge e : edges)
                {
                    String type = (String)e.getAttribute(typecol);
                    if(type.equals("not original/virtual") || type.equals("original/virtual")) {
                        double virtual_weight = (double)(Double)e.getAttribute(distcol);
                        // double transformed =1.0- 1.0/ (1 + Math.exp(- Math.pow((virtual_weight-mean_vw)/std_vw- level*1.0,2)/bandwidth) );
                        double transformed =5*(  1.0- 1.0/ (1 + Math.exp(- Math.pow((virtual_weight-mean_vw)/std_vw- level*1.0,2)/bandwidth) )  );

                        // transformed = 1*(usingcorrelation) + transformed*(1-usingcorrelation);
                        e.setAttribute("weight", (Double)transformed);
                        running_sum += transformed;
                        virtual_count++;
                    }
                    // else {
                    //     e.setAttribute("weight", (Double)ordinaryweight);
                    // }
                }
                // double mean = running_sum/virtual_count;
                // for(Edge e : edges) {
                //     String type = (String)e.getAttribute(typecol);
                //     if(type.equals("not original/virtual") || type.equals("original/virtual")) {
                //         double transformed = (double)(Double)e.getAttribute("weight");
                //         e.setAttribute("weight", 1+(Double)(transformed/mean));
                //     }
                //     else {
                //         e.setAttribute("weight", 1);
                //     }
                // }
            // } else {
            //     for(Edge e : edges) {
            //         String type = (String)e.getAttribute(typecol);
            //         double correlation = (double)(Double)e.getAttribute(corcol);
            //         if(type.equals("original/virtual") || type.equals("original/not virtual")) {
            //             double transformed =1.0- 1.0/ (1 + Math.exp(- ((correlation-mean_c)/std_c - level*1.0)/bandwidth) );
            //             e.setAttribute("weight", (Double)transformed);
            //         }
            //     }
            // }
        } catch(Exception e) {
            System.out.println(e);
        }
    }
}



// Integer timeslice = level;
// Integer synthesis = (Integer)node.getAttribute("SynthesisTime");
// Integer absorption = (Integer)node.getAttribute("AbsorptionTime");
// if(synthesis <= timeslice && (absorption > timeslice || absorption == -1)) {
//     return true;
// }
// return false;
