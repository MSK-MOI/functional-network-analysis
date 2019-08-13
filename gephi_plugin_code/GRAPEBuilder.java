package grapeplots;

import javax.swing.Icon;
import javax.swing.JPanel;
import org.gephi.filters.api.FilterLibrary;
import org.gephi.filters.spi.Category;
import org.gephi.filters.spi.Filter;
import org.gephi.filters.spi.FilterBuilder;
import org.openide.util.lookup.ServiceProvider;

import org.gephi.project.api.Workspace;

/**
 * Filter builder for the {@link GRAPE} filter.
 * <p>
 * This class configures how the filter should be integrated. It specifies it
 * belongs to the edge category, the name, icon and description.
 * <p>
 * This example doesn't have any user interface so the <code>getPanel()</code>
 * returns null.
 * 
 * @author Mathieu Bastian
 */
@ServiceProvider(service = FilterBuilder.class)
public class GRAPEBuilder implements FilterBuilder {

    @Override
    public Category getCategory() {
        return FilterLibrary.TOPOLOGY;
    }

    @Override
    public String getName() {
        return "GRAPE Plot";
    }

    @Override
    public Icon getIcon() {
        return null;
    }

    @Override
    public String getDescription() {
        return "GRAPE Plot";
    }

    @Override
    public Filter getFilter(Workspace ws) {
        return new GRAPE();
    }

    @Override
    public JPanel getPanel(Filter filter) {
        return new GRAPEPanel((GRAPE)filter);
    }

    @Override
    public void destroy(Filter filter) {
    }
}
