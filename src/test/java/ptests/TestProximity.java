package ptests;


import javax.media.jai.JAI;


import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import org.geotools.referencing.CRS;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
//import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.geotools.referencing.GeodeticCalculator;
//import org.geotools.geometry.jts.JTS;
//import org.geotools.geometry.jts.JTSFactoryFinder;
//import org.opengis.referencing.NoSuchAuthorityCodeException;
//import org.opengis.referencing.FactoryException;


/*
Latitude	Longitude
34.134782   -118.32387 - Hollywood Sign
34.13307	-118.344
33.995244	-118.454303
34.093942	-118.377523
34.0482681	-118.4627062
34.04852367	-118.3857939
34.02929123	-118.4554026
34.08537607	-118.3438039
 */

public class TestProximity {

    private static double[][] points = {
            {34.134782, 34.13307, 33.995244, 34.093942, 34.0482681, 34.04852367, 34.02929123, 34.08537607}, // lat
            {-118.32387, -118.344, -118.454303, -118.377523, -118.4627062, -118.3857939, -118.4554026, -118.3438039}  // lon
    };


    /**
     * The following will calculate the straigt line distance between two points in kilometers.
     */

    static final double R = 6372.8; // In kilometers

    /**
     * haversine Distance
     *
     * Reviewed against the following (use 6371.0 as mean earth radius)
     * https://www.vcalc.com/wiki/vCalc/Haversine+-+Distance
     *
     * @param lat1 The degrees latitude of source point.
     * @param lon1 The degrees longitude of source point.
     * @param lat2 The degrees latitude of destination point.
     * @param lon2 The degrees longitude of destination point.
     * @return The distance in kilometers.
     */
    static double haversine(double lat1, double lon1, double lat2, double lon2) {

        double dLat = Math.toRadians(lat2 - lat1);
        double dLon = Math.toRadians(lon2 - lon1);
        lat1 = Math.toRadians(lat1);
        lat2 = Math.toRadians(lat2);

        double a = Math.pow(Math.sin(dLat / 2), 2) + Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);
        double c = 2 * Math.asin(Math.sqrt(a));
        return R * c;
    }


    /**
     * For JTS tests
     *
     * Use https://www.darrinward.com/lat-long/ to map points in Google Map.
     *
     * GeodeticCalculator results compare well to
     *
     */
    @Test
    public void radialDistance() {

        System.out.println("*** Point to Point Distance ***\n");

        // NOTE:
        // EPSG:3785 - Google perfect sphere
        // EPSG:4326 - WGS84
        CoordinateReferenceSystem crs = null;
        try {
            crs = CRS.decode("EPSG:4326");
            String wkt = crs.toWKT();
            //System.out.println("wkt for EPSG:4326 ::");
            //System.out.println(wkt);
            //System.out.println();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        // NOTE: The default calculator is WGS84
        GeodeticCalculator gc = new GeodeticCalculator();
        //GeodeticCalculator gc = new GeodeticCalculator(crs);
        gc.setStartingGeographicPoint(points[1][0], points[0][0]);

        for (int i = 1; i < points[0].length; i++) {

            double distanceH = -1.0;
            double distanceO = -1.0;

            double x = points[0][i];
            double y = points[1][i];

            gc.setDestinationGeographicPoint(y, x);
            distanceO = gc.getOrthodromicDistance() / 1000;

            distanceH = haversine(points[0][0], points[1][0], x, y);

            double diff = Math.abs(distanceH - distanceO);

            NumberFormat p = new DecimalFormat("#0.0");
            NumberFormat k = new DecimalFormat("#0.0000");
            System.out.println(points[0][0] + "," + points[1][0] + "  ->  " + x + "," + y);
            // System.out.println("Distance - " + points[0][0] + "," + points[1][0] + ", " + points[0][i] + "," + points[1][i] + " :: " +
            System.out.println(
                    "haversine:  " + k.format(distanceH) + "\t" +
                    "orthodomic: " + k.format(distanceO) + "\tdiff(m): " + p.format(diff*1000) + "\n");


            //distance = JTS.orthodromicDistance( start, end, crs);
            //int totalmeters = (int) distance);
            //int km = totalmeters / 1000;
            //int meters = totalmeters - (km * 1000);
            //float remaining_cm = (float) (distanceO - totalmeters) * 10000;
            //remaining_cm = Math.round(remaining_cm);
            //float cm = remaining_cm / 100;
            //System.out.println("Distance = " + km + "km " + meters + "m " + cm + "cm");

        }
        System.out.println("");

        assertEquals(0, 0);
    }


    /**
     * Example of whether a point is in a polygon.
     *
     */
/*
    @Test
    public void pointInPolygon() {

        final GeometryFactory gf = new GeometryFactory();

        final ArrayList<Coordinate> points = new ArrayList<Coordinate>();
        points.add(new Coordinate(-10, -10));
        points.add(new Coordinate(-10, 10));
        points.add(new Coordinate(10, 10));
        points.add(new Coordinate(10, -10));
        points.add(new Coordinate(-10, -10));
        final Polygon polygon = gf.createPolygon(new LinearRing(new CoordinateArraySequence(points
                .toArray(new Coordinate[points.size()])), gf), null);

        final Coordinate coord = new Coordinate(0, 0);
        final Point point = gf.createPoint(coord);

        System.out.println(point.within(polygon));

    }
*/

}

/*


import java.util.ArrayList;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;

public class GeoTest {

  public static void main(final String[] args) {

    final GeometryFactory gf = new GeometryFactory();

    final ArrayList<Coordinate> points = new ArrayList<Coordinate>();
    points.add(new Coordinate(-10, -10));
    points.add(new Coordinate(-10, 10));
    points.add(new Coordinate(10, 10));
    points.add(new Coordinate(10, -10));
    points.add(new Coordinate(-10, -10));
    final Polygon polygon = gf.createPolygon(new LinearRing(new CoordinateArraySequence(points
        .toArray(new Coordinate[points.size()])), gf), null);

    final Coordinate coord = new Coordinate(0, 0);
    final Point point = gf.createPoint(coord);

    System.out.println(point.within(polygon));

  }

}

 */