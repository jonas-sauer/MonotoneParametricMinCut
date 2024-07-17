package OSM;

public class GeographicPoint implements IPoint<GeographicPoint> {

    public static final int EARTH_RADIUS_IN_CENTIMETRE = 637813700;

    private final double lat;
    private final double lon;

    public GeographicPoint(double lat, double lon) {
        this.lat = lat;
        this.lon = lon;
    }

    public double getLat() {
        return lat;
    }

    public double getLon() {
        return lon;
    }

    @Override
    public double getX() {
        return lat;
    }

    @Override
    public double getY() {
        return lon;
    }

    @Override
    public double getDistX(GeographicPoint p) {
        GeographicPoint temp = new GeographicPoint(p.getX(), this.getY());
        return getDist(temp);
    }

    @Override
    public double getDistY(GeographicPoint p) {
        GeographicPoint temp = new GeographicPoint(this.getX(), p.getY());
        return getDist(temp);
    }

    @Override
    public double getDist(GeographicPoint p) {
        return getDistance(this, p);
    }

    public static int getDistance(GeographicPoint from, GeographicPoint to) {
        double heightFrom = degreesToRadians(from.getLon());
        double heightTo = degreesToRadians(to.getLon());
        double widthFrom = degreesToRadians(from.getLat());
        double widthTo = degreesToRadians(to.getLat());
        double e = Math.acos(Math.sin(widthFrom) * Math.sin(widthTo) + Math.cos(widthFrom) * Math.cos(widthTo) * Math.cos(heightTo - heightFrom));
        return (int) Math.round(e * EARTH_RADIUS_IN_CENTIMETRE);
    }

    public GeographicPoint minimize(GeographicPoint p) {
        return new GeographicPoint(Math.min(this.getLat(), p.getLat()), Math.min(this.getLon(), p.getLon()));
    }

    public GeographicPoint maximize(GeographicPoint p) {
        return new GeographicPoint(Math.max(this.getLat(), p.getLat()), Math.max(this.getLon(), p.getLon()));
    }

    public static double degreesToRadians(double radius) {
        return radius / 180 * Math.PI;
    }

}
