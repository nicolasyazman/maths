package io.github.nicolasyazman.physics;

public class Point2D {

	private double[] coordinates;
	
	public Point2D(double x, double y) {
		this.coordinates = new double[2];
		this.coordinates[0] = x;
		this.coordinates[1] = y;
	}
	
	public double GetX() {
		return this.coordinates[0];
	}
	
	public double GetY() {
		return this.coordinates[1];
	}	
	
}
