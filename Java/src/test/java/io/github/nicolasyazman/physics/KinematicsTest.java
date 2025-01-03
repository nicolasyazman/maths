package io.github.nicolasyazman.physics;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import io.github.nicolasyazman.maths.Trigonometry;

public class KinematicsTest {

	@Test
	public void rotateAVector_X1_Y0_90DegreesCounterclockwiseShouldBeEqualToX0Y1() {
		Point2D p1 = new Point2D(1.0,0.0);
		Point2D res = Kinematics.RotateCounterclockwise(p1, Trigonometry.PI/2);
		assertEquals(0.0, res.GetX());
		assertEquals(1.0, res.GetY());
	}
	
	@Test
	public void rotateAVector_X1_Y0_180DegreesCounterclockwiseShouldBeEqualToXMinus1Y0() {
		Point2D p1 = new Point2D(1.0,0.0);
		Point2D res = Kinematics.RotateCounterclockwise(p1, Trigonometry.PI);
		assertEquals(-1.0, res.GetX());
		assertEquals( 0.0, res.GetY());
	}
	
	@Test
	public void rotateAVectorX1Y045DegreesCounterclockwiseShouldBeEqualToXMinus1Y0() {
		Point2D p1 = new Point2D(1.0,0.0);
		Point2D res = Kinematics.RotateClockwise(p1, Trigonometry.PI / 4);
		assertEquals(0.707106781186547, res.GetX(), Math.pow(10, -10));
		assertEquals(0.707106781186547, res.GetY(), Math.pow(10, -10));
	}
	
	
}
