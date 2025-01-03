package io.github.nicolasyazman.physics;

import com.aparapi.Kernel;
import com.aparapi.Range;

import io.github.nicolasyazman.maths.Matrix;
import io.github.nicolasyazman.maths.Trigonometry;

public class Kinematics {


	/**
	 * Gives the unitary rotation matrix corresponding to a rotation angle.
	 * Helper function.
	 * @param theta The rotation angle, in radians.
	 * @param clockwise If true, the 2D rotation matrix is clockwise, otherwise it is counterclockwise.
	 * @return A rotation matrix 
	 */
	public static double[][] rotationMatrix2D(double theta, boolean clockwise) {
		double cosTheta = Trigonometry.cos(theta);
		double sinTheta = Trigonometry.sin(theta);
		
		if (clockwise == true) {
			sinTheta = -sinTheta;
		}
		
		return new double[][] {{cosTheta, -sinTheta}, 
				{sinTheta, cosTheta}};
	}
	
	/**
	 * Performs a counterclockwise rotation of a 2D Point.
	 * Warning: The center of rotation is expected to be 0.0.
	 * If you want another center of rotation, you can use the function RotateCounterClockwise(Point2D, double theta, Point2D center) instead.
	 * @param p The original point.
	 * @param theta The angle of rotation (in radians).
	 * @return A new 2D Point with the rotation applied.
	 */
	public static Point2D RotateCounterclockwise(Point2D p, double theta) {
		double pointOriginalCoordinates[][] = {{p.GetX()}, {p.GetY()}};
		double[][] rotMat = rotationMatrix2D(theta, false);
		
		double[][] res = Matrix.MatMul(rotMat, pointOriginalCoordinates);
		
		return new Point2D(res[0][0], res[1][0]);
	}
	
	/**
	 *
	 * Example usage: You have an image and you want to rotate it around itself, not around the Pixel [0,0] of the screen.
	 * @param p
	 * @param theta
	 * @param centerOfRotation
	 * @return
	 */
	public static Point2D RotateCounterclockwise(Point2D p, double theta, Point2D centerOfRotation) {
		double pointOriginalCoordinates[][] = {{p.GetX() - centerOfRotation.GetX()}, {p.GetY() - centerOfRotation.GetY()}};
		double[][] rotMat = rotationMatrix2D(theta, false);
		
		double[][] res = Matrix.MatMul(rotMat, pointOriginalCoordinates);
		
		return new Point2D(res[0][0] + centerOfRotation.GetX(), res[1][0] + centerOfRotation.GetY());
	}
	
	
	/**
	 * Performs a clockwise rotation of a 2D Point.
	 * Warning: The center of rotation is expected to be 0.0.
	 * If you want another center of rotation, you can use the function RotateClockwise(Point2D, double theta, Point2D center) instead.
	 * @param p The original point.
	 * @param theta The angle of rotation (in radians).
	 * @return A new 2D Point with the rotation applied.
	 */
	public static Point2D RotateClockwise(Point2D p, double theta) {
		double pointOriginalCoordinates[][] = {{p.GetX()}, {p.GetY()}};
		double[][] rotMat = rotationMatrix2D(theta, true);
		
		double[][] res = Matrix.MatMul(rotMat, pointOriginalCoordinates);
		
		return new Point2D(res[0][0], res[1][0]);
	}
	
	public static void GPURotateImage(byte[] image, int rows, int cols, double theta) {
		
		final int size = rows * cols;
		
        final double[][] RotMat = rotationMatrix2D(theta, false);
        
        final double[] sum = new double[size];
        final double[] xRot = new double[size];
        final double[] yRot = new double[size];
        
		Kernel kernel = new Kernel(){
	           @Override public void run() {
	              int gid = getGlobalId();
	              double y = (int)(gid / cols);
	              double x = (int)(gid % cols);
	              double currentRotation = Trigonometry.atan2(y, x);
	              double newRotation = currentRotation+theta;
	              
	              xRot[gid] = x * RotMat[];
	           }
	        };
	        
	        kernel.execute(Range.create(size));
	}
	
}
