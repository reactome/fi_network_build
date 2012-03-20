/**
 * 
 */
package org.reactome.panther;

/**
 * Helper class to create a Box with x, y coordinates
 * of the origin, as well as width and height
 * @author andreash
 */
public class Box {
	private Integer x;
	private Integer y;
	private Integer w;
	private Integer h;

	/**
	 * Creates a Box and sets the coordinates for it
	 * @param x The x-value
	 * @param y The y-value
	 * @param w The width
	 * @param h The height
	 */
	public Box(Integer x, Integer y, Integer w, Integer h) {
		this.x = x;
		this.y = y;
		this.w = w;
		this.h = h;
	}
	
	/**
	 * Creates an empty Box
	 */
	public Box() {
	}
	
	/**
	 * Returns the String representation of a Box with x, y, width and height
	 * @return The String representation
	 */
	public String toString() {
		return "Box[x="+getX()+", y="+getY()+", w="+getW()+", h="+getH()+"]";
	}
	
	/**
	 * Returns the x-Value
	 * @return The x-Value
	 */
	public Integer getX() {
		return x;
	}

	/**
	 * Returns the y-Value
	 * @return The y-Value
	 */
	public Integer getY() {
		return y;
	}

	/**
	 * Returns the width
	 * @return The width
	 */
	public Integer getW() {
		return w;
	}

	/**
	 * Returns the height
	 * @return The height
	 */
	public Integer getH() {
		return h;
	}

	/**
	 * Sets the x-value
	 * @param x The x-value
	 */
	public void setX(int x) {
		this.x = x;
	}
	
	/**
	 * Sets the y-value
	 * @param y The y-value
	 */
	public void setY(int y) {
		this.y = y;
	}
	
	/**
	 * Sets the width
	 * @param w The width
	 */
	public void setW(int w) {
		this.w = w;
	}
	
	/**
	 * Sets the height
	 * @param h The height
	 */
	public void setH(int h) {
		this.h = h;
	}
}
