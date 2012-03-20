/**
 * 
 */
package org.reactome.panther;

/**
 * A helper class for source and target coordinates
 * @author andreash
 *
 */
public class Line {
	private int sourceX;
	private int sourceY;
	private int targetX;
	private int targetY;
	
	private int nrSource;
	private int nrTarget;
	
	/**
	 * Creates a Line object and initialises all values with <code>zero</code>, a point with the coordinates (0,0)
	 */
	public Line() {
		sourceX = 0;
		sourceY = 0;
		targetX = 0;
		targetY = 0;
		this.nrSource = 0;
		this.nrTarget = 0;
	}
	
	/**
	 * Creates a Line object and initialises the values with the given parameters
	 * @param sourceX The x-value of the source point
	 * @param sourceY The y-value of the source point
	 * @param targetX The x-value of the target point
	 * @param targetY The y-value of the target point
	 */
	public Line(int sourceX, int sourceY, int targetX, int targetY) {
		this.sourceX = sourceX;
		this.sourceY = sourceY;
		this.targetX = targetX;
		this.targetY = targetY;
		this.nrSource = 1;
		this.nrTarget = 1;
	}
	
	/**
	 * Adds a Point to the source coordinates
	 * @param sourceX The x-value of the source point
	 * @param sourceY The y-value of the source point
	 */
	public void addToSource(int sourceX, int sourceY) {
		this.sourceX += sourceX;
		this.sourceY += sourceY;
		nrSource++;
	}
	
	/**
	 * Adds a Point to the target coordinates
	 * @param targetX The x-value of the target point
	 * @param targetY The y-value of the target point
	 */
	public void addToTarget(int targetX, int targetY) {
		this.targetX += targetX;
		this.targetY += targetY;
		nrTarget++;
	}
	
	/**
	 * Returns the String representation of the Line, with the averaged source and target point
	 * @return The String representation of the Line object
	 */
	public String toString() {
		return "[source: {x:"+new Integer(getAvgX()).toString()+"|y:"+new Integer(getAvgY()).toString()+"}" +
				"target: {x:"+new Integer(getAvgX()).toString()+"|y:"+new Integer(getAvgY()).toString()+"}]";
	}
	
	/**
	 * Returns the averaged x-value of the source point
	 * @return The x-Value of the source point
	 */
	public int getSourceX() {
		return nrSource==0 ? -1 : sourceX/nrSource;
	}
	
	/**
	 * Returns the averaged y-value of the source point
	 * @return The y-Value of the source point
	 */
	public int getSourceY() {
		return nrSource==0 ? -1 : sourceY/nrSource;
	}
	
	/**
	 * Returns the averaged x-value of the target point
	 * @return The x-Value of the target point
	 */
	public int getTargetX() {
		return nrTarget==0 ? -1 : targetX/nrTarget;
	}
	
	/**
	 * Returns the averaged y-value of the target point
	 * @return The y-Value of the target point
	 */
	public int getTargetY() {
		return nrTarget==0 ? -1 : targetY/nrTarget;
	}
	
	/**
	 * Sets a new source x-value. the y value will then be set to the averaged value of the previous inputs.
	 * @param sourceX The x-value of the source point
	 */
	public void setSourceX(int sourceX) {
		this.sourceX = sourceX;
		this.sourceY = getSourceY();
		this.nrSource = 1;
	}
	
	/**
	 * Sets a new source y-value. the x value will then be set to the averaged value of the previous inputs.
	 * @param sourceY The y-value of the source point
	 */
	public void setSourceY(int sourceY) {
		this.sourceY = sourceY;
		this.sourceX = getSourceX();
		this.nrSource = 1;
	}
	
	/**
	 * Sets a new target x-value. the y value will then be set to the averaged value of the previous inputs.
	 * @param targetX The x-value of the target point
	 */
	public void setTargetX(int targetX) {
		this.targetX = targetX;
		this.targetY = getTargetY();
		this.nrTarget = 1;
	}

	/**
	 * Sets a new target y-value. the x value will then be set to the averaged value of the previous inputs.
	 * @param targetY The y-value of the target point
	 */
	public void setTargetY(int targetY) {
		this.targetY = targetY;
		this.targetX = getTargetX();
		this.nrTarget = 1;
	}
	
	/**
	 * Returns the x-value as an average of all used points, the middle point
	 * @return the averaged x-value
	 */
	public int getAvgX() {
		if(nrSource != 0 && nrTarget != 0) {
			int x1 = sourceX/nrSource;
			int x2 = targetX/nrTarget;
			
			return x1+Math.abs(x2-x1)/2;
		}
		return -1;
	}

	/**
	 * Returns the y-value as an average of all used points, the middle point
	 * @return The averaged y-value
	 */
	public int getAvgY() {
		if(nrSource != 0 && nrTarget != 0) {
			int y1 = sourceY/nrSource;
			int y2 = targetY/nrTarget;
			
			return y1+Math.abs(y2-y1)/2;
		}
		return -1;
	}
}
