/**
 * 
 */
package org.reactome.panther;

/**
 * An Exception that signals that the pathway that is intended to be converted was not found
 * @author andreash
 *
 */
public class PathwayNotFoundException extends Exception {
	private static final long serialVersionUID = 1L;

	/**
	 * Constructs an <code>PathwayNotFoundException</code>.
	 */
	public PathwayNotFoundException() {
		super("Pathway not found");
	}

	/**
	 * Constructs an <code>PathwayNotFoundException</code> with the 
    * specified detailed message.
	 * @param message the detailed message
	 */
	public PathwayNotFoundException(String message) {
		super("Pathway not found - " + message);
	}

	 /**
    * Constructs an <code>PathwayNotFoundException</code> with the 
    * specified cause. 
    * @param cause the cause.
    */
	public PathwayNotFoundException(Throwable cause) {
		super(cause);
	}

	/**
	 * Constructs an <code>PathwayNotFoundException</code> with the 
    * specified detailed message and cause.
	 * @param message the detailed message
	 * @param cause the cause
	 */
	public PathwayNotFoundException(String message, Throwable cause) {
		super("Pathway not found - " + message, cause);
	}
}
