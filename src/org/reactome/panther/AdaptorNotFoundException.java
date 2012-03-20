/**
 * 
 */
package org.reactome.panther;

/**
 * An Exception that should tell the user that an Adaptor has not been set or cannot be used.
 * @author andreash
 *
 */
public class AdaptorNotFoundException extends Exception {
	private static final long serialVersionUID = 1L;

	/**
	 * Constructs an <code>AdaptorNotFoundException</code>.
	 */
	public AdaptorNotFoundException() {
		super("Adaptor not found");
	}

	/**
	 * Constructs an <code>AdaptorNotFoundException</code> with the 
    * specified detailed message.
	 * @param message the detailed message
	 */
	public AdaptorNotFoundException(String message) {
		super(message);
	}

	/**
    * Constructs an <code>AdaptorNotFoundException</code> with the 
    * specified cause. 
    * @param cause the cause.
    */
	public AdaptorNotFoundException(Throwable cause) {
		super(cause);
	}

	/**
	 * Constructs an <code>AdaptorNotFoundException</code> with the 
    * specified detailed message and cause.
	 * @param message the detailed message
	 * @param cause the cause
	 */
	public AdaptorNotFoundException(String message, Throwable cause) {
		super(message, cause);
	}
}
