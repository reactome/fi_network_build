/*
 * Created on Sep 28, 2006
 *
 */
package org.reactome.funcInt;

/**
 * The types that indicate what classes are used to extract functional interctions. There
 * are only three types avaiable: COMPLEX, REACTION and INTERACTION (classes defined in 
 * Reactome model).
 * @author guanming
 *
 */
public enum ReactomeSourceType {
   COMPLEX,
   REACTION,
   INTERACTION,
   TARGETED_INTERACTION
}
