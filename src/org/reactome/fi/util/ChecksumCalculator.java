/*
 * Created on Mar 26, 2009
 *
 */
package org.reactome.fi.util;

import jonelo.jacksum.JacksumAPI;
import jonelo.jacksum.algorithm.AbstractChecksum;


/**
 * Basically a warpper to AbstractChecksum to make client code simpler.
 * @author wgm
 *
 */
public class ChecksumCalculator {
    // To calculate sequence check sum
    private AbstractChecksum checksum;
    
    public ChecksumCalculator() {
        try {
            checksum = JacksumAPI.getChecksumInstance("crc64");
            checksum.setEncoding(AbstractChecksum.HEX_UPPERCASE);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public String calculateChecksum(String sequence) {
        checksum.reset();
        checksum.update(sequence.getBytes());
        return checksum.getFormattedValue();
    }
}
