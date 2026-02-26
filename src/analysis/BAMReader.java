package analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import data.Gene;

public class BAMReader {

    public static void processChromosome(File bamFile, String chr,
                                         IntervalTree<Gene> geneTree) throws Exception {
        SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        
		/*
		 * SAMFileHeader header = reader.getFileHeader(); for (SAMSequenceRecord seq :
		 * header.getSequenceDictionary().getSequences()) {
		 * System.out.println(seq.getSequenceName()); }
		 */
        
        SAMSequenceRecord seq = reader.getFileHeader().getSequence(chr);
        int chrLength = seq.getSequenceLength();
        
        SAMRecordIterator iter = reader.query(chr, 1, chrLength, false);
        
        List<SAMRecord> plusReads = new ArrayList<>();
        List<SAMRecord> minusReads = new ArrayList<>();

        while(iter.hasNext()) {
            SAMRecord rec = iter.next();
            if(rec.getReadUnmappedFlag()) continue;

            int readStart = rec.getAlignmentStart();
            int readEnd   = rec.getAlignmentEnd();
            boolean readStrand = rec.getReadNegativeStrandFlag();

            Iterator<IntervalTree.Node<Gene>> overlappingGenes =
                    geneTree.overlappers(readStart, readEnd);

            while(overlappingGenes.hasNext()) {
                Gene gene = overlappingGenes.next().getValue();

                // check strand match
                if(readStrand == gene.getStrand()) {
                	if(readStrand) {
                		minusReads.add(rec);
                	}else {
                		plusReads.add(rec);
                	}
                    System.out.println("Read " + rec.getReadName() +
                                       " overlaps gene " + gene.getGeneId());
                }
            }
        }

        iter.close();
        reader.close();
        
        List<Integer> plusBeginnings = new ArrayList<>();
        List<Integer> minusBeginnings = new ArrayList<>();
        
        for(SAMRecord read: plusReads) {
        	plusBeginnings.add(getReadBeginning(read));
        }
        for(SAMRecord read: minusReads) {
        	minusBeginnings.add(getReadBeginning(read));
        }
        Collections.sort(plusBeginnings);
        int index = 0;
        int window = 10;
        for (int minusBegin : minusBeginnings) {
	        // Find the first plus read that could overlap
	        int firstIdx = Collections.binarySearch(plusBeginnings, minusBegin - window);
	        if (firstIdx < 0) firstIdx = -firstIdx - 1;
	
	        // Find the last plus read that could overlap
	        int lastIdx = Collections.binarySearch(plusBeginnings, minusBegin + window);
	        if (lastIdx < 0) lastIdx = -lastIdx - 2; // -2 because it returns insertion point -1
	
	        // Slice the list
	        List<Integer> candidatePlusReads = plusBeginnings.subList(firstIdx, lastIdx + 1);
	        
	        int idx = Collections.binarySearch(candidatePlusReads, minusBegin);
	        if (idx >= 0) {
	            System.out.println("Exact overlap at: " + minusBegin);
	        }
        }
    }
    
    private static int getReadBeginning(SAMRecord rec) {
        if (rec.getReadNegativeStrandFlag()) {
            return rec.getAlignmentEnd();   // minus strand
        } else {
            return rec.getAlignmentStart(); // plus strand
        }
    }
}