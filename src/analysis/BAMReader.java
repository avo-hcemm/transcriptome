package analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import data.Gene;
import data.Read;

public class BAMReader {

    public static void processChromosome(File bamFile, String chr, int geneCount,
                                         IntervalTree<Gene> geneTree) throws Exception {
        SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        System.out.println("Chr: " + chr);
		
		/*
		 * SAMFileHeader header = reader.getFileHeader(); for (SAMSequenceRecord seq :
		 * header.getSequenceDictionary().getSequences()) {
		 * System.out.println(seq.getSequenceName()); }
		 */
        
        SAMSequenceRecord seq = reader.getFileHeader().getSequence(chr);
        int chrLength = seq.getSequenceLength();
        
        SAMRecordIterator iter = reader.query(chr, 1, chrLength, false);
        
        List<Read> plusTargetReads = new ArrayList<>();
        List<Read> minusTargetReads = new ArrayList<>();

        int readsCount = 0;
        while(iter.hasNext()) {
        	readsCount ++;
            SAMRecord rec = iter.next();
            if(rec.getReadUnmappedFlag()) continue;

            int readStart = rec.getAlignmentStart();
            int readEnd   = rec.getAlignmentEnd();
            boolean readStrand = rec.getReadNegativeStrandFlag();

            Iterator<IntervalTree.Node<Gene>> overlappingGenes =
                    geneTree.overlappers(readStart, readEnd);

            while(overlappingGenes.hasNext()) {
                Gene gene = overlappingGenes.next().getValue();

                // check strand match & assign target reads to the correct list
                if(readStrand == gene.getStrand()) {
                	if(readStrand) {
                		plusTargetReads.add(initializeRead(rec));
                	}else {
                		minusTargetReads.add(initializeRead(rec));
                	}
                    System.out.println("Read " + rec.getReadName() +
                                       " overlaps gene " + gene.getGeneId());
                }
            }
        }
        iter.close();
        reader.close();
        System.out.println("Read count:" + readsCount );
        
        minusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
        plusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
		/*
		 * //nested loops for (Read minusRead : minusTargetReads) { int minusBegin =
		 * minusRead.getReadStart(); for(Read plusRead: plusTargetReads){ int plusBegin
		 * = plusRead.getReadStart(); int plusLength = plusRead.getReadLength();
		 * if(minusBegin >= plusBegin && minusBegin <= plusLength)
		 * System.out.println("Exact overlap at: " + minusBegin); } }
		 */
		  
		  // advance indeces
		  int i = 0, j = 0;
		  int overlapCount = 0;
		  System.out.println("Overlapping reads:");
		  while (i < minusTargetReads.size() && j < plusTargetReads.size()) {
			  Read minusTargetRead = minusTargetReads.get(i);
			  Read plusTargetRead = plusTargetReads.get(j);
			  
			  int minusBegin = minusTargetRead.getReadStart();
			  int minusEnd = minusBegin + minusTargetRead.getReadLength();
			  int plusBegin = plusTargetRead.getReadStart(); 
			  int plusEnd = plusBegin + plusTargetRead.getReadLength();
			  
			  if(minusBegin < plusEnd && minusEnd > plusBegin) { 
				  int jj = j;
				  while(jj < plusTargetReads.size() && plusTargetReads.get(jj).getReadStart() < minusEnd) {
					  overlapCount ++; //overlap region
					  System.out.println("- strand: " + minusBegin + ", + strand: " + plusBegin);
					  jj++;
				  }
				  i++;
			  }else if(minusEnd < plusBegin){ 
				  i ++; // next read on - strand
			  }else if(minusBegin > plusEnd) {
				  j++; // next read on + strand
			  }
		  }
		  System.out.println("Chromosome: "+chr+ " - Number of protein coding genes: "+geneCount+" - Number of overlaps: "+overlapCount);
		  System.out.println("The program has terminated");
		 
    }
    
    private static Read initializeRead(SAMRecord rec) {
    	int start = rec.getAlignmentStart();
    	int length = rec.getReadLength();
        return new Read(start, length);
    }
    
}