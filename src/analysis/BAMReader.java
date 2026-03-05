package analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import data.Gene;
import data.Read;

public class BAMReader {

    public static int[] processChromosome(File bamFile, String chr,
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
        int readsUnmappedCount = 0;
        int readsSecAlignCount = 0;
        int readsLowMQCount = 0;
        int readsProperlyMappedCount = 0;
        int overlappingGenesCount = 0;
        Set<Gene> overlappingGenesSet = new HashSet<>();
        
        while(iter.hasNext()) {
        	readsCount ++;
            SAMRecord rec = iter.next();
            if(rec.getReadUnmappedFlag()) {readsUnmappedCount ++;continue;}
            if(rec.isSecondaryOrSupplementary()) {readsSecAlignCount ++;continue;}
            if (rec.getMappingQuality() < 20) { readsLowMQCount ++;continue;}
            readsProperlyMappedCount ++;
            int readStart = rec.getAlignmentStart();
            int readEnd   = rec.getAlignmentEnd();
            boolean readStrand = rec.getReadNegativeStrandFlag();// true for reverse strand, false for forward strand

            Iterator<IntervalTree.Node<Gene>> readsOverlappingGenes =
                    geneTree.overlappers(readStart, readEnd);
            
            while(readsOverlappingGenes.hasNext()) {
            	
                Gene gene = readsOverlappingGenes.next().getValue();

                // check strand match & assign target reads to the correct list
                if(readStrand == gene.getStrand()) {
                	overlappingGenesSet.add(gene);
	            	overlappingGenesCount ++;
	            	if(readStrand) {
	            		minusTargetReads.add(initializeRead(rec, gene.getGeneName()));
	            	}else {
	            		plusTargetReads.add(initializeRead(rec,gene.getGeneName()));
	            	}
	                System.out.println("Read " + rec.getReadName() +
                                       " overlaps gene " + gene.getGeneId() +" : "+ gene.getGeneName());
                }
            }
        }
        iter.close();
        reader.close();
        
        minusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
        plusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
		
		 // advance indeces
		 int i = 0, j = 0;
		 int overlappingReadsCount = 0;
		 System.out.println("Overlapping reads:");
		 while (i < minusTargetReads.size() && j < plusTargetReads.size()) {
			 
			 Read minusTargetRead = minusTargetReads.get(i);
			 int minusBegin = minusTargetRead.getReadStart();
			 int minusEnd = minusBegin + minusTargetRead.getReadLength();

			 Read plusTargetRead = plusTargetReads.get(j);
			 int plusBegin = plusTargetRead.getReadStart(); 
			 int plusEnd = plusBegin + plusTargetRead.getReadLength();
			  
			 if(minusBegin < plusEnd && minusEnd > plusBegin) { 
				 int jj = j;
				 while(jj < plusTargetReads.size() && plusTargetReads.get(jj).getReadStart() < minusEnd) {
					 overlappingReadsCount ++; //overlap region
				 Read targetRead = plusTargetReads.get(jj);
//				 System.out.println("- strand:"+minusTargetReadName+" (" + minusBegin +","+ minusEnd+ "), + strand:"+plusTargetReadName+" (" + plusBegin+","+ plusEnd+ ")");
				 System.out.println("- "+minusTargetRead.toString()+", + "+targetRead.toString());
				 jj++;
			 }
			 i++;
		 }else if(minusEnd < plusBegin){ 
			 i ++; // next read on - strand
		 }else if(minusBegin > plusEnd) {
			 j++; // next read on + strand
			 }
		 }
		 double prop = (double)readsProperlyMappedCount/(readsCount)*100.0;
		 double prop2 = (double)readsUnmappedCount/(readsCount)*100.0;
		 double prop3 = (double)readsSecAlignCount/(readsCount)*100.0;
		 double prop4 = (double)readsLowMQCount/(readsCount)*100.0;
		 int overlapGeneSetSize = overlappingGenesSet.size();
	     System.out.printf("Percentage of unmapped reads: %.2f%% | Percentage of secondary alignments: %.2f%% | Percentage of alignments with low mapping quality: %.2f%%%n",
	    		 prop2,prop3,prop4); 
	     System.out.printf("Total reads: %1d -> Percentage of properly mapped reads: %.2f%% %n", readsCount ,prop);
	     System.out.println("Reads overlapping protein-coding genes:" + overlappingGenesCount+" -> Unique overlapped genes:" +overlapGeneSetSize);
	        
		 return new int[] {overlapGeneSetSize, overlappingReadsCount}; 
    }
    
    private static Read initializeRead(SAMRecord rec, String geneName) {
    	int start = rec.getAlignmentStart();
    	int length = rec.getReadLength();
    	String readName = rec.getReadName();
        return new Read(start, length, readName, geneName);
    }
    
}