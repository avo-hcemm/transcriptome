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
        int readsProperlyMappedCount = 0;
        int overlappingGenesCount = 0;
        Set<Gene> overlappingGenesSet = new HashSet<>();
        
        while(iter.hasNext()) {
        	readsCount ++;
            SAMRecord rec = iter.next();
            if(rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()) continue;
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
        System.out.println("Read count:" + readsCount + " -> Properly mapped reads count:" + readsProperlyMappedCount );
        System.out.println("Mapped genes count:" + overlappingGenesCount+" -> Overlapping genes count:" +overlappingGenesSet.size());
        
        minusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
        plusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
		
		  
		 // advance indeces
		 int i = 0, j = 0;
		 int overlapCount = 0;
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
					 overlapCount ++; //overlap region
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
		 System.out.println("Chromosome: "+chr+ " - Number of protein coding genes: "+geneCount+" - Number of overlaps: "+overlapCount);
		 System.out.println("The program has terminated"); 
    }
    
    private static Read initializeRead(SAMRecord rec, String geneName) {
    	int start = rec.getAlignmentStart();
    	int length = rec.getReadLength();
    	String readName = rec.getReadName();
        return new Read(start, length, readName, geneName);
    }
    
}