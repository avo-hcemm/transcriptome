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
import data.Pair;
import data.Read;

public class BAMReader {
final static private int MQuality = 10;

    public static int[] processChromosome(File bamFile, String chr,
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
            if (rec.getMappingQuality() < MQuality) { readsLowMQCount ++;continue;}
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
        
        HashSet<Pair<String, String>> pairedGenesSet = new  HashSet<Pair<String, String>>();
		
		 // advance indeces
		 int i = 0, j = 0;
		 int overlappingReadsCount = 0;
		 int uniqueOverlappingReadsCount = 0;
		 System.out.println("Overlapping reads:");
		 while (i < minusTargetReads.size() && j < plusTargetReads.size()) {
			 
			 Read minusTargetRead = minusTargetReads.get(i);
			 int minusBegin = minusTargetRead.getReadStart();
			 int minusEnd = minusBegin + minusTargetRead.getReadLength();

			 Read plusTargetRead = plusTargetReads.get(j);
			 int plusBegin = plusTargetRead.getReadStart(); 
			 int plusEnd = plusBegin + plusTargetRead.getReadLength();
			  
			 if(minusBegin < plusEnd && minusEnd > plusBegin) { 
				 
				 //long table of overlapping read pairs
				 int jj = j;
				 /* Overlapping region scan */
				 while(jj < plusTargetReads.size() && plusTargetReads.get(jj).getReadStart() < minusEnd) {
					 Read mateTargetRead = plusTargetReads.get(jj);
					 int mateBegin = mateTargetRead.getReadStart();
					 int mateEnd = mateBegin + mateTargetRead.getReadLength();
					 int overlapStart = Math.max(minusTargetRead.getReadStart(), mateBegin);
					 int overlapEnd = Math.min(minusTargetRead.getReadStart() + minusTargetRead.getReadLength(), mateEnd);
					 overlappingReadsCount ++; //overlap region
					 System.out.println("- "+minusTargetRead.toString()+", + "+mateTargetRead.toString());
					 //gene pairs
					 Iterator<IntervalTree.Node<Gene>> newReadsOverlappingGenes =
			                    geneTree.overlappers(overlapStart, overlapEnd);
					 
					  initializeNewPairs(newReadsOverlappingGenes, pairedGenesSet);
					 jj++;
				 }
				 /* Overlapping region end */
			 i++;
			 }else if(minusEnd <= plusBegin){ 
				 i ++; // next read on - strand
			 }else{
				 j++; // next read on + strand
			 }
		 }
		 System.out.printf("Paired read-overlapping genes: %n");
		 for(Pair<String,String> pair: pairedGenesSet) {
			 uniqueOverlappingReadsCount ++;
			 System.out.println("- "+pair.getFirst()+", + "+pair.getSecond());
		 }
		 System.out.println("Number of unique reads overlaps: "+ uniqueOverlappingReadsCount);
		 double prop = (double)readsProperlyMappedCount/(readsCount)*100.0;
		 double prop2 = (double)readsUnmappedCount/(readsCount)*100.0;
		 double prop3 = (double)readsSecAlignCount/(readsCount)*100.0;
		 double prop4 = (double)readsLowMQCount/(readsCount)*100.0;
		 int overlapGeneSetSize = overlappingGenesSet.size();
	     System.out.printf("Unmapped reads: %.2f%% %nSecondary/alternative alignments: %.2f%% %nAlignments with low mapping quality: %.2f%%%n",
	    		 prop2,prop3,prop4); 
	     System.out.printf("Total reads: %1d -> Percentage of properly mapped reads: %.2f%% %n", readsCount ,prop);
	     System.out.println("Protein-coding genes overlapping reads :" + overlappingGenesCount+" -> Unique overlapped genes:" +overlapGeneSetSize);
	        
		 return new int[] {overlapGeneSetSize, overlappingReadsCount}; 
    }
    
    private static Read initializeRead(SAMRecord rec, String geneName) {
    	int start = rec.getAlignmentStart();
    	int length = rec.getReadLength();
    	String readName = rec.getReadName();
        return new Read(start, length, readName, geneName);
    }
    
//    private static void initializePairs(Iterator<IntervalTree.Node<Gene>> Genes, HashSet<Pair<String, String>> processedPairs) {
//    	Pair<String, String> pair = new Pair<String, String>();
//    	while(Genes.hasNext()) {
//    		Gene gene = Genes.next().getValue();
//	    	if(!pair.isComplete()){
//	    		if(!gene.getStrand())
//	            	pair.setSecond(gene.getGeneName());// + strand
//	    		 else
//	    			 pair.setFirst(gene.getGeneName());// - strand
//			 }else{
//				 processedPairs.add(pair);
//				 pair = new Pair<String, String>();
//				 if(!gene.getStrand())
//				    	pair.setSecond(gene.getGeneName());// + strand
//				 else
//					 pair.setFirst(gene.getGeneName());// - strand
//			 }
//    	}
//    	 if(pair.isComplete())
//    		 processedPairs.add(pair);
//    }
    
    private static void initializeNewPairs(Iterator<IntervalTree.Node<Gene>> Genes, HashSet<Pair<String, String>> pairSet) {
    	ArrayList<String> minusGenes = new ArrayList<>();
    	ArrayList<String> plusGenes = new ArrayList<>();
    	while(Genes.hasNext()) {
    		Gene gene = Genes.next().getValue();
    		boolean strand = gene.getStrand();
    		String geneName = gene.getGeneName();
    		if(strand)
    			minusGenes.add(geneName);
    		else
    			plusGenes.add(geneName);
    	}
    	
    	if(minusGenes.isEmpty() || plusGenes.isEmpty()) return; 
    	
//    	int minSize = Math.min(minusGenes.size(),plusGenes.size());
//    	int maxSize = Math.max(minusGenes.size(),plusGenes.size());
//    	for(int i=0; i < minSize; i++ ) {
//    		for(int ii=0; ii < minSize; ii++) {
//	    		Pair<String, String> pair = new Pair<String, String>();
//	    		pair.setFirst(minusGenes.get(i));
//	    		pair.setSecond(plusGenes.get(ii));
//	    		if(pair.isComplete()) pairSet.add(pair);
//    		}
//    	}
//    	
//    	for(int j=minSize; j< maxSize; j++) {
//    		for(int jj=0; jj < minSize; jj++) {
//	    		Pair<String, String> pair = new Pair<String, String>();
//	    		if(minusGenes.size()<plusGenes.size()) {
//	    			pair.setFirst(minusGenes.get(jj));
//	    			pair.setSecond(plusGenes.get(j));
//	    		}else {
//	    			pair.setFirst(minusGenes.get(j));
//	    			pair.setSecond(plusGenes.get(jj));
//	    		}
//	    		if(pair.isComplete()) pairSet.add(pair);
//    		}
//    	}
    	for(String minus : minusGenes){
	        for(String plus : plusGenes){
	            Pair<String,String> pair = new Pair<>();
	            pair.setFirst(minus);
	            pair.setSecond(plus);
	            pairSet.add(pair);
	        }
    	}
    }
}