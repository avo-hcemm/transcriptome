package analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
            
            int genesCount = 0;
            Read read = null;
            
            while(readsOverlappingGenes.hasNext()) {
           
                Gene gene = readsOverlappingGenes.next().getValue();
                genesCount ++;

                // check strand match & assign target reads to the correct list
                if(readStrand == gene.getStrand()) {
                	overlappingGenesSet.add(gene);
	            	overlappingGenesCount ++;
	            	if(read == null ) {
	            		read = initializeRead(chr,rec, gene.getGeneName());
		            	if(readStrand) minusTargetReads.add(read);
		            	else plusTargetReads.add(read);			            		
	            	}else {    		
		            	if(readStrand) minusTargetReads.add(updateRead(read,gene.getGeneName()));
		            	else plusTargetReads.add(updateRead(read,gene.getGeneName()));	
	            		}
	            }
//	                System.out.println("Read " + rec.getReadName() +
//                                     " overlaps gene " + gene.getGeneId() +" : "+ gene.getGeneName());
            }
        }
        iter.close();
        reader.close();
        
        minusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
        plusTargetReads.sort(Comparator.comparingInt(Read::getReadStart));
        
        HashSet<Pair<String, String>> pairedGenesSet = new  HashSet<Pair<String, String>>();
        boolean pairedGenesFound = false;
		
		 // advance indeces
		 int i = 0, j = 0;
		 int overlappingReadsCount = 0;
		 int uniqueOverlappingReadsCount = 0;
		 System.out.println("Overlapping reads: [OverlapStart, OverlapEnd], Overlap type, Strand, Read{Start-End, Name, Genes}");
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
					 int overlapStart = Math.max(minusBegin, mateBegin);
					 int overlapEnd = Math.min(minusEnd, mateEnd);
					 String type = null;
					 if(overlapStart == mateBegin) {
						 if(overlapEnd == minusEnd)
							 type = "5'";
						 else{
							 type = "N";
						 }
					 }else if(overlapStart == minusBegin) {
						 if(overlapEnd == mateEnd)
							 type = "3'";
						 else{
							 type = "N";
						 }
					 }
					 overlappingReadsCount ++; //overlap region
					
					 System.out.println("[ "+overlapStart+", "+overlapEnd+"], OverlapType:"+type+" Strand:- "+minusTargetRead.toString()+", Strand:+ "+mateTargetRead.toString());
					 //validating the actual gene pairs
					 Iterator<IntervalTree.Node<Gene>> newReadsOverlappingGenes =
			                    geneTree.overlappers(overlapStart, overlapEnd);
					 
					 pairedGenesFound = initializeNewPairs(newReadsOverlappingGenes, pairedGenesSet);
//					 if(pairedGenesFound) printNewPairs(overlapStart, overlapEnd, pairedGenesSet);
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
		 System.out.printf("Read-covered overlapping genes: %n");
		 for(Pair<String,String> pair: pairedGenesSet) {
			 uniqueOverlappingReadsCount ++;
			 System.out.println("- "+pair.getFirst()+", + "+pair.getSecond());
		 }
		 System.out.println("Number of unique reads read-overlapping genes: "+ uniqueOverlappingReadsCount);
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
    
    private static Read initializeRead(String chr, SAMRecord rec, String geneName) {
    	int start = rec.getAlignmentStart();
    	int length = rec.getReadLength();
    	String readName = rec.getReadName();
        return new Read(chr, start, length, readName, geneName);
    }
    
    private static Read updateRead(Read read, String geneName) {
    	read.addGene(geneName);
    	return  read;
    }
    
    private static boolean initializeNewPairs(Iterator<IntervalTree.Node<Gene>> Genes, HashSet<Pair<String, String>> overlapGenePairs) {
    	ArrayList<Pair<String, Integer>> minusGenes = new ArrayList<Pair<String, Integer>>();
    	ArrayList<Pair<String, Integer>> plusGenes = new ArrayList<Pair<String, Integer>>();
    	while(Genes.hasNext()) {
    		Gene gene = Genes.next().getValue();
    		boolean strand = gene.getStrand();
    		String geneName = gene.getGeneName();
    		int geneStart = gene.getStart();
    		Pair<String, Integer> pair = new Pair<String, Integer>();
    		pair.setFirst(geneName);
    		pair.setSecond(geneStart);
    		if(strand)
    			minusGenes.add(pair);
    		else
    			plusGenes.add(pair);
    	}
    	
    	if(minusGenes.isEmpty() || plusGenes.isEmpty()) return false; 
    	

    	for(Pair<String, Integer> minusPair : minusGenes){
    		String minusName = minusPair.getFirst();
	        for(Pair<String, Integer> plusPair : plusGenes){
	        	String plusName = plusPair.getFirst();
	            Pair<String,String> pair = new Pair<>();
	            pair.setFirst(minusName);
	            pair.setSecond(plusName);
	            overlapGenePairs.add(pair);
	        }
    	}
    	return true;
    }
    private static void printNewPairs(int intervalStart, int intervalEnd, HashSet<Pair<String, String>> pairSet) {
    	System.out.println("[ "+intervalStart+", "+intervalEnd+"]");
    	for(Pair<String, String> p : pairSet) {
    		System.out.println("- "+p.getFirst()+", + "+p.getSecond());
    	}
    	
    }
}