package analysis;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import data.Gene;
import htsjdk.samtools.util.IntervalTree;

public class GeneOverlap {

	
	public static HashMap<String, HashSet<String>> getOverlaps(IntervalTree<Gene> geneTree) {
		HashMap<String, HashSet<String>> overlappingGenesMap = new HashMap<String, HashSet<String>>();
		
		Iterator<IntervalTree.Node<Gene>> iter = geneTree.iterator();
		
		HashSet<String> processedGeneSet = new HashSet<>();
		while(iter.hasNext()) {
			Gene gene = iter.next().getValue();
			int geneStart = gene.getStart();
			int geneEnd = gene.getEnd();
			String geneName = gene.getGeneName();
			
			if(processedGeneSet.contains(geneName))continue;
			
			Iterator<IntervalTree.Node<Gene>> overlappingGenes =
	                 geneTree.overlappers(geneStart, geneEnd);
			
			HashSet<String> geneSet = new HashSet<>();
			 while(overlappingGenes.hasNext()) {
				 Gene ogene = overlappingGenes.next().getValue();
				 String ogeneName = ogene.getGeneName();
				 if (!geneName.equals(ogeneName) 
		                    && gene.getStrand() != ogene.getStrand()) {
					 	processedGeneSet.add(ogeneName);
		                geneSet.add(ogeneName);
		            }
			 }
			 overlappingGenesMap.put(geneName,geneSet);
		}
		return overlappingGenesMap;
		 
	}
	
	public static void printOverlaps(HashMap<String, HashSet<String>> overlappingGenesMap) {
		int count = 0;
		for (HashMap.Entry<String, HashSet<String>> entry : overlappingGenesMap.entrySet()) {
	        String gene = entry.getKey();
	        HashSet<String> overlapSet = entry.getValue();
	        System.out.println(gene + " " + overlapSet);
	        if(overlapSet.size() != 0) count = count + overlapSet.size();
	    }
		System.out.println("Number of overlapping protein-coding genes:"+ count);
	}
	
}
