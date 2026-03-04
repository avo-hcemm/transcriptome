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
		
		while(iter.hasNext()) {
			Gene gene = iter.next().getValue();
			int geneStart = gene.getStart();
			int geneEnd = gene.getEnd();
			String geneName = gene.getGeneName();
			
			Iterator<IntervalTree.Node<Gene>> overlappingGenes =
	                 geneTree.overlappers(geneStart, geneEnd);
			
			HashSet<String> geneSet = new HashSet<>();
			 while(overlappingGenes.hasNext()) {
				 Gene ogene = overlappingGenes.next().getValue();
				 if(gene.getStrand() != ogene.getStrand())
					 geneSet.add(ogene.getGeneName());
			 }
			 overlappingGenesMap.put(geneName,geneSet );
		}
		return overlappingGenesMap;
		 
	}
	
	public static void printOverlaps(HashMap<String, HashSet<String>> overlappingGenesMap) {
		overlappingGenesMap.forEach((g,o) -> System.out.println(g+" "+o));
	}
	
}
