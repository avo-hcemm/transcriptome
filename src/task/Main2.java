package task;

import htsjdk.samtools.util.IntervalTree;
import java.io.File; 

import analysis.GeneOverlap;
import analysis.GTFIntervalTreeBuilder;
import data.Gene;


public class Main2 {
    public static void main(String[] args) throws Exception {
    	
    	String gtfFile = args[0];
        String chr = args[1];
        
       // List of chromosomes to process
        String[] chromosomes = {
            "1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20","21","22","X","Y"
        };

        
    	GTFIntervalTreeBuilder builder = new GTFIntervalTreeBuilder();
        builder.parseGTF(new File(gtfFile));

        // get interval tree for input chromosome 
        IntervalTree<Gene> chrTree = builder.getTreeForChromosome(chr);
        
        if(chrTree == null) {
        	System.out.println("Null gene tree. The program will exit");
        	System.exit(0);
        }
        
        GeneOverlap.printOverlaps(GeneOverlap.getOverlaps(chrTree));
        
    }
}
