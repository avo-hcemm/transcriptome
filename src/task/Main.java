package task;
import analysis.BAMReader;
import data.Gene;

import htsjdk.samtools.util.IntervalTree;
import java.io.File;

public class Main {
    public static void main(String[] args) throws Exception {

        // example: build geneTree for chr1
        IntervalTree<Gene> geneTree = new IntervalTree<>();
        
        /* * SOD1 Gene (GRCh38 coordinates)
         * Start: 31,659,666
         * End: 31,668,931
         */
        int start = 31659666;
        int end = 31668931;
        
        geneTree.put(start, end, new Gene("SOD1", start, end, false));
        
        // Adding RUNX1 to your gene tree
        start = 34787801;
        end = 35049334;
        geneTree.put(start, end, new Gene("RUNX1", start, end, false));

        String exampleBamFile = args[0];
        String chr = args[1];
        File bamFile = new File(exampleBamFile);

        BAMReader.processChromosome(bamFile, chr, geneTree);
    }
}
