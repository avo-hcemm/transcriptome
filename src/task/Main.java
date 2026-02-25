package task;
import analysis.BAMReader;
import data.Gene;

import htsjdk.samtools.util.IntervalTree;
import java.io.File;

public class Main {
    public static void main(String[] args) throws Exception {

        // example: build geneTree for chr1
        IntervalTree<Gene> geneTree = new IntervalTree<>();
        geneTree.put(1000, 5000, new Gene("BRCA1", 1000, 5000, false));

        File bamFile = new File("example.bam");

        BAMReader.processChromosome(bamFile, "chr1", geneTree);
    }
}
