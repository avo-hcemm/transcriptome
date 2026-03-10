package task;

import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import analysis.BAMReader;
import analysis.GTFIntervalTreeBuilder;
import data.Gene;


public class Main {
    public static void main(String[] args) throws Exception {
    	
    	String gtfFile = args[0];
    	String exampleBamFile = args[1];
        String chr = args[2];
        
        List<String> chrList = new ArrayList<>();

        // Default list of human chromosomes if no arguments
        List<String> allHumanChromosomes = Arrays.asList(
            "1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20",
            "21","22","X","Y","MT"
        );

        if (args.length <= 2) {
            // Use default chromosomes if no additional arguments
            chrList.addAll(allHumanChromosomes);
        } else {
            // From the 3rd argument onward, add chromosome names to the list
            for (int i = 2; i < args.length; i++) {
                chrList.add(args[i]);
            }
        }
        
    	GTFIntervalTreeBuilder builder = new GTFIntervalTreeBuilder();
        builder.parseGTF(new File(gtfFile));

        // get interval tree for input chromosome
        IntervalTree<Gene> chrTree = builder.getTreeForChromosome(chr);
        
        if(chrTree == null) {
        	System.out.println("Null gene tree. The program will exit");
        	System.exit(0);
        }
        int geneCount = builder.getProteinCodingGenes(chr);
        
        File bamFile = new File(exampleBamFile);
        int[] overlapCount = BAMReader.processChromosome(bamFile, chr, chrTree);
        System.out.println("Chromosome "+chr+" summary:");
        double prop = (double)overlapCount[0]/geneCount*100;
        System.out.printf("Protein coding genes:%d | Percentage of overlapped genes over protein coding genes:%.2f%%"
        		+ " | Total number of overlapping reads: %d %n",+geneCount,prop,overlapCount[1]);
		System.out.println("End of the program");
        /*
         * // Fixed thread pool: 6 threads for your 6-core CPU
        ExecutorService executor = Executors.newFixedThreadPool(6);

        for (String chr : chromosomes) {
            IntervalTree<Gene> tree = builder.getTreeForChromosome(chr);
            if (tree == null) continue; // skip chromosomes not in GTF

            executor.submit(() -> {
                try {
                    BAMReader.processChromosome(bamFile, chr, tree);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
        }

        executor.shutdown();
        while (!executor.isTerminated()) {
            Thread.sleep(1000); // wait for all threads to finish
        }

        System.out.println("All chromosomes processed.");
    }
         * */
    }
}
