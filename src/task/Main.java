package task;

import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import analysis.BAMReader;
import analysis.GTFIntervalTreeBuilder;
import data.Gene;


public class Main {
    public static void main(String[] args) throws Exception {
    	
    	String gtfFile = args[0];
    	String exampleBamFile = args[1];
//        String chr = args[2];
        
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
        	chrList.addAll(Arrays.asList(args).subList(2, args.length));
        }
        
        File bamFile = new File(exampleBamFile);
        GTFIntervalTreeBuilder builder = new GTFIntervalTreeBuilder();
        builder.parseGTF(new File(gtfFile));

         // Fixed thread pool: 6 threads for my Mac 6-core CPU
        int threads = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(threads);
//        int[] overlapCount = new int[2];
        for (String chr : chrList) {
            IntervalTree<Gene> chrTree = builder.getTreeForChromosome(chr);
            if (chrTree == null) continue; // skip chromosomes not in GTF

            executor.submit(() -> {
                try {
                	int geneCount = builder.getProteinCodingGenes(chr);
                	int[] overlapCount = BAMReader.processChromosome(bamFile, chr, chrTree);
                	System.out.println("Chromosome "+chr+" summary:");
                    double prop = geneCount > 0 ?(double)overlapCount[0]/geneCount*100: 0;
                    System.out.printf("Protein coding genes:%d | Percentage of overlapped genes over protein coding genes:%.2f%%"
                    		+ " | Total number of overlapping reads: %d %n",geneCount,prop,overlapCount[1]);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
        }
        executor.shutdown();
        executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);

        System.out.println("All chromosomes processed.");
        System.out.println("End of the program");
    }
}
