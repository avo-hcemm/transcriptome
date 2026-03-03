package analysis;

import htsjdk.samtools.util.IntervalTree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import data.Gene;

public class GTFIntervalTreeBuilder {

    // Map chromosome -> interval tree
    private final Map<String, IntervalTree<Gene>> chromosomeTrees = new HashMap<>();

    public void parseGTF(File gtfFile) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(gtfFile),16 * 1024 * 1024)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue; // skip header

                int idx1 = line.indexOf('\t');       // end of chr
                String chr = line.substring(0,idx1);
//                if(!chr.equals("21")) continue; // TO BE DELETED
                
                int idx2 = line.indexOf('\t', idx1+1); 
                int idx3 = line.indexOf('\t', idx2+1); 
                int idx4 = line.indexOf('\t', idx3+1); 
                int idx5 = line.indexOf('\t', idx4+1); 
                int idx6 = line.indexOf('\t', idx5+1); 
                int idx7 = line.indexOf('\t', idx6+1); 
                int idx8 = line.indexOf('\t', idx7+1); 
               
                String featureType = line.substring(idx2+1, idx3);
//                if (!featureType.equals("transcript")) continue;// Keep only protein-coding transcripts
                if (!featureType.equals("gene")) continue;
                
                String attributes = line.substring(idx8+1);
//                if (!attributes.contains("transcript_biotype \"protein_coding\"")) continue;
                if (!attributes.contains("gene_biotype \"protein_coding\"")) continue;
                
                
                int start= Integer.parseInt(line.substring(idx3+1, idx4));
                
                
                int end= Integer.parseInt(line.substring(idx4+1, idx5));
                
                
                String strandStr =line.substring(idx6+1, idx7);
                Boolean strand = null;
               
                if (strandStr.equals("+")) {
                    strand = true;
                } else if (strandStr.equals("-")) {
                    strand = false;
                }
                
                // extract gene_id
                String geneId = "unknown";
                
                int geneIdIndex = attributes.indexOf("gene_id \"");
                if (geneIdIndex >= 0) {
                    int startIndex = geneIdIndex + 9;
                    int endIndex = attributes.indexOf("\"", startIndex);
                    geneId = attributes.substring(startIndex, endIndex);
                }
                
                // extract gene_name
                String geneName = "unknown";
 
                int geneNameIndex = attributes.indexOf("gene_name \"");
                if (geneNameIndex >= 0) {
                    int startIndex = geneNameIndex + 11;
                    int endIndex = attributes.indexOf("\"", startIndex);
                    geneName = attributes.substring(startIndex, endIndex);
                }

                // initialize interval tree if absent
                chromosomeTrees.putIfAbsent(chr, new IntervalTree<>());
                IntervalTree<Gene> tree = chromosomeTrees.get(chr);

                // add interval
                tree.put(start, end, new Gene(geneId, geneName, start, end, strand));
            }
        }
        System.out.println("GeneTree parsed");
    }

    public IntervalTree<Gene> getTreeForChromosome(String chr) {
        return chromosomeTrees.get(chr);
    }
    
    public int getProteinCodingGenes(String chr) {
    	IntervalTree<Gene> tree = chromosomeTrees.get(chr);
    	if(tree == null) return 0;
        return chromosomeTrees.get(chr).size();
    }
}