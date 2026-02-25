package analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalTree;
import java.io.File;
import java.util.Iterator;

public class BAMReader {

    public static void processChromosome(File bamFile, String chr,
                                         IntervalTree<Gene> geneTree) throws Exception {
        SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        SAMRecordIterator iter = reader.query(chr, 0, 0, false);

        while(iter.hasNext()) {
            SAMRecord rec = iter.next();
            if(rec.getReadUnmappedFlag()) continue;

            int readStart = rec.getAlignmentStart();
            int readEnd   = rec.getAlignmentEnd();
            boolean readStrand = rec.getReadNegativeStrandFlag();

            Iterator<IntervalTree.Node<Gene>> overlaps =
                    geneTree.overlappers(readStart, readEnd);

            while(overlaps.hasNext()) {
                Gene gene = overlaps.next().getValue();

                // check strand match
                if(readStrand == gene.getStrand()) {
                    System.out.println("Read " + rec.getReadName() +
                                       " overlaps gene " + gene.getGeneId());
                }
            }
        }

        iter.close();
        reader.close();
    }
}