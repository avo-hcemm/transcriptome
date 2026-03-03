package data;

public class Gene {
    private String geneId;
    private String geneName;
    private int start;
    private int end;
    private boolean strand; // true for - strand, false for + strand

    public Gene(String geneId, String geneName, int start, int end, boolean strand) {
    	this.geneName = geneName;
        this.geneId = geneId;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    // getters
    public String getGeneId() { return geneId; }
    public String getGeneName() { return geneName; }
    public int getStart() { return start; }
    public int getEnd() { return end; }
    public boolean getStrand() { return strand; }
}