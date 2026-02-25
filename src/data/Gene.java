package data;

public class Gene {
    private String geneId;
    private int start;
    private int end;
    private boolean strand; // false = +, true = -

    public Gene(String geneId, int start, int end, boolean strand) {
        this.geneId = geneId;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    // getters
    public String getGeneId() { return geneId; }
    public int getStart() { return start; }
    public int getEnd() { return end; }
    public boolean getStrand() { return strand; }
}