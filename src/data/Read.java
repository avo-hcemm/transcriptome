package data;

import java.util.HashSet;

public class Read {
	private final String chromosome;
    private final int readStart;
    private final int readLength;
    private final String readName;
    private final HashSet<String> overlappingGenes = new HashSet<>();

    public Read(String chromosome, int readStart, int readLength, String readName, String geneName) {
    	this.chromosome = chromosome;
    	this.overlappingGenes.add(geneName);
    	this.readName = readName;
        this.readStart = readStart;
        this.readLength = readLength;
    }
    
    public Read(String chromosome, int readStart, int readLength, String readName) {
    	this.chromosome = chromosome;
    	this.readName = readName;
        this.readStart = readStart;
        this.readLength = readLength;
    }

    public int getReadStart() {
        return readStart;
    }

    public int getReadLength() {
        return readLength;
    }

    public String getReadName() {
        return readName;
    }
    
    public HashSet<String> getGeneName() {
        return overlappingGenes;
    }
    
    public void addGene(String geneName) {
    	this.overlappingGenes.add(geneName);
    }
    
    public String getChromosome() {
    	return this.chromosome;
    }
    
    public String toStringShort() {
    	int readEnd = this.readStart + this.readLength;
        return "Read{" +
                "coordinates=(" + this.readStart +
                ", " + readEnd +
                "), name=" + this.readName +
                "), genes=" + this.overlappingGenes +
                '}';
    }
    @Override
    public String toString() {
    	int readEnd = this.readStart + this.readLength;
        return "Read{" +
                "coordinates=(" + this.readStart +
                ", " + readEnd +
                "), name=" + this.readName +
                "), genes=" + this.overlappingGenes +
                '}';
    }
}
