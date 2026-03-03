package data;

public class Read {
    private final int readStart;
    private final int readLength;
    private final String readName;
    private final String geneName;

    public Read(int readStart, int readLength, String readName, String geneName) {
    	this.geneName = geneName;
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
    
    public String getGeneName() {
        return geneName;
    }
    @Override
    public String toString() {
    	int readEnd = readStart + readLength;
        return "Read{" +
                "coordinates=(" + readStart +
                ", " + readEnd +
                "), name=" + readName +
                "), gene=" + geneName +
                '}';
    }
}
