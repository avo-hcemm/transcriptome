package data;

import java.util.Objects;

public class Pair<A, B> {
    private A first;
    private B second;

    // 1. Empty Constructor (Allows: new Pair<>())
    public Pair() {}

    // 2. Setters (Allows: pair.setFirst("aa"))
    public void setFirst(A first) { this.first = first; }
    public void setSecond(B second) { this.second = second; }

    // 3. Getters
    public A getFirst() { return first; }
    public B getSecond() { return second; }
    
    // 4. Complete object
    public boolean isComplete() { return (getFirst() != null && getSecond() != null); }

    // 5. ESSENTIAL: equals and hashCode for the "Unique" List (Set)
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Pair<?, ?> pair = (Pair<?, ?>) o;
        return Objects.equals(first, pair.first) && 
               Objects.equals(second, pair.second);
    }

    @Override
    public int hashCode() {
        return Objects.hash(first, second);
    }

    @Override
    public String toString() {
        return "(" + first + ", " + second + ")";
    }
}