package bioinformatics;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class SubstitutionMatrix {

    /**
     * Substitution matrix class is a 2d array in the form of a map of maps.
     * The map contains keys of substitution matrix row headers (i.e. amino acids/placeholders)
     * paired to maps of substitution matrix column headers and their corresponding values.
     * A map object is used for convenient lookup.
     */
    private Map<String, Map<String, Integer>> substitutionMatrix;

    public SubstitutionMatrix(String filename) throws FileNotFoundException {
        File file = new File(filename);
        initializeMatrixFromFile(file);
    }

    public SubstitutionMatrix(File file) throws FileNotFoundException {
        initializeMatrixFromFile(file);
    }

    /**
     * Reads a formatted substitution matrix file into the substitution matrix
     * map object.
     * @param file the substitution matrix file, File object
     * @throws FileNotFoundException
     */
    private void initializeMatrixFromFile(File file) throws FileNotFoundException {
        Scanner scanner = new Scanner(file);
        substitutionMatrix = new HashMap<String, Map<String, Integer>>();

        String[] keys = scanner.nextLine().trim().split("  ");
        while (scanner.hasNext()) {
            String key = scanner.next().trim();
            Map<String, Integer> tempMap = new HashMap<String, Integer>();
            for (String s : keys) {
                tempMap.put(s, Integer.parseInt(scanner.next()));
            }
            substitutionMatrix.put(key, tempMap);
        }
        scanner.close();
    }

    /**
     * Gets score of two matched keys
     * @param a,b String, keys to match
     * @return integer, score value
     */
    public int getScore(String a, String b) {
        return substitutionMatrix.get(a.toUpperCase()).get(b.toUpperCase());
    }

}
