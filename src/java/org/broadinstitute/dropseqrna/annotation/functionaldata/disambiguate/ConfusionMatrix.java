package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;


import java.io.PrintStream;


import java.io.PrintStream;

import java.io.PrintStream;

public class ConfusionMatrix<T extends Enum<T>> {
    private Class<T> enumClass;
    private int[][] matrix;

    public ConfusionMatrix(Class<T> enumClass) {
        if (enumClass == null || !enumClass.isEnum()) {
            throw new IllegalArgumentException("Invalid enum class.");
        }

        this.enumClass = enumClass;
        int numClasses = enumClass.getEnumConstants().length;
        matrix = new int[numClasses][numClasses];
    }

    public void update(T actualClass, T predictedClass) {
        if (actualClass == null || predictedClass == null) {
            throw new IllegalArgumentException("Class values cannot be null.");
        }

        int actualIndex = actualClass.ordinal();
        int predictedIndex = predictedClass.ordinal();
        matrix[actualIndex][predictedIndex]++;
    }

    public void writeFile(PrintStream printStream, char delimiter, String comment) {
        T[] enumConstants = enumClass.getEnumConstants();

        if (comment!=null)
            printStream.println("# "+comment);

        // Print header row
        for (T constant : enumConstants) {
            printStream.printf("%c%s", delimiter, constant.toString());
        }
        printStream.println();

        for (int i = 0; i < matrix.length; i++) {
            // Print row header
            printStream.print(enumConstants[i].toString());

            // Print matrix values
            for (int j = 0; j < matrix[i].length; j++) {
                printStream.printf("%c%d", delimiter, matrix[i][j]);
            }
            printStream.println();
        }
    }

}
