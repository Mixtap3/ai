/**
 * Created by Leonard on 2017-09-07.
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Scanner;
import java.util.Locale;

public class Main {

    Scanner sc;
    Double[][] A;
    Double[][] B;
    Double[][] pi;


    Main(Scanner sc) {
        this.sc = sc;
        this.createMatrix();
    }

    void createMatrix() {

        int arow = sc.nextInt();
        int acol = sc.nextInt();
        this.A = new Double[arow][acol];
        for (int i = 0; i < arow; i++) {
            for (int j = 0; j < acol; j++) {
                this.A[i][j] = sc.nextDouble();
            }
        }

        int brow = sc.nextInt();
        int bcol = sc.nextInt();
        this.B = new Double[brow][bcol];
        for (int i = 0; i < brow; i++) {
            for (int j = 0; j < bcol; j++) {
                this.B[i][j] = sc.nextDouble();
            }
        }

        int pirow = sc.nextInt();
        int picol = sc.nextInt();
        this.pi = new Double[pirow][picol];
        for (int i = 0; i < pirow; i++) {
            for (int j = 0; j < picol; j++) {
                this.pi[i][j] = sc.nextDouble();
            }
        }

        this.multiply();
    }

    void multiply() {
        Double[][] newpi = multiplicar(this.pi, this.A);
        Double[][] observations = multiplicar(newpi, this.B);

        this.writeOut(observations);
    }

    void writeOut(Double[][] vec) {
        int cols = vec[0].length;

        StringBuilder sb = new StringBuilder();
        sb.append("1 " + cols);
        for (int i = 0; i < cols; i++) {
            double d = (double)Math.round(vec[0][i] * 100d) / 100d;
            sb.append(" " + d);
        }
        System.out.println(sb.toString());
    }

    public static Double[][] multiplicar(Double[][] A, Double[][] B) {

        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;

        if (aColumns != bRows) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        Double[][] C = new Double[aRows][bColumns];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                C[i][j] = 0.00000;
            }
        }

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }

    public static void main(String[] args) throws IOException {
//        File file = new File("testCase.txt");
//        System.setIn(new FileInputStream(file));

        Scanner sc = new Scanner(System.in).useLocale(Locale.US);
        Main m = new Main(sc);
    }
}
