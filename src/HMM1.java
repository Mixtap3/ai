import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Locale;
import java.util.Scanner;

/**
 * Created by Leonard on 2017-09-08.
 */

public class HMM1 {
    Scanner sc;
    Double[][] A;
    Double[][] B;
    Double[][] pi;
    int[] obs;
    Double[][] alpha;
    int piLength;


    HMM1(Scanner sc) {
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

        obs = new int[sc.nextInt()];

        for (int i = 0; i < obs.length; i++) {
            obs[i] = sc.nextInt();
        }
        piLength = pi[0].length;
        this.forward();
    }

    void forward(){
        alpha = new Double[1][piLength];
        this.initAlpha();

        for (int i = 1; i < obs.length; i++) {
            Double[][] newAlpha = multiplicar(alpha, A);
            for (int j = 0; j < piLength; j++) {
                alpha[0][j] = newAlpha[0][j] * B[j][obs[i]];
            }
        }

        Double result = 0.0;
        for(Double d : alpha[0]){
            result += d;
        }
        double r = (double)Math.round(result * 1000000d) / 1000000d;
        System.out.println(r);
    }

    void initAlpha(){
        for (int i = 0; i < this.piLength; i++) {
            alpha[0][i] = pi[0][i] * B[i][obs[0]];
        }
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
//        File file = new File("testHmm1.txt");
//        System.setIn(new FileInputStream(file));

        Scanner sc = new Scanner(System.in).useLocale(Locale.US);
        HMM1 m = new HMM1(sc);
    }
}
