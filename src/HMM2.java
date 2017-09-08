import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;

/**
 * Created by Leonard on 2017-09-08.
 */
public class HMM2 {
    Scanner sc;
    Double[][] A;
    Double[][] B;
    Double[][] pi;
    int[] obs;
    Double[][] delta;
    int[][] deltaIdx;
    int stateLength;


    HMM2(Scanner sc) {
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

        stateLength = pi[0].length;
        deltaIdx = new int[stateLength][obs.length];

        this.forward();
    }

    void forward(){
        delta = new Double[1][stateLength];
        this.initAlpha();

        for (int i = 1; i < obs.length; i++) {
            Double[][] newDelta = multiplicar(delta, A, i);
            for (int j = 0; j < stateLength; j++) {
                delta[0][j] = newDelta[0][j] * B[j][obs[i]];
            }
        }


        this.backTrack();
    }

    void backTrack(){
        int[] result = new int[obs.length];

        Double hi = 0.0;
        int index = 0;
        for (int i = 0; i < delta[0].length; i++) {
            if (delta[0][i] > hi){
                hi = delta[0][i];
                index = i;
            }
        }
        result[obs.length-1] = index;

        for (int i = obs.length-1; i > 0; i--) {
            int val = result[i];
            result[i-1] = deltaIdx[val][i];
        }

        StringBuilder sb = new StringBuilder();
        sb.append(result[0]);
        for (int i = 1; i < result.length; i++) {
            sb.append(" " + result[i]);
        }
        System.out.println(sb.toString());
    }

    void initAlpha(){
        for (int i = 0; i < this.stateLength; i++) {
            delta[0][i] = pi[0][i] * B[i][obs[0]];
        }
    }

    public Double[][] multiplicar(Double[][] A, Double[][] B, int index) {

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
                double hi = 0;
                for (int k = 0; k < aColumns; k++) { // aColumn
                    double prod = A[i][k] * B[k][j];
                    if(prod > hi) {
                        C[i][j] += A[i][k] * B[k][j];
                        //lägg in index
                        deltaIdx[j][index] = k;
                        hi = prod;
                    }
                }
            }
        }

        return C;
    }

    public static void main(String[] args) throws IOException {
//        File file = new File("testCaseHmm2.txt");
//        System.setIn(new FileInputStream(file));

        Scanner sc = new Scanner(System.in).useLocale(Locale.US);
        HMM2 m = new HMM2(sc);
    }
}
