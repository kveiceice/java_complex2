import java.util.Scanner;

class Complex{
    public final double real;
    public final double imaginary;

    public Complex(double re, double im){
        this.real = (re == -0) ? 0 : re;
        this.imaginary = (im == -0) ? 0 : im;
    }

    public Complex add(Complex other){
        Complex new_c = new Complex(this.real + other.real, this.imaginary + other.imaginary);
        return new_c;
    }

    public Complex subtract(Complex other){
        Complex new_c = new Complex(this.real - other.real, this.imaginary - other.imaginary);
        return new_c;
    }

    public Complex multiply(Complex other){
        // (a+bi)(c+di) = (ac-bd) + (ad+bc)i
        Complex new_c = new Complex((this.real * other.real) - (this.imaginary * other.imaginary), (this.real * other.imaginary) + (this.imaginary * other.real));
        return new_c;
    }

    public Complex divide(Complex other){
        //A+BI
        //C+DI
        double denom = (other.real)*(other.real) + (other.imaginary)*(other.imaginary);
        if (denom == 0){
            throw new IllegalArgumentException("Второе комплексное число равно 0!");
        }
        Complex new_c = new Complex((this.real*other.real + this.imaginary*other.imaginary)/denom, (this.imaginary * other.real - this.real * other.imaginary)/denom);
        return new_c;
    }

    public double abs(){
        double abs_c = Math.sqrt(real*real+imaginary*imaginary);
        return abs_c;
    }

    public String ToString(){
        if (imaginary >= 0){
            String new_c = String.format("%.2f + %.2fi", real, imaginary);
            return new_c;
        }
        String new_c = String.format("%.2f - %.2fi", real, -imaginary);
        return new_c;
    }
}

class ComplexMatrix{
    private final int rows;
    private final int cols;
    private final Complex[][] matrix;

    public ComplexMatrix(int rows, int cols){
        this.rows = rows;
        this.cols = cols;
        this.matrix = new Complex[rows][cols];
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                matrix[i][j] = new Complex(0, 0);
            }
        }
    }

    public void put(int i, int j, Complex value){
        matrix[i][j] = value;
        return;
    }

    public ComplexMatrix add(ComplexMatrix other){
        if(rows != other.rows || cols != other.cols){
            throw new IllegalArgumentException("Матрицы нельзя сложить(Не одинаковое кол-во столбцов или строк.)");
        }
        ComplexMatrix res = new ComplexMatrix(rows, cols);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; ++j){
               res.matrix[i][j] = matrix[i][j].add(other.matrix[i][j]);
            }
        }
        return res;
    }

    public ComplexMatrix subtract(ComplexMatrix other){
        if(rows != other.rows || cols != other.cols){
            throw new IllegalArgumentException("Матрицы нельзя вычесть(Не одинаковое кол-во столбцов или строк.)");
        }
        ComplexMatrix res = new ComplexMatrix(rows, cols);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; ++j){
                res.matrix[i][j] = matrix[i][j].subtract(other.matrix[i][j]);
            }
        }
        return res;
    }

    public ComplexMatrix multiply(ComplexMatrix other){
        if(cols != other.rows){
            throw new IllegalArgumentException("Матрицы нельзя умножить)");
        }
        ComplexMatrix res = new ComplexMatrix(this.rows, other.cols);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < other.cols; ++j){
                Complex sum = new Complex(0, 0);
                for(int k = 0; k < cols; ++k){

                    sum = sum.add(matrix[i][k].multiply(other.matrix[k][j]));


                }
                res.put(i, j, sum);
            }
        }
        return res;
    }

    public ComplexMatrix Transponse(){
        ComplexMatrix t_matrix = new ComplexMatrix(this.cols, this.rows);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; ++j){
                t_matrix.put(j, i, matrix[i][j]);
            }
        }
        return t_matrix;
    }

    public Complex det(){
        if(rows != cols){
            throw new IllegalArgumentException("Матрица не квадратичная, нельзя посчитать определитель");
        }

        if (rows == 1) {
            return matrix[0][0];
        }

        if (rows == 2) {
            return matrix[0][0].multiply(matrix[1][1]).subtract(matrix[0][1].multiply(matrix[1][0]));
        }

        Complex determinant = new Complex(0, 0);

        for (int column = 0; column < cols; column++) {
            ComplexMatrix minorMatrix = getMinorMatrix(0, column);

            Complex minorDeterminant = minorMatrix.det();

            Complex sign;
            if (column % 2 == 0) {
                sign = new Complex(1, 0);  // +1
            } else {
                sign = new Complex(-1, 0); // -1
            }

            Complex element = matrix[0][column];

            Complex term = sign.multiply(element).multiply(minorDeterminant);
            determinant = determinant.add(term);
        }



        return determinant;
    }

    private ComplexMatrix getMinorMatrix(int skipRow, int skipCol) {

        int minorRows = rows - 1;
        int minorCols = cols - 1;

        ComplexMatrix minor = new ComplexMatrix(minorRows, minorCols);

        int minorRowIndex = 0;

        for (int row = 0; row < rows; row++) {
            if (row == skipRow) {
                continue;
            }

            int minorColIndex = 0;

            for (int col = 0; col < cols; col++) {
                if (col == skipCol) {
                    continue;
                }

                minor.matrix[minorRowIndex][minorColIndex] = matrix[row][col];

                minorColIndex++;
            }

            minorRowIndex++;
        }
        return minor;
    }

    public ComplexMatrix inverse(){
        if (rows != cols)
            throw new IllegalArgumentException("Обратная только для квадратных");
        Complex det = this.det();
        if (det.abs() < 1e-10)
            throw new ArithmeticException("Матрица вырожденная, обратной не существует");

        ComplexMatrix cofactors = new ComplexMatrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                ComplexMatrix minor = new ComplexMatrix(rows - 1, cols - 1);
                int mi = 0;
                for (int r = 0; r < rows; r++) {
                    if (r == i) continue;
                    int mj = 0;
                    for (int c = 0; c < cols; c++) {
                        if (c == j) continue;
                        minor.matrix[mi][mj] = matrix[r][c];
                        mj++;
                    }
                    mi++;
                }
                Complex sign = ((i + j) % 2 == 0) ? new Complex(1, 0) : new Complex(-1, 0);
                cofactors.matrix[i][j] = sign.multiply(minor.det());
            }
        }
        ComplexMatrix adjugate = cofactors.Transponse();
        ComplexMatrix result = new ComplexMatrix(rows, cols);
        Complex detInv = new Complex(1, 0).divide(det);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                result.matrix[i][j] = adjugate.matrix[i][j].multiply(detInv);
        return result;
    }

    public ComplexMatrix divide(ComplexMatrix other){
        ComplexMatrix invers = other.inverse();
        ComplexMatrix res = this.multiply(invers);
        return res;
    }

    public void print(){
        for(int k = 0; k < (21*cols + 2 * cols); k++){
            System.out.print("-");
        }
        System.out.println();
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                System.out.printf("|   %15s   |", matrix[i][j].ToString());
            }
            System.out.println();
            for(int k = 0; k < (21*cols + 2 * cols); k++){
                System.out.print("-");
            }
            System.out.println();
        }
    }
}

public class Main{
    public static ComplexMatrix inputMatrix(Scanner sc){
        System.out.print("Введите количество строк: ");
        int rows = sc.nextInt();
        System.out.print("Введите количество столбцов: ");
        int coll = sc.nextInt();
        ComplexMatrix Matrix = new ComplexMatrix(rows, coll);
        System.out.printf("Вводите комплексные числа через пробел, всего чисел ввести - %d\n", coll*rows);
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < coll; j++){
                double real = sc.nextDouble();
                double imaginary = sc.nextDouble();
                Matrix.put(i, j, new Complex(real, imaginary));

            }
        }
        return Matrix;
    }

    private static ComplexMatrix get_matrix_from_name(String name, ComplexMatrix A, ComplexMatrix B, ComplexMatrix C){
        switch(name.toLowerCase()){
            case "a": return A;
            case "b": return B;
            case "c": return C;
            default: return null;
        }
    }


    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        ComplexMatrix A = null;
        ComplexMatrix B = null;
        ComplexMatrix C = null;
        ComplexMatrix option1 = null;
        ComplexMatrix option2 = null;
        ComplexMatrix res = null;

        System.out.println("----СКРИПТ ДЛЯ РАБОТЫ С МАТРИЦАМИ КОМПЛЕКСНЫХ ЧИСЕЛ----");
        System.out.println("Матрица C является временной хранящей последнее существенное действие");
        System.out.println("Любой результат последующих операций(кроме вывода и подсчета определителя и задания A/B) изменит матрицу C");
        System.out.println("Комплексное число вида 2+2i записывать в виде 2 2");
        System.out.println("Буквы вводить на английском");

        System.out.println();
        while (true){
            System.out.println("1 - Задать/Изменить матрицу A");
            System.out.println("2 - Задать/Изменить матрицу B");
            System.out.println("3 - Сложить");
            System.out.println("4 - Вычесть");
            System.out.println("5 - Умножить");
            System.out.println("6 - Поделить");
            System.out.println("7 - Транспонировать");
            System.out.println("8 - Определитель");
            System.out.println("9 - Обратная");
            System.out.println("10 - Вывести");
            System.out.println("-1 - завершить работу");
            System.out.println("Введите действие: ");
            int option = sc.nextInt();

            try{
                switch(option){
                    case 1: A = inputMatrix(sc); break;
                    case 2: B = inputMatrix(sc); break;

                    case 3:
                        System.out.print("Введите первую матрицу для сложения(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        System.out.print("Введите вторую матрицу для сложения(A/B/C): ");
                        option2 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null || option2 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицы(или же была введена несуществующая матрица)");
                        }
                        res = option1.add(option2);
                        C = res;
                        res.print();
                        break;
                    case 4:
                        System.out.print("Введите первую матрицу для вычитания(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        System.out.print("Введите вторую матрицу для вычитания(A/B/C): ");
                        option2 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null || option2 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицы(или же была введена несуществующая матрица)");
                        }
                        res = option1.subtract(option2);
                        C = res;
                        res.print();
                        break;
                    case 5:
                        System.out.print("Введите первую матрицу для умножения(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        System.out.print("Введите вторую матрицу для умножения(A/B/C): ");
                        option2 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null || option2 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицы(или же была введена несуществующая матрица)");
                        }
                        res = option1.multiply(option2);
                        C = res;
                        res.print();
                        break;
                    case 6:
                        System.out.print("Введите первую матрицу для деления(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        System.out.print("Введите вторую матрицу для деления(A/B/C): ");
                        option2 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null || option2 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицы(или же была введена несуществующая матрица)");
                        }
                        res = option1.divide(option2);
                        C = res;
                        res.print();
                        break;
                    case 7:
                        System.out.print("Введите матрицу для транспонирования(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицу(или была введена несуществующая матрица)");
                        }
                        res = option1.Transponse();
                        C = res;
                        res.print();
                        break;

                    case 8:
                        System.out.print("Введите матрицу для вычисления определителя(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицу(или была введена несуществующая матрица)");
                        }
                        String matr_name;
                        if(option1 == A){
                            matr_name = "A";
                        }
                        else if(option1 == B){
                            matr_name = "B";
                        }
                        else if(option1 == C){
                            matr_name = "C";
                        }
                        else{
                            matr_name = "Неизвестная матрица";
                        }
                        Complex res_c = option1.det();
                        System.out.print("Детерминант матрицы " + matr_name + " = ");
                        System.out.println(res_c.ToString());
                        break;
                    case 9:
                        System.out.print("Введите матрицу для подсчета обратной(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицу(или была введена несуществующая матрица)");
                        }
                        res = option1.inverse();
                        C = res;
                        res.print();
                        break;
                    case 10:
                        System.out.print("Введите матрицу для вывода(A/B/C): ");
                        sc.nextLine();
                        option1 = get_matrix_from_name(sc.nextLine(), A, B, C);
                        if(option1 == null){
                            throw new IllegalArgumentException("Сначала создайте матрицу(или была введена несуществующая матрица)");
                        }
                        option1.print();
                        break;


                    case -1:
                        return;
                    default: System.out.println("Нет такой опции");
                }
            } catch(Exception exc){
                System.out.println("Произошла ошибка - " + exc.getMessage());
            }
        }
    }
}