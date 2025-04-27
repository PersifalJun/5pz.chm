package org.entity;

public class SVD {

    static final double EPS = 1e-20;

    public static class SVDResult {
        public double[][] U;
        public double[][] Sigma;
        public double[][] V;
    }

    // Основной метод: SVD-разложение матрицы A
    public static SVDResult decompose(double[][] A) {
        int N = A.length;
        double[][] U = identity(N); // Матрица U — начальная единичная
        double[][] Sigma = copy(A); // Сигма — копия A
        double[][] V = identity(N); // Матрица V — начальная единичная

        // Этап I: бидиагонализация
        for (int i = 0; i < N - 1; i++) {
            columnTransformation(Sigma, U, i, i);   // Обнуление ниже диагонали в столбце
            rowTransformation(Sigma, V, i, i + 1);  // Обнуление справа от диагонали в строке
        }

        // Этап II: преследование диагонали
        double[] up = new double[N - 1];
        double[] down = new double[N - 1];
        int countUpElements;

        do {
            countUpElements = 0;
            for (int i = 0; i < N - 1; i++) {
                if (Math.abs(Sigma[i][i + 1] - up[i]) > EPS) {
                    up[i] = Sigma[i][i + 1];
                    deleteUpperTriangle(Sigma, V, i, i + 1);
                } else {
                    countUpElements++;
                }
            }

            for (int i = 0; i < N - 1; i++) {
                if (Math.abs(Sigma[i + 1][i] - down[i]) > EPS) {
                    down[i] = Sigma[i + 1][i];
                    deleteLowerTriangle(Sigma, U, i + 1, i);
                }
            }
        } while (countUpElements != N - 1);

        // Убираем отрицательные значения на диагонали
        absSingularValues(Sigma, U);

        // Сортируем сингуляры по убыванию
        sortSingularValues(Sigma, U, V);

        // Возвращаем результат
        SVDResult result = new SVDResult();
        result.U = U;
        result.Sigma = Sigma;
        result.V = V;
        return result;
    }

    // Создание единичной матрицы
    private static double[][] identity(int N) {
        double[][] id = new double[N][N];
        for (int i = 0; i < N; i++) id[i][i] = 1.0;
        return id;
    }

    // Копирование матрицы
    private static double[][] copy(double[][] A) {
        int N = A.length;
        double[][] copy = new double[N][N];
        for (int i = 0; i < N; i++)
            System.arraycopy(A[i], 0, copy[i], 0, N);
        return copy;
    }

    // Преобразование Хаусхолдера для столбца
    private static void columnTransformation(double[][] A, double[][] U, int i, int j) {
        int N = A.length;
        double[] p = new double[N];
        double s = 0;
        for (int I = j; I < N; I++) s += A[I][i] * A[I][i];

        if (Math.sqrt(Math.abs(s - A[j][i] * A[j][i])) > EPS) {
            double beta = A[j][i] < 0 ? Math.sqrt(s) : -Math.sqrt(s);
            double mu = 1.0 / (beta * (beta - A[j][i]));

            for (int I = 0; I < N; I++) {
                if (I >= j) p[I] = A[I][i];
            }
            p[j] -= beta;

            // Обновляем A
            for (int m = 0; m < N; m++) {
                double temp = 0;
                for (int n = j; n < N; n++) temp += A[n][m] * p[n];
                temp *= mu;
                for (int n = j; n < N; n++) A[n][m] -= temp * p[n];
            }

            // Обновляем U
            for (int m = 0; m < N; m++) {
                double temp = 0;
                for (int n = j; n < N; n++) temp += U[m][n] * p[n];
                temp *= mu;
                for (int n = j; n < N; n++) U[m][n] -= temp * p[n];
            }
        }
    }

    // Преобразование Хаусхолдера для строки
    private static void rowTransformation(double[][] A, double[][] V, int i, int j) {
        int N = A.length;
        double[] p = new double[N];
        double s = 0;
        for (int I = j; I < N; I++) s += A[i][I] * A[i][I];

        if (Math.sqrt(Math.abs(s - A[i][j] * A[i][j])) > EPS) {
            double beta = A[i][j] < 0 ? Math.sqrt(s) : -Math.sqrt(s);
            double mu = 1.0 / (beta * (beta - A[i][j]));

            for (int I = 0; I < N; I++) {
                if (I >= j) p[I] = A[i][I];
            }
            p[j] -= beta;

            // Обновляем A
            for (int m = 0; m < N; m++) {
                double temp = 0;
                for (int n = j; n < N; n++) temp += A[m][n] * p[n];
                temp *= mu;
                for (int n = j; n < N; n++) A[m][n] -= temp * p[n];
            }

            // Обновляем V
            for (int m = 0; m < N; m++) {
                double temp = 0;
                for (int n = j; n < N; n++) temp += V[m][n] * p[n];
                temp *= mu;
                for (int n = j; n < N; n++) V[m][n] -= temp * p[n];
            }
        }
    }

    // Зануление элемента ниже диагонали (нижний треугольник)
    private static void deleteLowerTriangle(double[][] A, double[][] U, int I, int J) {
        int N = A.length;
        if (Math.abs(A[I][J]) > EPS) {
            double r = Math.sqrt(A[I][J] * A[I][J] + A[J][J] * A[J][J]);
            double c = A[J][J] / r;
            double s = A[I][J] / r;
            for (int k = 0; k < N; k++) {
                double temp1 = c * A[J][k] + s * A[I][k];
                double temp2 = c * A[I][k] - s * A[J][k];
                A[J][k] = temp1;
                A[I][k] = temp2;

                temp1 = c * U[k][J] + s * U[k][I];
                temp2 = c * U[k][I] - s * U[k][J];
                U[k][J] = temp1;
                U[k][I] = temp2;
            }
        }
        A[I][J] = 0;
    }

    // Зануление элемента выше диагонали (верхний треугольник)
    private static void deleteUpperTriangle(double[][] A, double[][] V, int I, int J) {
        int N = A.length;
        if (Math.abs(A[I][J]) > EPS) {
            double r = Math.sqrt(A[I][J] * A[I][J] + A[I][I] * A[I][I]);
            double c = A[I][I] / r;
            double s = -A[I][J] / r;
            for (int k = 0; k < N; k++) {
                double temp1 = c * A[k][I] - s * A[k][J];
                double temp2 = c * A[k][J] + s * A[k][I];
                A[k][I] = temp1;
                A[k][J] = temp2;

                temp1 = c * V[k][I] - s * V[k][J];
                temp2 = c * V[k][J] + s * V[k][I];
                V[k][I] = temp1;
                V[k][J] = temp2;
            }
        }
    }

    // Приведение всех сингулярных чисел к положительным
    private static void absSingularValues(double[][] Sigma, double[][] U) {
        int N = Sigma.length;
        for (int i = 0; i < N; i++) {
            if (Sigma[i][i] < 0) {
                Sigma[i][i] = -Sigma[i][i];
                for (int j = 0; j < N; j++) {
                    U[j][i] = -U[j][i];
                }
            }
        }
    }

    // Сортировка сингулярных чисел по убыванию
    private static void sortSingularValues(double[][] Sigma, double[][] U, double[][] V) {
        int N = Sigma.length;
        for (int I = 0; I < N; I++) {
            double max = Sigma[I][I];
            int idx = I;
            for (int j = I + 1; j < N; j++) {
                if (Sigma[j][j] > max) {
                    max = Sigma[j][j];
                    idx = j;
                }
            }
            if (I != idx) {
                double temp = Sigma[I][I];
                Sigma[I][I] = Sigma[idx][idx];
                Sigma[idx][idx] = temp;

                swapColumns(U, I, idx);
                swapColumns(V, I, idx);
            }
        }
    }

    // Перестановка столбцов в матрице
    private static void swapColumns(double[][] M, int col1, int col2) {
        int N = M.length;
        for (int i = 0; i < N; i++) {
            double temp = M[i][col1];
            M[i][col1] = M[i][col2];
            M[i][col2] = temp;
        }
    }
}
