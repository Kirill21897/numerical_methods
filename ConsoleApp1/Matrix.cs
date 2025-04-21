using System;
using System.Collections.Generic;
using System.Collections;
using System.Text;
using System.Linq;
using System.Numerics;
using System.Reflection.Metadata;
using System.Xml.Linq;
using System.Collections.Immutable;
using System.ComponentModel.Design;
using System.Xml.XPath;
using System.Diagnostics.Tracing;

namespace ConsoleChMethod
{
   class Matrix
    {
        protected int rows, columns;
        protected double[,] data;
        private const double Eps = 0.000001;
        public Matrix(int r, int c)
        {
            this.rows = r; this.columns = c;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) data[i, j] = 0;
        }
        public Matrix(double[,] mm)
        {
            this.rows = mm.GetLength(0); this.columns = mm.GetLength(1);
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = mm[i, j];
        }
        public Matrix(Matrix M)
        {
            this.rows = M.rows; this.columns = M.columns;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = M[i, j];
        }
        public int Rows { get { return rows; } }
        public int Columns { get { return columns; } }

        public double this[int i, int j]
        {
            get
            {
                if (i < 0 || j < 0 || i >= rows || j >= columns)
                {
                    // Console.WriteLine(" Индексы вышли за пределы матрицы ");
                    return Double.NaN;
                }
                else
                    return data[i, j];
            }
            set
            {
                if (i < 0 || j < 0 || i >= rows || j >= columns)
                {
                    //Console.WriteLine(" Индексы вышли за пределы матрицы ");
                }
                else
                    data[i, j] = value;
            }
        }
        public Vector GetRow(int r)
        {
            if (r >= 0 && r < rows)
            {
                Vector row = new Vector(columns);
                for (int j = 0; j < columns; j++) row[j] = data[r, j];
                return row;
            }
            return null!;
        }
        public Vector GetColumn(int c)
        {
            if (c >= 0 && c < columns)
            {
                Vector column = new Vector(rows);
                for (int i = 0; i < rows; i++) column[i] = data[i, c];
                return column;
            }
            return null!;
        }
        public bool SetRow(int index, Vector r)
        {
            if (index < 0 || index > rows) return false;
            if (r.Size != columns) return false;
            for (int k = 0; k < columns; k++) data[index, k] = r[k];
            return true;
        }
        public bool SetColumn(int index, Vector c)
        {
            if (index < 0 || index > columns) return false;
            if (c.Size != rows) return false;
            for (int k = 0; k < rows; k++) data[k, index] = c[k];
            return true;
        }
        public double Norma1()
        {
            double s = 0;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    s += data[i, j] * data[i, j];
            return Math.Sqrt(s);
        }
        public double Norma2()
        {
            double max = 0, s = 0;
            for (int i = 0; i < rows; i++)
            {
                s = 0;
                for (int j = 0; j < columns; j++)
                    s += Math.Abs(data[i, j]);
                if (s > max) max = s;
            }
            return max;
        }
        public double Norma3()
        {
            double max = 0, s = 0;
            for (int j = 0; j < columns; j++)
            {
                s = 0;
                for (int i = 0; i < rows; i++)
                    s += Math.Abs(data[i, j]);
                if (s > max) max = s;
            }
            return max;
        }
        //умножение матрицы на вектор
        public static Vector operator *(Matrix a, Vector b)
        {
            if (a.columns != b.Size) return null!;
            Vector r = new Vector(a.rows);
            for (int i = 0; i < a.rows; i++)
            {
                r[i] = a.GetRow(i) * b;
            }
            return r;
        }
        //печать
        public void Print()
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    Console.Write($"{data[i, j]} \t");             //t - горизонтальная табуляция
                }
                Console.WriteLine();
            }
            Console.WriteLine("\n");
        }
        public override string ToString()
        {
            string s = "{\n";
            for (int i = 0; i < rows; i++)
                s += GetRow(i).ToString() + "\n";
            s += "}";
            return @s;
        }
        //Транспонированние
        public Matrix Transpose()
        {
            Matrix transposeMatrix = new Matrix(columns, rows);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    transposeMatrix.data[j, i] = data[i, j];
                }
            }
            return transposeMatrix;
        }
        
        //Умножение матрицы на чило
        public static Matrix MultByNum(Matrix m, double c)   //Умножаем матрицу на число
        {
            Matrix result = new Matrix(m.rows, m.columns);
            for (int i = 0; i < m.rows; i++)
            {
                for (int j = 0; j < m.columns; j++)
                {
                    result[i, j] = m[i, j] * c;
                }
            }
            return result;
        }
        public static Matrix operator *(Matrix m, double c)
        {
            return Matrix.MultByNum(m, c);
        }

        public static Matrix operator *(double c, Matrix m)
        {
            return Matrix.MultByNum(m, c);
        }
        //Умножение матриц
        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            if (m1.columns != m2.rows)
            {
                throw new Exception("Количество столбцов первой матрицы не равно количеству строк второй");
            }
            Matrix result = new Matrix(m1.rows, m2.columns);
            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m2.columns; j++)
                {
                    result[i, j] = 0;
                    for (int k = 0; k < m1.columns; k++)
                    {
                        result[i, j] += m1[i, k] * m2[k, j];
                    }
                }
            }
            return result;
        }
       
        //Сложение матриц
        public static Matrix operator +(Matrix m1, Matrix m2)     //Сложение матриц
        {
            if (m1.rows != m2.rows || m1.columns != m2.columns)
            {
                throw new Exception("Матрицы не совпадают по размерности");
            }
            Matrix result = new Matrix(m1.rows, m1.columns);

            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m2.columns; j++)
                {
                    result[i, j] = m1[i, j] + m2[i, j];
                }
            }
            return result;
        }
        //Вычитание матриц
        public static Matrix operator -(Matrix m1, Matrix m2)     //Вычитание матриц
        {
            if (m1.rows != m2.rows || m1.columns != m2.columns)
            {
                throw new Exception("Матрицы не совпадают по размерности");
            }
            Matrix result = new Matrix(m1.rows, m2.columns);

            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m2.columns; j++)
                {
                    result[i, j] = m1[i, j] - m2[i, j];
                }
            }
            return result;
        }
        public static Matrix operator -(Matrix m1)     //Отрицание матриц
        {

            Matrix result = new Matrix(m1.rows, m1.columns);

            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m1.columns; j++)
                {
                    result[i, j] = -m1[i, j];
                }
            }
            return result;
        }
        public static Matrix EdMatrix(int size)
        {
            Matrix tm = new Matrix(size, size);
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                    tm.data[i, j] = 0.0;
                tm.data[i, i] = 1.0;
            }
            return tm;
        }
        public void SwapRows(int r1, int r2)
        {
            if (r1 < 0 || r2 < 0 || r1 >= rows || r2 >= rows || (r1 == r2)) return;
            Vector v1 = GetRow(r1);
            Vector v2 = GetRow(r2);
            SetRow(r2, v1);
            SetRow(r1, v2);
        }
        public void SwapColumns(int c1, int c2)
        {
            if (c1 < 0 || c2 < 0 || c1 >= columns || c2 >= columns || (c1 == c2)) return;
            Vector v1 = GetColumn(c1);
            Vector v2 = GetColumn(c2);
            SetColumn(c2, v1);
            SetColumn(c1, v2);
        }
              public static Vector Solve_LU_DOWN_Treug(Matrix a, Vector b)
        {
            int rows = a.rows; int columns = a.columns;
            if (columns != rows || rows != b.Size) return null!;
            for (int i = 0; i < rows; i++)
            {
                if (a.data[i, i] == 0) return null!;
                for (int j = i + 1; j < rows; j++)
                    if (Math.Abs(a.data[i, j] )> Eps) return null!;
            }
            Vector x = new Vector(rows);
            x[0] = b[0] / a.data[0, 0];
            for (int i = 1; i < rows; i++)
            {
                double s = 0;
                for (int k = 0; k < i; k++)
                    s += a.data[i, k] * x[k];
                x[i] = (b[i] - s) / a.data[i, i];
            }
            return x;
        }
        public Matrix Copy()
        {
            Matrix r = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) r[i, j] = data[i, j];
            return r;
        }
/*
        public double[] SolveUpperTriangular(double[,] matrix, double[] b)
        {
            int n = b.Length;
            double[] x = new double[n];

            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = b[i];
                for (int j = i + 1; j < n; j++)
                {
                    x[i] -= matrix[i, j] * x[j];
                }
                x[i] /= matrix[i, i];
            }

            return x;
        }
*/
        /*public static double ScalarDot (Vector a, Vector b){
            double scalarDot = 0;
            if (a.Size == b.Size){
                for (int i = 0; i < a.Size; i++){
                    scalarDot += (a[i] * b[i]);
                }
            }
            return scalarDot;
        }*/
        public static Vector SolveUpperTriangular(Matrix matrix, Vector b)
        {
            int n = b.Size;
            Vector x = new Vector(n);

            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = b[i];
                for (int j = i + 1; j < n; j++)
                {
                    x[i] -= matrix[i, j] * x[j];
                }
                x[i] /= matrix[i, i];
            }

            return x;
        }

        public static Matrix ObrSolveUpperTriangular(Matrix matrix)
        {
            int n = matrix.Rows;
            Matrix obr = new Matrix(matrix.Rows, matrix.Columns);
            for (int i = 0; i < n; i++)
            {
                Vector x = new Vector(n);
                Vector b = new Vector(n); 
                b[i] = 1;
                x = SolveUpperTriangular(matrix, b);
                obr.SetColumn(i, x);
            }        
            return obr;
        }

        public static Vector SolveLowerTriangular(Matrix matrix, Vector b)
        {
            int n = b.Size;
            Vector x = new Vector(n);

            for (int i = 0; i < n; i++)
            {
                x[i] = b[i];
                for (int j = 0; j < i; j++)
                {
                    x[i] -= matrix[i, j] * x[j];
                }
                x[i] /= matrix[i, i];
            }

            return x;
        }

        public static Matrix ObrSolveLowerTriangular(Matrix matrix)
        {
            int n = matrix.Rows;
            Matrix Obr = new Matrix(matrix.Rows, matrix.Columns);
            for (int i = 0; i < n; i++)
            {
                Vector x = new Vector(n);
                Vector b = new Vector(n);
                b[i] = 1;
                x = SolveLowerTriangular(matrix, b);
                Obr.SetColumn(i, x); 
            }
            return Obr;
        }
        

        public static Vector SolveGaussian(Matrix augmentedMatrix)
        {
            int n = augmentedMatrix.Rows;
            Vector x = new Vector(n);

            // Прямой ход
            for (int i = 0; i < n; i++)
            {
                // Поиск максимального элемента в столбце
                double maxEl = Math.Abs(augmentedMatrix[i, i]);
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmentedMatrix[k, i]) > maxEl)
                    {
                        maxEl = Math.Abs(augmentedMatrix[k, i]);
                        maxRow = k;
                    }
                }

                // Меняем местами текущую строку и строку с максимальным элементом
                for (int k = i; k < n + 1; k++)
                {
                    double tmp = augmentedMatrix[maxRow, k];
                    augmentedMatrix[maxRow, k] = augmentedMatrix[i, k];
                    augmentedMatrix[i, k] = tmp;
                }

                // Обнуляем элементы ниже главной диагонали
                for (int k = i + 1; k < n; k++)
                {
                    double c = -augmentedMatrix[k, i] / augmentedMatrix[i, i];
                    for (int j = i; j < n + 1; j++)
                    {
                        if (i == j)
                        {
                            augmentedMatrix[k, j] = 0;
                        }
                        else
                        {
                            augmentedMatrix[k, j] += c * augmentedMatrix[i, j];
                        }
                    }
                }
            }

            // Обратный ход
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = augmentedMatrix[i, n] / augmentedMatrix[i, i];
                for (int j = i - 1; j >= 0; j--)
                {
                    augmentedMatrix[j, n] -= augmentedMatrix[j, i] * x[i];
                }
            }

            return x;
        }

        public static Matrix obrSolveGaussian(Matrix matrix){
            int n = matrix.Rows;
            // Создаем расширенную матрицу, состоящую из исходной матрицы и единичной матрицы
            Matrix augmentedMatrix = new Matrix(n, 2 * n);
            // Заполняем левую часть (исходная матрица) и правую часть (единичная матрица)
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    augmentedMatrix[i, j] = matrix[i, j];
                }
                augmentedMatrix[i, i + n] = 1; // Единичная матрица
            }

            // Прямой ход
            for (int i = 0; i < n; i++)
            {
                // Поиск максимального элемента в столбце
                double maxEl = Math.Abs(augmentedMatrix[i, i]);
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmentedMatrix[k, i]) > maxEl)
                    {
                        maxEl = Math.Abs(augmentedMatrix[k, i]);
                        maxRow = k;
                    }
                }

                // Меняем местами текущую строку и строку с максимальным элементом
                for (int k = 0; k < 2 * n; k++)
                {
                    double tmp = augmentedMatrix[maxRow, k];
                    augmentedMatrix[maxRow, k] = augmentedMatrix[i, k];
                    augmentedMatrix[i, k] = tmp;
                }

                // Обнуляем элементы ниже главной диагонали
                for (int k = i + 1; k < n; k++)
                {
                    double c = -augmentedMatrix[k, i] / augmentedMatrix[i, i];
                    for (int j = 0; j < 2 * n; j++)
                    {
                        if (i == j)
                        {
                            augmentedMatrix[k, j] = 0;
                        }
                        else
                        {
                            augmentedMatrix[k, j] += c * augmentedMatrix[i, j];
                        }
                    }
                }
            }

            // Обратный ход
            for (int i = n - 1; i >= 0; i--)
            {
                // Нормализуем строку
                double diagElement = augmentedMatrix[i, i];
                for (int j = 0; j < 2 * n; j++)
                {
                    augmentedMatrix[i, j] /= diagElement;
                }

                for (int j = i - 1; j >= 0; j--)
                {
                    double factor = augmentedMatrix[j, i];
                    for (int k = 0; k < 2 * n; k++)
                    {
                        augmentedMatrix[j, k] -= factor * augmentedMatrix[i, k];
                    }
                }
            }

            // Извлекаем обратную матрицу из расширенной матрицы
            Matrix inverseMatrix = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    inverseMatrix[i, j] = augmentedMatrix[i, j + n]; // Правая часть
                }
            }

            return inverseMatrix;
        }

        public static (Matrix Q, Matrix R) QR_decomposition_by_Gram_Schmidt_method(Matrix input_matrix)
        {
            int n = input_matrix.Rows; 
            int m = input_matrix.Columns;
            if (n != m){
                return (null, null)!;
            }
            Matrix Q_matrix = new Matrix(n, m); 
            Matrix R_matrix = new Matrix(m, m); 

            Q_matrix.SetColumn(0, input_matrix.GetColumn(0));

            for (int i = 0; i < m; i++)
            {
                Vector column_to_Q_matrix = input_matrix.GetColumn(i);

                for (int j = 0; j < i; j++)
                {
                    double dotProduct = Q_matrix.GetColumn(j).ScalarMultiply(column_to_Q_matrix);
                    R_matrix[j, i] = dotProduct; 

                    Vector projection = Q_matrix.GetColumn(j).MultiplyScalar(dotProduct);
                    column_to_Q_matrix = column_to_Q_matrix.Minus(projection);
                }

                double norm = column_to_Q_matrix.Norma1(); 
                if (norm < 1e-10) 
                {
                    throw new InvalidOperationException("The input matrix has linearly dependent columns.");
                }

                R_matrix[i, i] = norm; 
                Q_matrix.SetColumn(i, column_to_Q_matrix.MultiplyScalar(1.0 / norm));
            }

            return (Q_matrix, R_matrix);
        }

        public static Vector Solve_LU_OrtGS (Matrix A, Vector b){
            var (Q, R) = QR_decomposition_by_Gram_Schmidt_method(A);
            Matrix QT = Q.Transpose();
            Vector b_res = QT * b;
            // Rx = D_obr * D * b
            Vector res = SolveUpperTriangular(R, b_res); 
            return res;
        }

         public static Vector MPP(Matrix A, Vector B)
        {
            if (A.Rows != A.Columns) return null!;
            if (A.Rows != B.Size) return null!;

            int n = A.Rows;
            double eps = 0.00000001;
            int max;
            double tmp;

            for (int j = 0; j < n; j++)
            {
                max = j;
                for (int i = j + 1; i < n; i++)
                {
                    if (Math.Abs(A[i, j]) > Math.Abs(A[max, j])) { max = i; };
                }

                if (max != j)
                {
                    Vector temp = A.GetRow(max); A.SetRow(max, A.GetRow(j)); A.SetRow(j, temp);
                    tmp = B[max]; B[max] = B[j]; B[j] = tmp;

                }
                if (Math.Abs(A[j, j]) < eps) return null!;

            }

            Vector beta = new Vector(n);
            for (int i = 0; i < beta.Size; i++)
            {
                beta[i] = B[i] / A[i, i];
            }

            Matrix alpha = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i != j) alpha[i, j] = A[i, j] / A[i, i];
                    else alpha[i, j] = 0;
                }
                beta[i] = B[i] / A[i, i];
            }

            Vector prev_x = beta;
            Vector current_x;
            Vector delta;

            do
            {
                current_x = beta - alpha * prev_x;
                delta = current_x - prev_x;
                prev_x = current_x;
            }
            while (delta.Norma1() > eps);

            return prev_x;
        }


        public static Vector GradientDescentSolver(Matrix a, Vector b, double step, double eps)
        {
            int n = b.Size;
            Vector x = new Vector(n);
            int count = 0;

            while (true)
            {
                Vector gradient = a.Transpose() * ((a * x).Minus(b));
                Vector x_new = x.Minus(gradient.MultiplyScalar(step));

                if ((x_new.Minus(x)).Norma1() < eps)
                {
                    Console.WriteLine($"Сходимость достигнута на итерации {count + 1}");
                    break;
                }

                x = x_new; 
                count += 1;
            }

            return x;
        }

        public static Matrix Pseudo_inverse_matrix(Matrix AA)
        {
            Matrix A = AA.Copy();
            Double deter = Matrix.Determinant(AA);
            if (deter == 0 || A.Rows != A.Columns)
            {
                return null!;
            }
            Matrix c = new Matrix(A.Rows, A.Columns);
            Matrix X = new Matrix(A.Rows, A.Columns);
            int n = A.Rows;
            for (int i = 0; i < n; i++)
            {
                c[i, i] = 1;
            }
            for (int k = 0; k < n; k++)
            {
                double maxEl = Math.Abs(A[k, k]);
                int maxRow = k;
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(A[i, k]) > maxEl)
                    {
                        maxEl = Math.Abs(A[i, k]);
                        maxRow = i;
                    }
                }

                for (int j = 0; j < n; j++) // Меняем местами k-ую и максимальную строки
                {
                    double tmp = A[maxRow, j];
                    double tmpc = c[maxRow, j];
                    A[maxRow, j] = A[k, j];
                    c[maxRow, j] = c[k, j];
                    A[k, j] = tmp;
                    c[k, j] = tmpc;
                }
                for (int i = k; i < n; i++) // Проходим по всем строкам 
                {
                    if (i == k) continue; // Пропускаем k-ую строку
                    double mull = A[i, k] / A[k, k]; // Делим элементы строки на максимальный элемент
                    for (int j = k; j < n; j++)
                    { // Вычитаем из каждого элемента строки i произведение элемента k-ой строки на mull
                        A[i, j] -= A[k, j] * mull;
                        c[i, j] -= c[k, j] * mull;
                    }
                }
            }
            // obr 
            for (int i = 0; i < n; i++)
            {
                Vector tmp = new Vector(A.Rows);
                for (int j = 0; j < n; j++)
                {
                    tmp[j] = c[j, i];
                }
                Vector tmp2 = Solve_UP_Down_Tri(A, tmp);
                for (int j = 0; j < n; j++)
                {
                    X[j, i] = tmp2[j];
                }
            }
            return X;
        }

        public static Vector Jacoby_method(Matrix A, Vector b, int maxIterations, double eps)
        {
            int n = b.Size;
            Vector x = new Vector(n);
            Vector xNew = new Vector(n);

            for (int k = 0; k < maxIterations; k++)
            {
                for (int i = 0; i < n; i++)
                {
                    xNew[i] = b[i];

                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                        {
                            xNew[i] -= A[i, j] * x[j];
                        }
                    }

                    xNew[i] /= A[i, i];
                }

                double maxDiff = 0;
                for (int i = 0; i < n; i++)
                {
                    maxDiff = Math.Max(maxDiff, Math.Abs(xNew[i] - x[i]));
                }

                x = xNew.Copy();

                if (maxDiff < eps)
                {
                    Console.WriteLine($"Сошлось за {k + 1} итераций.");
                    break;
                }
            }

            return x;
        }

        public static double Determinant(Matrix AA)
        {
            Matrix A = AA.Copy();
            int n = A.Rows;
            for (int k = 0; k < n; k++)
            {
                double maxEl = Math.Abs(A[k, k]);
                int maxRow = k;
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(A[i, k]) > maxEl)
                    {
                        maxEl = Math.Abs(A[i, k]);
                        maxRow = i;
                    }
                }

                for (int j = k; j < n; j++)
                {
                    double tmp = A[maxRow, j];
                    A[maxRow, j] = A[k, j];
                    A[k, j] = tmp;
                }
                for (int i = 0; i < n; i++)
                {
                    if (i == k) continue;
                    double mull = A[i, k] / A[k, k];
                    for (int j = k; j < n; j++)
                        A[i, j] -= A[k, j] * mull;
                }
            }
            double answer = 1;
            for (int i = 0; i < n; i++)
            {
                answer *= A[i, i];
            }
            return answer;
        }

        public static Vector Solve_UP_Down_Tri(Matrix a, Vector b)
        {
            int rows = a.rows; int columns = a.columns;
            if (columns != rows || rows != b.Size) return null!;
            for (int i = rows - 1; i >= 0; i--)
            {
                if (a.data[i, i] == 0) return null!;
                for (int j = 0; j < i; j++)
                    if (Math.Abs(a.data[i, j]) > 0.000000001) return null!;
            }
            Vector x = new Vector(rows);
            x[rows - 1] = b[rows - 1] / a.data[rows - 1, rows - 1];
            for (int i = rows - 2; i >= 0; i--)
            {
                double s = 0;
                for (int k = rows - 1; k > i; k--)
                    s += a.data[i, k] * x[k];
                x[i] = (b[i] - s) / a.data[i, i];
            }
            return x;
        }

    }
}
