using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleChMethod
{
    class Matrix
    {
        protected int rows, columns;           // Количество строк и столбцов
        protected double[,] data;              // Данные матрицы
        private const double Eps = 0.000001;   // Точность для сравнений

        // Конструктор: создаёт матрицу заданного размера (r x c), заполненную нулями
        public Matrix(int r, int c)
        {
            this.rows = r;
            this.columns = c;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) data[i, j] = 0;
        }

        // Конструктор: создаёт матрицу на основе двумерного массива
        public Matrix(double[,] mm)
        {
            this.rows = mm.GetLength(0);
            this.columns = mm.GetLength(1);
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = mm[i, j];
        }

        // Конструктор копирования
        public Matrix(Matrix M)
        {
            this.rows = M.rows;
            this.columns = M.columns;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = M[i, j];
        }

        // Свойства для получения размеров матрицы
        public int Rows { get { return rows; } }
        public int Columns { get { return columns; } }

        // Индексатор: доступ к элементам матрицы через [i,j]
        public double this[int i, int j]
        {
            get
            {
                if (i < 0 || j < 0 || i >= rows || j >= columns)
                    return Double.NaN;
                else
                    return data[i, j];
            }
            set
            {
                if (i < 0 || j < 0 || i >= rows || j >= columns)
                    return;
                else
                    data[i, j] = value;
            }
        }

        // Получить строку матрицы как Vector
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

        // Получить столбец матрицы как Vector
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

        // Установить строку матрицы
        public bool SetRow(int index, Vector r)
        {
            if (index < 0 || index > rows) return false;
            if (r.Size != columns) return false;
            for (int k = 0; k < columns; k++) data[index, k] = r[k];
            return true;
        }

        // Установить столбец матрицы
        public bool SetColumn(int index, Vector c)
        {
            if (index < 0 || index > columns) return false;
            if (c.Size != rows) return false;
            for (int k = 0; k < rows; k++) data[k, index] = c[k];
            return true;
        }

        // Норма матрицы (евклидова норма)
        public double Norma1()
        {
            double s = 0;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    s += data[i, j] * data[i, j];
            return Math.Sqrt(s);
        }

        // Максимум сумм по строкам
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

        // Максимум сумм по столбцам
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

        // Умножение матрицы на вектор
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

        // Печать матрицы в консоль
        public void Print()
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    Console.Write($"{data[i, j]} \t");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        // Преобразование в строку
        public override string ToString()
        {
            string s = "{\n";
            for (int i = 0; i < rows; i++)
                s += GetRow(i).ToString() + "\n";
            s += "}";
            return s;
        }

        // Транспонирование матрицы
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

        // Умножение матрицы на число
        public static Matrix MultByNum(Matrix m, double c)
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

        public static Matrix operator *(Matrix m, double c) => MultByNum(m, c);
        public static Matrix operator *(double c, Matrix m) => MultByNum(m, c);

        // Умножение двух матриц
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

        // Сложение матриц
        public static Matrix operator +(Matrix m1, Matrix m2)
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

        // Вычитание матриц
        public static Matrix operator -(Matrix m1, Matrix m2)
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

        // Отрицание матрицы
        public static Matrix operator -(Matrix m1)
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

        // Создание единичной матрицы
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

        // Перестановка строк
        public void SwapRows(int r1, int r2)
        {
            if (r1 < 0 || r2 < 0 || r1 >= rows || r2 >= rows || (r1 == r2)) return;
            Vector v1 = GetRow(r1);
            Vector v2 = GetRow(r2);
            SetRow(r2, v1);
            SetRow(r1, v2);
        }

        // Перестановка столбцов
        public void SwapColumns(int c1, int c2)
        {
            if (c1 < 0 || c2 < 0 || c1 >= columns || c2 >= columns || (c1 == c2)) return;
            Vector v1 = GetColumn(c1);
            Vector v2 = GetColumn(c2);
            SetColumn(c2, v1);
            SetColumn(c1, v2);
        }

        // Решение системы с нижнетреугольной матрицей
        public static Vector Solve_LU_DOWN_Treug(Matrix a, Vector b)
        {
            int rows = a.rows;
            int columns = a.columns;
            if (columns != rows || rows != b.Size) return null!;
            for (int i = 0; i < rows; i++)
            {
                if (a.data[i, i] == 0) return null!;
                for (int j = i + 1; j < rows; j++)
                    if (Math.Abs(a.data[i, j]) > Eps) return null!;
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

        // Копирование матрицы
        public Matrix Copy()
        {
            Matrix r = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) r[i, j] = data[i, j];
            return r;
        }

        // Решение системы с верхнетреугольной матрицей
        public static Vector SolveUpperTriangular(Matrix matrix, Vector b)
        {
            int n = b.Size;
            Vector x = new Vector(n);
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = b[i];
                for (int j = i + 1; j < n; j++)
                    x[i] -= matrix[i, j] * x[j];
                x[i] /= matrix[i, i];
            }
            return x;
        }

        // Обращение верхнетреугольной матрицы
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

        // Решение системы с нижнетреугольной матрицей
        public static Vector SolveLowerTriangular(Matrix matrix, Vector b)
        {
            int n = b.Size;
            Vector x = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                x[i] = b[i];
                for (int j = 0; j < i; j++)
                    x[i] -= matrix[i, j] * x[j];
                x[i] /= matrix[i, i];
            }
            return x;
        }

        // Обращение нижнетреугольной матрицы
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

        // Решение СЛАУ методом Гаусса
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
                // Перестановка строк
                for (int k = i; k < n + 1; k++)
                {
                    double tmp = augmentedMatrix[maxRow, k];
                    augmentedMatrix[maxRow, k] = augmentedMatrix[i, k];
                    augmentedMatrix[i, k] = tmp;
                }
                // Обнуление ниже диагонали
                for (int k = i + 1; k < n; k++)
                {
                    double c = -augmentedMatrix[k, i] / augmentedMatrix[i, i];
                    for (int j = i; j < n + 1; j++)
                        augmentedMatrix[k, j] += c * augmentedMatrix[i, j];
                }
            }

            // Обратный ход
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = augmentedMatrix[i, n] / augmentedMatrix[i, i];
                for (int j = i - 1; j >= 0; j--)
                    augmentedMatrix[j, n] -= augmentedMatrix[j, i] * x[i];
            }

            return x;
        }

        // Обращение матрицы методом Гаусса
        public static Matrix obrSolveGaussian(Matrix matrix)
        {
            int n = matrix.Rows;
            Matrix augmentedMatrix = new Matrix(n, 2 * n);

            // Инициализация расширенной матрицы
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    augmentedMatrix[i, j] = matrix[i, j];
                augmentedMatrix[i, i + n] = 1;
            }

            // Прямой ход
            for (int i = 0; i < n; i++)
            {
                // Поиск максимального элемента
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

                // Перестановка строк
                for (int k = 0; k < 2 * n; k++)
                {
                    double tmp = augmentedMatrix[maxRow, k];
                    augmentedMatrix[maxRow, k] = augmentedMatrix[i, k];
                    augmentedMatrix[i, k] = tmp;
                }

                // Обнуление ниже главного элемента
                for (int k = i + 1; k < n; k++)
                {
                    double c = -augmentedMatrix[k, i] / augmentedMatrix[i, i];
                    for (int j = 0; j < 2 * n; j++)
                        augmentedMatrix[k, j] += c * augmentedMatrix[i, j];
                }
            }

            // Обратный ход
            for (int i = n - 1; i >= 0; i--)
            {
                double diagElement = augmentedMatrix[i, i];
                for (int j = 0; j < 2 * n; j++)
                    augmentedMatrix[i, j] /= diagElement;

                for (int j = i - 1; j >= 0; j--)
                {
                    double factor = augmentedMatrix[j, i];
                    for (int k = 0; k < 2 * n; k++)
                        augmentedMatrix[j, k] -= factor * augmentedMatrix[i, k];
                }
            }

            // Извлечение обратной матрицы
            Matrix inverseMatrix = new Matrix(n, n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    inverseMatrix[i, j] = augmentedMatrix[i, j + n];

            return inverseMatrix;
        }

        // QR-разложение методом Грама–Шмидта
        public static (Matrix Q, Matrix R) QR_decomposition_by_Gram_Schmidt_method(Matrix input_matrix)
        {
            int n = input_matrix.Rows;
            int m = input_matrix.Columns;
            if (n != m) return (null!, null!);

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
                    column_to_Q_matrix = column_to_Q_matrix.Minus(Q_matrix.GetColumn(j).MultiplyScalar(dotProduct));
                }

                double norm = column_to_Q_matrix.Norma1();
                if (norm < 1e-10)
                    throw new InvalidOperationException("Столбцы матрицы линейно зависимы.");
                R_matrix[i, i] = norm;
                Q_matrix.SetColumn(i, column_to_Q_matrix.MultiplyScalar(1.0 / norm));
            }

            return (Q_matrix, R_matrix);
        }

        // Решение СЛАУ через QR-разложение
        public static Vector Solve_LU_OrtGS(Matrix A, Vector b)
        {
            var (Q, R) = QR_decomposition_by_Gram_Schmidt_method(A);
            Matrix QT = Q.Transpose();
            Vector b_res = QT * b;
            return SolveUpperTriangular(R, b_res);
        }

        // Метод последовательных приближений
        public static Vector MPP(Matrix A, Vector B)
        {
            if (A.Rows != A.Columns || A.Rows != B.Size) return null!;
            int n = A.Rows;
            double eps = 1e-8;
            int maxRow;
            //double tmp;

            for (int j = 0; j < n; j++)
            {
                maxRow = j;
                for (int i = j + 1; i < n; i++)
                    if (Math.Abs(A[i, j]) > Math.Abs(A[maxRow, j]))
                        maxRow = i;

                if (maxRow != j)
                {
                    A.SwapRows(maxRow, j);
                    (B[maxRow], B[j]) = (B[j], B[maxRow]);
                }

                if (Math.Abs(A[j, j]) < eps) return null!;
            }

            Vector beta = new Vector(n);
            for (int i = 0; i < n; i++)
                beta[i] = B[i] / A[i, i];

            Matrix alpha = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    alpha[i, j] = i != j ? A[i, j] / A[i, i] : 0;
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
            } while (delta.Norma1() > eps);

            return prev_x;
        }

        // Градиентный спуск
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
                count++;
            }
            return x;
        }

        // Псевдообратная матрица (через Гаусса)
        public static Matrix Pseudo_inverse_matrix(Matrix AA)
        {
            Matrix A = AA.Copy();
            Double deter = Determinant(AA);
            if (deter == 0 || A.Rows != A.Columns) return null!;

            int n = A.Rows;
            Matrix c = new Matrix(n, n);
            Matrix X = new Matrix(n, n);
            for (int i = 0; i < n; i++) c[i, i] = 1;

            for (int k = 0; k < n; k++)
            {
                double maxEl = Math.Abs(A[k, k]);
                int maxRow = k;
                for (int i = k + 1; i < n; i++)
                    if (Math.Abs(A[i, k]) > maxEl)
                    {
                        maxEl = Math.Abs(A[i, k]);
                        maxRow = i;
                    }

                for (int j = 0; j < n; j++)
                {
                    (A[maxRow, j], A[k, j]) = (A[k, j], A[maxRow, j]);
                    (c[maxRow, j], c[k, j]) = (c[k, j], c[maxRow, j]);
                }

                for (int i = k; i < n; i++)
                {
                    if (i == k) continue;
                    double mull = A[i, k] / A[k, k];
                    for (int j = k; j < n; j++)
                    {
                        A[i, j] -= A[k, j] * mull;
                        c[i, j] -= c[k, j] * mull;
                    }
                }
            }

            for (int i = 0; i < n; i++)
            {
                Vector tmp = new Vector(n);
                for (int j = 0; j < n; j++)
                    tmp[j] = c[j, i];
                Vector tmp2 = Solve_UP_Down_Tri(A, tmp);
                for (int j = 0; j < n; j++)
                    X[j, i] = tmp2[j];
            }

            return X;
        }

        // Метод Якоби
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
                        if (i != j)
                            xNew[i] -= A[i, j] * x[j];
                    xNew[i] /= A[i, i];
                }

                double maxDiff = 0;
                for (int i = 0; i < n; i++)
                    maxDiff = Math.Max(maxDiff, Math.Abs(xNew[i] - x[i]));

                x = xNew.Copy();
                if (maxDiff < eps)
                {
                    Console.WriteLine($"Сошлось за {k + 1} итераций.");
                    break;
                }
            }

            return x;
        }

        // Определитель матрицы
        public static double Determinant(Matrix AA)
        {
            Matrix A = AA.Copy();
            int n = A.Rows;
            for (int k = 0; k < n; k++)
            {
                double maxEl = Math.Abs(A[k, k]);
                int maxRow = k;
                for (int i = k + 1; i < n; i++)
                    if (Math.Abs(A[i, k]) > maxEl)
                    {
                        maxEl = Math.Abs(A[i, k]);
                        maxRow = i;
                    }

                for (int j = k; j < n; j++)
                {
                    (A[maxRow, j], A[k, j]) = (A[k, j], A[maxRow, j]);
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
                answer *= A[i, i];
            return answer;
        }

        // Вспомогательный метод для псевдообратной матрицы
        public static Vector Solve_UP_Down_Tri(Matrix a, Vector b)
        {
            int rows = a.rows;
            int columns = a.columns;
            if (columns != rows || rows != b.Size) return null!;

            for (int i = rows - 1; i >= 0; i--)
            {
                if (a.data[i, i] == 0) return null!;
                for (int j = 0; j < i; j++)
                    if (Math.Abs(a.data[i, j]) > Eps) return null!;
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

        // Решение трёхдиагональной системы методом прогонки
        public static double[] TridiagonalMatrixSolver(double[] a, double[] b, double[] c, double[] d)
        {
            int n = d.Length;
            double[] p = new double[n];
            double[] q = new double[n];
            double[] solution = new double[n];

            p[0] = -c[0] / b[0];
            q[0] = d[0] / b[0];

            for (int i = 1; i < n; i++)
            {
                double denom = b[i] + a[i - 1] * p[i - 1];
                p[i] = -c[i] / denom;
                q[i] = (d[i] - a[i - 1] * q[i - 1]) / denom;
            }

            solution[n - 1] = q[n - 1];
            for (int i = n - 2; i >= 0; i--)
                solution[i] = p[i] * solution[i + 1] + q[i];

            return solution;
        }
    }
}