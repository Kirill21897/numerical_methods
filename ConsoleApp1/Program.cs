using System;
using System.ComponentModel.DataAnnotations;
using System.Security.AccessControl;
using ConsoleChMethod;

// Объявление пространства имён проекта
namespace ConsoleApp1;

// Основной класс программы
class Program
{
    // Точка входа в приложение
    static void Main(string[] agrs)
    {
        // Флаг для работы цикла меню
        bool flag = true;
        
        while (flag)
        {
            // Вывод меню выбора теста
            Console.WriteLine("Выберете тест:\n" +
                              "1) Нахождение корня\n" +
                              "2) Вектора\n" +
                              "3) Матрицы\n" +
                              "4) Сплайн\n" +
                              "5) Метод наименьших квадратов\n" +
                              "6) Дифференциальные уравнения\n" +
                              "7) Интегралы\n" +
                              "8) Выход");

            // Чтение и преобразование пользовательского ввода
            string ex_string = Console.ReadLine()!;
            int ex_number = int.Parse(ex_string);

            // Обработка выбора пользователя через switch-case
            switch(ex_number)
            {
                case 1:
                    Console.WriteLine("\n-----Нахождение корня-----");
                    try
                    {
                        // Пример функции: f(x) = cos(x), для которой ищем корень
                        Func<double, double> func = x => Math.Cos(x);

                        // Параметры методов поиска корней
                        double leftBorder = 0.9;   // Левая граница интервала
                        double rightBorder = 2.7;  // Правая граница интервала
                        double eps = 1e-6;         // Точность

                        // Метод деления пополам (бисекции)
                        double rootBisection = EquationSolution.BisectionMethod(leftBorder, rightBorder, eps, func);
                        Console.WriteLine($"Корень методом деления пополам: {rootBisection}");

                        // Метод Ньютона (итерационный, требует начального приближения)
                        double rootNewton = EquationSolution.NewtonMethod(1, eps, func);
                        Console.WriteLine($"Корень методом Ньютона: {rootNewton}");

                        // Метод простых итераций
                        double rootIter = EquationSolution.SimpleIterations(1, eps, func);
                        Console.WriteLine($"Корень методом SimpleIterations: {rootIter}");

                        // Метод Стеффенсена (ускоряет сходимость итераций)
                        double rootStef = EquationSolution.SteffMethod(1, eps, func);
                        Console.WriteLine($"Корень методом Стеффенсена: {rootStef}");

                        Console.WriteLine("--------------------------");
                    }
                    catch (InvalidOperationException ex)
                    {
                        Console.WriteLine($"\n{ex}");
                        Console.WriteLine("--------------------------");
                    }

                    flag = false;
                    break;

                case 2:
                    Console.WriteLine("\n-----Вектора-----");

                    // Создание двух векторов
                    Vector v1 = new Vector(new double[] { 1.0, 2.0, 3.0 });
                    Vector v2 = new Vector(new double[] { 4.0, 5.0, 6.0 });

                    // Вывод векторов
                    Console.WriteLine("Vector v1: " + v1);
                    Console.WriteLine("Vector v2: " + v2);

                    // Сложение векторов
                    Vector v3 = v1 + v2;
                    Console.WriteLine("v1 + v2: " + v3);

                    // Разность векторов
                    Vector v4 = v1 - v2;
                    Console.WriteLine("v1 - v2: " + v4);

                    // Умножение вектора на скаляр
                    double scalar = 2.0;
                    Vector v5 = v1 * scalar;
                    Console.WriteLine("v1 * " + scalar + ": " + v5);

                    // Скалярное произведение векторов
                    double dotProduct = v1 * v2;
                    Console.WriteLine("v1 . v2: " + dotProduct);

                    // Нормализация вектора
                    Vector v6 = v1.Normalize();
                    Console.WriteLine("Normalized v1: " + v6);

                    // Генерация случайного нормированного вектора
                    Vector v7 = Vector.NormalizeRandom(3);
                    Console.WriteLine("Random normalized vector: " + v7);

                    Console.WriteLine("-----------------");
                    flag = false;
                    break;

                case 3:
                    Console.WriteLine("\n-----Матрицы-----");

                    // Верхняя треугольная матрица и правая часть
                    double[,] upperMatrixData = {
                        { 2, -1, 1 },
                        { 0, 3, -1 },
                        { 0, 0, 1 }
                    };
                    double[] bUpper = { 3, 7, 2 };

                    Matrix upperTriangularMatrix = new Matrix(upperMatrixData);
                    Vector BU = new Vector(bUpper);

                    Console.WriteLine("(Решение верхней треугольной матрицы)");
                    Console.WriteLine("Верхняя треугольная матрица:");
                    upperTriangularMatrix.Print();

                    // Решение системы с верхней треугольной матрицей
                    Vector solutionUpper = Matrix.SolveUpperTriangular(upperTriangularMatrix, BU);
                    Console.WriteLine($"Ответ: {solutionUpper}");
                    Console.WriteLine($"--------");

                    // Нахождение обратной матрицы к верхней треугольной
                    Console.WriteLine("(Нахождение обратной матрицы к верхней треугольной)");
                    Matrix ObrUp = Matrix.ObrSolveUpperTriangular(upperTriangularMatrix);
                    ObrUp.Print();
                    Console.WriteLine($"--------");

                    // Аналогично для нижней треугольной матрицы
                    double[,] lowerMatrixData = {
                        { 2, 0, 0 },
                        { 1, 3, 0 },
                        { 1, -1, 1 }
                    };
                    double[] bLower = { 2, 5, 3 };

                    Matrix lowerTriangularMatrix = new Matrix(lowerMatrixData);
                    Vector BL = new Vector(bLower);

                    Console.WriteLine("(Решение нижней треугольной матрицы)");
                    Console.WriteLine("Нижняя треугольная матрица:");
                    lowerTriangularMatrix.Print();

                    Vector solutionLower = Matrix.SolveLowerTriangular(lowerTriangularMatrix, BL);
                    Console.WriteLine($"Ответ: {solutionLower}");
                    Console.WriteLine($"--------");

                    Console.WriteLine("(Нахождение обратной матрицы к нижней треугольной)");
                    Matrix ObrLow = Matrix.ObrSolveLowerTriangular(lowerTriangularMatrix);
                    ObrLow.Print();
                    Console.WriteLine($"--------");

                    // Решение СЛАУ методом Гаусса
                    double[,] augmentedMatrix = {
                        { 3, 2, -5, -1},
                        { 2, -1, 3, 13},
                        { 1, 2, -1, 9}
                    };

                    Matrix matrix = new Matrix(augmentedMatrix);
                    Console.WriteLine("(Метод Гаусса)");
                    Console.WriteLine("Расширенная матрица:");
                    matrix.Print();

                    Vector solutionGauss = Matrix.SolveGaussian(matrix);
                    Console.WriteLine($"Ответ: {solutionGauss}");
                    Console.WriteLine($"--------");

                    // Нахождение обратной матрицы методом Гаусса
                    Matrix ObrGauss = Matrix.obrSolveGaussian(lowerTriangularMatrix);
                    Console.WriteLine("Нахождение обратной матрицы методом Гаусса");
                    Console.WriteLine("Исходная матрица:");
                    lowerTriangularMatrix.Print();
                    Console.WriteLine("Обратная матрица");
                    ObrGauss.Print();
                    Console.WriteLine($"--------");

                    // QR-разложение методом Грама–Шмидта
                    Console.WriteLine("QR - разложение методом Грама - Шмидта\n");
                    Console.WriteLine("Начальная матрица:");
                    lowerTriangularMatrix.Print();

                    var (Q, R) = Matrix.QR_decomposition_by_Gram_Schmidt_method(lowerTriangularMatrix);
                    Console.WriteLine("Матрица Q:");
                    Q.Print();
                    Console.WriteLine("Матрица R:");
                    R.Print();
                    Console.WriteLine("Произведение матриц Q * R:");
                    (Q * R).Print();
                    Console.WriteLine($"--------");

                    // Решение СЛАУ методом LU с ортогонализацией
                    double[,] augmentedMatrix_ = {
                        { 3, 2, -5 },
                        { 2, -1, 3 },
                        { 1, 2, -1 }
                    };
                    double[] b_ = { 2, 5, 3 };
                    Vector B_ = new Vector(b_);
                    Matrix matrix_ = new Matrix(augmentedMatrix_);

                    Console.WriteLine("Solve_LU_OrtGS");
                    Console.Write("Исходная матрица:\n");
                    matrix_.Print();
                    Console.WriteLine("Вектор:");
                    Console.WriteLine($"{B_}\n");

                    Vector res = Matrix.Solve_LU_OrtGS(matrix_, B_);
                    Console.WriteLine($"Результат решения\n{res}");
                    Console.WriteLine($"--------");

                    // Метод последовательных приближений
                    Console.WriteLine("Метод последовательных приближений");
                    double[,] augmentedMatrix_mpp = {
                        { 4, -1, 1 },
                        { 2, 3, -2 },
                        { -1, 2, 4 }
                    };
                    double[] b_mpp = { 1, 2, 3 };
                    Vector B_mpp = new Vector(b_mpp);
                    Matrix matrix_mpp = new Matrix(augmentedMatrix_mpp);

                    Console.WriteLine($"Исходная матрица:\n");
                    matrix_mpp.Print();

                    Vector solution_mpp = Matrix.MPP(matrix_mpp, B_mpp);
                    Console.WriteLine($"Результат:\n{solution_mpp}");
                    Console.WriteLine("-----------------");

                    // Метод градиентного спуска
                    Console.WriteLine("Метод наименьшего спуска (Метод градиентного спуска)");
                    double[,] augmentedMatrix_grad = {
                        { 4, -1, 1 },
                        { 2, 3, -2 },
                        { -1, 2, 4 }
                    };
                    double[] b_grad = {1, 2, 3};
                    Vector B_garad = new Vector(b_grad);
                    Matrix matrix_grad = new Matrix(augmentedMatrix_grad);

                    Console.WriteLine($"Исходная матрица:\n");
                    matrix_mpp.Print();

                    double eps_grad = 1e-6;
                    double step = 0.01;
                    Vector solution_grad = Matrix.GradientDescentSolver(matrix_grad, B_garad, step, eps_grad);
                    Console.WriteLine($"Результат:\n{solution_grad}");
                    Console.WriteLine("-----------------");

                    // Псевдообратная матрица методом Гревиля
                    Console.WriteLine("Нахождение псевдообратной матрицы методом Гревиля");
                    double[,] augmentedMatrix_Grevile = {
                        { 1, 2 },
                        { 3, 4 },
                        { 5, 6 }
                    };
                    Matrix matrix_Grevile = new Matrix(augmentedMatrix_Grevile);

                    Console.WriteLine($"Исходная матрица:\n");
                    matrix_Grevile.Print();

                    Matrix pseudo_inverse_matrix_Grevile = Matrix.Pseudo_inverse_matrix(matrix_Grevile);
                    Console.WriteLine("Результат:");
                    pseudo_inverse_matrix_Grevile.Print();
                    Console.WriteLine("-----------------");

                    // Метод Якоби
                    Console.WriteLine("Решение СЛАУ методом Якоби");
                    double[,] augmentedMatrix_jacoby = {
                        { 4, -1, 1 },
                        { 2, 3, -2 },
                        { -1, 2, 4 }
                    };
                    double[] b_jacoby = {1, 2, 3};
                    Vector B_jacoby = new Vector(b_jacoby);
                    Matrix matrix_jacoby = new Matrix(augmentedMatrix_jacoby);

                    Console.WriteLine($"Исходная матрица:\n");
                    matrix_jacoby.Print();

                    double eps_jacoby = 1e-6;
                    int mxiterations_jacoby = 1000;
                    Vector solution_Jacoby = Matrix.Jacoby_method(matrix_jacoby, B_jacoby, mxiterations_jacoby, eps_jacoby);
                    Console.WriteLine($"Результат:\n{solution_Jacoby}");
                    Console.WriteLine("-----------------");

                    flag = false;
                    break;

                case 4:
                    Console.WriteLine("\n-----Сплайн-----");

                    // Задаем точки для построения сплайна
                    double[] arr_x = { 0, 1, 2, 3, 4 };
                    double[] arr_y = { 0, 1, 4, 9, 16 };

                    Vector x = new Vector(arr_x);
                    Vector y = new Vector(arr_y);

                    // Создаем сплайн
                    var spline = new Spline(x, y);

                    // Интерполяция значений
                    for (double xi = 0; xi <= 4; xi += 0.5)
                    {
                        Console.WriteLine($"x = {xi:F2}, y = {spline.Interpolate(xi):F4}");
                    }

                    Console.WriteLine("-----------------");
                    flag = false;
                    break;

                case 5:
                    Console.WriteLine("\n-----МНК-----");

                    // Исходные данные
                    Vector x_ = new Vector(new double[] { 0.0, 1.0, 2.0, 3.0, 4.0 });
                    Vector y_ = new Vector(new double[] { 1.1, 2.9, 5.0, 8.8, 16.2 });

                    // Базисные функции для МНК
                    NaimQuadrat.PsiFunction[] basisFunctions = new NaimQuadrat.PsiFunction[]
                    {
                        val => 1,           // Константа
                        val => val,         // Линейная
                        val => val * val    // Квадратичная
                    };

                    // Создаем объект МНК
                    NaimQuadrat nq = new NaimQuadrat(x_, y_, basisFunctions);

                    // Вывод результатов
                    Console.WriteLine("Найденные параметры (коэффициенты):");
                    Console.WriteLine(nq.parameters);
                    Console.WriteLine("Значение функционала (невязка):");
                    Console.WriteLine(nq.GetCriterion());

                    Console.WriteLine("Предсказанные значения:");
                    for (int i = 0; i < x_.Size; i++)
                    {
                        double predicted = nq.parameters * nq.GetFunctionValues(x_[i]);
                        Console.WriteLine($"x = {x_[i]} -> y_pred = {predicted:F4}");
                    }

                    Console.WriteLine("-----------------");
                    flag = false;
                    break;

                case 6:
                    Console.WriteLine("\n-----Дифференциальные уравнения-----");

                    // Начальные условия для гармонического осциллятора
                    Vector initialState = new Vector(new double[] { 1, 0 }); // x0 = 1, v0 = 0
                    double t0 = 0;
                    double tEnd = 10;
                    int steps = 100;

                    Console.WriteLine("Рассмотрим движение гармонического осциллятора — это как маятник или пружина.");
                    Console.WriteLine("Уравнение: x'' = -x\n");
                    Console.WriteLine("Начальные условия:");
                    Console.WriteLine("x(0) = 1   (начальное положение)");
                    Console.WriteLine("v(0) = 0   (начальная скорость)\n");
                    Console.WriteLine("Мы будем решать это уравнение на интервале от t=0 до t=10,");
                    Console.WriteLine("с шагом по времени и с использованием нескольких численных методов.\n");

                    // Решение методом Рунге–Кутты 4-го порядка
                    Matrix resultRK4 = EquationSolver.SolveRungeKutta4(t0, tEnd, initialState, steps, EquationSolver.Oscillator);
                    Console.WriteLine("Метод Рунге–Кутты 4-го порядка — самый точный из используемых.");
                    Console.WriteLine("t       | Положение (x)    | Скорость (v)");
                    Console.WriteLine("--------------------------------------------------");
                    for (int i = 0; i <= steps; i += steps / 10)
                    {
                        Vector col = resultRK4.GetColumn(i);
                        Console.WriteLine($"{col[0]:F2}     | {col[1]:F6}        | {col[2]:F6}");
                    }
                    Console.WriteLine("\n...");

                    // Метод Эйлера
                    Matrix resultEuler = EquationSolver.SolveEuler(t0, tEnd, initialState, steps, EquationSolver.Oscillator);
                    Console.WriteLine("\nМетод Эйлера — простой, но менее точный.");
                    Console.WriteLine("t       | Положение (x)    | Скорость (v)");
                    Console.WriteLine("--------------------------------------------------");
                    for (int i = 0; i <= steps; i += steps / 10)
                    {
                        Vector col = resultEuler.GetColumn(i);
                        Console.WriteLine($"{col[0]:F2}     | {col[1]:F6}        | {col[2]:F6}");
                    }
                    Console.WriteLine("\n...");

                    // Метод Адамса
                    Matrix resultAdams = EquationSolver.SolveAdams(t0, tEnd, initialState, steps, EquationSolver.Oscillator);
                    Console.WriteLine("\nМетод Адамса — многошаговый, требует начального приближения.");
                    Console.WriteLine("t       | Положение (x)    | Скорость (v)");
                    Console.WriteLine("--------------------------------------------------");
                    for (int i = 0; i <= steps; i += steps / 10)
                    {
                        Vector col = resultAdams.GetColumn(i);
                        Console.WriteLine($"{col[0]:F2}     | {col[1]:F6}        | {col[2]:F6}");
                    }

                    // Модифицированный метод Рунге–Кутты 2-го порядка
                    Matrix resultModifiedRK2 = EquationSolver.SolveRungeKutta2_Modified(t0, tEnd, initialState, steps, EquationSolver.Oscillator);
                    Console.WriteLine("\nМодифицированный метод Рунге–Кутты 2-го порядка — улучшенная версия классического RK2.");
                    Console.WriteLine("Использует предиктор-корректор для повышения точности.");
                    Console.WriteLine("t       | Положение (x)    | Скорость (v)");
                    Console.WriteLine("--------------------------------------------------");
                    for (int i = 0; i <= steps; i += steps / 10)
                    {
                        Vector col = resultModifiedRK2.GetColumn(i);
                        Console.WriteLine($"{col[0]:F2}     | {col[1]:F6}        | {col[2]:F6}");
                    }

                    Console.WriteLine("-----------------");
                    flag = false;
                    break;

                case 7:
                    Console.WriteLine("\n-----Интергралы-----");

                    // Целевая функция: f(x) = sin(x)
                    static double Function(double x) => Math.Sin(x);

                    // Точное значение интеграла ∫sin(x)dx от a до b
                    static double ExactIntegral(double a, double b) => -Math.Cos(b) + Math.Cos(a);

                    double a = 0;
                    double b = Math.PI;

                    // Точное значение
                    double exact = ExactIntegral(a, b);

                    // Численные методы интегрирования
                    double rect = IntegralMethods.MidpointIntegral(a, b, 1e-6, Function);
                    double trap = IntegralMethods.TrapezoidIntegral(a, b, 1e-6, Function);
                    double simp = IntegralMethods.SimpsonIntegral(a, b, 1e-6, Function);
                    double cheb = IntegralMethods.Chebyshev5(Function, a, b);

                    // Вывод результатов
                    Console.WriteLine("=== Численное интегрирование функции f(x) = sin(x) ===");
                    Console.WriteLine($"Интеграл на интервале [{a:F4}, {b:F4}]");
                    Console.WriteLine();
                    Console.WriteLine($"Метод прямоугольников: {rect:F8}");
                    Console.WriteLine($"Метод трапеций:       {trap:F8}");
                    Console.WriteLine($"Метод Симпсона:        {simp:F8}");
                    Console.WriteLine($"Метод Чебышева (n=5):  {cheb:F8}");
                    Console.WriteLine($"Точное значение:        {exact:F8}");
                    Console.WriteLine($"Погрешность (Симпсон):  {Math.Abs(simp - exact):F8}");
                    Console.WriteLine($"Погрешность (Чебышев):  {Math.Abs(cheb - exact):F8}");
                    Console.WriteLine("-----------------");

                    flag = false;
                    break;

                case 8:
                    flag = false;
                    break;

                default:
                    Console.WriteLine("Ошибка! Введите верный номер теста");
                    break;
            }
        }
    }
}