using ConsoleChMethod;

public class Spline
{
    // Узлы сетки по X и Y
    readonly Vector x;
    readonly Vector y;

    // Количество интервалов между узлами
    int intervalCount;

    // Коэффициенты кубического сплайна
    Vector coeffA = new Vector(0);  // A
    Vector coeffB = new Vector(0);  // B
    Vector coeffC = new Vector(0);  // C
    Vector coeffD = new Vector(0);  // D

    // Конструктор: принимает точки (x, y) и строит сплайн
    public Spline(Vector x, Vector y)
    {
        this.x = x.Copy();
        this.y = y.Copy();

        if (x.GetSize() != y.GetSize())
            throw new ArgumentException("Количество аргументов не совпадает с количеством значений функции");

        intervalCount = x.GetSize() - 1;  // Число интервалов на 1 меньше числа точек

        CreateSpline();  // Построение сплайна
    }

    // Метод построения кубического сплайна
    private void CreateSpline()
    {
        // Вычисляем длины интервалов между соседними точками
        double[] intervalLengths = new double[intervalCount];
        for (int i = 0; i < intervalCount; i++)
            intervalLengths[i] = x[i + 1] - x[i];

        // Инициализируем коэффициенты A: A[i] = y[i+1]
        coeffA = new Vector(intervalCount);
        for (int i = 0; i < intervalCount; i++)
            coeffA[i] = y[i + 1];

        // diagA — нижняя диагональ трехдиагональной матрицы
        double[] diagA = new double[intervalCount - 1];
        diagA[0] = 0;
        for (int i = 1; i < intervalCount - 1; i++)
            diagA[i] = intervalLengths[i];

        // diagC — главная диагональ
        double[] diagC = new double[intervalCount - 1];
        for (int i = 0; i < intervalCount - 1; i++)
            diagC[i] = 2 * (intervalLengths[i] + intervalLengths[i + 1]);

        // diagB — верхняя диагональ
        double[] diagB = new double[intervalCount - 1];
        for (int i = 0; i < intervalCount - 2; i++)
            diagB[i] = intervalLengths[i + 1];
        diagB[intervalCount - 2] = 0;

        // rhsVector — правая часть системы уравнений
        double[] rhsVector = new double[intervalCount - 1];
        for (int i = 0; i < intervalCount - 1; i++)
        {
            double slopeNext = (y[i + 2] - y[i + 1]) / intervalLengths[i + 1];  // Наклон следующего отрезка
            double slopePrev = (y[i + 1] - y[i]) / intervalLengths[i];         // Наклон предыдущего отрезка
            rhsVector[i] = 6 * (slopeNext - slopePrev);                        // Формула для С
        }

        // Решаем систему методом прогонки
        double[] resultVector = Matrix.TridiagonalMatrixSolver(diagA, diagC, diagB, rhsVector);

        // Заполняем коэффициенты C
        coeffC = new Vector(intervalCount + 1);
        for (int i = 0; i < intervalCount - 1; i++)
            coeffC[i + 1] = resultVector[i];

        coeffC[0] = 0;                     // Граничное условие: C[0] = 0
        coeffC[intervalCount] = 0;         // Граничное условие: C[n] = 0

        // Заполняем коэффициенты D
        coeffD = new Vector(intervalCount);
        for (int i = 0; i < intervalCount; i++)
        {
            if (i == 0)
                coeffD[i] = coeffC[i] / intervalLengths[i];
            else
                coeffD[i] = (coeffC[i] - coeffC[i - 1]) / intervalLengths[i];
        }

        // Заполняем коэффициенты B
        coeffB = new Vector(intervalCount);
        for (int i = 0; i < intervalCount; i++)
        {
            coeffB[i] = (intervalLengths[i] / 2) * coeffC[i + 1]
                        - (Math.Pow(intervalLengths[i], 2) / 6) * coeffD[i]
                        + (y[i + 1] - y[i]) / intervalLengths[i];
        }
    }

    // Метод интерполяции в точке xt
    public double Interpolate(double xt)
    {
        // Если точка вне диапазона — возвращаем NaN
        if (xt < x[0] || xt > x[intervalCount])
            return double.NaN;

        int interval = 0;

        // Находим нужный интервал
        for (int i = 0; i < intervalCount; i++)
        {
            if (xt >= x[i] && xt <= x[i + 1])
            {
                interval = i;
                break;
            }
        }

        // Вычисляем значение сплайна в точке xt
        double dx = xt - x[interval];
        double value = y[interval] +
                       coeffB[interval] * dx +
                       (coeffC[interval] / 2) * Math.Pow(dx, 2) +
                       (coeffD[interval] / 6) * Math.Pow(dx, 3);

        return value;
    }
}