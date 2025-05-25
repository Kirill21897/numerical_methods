public static class IntegralMethods
{
    // === Метод средних прямоугольников ===
    // Вычисляет определённый интеграл функции f на интервале [a, b]
    // с заданной точностью eps методом средних прямоугольников.
    public static double MidpointIntegral(double a, double b, double eps, Func<double, double> f)
    {
        int n = 1;           // Начальное количество разбиений
        double prev = 0;     // Предыдущее значение интеграла для сравнения

        do
        {
            n *= 2;          // Увеличиваем количество разбиений вдвое
            double h = (b - a) / n;  // Шаг между узлами
            double sum = 0;

            // Суммирование значений функции в серединах отрезков
            for (int i = 0; i < n; i++)
            {
                double mid = a + h * (i + 0.5);  // Середина текущего подынтервала
                sum += f(mid);
            }

            double current = sum * h;  // Текущее приближение интеграла

            // Если достигнута нужная точность — возвращаем результат
            if (Math.Abs(current - prev) < eps)
                return current;

            prev = current;  // Сохраняем текущее значение для следующей итерации

        } while (true);  // Цикл продолжается до достижения точности
    }

    // === Метод трапеций ===
    // Вычисляет определённый интеграл методом трапеций.
    public static double TrapezoidIntegral(double a, double b, double eps, Func<double, double> f)
    {
        int n = 1;
        double prev = 0;

        do
        {
            n *= 2;
            double h = (b - a) / n;
            double sum = 0.5 * (f(a) + f(b));  // Граничные значения делим пополам

            // Суммируем промежуточные значения
            for (int i = 1; i < n; i++)
            {
                sum += f(a + i * h);
            }

            double current = sum * h;

            if (Math.Abs(current - prev) < eps)
                return current;

            prev = current;

        } while (true);
    }

    // === Метод Симпсона ===
    // Используется для численного интегрирования с высокой точностью.
    public static double SimpsonIntegral(double a, double b, double eps, Func<double, double> f)
    {
        int n = 2;  // Должно быть чётным
        double prev = 0;

        do
        {
            n *= 2;
            double h = (b - a) / n;
            double sumOdd = 0, sumEven = 0;

            // Разделение суммы на "чётные" и "нечётные" точки
            for (int i = 1; i < n; i++)
            {
                double x = a + i * h;
                if (i % 2 == 1)
                    sumOdd += f(x);
                else
                    sumEven += f(x);
            }

            // Формула Симпсона: (h/3) * [f(a) + f(b) + 4*sumOdd + 2*sumEven]
            double current = (h / 3) * (f(a) + f(b) + 4 * sumOdd + 2 * sumEven);

            if (Math.Abs(current - prev) < eps)
                return current;

            prev = current;

        } while (true);
    }

    // === Двойной интеграл методом Симпсона по обоим направлениям ===
    // Вычисляет двойной интеграл ∫∫f(x,y) dx dy по прямоугольной области [x1,x2]×[y1,y2]
    public static double DoubleIntegral(
        double x1, double x2,
        double y1, double y2,
        double eps,
        Func<double, double, double> f)
    {
        int nx = 2, ny = 2;
        double prev = 0;

        while (true)
        {
            double hx = (x2 - x1) / nx;
            double hy = (y2 - y1) / ny;
            double total = 0;

            // Перебор всех узлов сетки
            for (int i = 0; i <= nx; i++)
            {
                double x = x1 + i * hx;
                for (int j = 0; j <= ny; j++)
                {
                    double y = y1 + j * hy;

                    // Весовые коэффициенты по правилам Симпсона
                    int wx = (i == 0 || i == nx) ? 1 : (i % 2 == 0) ? 2 : 4;
                    int wy = (j == 0 || j == ny) ? 1 : (j % 2 == 0) ? 2 : 4;

                    total += wx * wy * f(x, y);
                }
            }

            // Общая формула двойного интеграла методом Симпсона
            double current = total * hx * hy / 9;

            if (Math.Abs(current - prev) < eps)
                return current;

            prev = current;
            nx *= 2;
            ny *= 2;  // Уточняем сетку
        }
    }

    // === Метод Чебышева для интегрирования (n=5 узлов) ===
    // Использует специальные узлы и веса, подходящие для полиномов до 9-й степени.
    public static double Chebyshev5(Func<double, double> f, double a, double b)
    {
        // Узлы на стандартном интервале [-1, 1]
        List<double> nodes = new List<double>
        {
            -Math.Sqrt((7 + 2 * Math.Sqrt(7)) / 15),
            -Math.Sqrt((7 - 2 * Math.Sqrt(7)) / 15),
            0,
            Math.Sqrt((7 - 2 * Math.Sqrt(7)) / 15),
            Math.Sqrt((7 + 2 * Math.Sqrt(7)) / 15)
        };

        // Соответствующие веса
        List<double> weights = new List<double>
        {
            128.0 / 225.0,
            (32 * (14 + Math.Sqrt(7))) / 225.0,
            16.0 / 225.0,
            (32 * (14 + Math.Sqrt(7))) / 225.0,
            128.0 / 225.0
        };

        // Преобразование от [a, b] к [-1, 1], чтобы использовать узлы Чебышева
        double c1 = (b - a) / 2;   // Коэффициент масштабирования
        double c2 = (a + b) / 2;   // Коэффициент сдвига

        double integral = 0;
        for (int i = 0; i < nodes.Count; i++)
        {
            double x_i = c1 * nodes[i] + c2;  // Преобразование узла из [-1, 1] в [a, b]
            integral += weights[i] * f(x_i);  // Сумма с учётом веса
        }

        return integral * c1;  // Умножение на длину интервала для корректного масштабирования
    }
}