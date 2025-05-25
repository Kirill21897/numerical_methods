using ConsoleChMethod;
namespace ConsoleApp1;

class NaimQuadrat
{
    // Входные значения X (аргументы функции)
    private readonly Vector argumentX;

    // Выходные значения Y (значения целевой функции)
    private readonly Vector argumentY;

    // Рассчитанные параметры модели (коэффициенты регрессии)
    public Vector parameters = new Vector(0);

    // Делегат для базисной функции
    public delegate double PsiFunction(double value);

    // Массив базисных функций, например: 1, x, x^2 и т.д.
    public PsiFunction[] psiFunctions;

    // Количество точек данных
    public int argumentCount;

    // Количество базисных функций
    public int psiFunctionCount;


    // Конструктор
    public NaimQuadrat(Vector x, Vector y, PsiFunction[] functions)
    {
        if (x.Size != y.Size)
            throw new Exception("Количество аргументов не совпадает с количеством значений функции");

        this.argumentX = x;
        this.argumentY = y;
        this.argumentCount = x.Size;

        this.psiFunctions = functions;
        this.psiFunctionCount = functions.Length;

        CalculateParameters(); // Расчёт параметров методом наименьших квадратов
    }

    // Основной метод — вычисляет параметры модели по методу наименьших квадратов
    private void CalculateParameters()
    {
        // Создаем матрицу H, где каждая строка — значения базисных функций для x_i
        Matrix H = new Matrix(argumentCount, psiFunctionCount);
        for (int i = 0; i < argumentCount; i++)
        {
            Vector row = GetFunctionValues(argumentX[i]); // Получаем значения базисных функций для x_i
            H.SetRow(i, row); // Устанавливаем строку матрицы
        }

        // Транспонируем матрицу H
        Matrix Ht = H.Transpose();

        // D = Ht * H — матрица системы нормальных уравнений
        Matrix D = Ht * H;

        // Находим псевдообратную матрицу вместо обычной обратной
        Matrix D_PseudoInverse = Matrix.Pseudo_inverse_matrix(D);

        // Q = D^+ * Ht — матрица преобразования
        Matrix Q = D_PseudoInverse * Ht;

        // parameters = Q * Y — финальный расчёт параметров
        parameters = Q * argumentY;
    }

    // Возвращает вектор значений базисных функций для заданного x
    public Vector GetFunctionValues(double x)
    {
        Vector result = new Vector(psiFunctionCount);
        for (int i = 0; i < psiFunctionCount; i++)
        {
            result[i] = psiFunctions[i](x); // Применяем каждую функцию к x
        }
        return result;
    }

    // Считает критерий качества — L1-норму остатков
    public double GetCriterion()
    {
        Vector residuals = new Vector(argumentCount);
        for (int i = 0; i < argumentCount; i++)
        {
            Vector phi = GetFunctionValues(argumentX[i]); // Значения базисных функций для x_i
            residuals[i] = argumentY[i] - (parameters * phi); // Разница между реальным и предсказанным y_i
        }
        return residuals.Norma1(); // Возвращаем L1-норму вектора остатков
    }
}