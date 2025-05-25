using System;

public class Vector
{
    protected int size;         // Размерность вектора
    protected double[] data;    // Массив, хранящий элементы вектора
    static Random rnd = new Random();  // Для генерации случайных чисел

    // Свойство для получения размерности вектора
    public int Size { get { return size; } }

    // Конструктор: создаёт вектор заданного размера, заполненный нулями
    public Vector(int size)
    {
        this.size = size;
        data = new double[size];
    }

    // Получить массив данных вектора
    public double[] GetElements()
    {
        return data;
    }

    // Конструктор копирования: создаёт вектор из массива
    public Vector(double[] v)
    {
        this.size = v.Length;
        data = new double[size];
        for (int i = 0; i < size; i++) data[i] = v[i];
    }

    // Индексатор: доступ к элементам через [i]
    public double this[int index]
    {
        get { return data[index]; }
        set { data[index] = value; }
    }

    // Возвращает размер вектора
    public int GetSize() { return size; }

    // Устанавливает значение элемента по индексу
    public bool SetElement(double el, int index)
    {
        if (index < 0 || index >= size) return false;
        data[index] = el;
        return true;
    }

    // Получает значение элемента по индексу
    public double GetElement(int index)
    {
        if (index < 0 || index >= size) return default(double);
        return data[index];
    }

    // Создаёт копию вектора
    public Vector Copy()
    {
        Vector rez = new Vector(data);
        return rez;
    }

    // Преобразование в строку (для вывода в консоль)
    public override string ToString()
    {
        return $"({string.Join(", ", data)})";
    }

    // Евклидова норма вектора
    public double Norma1()
    {
        double s = 0;
        for (int i = 0; i < size; i++)
            s += data[i] * data[i];
        return Math.Sqrt(s);
    }

    // Скалярное произведение с другим вектором
    public double ScalarMultiply(Vector b)
    {
        if (size != b.size) return 0;
        double s = 0;
        for (int i = 0; i < size; i++)
            s += data[i] * b.data[i];
        return s;
    }

    // Умножение вектора на скаляр
    public Vector MultiplyScalar(double c)
    {
        Vector rez = new Vector(size);
        for (int i = 0; i < size; i++) rez.data[i] = data[i] * c;
        return rez;
    }

    // Нормализация вектора (деление на его длину)
    public Vector Normalize()
    {
        Vector rez = new Vector(size);
        double d = Norma1();
        for (int i = 0; i < size; i++)
            rez.data[i] = (d != 0) ? data[i] / d : 0;
        return rez;
    }

    // Генерирует случайный нормализованный вектор
    public static Vector NormalizeRandom(int size)
    {
        Vector rez = new Vector(size);
        for (int i = 0; i < size; i++)
            rez.data[i] = (rnd.NextDouble() - 0.5) * 2.0;
        return rez.Normalize();
    }

    // Отрицание вектора (-vector)
    public Vector UMinus()
    {
        Vector rez = new Vector(size);
        for (int i = 0; i < size; i++) rez.data[i] = -data[i];
        return rez;
    }

    // Сложение векторов
    public Vector Add(Vector c)
    {
        Vector rez = new Vector(size);
        for (int i = 0; i < size; i++) rez.data[i] = data[i] + c.data[i];
        return rez;
    }

    // Вычитание векторов
    public Vector Minus(Vector c)
    {
        Vector rez = new Vector(size);
        for (int i = 0; i < size; i++) rez.data[i] = data[i] - c.data[i];
        return rez;
    }

    // Перегрузка оператора +: сложение двух векторов
    public static Vector operator +(Vector a, Vector b)
    {
        if (a.size == b.size)
        {
            Vector c = new Vector(a.size);
            for (int i = 0; i < a.size; i++)
                c[i] = a[i] + b[i];
            return c;
        }
        throw new ArgumentException("Vectors must be of the same size.");
    }

    // Перегрузка оператора -: вычитание двух векторов
    public static Vector operator -(Vector a, Vector b)
    {
        if (a.size == b.size)
        {
            Vector c = new Vector(a.size);
            for (int i = 0; i < a.size; i++)
                c[i] = a[i] - b[i];
            return c;
        }
        throw new ArgumentException("Vectors must be of the same size.");
    }

    // Перегрузка оператора *: умножение вектора на число
    public static Vector operator *(Vector a, double c)
    {
        Vector r = new Vector(a.size);
        for (int i = 0; i < a.size; i++)
            r[i] = a[i] * c;
        return r;
    }

    // Перегрузка оператора *: умножение числа на вектор
    public static Vector operator *(double c, Vector a)
    {
        Vector r = new Vector(a.size);
        for (int i = 0; i < a.size; i++)
            r[i] = a[i] * c;
        return r;
    }

    // Перегрузка оператора *: скалярное произведение двух векторов
    public static double operator *(Vector a, Vector b)
    {
        if (a.size == b.size)
        {
            double s = 0.0;
            for (int i = 0; i < a.size; i++)
                s += a[i] * b[i];
            return s;
        }
        throw new ArgumentException("Vectors must be of the same size.");
    }

    // Перегрузка оператора /: деление вектора на число
    public static Vector operator /(Vector a, double c)
    {
        Vector result = new Vector(a.Size);
        for (int i = 0; i < a.Size; i++)
            result[i] = a[i] / c;
        return result;
    }

    // Перегрузка оператора /: деление вектора на целое число
    public static Vector operator /(Vector a, int c)
    {
        return a / (double)c;
    }
}