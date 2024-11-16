#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <getopt.h>
using namespace std;

// Перечисление типов элементов
enum element_type { unknown, linear, cubic = 3 };

// Перечисление степеней ограничений
enum restrict_grade { first, second };

// Структура для передачи опций программе
typedef struct {
    element_type type;  // Тип элемента (линейный или кубический)
    int elemAmount;     // Количество конечных элементов
    int help = 0;       // Флаг отображения справки
} options_struct;

// Структура для описания граничных условий
typedef struct {
    restrict_grade grade; // Степень ограничения (первый или второй род)
    double pos;           // Позиция границы
    double val;           // Значение ограничения
} restriction;

// Прототипы функций
void printUsage();
double realSolve(double x);
double realSolvePrime(double x);
void printMatrix(vector<double> &matrix, int len);
options_struct getOptions(const int argc, char *argv[]);
int solveSLAU(vector<double> &retXVector, vector<double> inAMatrix, vector<double> inBVector, int size);
int restrictType2(vector<double> &matrix, vector<double> &vector, int a, int len, int line, double value);
int restrictType1(vector<double> &matrix, vector<double> &vector, int a, int len, int line, double value);
int makeLinearSLAU(vector<double> &resultMatrix, vector<double> &resultVector, int size, double L, double a, double b, double c, double d);
int makeCubicSLAU(vector<double> &resultMatrix, vector<double> &resultVector, int size, double L, double a, double b, double c, double d);
void exportStiffnessMatrixAndLoadVector(const vector<double> &stiffnessMatrix, const vector<double> &loadVector, int size,
                                        const string &matrixFilename = "stiffness_matrix.txt",
                                        const string &vectorFilename = "load_vector.txt");
void export_data(const int size, options_struct settings, vector<pair<double, double>> plotMKE,
                 vector<pair<double, double>> plotReal, vector<pair<double, double>> plotAbsolutelyReal,
                 vector<double> displacements, vector<double> errors);

// Новая функция для решения расширенной СЛАУ
int solveExtendedSLAU(vector<double> &retXVector, vector<double> &retDerivatives,
                      vector<double> inAMatrix, vector<double> inBVector, int size) {
    // Прямой ход метода Гаусса
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            double d = inAMatrix[j * size + i] / inAMatrix[i * size + i];
            for (int k = i; k < size; k++)
                inAMatrix[j * size + k] -= d * inAMatrix[i * size + k];
            inBVector[j] -= d * inBVector[i];
        }
    }

    // Обратный ход метода Гаусса для определения производных
    for (int i = size - 1; i >= 0; --i) {
        double rSum = 0;
        for (int j = i + 1; j < size; ++j)
            rSum += retXVector[j] * inAMatrix[i * size + j];
        retXVector[i] = (inBVector[i] - rSum) / inAMatrix[i * size + i];
    }

    // Сохранение производной в конце для повторного использования
    retDerivatives[0] = retXVector[0];  // Производная слева
    retDerivatives[1] = retXVector[size - 1];  // Производная справа

    return 0;
}

// Основная программа
int main(int argc, char **argv) {
    // Получаем настройки из командной строки
    const options_struct settings = getOptions(argc, argv);

    // Данные по задаче
    const double a = 7, b = 0, c = -15, d = 60;
    restriction lower = {second, -7, 0};  // Начальное условие второго рода
    restriction upper = {second, 4, 4};   // Граничное условие второго рода

    // Размеры сетки и шаг
    int size = settings.elemAmount * settings.type + 1;  // Количество узлов
    double step = (upper.pos - lower.pos) / settings.elemAmount;

    // Инициализация векторов и матриц
    vector<double> stiffnessMatrix(size * size, 0);  // Матрица жесткости
    vector<double> loadVector(size, 0);              // Вектор нагрузок
    vector<double> displacements(size, 0);           // Вектор смещений
    vector<double> derivatives(2, 0);                // Производные в крайних точках

    // Построение матрицы жесткости и вектора нагрузок
    if (settings.type == linear) {
        makeLinearSLAU(stiffnessMatrix, loadVector, size, step, a, b, c, d);
    } else if (settings.type == cubic) {
        makeCubicSLAU(stiffnessMatrix, loadVector, size, step, a, b, c, d);
    } else {
        cout << "Ошибка: не задан тип конечных элементов." << endl;
        return -1;
    }

    // Первый этап: решение расширенной СЛАУ для нахождения производной слева
    solveExtendedSLAU(displacements, derivatives, stiffnessMatrix, loadVector, size);

    // Используем вычисленное значение производной для обновления граничного условия
    lower.val = derivatives[0];  // Устанавливаем производную как значение граничного условия
    restrictType2(stiffnessMatrix, loadVector, a, size, 0, lower.val);  // Применяем новое граничное условие
    restrictType2(stiffnessMatrix, loadVector, a, size, size - 1, upper.val);  // Условие второго рода на правой границе

    // Второй этап: решаем задачу заново с новым граничным условием
    solveSLAU(displacements, stiffnessMatrix, loadVector, size);

    // Вычисление погрешностей и формирование данных для графиков
    vector<double> errors;
    double maxError = 0;
    double node = lower.pos;
    vector<pair<double, double>> plotMKE;
    vector<pair<double, double>> plotReal;

    for (int i = 0; i < size; i++) {
        double realU = realSolve(node);
        double error = fabs(displacements[i] - realU);
        if (error > maxError) maxError = error;
        errors.push_back(error);
        plotMKE.push_back(make_pair(node, displacements[i]));
        plotReal.push_back(make_pair(node, realU));
        node += step / settings.type;
    }

    // Вывод результатов
    cout << "Производная слева (после расчета): " << derivatives[0] << endl;
    cout << "Решение (смещения):" << endl;
    for (int i = 0; i < size; ++i) {
        cout << "u(" << lower.pos + i * step / settings.type << ") = " << displacements[i] << endl;
    }
    cout << "Максимальная ошибка: " << maxError << endl;

    // Экспорт результатов
    export_data(size, settings, plotMKE, plotReal, plotReal, displacements, errors);
    exportStiffnessMatrixAndLoadVector(stiffnessMatrix, loadVector, size);

    return 0;
}

// Вспомогательные функции
void printMatrix(vector<double> &matrix, int len) {
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            cout << matrix[i * len + j] << " ";
        }
        cout << endl;
    }
}
