#include <fstream>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <vector>
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
void printMatrix(vector<double> &matrix, int len);
options_struct getOptions(const int argc, char *argv[]);
int solveSLAU(vector<double> &retXVector, vector<double> inAMatrix, vector<double> inBVector, int size);
int restrictType2(vector<double> &matrix, vector<double> &vector, int a, int len, int line, double value);
int restrictType1(vector<double> &matrix, vector<double> &vector, int a, int len, int line, double value);
int makeLinearSLAU(vector<double> &resultMatrix, vector<double> &resultVector, int size, double L, double a, double b, double c, double d);
int makeCubicSLAU(vector<double> &resultMatrix, vector<double> &resultVector, int size, double L, double a, double b, double c, double d);

double realSolvePrime(double x) {
    double C2 = -0.000141861924;
    double C1 = 0.007826998720;
    double sqrtVal = sqrt(15.0 / 7.0);
    return C1 * sqrtVal * exp(sqrtVal * x) - C2 * sqrtVal * exp(-sqrtVal * x);
}
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

void exportStiffnessMatrixAndLoadVector(const vector<double>& stiffnessMatrix, 
                                        const vector<double>& loadVector, 
                                        int size, 
                                        const string& matrixFilename = "stiffness_matrix.txt", 
                                        const string& vectorFilename = "load_vector.txt") {
    // Открываем файл для матрицы жесткости
    ofstream matrixFile(matrixFilename);
    if (matrixFile.is_open()) {
        matrixFile << "Stiffness Matrix:\n";
        matrixFile << setprecision(4) << fixed;  // Устанавливаем формат вывода для чисел
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                matrixFile << setw(10) << stiffnessMatrix[i * size + j] << " ";
            }
            matrixFile << "\n";
        }
        matrixFile.close();
        cout << "Stiffness matrix exported to " << matrixFilename << endl;
    } else {
        cerr << "Unable to open file for stiffness matrix export." << endl;
    }

    // Открываем файл для вектора нагрузок
    ofstream vectorFile(vectorFilename);
    if (vectorFile.is_open()) {
        vectorFile << "Load Vector:\n";
        vectorFile << setprecision(4) << fixed;  // Устанавливаем формат вывода для чисел
        for (int i = 0; i < size; ++i) {
            vectorFile << "Node " << setw(3) << i << ": " << setw(10) << loadVector[i] << "\n";
        }
        vectorFile.close();
        cout << "Load vector exported to " << vectorFilename << endl;
    } else {
        cerr << "Unable to open file for load vector export." << endl;
    }
}

// addBoundaryConditionWithoutPrecomputedDerivative(stiffnessMatrix, loadVector, size, step, 0);
void addBoundaryConditionWithoutPrecomputedDerivative(vector<double>& matrix, vector<double>& loadVector, int& size, double step, int boundaryIndex) {
    // Увеличиваем размерность системы для включения производной как переменной
    int newSize = size + 1;
    vector<double> newMatrix(newSize * newSize, 0.0); // Создаем увеличенную матрицу жесткости
    vector<double> newLoadVector(newSize, 0.0);       // Создаем увеличенный вектор нагрузок

    // Копируем старую матрицу жесткости в увеличенную
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            newMatrix[i * newSize + j] = matrix[i * size + j];
        }
        newLoadVector[i] = loadVector[i]; // Копируем вектор нагрузок
    }

    // Добавляем новое уравнение для граничного условия второго рода
    // Связываем новую переменную (производную) с граничным узлом и следующим узлом
    newMatrix[boundaryIndex * newSize + newSize - 1] = -1.0;        // Влияние производной на граничный узел
    newMatrix[(newSize - 1) * newSize + boundaryIndex] = -1.0;      // Симметричный элемент
    newMatrix[(newSize - 1) * newSize + boundaryIndex + 1] = 1.0;   // Связь с соседним узлом
    newMatrix[(newSize - 1) * newSize + newSize - 1] = -step;       // Учет конечной разности (производная)

    // Вектор нагрузки для нового уравнения ставим в ноль (нулевой поток на границе)
    newLoadVector[newSize - 1] = 0.0;

    // Обновляем исходные matrix и loadVector, заменяя их увеличенными версиями
    matrix = move(newMatrix);
    loadVector = move(newLoadVector);
    size = newSize;

    cout << "Boundary condition added without precomputed derivative at index " << newSize - 1 << endl;
}



// Массив указателей на функции для применения граничных условий
int (*restrict[2]) (vector<double> &matrix, vector<double> &vector,
     int a, int len, int line,
     double value) = {restrictType1, restrictType2};

int main(int argc, char **argv) {
    // Обработка опций программы
    const options_struct settings = getOptions(argc, argv);
    if (settings.help)
        return 0;

    // Данные по варианту задачи
    const double a = 7, b = 0, c = -15, d = 60;
    restriction lower = {first, -7, 0}, upper = {second, 4, 4};
    // lower = {second, -7, realSolvePrime(-7)};
    cout << "Аналитическая производная в точке -7 = " << realSolvePrime(-7) << endl;
    // restriction lower = {second, -7, 0}, upper = {second, 4, 4};


    // Вычисление необходимых постоянных
    int size = settings.elemAmount * settings.type + 1; // Размер матрицы и векторов
    const double step = (upper.pos - lower.pos) / settings.elemAmount; // Шаг сетки

    // Используемые структуры данных
    vector<double> stiffnessMatrix(size * size, 0); // Матрица жесткости
    vector<double> loadVector(size, 0);             // Вектор нагрузок
    vector<double> displacements(size, 0);          // Вектор смещений (решение)
    vector<double> errors;                           // Вектор погрешностей
    vector<pair<double, double>> plotMKE;            // Данные для графика МКЭ
    vector<pair<double, double>> plotReal;           // Данные для графика точного решения
    vector<pair<double, double>> plotAbsolutelyReal;           // Данные для графика точного решения
    
    // Заполняем матрицу жесткости и вектор нагрузок в зависимости от типа элементов
    if (settings.type == linear)
        makeLinearSLAU(stiffnessMatrix, loadVector, size, step, a, b, c, d);
    else if (settings.type == cubic)
        makeCubicSLAU(stiffnessMatrix, loadVector, size, step, a, b, c, d);
    else {
        cout << "Произошла ошибка: тип уравнения не задан." << endl;
        return -1;
    }
    
    // Применяем граничные условия
    restrict[lower.grade](stiffnessMatrix, loadVector, a, size, 0, lower.val);
    // restrict[lower.grade](stiffnessMatrix, loadVector, a, size, 0, lower.val);
    addBoundaryConditionWithoutPrecomputedDerivative(stiffnessMatrix, loadVector, size, step, 0);
    // restrict[middle.grade](stiffnessMatrix, loadVector, a, size, 0.635 * size, middle.val);  // Доп. Задание: ограничение первого рода в точке 0
    // restrict[upper.grade](stiffnessMatrix, loadVector, a, size, size - 1, upper.val);

    // Решаем систему линейных алгебраических уравнений методом Гаусса
    solveSLAU(displacements, stiffnessMatrix, loadVector, size);

    // Вычисляем погрешности и заполняем массивы для вывода графиков
    double maxError = 0;
    double node = lower.pos;
    for (int i = 0; i < size; i++) {
        double realU = realSolve(node); // Точное решение в узле
        double error = fabs(displacements[i] - realU); // Погрешность
        if (error > maxError)
            maxError = error;
        errors.push_back(error);

        plotMKE.push_back(make_pair(node, displacements[i])); // Данные МКЭ
        plotReal.push_back(make_pair(node, realU));           // Данные точного решения

        // cout << "u(" << node << ") = " << displacements[i] << endl;
        node += step / settings.type; // Переход к следующему узлу
    }
    double x = lower.pos;   // Для дебага: рисую настоящий график решения.
    while (x < upper.pos) {
        double realU = realSolve(x); // Точное решение в узле
        plotAbsolutelyReal.push_back(make_pair(x, realU));           // Данные точного решения
        x += 0.1; // Переход к следующему узлу
    }

    // Выводим наибольшую погрешность
    cout << endl << "\033[4m" << "Ошибка не превышает " << maxError << "\033[0m" << "\n\n";

    // Экспортируем данные в файлы
    void export_data(
                    const int size, 
                    options_struct settings, 
                    vector<pair<double, double>> plotMKE,
                    vector<pair<double, double>> plotReal,
                    vector<pair<double, double>> plotAbsolutelyReal,
                    vector<double> displacements,
                    vector<double> errors
            );
    export_data(size, settings, plotMKE, plotReal, plotAbsolutelyReal, displacements, errors);
    exportStiffnessMatrixAndLoadVector(stiffnessMatrix, loadVector, size);

    return 0;
}

// Заполнение СЛАУ линейными конечными элементами
int makeLinearSLAU(vector<double> &resultMatrix,
                                     vector<double> &resultVector, int size, double L,
                                     double a, double b, double c, double d) {
    for (int i = 0; i < size - 1; ++i) {
        // Заполнение матрицы жесткости для линейных КЭ
        resultMatrix[i * size + i] += -a / L - b / 2. + c * L / 3.;
        resultMatrix[i * size + i + 1] += a / L + b / 2. + c * L / 6.;
        resultMatrix[(i + 1) * size + i] += a / L - b / 2. + c * L / 6.;
        resultMatrix[(i + 1) * size + i + 1] += -a / L + b / 2. + c * L / 3.;
        // Заполнение вектора нагрузок
        resultVector[i] += -d * L / 2.;
        resultVector[i + 1] += -d * L / 2.;
    }
    return 0;
}

// Заполнение СЛАУ кубическими конечными элементами
int makeCubicSLAU(vector<double> &resultMatrix,
                                    vector<double> &resultVector, int size, double L,
                                    double a, double b, double c, double d) {
    for (int i = 0; i < size - 3; i += 3) {
        // Заполнение матрицы жесткости для кубических КЭ
        resultMatrix[i * size + i] +=
                -a * 37 / (10 * L) - b / 2. + c * 8 * L / 105.;
        resultMatrix[i * size + i + 1] +=
                a * 189 / (40 * L) + b * 57 / 80. + c * 33 * L / 560.;
        resultMatrix[i * size + i + 2] +=
                -a * 27 / (20 * L) - b * 3 / 10. - c * 3 * L / 140.;
        resultMatrix[i * size + i + 3] +=
                a * 13 / (40 * L) + b * 7 / 80. + c * 19 * L / 1680.;

        resultMatrix[(i + 1) * size + i] +=
                a * 189 / (40 * L) - b * 57 / 80. + c * 33 * L / 560.;
        resultMatrix[(i + 1) * size + i + 1] +=
                -a * 54 / (5 * L) + c * 27 * L / 70.;
        resultMatrix[(i + 1) * size + i + 2] +=
                a * 297 / (40 * L) + b * 81 / 80. - c * 27 * L / 560.;
        resultMatrix[(i + 1) * size + i + 3] +=
                -a * 27 / (20 * L) - b * 3 / 10. - c * 3 * L / 140.;

        resultMatrix[(i + 2) * size + i] +=
                -a * 27 / (20 * L) + b * 3 / 10. - c * 3 * L / 140.;
        resultMatrix[(i + 2) * size + i + 1] +=
                a * 297 / (40 * L) - b * 81 / 80. - c * 27 * L / 560.;
        resultMatrix[(i + 2) * size + i + 2] +=
                -a * 54 / (5 * L) + c * 27 * L / 70.;
        resultMatrix[(i + 2) * size + i + 3] +=
                a * 189 / (40 * L) + b * 57 / 80. + c * 33 * L / 560.;

        resultMatrix[(i + 3) * size + i] +=
                a * 13 / (40 * L) - b * 7 / 80. + c * 19 * L / 1680.;
        resultMatrix[(i + 3) * size + i + 1] +=
                -a * 27 / (20 * L) + b * 3 / 10. - c * 3 * L / 140.;
        resultMatrix[(i + 3) * size + i + 2] +=
                a * 189 / (40 * L) - b * 57 / 80. + c * 33 * L / 560.;
        resultMatrix[(i + 3) * size + i + 3] +=
                -a * 37 / (10 * L) + b / 2. + c * 8 * L / 105.;

        // Заполнение вектора нагрузок для кубических КЭ
        resultVector[i] += -d * L / 8.;
        resultVector[i + 1] += -d * 3 * L / 8.;
        resultVector[i + 2] += -d * 3 * L / 8.;
        resultVector[i + 3] += -d * L / 8.;
    }
    return 0;
}

// Ограничение первого рода (задание значения функции в узле)
int restrictType1(vector<double> &matrix, vector<double> &vector,
                                    int a, int len, int line, double value) {
    // Обнуляем соответствующую строку матрицы жесткости
    for (auto iter{matrix.begin() + len * line}; iter != matrix.begin() + len * (line + 1); iter++) {
        *iter = 0;
    }
    // Устанавливаем диагональный элемент матрицы равным 1
    matrix[len * line + line] = 1;
    // Устанавливаем значение вектора нагрузок
    vector[line] = value;
    return 0;
}

// Ограничение второго рода (задание производной функции в узле)
int restrictType2(vector<double> &matrix, vector<double> &vector, int a, int len, int line, double value) {
    // Модифицируем вектор нагрузок в зависимости от положения узла
    vector[line] += a * (line ? -value : value);
    return 0;
}

// Простейший вывод матрицы (не используется, написан для дебага)
void printMatrix(vector<double> &matrix, int len) {
    int i = 0;
    cout << endl;
    for (auto node : matrix) {
        cout << node << '\t';
        i++;
        if (i == len) {
            i = 0;
            cout << endl;
        }
    }
    cout << endl;
}

// Решение СЛАУ методом Гаусса
int solveSLAU(vector<double> &retXVector, vector<double> inAMatrix,
                            vector<double> inBVector, int size) {
    // Прямой ход метода Гаусса
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            double d = inAMatrix[j * size + i] / inAMatrix[i * size + i];
            for (int k = i; k < size; k++)
                inAMatrix[j * size + k] -= d * inAMatrix[i * size + k];
            inBVector[j] -= d * inBVector[i];
        }
    }

    // Обратный ход метода Гаусса
    for (int i = size - 1; i >= 0; --i) {
        double rSum = 0;
        for (int j = i + 1; j < size; ++j)
            rSum += retXVector[j] * inAMatrix[i * size + j];
        retXVector[i] = (inBVector[i] - rSum) / inAMatrix[i * size + i];
    }
    return 0;
}

// Точное решение задачи
double realSolve(double x) {
    double C2 = -0.000141861924;
    double C1 = 0.007826998720;
    return C1 * exp(sqrt(15. / 7) * x) + C2 * exp(-sqrt(15. / 7) * x) + 4;
}

// Экспорт данных в файлы
void export_data(
                        const int size, 
                        options_struct settings, 
                        vector<pair<double, double>> plotMKE,
                        vector<pair<double, double>> plotReal,
                        vector<pair<double, double>> plotAbsolutelyReal,
                        vector<double> displacements,
                        vector<double> errors
            ) {
    fstream out;

    // Создаём файл с результатами рассчетов в формате таблицы LaTeX
    out.open("out.txt", ios::out);
    if (out.is_open()) {
        if (size < 21) {
            // Форматирование таблицы для небольшого количества элементов
            out << "\\begin{table}[H]\n\t\\centering\n\t\\begin{tabular}{|c|c|c|c|}"
                         "\n";
            out << "\t\t\\hline x \t& Аналитическое решение \t& Решение МКЭ \t& "
                         "Погрешность\t\\\\\n";
            for (int i = 0; i < size; i++)
                out << "\t\t\\hline " << plotReal[i].first << "\t & "
                        << plotReal[i].second << "\t & " << displacements[i] << "\t & "
                        << errors[i] << "\t\\\\\n";
            out << "\t\t\\hline\n\t\\end{tabular}\n\t\\caption{"
                    << settings.elemAmount << " КЭ}\n\t\\label{tab:"
                    << "lc"[settings.type / 2] << '_' << settings.elemAmount
                    << "}\n\\end{table}\n";
        } else {
            // Форматирование таблицы для большого количества элементов с использованием longtable
            out << "\\begin{longtable}[c]{|c|c|c|c|}\n";
            out << "\t\\hline x 	& Аналитическое решение 	& Решение МКЭ "
                         "	& Погрешность	\\\\ \\endhead\n";
            out << "\t\\caption{" << settings.elemAmount
                    << " КЭ}  \\endfoot\n\t\\caption{" << settings.elemAmount
                    << " КЭ\\label{tab:"
                    << "lc"[settings.type / 2] << '_' << settings.elemAmount
                    << "}} \\endlastfoot\n";
            for (int i = 0; i < size; i++)
                out << "\t\\hline " << plotReal[i].first << "\t & "
                        << plotReal[i].second << "\t & " << displacements[i] << "\t & "
                        << errors[i] << "\t\\\\\n";
            out << "\t\\hline\n\\end{longtable}\n";
        }
        out.close();
    } else {
        cerr << "Не удалось открыть выходной файл." << endl;
    }

    // Экспорт данных МКЭ в файл MKE.txt
    out.open("MKE.txt", ios::out);
    if (out.is_open()) {
        for (auto &vector : plotMKE) {
            out << vector.first << " " << vector.second << endl;
        }
        out.close();
    } else {
        cerr << "Не удалось открыть файл MKE.txt." << endl;
    }

    // Экспорт данных точного решения в файл Real.txt
    out.open("Real.txt", ios::out);
    if (out.is_open()) {
        for (auto &vector : plotReal) {
            out << vector.first << " " << vector.second << endl;
        }
        out.close();
    } else {
        cerr << "Не удалось открыть файл Real.txt." << endl;
    }
    out.close();

    // Экспорт данных точного решения в файл AbsolutelyReal.txt
    out.open("AbsolutelyReal.txt", ios::out);
    if (out.is_open()) {
        for (auto &vector : plotAbsolutelyReal) {
            out << vector.first << " " << vector.second << endl;
        }
        out.close();
    } else {
        cerr << "Не удалось открыть файл AbsolutelyReal.txt." << endl;
    }
    out.close();
}

// Функция отображения справки по использованию программы
void printUsage() {
    std::cout << std::endl;
    std::cout << "Использование:" << std::endl;
    std::cout << "./mke [-l | -c] [-s <SIZE>]" << std::endl;
    std::cout << std::endl;
    std::cout << "Опции:" << std::endl;
    std::cout << "\033[1m" << "\t-h, --help" << "\033[0m" << std::endl;
    std::cout << "\t\tпоказать справку." << std::endl;
    std::cout << std::endl;
    std::cout << "\033[1m" << "\t-s, --size" << "\033[0m" << std::endl;
    std::cout << "\t\tустановить количество элементов (по умолчанию 20)." << std::endl;
    std::cout << std::endl;
    std::cout << "\033[1m" << "\t-l, --linear" << "\033[0m" << std::endl;
    std::cout << "\t\tиспользовать линейные уравнения для элементов (используется по умолчанию)." << std::endl;
    std::cout << std::endl;
    std::cout << "\033[1m" << "\t-c, --cubic" << "\033[0m" << std::endl;
    std::cout << "\t\tиспользовать кубические уравнения для элементов." << std::endl;
    std::cout << std::endl;
}

// Обработка опций программы
options_struct getOptions(const int argc, char *argv[]) {
    void printUsage();
    options_struct options;
    options.type = unknown;
    options.elemAmount = 0;

    // Определение длинных опций для getopt_long
    const struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"linear", no_argument, NULL, 'l'},
        {"cubic", no_argument, NULL, 'c'},
        {"size", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}
    };
    int long_optind;
    char c;

    // Цикл обработки опций
    while ((c = getopt_long(argc, argv, "hlcs:", long_options, &long_optind)) != -1) {
        switch (c) {
        case 'h':
            printUsage(); // Вывод справки
            options.help = 1;
            return options;
        case 'l':
            options.type = linear; // Установка линейного типа элементов
            break;
        case 'c':
            options.type = cubic; // Установка кубического типа элементов
            break;
        case 's':
            options.elemAmount = atoi(optarg); // Установка количества элементов
            break;
        default:
            break;
        }
    }

    // Установка значения по умолчанию для количества элементов, если не указано
    if (!options.elemAmount) {
        cout << "Количество элементов - (" << "\033[1;4;37;41m" << 20 << "\033[0m" << ")" << endl;
        options.elemAmount = 20;
    }

    // Установка типа элементов по умолчанию, если не указано
    if (options.type == unknown) {
        cout << "Тип функции - (" << "\033[1;4;43m" << "linear" << "\033[0m" << ")" << endl;
        options.type = linear;
    }

    return options;
}
