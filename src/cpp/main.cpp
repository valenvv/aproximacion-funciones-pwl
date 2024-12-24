#include <string>
#include <iostream>
#include <fstream>
#include "include/json.hpp"
#include <vector>
#include <utility> // pair
#include <cmath>
#include <algorithm>

// Para libreria de JSON.
using namespace nlohmann;
using namespace std;

// variables globales
const double BIG_NUMBER = 1e10;
string instance_name;
json instance;
int K;
int m1;
int m2;
vector<pair<int, int>> bestSol; // Error total de la mejor solución hasta el momento
double bestObj = BIG_NUMBER;
vector<double> grid_x;
vector<double> grid_y;
vector<vector<vector<double>>> C1;
vector<vector<vector<double>>> C2;

// auxiliares
bool esFactible(const vector<pair<int, int>>& pwl, int m1);
double errorTotal(const vector<pair<int,int>> & pwl);
double errorPorPieza(const pair<double, double>& p1, const pair<double, double>& p2);
double recta(const pair<double, double>& p1, const pair<double, double>& p2, double x);
vector<vector<vector<double>>> crearCubo(int x, int y, int z);
vector<double> discretizacion(double v1, double vn, int m);
void recorrerOrdenadas(int abscisa, vector<pair<int, int>>& actual);

// algoritmos
void fuerzaBruta(vector<pair<int, int>>& actual);
void backtracking(std::vector<std::pair<int, int>>& actual);
double progDinamicaTD(int M, int i, int j);
void progDinamicaBU(int M);

// auxiliares de dinámica
void reconstruirSolucion(vector<vector<vector<double>>>& C);
void solucionOptima();

// funciones auxiliares
bool esFactible(const vector<pair<int, int>>& pwl) {
    /*
    Devuelve true si la función continua pwl pasada por parámetro es factible.
    pwl es factible si el primer breakpoint se encuentra en la primera columna de la grilla
    de la discretización (i = 0) y el último breakpoint está en la última columna (i = m1 - 1).
    */
    return pwl[0].first == 0 && pwl.back().first == m1 - 1;
}

// grid 
vector<double> discretizacion(double v1, double vn, int m) {
    vector<double> res;
    double delta = (vn - v1) / (m - 1);
    for (int i = 0; i < m; i++) {
        res.push_back(v1 + i * delta);
    }
    return res;
}

double errorTotal(const vector<pair<int,int>> & pwl) {
    /*
    Calcula el error de aproximación de una función continua pwl, el cual está dado por 
	la suma de los errores de cada una de sus piezas.
    */

   double res = 0;
   if (!pwl.empty()) { // Verificar si el vector pwl está vacío, en la primera ejecucion de backtracking lo está
       for (int i = 0; i < pwl.size() - 1; i++) {
            pair<double, double> p1 = make_pair(grid_x[pwl[i].first], grid_y[pwl[i].second]);
            pair<double, double> p2 = make_pair(grid_x[pwl[i + 1].first], grid_y[pwl[i + 1].second]);
            res += errorPorPieza(p1, p2);
        }
    }
    return res;
}

double errorPorPieza(const pair<double, double>& p1, const pair<double, double>& p2) {
    /*
    Calcula el error de la aproximación de la pieza definida por los breakpoints p1 y p2, teniendo en
	cuenta los puntos comprendidos entre p1 (sin incluir) y p2 (incluido). En el caso de ser la primera
	pieza de una función, sí se incluye p1 a fin de considerar todos los puntos.
    */
   double res = 0;
   int i = 0;

    if (p1.first == grid_x[0]) {
        res += abs(instance["y"][0].get<double>() - recta(p1, p2, instance["x"][0].get<double>()));
    }

    while (i < instance["n"] && instance["x"][i] <= p2.first) {
        double x = instance["x"][i];
        double y = instance["y"][i];

        if (p1.first < x) {
            res += abs(y - recta(p1, p2, x));
        }

        i++;
    }

    return res;
}


double recta(const pair<double, double>& p1, const pair<double, double>& p2, double x) {
    /*
    Calcula el valor de y en la recta definida por los puntos p1 y p2 para un valor de x dado.
    */
    return (p2.second - p1.second) / (p2.first - p1.first) * (x - p1.first) + p1.second;
}

void recorrerOrdenadas(int abscisa, vector<pair<int, int>>& actual) {
    /*
    Explora todas las combinaciones posibles de breakpoints en la discretización,
    recorriendo todas las ordenadas que se encuentran en la abscisa recibida por parámetro
    y utilizando backtracking.
    */
    for (int ordenada = 0; ordenada < m2; ++ordenada) {
        actual.push_back(make_pair(abscisa, ordenada));
        backtracking(actual);
        actual.pop_back();
    }
}

// algoritmos pwl
void fuerzaBruta(vector<pair<int, int>> &actual) {

    if (actual.size() == K) {
        if (esFactible(actual) && errorTotal(actual) < bestObj) {
            bestSol = actual;
            bestObj = errorTotal(actual);
        }
    } else {
        int siguienteAbscisa;
        if (actual.empty()) {
            siguienteAbscisa = 0;
        } else {
            siguienteAbscisa = actual.back().first + 1;
        }

        for (int abscisa = siguienteAbscisa; abscisa < m1; ++abscisa) {
            for (int ordenada = 0; ordenada < m2; ++ordenada) {
                actual.push_back(make_pair(abscisa, ordenada));
                fuerzaBruta(actual);
                actual.pop_back();

            }
        }
    }
}

void backtracking(std::vector<std::pair<int, int>>& actual) {
    if (actual.size() == K) { // caso base
        if (esFactible(actual) && errorTotal(actual) < bestObj) {
            bestSol = actual; // Guardar la mejor solución hasta el momento
            bestObj = errorTotal(actual);
        }
        
    } else if (errorTotal(actual) < bestObj) { // poda optimalidad
        int siguienteAbscisa;
        if (actual.empty()) { // poda factibilidad
            recorrerOrdenadas(0, actual);
        } else if (actual.size() == K - 1) { // poda factibilidad
            recorrerOrdenadas(m1 - 1, actual);
        } else {
            siguienteAbscisa = actual.back().first + 1;

            for (int abscisa = siguienteAbscisa; abscisa < m1 - (K - actual.size()) + 1; ++abscisa) { // poda factibilidad
                recorrerOrdenadas(abscisa, actual);
            }
        }
    }
}

vector<vector<vector<double>>> crearCubo(int x, int y, int z) {
    vector<vector<vector<double>>> cubo;
    for (int k = 0; k < x; ++k) { //dimension x del cubo
        vector<vector<double>> plano_y; // representa una fila en la dimension y del cubo

        for (int j = 0; j < y; ++j) {
            vector<double> fila_z; // representa una columna en la dimensión z del cubo
            
            for (int i = 0; i < z; ++i) {
                fila_z.push_back(BIG_NUMBER);
            }
            plano_y.push_back(fila_z);
        }
        cubo.push_back(plano_y);
    }
        return cubo;
}

double progDinamicaTD(int M, int i, int j) {
     
    if (C1[i][j][M] != BIG_NUMBER) {
        return C1[i][j][M];
    }
    
    double minimo = BIG_NUMBER;
    pair <double, double> ultBp = make_pair(grid_x[i], grid_y[j]); 

    if (M == 1) { // Caso base
        for (int ordenada = 0; ordenada < m2; ++ordenada) {
            pair<double, double> p1 = make_pair(grid_x[0], grid_y[ordenada]);
            double error = errorPorPieza(p1, ultBp);
            minimo = min(minimo, error);
        }
    } else {
        for (int abscisa = M - 1; abscisa < i; ++abscisa) { // Considera espacio necesario para piezas restantes
            for (int ordenada = 0; ordenada < m2; ++ordenada) {
                pair<double, double> p1= make_pair(grid_x[abscisa], grid_y[ordenada]);
                double error = errorPorPieza(p1, ultBp) + progDinamicaTD(M - 1, abscisa, ordenada);
                minimo = min(minimo, error);
            }
        }
    }

    C1[i][j][M] = minimo; 
    return C1[i][j][M];
}

void solucionOptima() {
    double res = BIG_NUMBER;
    for (int j = 0; j < m2; ++j) {
        res = min(res, progDinamicaTD(K - 1, m1 - 1, j)); // devuelve mínimo entre valores de la última columna
    }
   bestObj = res;
}

void progDinamicaBU(int M) {

    for (int abscisa = 1; abscisa < m1 - (M - 1); ++abscisa) {
        for (int ordenada = 0; ordenada < m2; ++ordenada) {
            pair<double, double> p2 = make_pair(grid_x[abscisa], grid_y[ordenada]);
            double minimo = BIG_NUMBER;
            for (int l = 0; l < m2; ++l) {
               pair<double, double> p1 = make_pair(grid_x[0], grid_y[l]);
                double error = errorPorPieza(p1, p2);
                minimo = min(minimo, error);
            }
            C2[abscisa][ordenada][1] = minimo;
        }
    }

    for (int m = 2; m < M+1 ; ++m) { //recorre cant piezas
        for (int abscisa = m; abscisa < m1 - (M - m); ++abscisa) { //considera espacio necesario para piezas restantes.
            for (int ordenada = 0; ordenada < m2; ++ordenada) {
                pair<double, double> p2 = make_pair(grid_x[abscisa], grid_y[ordenada]);
                double minimo = BIG_NUMBER;
                for (int k = m - 1; k < abscisa; ++k) { //recorre abscisas que esten antes de p2
                    for (int l = 0; l < m2; ++l) {
                        pair<double, double> p1 = make_pair(grid_x[k], grid_y[l]);
                        double error = errorPorPieza(p1, p2) + C2[k][l][m - 1];
                        minimo = min(minimo, error);
                    }
                }
                C2[abscisa][ordenada][m] = minimo;
            }
        }
    }

    double min_obj = BIG_NUMBER;
    for (int o = 0; o < m2; ++o) {
        min_obj = min(min_obj, C2[m1 - 1][o][M]);
    }
  bestObj = min_obj;
}

// auxiliares de dinámica
void reconstruirSolucion(vector<vector<vector<double>>>& C) {
    const double EPSILON = 1e-6;
    vector<pair<int, int>> solucion; 

    int M = K - 1;  // Número de piezas
    int i = m1 - 1;
    int j = 0;
    double minimo =BIG_NUMBER;
    for (int k = 0; k < m2; k++) { // Buscar el indice del ult breakpoint
        if (C[i][k][M] < minimo) {
            minimo = C[i][k][M];
            j = k;
        }
    }
    solucion.push_back(make_pair(i, j));

    while (M > 1) {
        bool encontramosBp = false;
        pair<double, double> ultBp = make_pair(grid_x[i], grid_y[j]);

        for (int x = i - 1; x >= M - 1; x--) { // Retroceder en la matriz C para encontrar el próximo breakpoint
            // Buscamos en las celdas anteriores en la matriz C para encontrar la celda que se utilizó para calcular el valor de la celda actual.
            // Esta celda lleva a la posición del próximo breakpoint.
            for (int y = 0; y < m2; y++) {
                pair<double, double> p1 = make_pair(grid_x[x], grid_y[y]);
                if (fabs(C[i][j][M] - (errorPorPieza(p1, ultBp) + C[x][y][M - 1])) < EPSILON) { // Significa que el breakpoint (x, y) es parte de la solución
                    i = x;
                    j = y;
                    encontramosBp = true;
                    break;
                }
            }
            if (encontramosBp) {
                solucion.push_back(make_pair(i, j)); // Agregamos el breakpoint encontrado
                M--; // Disminuimos la cant de piezas por agregar
                break;
            }
        }
    }
    
    // Cuando quede una pieza (M = 1), el breakpoint debe estar en la primera columna
    pair<double, double> ultBp = make_pair(grid_x[i], grid_y[j]);
    for (int y = 0; y < m2; y++) {
        pair<double, double> p1 = make_pair(grid_x[0], grid_y[y]);
        if (fabs(C[i][j][M] - errorPorPieza(p1, ultBp)) < EPSILON) {
            solucion.push_back(make_pair(0, y));
            break;
        }
    }

    reverse(solucion.begin(), solucion.end()); // Devolver la solución invertida para que esté en orden creciente de índices de abscisas
    bestSol = solucion;
}

int main() {
    vector<string> instancias = {"aspen_simulation", "ethanol_water_vle", "optimistic_instance", "titanium", "toy_instance"};

    cout << "Elija una instancia:" << endl;
    for (int i = 0; i < instancias.size(); i++) {
        cout << i + 1 << ". "<< instancias[i] << endl;
    }

    int opcion1;
    cout << "Ingrese el número de la instancia que desea ejecutar: ";
    cin >> opcion1;

    if (opcion1 < 1 || opcion1 > 5) {
        cout << "Opción no válida. Por favor, elija un número del 1 al 5." << endl;
        exit(0);
    }

    instance_name = instancias[opcion1 -1] + ".json";
    string file_name = "../../data/" + instance_name;
    ifstream input(file_name);

    input >> instance;
    input.close();

    cout << "Elija cantidad de breakpoints K, m1 y m2:" << endl;

    cout << "Ingrese la cantidad de breakpoints (> 1): ";
    cin >> K;
    cout << "Ingrese m1: ";
    cin >> m1;
    cout << "Ingrese m2: ";
    cin >> m2;

    grid_x = discretizacion(*min_element(instance["x"].begin(),instance["x"].end()), *max_element(instance["x"].begin(),instance["x"].end()), m1);
    grid_y = discretizacion(*min_element(instance["y"].begin(),instance["y"].end()), *max_element(instance["y"].begin(),instance["y"].end()), m2);
    vector<pair<int, int>> actual;
    
    vector<string> metodos = {"Fuerza Bruta", "Backtracking", "Programacion Dinamica (Top-Down)", "Programacion Dinamica (Bottom-Up)"};

    cout << "Elija un método:" << endl;
    for (int i = 0; i < metodos.size(); i++) {
        cout << i + 1 << ". "<< metodos[i] << endl;
    }
    int opcion2;
    cout << "Ingrese el número del método que desea ejecutar: ";
    cin >> opcion2;

    if (opcion2 == 1) { // fuerza bruta
        fuerzaBruta(actual);
    } else if (opcion2 == 2) { // backtracking
        backtracking(actual);
    } else if (opcion2 == 3) { // top-down
        C1 = crearCubo(m1, m2, K - 1 + 1);
        solucionOptima();
        reconstruirSolucion(C1);
    } else if (opcion2 == 4) { // bottom-up
        C2 = crearCubo(m1, m2, K - 1 + 1);
        progDinamicaBU(K-1);
        reconstruirSolucion(C2);
    } else {
        cout << "Opción no válida. Por favor, elija un número del 1 al 4." << endl;
        exit(0);
    }

    cout << "Índices de los breakpoints de la función continua PWL que minimizan el error total:" << endl;
    for (const auto& p : bestSol) {
        cout << "(" << p.first << ", " << p.second << ") " ;
    }
    cout << endl << "Coordenadas de los breakpoints de la función continua PWL que minimizan el error total:" << endl;
    for (const auto& bp : bestSol) {
        double x = grid_x[bp.first];
        double y = grid_y[bp.second];

        cout << "(" << x << ", " << y << ") ";
    }
    cout << endl << "Mejor error total encontrado: " << bestObj << endl;

    // Representamos la solución con un diccionario que indica:
    // - n: cantidad de breakpoints
    // - x: lista con las coordenadas de la abscisa para cada breakpoint
    // - y: lista con las coordenadas de la ordenada para cada breakpoint
    // - obj: error total de la solución
	json solution;
	solution["n"] = bestSol.size();
    solution["x"] = json::array();
	solution["y"] = json::array();
	solution["obj"] = bestObj;

    for (const auto& bp : bestSol) {
        solution["x"].push_back(grid_x[bp.first]);
        solution["y"].push_back(grid_y[bp.second]);
    }

    string archivo_solucion = "solution_" + instance_name;
    ofstream output(archivo_solucion);
    output << solution;
    output.close();

    return 0;
}