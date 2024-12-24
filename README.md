# **Aproximación de datos vía funciones lineales continuas a trozos**

## Descripción
Este programa implementa diferentes algoritmos para encontrar una aproximación de una función continua Piecewise Linear (PWL) que minimice el error total de aproximación, dados ciertos puntos de datos. La aproximación se realiza colocando una cantidad determinada de breakpoints en una grilla de discretización y ajustando segmentos lineales entre ellos.

## Requisitos
### Python
- Python 3.x
- Numpy
- Matplotlib

### C++
- Compilador de C++
- nlohmann/json.hpp

## Instancias de Prueba
El programa viene con cinco instancias de prueba ubicadas en la carpeta data. Estas instancias contienen conjuntos de puntos de datos para los cuales se desea encontrar una aproximación PWL.

## Métodos Disponibles
- Fuerza Bruta
- Backtracking
- Programación Dinámica (Top-Down)
- Programación Dinámica (Bottom-Up)
- Programación Dinámica Refinada (Top-Down) - PYTHON: Implementa un enfoque de programación dinámica optimizado con un grafo para almacenar cálculos intermedios.

## Uso

1. Clona este repositorio o descarga los archivos del mismo.
2. Asegúrate de tener las dependencias instaladas.
3. Ejecuta el archivo main.py o compila y ejecuta el archivo main.cpp mediante la terminal, dependiendo su preferencia.
4. Sigue las instrucciones en la terminal para seleccionar la instancia, los parámetros necesarios y el método de optimización deseado. Introduce los valores solicitados cuando se te indique y presiona Enter para confirmar.

## Resultados
Una vez completada la ejecución, se imprimirán por pantalla los índices y coordenadas de los breakpoints de la función continua PWL y el error total de la aproximación. Además, se generará un archivo solution_instance.json que contiene la solución encontrada.
