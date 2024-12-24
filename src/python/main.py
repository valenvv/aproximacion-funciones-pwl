import json
import numpy as np
import matplotlib.pyplot as plt

import sys # Para finalizar la ejecución del programa en caso de que el usuario ingrese algo incorrecto
from typing import List, Tuple # Anotaciones de tipo

from Grafo import Grafo # Implementación extra

BIG_NUMBER = 1e10

def esFactible(pwl: List[Tuple[int, int]]) -> bool:
	'''
	Devuelve True si la función continua pwl pasada por parámetro es factible. 
	pwl es factible si el primer breakpoint se encuentra en la primera columna de la grilla
	de la discretización (i = 0) y el último breakpoint está en la última columna (i = m1)
	'''
	return pwl[0][0] == 0 and pwl[-1][0] == m1 - 1

def errorTotal(pwl: List[Tuple[int, int]]) -> float:
	'''
	Calcula el error de aproximación de una función continua pwl, el cual está dado por 
	la suma de los errores de cada una de sus piezas.
	'''
	res: float = 0
	for i in range(len(pwl) - 1):
		p1: Tuple[float,float] = grid_x[pwl[i][0]], grid_y[pwl[i][1]]
		p2: Tuple[float,float] = grid_x[pwl[i + 1][0]], grid_y[pwl[i + 1][1]]
		res += errorPorPieza(p1, p2)
	return res

def errorPorPieza(p1: Tuple[float, float], p2: Tuple[float, float]) -> float:
	'''
	Calcula el error de la aproximación de la pieza definida por los breakpoints p1 y p2, teniendo en
	cuenta los puntos comprendidos entre p1 (sin incluir) y p2 (incluido). En el caso de ser la primera
	pieza de una función, sí se incluye p1 a fin de considerar todos los puntos.
	'''
	res: float = 0
	i: int = 0

	if p1[0] == grid_x[0]: # primera pieza
		res += abs(instance['y'][0] - recta(p1, p2, instance['x'][0]))

	while i < instance['n'] and instance['x'][i] <= p2[0]:
		x: float = instance['x'][i]
		y: float = instance['y'][i]
		
		if p1[0] < x:
			res += abs(y - recta(p1, p2, x))
		
		i += 1

	return res

def recta(p1: Tuple[float, float], p2: Tuple[float, float], x: float) -> float:
	'''
	Calcula el valor de y en la recta definida por los puntos p1 y p2 para un valor de x dado.
	'''
	return (p2[1] - p1[1]) / (p2[0] - p1[0]) * (x - p1[0]) + p1[1]

def recorrerOrdenadas(abscisa: int, actual: List[Tuple[int, int]]) -> None:
	'''
	Explora todas las combinaciones posibles de breakpoints en la discretización, recorriendo 
	todas las ordenadas que se encuentran en la abscisa recibida por parámetro y utilizando backtracking.
	'''
	for ordenada in range(m2):
		actual.append((abscisa, ordenada))
		backtracking(actual)
		actual.pop()

def crearCubo(x: int, y: int, z: int) -> List[List[List[None]]]:
	'''
	Genera un cubo con dimensiones x, y y z e inicializa con el valor None en cada una de sus celdas.
	'''
	cubo: List[List[List[None]]] = [[[None for i in range(z)] for j in range(y)] for k in range(x)]
	return cubo

def regrLineal():
	x = np.array(instance["x"])
	y = np.array(instance["y"])

	# Calcular la pendiente y la ordenada al origen
	pendiente, ordOrigen = np.polyfit(x, y, 1)

	predicciones = pendiente * x + ordOrigen

	# Calcular el error 
	p1 = grid_x[0], predicciones[0]
	p2 = grid_x[-1], predicciones[-1]
	error = errorPorPieza(p1, p2)
	# Graficar los puntos de datos y la recta de predicciones
	plt.scatter(x, y, label='Datos originales')
	plt.plot(x, predicciones, color='red', label='Recta de predicciones')
	plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
	plt.title('Regresión Lineal')
	plt.show()

	return error

##### ALGORITMOS #####

def fuerzaBruta(actual: List[Tuple[int, int]]) -> None:
	if len(actual) == K:	#caso base
		if esFactible(actual) and errorTotal(actual) < best['obj']:
			best['sol'] = actual.copy()
			best['obj'] = errorTotal(actual)

	else:
		if len(actual) == 0:		# si no hay breakpoints, iteramos desde primer columna
			siguienteAbscisa = 0			
		else:						# para no generar breakpoints en la misma columna			
			siguienteAbscisa = actual[-1][0] + 1

		for abscisa in range(siguienteAbscisa, m1):	
			for ordenada in range(m2):
				actual.append((abscisa, ordenada))
				fuerzaBruta(actual)
				actual.pop()

def backtracking(actual: List[Tuple[int, int]]) -> None:
	if len(actual) == K:		# caso base
		if errorTotal(actual) < best['obj']:
			best['sol'] = actual.copy()
			best['obj'] = errorTotal(actual)

	elif errorTotal(actual) < best['obj']: # poda optimalidad
		if len(actual) == 0: # poda factibilidad
			recorrerOrdenadas(0, actual)
		elif len(actual) == K - 1: # poda factibilidad
			recorrerOrdenadas(m1 - 1, actual)
		else:
			siguienteAbscisa = actual[-1][0] + 1

			for abscisa in range(siguienteAbscisa, m1 - (K - len(actual)) + 1): # poda factibildad
				recorrerOrdenadas(abscisa, actual)

def progDinamicaTD(M: int, i: int, j: int) -> float: #Top down
	if C1[i][j][M] is not None:
		return C1[i][j][M]
	
	minimo = BIG_NUMBER
	ultBp = grid_x[i], grid_y[j]

	if M == 1:			# caso base
		for ordenada in range(m2):
			p1 = grid_x[0], grid_y[ordenada]
			error = errorPorPieza(p1, ultBp)
			minimo = min(minimo, error)

	else:
		for abscisa in range(M - 1, i): #considera espacio necesario para piezas restantes.
			for ordenada in range(m2):
				p1 = grid_x[abscisa], grid_y[ordenada]
				error = errorPorPieza(p1, ultBp) + progDinamicaTD(M - 1, abscisa, ordenada)
				minimo = min(minimo, error)

	C1[i][j][M] = minimo
	return C1[i][j][M]

def progDinamicaBU(M:int) -> None: #Bottom Up
	for abscisa in range(1, m1 - (M - 1)):  # recorre para M = 1 y considera espacio necesario para piezas restantes.
		for ordenada in range(m2):
			p2 = grid_x[abscisa], grid_y[ordenada]
			C2[abscisa][ordenada][1] = min([errorPorPieza((grid_x[0], grid_y[l]), p2) for l in range(m2)])
	
	for m in range(2, M + 1):	# recorre cant piezas
		for abscisa in range(m, m1 - (M - m)): #considera espacio necesario para piezas restantes.
			for ordenada in range(m2):
				p2 = grid_x[abscisa], grid_y[ordenada]
				minimo = BIG_NUMBER
				for k in range(m - 1, abscisa):	# recorre abscisas que esten antes de p2
					for l in range(m2):
						p1 = grid_x[k], grid_y[l]
						error = errorPorPieza(p1, p2) + C2[k][l][m - 1]
						minimo = min(minimo, error)
				C2[abscisa][ordenada][m] = minimo
	best['obj'] = min([C2[m1 - 1][o][M] for o in range(m2)])
		
def solucionOptima(refinado: bool) -> None: 
	res = BIG_NUMBER
	
	for j in range(m2):
		if not refinado: # si no estamos usando el algoritmo refinado
			res = min(res, progDinamicaTD(K - 1, m1 - 1, j))  #dev min entre valores de ult columna
		else: # algoritmo refinado con grafo (implementación extra)
			res = min(res, progDinamicaRefinada(K - 1, m1 - 1, j))

	best['obj'] = res

def reconstruirSolucion(C: List[List[List[float]]]) -> None:
	solucion = []

	M = K - 1  # Número de piezas
	i = m1 - 1
	j = np.argmin([C[i][k][M] for k in range(m2)])  # Buscar el indice del ult breakpoint
	solucion.append((i, j))

	while M > 1:
		encontramosBp = False
		ultBp = grid_x[i], grid_y[j]
		for x in range(i - 1, M - 2, -1):  # Retroceder en la matriz C para encontrar el próximo breakpoint
			# Buscamos en las celdas anteriores en la matriz C para encontrar la celda que se utilizó para calcular el valor de la celda actual. Esta celda lleva a la posición del próximo breakpoint.
			for y in range(m2): 
				p1 = grid_x[x], grid_y[y]
				if C[i][j][M] == errorPorPieza(p1, ultBp) + C[x][y][M - 1]: # Significa que el breakpoint (x, y) es parte de la solución
					i, j = x, y
					encontramosBp = True
					break
			if encontramosBp:
				solucion.append((i, j)) # Agregamos el breakpoint encontrado
				M -= 1 # Disminuimos la cant de piezas por agregar
				break

	# Cuando quede una pieza (M = 1), el breakpoint debe estar en la primera columna
	ultBp = grid_x[i], grid_y[j]
	for y in range(m2):
		p1 = grid_x[0], grid_y[y]
		if C[i][j][M] == errorPorPieza(p1, ultBp):
			i, j = 0, y
			solucion.append((i, j))
			break

	best['sol'] = solucion[::-1]  # Devolver la solución invertida para que esté en orden creciente de índices de abscisas

### IMPLEMENTACION EXTRA ###
def progDinamicaRefinada(M: int, i: int, j: int) -> float: # Top-down refinado
	if C3[i][j][M] is not None:
		return C3[i][j][M]

	minimo = BIG_NUMBER
	ultBp = grid_x[i], grid_y[j]

	if M == 1:			# caso base
		for ordenada in range(m2):
			if grafoErrores.existeArista((0, ordenada), (i, j)):
				errorPieza = grafoErrores.obtenerArista((0, ordenada), (i, j))
			else:
				p1 = grid_x[0], grid_y[ordenada]
				errorPieza = errorPorPieza(p1, ultBp)
				grafoErrores.agregarArista((0, ordenada), (i, j), errorPieza)
			
			minimo = min(minimo, errorPieza)

	else:
		for abscisa in range(M - 1, i): #considera espacio necesario para piezas restantes.
			for ordenada in range(m2):
				if grafoErrores.existeArista((abscisa, ordenada), (i, j)):
					errorPieza = grafoErrores.obtenerArista((abscisa, ordenada), (i, j))
				else:
					p1 = grid_x[abscisa], grid_y[ordenada]
					errorPieza = errorPorPieza(p1, ultBp)
					grafoErrores.agregarArista((abscisa, ordenada), (i, j), errorPieza)

				error = errorPieza + progDinamicaRefinada(M - 1, abscisa, ordenada)
				minimo = min(minimo, error)

	C3[i][j][M] = minimo
	return C3[i][j][M]

def reconstruirSolucionRefinado(C: List[List[List[float]]]) -> None:
	solucion = []

	M = K - 1  # Número de piezas
	i = m1 - 1
	j = np.argmin([C[i][k][M] for k in range(m2)])  # Buscar el indice del ult breakpoint
	solucion.append((i, j))

	while M > 1:
		encontramosBp = False
		for x in range(i - 1, M - 2, -1):  # Retroceder en la matriz C para encontrar el próximo breakpoint
			# Buscamos en las celdas anteriores en la matriz C para encontrar la celda que se utilizó para calcular el valor de la celda actual. Esta celda lleva a la posición del próximo breakpoint.
			for y in range(m2): 
				if C[i][j][M] == grafoErrores.obtenerArista((x, y), (i, j)) + C[x][y][M - 1]: # Significa que el breakpoint (x, y) es parte de la solución
					i, j = x, y
					encontramosBp = True
					break
			if encontramosBp:
				solucion.append((i, j)) # Agregamos el breakpoint encontrado
				M -= 1 # Disminuimos la cant de piezas por agregar
				break

	# Cuando quede una pieza (M = 1), el breakpoint debe estar en la primera columna
	for y in range(m2):
		if C[i][j][M] == grafoErrores.obtenerArista((0, y), (i, j)):
			i, j = 0, y
			solucion.append((i, j))
			break

	best['sol'] = solucion[::-1]  # Devolver la solución invertida para que esté en orden creciente de índices de abscisas

def main() -> None:

	instancias = ["aspen_simulation", "ethanol_water_vle", "optimistic_instance", "titanium", "toy_instance"]
	print("Elija una instancia:")
	for i in range(len(instancias)):
		print(f"{i + 1}. {instancias[i]}")

	global instance_name

	opcion1 = int(input("Ingrese el número de la instancia que desea ejecutar: "))

	if 1 > opcion1 or opcion1 > 5:
		print("Opción no válida. Por favor, elija un número del 1 al 5.")
		sys.exit()
	
	instance_name = instancias[opcion1 - 1] + ".json"
	filename = "../../data/" + instance_name
	with open(filename) as f:
		global instance
		instance = json.load(f)

	print("Elija cantidad de breakpoints K, m1 y m2:")

	global K
	global m1
	global m2
	global grid_x
	global grid_y
	global best
	
	K = int(input("Ingrese la cantidad de breakpoints (> 1): "))
	m1 = int(input("Ingrese m1: ")) # Cantidad de puntos en la discretización para las abscisas
	m2 = int(input("Ingrese m2: ")) # Cantidad de puntos en la discretización para las ordenadas

	# Ejemplo para definir una grilla de m x n.
	grid_x = np.linspace(min(instance["x"]), max(instance["x"]), num=m1, endpoint=True)
	grid_y = np.linspace(min(instance["y"]), max(instance["y"]), num=m2, endpoint=True)
 
	best= {} # Mejor solución encontrada hasta el momento.
	best['sol'] = [None]*(K) # Lista de tuplas con indices de los breakpoints
	best['obj'] = BIG_NUMBER # Error total de la mejor solución hasta el momento
	
	print("Elija un método:")
	print("1. Fuerza Bruta")
	print("2. Backtracking")
	print("3. Programacion Dinamica (Top-Down)")
	print("4. Programacion Dinamica (Bottom-Up)")
	print("5. Programacion Dinamica refinada (Top-Down)")

	opcion2 = input("Ingrese el número del metodo que desea ejecutar: ")

	if opcion2 == "1": # Fuerza Bruta
		fuerzaBruta(list())
	elif opcion2 == "2": # Backtracking
		backtracking(list())
	elif opcion2 == "3": # Top-Down
		global C1
		C1 = crearCubo(m1, m2, K - 1 + 1) # K es cantidad de breakpoints, K - 1 es cant de piezas
		solucionOptima(False)
		reconstruirSolucion(C1) 
	elif opcion2 == "4": # Bottom-Up
		global C2
		C2 = crearCubo(m1, m2, K - 1 + 1) # K es cantidad de breakpoints, K - 1 es cant de piezas
		progDinamicaBU(K - 1)
		reconstruirSolucion(C2)
	elif opcion2 == "5": # Implementación extra - Top-down refinado
		global C3
		global grafoErrores
		C3 = crearCubo(m1, m2, K - 1 + 1) # K es cantidad de breakpoints, K - 1 es cant de piezas
		grafoErrores = Grafo(m1, m2)
		solucionOptima(True) # True para que implemente el algoritmo refinado
		reconstruirSolucionRefinado(C3)
	else:
		print("Opción no válida. Por favor, elija un número del 1 al 5.")
		sys.exit()

	# Representamos la solución con un diccionario que indica:
	# - n: cantidad de breakpoints
	# - x: lista con las coordenadas de la abscisa para cada breakpoint
	# - y: lista con las coordenadas de la ordenada para cada breakpoint
	# - obj: error total de la solución
	solution = {}
	solution['n'] = len(best['sol'])
	solution['x'] = [grid_x[x[0]] for x in best['sol']]
	solution['y'] = [grid_y[x[1]] for x in best['sol']]
	solution['obj'] = best['obj']

	print("Índices de los breakpoints de la función continua PWL que minimizan el error total:")
	print(best['sol'])
	print("Coordenadas de los breakpoints de la función continua PWL que minimizan el error total:")
	print([(grid_x[bp[0]], grid_y[bp[1]]) for bp in best['sol']])
	print("Error total:")
	print(solution['obj'])

	# Se guarda el archivo en formato JSON
	with open('solution_' + instance_name, 'w') as f:
		json.dump(solution, f)

if __name__ == "__main__":
	main()