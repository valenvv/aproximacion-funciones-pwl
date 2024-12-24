from typing import Tuple

class Grafo:
    def __init__(self, m1: int, m2: int) -> None:
        self._m1 = m1
        self._m2 = m2
        n = m1 * m2
        self._matriz = [[-1 for i in range(n)] for j in range(n)] # matriz adyacencia

    def agregarArista(self, p1: Tuple[int, int], p2: Tuple[int, int], peso: float) -> None:
        x, y = self.obtenerIndices(p1, p2)
        self._matriz[x][y] = self._matriz[y][x] = peso
        
    def eliminarArista(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> None:
        x, y = self.obtenerIndices(p1, p2)
        self._matriz[x][y] = self._matriz[y][x] = -1

    def existeArista(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> bool:
        x, y = self.obtenerIndices(p1, p2)
        return self._matriz[x][y] != -1 # si no es -1, ya fue calculada
    
    def obtenerArista(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> float:
        x, y = self.obtenerIndices(p1, p2)
        return self._matriz[x][y]
    
    def tamanio(self) -> int: # cantidad de nodos
        return self._m1 * self._m2
    
    def obtenerIndices(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> Tuple[int, int]:
        x = p1[0] * self._m2 + p1[1] # indice de p1 en la matriz
        y = p2[0] * self._m2 + p2[1] # indice de p2 en la matriz
        return x, y
    
# no agregamos ni eliminamos vértices porque la grilla no cambia
# -1 significa que no tiene arista