# Task02, only python version

I am sorry, I will put here the task 02 only in python because I had not enought time to complete it even on C++

# Problem 1
```python
import numpy as np

# Definizione delle dimensioni del vettore
N = 10  # Per la task dopo uso anche 10**6 e 10**8
a = 3
# Creazione di un vettore di dimensione N con tutti gli elementi uguali a 0.1
x = np.full(N, 0.1)
# Creazione di un vettore di dimensione N con tutti gli elementi uguali a 7.1
y = np.full(N, 7.1)
# Calcolo del vettore risultante dalla somma ponderata di x e y moltiplicati per a
d = a * x + y
print (d)
# Stampa della lunghezza di d per verificare che sia corretta
print(len(d))  
# Verifica che tutti gli elementi di d siano vicini a 7.4
print(np.allclose(d, [7.4] * len(d)))

```
# Problem 2

```python
import numpy as np

# Definisco la matrice
for N in [10]:
    # Creazione di una matrice NxN con tutti gli elementi uguali a 3
    A = np.full((N, N), 3.0)
    # Creazione di una matrice NxN con tutti gli elementi uguali a 7.1
    B = np.full((N, N), 7.1)
    # Moltiplicazione delle matrici
    C = np.dot(A, B)
    
    # Verifica che tutti gli elementi di C siano 21.3 (questo controllo è teorico e mostrerà un errore nel contesto reale poiché il risultato atteso è 21.3*N)
    check = np.allclose(C, 21.3 * np.ones((N, N)))

    # Stampa della matrice risultante e del risultato della verifica
    print(C)
    print(f"Dimension {N}x{N}: All elements are 21.3:", check)
    # Nota: non può essere uguale a 21.3 in quanto prodotto tra matrici, sarà uguale a 21.3*N
