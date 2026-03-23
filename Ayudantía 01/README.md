# Ayudantía 01 — Semivariograma Esférico

## Contenidos

- Introducción al semivariograma como herramienta de la geoestadística
- **Modelo esférico**: definición, fórmula y propiedades
- Parámetros del semivariograma:
  - **Nugget** ($\alpha$): discontinuidad en el origen, representa variabilidad a escala muy pequeña
  - **Partial Sill** ($\beta$): varianza estructural del proceso
  - **Sill** ($\alpha + \beta$): varianza total asintótica
  - **Range**: distancia a la cual el semivariograma alcanza el sill
- Visualización del semivariograma con `ggplot2`

## Modelo esférico

$$
\gamma(t) = \begin{cases}
0 & t = 0 \\
\alpha + \beta \left(\dfrac{3t}{2} - \dfrac{t^3}{2}\right) & 0 < t < 1 \\
\alpha + \beta & t \geq 1
\end{cases}
$$

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `Ay-01.pdf` | Enunciado de la ayudantía |
| `Ay_01_Sol.pdf` | Solución detallada |
| `Ay_01_Espacial_R.ipynb` | Notebook R: implementación y visualización del semivariograma esférico |

## Librerías R utilizadas

```r
library(ggplot2)  # visualización
library(dplyr)    # manipulación de datos
```
