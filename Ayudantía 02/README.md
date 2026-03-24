# Ayudantía 02 — Análisis Variográfico: Acuífero Wolfcamp

## Descripción

En esta sesión se realiza un análisis geoestadístico completo sobre el dataset `wolfcamp` (paquete `geoR`), que contiene mediciones de carga piezométrica (presión del agua subterránea, en pies) en 85 pozos del acuífero Wolfcamp, Texas.

El foco está en el análisis del variograma empírico: su estimación, sensibilidad a parámetros, detección de anisotropía y ajuste de modelos teóricos, considerando la presencia de tendencia espacial de primer orden.

---

## Contenidos

| Parte | Descripción |
|-------|-------------|
| **a** | Visualización espacial y diagnóstico de estacionariedad de primer orden |
| **b** | Sensibilidad del variograma empírico al ancho de lag (`width`) |
| **c** | Variogramas direccionales y detección de anisotropía |
| **d** | Ajuste de modelos teóricos: Exponencial, Esférico y Gaussiano |
| **e** | Remoción de tendencia lineal y análisis de residuos |

---

## Conceptos clave

- **Estacionariedad de primer orden**: E[Z(s)] = μ constante en todo el dominio.
- **Estimador de Matheron**: γ̂(hₖ) = 1/(2|N(hₖ)|) · Σ [Z(sᵢ) − Z(sⱼ)]²
- **Parámetros del variograma**: nugget (c₀), sill (c₀ + c), range (a).
- **Modelos teóricos**:
  - Exponencial: γ(h) = c₀ + c·[1 − exp(−h/a)]
  - Esférico: γ(h) = c₀ + c·[3h/(2a) − h³/(2a³)] para h ≤ a
  - Gaussiano: γ(h) = c₀ + c·[1 − exp(−h²/a²)]
- **Modelo con tendencia**: Z(s) = μ(s) + ε(s), con μ(s) = β₀ + β₁X + β₂Y.

---

## Resultados principales

- Se detecta un gradiente espacial claro en dirección SW → NE, invalidando el supuesto de media constante (R² = 0.891 para la tendencia lineal).
- El variograma original crece sin sill hasta ~25.000 ft² — efecto directo de la tendencia.
- Los variogramas direccionales revelan anisotropía: mayor continuidad espacial en la dirección 135° (NW-SE), perpendicular al gradiente.
- Tras remover la tendencia, el variograma de residuos se estabiliza en ~4.000–4.500 ft² con sill definido, recuperando la estacionariedad de segundo orden.
- El modelo **Esférico** presenta el mejor ajuste (SSR ≈ 19.097) frente al Exponencial (SSR ≈ 29.924) y el Gaussiano (SSR ≈ 385.643).

---

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `Ay_02.pdf` | Enunciado de la ayudantía |
| `Ay_02_Espacial.R` | Script R comentado con el análisis completo |

---

## Librerías R utilizadas

```r
library(geoR)      # dataset wolfcamp y herramientas geoestadísticas
library(gstat)     # variogram(), fit.variogram(), vgm()
library(tidyverse) # manipulación de datos y visualización (ggplot2)
library(patchwork) # composición de múltiples gráficos
```
