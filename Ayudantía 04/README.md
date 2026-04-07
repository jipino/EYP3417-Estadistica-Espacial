# Ayudantía 04 — Kriging Simple y Ordinario

## Tema

Generación de mapas predictivos sobre el dominio del acuífero **Wolfcamp** mediante técnicas de Kriging, usando los parámetros θ = (σ², φ, τ²) estimados por máxima verosimilitud en la Ayudantía 03. Se comparan Kriging Simple (KS) y Kriging Ordinario (KO), se analizan sus propiedades teóricas y se discute el costo computacional.

---

## Contenidos

| Parte | Descripción |
|-------|-------------|
| **a** | Configuración del espacio de predicción: grilla regular 50×50 = 2.500 puntos sobre el dominio. Extracción del vector θ = (σ² = 3.782 ft², φ = 37.28 km, τ² = 672.75 ft²) desde `fit_ml`. |
| **b** | **Kriging Simple (KS)**: media constante y conocida E[Z(x)] = μ. Predicción con `krige.conv` (type.krige = "SK"). Verificación analítica de λ₀ = μ(1 − Σλᵢ) en el centroide geográfico. |
| **c** | **Kriging Ordinario (KO)**: media desconocida, restricción Σλᵢ = 1 via multiplicadores de Lagrange. Verificación del sistema aumentado (n+1)×(n+1). Comparación de mapas KS vs KO y mapa de diferencias. |
| **d** | **Análisis de varianza e interpolador exacto**: mapas de varianza de predicción para KS y KO. Demostración empírica de la propiedad de interpolador exacto con τ² = 0: error máximo = 0 ft y varianza máxima = 0 ft² sobre los 85 sitios. |
| **e** | **Costo computacional O(n³)**: fundamento teórico (factorización de Cholesky) y demostración empírica midiendo `solve()` para n ∈ {50, 85, 200, 500, 1000}. Extrapolación y motivación de métodos aproximados (NNGP, low-rank, SPDE/INLA, tapering). |

---

## Conceptos clave

- **Kriging Simple**: predictor Z*(x₀) = λ₀ + Σλᵢ Z(xᵢ) con λ₀ = μ(1 − Σλᵢ). En zonas sin datos, el predictor converge a la media conocida μ.
- **Kriging Ordinario**: mismo predictor con λ₀ = 0 y restricción Σλᵢ = 1. La media se estima implícitamente — no se requiere conocerla.
- **Sistema de Lagrange**: [Σ 1; 1ᵀ 0][λ; η] = [k; 1]. El multiplicador η puede ser negativo.
- **Varianza de kriging**: depende únicamente del diseño de muestreo y del modelo de covarianza, no de los valores observados. Mínima en datos, máxima en bordes.
- **Var(KO) ≥ Var(KS)**: KO paga el costo de no conocer μ con mayor incertidumbre de predicción.
- **Interpolador exacto**: con τ² = 0, Z*(xᵢ) = Z(xᵢ) y σ²ₖ(xᵢ) = 0 en todos los sitios observados. Con τ² > 0 el kriging es un smoother.
- **KS ≈ KO cuando μ_KS = β̂_GLS**: si la media "conocida" de KS es el estimador GLS óptimo (como devuelve `likfit`), ambos predictores son teóricamente equivalentes.
- **O(n³)**: costo de la factorización de Cholesky. Kriging exacto inviable para n > 10⁴.

---

## Resultados principales

**Vector θ estimado por MV (Ay_03):**

| Parámetro | Valor | Interpretación |
|-----------|-------|----------------|
| σ² | 3.782 ft² | Varianza del proceso (sill parcial) |
| φ | 37.28 km | Escala; rango práctico 3φ = 111.8 km |
| τ² | 672.75 ft² | Nugget (~15% del sill total) |
| μ (GLS) | 14.02 ft | Media conocida usada en KS |

**Verificación analítica en x₀ = (−30.2, −7.4) km (centroide):**

| Métrica | KS | KO |
|---------|----|----|
| Σλᵢ | 0.9399 | 1.0000 (exacto) |
| λ₀ | 0.843 ft | 0 |
| η (Lagrange) | — | −15.77 |
| Predictor | 51.16 ft | 51.16 ft |
| ECM | — | 2.511.97 ft² |

**Interpolador exacto (τ² = 0, KO, 85 sitios):**

| Métrica | Valor |
|---------|-------|
| Error máximo | 0 ft |
| Varianza máxima | 0 ft² |

**Costo computacional:**

| n | Tiempo solve() | Factor vs n/2 |
|---|----------------|---------------|
| 500 | ~0.03 s | — |
| 1.000 | ~0.20 s | ×6.7 (teórico: ×8) |
| 10.000 (extrapolado) | ~200 s | — |
| 50.000 (extrapolado) | ~3 días | — |

---

## Librerías R

```r
library(geoR)      # krige.conv, krige.control, likfit
library(gstat)     # vgm, variogramLine
library(tidyverse) # dplyr, ggplot2, tibble
library(patchwork) # composición de gráficos
```

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `Ay_04.pdf` | Enunciado de la ayudantía |
| `Ay_04_Espacial.R` | Script R comentado con toda la sesión |
