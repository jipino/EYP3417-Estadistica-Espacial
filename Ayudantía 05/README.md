# Ayudantía 05 — Pipeline Geoestadístico Completo: Isla Rongelap

## Tema

Recorrido del pipeline geoestadístico completo sobre el dataset **rongelap** (paquete geoR): 157 mediciones de radiación gamma en la isla Rongelap (Pacífico), contaminada por pruebas nucleares estadounidenses en los años 1950. La variable de interés es la tasa de radiación r(xᵢ) = data / units.m (conteos por segundo). Se integran todos los contenidos de las ayudantías 01–04 en un análisis unificado, más una pregunta teórica sobre las propiedades del Kriging Ordinario.

---

## Contenidos

| Parte | Descripción |
|-------|-------------|
| **1a** | **Derivación del ECM de KO**: demostración de σ²_KO(x₀) = C(0) − λᵀk + η a partir del sistema aumentado de Lagrange. |
| **1b** | **Interpretación de η y desigualdad σ²_KO ≥ σ²_KS**: el multiplicador η como "costo" de no conocer μ; argumento formal de la desigualdad. |
| **1c** | **Condiciones geográficas para |η| grande**: zonas alejadas de los datos vs zonas bien rodeadas; conexión con el comportamiento de KS cuando μ es incorrecta. |
| **2a** | **EDA y tendencia espacial**: scatter coloreado por r(xᵢ), modelo lm(rate ~ x + y), R² = 0.149, decisión de trabajar con residuos. |
| **2b** | **Variograma experimental y ajuste MCO/MCP**: variograma omnidireccional, sensibilidad al width y cutoff, variogramas direccionales, ajuste Exponencial por MCO y MCP (cutoff = 1.500 m). |
| **2c** | **Estimación por MV y selección de modelo**: likfit() Exponencial vs Esférico, comparación por AIC, reporte de θ̂ = (σ², φ, τ²) y comparación MCO/MCP/MV. |
| **2d** | **Mapas de predicción KS y KO**: grilla 100×100, kriging con media aritmética (μ = 0) y KO, mapas de varianza de predicción para ambos métodos. |
| **2e** | **Comparación de predictores**: mapa de diferencias Z*_KO − Z*_KS; análisis de dónde y por qué divergen; discusión del caso μ_KS = μ_GLS. |
| **2f** | **Validación cruzada LOO**: xvalid() para KO, cálculo de RMSE y MSSE, gráfico de diagnóstico errores vs predichos. |

---

## Conceptos clave

- **Efecto hoyo**: variograma no monótono que sube hasta h ≈ 882 m y luego desciende. Refleja la geometría irregular de la isla y la heterogeneidad espacial compleja de la contaminación.
- **Anisotropía fuerte**: los variogramas direccionales muestran sills muy distintos (2.7 vs 6.4 ct²/s²) según la dirección. La dirección E-W es degenerada por la falta de pares.
- **Nugget dominante**: τ²/(σ²+τ²) ≈ 72% del sill total. Alta variabilidad a escala sub-muestral; el kriging nunca puede ser muy preciso.
- **Rango corto**: φ = 249 m (Esférico). A más de 250 m de cualquier sitio, la correlación es 0 y la predicción colapsa a la media.
- **μ_KS = 0 ≠ μ_GLS = −0.248**: la diferencia entre la media aritmética de los residuos OLS (exactamente 0) y la estimada por GLS genera el mapa de diferencias KO − KS constante e igual a −0.248 en el mar.
- **MSSE ≈ 1**: modelo bien calibrado — la varianza de kriging predice correctamente la magnitud del error de predicción.

---

## Resultados principales

**Modelo de tendencia lineal:**

| Coeficiente | Estimación | p-valor | Interpretación |
|-------------|-----------|---------|----------------|
| β₀ | 6.40 ct/s | < 0.001 | intercepto |
| β₁ (x) | −4.5×10⁻⁴ | 0.071 | leve gradiente E-W, marginal |
| β₂ (y) | 2.8×10⁻⁵ | 0.967 | sin gradiente N-S |
| R² | 0.149 | — | solo 15% de variabilidad explicada |

**Parámetros MCO / MCP / MV (modelo Exponencial o Esférico según columna):**

| Método | Modelo | σ² | φ (m) | τ² |
|--------|--------|----|-------|----|
| MCO | Exponencial | ~4.79 | ~178 | ~2.65 |
| MCP | Exponencial | ~4.37 | ~250 | ~3.30 |
| **MV** | **Esférico** | **1.83** | **249** | **4.72** |

**Selección de modelo por AIC:**

| Modelo | AIC |
|--------|-----|
| Exponencial | 738.86 |
| **Esférico** | **738.42** ← seleccionado |

**Mapas de predicción (residuos, ct/s):**

| Métrica | KS | KO |
|---------|----|----|
| Predicción mín. | −2.86 | −3.01 |
| Predicción máx. | +1.41 | +1.38 |
| Varianza mín. | 5.38 | 5.38 |
| Varianza máx. | 6.55 | 6.63 |
| Diferencia KO − KS | — | [−0.248, −0.014] |

**Validación cruzada LOO (KO):**

| Métrica | Valor | Interpretación |
|---------|-------|----------------|
| RMSE | 2.49 ct/s | error de predicción promedio |
| MSSE | 1.019 | ≈ 1 → modelo bien calibrado |

---

## Librerías R

```r
library(geoR)      # variog, variog4, variofit, likfit, krige.conv, xvalid
library(gstat)     # (cargada por compatibilidad)
library(tidyverse) # dplyr, ggplot2, tibble
library(patchwork) # composición de gráficos
```

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `Ay_05.pdf` | Enunciado de la ayudantía |
| `Ay_05_Espacial.R` | Script R comentado con toda la sesión |
