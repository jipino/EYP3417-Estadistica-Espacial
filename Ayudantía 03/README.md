# Ayudantía 03 — Inferencia Paramétrica

## Tema

Estimación de parámetros de covarianza en campos aleatorios gaussianos. Se trabaja sobre el dataset **Wolfcamp** (acuífero, 85 pozos) extendiendo el análisis variográfico de la Ayudantía 02 hacia la estimación formal por máxima verosimilitud y la detección de anisotropía geométrica.

---

## Contenidos

| Parte | Descripción |
|-------|-------------|
| **a** | Modelo de tendencia lineal, extracción de residuos y variograma experimental de los residuos. Justificación de la necesidad de estacionariedad para estimar θ = (σ², φ, τ²). |
| **b** | Ajuste del variograma por **MCO** (fit.method = 6) y **MCP** (fit.method = 7). Diferencias en los pesos, impacto en el nugget y el range estimado. |
| **c** | Estimación por **Máxima Verosimilitud** con `likfit` (geoR). Log-verosimilitud gaussiana, ventajas de MV sobre mínimos cuadrados, comparación de los tres métodos (MCO / MCP / MV). |
| **d** | Familia **Matérn**: comparación entre ν = 0.5 (Exponencial) y ν = 1.5 vía AIC. Interpretación del parámetro de suavidad ν y diferenciabilidad media cuadrática. |
| **e** | **Anisotropía geométrica** vía verosimilitud perfilada. Búsqueda en grilla sobre ángulos {0°, 45°, 90°, 135°} y razones {0.3, 0.5, 0.7, 0.9}. Transformación manual de coordenadas como workaround al bug de geoR en R 4.x. |

---

## Conceptos clave

- **Estacionariedad de segundo orden**: necesaria para que γ(h) dependa solo de h y no de la ubicación s.
- **MCO vs MCP**: los pesos wₖ determinan qué parte del variograma empírico se ajusta con mayor precisión. MCP pondera por Nₖ / [γθ(hₖ)]², favoreciendo los lags cortos.
- **MV en geostatística**: θ̂ = argmax ℓ(θ|Z), con ℓ la log-verosimilitud gaussiana. Costo O(n³), asintóticamente eficiente bajo gaussianidad.
- **AIC = −2·loglik + 2·p**: criterio de selección de modelos que penaliza la complejidad. ΔAIC < 2 → soporte empírico similar.
- **Familia Matérn**: ν controla la diferenciabilidad del campo. ν = 0.5 ↔ Exponencial (continuo, no diferenciable); ν = 1.5 → una vez diferenciable.
- **Anisotropía geométrica**: la correlación decae a distinto ritmo según la dirección. Se parametriza con psi.A (ángulo del eje de mayor correlación) y psi.R (razón rango_menor / rango_mayor).

---

## Resultados principales

| Método | Nugget (ft²) | PSill (ft²) | Rango práctico (km) |
|--------|-------------|-------------|---------------------|
| MCO    | 722         | 4034        | ~59                 |
| MCP    | 1228        | 5074        | ~130                |
| MV (Exp) | 673       | 3782        | ~112                |

**Matérn por AIC:**

| Modelo | LogLik | AIC |
|--------|--------|-----|
| Matérn ν = 0.5 (Exp) | −459.02 | 926.03 |
| Matérn ν = 1.5 | −459.53 | 927.07 |

Preferencia débil por ν = 0.5 (ΔAIC ≈ 1.0, modelos esencialmente equivalentes).

**Anisotropía:**

| Modelo | LogLik | AIC |
|--------|--------|-----|
| Isotrópico | −459.02 | 926.03 |
| Anisotrópico (135°, ratio = 0.7) | −458.75 | 925.51 |

ΔAIC = 0.52. La dirección de mayor correlación (135°, NW-SE) es robusta y físicamente interpretable: perpendicular al gradiente de flujo SW → NE del acuífero.

---

## Librerías R

```r
library(geoR)      # likfit, as.geodata, geodata
library(gstat)     # variogram, fit.variogram, variogramLine, vgm
library(tidyverse) # dplyr, ggplot2, tibble
library(patchwork) # composición de gráficos
```

## Archivos

| Archivo | Descripción |
|---------|-------------|
| `Ay_03.pdf` | Enunciado de la ayudantía |
| `Ay_03_Espacial.R` | Script R comentado con toda la sesión |
