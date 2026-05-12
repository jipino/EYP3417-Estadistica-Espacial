# Ayudantía 08: Evaluación Crítica de Dependencia Espacial en Delincuencia Comunal

**Curso:** EYP3417 Estadística Espacial  
**Profesor:** Alfredo Alegría  
**Ayudante:** Juan Pino  
**Período:** 2026  

---

## Descripción General

Esta ayudantía es un análisis **crítico y comparativo** de la dependencia espacial en tasas de delincuencia por comuna de la Región Metropolitana de Santiago (2018–2024). El objetivo es:

1. **Detectar** qué tipo de dependencia espacial existe (global y local)
2. **Validar** supuestos específicos de tres modelos competidores (OLS, CAR, SAR)
3. **Comparar** críticamente qué modelo es más defendible según diagnósticos
4. **Evaluar robustez** ante cambios de covariables
5. **Argumentar** una conclusión basada en evidencia, teoría y limitaciones

A diferencia de Ayudantía 06 (que implementó modelos), esta ayudantía enfatiza **evaluación de supuestos, diagnósticos espaciales y decisión de modelo** bajo incertidumbre.

---

## Estructura del Script

El script `Ay_08_Diagnostico_Espacial.R` sigue cinco problemas propuestos:

### Problema A: Diagnóstico Preliminar

**Objetivo:** Detectar si existe autocorrelación espacial significativa.

**Componentes:**

1. **Matriz de Vecindad W (Queen Adjacency):**
   - Se construye usando `poly2nb()` (convexidad de primer orden)
   - Se estandariza por fila con `nb2listw(style = "W")`
   - Resumen: Número de comunas, promedio de vecinos, rango

2. **Test de Moran I Global:**
   - Hipótesis nula: $I = 0$ (no hay autocorrelación espacial)
   - Estadístico: $I = \frac{n \sum_{i,j} w_{ij} z_i z_j}{\sum_i z_i^2}$
   - Interpretación:
     - $I > E[I]$ y $p < 0.05$ → **Clustering positivo significativo** (hot/cold spots)
     - $I < E[I]$ y $p < 0.05$ → **Dispersión significativa** (raro en delincuencia)
     - $p > 0.05$ → Compatible con aleatoriedad espacial

3. **Gráfico de Moran (Bivariado Estandarizado):**
   - Eje X: Variable estandarizada $Z_i = \frac{y_i - \bar{y}}{s}$
   - Eje Y: Lag espacial estandarizado $W·Z$
   - Cuadrantes:
     - **HH:** Hot spots (alta tasa, vecinos altos) → acción prioritaria
     - **LL:** Cold spots (baja tasa, vecinos bajos) → replicar prácticas
     - **HL:** Outliers negativos (alta tasa, vecinos bajos) → anomalía local
     - **LH:** Outliers positivos (baja tasa, vecinos altos) → isla de paz

**Interpretación Conjunta:**
- Si Moran I es significativo y hay concentración HH+LL → Clustering fuerte
- Identifica comunas específicas para intervención diferenciada

---

### Problema B: Test de Especificación

**Objetivo:** Evaluar si OLS es suficiente o si hay estructura espacial omitida.

**Componentes:**

1. **Regresión OLS Base:**
   - Modelo: $Y_i = \beta_0 + \beta_1 X_i + \varepsilon_i$ (sin estructura espacial)
   - Variables: `tasa_delitos ~ ingreso_pc`
   - Salida: Coeficientes, $R^2$, p-values
   - Interpretación: ¿Ingreso es predictor significativo? ¿Qué magnitud?

2. **Test de Moran I en Residuos:**
   - **Hipótesis crítica:** Después de controlar ingreso, ¿quedan residuos con autocorrelación?
   - Si $I_{residuos} > 0$ y $p < 0.05$: **OLS OMITE estructura espacial importante**
   - Si $I_{residuos} \approx 0$ y $p > 0.05$: La dependencia espacial se explica por covariables
   - Interpretación: El test de Moran en residuos es el **puente entre diagnóstico y modelado**

3. **Conclusión sobre Necesidad de Modelo Espacial:**
   - Si Moran I residuos significativo: Modelo espacial (CAR/SAR) es necesario
   - Si no significativo: OLS es aceptable, pero CAR/SAR pueden mejorar por eficiencia

---

### Problema C: Validación de Supuestos Específicos

**Objetivo:** Comparar supuestos y diagnósticos de CAR vs SAR.

#### C1: Modelo CAR (Conditional Autoregressive)

**Especificación:**
$$E[Y_i | Y_{-i}, X] = \beta_0 + \beta_1 X_i + \lambda \sum_{j \sim i} w_{ij} Y_j$$
$$\text{Var}(Y_i | Y_{-i}, X) = \sigma^2$$

**Supuestos:**
- Distribuciones condicionales gaussianas locales
- Markov local: $Y_i | Y_{-i}$ depende solo de vecinos
- Simetría de la matriz de precisión

**Parámetro $\lambda$:**
- Interpretación: Dependencia CONDICIONAL espacial
- Si $\lambda > 0$ y significativo: Clustering; si $\lambda < 0$: Dispersión
- Rango típico: $[-1, 1]$ (para matriz de pesos estandarizada)

**Diagnósticos:**
- Normalidad de residuos (Shapiro-Wilk test)
- Q-Q plot de residuos
- Gráfico de valores ajustados vs residuos (heterocedasticidad)

**Ventajas:**
- Más flexible: no asume causalidad simultánea
- Interpretación local: cada región afecta su vecindario

#### C2: Modelo SAR (Spatial Autoregressive / Lag Espacial)

**Especificación:**
$$Y = \rho W Y + X \beta + \varepsilon$$
$$\text{Donde } \varepsilon \sim N(0, \sigma^2 I)$$

**Supuestos:**
- Lag espacial es predictor válido (no endógeno)
- **Simultaneidad:** $Y_i$ afecta $Y_j$ y $Y_j$ afecta $Y_i$ a la vez (two-way feedback)
- No hay dependencia en residuos (errores independientes)

**Parámetro $\rho$:**
- Interpretación: Parámetro de lag espacial (autoregresión espacial)
- Si $\rho > 0$ y significativo: Spillover positivo (interdependencia)
- Rango: $(\lambda_{min}^{-1}, \lambda_{max}^{-1})$ (valores propios de W)

**Multiplicadores Espaciales:**
- **Efecto directo:** Cambio en $X_i$ sobre $Y_i$ (parcial)
- **Efecto indirecto:** Propagación a través de vecinos
- **Efecto total:** Suma de ambos (puede ser > efecto directo)
- Interpretación: Si $\rho$ es grande, hay "contagio" importante entre comunas

**Diagnósticos:**
- Normalidad de residuos (igual que CAR)
- Two-way feedback es realista? ¿Merece modelo SAR?

#### C3: Comparación de Modelos

| Aspecto | OLS | CAR | SAR |
|---------|-----|-----|-----|
| **Supuesto independencia** | SÍ (estricto) | NO | NO |
| **Markov local** | — | SÍ | — |
| **Simultaneidad bidireccional** | NO | NO | SÍ |
| **Mecanismo causal** | Directo solo | Contagio condicional | Retroalimentación mutua |
| **Interpretabilidad** | Fácil | Media | Compleja (multiplicadores) |
| **Requisitos datos** | Mínimos | Moderados | Altos (exogeneidad de W·Y) |

---

### Problema D: Robustez y Sensibilidad

**Objetivo:** Evaluar si la dependencia espacial es robusta a cambios de covariables.

**Estrategia:**
1. Reajustar OLS, CAR, SAR usando **densidad poblacional** en lugar de ingreso
2. Comparar cambios en parámetros espaciales ($\lambda$, $\rho$)

**Indicadores de Robustez:**
- $|\Delta \lambda / \lambda| < 20\%$: Cambio pequeño → robustez
- $|\Delta \rho / \rho| < 20\%$: Estable
- Persistencia de significancia: ¿Se mantiene el p-value < 0.05?

**Interpretación:**
- **Cambios sustanciales** → Dependencia espacial es **confusión de covariable** (no genuina)
- **Estabilidad** → Dependencia espacial es **fenómeno genuino**, no artefacto

---

### Problema E: Conclusión Integral

**Estructura de Argumentación:**

1. **Supuestos Teóricos Cumplidos:**
   - CAR: ¿Normalidad de residuos? ¿Es realista Markov local?
   - SAR: ¿Normalidad? ¿Es realista two-way feedback?
   - OLS: Si hay autocorrelación residual, supuesto de independencia es violado

2. **Interpretación Causal:**
   - ¿Qué significa $\lambda$ o $\rho$ en contexto de delincuencia?
   - ¿Son mecanismos plausibles? (contagio, desplazamiento, operaciones conjuntas)

3. **Diagnósticos:**
   - Moran I global: ¿Hay clustering inicial?
   - Moran I residuos OLS: ¿Hay estructura omitida?
   - Normalidad: ¿Qué modelo tiene mejores residuos?
   - Robustez: ¿Cambian los parámetros con covariables?

4. **Limitaciones de Datos/Contexto:**
   - Tamaño muestral: $n = 51$ comunas (pequeño)
   - Matriz de pesos: Queen es discrecional, ¿qué ocurre con Rook o distancia?
   - Omisión de variables: Educación, policía, capital social omitidas
   - Heterogeneidad espacial: ¿Central vs periférico tienen estructuras diferentes?
   - Endogeneidad: En SAR, ¿es W·Y realmente exógeno?

**Conclusión Defendible:**
Después de toda la evidencia, argumentar **cuál modelo (OLS, CAR, o SAR) es más defendible**, reconociendo que ninguno es perfecto y cada uno tiene limitaciones.

---

## Datos Utilizados

- **Delincuencia:** `cead_delincuencia_chile.parquet` (Centro de Estadísticas y Análisis de Delito, 2018–2024)
- **Comunas:** Shapefile de Chile, 52 comunas de la RM (51 tras filtrado de NAs)
- **Covariables:**
  - Ingreso per cápita (aproximación CASEN 2022, miles de pesos mensuales)
  - Densidad poblacional (hab/km², aproximada)

---

## Librerías Requeridas

```r
library(tidyverse)      # dplyr, ggplot2, tidyr
library(lubridate)      # Manejo de fechas
library(sf)             # Objetos espaciales (shapefiles)
library(spdep)          # Análisis de dependencia espacial (moran.test, poly2nb)
library(spatialreg)     # Regresión espacial (spautolm para CAR, lagsarlm para SAR)
library(arrow)          # Lectura de parquet
library(ggplot2)        # Visualización
library(patchwork)      # Composición de gráficos
```

---

## Cómo Ejecutar

```bash
# En RStudio o terminal R:
source("Ay_08_Diagnostico_Espacial.R")
```

El script ejecuta todos los problemas secuencialmente y genera:
- Tablas de resultados en consola
- Gráficos de Moran
- Diagnósticos de normalidad (histogramas, Q-Q plots)
- Síntesis comparativa y recomendación final

---

## Interpretaciones Principales

### ¿Hay Dependencia Espacial?

**Pregunta:** ¿Las tasas de delincuencia de comunas vecinas están relacionadas?

- **Sí, significativa (Moran I p < 0.05):** Hay clustering. Necesario modelo espacial.
- **No (Moran I p > 0.05):** Compatible con aleatoriedad. OLS podría ser suficiente.

### ¿OLS es Suficiente?

**Pregunta:** Después de controlar ingreso, ¿desaparece la autocorrelación espacial?

- **Sí (Moran I residuos p > 0.05):** OLS captura esencialmente la estructura.
- **No (Moran I residuos p < 0.05):** Hay estructura omitida. CAR/SAR son necesarios.

### ¿CAR o SAR?

**Pregunta:** ¿Cuál mecanismo es más realista?

- **CAR:** Dependencia condicional (tasas de una comuna condicionadas en vecinas)
- **SAR:** Retroalimentación simultánea (efectos bidireccionales)

**Decisión:** Basado en normalidad, significancia de parámetros y plausibilidad causal.

### ¿Es Robusto?

**Pregunta:** ¿Los parámetros cambian mucho con densidad en lugar de ingreso?

- **Cambios < 20%:** Robusto. El fenómeno espacial es genuino.
- **Cambios > 20% o p-value pierde significancia:** Sensible a covariables. Posible confusión.

---

## Lecciones Clave

1. **Autocorrelación Global ≠ Dependencia Causal:** Moran I significativo no prueba que SAR/CAR sean necesarios; podría ser confusión de covariables.

2. **Test de Moran en Residuos es Crítico:** Es el puente entre diagnóstico exploratorio (Prob. A) y decisión de modelo (Prob. C).

3. **CAR vs SAR:** CAR es más flexible (no asume simultaneidad); SAR es más exigente pero puede capturar mecanismos de retroalimentación mutua.

4. **Supuestos Importan:** Normalidad de residuos, exogeneidad de covariables, especificación correcta de W.

5. **Contexto Domina:** Sin teoría sobre *por qué* las comunas están relacionadas, cualquier modelo es especulativo.

6. **Humildad Estadística:** $n = 51$ es pequeño; inferencia es incierta; múltiples modelos compiten.

---

## Próximas Investigaciones

- Incluir educación, policía, desigualdad como covariables
- Explorar estructuras alternativas de W (Rook, distancia, k-vecinos)
- Modelos dinámicos (cambios 2018–2024)
- Heterogeneidad espacial (regímenes locales)
- Análisis de LISA (Local Indicators of Spatial Association)
- Validación cruzada espacial (train/test sets)
- Distribuciones alternativas (Poisson, Binomial Negativa para conteos)

---

**Contacto:** jipinov95@gmail.com  
**Última actualización:** Mayo 2026  
**Tipo de sesión:** Práctica con énfasis en evaluación crítica de supuestos y decisión de modelo
