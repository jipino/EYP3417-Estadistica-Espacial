# Ayudantía 06: Datos Areales — Modelos Ising y CAR

**Curso:** EYP3417 Estadística Espacial  
**Profesor:** Alfredo Alegría  
**Ayudante:** Juan Pino  
**Período:** 2026  

## Descripción

Esta ayudantía explora el modelamiento de datos areales (comunas de la Región Metropolitana de Santiago) usando dos enfoques clave en estadística espacial:

- **Modelo de Ising**: Captura la interacción espacial binaria (presencia/ausencia de delincuencia alta)
- **Modelo CAR (Conditional Autoregressive)**: Modela tasas continuas con dependencia espacial condicional

El análisis integra datos de delincuencia (CEAD 2018–2024) con información socioeconómica (ingreso per cápita), revelando la estructura espacial y los determinantes de la criminalidad en Santiago.

## Contenido

### Problema 1: Visualización de Datos Areales y Estructura de Red

**Objetivos:**
- Cargar y explorar shapefiles (comunas RM)
- Construir matriz de pesos espaciales (Queen contiguity)
- Visualizar tasas de delincuencia cruda y transformada
- Aplicar transformación Freeman-Tukey para estabilizar varianza
- Graficar la estructura de vecindad (grafo de dependencia)

**Secciones:**
- **1a)** Shapefile y matriz de pesos espaciales W (poly2nb, nb2listw)
- **1b)** Datos de delincuencia (CEAD) y mapa coroplético de tasa cruda
- **1c)** Grafo de vecindades superpuesto sobre el mapa
- **Freeman-Tukey:** Transformación √(x/n) + √((x+1)/n) para varianza estable
- **Evolución anual:** Facet-wrap con 8 años (2018–2025, parcial)
- **Comparativa:** Mapa cruda vs Freeman-Tukey lado a lado

**Resultados clave:**
- Gradiente oriente-occidente claro: oriente (Las Condes, Providencia) con tasas bajas; norte/sur-poniente (Colina, Melipilla) con tasas altas
- Clustering espacial evidente: comunas vecinas tienden a compartir niveles similares
- Impacto COVID-19 visible en 2020: transición de tasas altas a tasas más bajas
- Transformación Freeman-Tukey modera outliers de comunas pequeñas (Tiltil, Alhué) sin distorsionar patrón urbano

### Problema 2: Modelo de Ising

**Objetivos:**
- Binarizar tasas en la mediana regional
- Estimar parámetros de Ising usando pseudoverosimilitud (glm)
- Interpretar atracción/repulsión espacial mediante β
- Calcular probabilidades condicionales

**Especificación matemática:**
```
π(x) ∝ exp(α · Σ_i x_i + β · Σ_{i~j} x_i·x_j)

Condicional logística:
π_i(1 | x_{-i}) = exp(α + β·Σ_{j~i} x_j) / (1 + exp(α + β·Σ_{j~i} x_j))
```

**Secciones:**
- **2a)** Binarización: X_i = 1 si tasa_ft > mediana regional
- **2b)** Gráfico de sensibilidad: variación de β sobre dos escenarios
- **2c)** Estimación por pseudoverosimilitud (glm, family=binomial)
- **2d)** Probabilidades condicionales estimadas

**Resultados clave:**
- **α̂ = -0.5992** (p=0.367, NO SIG.): prevalencia global baja (~35%)
- **β̂ = 0.2085** (p=0.324, NO SIG.): atracción espacial débil, OR=1.23
- **Clustering visual vs no-significancia:** Mapa binario muestra clustering claro (rojo-rojo, azul-azul), pero estimación por pseudoverosimilitud es débil
- **Probabilidades condicionales:** P(X_i=1|todos vecinos=1)=0.7444 vs P(X_i=1|todos vecinos=0)=0.3545 → efecto neto de 39 pp sustancial pero β no significativo

### Problema 3: Modelo CAR (Conditional Autoregressive)

**Objetivos:**
- Ajustar modelo CAR con covariable (ingreso per cápita)
- Estimar parámetro de dependencia espacial λ (ρ)
- Comparar OLS (sin estructura espacial) vs CAR (con dependencia)
- Diagnosticar residuos y detectar autocorrelación

**Especificación matemática:**
```
X ~ N(0, (I - B)^{-1} K)
donde B = λ·N (N = matriz binaria de vecinos)
      K = σ² · I

E(X_i | X_{-i}) = λ · Σ_{j~i} x_j / n_i
Var(X_i | X_{-i}) = σ²
```

**Secciones:**
- **3a)** Especificación del modelo CAR con covariables
- **3b)** Ajuste mediante spautolm (máxima verosimilitud)
- **3c)** Medias condicionales del CAR (lag espacial)
- **3d)** Comparación de residuos OLS vs CAR
- **3e)** Test de Moran's I sobre residuos (futuro)

**Resultados clave:**
- **λ̂ = -0.0310** (p=0.9343, NO SIG.): dependencia espacial prácticamente nula
- **β̂₀ = 0.4795** (p<2e-16, SIG.): intercepto significativo
- **β̂_ingreso = -0.000035** (p=0.4552, NO SIG.): efecto débil e insignificante
- **R² = 0.00918** (OLS): ingreso explica menos del 1% de varianza
- **Residuos OLS vs CAR:** Mapas prácticamente idénticos, confirma que λ̂≈0 no absorbió estructura espacial
- **Conclusión:** Ingreso es mal predictor de delincuencia; factores no observados (educación, densidad, desigualdad, policía) son probablemente determinantes

## Datos

- **comunas.shp**: Shapefile de comunas de Chile (RM filtrada: 52 comunas)
- **cead_delincuencia_chile.parquet**: Datos del Centro de Estadísticas y Análisis de Delito (2018–2024)
  - Fuente: https://github.com/bastianolea/delincuencia_chile
- **Ingreso per cápita**: Elaboración aproximada basada en CASEN 2022 (miles de pesos mensuales)

## Librerías requeridas

```r
library(tidyverse)      # dplyr, ggplot2, tidyr
library(lubridate)      # year(), date handling
library(sf)             # Spatial objects, st_read(), st_geometry()
library(spdep)          # poly2nb(), nb2listw(), moran.test()
library(spatialreg)     # spautolm() para CAR
library(arrow)          # read_parquet()
library(ggplot2)        # geom_sf(), facet_wrap(), scales
library(mapview)        # Interactive maps
library(patchwork)      # Plot composition (+)
```

## Cómo ejecutar

1. **Descargar archivos:** Asegúrate que `comunas.shp` y `cead_delincuencia_chile.parquet` estén en el mismo directorio que el script
2. **Ejecutar script:** `source("Ay_06_Espacial.R")` en R
3. **Resultados:** Se generarán mapas (PNG/PDF automáticamente) y salida en consola con estadísticas

## Interpretaciones principales

### Pregunta central: ¿La delincuencia en la RM depende del ingreso y la ubicación?

**Respuesta:**
- **Visualmente:** Sí (Problema 1): patrón claro oriente-occidente, clustering evidente
- **Estadísticamente:** Débilmente o no (Problemas 2–3):
  - Ising: β̂ positivo pero no significativo (p=0.324)
  - CAR: λ̂ ≈ 0 (p=0.934), ingreso no significativo (p=0.455)
  - R² muy bajo (0.009): ingreso explica < 1%

**Conclusión:**
La delincuencia en Santiago está **determinada por factores más complejos** que solo ingreso e interdependencia espacial simple. Candidatos:
- Educación (nivel sociocultural)
- Densidad poblacional y urbanización
- Desigualdad interna dentro de comunas
- Capacidad y presencia policial
- Infraestructura de vigilancia (cámaras, iluminación)
- Capital social y cohesión comunitaria
- Dinámicas criminales históricas

### Lecciones sobre estadística espacial

1. **Clustering visual ≠ Dependencia espacial significativa:** Problema 1 muestra clustering claro, pero Ising/CAR no lo capturan como significativo. Puede deberse a covariables que generan similitud espacial sin dependencia directa.

2. **Freeman-Tukey es crucial para datos areales:** Estabiliza varianza en comunas pequeñas, mejora legibilidad de mapas sin perder información.

3. **OLS vs CAR:** Cuando λ̂ ≈ 0, no hay mejora CAR. Esto ocurre cuando las covariables explican la estructura espacial aparente.

4. **R² bajo es una advertencia:** Indica que el modelo omite determinantes importantes. Expandir covariables, considerar no-linealidades, o aceptar heterogeneidad irreducible.

## Próximas sesiones

- Incluir covariables adicionales (educación, densidad poblacional)
- Explorar relaciones no-lineales (spline espacial)
- Modelo SAR (Spatial Autoregressive) como alternativa a CAR
- Análisis de cambio temporal con modelos dinámicos

## Referencias

- Alegría, A. (2025). *Estadística Espacial*. Apuntes de cátedra.
- Cressie, N. (1993). *Statistics for Spatial Data* (Revised Edition). Wiley.
- Waller, L. A., & Gotway, C. A. (2004). *Applied Spatial Statistics for Public Health Data*. Wiley.
- Bastiaan Olea. *Delincuencia en Chile* [Repositorio GitHub]. https://github.com/bastianolea/delincuencia_chile

---

**Contacto:** jipinov95@gmail.com  
**Última actualización:** Abril 2026
