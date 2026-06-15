# Ayudantía 10: Tests de CSR — Funciones G, F, K y L con Datos Reales

**Curso:** EYP3417 Estadística Espacial  
**Profesor:** Alfredo Alegría  
**Ayudante:** Juan Pino  
**Período:** 2026  

---

## Descripción General

Esta ayudantía aplica los tests estadísticos de **Aleatoriedad Espacial Completa (CSR)** sobre patrones de puntos reales, utilizando cuatro funciones clásicas de `spatstat`: **G**, **F**, **K** y **L**. El enfoque es práctico: se comparan patrones con estructuras espaciales conocidas (aleatorio, inhibido, agrupado) y se analiza un caso real de localización de tiendas comerciales.

Los dos ejercicios permiten:
1. **Calibrar la lectura** de los tests G y F sobre patrones de referencia con comportamiento espacial conocido (CSR, repulsión, clustering).
2. **Aplicar** los tests K y L sobre un patrón real (Starbucks en Massachusetts), razonando sobre causas geográficas y socioeconómicas del agrupamiento.

---

## Estructura del Script

### Ejercicio 1: Tests G y F sobre patrones de referencia

Se cargan tres conjuntos de datos clásicos de `spatstat`:

| Dataset | n | Patrón esperado |
|---|---|---|
| `japanesepines` | 65 | Aleatorio (compatible con CSR) |
| `cells` | 42 | Inhibido / repulsivo (hard-core) |
| `redwood` | 62 | Agrupado (clustering espacial) |

Todos están reescalados a la ventana unitaria $[0,1]^2$.

#### Función G — distancia al vecino más cercano

$G(r)$ mide la proporción de puntos cuyo vecino más cercano está a distancia $\leq r$:

$$G(r) = P(d_{\min} \leq r)$$

Bajo CSR, $G(r) = 1 - e^{-\lambda \pi r^2}$. La comparación con las bandas de confianza de Monte Carlo (100 simulaciones, $\alpha = 0.05$, nivel efectivo $\approx 0.0198$) permite detectar:

- Curva observada **sobre** la banda → $G$ crece más rápido de lo esperado → vecinos más **cercanos** de lo aleatorio → **clustering**.
- Curva observada **bajo** la banda → $G$ crece más lento → vecinos más **lejanos** de lo aleatorio → **repulsión / inhibición**.

**Resultados:**
- `japanesepines`: curva dentro de la banda → **no se rechaza CSR**.
- `cells`: curva bajo la banda en la mayor parte del rango → se rechaza CSR, **repulsión** confirmada.
- `redwood`: curva sobre la banda → se rechaza CSR, **clustering** confirmado.

#### Función F — función de espacio vacío

$F(r)$ mide la proporción de puntos de una grilla aleatoria (no del patrón) cuyo punto más cercano del patrón está a distancia $\leq r$:

$$F(r) = P(d(u, \mathbf{X}) \leq r), \quad u \text{ uniformemente distribuido en } S$$

Bajo CSR, $F(r) = G(r)$. La función F es sensible a la *cantidad de espacio vacío* y complementa a G:

- Curva observada **sobre** la banda → hay menos espacio vacío de lo esperado → puntos distribuidos más uniformemente → **inhibición**.
- Curva observada **bajo** la banda → hay más espacio vacío de lo esperado → puntos concentrados → **clustering**.

**Resultados:**
- `japanesepines`: curva dentro de la banda → **no se rechaza CSR**.
- `cells`: curva sobre la banda → se rechaza CSR, **inhibición** confirmada.
- `redwood`: curva bajo la banda → se rechaza CSR, **clustering** confirmado.

---

### Ejercicio 2: Tests K y L sobre distribución de Starbucks en Massachusetts

Se carga el dataset `starbucks` ($n = 171$ tiendas) junto a la ventana irregular `ma` (estado de Massachusetts, coordenadas UTM). La visualización inicial revela una fuerte concentración de tiendas en la región noreste (Greater Boston), con zonas vacías en el interior del estado — consistente con un patrón agrupado.

#### Función K — conteo acumulado de vecinos

$K(h)$ cuenta el número esperado de puntos adicionales dentro de un radio $h$ de un punto típico, normalizado por $\lambda$. Bajo CSR, $K(h) = \pi h^2$. La curva observada **por encima** de la banda indica más vecinos que bajo CSR → **clustering**.

**Resultado:** la curva observada supera las bandas de confianza a lo largo de todo el rango $[0, 38488]$ m → se rechaza CSR con fuerte evidencia de agrupamiento.

#### Función L — transformación estabilizadora de K

$L(h) = \sqrt{K(h)/\pi}$, que satisface $L(h) = h$ bajo CSR, facilitando la comparación visual. La curva observada **consistentemente por encima** de la banda confirma que la distancia acumulada entre tiendas es significativamente menor a la esperada bajo aleatoriedad.

**Resultado:** se rechaza CSR. Existe evidencia muy fuerte de **agrupamiento espacial** en la distribución de tiendas Starbucks de Massachusetts, consistente con la localización preferencial en zonas urbanas densamente pobladas.

---

## Datos Utilizados

- **`japanesepines`, `cells`, `redwood`:** incluidos en `spatstat` (Ripley, Diggle).
- **`starbucks`, `ma`:** cargados desde el repositorio de Gimond (`mgimond/Spatial` en GitHub) — tiendas Starbucks y polígono del estado de Massachusetts en coordenadas UTM.

---

## Librerías Requeridas

```r
library(spatstat)
library(spatstat.geom)
library(spatstat.explore)
```

---

## Cómo Ejecutar

```r
source("Ay_10_Espacial.R")
```

El script ejecuta los dos ejercicios secuencialmente y genera:
- Panel comparativo de los tres patrones de referencia
- Bandas de confianza de la función G para los tres patrones
- Bandas de confianza de la función F para los tres patrones
- Gráfico de distribución de Starbucks en Massachusetts
- Bandas de confianza de la función K para Starbucks
- Bandas de confianza de la función L para Starbucks

---

## Lecciones Clave

1. **G y F son complementarias:** G mide distancias *entre puntos del patrón*; F mide distancias *desde el espacio vacío al patrón más cercano*. Clustering y repulsión producen desviaciones en direcciones opuestas en ambas funciones, lo que permite un diagnóstico robusto al usarlas en conjunto.

2. **K y L amplían el análisis a todas las escalas:** mientras G y F resumen en una sola estadística (el vecino más cercano), K y L acumulan información a múltiples radios $r$, detectando estructuras de clustering o repulsión que operan a distintas escalas.

3. **Significancia Monte Carlo:** con 100 simulaciones, el nivel efectivo del test es $2/101 \approx 0.0198$ — más conservador que el $\alpha = 0.05$ nominal, reduciendo falsos positivos.

4. **La geometría de la ventana importa:** en el ejercicio de Starbucks, la ventana irregular de Massachusetts requiere correcciones de borde (isotrópica en K y L) para no confundir efectos de frontera con clustering real.

5. **El patrón refleja el fenómeno subyacente:** el agrupamiento de Starbucks no es aleatorio — responde a densidad poblacional, infraestructura vial y poder adquisitivo, variables que el análisis espacial detecta pero no explica por sí solo.

---

**Contacto:** jipinov95@gmail.com  
**Última actualización:** Junio 2026  
**Tipo de sesión:** Práctica — tests de CSR con funciones G, F, K y L sobre patrones reales
