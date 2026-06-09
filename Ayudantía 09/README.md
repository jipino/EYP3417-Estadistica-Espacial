# Ayudantía 09: Procesos Puntuales — Intensidad, PPP, Funciones K/L y Modelos Log-lineales

**Curso:** EYP3417 Estadística Espacial  
**Profesor:** Alfredo Alegría  
**Ayudante:** Juan Pino  
**Período:** 2026  

---

## Descripción General

Esta ayudantía introduce el análisis de **patrones de puntos espaciales** mediante el Proceso Puntual de Poisson (PPP) — el modelo nulo de referencia en la disciplina — y las herramientas básicas de `spatstat` para:

1. **Caracterizar** un PPP y distinguir entre intensidad homogénea y no homogénea (estructura de **primer orden**)
2. **Estimar** la función de intensidad de forma no paramétrica mediante suavizamiento kernel
3. **Detectar** dependencia de **segundo orden** (atracción/repulsión) usando las funciones K y L de Ripley
4. **Ajustar** modelos paramétricos de intensidad log-lineal y **validar** su adecuación mediante bandas de confianza de Monte Carlo

A diferencia de las ayudantías de datos areales (06–08), aquí los datos son **realizaciones de procesos puntuales** (ubicaciones de eventos), y el foco pasa de matrices de vecindad a funciones de intensidad y estadísticos de segundo orden.

---

## Estructura del Script

El script `Ay_09_Espacial.R` sigue los tres problemas propuestos en el enunciado.

### Problema 1: Intensidad Espacial y el Proceso de Poisson (PPP)

**a) Definición y propiedades de un PPP**

Un PPP en una región $S$ queda caracterizado por dos propiedades:

- **Distribución de Poisson del conteo:** para toda subregión acotada $B \subseteq S$, $N(B) \sim \text{Poisson}(\mu(B))$ con $\mu(B) = \int_B \lambda(u)\, du$.
- **Independencia en regiones disjuntas:** los conteos $N(B_1), \dots, N(B_k)$ en subconjuntos disjuntos son independientes.

La diferencia entre un PPP **homogéneo** ($\lambda(u) = \lambda$ constante) y uno **no homogéneo** ($\lambda(u) = \lambda(x,y)$ variable) es puramente de **primer orden** (media): en ningún caso hay interacción entre puntos.

**b) Simulación con `rpoispp()`**

Se simulan, en la ventana unitaria $[0,1]^2$:
- Un PPP homogéneo con $\lambda = 100$.
- Un PPP no homogéneo con $\lambda(x,y) = 200 \cdot e^{-3x}$.

El PPP no homogéneo concentra sus puntos cerca de $x = 0$ (donde $\lambda$ es máxima, $\lambda(0,y) = 200$) y se diluye hacia $x = 1$ ($\lambda(1,y) = 200 e^{-3} \approx 9.96$), reflejando visualmente el decaimiento exponencial de la intensidad.

**c) Estimación kernel de la intensidad y el rol del bandwidth**

Se estima $\lambda(x,y)$ de forma no paramétrica con `density.ppp()` (kernel gaussiano), comparando un bandwidth seleccionado por la regla de Diggle (`bw.diggle()`) contra versiones $5$ veces más pequeña y más grande:

- **Bandwidth muy pequeño:** sobreajuste — el mapa retrata casi literalmente la nube de puntos individuales (alta varianza, bajo sesgo).
- **Bandwidth muy grande:** sobre-suavizado — el patrón de decaimiento exponencial se diluye hacia una superficie casi constante (alto sesgo, baja varianza).
- **Bandwidth óptimo:** balance sesgo-varianza que recupera visualmente la estructura $\lambda(x,y) = 200 e^{-3x}$.

---

### Problema 2: Análisis de Dependencia de Segundo Orden — Funciones K y L

**a) La función K de Ripley y su relación con L**

$K(h)$ mide el número esperado de puntos adicionales dentro de un radio $h$ de un punto típico, normalizado por $\lambda$:
$$K(h) = \frac{1}{\lambda}\, \mathbb{E}[\#\{\text{otros puntos a distancia} \leq h\}]$$

Bajo Aleatoriedad Espacial Completa (CSR), $K(h) = \pi h^2$, y la transformación estabilizadora $L(h) = \sqrt{K(h)/\pi}$ satisface $L(h) = h$. Por ello se grafica $L(h) - h$, cuyo valor esperado bajo CSR es $0$:

- $\widehat{L}(h) - h > 0$ → **agregación / atracción** (clustering) a esa escala.
- $\widehat{L}(h) - h < 0$ → **regularidad / repulsión** (inhibición, p. ej. hard-core) a esa escala.
- $\widehat{L}(h) - h \approx 0$ → compatible con CSR.

**b) Estimación empírica para el patrón `cells`**

Usando el conjunto `cells` de `spatstat` ($n = 42$ centros celulares), se estima $\widehat{L}(h) - h$ con `Lest()`. La curva resultante **no es monótona**: cae a un mínimo pronunciado ($\approx -0.09$ cerca de $h \approx 0.10$), cruza la referencia CSR cerca de $h \approx 0.17$–$0.18$, alcanza un leve máximo positivo ($\approx +0.005$ en $h \approx 0.20$) y vuelve a oscilar levemente bajo cero hacia $h = 0.25$.

**c) ¿Es razonable un modelo hard-core?**

Sí. El "valle profundo seguido de un pequeño rebote" es la firma típica de un proceso con **inhibición de corto alcance**: existe una distancia mínima aproximada (señalada por el mínimo de la curva, $\approx 0.10$) por debajo de la cual casi no se observan pares de células, y luego los centros se "acomodan" a una distancia algo mayor y más regular. Esto es coherente con la naturaleza física del fenómeno — las células ocupan espacio y sus núcleos no pueden solaparse —, lo que justifica un modelo de tipo **hard-core** (núcleo duro).

---

### Problema 3: Modelamiento Paramétrico y Validación mediante Monte Carlo

**a) Log-verosimilitud del modelo log-lineal**

Para $\log \lambda(u; \theta) = \theta_1 + \theta_2 x(u)$, la log-verosimilitud del PPP en una ventana $A$ es:
$$\ell_A(\theta) = \sum_{u_i \in A} \log \lambda(u_i; \theta) - \int_A \lambda(u; \theta)\, du = n\theta_1 + \theta_2 \sum_{i=1}^n x_i - \int_A e^{\theta_1 + \theta_2 x(u)}\, du$$

El primer término favorece intensidades altas donde se observaron puntos; el segundo penaliza modelos que sobreestiman la intensidad total. `ppm()` maximiza $\ell_A(\theta)$ numéricamente.

**b) Ajuste con `ppm()` e interpretación de $\hat\theta_2$**

Se ajusta `ppm(ppp_no_homogeneo ~ x)` al patrón simulado en 1(b). El estimador resulta $\hat\theta_2 < 0$ (su signo y orden de magnitud son coherentes con el valor poblacional $\theta_2 = -3$ de $\lambda(x,y) = 200 e^{-3x} = \exp(\log 200 - 3x)$, y su intervalo de confianza al 95% contiene a $-3$), confirmando que la intensidad decrece con $x$. El test de Wald asociado rechaza contundentemente $H_0: \theta_2 = 0$, evidenciando que la dependencia de la intensidad respecto a $x$ es estadísticamente muy significativa.

**c) Validación con bandas de Monte Carlo para $L(h) - h$**

Se construyen bandas de confianza al 95% (`envelope()`, $99$ simulaciones) bajo el modelo de Poisson log-lineal ajustado. La curva observada **sigue de cerca** la media simulada bajo $H_0$ y se mantiene **siempre dentro** de las bandas grises. Cabe destacar que tanto la curva observada como el envelope muestran una **tendencia creciente** con $h$ — esto **no** es evidencia de agregación, sino un artefacto esperado de la **heterogeneidad de primer orden** (que el envelope reproduce correctamente bajo $H_0$). La conclusión es que el modelo log-lineal ajustado captura adecuadamente toda la estructura del patrón, sin necesidad de invocar mecanismos de interacción de segundo orden — exactamente lo esperado, ya que el patrón fue generado como un PPP.

---

## Datos Utilizados

- **Simulados:** PPP homogéneo y no homogéneo en $[0,1]^2$ generados con `rpoispp()` (semilla `set.seed(2026)`).
- **`cells`:** conjunto de datos incluido en `spatstat` — centros de 42 células en un corte histológico microscópico.

---

## Librerías Requeridas

```r
library(spatstat)      # ppp, owin, rpoispp, density.ppp, Kest, Lest, ppm, envelope
library(ggplot2)       # Visualización
library(patchwork)     # Composición de gráficos
library(viridisLite)   # Paletas de color (plasma) para mapas de intensidad
```

---

## Cómo Ejecutar

```r
# En RStudio o terminal R:
source("Ay_09_Espacial.R")
```

El script ejecuta los tres problemas secuencialmente y genera:
- Gráficos comparativos de los PPP simulados (homogéneo vs. no homogéneo)
- Mapas de intensidad estimada por kernel para distintos bandwidths
- Gráfico de $L(h) - h$ para el patrón `cells`
- Tabla de coeficientes del modelo log-lineal ajustado y bandas de Monte Carlo para $L(h) - h$

---

## Lecciones Clave

1. **Primer orden vs. segundo orden:** la intensidad $\lambda(u)$ describe *dónde* es más probable observar puntos (primer orden); las funciones K/L describen si los puntos *interactúan* entre sí más allá de lo que dicta esa intensidad (segundo orden). Son preguntas distintas y requieren herramientas distintas.

2. **El bandwidth no es un detalle técnico:** la elección del parámetro de suavizado determina si se observa ruido muestral, una superficie sobre-simplificada, o la verdadera estructura de intensidad subyacente.

3. **$L(h) - h$ como diagnóstico visual:** transformar $K(h)$ a $L(h) - h$ permite comparar directamente contra una referencia horizontal en cero (CSR), facilitando detectar agregación (curva sobre cero) o repulsión (curva bajo cero) a distintas escalas $h$.

4. **No confundir heterogeneidad con interacción:** una tendencia no nula en $L(h) - h$ puede deberse enteramente a que la intensidad no es constante (efecto de primer orden), no a una verdadera atracción o repulsión entre puntos (efecto de segundo orden) — el envelope de Monte Carlo bajo el modelo ajustado permite distinguir ambos escenarios correctamente.

5. **Validación por simulación:** los envelopes de Monte Carlo son una herramienta general para contrastar si un modelo ajustado "explica" toda la estructura observada, comparando la curva empírica contra el rango de curvas que el propio modelo generaría.

---

**Contacto:** jipinov95@gmail.com  
**Última actualización:** Junio 2026  
**Tipo de sesión:** Introducción a procesos puntuales — intensidad, dependencia de segundo orden y modelamiento paramétrico
