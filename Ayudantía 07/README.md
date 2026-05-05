# Ayudantía 07: Campos Aleatorios de Markov, Modelos CAR y Procesos de Poisson Espaciales

**Curso:** EYP3417 Estadística Espacial  
**Profesor:** Alfredo Alegría  
**Ayudante:** Juan Pino  
**Período:** 2026  

---

## Descripción

Esta ayudantía teórica profundiza en los fundamentos matemáticos de tres temas fundamentales en estadística espacial:

1. **Campos Aleatorios de Markov (MRF):** Propiedades de independencia condicional en estructuras bipartitas
2. **Modelos CAR (Conditional Autoregressive):** Especificación de distribuciones condicionales gaussianas con dependencia espacial
3. **Procesos de Poisson Espaciales:** Condiciones de existencia y estabilidad en autorregresiones de conteos

El foco está en **demostraciones rigurosas**, **teoremas fundamentales** (Hammersley-Clifford, Brook/Besag) y el **análisis de consistencia** de modelos espaciales.

---

## Contenido por Problema

### Problema 1: Campos Aleatorios de Markov en Grillas Bipartitas

**Contexto:** Sea $T \subset \mathbb{Z}^2$ una grilla finita particionada en dos conjuntos disjuntos $S_1$ y $S_2$ en patrón de tablero de ajedrez (ningún par dentro de $S_1$ o $S_2$ es vecino).

#### Parte a) Independencia Condicional en MRF

**Objetivo:** Demostrar que $\{X_i : i \in S_1\}$ son condicionalmente independientes dado $\{X_j : j \in S_2\}$.

**Estrategia:**
- Aplicar la propiedad de Markov local: $\pi(x_i | x_{T\setminus\{i\}}) = \pi(x_i | x_{\partial i})$
- Observar que todo vecino de $i \in S_1$ pertenece a $S_2$: $\partial i \subseteq S_2$
- Demostrar mediante inducción finita que la distribución conjunta factoriza:
$$\pi(x_{S_1} | x_{S_2}) = \prod_{i \in S_1} \pi(x_i | x_{S_2})$$

**Conceptos clave:**
- Estructura bipartita del grafo de vecindad
- Ausencia de "path" de comunicación directa entre sitios en $S_1$
- Cálculo mediante regla de la cadena e independencia condicional

#### Parte b) Teorema de Hammersley-Clifford y Cliques

**Objetivo:** Justificar la simplificación de la densidad conjunta identificando los cliques.

**Teorema (Hammersley-Clifford):** Para un MRF positivo,
$$\pi(x) = Z^{-1} \exp\left(\sum_{C \in \mathcal{C}} \varphi_C(x_C)\right)$$

donde $\mathcal{C}$ es la colección de cliques maximos del grafo de vecindad.

**Identificación de cliques en estructura de ajedrez:**

La bipartición implica que no existen cliques de tamaño $\geq 3$ (no hay triángulos). Los únicos cliques posibles son:

| Tipo | Forma | Descripción |
|------|-------|-------------|
| **Singletons** | $\{i\}$, $i \in T$ | Cada sitio por sí solo |
| **Pares** | $\{i,j\}$, $i \sim j$ | Aristas del grafo (siempre: $i \in S_1, j \in S_2$) |

**Simplificación resultante:**
$$\pi(x) = Z^{-1} \exp\left(\sum_{i \in T} \varphi_{\{i\}}(x_i) + \sum_{i \sim j} \varphi_{\{i,j\}}(x_i, x_j)\right)$$

**Interpretación:** La estructura bipartita **elimina términos de interacción de orden 3 o superior**, reduciendo la complejidad computacional y reflejando que la única comunicación entre $S_1$ y $S_2$ pasa a través de aristas directas.

---

### Problema 2: Modelos de Autorregresi´on Condicional (CAR)

**Contexto:** Especificación condicional gaussiana:
$$X_i | x_{T\setminus\{i\}} \sim N\left(\sum_{j \neq i} b_{ij} x_j, \kappa_i\right)$$

donde $B = (b_{ij})$ es una matriz de pesos (diagonal 0) y $K = \text{diag}(\kappa_i)$ una matriz diagonal de varianzas.

#### Parte a) Propiedad de Markov Espacial

**Objetivo:** Demostrar que si $b_{ij} = 0$ para $i \not\sim j$, entonces la distribución gaussiana resultante es un MRF.

**Demostración:**
- La media condicional $\mu_i = \sum_{j \in \partial i} b_{ij} x_j$ (solo depende de vecinos)
- La varianza $\kappa_i$ **no depende** de $x_{T\setminus\{i\}}$
- Por tanto: $\pi(x_i | x_{T\setminus\{i\}}) = \pi(x_i | x_{\partial i})$ ✓ (Markov local)

**Consecuencia:** El campo aleatorio gaussiano es un MRF con respecto a la vecindad $\sim$.

#### Parte b) Funciones de Interacción de Clique

La densidad conjunta gaussiana es:
$$\pi(x) \propto \exp\left(-\frac{1}{2}x^\top Q x\right), \quad Q = K^{-1}(I - B)$$

con entradas $Q_{ii} = 1/\kappa_i$ y $Q_{ij} = -b_{ij}/\kappa_i$ (para $i \neq j$).

Expandiendo la forma cuadrática:
$$-\frac{1}{2}x^\top Q x = -\frac{1}{2}\sum_{i \in T} \frac{x_i^2}{\kappa_i} + \sum_{\substack{i,j \in T \\ i \sim j, i < j}} \frac{b_{ij}}{\kappa_i} x_i x_j$$

**Cliques y potenciales:**

| Clique | Función | Interpretación |
|--------|---------|-----------------|
| $\{i\}$ | $\varphi_{\{i\}}(x_i) = -\frac{x_i^2}{2\kappa_i}$ | Penalización cuadrática; escala proporcional a $\kappa_i$ |
| $\{i,j\}, i \sim j$ | $\varphi_{\{i,j\}}(x_i,x_j) = \frac{b_{ij}}{\kappa_i} x_i x_j$ | Dependencia pairwise; $b_{ij} > 0 \Rightarrow$ suavidad; $b_{ij} < 0 \Rightarrow$ repulsión |

#### Parte c) Importancia de la Simetría de Q

**Pregunta central:** ¿Por qué es crítico que $Q = K^{-1}(I - B)$ sea simétrica y definida positiva?

**Respuesta:**

1. **Simetría** $\Rightarrow$ compatibilidad de condicionales:
   - La simetría $Q_{ij} = Q_{ji}$ impone la **condición de reciprocidad**:
   $$\frac{b_{ij}}{\kappa_i} = \frac{b_{ji}}{\kappa_j}$$
   - Esto asegura que existe una densidad conjunta $\pi(x)$ que genera simultáneamente **todas** las condicionales especificadas (coherencia del modelo).
   - Sin simetría, los pesos $b_{ij}$ y $b_{ji}$ son incompatibles: no existe $\pi$ que los reproduzca.

2. **Definida positiva** $\Rightarrow$ normalizabilidad:
   - Garantiza que $Z = \int \exp(-\frac{1}{2}x^\top Q x) \, dx < \infty$ (la integral converge).
   - Sin definida positiva, $Z = +\infty$ y $\pi$ no es una distribución válida.

**Conclusión:** Simetría + definida positiva son condiciones **necesarias y suficientes** para que el CAR sea un MRF gaussiano bien definido. El enunciado asume que $(I-B)^{-1}K$ posee ambas propiedades, lo que garantiza la consistencia.

---

### Problema 3: Autorregresión de Poisson Espacial

**Contexto:** Intento de definir un modelo con condicionales Poisson:
$$\pi_i(y_i | y_{T\setminus\{i\}}) = \frac{e^{-\mu_i} \mu_i^{y_i}}{y_i!}, \quad \log \mu_i = \theta \sum_{j \sim i} y_j$$

donde $\theta$ es un parámetro de interacción espacial.

#### Parte a) Lema de Brook (Teorema de Factorización de Besag)

**Objetivo:** Derivar la forma funcional de la densidad conjunta.

**Lema (Brook):** Para cualquier distribución conjunta positiva y ordenamiento de sitios:
$$\frac{\pi(y)}{\pi(0)} = \prod_{i=1}^n \frac{\pi_i(y_i | y_1, \ldots, y_{i-1}, 0, \ldots, 0)}{\pi_i(0 | y_1, \ldots, y_{i-1}, 0, \ldots, 0)}$$

**Cálculo paso a paso:**

En el paso $i$, los sitios $j < i$ tienen valor $y_j$ y los sitios $j > i$ tienen valor 0. La media condicional es:
$$\mu_i^{(i)} = \exp\left(\theta \sum_{j \sim i, j < i} y_j\right)$$

El cociente del factor $i$ es:
$$\frac{\pi_i(y_i | \ldots)}{\pi_i(0 | \ldots)} = \left(\mu_i^{(i)}\right)^{y_i} / y_i! = \exp\left(y_i \cdot \theta \sum_{j \sim i, j < i} y_j\right) / y_i!$$

**Producto sobre todos los sitios:**
$$\frac{\pi(y)}{\pi(0)} = \frac{1}{\prod_i y_i!} \exp\left(\theta \sum_{\substack{i,j \in T \\ i \sim j, i < j}} y_i y_j\right)$$

Con $\pi(0) = e^{-|T|}$, obtenemos:

$$\boxed{\pi(y) \propto \frac{\exp\left(\theta \sum_{i \sim j} y_i y_j\right)}{\prod_{i \in T} y_i!}}$$

**Estructura de cliques:** Distribución de Gibbs con potencial de clique de tamaño 2 (pares vecinos).

#### Parte b) Condición de Sumabilidad

**Objetivo:** Demostrar que la distribución está bien definida ($Z < \infty$) si y solo si $\theta \leq 0$.

**Caso $\theta \leq 0$ (existencia):**

Como $y_i y_j \geq 0$ y $\theta \leq 0$, tenemos:
$$\exp\left(\theta \sum_{i \sim j} y_i y_j\right) \leq 1$$

Por tanto:
$$Z(\theta) = \sum_{y \in \mathbb{N}_0^T} \frac{\exp(\theta \sum_{i \sim j} y_i y_j)}{\prod_i y_i!} \leq \sum_{y} \frac{1}{\prod_i y_i!} = \prod_{i=1}^{|T|} e = e^{|T|} < \infty$$

✓ La distribución existe y está normalizada.

**Caso $\theta > 0$ (no existencia):**

Considérese la familia de configuraciones uniformes $y^{(\lambda)} = (\lambda, \lambda, \ldots, \lambda)$ para $\lambda \in \mathbb{N}$.

El número de pares vecinos es $|E|$, así que:
$$Z(\theta) \geq \sum_{\lambda=0}^\infty \frac{\exp(\theta |E| \lambda^2)}{(\lambda!)^{|T|}}$$

Por la fórmula de Stirling, $\lambda! \approx \sqrt{2\pi\lambda}(\lambda/e)^\lambda$, así:
$$(\lambda!)^{|T|} \approx \exp(|T| \lambda \log \lambda)$$

El numerador crece como $\exp(c\lambda^2)$ con $c = \theta|E| > 0$ (cuadrático en $\lambda$), mientras que el denominador crece solo como $\exp(c' \lambda \log \lambda)$ (subexponencial cuadrático).

Por tanto:
$$\frac{\exp(\theta |E| \lambda^2)}{(\lambda!)^{|T|}} \to +\infty \quad \text{cuando} \quad \lambda \to \infty$$

Los términos de la suma **no convergen a cero**, luego $Z(\theta) = +\infty$ y la distribución **no está normalizada**.

**Conclusión:**
$$\boxed{\pi(y) \text{ es una distribución válida} \iff \theta \leq 0}$$

#### Parte c) Interpretación Física para $\theta > 0$

**Mecanismo de inestabilidad:**

Cuando $\theta > 0$, la media condicional $\mu_i = \exp(\theta \sum_{j \sim i} y_j)$ es **estrictamente creciente** en los valores vecinos. Esto genera:

1. **Retroalimentación positiva:** Un valor alto en $y_i$ eleva $\mu_j$ para todos los $j \sim i$, incrementando la probabilidad de que $y_j$ sea grande, lo que a su vez incrementa $\mu_i$, sin ningún mecanismo de amortiguación.

2. **Regímenes según $\theta$:**

| Parámetro | Distribución | Comportamiento Espacial |
|-----------|--------------|------------------------|
| $\theta < 0$ | Bien definida | Competencia, suavidad, variación local (valores altos evitan sitios vecinos) |
| $\theta = 0$ | Bien definida (producto) | Independencia espacial: $Y_i \sim \text{Poisson}(1)$ |
| $\theta > 0$ | **No existe** | Retroalimentación positiva, $Z = +\infty$ |

3. **Formación de "clumps" (aglomeraciones):**

Intuitivamente, cuando $\theta > 0$, el sistema **favorece la coexistencia de valores simultáneamente altos en sitios vecinos**, generando grandes bloques de sitios adyacentes con conteos muy altos. Como los conteos no están acotados (a diferencia del modelo de Ising binario), la retroalimentación positiva es **matemáticamente irreparable**:

- La función de partición diverge porque la contribución de las configuraciones uniformes $y^{(\lambda)} = (\lambda, \ldots, \lambda)$ crece sin límite.
- El sistema **no tiene punto de equilibrio**: las aglomeraciones tienen incentivo a crecer indefinidamente, escapando al infinito.
- La distribución conjunta **no puede asignarse de manera coherente** y el modelo es inválido.

4. **Analogía con Ising ferromagnético:**

El comportamiento es análogo al modelo de Ising con interacciones ferromagnéticas ($\theta > 0$ en versión binaria): el sistema favorece la alineación completa. Pero en la versión Poisson, al no haber cota superior en los valores, la inestabilidad es **absoluta y matemática**, no meramente probabilística. La distribución conjunta no existe para ningún $\theta > 0$.

---

## Teoremas y Definiciones Clave

### Definiciones

- **MRF (Campo Aleatorio de Markov):** Un campo aleatorio que satisface la propiedad de Markov local.
- **Clique:** Subconjunto de sitios mutuamente vecinos.
- **CAR (Conditional Autoregressive):** Modelo especificado mediante distribuciones condicionales.
- **Matriz de Precisión:** $Q = K^{-1}(I-B)$ en modelos CAR gaussianos.

### Teoremas Fundamentales

1. **Hammersley-Clifford:** Equivalencia entre MRF y representación de Gibbs vía cliques.
2. **Lema de Brook (Factorización de Besag):** Cálculo de densidades conjuntas desde condicionales.
3. **Reciprocidad CAR:** Condición de simetría que asegura compatibilidad de condicionales gaussianas.

---

## Estructura de la Solución

El archivo `Ay_07_Sol.pdf` contiene:

- **Demostraciones rigurosas** de cada resultado
- **Cálculos detallados** de derivadas y expansiones algebraicas
- **Interpretaciones intuitivas** de conceptos matemáticos
- **Tablas comparativas** de comportamientos según parámetros

---

## Lecturas Complementarias

- Cressie, N. (1993). *Statistics for Spatial Data* (Revised Edition). Wiley. — Caps. 2–3.
- Waller, L. A., & Gotway, C. A. (2004). *Applied Spatial Statistics for Public Health Data*. Wiley. — Cap. 5.
- Besag, J. (1974). "Spatial Interaction and the Statistical Analysis of Lattice Systems." *Journal of the Royal Statistical Society: Series B*, 36(2).
- Hammersley, J. M., & Clifford, P. (1971). "Markov Fields on Finite Graphs and Lattices."

---

## Observaciones Pedagógicas

1. **Importancia de las condiciones de consistencia:** Los Problemas 2c y 3b demuestran que pequeñas variaciones paramétricas ($\theta$ de signo) pueden cambiar de forma **cualitativa** la validez del modelo.

2. **Estructura bipartita como simplificación:** El Problema 1 muestra cómo una disposición espacial particular (ajedrez) reduce la complejidad de cliques, ilustrando la interacción entre **topología** y **factorización probabilística**.

3. **De Poisson a Ising:** El Problema 3 justifica por qué la versión binaria del modelo es más robusta (acotada) que la versión de conteos.

---

**Contacto:** jipinov95@gmail.com  
**Última actualización:** Mayo 2026  
**Tipo de sesión:** Teórica (demostración y análisis)
