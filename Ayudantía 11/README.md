# Ayudantía 11: Procesos Puntuales — Teoría y Análisis Aplicado

**Curso:** EYP3417 Estadística Espacial  
**Profesor:** Alfredo Alegría  
**Ayudante:** Juan Pino  
**Período:** 2026  

---

## Descripción General

Esta ayudantía integra teoría rigurosa y análisis aplicado de **procesos puntuales de Poisson (PPP)**. Los dos primeros ejercicios desarrollan resultados matemáticos fundamentales — log-verosimilitud, distribuciones Palm, función K de Ripley y modelos hard-core — mientras el tercero implementa el pipeline completo de análisis sobre datos de trampas cámara simulados.

---

## Estructura

### Ejercicio 1 — PPP Inhomogéneo con $\lambda(x,y) = 1 + 4xy$

Sea $X$ un PPP inhomogéneo sobre $S = [0,1]^2$ con función de intensidad $\lambda(x,y) = 1 + 4xy$.

**a) Medida de intensidad**

$$\mu(S) = \int_0^1\int_0^1 (1+4xy)\,dx\,dy = 2$$

Se esperan en promedio **2 puntos** en $S$.

**b) Distribución de $N(S)$ y densidad condicional**

A partir de la densidad de Janossy $j_n(u_1,\ldots,u_n) = e^{-\mu(S)}\prod_i\lambda(u_i)$, se muestra que:

$$P(N(S)=n) = e^{-2}\frac{2^n}{n!} \quad \Rightarrow \quad N(S)\sim\text{Poisson}(2)$$

Dado $N(S)=n$, los puntos son i.i.d. con densidad $f(x,y) = \lambda(x,y)/\mu(S) = (1+4xy)/2$.

**c) Log-verosimilitud y estadísticos suficientes**

Para el modelo $\log\lambda(x,y;\theta) = \theta_1 + \theta_2 x + \theta_3 y + \theta_4 xy$:

$$\ell(\theta) = n\theta_1 + \theta_2\sum x_i + \theta_3\sum y_i + \theta_4\sum x_iy_i - \int_S e^{\theta_1+\theta_2 x+\theta_3 y+\theta_4 xy}\,dx\,dy$$

El vector suficiente es $\mathbf{T} = \bigl(n,\,\sum x_i,\,\sum y_i,\,\sum x_iy_i\bigr)$.

---

### Ejercicio 2 — Modelo de Matérn-I: Hard-Core y Función $K$

Sea $X_0$ un PPP($\lambda_0$) sobre $\mathbb{R}^2$. Se construye $X_1$ eliminando todo punto con al menos un vecino a distancia $\leq r$.

**a) Intensidad de $X_1$**

Por el teorema de Slivnyak-Mecke, dado $\xi\in X_0$, el proceso restante sigue siendo PPP($\lambda_0$). La probabilidad de supervivencia al thinning es:

$$p_{\text{surv}} = P(N(B(\xi,r)\setminus\{\xi\})=0) = e^{-\lambda_0\pi r^2}$$

$$\Rightarrow \quad \lambda_1 = \lambda_0\,e^{-\lambda_0\pi r^2}$$

**b) Función $K$ de Ripley para $X_1$**

Usando $\lambda_1^2 K(h) = 2\pi\int_0^h t\,\lambda_2(t)\,dt$ con $\lambda_2(h) = \lambda_0^2\exp\{-\lambda_0 V_r(h)\}\cdot\mathbf{1}(h>r)$:

$$K(h) = \begin{cases} 0 & h \leq r \\ 2\pi e^{2\lambda_0\pi r^2}\displaystyle\int_r^h t\,e^{-\lambda_0 V_r(t)}\,dt & r < h \leq 2r \\ \pi h^2 - 4\pi r^2 + C_r & h > 2r \end{cases}$$

donde $V_r(h) = 2\pi r^2 - 2r^2\arccos\!\tfrac{h}{2r} + \tfrac{h}{2}\sqrt{4r^2-h^2}$ y $C_r > 0$ es una constante.

**c) Comportamiento cualitativo de $L(h)-h$**

| Rango | $L(h)-h$ | Interpretación |
|---|---|---|
| $h < r$ | $= -h$ (mínimo en $h=r$: $-r$) | Ausencia total de pares — zona hard-core |
| $r < h \leq 2r$ | $< 0$, creciente | Inicio de pares, aún bajo el nivel CSR |
| $h \gg r$ | $\to 0^-$ | Recuperación asintótica al nivel de CSR |

La curva es **siempre negativa**, firma característica de un proceso con **inhibición de corto alcance** más regular que un PPP.

---

### Ejercicio 3 — Análisis Aplicado: Trampas Cámara en Bosque Templado

Se simulan $n = 216$ avistamientos de fauna en una ventana de $2\times 2$ km bajo el modelo $\log\lambda(x,y) = -9.0 - 0.001\cdot x$, que concentra la actividad faunística hacia el borde oeste.

#### a) Patrón de puntos

El patrón de $n=216$ puntos muestra concentración visible hacia el borde oeste, con densidad decreciente hacia el este, anticipando un PPP no homogéneo con tendencia en $x$.

#### b) Estimación kernel ($\sigma_{\text{Diggle}} = 121.8$ m)

| Bandwidth | Efecto |
|---|---|
| $\sigma/5 = 24.4$ m | Sobreajuste — picos aislados en cada punto (alta varianza) |
| $\sigma = 121.8$ m | Balance óptimo — recupera el gradiente oeste-este |
| $5\sigma = 609.1$ m | Sobresuavizado — gradiente diluido, superficie casi constante |

#### c) Tests G y F (99 simulaciones CSR, nivel efectivo $= 0.02$)

Ambas curvas permanecen dentro de las bandas de confianza. No se rechaza CSR a escala local. Esto es esperable para un PPP no homogéneo: G y F bajo CSR homogéneo no son sensibles a la heterogeneidad de primer orden a las escalas de vecindad más próximas.

#### d) Función $L(h)-h$ bajo CSR

La curva está dentro de las bandas para $h \lesssim 100$ m. A partir de $h \approx 150$ m sale progresivamente, alcanzando $L(h)-h \approx 65$ en $h = 500$ m. Esta desviación es un **artefacto de primer orden**: el envelope bajo CSR homogéneo no contempla la variación de $\lambda(x,y)$, sobreestimando los pares esperados a gran escala.

#### e) Modelo log-lineal y validación Monte Carlo

| Parámetro | $\hat\theta$ | Valor real | SE | IC 95% | Wald |
|---|---|---|---|---|---|
| $\theta_1$ | $-8.9285$ | $-9.0000$ | $0.1107$ | $[-9.146,\,-8.711]$ | $***$ |
| $\theta_2$ | $-0.0011$ | $-0.0010$ | $0.0001$ | $[-0.00135,\,-0.00083]$ | $***$ |

Ambos estimadores recuperan los valores poblacionales. $\hat\theta_2 < 0$ significativo (Z = -8.26, p < 0.001): cada 100 m adicionales en $x$ reducen $\lambda$ en un factor $\exp(100\hat\theta_2) \approx 0.897$ (~10% menos avistamientos).

Las bandas de Monte Carlo bajo el modelo ajustado muestran la curva observada siguiendo la media simulada y permaneciendo dentro de las bandas en todo el rango. **No se rechaza el modelo log-lineal**: la estructura espacial queda completamente explicada por la heterogeneidad de primer orden en $x$, sin necesidad de interacción de segundo orden.

---

## Archivos

| Archivo | Contenido |
|---|---|
| [`Ay_11.pdf`](./Ay_11.pdf) | Enunciado de la ayudantía |
| [`Ay_11_Espacial.R`](./Ay_11_Espacial.R) | Script R — Ejercicio 3 completo |
| [`Ayudantía 11 - Solución ejercicios 1 y 2.md`](./Ayudantía%2011%20-%20Solución%20ejercicios%201%20y%202.md) | Solución teórica detallada — Ejercicios 1 y 2 |

---

## Librerías Requeridas

```r
library(spatstat)     # ppp, owin, rpoispp, density.ppp, Lest, Gest, Fest, ppm, envelope
library(viridisLite)  # Paleta plasma para mapas de intensidad
```

---

## Lecciones Clave

1. **La densidad de Janossy conecta la definición del PPP con la inferencia:** la log-verosimilitud y la suficiencia emergen directamente de la factorización $j_n = e^{-\mu}\prod\lambda(u_i)$.

2. **Slivnyak-Mecke es la herramienta central para el PPP:** permite calcular probabilidades condicionadas a la presencia de un punto sin alterar la distribución del resto del proceso.

3. **El modelo Matérn-I produce $L(h)-h$ siempre negativa:** la inhibición hard-core se firma con un mínimo profundo en $h=r$ y recuperación asintótica a cero — patrón opuesto al clustering.

4. **G y F no detectan heterogeneidad de primer orden:** al comparar contra CSR homogéneo, pueden no rechazar la hipótesis nula aunque el patrón sea claramente no homogéneo. La función $L(h)-h$ bajo CSR sí la detecta, pero de forma espuria.

5. **El envelope bajo el modelo ajustado es el diagnóstico definitivo:** si la curva observada entra en las bandas del modelo log-lineal, la heterogeneidad de primer orden explica todo y no hay evidencia de interacción entre puntos.

---

**Contacto:** jipinov95@gmail.com  
**Última actualización:** Junio 2026  
**Tipo de sesión:** Teórico-práctica — PPP inhomogéneo, modelo Matérn-I y análisis de trampas cámara
