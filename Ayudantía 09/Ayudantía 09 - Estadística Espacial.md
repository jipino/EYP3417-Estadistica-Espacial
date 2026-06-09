
### **Solución Detallada de Ejercicios Teóricos**

#### **Ejercicio 1a: Definición y propiedades teóricas del Proceso de Poisson**

De acuerdo con el libro guía, si definimos formalmente un PPP (Proceso Puntual de Poisson), primero se define $\lambda$ como una medida positiva sobre los conjuntos de Borel $\mathcal{B}(S)$ con densidad $\rho$.

Un PPP con medida de intensidad $\lambda(\cdot)$ > 0 e intensidad $\rho(\cdot)$ se caracteriza por las siguientes propiedades fundamentales:

- Distribución de conteos: Para cualquier conjunto de Borel acotado $A \in \mathcal{B}_b(S)$ cuya medida cumpla con que $0 < \lambda(A) < \infty$, la V.A de conteo $N(A)$ tiene una distribución de Poisson con parámetro $\lambda(A)$ 
- Independencia y distribución espacial: Condicional a que el número de eventos $N(A) = n$, los puntos que componen la configuración $x \cap A$ son i.i.d, con una densidad de probabilidad que es proporcional a la función de intensidad $\rho(\xi)$ para cada $\xi \in A$. Esta propiedad implica que para conjuntos de Borel disjuntos, los conteos son variables aleatorias independientes. 

Diferencia entre PPP Homogéneo y no homogéneo:

La diferencia radica netamente sobre la estructura de su función de intensidad $\rho(\xi)$ 

- Se dice que un PPP es homogéneo si su intensidad es constante, es decir, su medida de intensidad es proporcional a la medida de Lebesgue $\nu$, obteniéndose $\lambda(\cdot) = \rho \nu(\cdot)$. En este caso, la uniformidad de la distribución espacial se suma a la independencia.
- Por otra parte, un PPP es no homogéneo si la densidad $\rho(\xi)$ no es constante y varía en el espacio geográfico. Un ejemplo típico, según el texto, es modelar $\rho(\xi)$ mediante un modelo log-lineal que depende de covariables locales  $z(\xi)$, tal que $\log \rho(\xi) = z(\xi)^T \beta$ . Esto introduce heterogeneidad o tendencia de primer orden en el patrón de puntos.

---

#### **Ejercicio 2a: Análisis de Dependencia de Segundo Orden**

La función $K$ evalúa la dependencia espacial empírica estandarizando la medida del momento factorial de segundo orden. Conceptualmente, para un PPP homogéneo, "$K(h)$ es la esperanza del número de otros puntos del proceso en una bola de radio $h$ centrada en un punto de $X$, normalizada por la intensidad $\rho$".

**Nota**: Recordemos que CSR representa a una hipótesis espacial en procesos puntuales, la cual significa que los puntos se distribuyen de manera independiente unos de otros en el espacio, aunque no necesariamente de manera uniforme.

La varianza de la función $K(h)$ tiende a volverse inestable para radios $h$ grandes. Para corregir esto y linealizar la relación teórica, se introduce la función $L$, la cual es una estadística de resumen útil para verificar la hipótesis de CSR (Completa Aleatoriedad Espacial). Para un espacio de dimensión $d$ se define como $L(h) = \{K(h)/b_d\}^{1/d}$ . Dado que estamos en $\mathbb{R}^2$, el volumen de la esfera de $r=1$ es $b_2 = \pi$, por lo que nuestra ecuación se reduce a $L(h) = \sqrt{K(h)/\pi}$

Bajo el escenario nulo de CSR (un PPP homogéneo), se cumple teóricamente que $K(h) = \pi h^2$. De esta forma, bajo CSR, el gráfico de $h \to L(h) = h$  es una línea recta, lo que equivale a decir que $\hat{L}(h) - h = 0$

Basándonos en las propiedades geométricas de esta curva, los criterios empíricos que justifican el rechazo de la independencia espacial son:

- Clustering: Si la función $L$ es cóncava y se ubica por encima de la recta teórica ($\hat{L}(h) - h > 0$) esto indica patrones con agregados. Hay un exceso de pares de puntos a esa distancia en comparación al azar.
- Regularidad (Repulsión): Si la función $L$ es convexa y se ubica por debajo de la recta ($\hat{L}(h) - h < 0$), esto indica patrones más regulares que los de un PPP. Existe un déficit de puntos a esa distancia forzando distanciamiento.

---

#### **Ejercicio 3a: Formulación de la Log-Verosimilitud Paramétrica**

Si un proceso $X$ es un PPP no homogéneo con una intensidad $\rho(\cdot;\theta)$, parametrizada por un vector $\theta$, el EMV se obtiene maximizando la log-densidad espacial sobre la ventana o dominio de observación $A$

La log-verosimilitud observada para los parámetros de la intensidad toma de forma explícita la siguiente forma:
$$l_A(\theta) = \sum_{\xi \in x \cap A} \log \rho(\xi;\theta) - \int_A \rho(\eta;\theta) d\eta$$

**Nota**: Computacionalmente, como el término aditivo $\int_A 1 d\eta = \nu(A)$ es una constante, que representa el área de la ventana y no depende de $\theta$ suele eliminarse al derivar o maximizar, pero teóricamente la ecuación incorpora este término.

El problema establece un modelo log-lineal donde la intensidad depende de una sola covariable espacial (la coordenda $x$), dada por la ecuación:
$$\log \rho(\xi;\theta) = \theta_1 + \theta_2 x_{\xi}$$
Despejando, la función de intensidad directa es:
$$\rho(\xi;\theta) = \exp(\theta_1 + \theta_2 x_{\xi})$$

Luego, sustituyendo ambas expresiones directamente en la ecuación de la log-verosimilitud, la justificación matemática de la solución teórica exacta es:

$$l_A(\theta) = \sum_{\xi \in x \cap A} (\theta_1 + \theta_2 x_{\xi}) - \int_A \exp(\theta_1 + \theta_2 \eta_x) d\eta$$
Esta expresión es la función objetivo que los algoritmos de optimización maximizan para encontrar $\hat{\theta}_1$ y $\hat{\theta}_2$.