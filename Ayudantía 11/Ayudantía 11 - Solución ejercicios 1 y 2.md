
## Ejercicio 1 — PPP Inhomogéneo con $\lambda(x,y) = 1 + 4xy$

## Ejercicio 1

Sea $X$ un PPP no homogéneo sobre la ventana $S = [0,1]^2$ con función de intensidad $\lambda(x,y) = 1 + 4xy$.

### Parte a

Calcule $\mu(S) = \int_S \lambda(x,y)\,dx\,dy$. ¿Cuántos puntos se esperan en promedio en $S$?

**Solución**: 

La medida de intensidad sobre la ventana completa es:

$$\mu(S) = \int_0^1 \int_0^1 (1 + 4xy)\,dx\,dy$$

Integrando primero respecto a $x$ :  

$$\int_0^1 (1 + 4xy)\,dx = \left[x + 2x^2 y\right]_0^1 = 1 + 2y$$

Luego, sustituyendo en la integral anterior se tiene:

$$\mu(S) = \int_0^1 (1 + 2y)\,dy = \left[y + y^2\right]_0^1 = 1 + 1 = 2$$
Por la propiedad fundamental del PPP, tenemos que  $\mathbb{E}[N(S)] = \mu(S)$, por lo tanto se esperan en promedio **2 puntos** en $S$. $\blacksquare$


### Parte b

 Verifique que $N(S) \sim \text{Poisson}(\mu(S))$ y que los puntos, dado $N(S) = n$, son i.i.d. con densidad $f(x,y) = \lambda(x,y)/\mu(S)$.

**Solución**:

Para un PPP no homogéneo con intensidad $\lambda(\cdot)$, la densidad de Janossy (densidad conjunta de la configuración no ordenada de $n$ puntos) respecto a la medida de Lebesgue $n$-dimensional es:

$$j_n(u_1, \ldots, u_n) = \displaystyle e^{-\mu(S)} \prod_{i=1}^n \lambda(u_i)$$

Esta expresión viene directamente de la definición del PPP, la probabilidad de observar exactamente un punto en cada región infinitesimal $du$ y ningún punto en el resto se factoriza en el producto de las intensidades locales, pesado por la función de vacío $e^{-\mu(S)}$.

Luego, integrando la densidad de Janossy sobre todas las configuraciones de $n$ puntos en $S$:

$$P(N(S) = n) = \int_{S^n} j_n(u_1,\ldots,u_n)\,du_1\cdots du_n = e^{-\mu(S)} \int_{S^n} \prod_{i=1}^n \lambda(u_i)\,du_1\cdots du_n$$
Por independencia de las integrales:

$$P(N(S) = n) = e^{-\mu(S)} \left(\int_S \lambda(u)\,du\right)^n = e^{-\mu(S)} \frac{\mu(S)^n}{n!} \cdot n!$$

Dividiendo por $n!$ para corregir la simetría (puntos no ordenados):

  

$$P(N(S) = n) = e^{-\mu(S)} \frac{\mu(S)^n}{n!}$$

que es precisamente la distribución $\text{Poisson}(\mu(S))$. Con $\mu(S) = 2$, se tiene $N(S) \sim \text{Poisson}(2)$. $\blacksquare$

Ahora, probemos que son i.i.d. Para ello veamos que la densidad condicional de la configuración no ordenada $(u_1, \ldots, u_n)$ dado $N(S) = n$ es:

$$f(u_1, \ldots, u_n \mid N(S) = n) = \frac{j_n(u_1, \ldots, u_n)}{P(N(S) = n)} = \frac{e^{-\mu(S)} \prod_{i=1}^n \lambda(u_i)}{e^{-\mu(S)}\,\mu(S)^n / n!} \cdot \frac{1}{n!}$$

Simplificando tenemos:
$$f(u_1, \ldots, u_n \mid N(S) = n) = \frac{\prod_{i=1}^n \lambda(u_i)}{\mu(S)^n} =
\prod_{i=1}^n \frac{\lambda(u_i)}{\mu(S)} = \prod_{i=1}^n f(u_i)$$

donde $f(u) = \lambda(u)/\mu(S) = (1 + 4xy)/2$. Dado que la densidad conjunta se factoriza como producto de densidades marginales idénticas, los $n$ puntos son **i.i.d.** con densidad $f(x,y) = (1+4xy)/2$. $\blacksquare$


### Parte c

Suponga que se observan $n$ puntos $\{(x_i, y_i)\}_{i=1}^n$. Escriba la log-verosimilitud del modelo $\log\lambda(x,y;\theta) = \theta_1 + \theta_2 x + \theta_3 y + \theta_4 xy$. Identifique los estadísticos suficientes.

**Solución**:

Para un PPP no homogéneo con intensidad $\lambda(u;\theta)$ observado en la ventana $S$, la log-verosimilitud es: 

$$\ell(\theta) = \sum_{i=1}^n \log \lambda(u_i;\theta) - \int_S \lambda(u;\theta)\,du$$
El primer término contribuye de manera positiva donde se observaron puntos (favorece intensidades altas en los datos), mientras que el segundo penaliza modelos que sobreestiman la intensidad total

Con $\log \lambda(x,y;\theta) = \theta_1 + \theta_2 x + \theta_3 y + \theta_4 xy$, tenemos:

$$\sum_{i=1}^n \log \lambda(x_i, y_i;\theta) = \sum_{i=1}^n (\theta_1 + \theta_2 x_i + \theta_3 y_i + \theta_4 x_i y_i)$$
$$\sum_{i=1}^n \log \lambda(x_i, y_i;\theta) = n\theta_1 + \theta_2 \sum_{i=1}^n x_i + \theta_3 \sum_{i=1}^n y_i + \theta_4 \sum_{i=1}^n x_i y_i$$

Ahora, desarrollando el término de la integral con $\lambda(x,y;\theta) = \exp(\theta_1 + \theta_2 x + \theta_3 y + \theta_4 xy)$:

$$\int_S \lambda(x,y;\theta)\,dx\,dy = e^{\theta_1} \int_0^1 e^{\theta_3 y} \left(\int_0^1 e^{(\theta_2 + \theta_4 y)x}\,dx\right) dy$$

  Para $\theta_2 + \theta_4 y \neq 0$, la integral interior es:

$$\int_0^1 e^{(\theta_2 + \theta_4 y)x}\,dx = \frac{e^{\theta_2 + \theta_4 y} - 1}{\theta_2 + \theta_4 y}$$

Por lo que el término integral queda:

$$\int_S \lambda(x,y;\theta)\,dx\,dy = e^{\theta_1} \int_0^1 e^{\theta_3 y} \cdot \frac{e^{\theta_2 + \theta_4 y} - 1}{\theta_2 + \theta_4 y}\,dy$$


Esta integral no admite forma cerrada elemental en general (depende de $\theta_4 y$ de manera no lineal), y es evaluada numéricamente por `ppm()`.

Uniendo todo, tenemos:

$$\boxed{\ell(\theta) = n\theta_1 + \theta_2 \sum_{i=1}^n x_i + \theta_3 \sum_{i=1}^n y_i + \theta_4 \sum_{i=1}^n x_i y_i - e^{\theta_1}\int_0^1 e^{\theta_3 y}\,\frac{e^{\theta_2 + \theta_4 y} - 1}{\theta_2 + \theta_4 y}\,dy}$$

El primer término de $\ell(\theta)$ tiene la forma $\theta^\top T(\mathbf{x})$, característica de la familia exponencial. Los parámetros $\theta = (\theta_1, \theta_2, \theta_3, \theta_4)$ aparecen linealmente en:


$$\ell(\theta) = \underbrace{n}_{T_1}\theta_1 + \underbrace{\textstyle\sum x_i}_{T_2}\theta_2 + \underbrace{\textstyle\sum y_i}_{T_3}\theta_3 + \underbrace{\textstyle\sum x_i y_i}_{T_4}\theta_4 - C(\theta)$$

donde $C(\theta) = \displaystyle \int_S e^{\theta_1+\theta_2 x+\theta_3 y+\theta_4 xy}\,dx\,dy$ depende solo de $\theta$. Por el criterio de factorización de Fisher-Neyman, el vector

$$\mathbf{T} = \left(n,\;\sum_{i=1}^n x_i,\;\sum_{i=1}^n y_i,\;\sum_{i=1}^n x_i y_i\right)$$

es **suficiente** para $\theta$. Cada componente tiene una interpretación directa: $T_1 = n$ determina el "nivel" de la intensidad, $T_2$ y $T_3$ capturan la tendencia lineal en $x$ e $y$, mientras que $T_4 = \sum x_i y_i$ captura la interacción entre coordenadas. $\blacksquare$

  

---
## Ejercicio 2 

Sea $X_0$ un PPP homogéneo con intensidad $\lambda_0$ sobre $\mathbb{R}^2$ y sea $r > 0$. Se construye el proceso $X_1$ eliminando todo punto de $X_0$ que tenga al menos un vecino a distancia $\leq r$.

### Parte a
 
Muestre que $\lambda_1 = \lambda_0\exp\{-\lambda_0\pi r^2\}$.

**Solución**:  Un punto $\xi \in X_0$ sobrevive al thinning si y solo si ningún otro punto de $X_0$ cae en la bola $B(\xi, r)$. La probabilidad de supervivencia es:

$$p_{\text{surv}} = P\bigl(N(B(\xi, r) \setminus \{\xi\}) = 0 \;\big|\; \xi \in X_0\bigr)$$

Por el **teorema de Slivnyak-Mecke** para el PPP. Dado que $\xi \in X_0$, el proceso restante $X_0 \setminus \{\xi\}$ sigue siendo un PPP$(\lambda_0)$ independiente de la posición de $\xi$. Por lo tanto:

$$N\bigl(B(\xi, r) \setminus \{\xi\}\bigr) \;\big|\; \xi \in X_0 \;\sim\; \text{Poisson}\bigl(\lambda_0 \cdot |B(\xi,r)|\bigr) = \text{Poisson}(\lambda_0 \pi r^2)$$
De lo anterior, la probabilidad de vacío es:
  
$$p_{\text{surv}} = P\bigl(N(B(\xi,r)\setminus\{\xi\}) = 0\bigr) = e^{-\lambda_0 \pi r^2}$$

El proceso $X_1$ es un **thinning independiente** de $X_0$: cada punto se retiene con probabilidad $p_{\text{surv}}$ que no depende de $\xi$ (homogeneidad de $X_0$). Por la propiedad de superposición del PPP bajo thinning independiente homogéneo, $X_1$ es un proceso puntual estacionario con intensidad:

$$\boxed{\lambda_1 = \lambda_0 \cdot p_{\text{surv}} = \lambda_0\, e^{-\lambda_0 \pi r^2}}$$  

**Observación:** A diferencia de un PPP, $X_1$ **no es** un PPP (sus puntos no son independientes), pero la intensidad de primer orden es bien definida y viene dada por la expresión anterior. La función $\lambda_0 \mapsto \lambda_1(\lambda_0)$ es no monótona: alcanza su máximo en $\lambda_0 = 1/(\pi r^2)$, con $\lambda_1^{\max} = 1/(e\pi r^2)$. $\blacksquare$

### Parte b 

Deduzca la función $K$ de Ripley para $X_1$ a partir de la intensidad de segundo orden dada.

**Solución**: 

Para un proceso puntual estacionario e isótropo en $\mathbb{R}^2$ con intensidad $\lambda$ e intensidad de segundo orden $\lambda_2(h)$ (función solo de la distancia $h = \|u - v\|$), la función $K$ de Ripley satisface:  

$$\lambda^2 K(h) = 2\pi \int_0^h t\,\lambda_2(t)\,dt$$

o equivalentemente, diferenciando:

$$\lambda^2 K'(h) = 2\pi h\,\lambda_2(h)$$
Se tiene por enunciado:

$$\lambda_2(h) = \lambda_0^2\,\exp\{-\lambda_0 V_r(h)\} \cdot \mathbf{1}(h > r)$$
donde $V_r(h)$ es el área de la unión de dos discos de radio $r$ separados una distancia $h$:

$$V_r(h) = \begin{cases} 2\pi r^2 - 2r^2 \arccos\!\left(\dfrac{h}{2r}\right) + \dfrac{h}{2}\sqrt{4r^2 - h^2}, & 0 \leq h \leq 2r \\[6pt] 2\pi r^2, & h > 2r \end{cases}$$

  La restricción $\mathbf{1}(h > r)$ refleja el hard-core: no pueden existir pares de puntos en $X_1$ a distancia menor que $r$.

Con $\lambda = \lambda_1 = \lambda_0 e^{-\lambda_0\pi r^2}$, se tiene $\lambda_1^2 = \lambda_0^2\, e^{-2\lambda_0\pi r^2}$. Aplicando la fórmula del primer paso tenemos:

  

$$K(h) = \frac{2\pi}{\lambda_1^2}\int_0^h t\,\lambda_2(t)\,dt = \frac{2\pi}{\lambda_0^2 e^{-2\lambda_0\pi r^2}}\int_0^h t\,\lambda_0^2 e^{-\lambda_0 V_r(t)}\,\mathbf{1}(t>r)\,dt$$

$$= 2\pi\, e^{2\lambda_0\pi r^2} \int_r^h t\,e^{-\lambda_0 V_r(t)}\,dt \qquad \text{(para } h > r\text{)}$$

**Tramo 1: $0 < h \leq r$.**

No existen pares a distancia $\leq r$, luego $\lambda_2(t) = 0$ para todo $t \in (0, r]$:  

$$\boxed{K(h) = 0, \qquad 0 < h \leq r}$$

**Tramo 2: $r < h \leq 2r$.**  

$$K(h) = 2\pi\, e^{2\lambda_0\pi r^2} \int_r^h t\,\exp\!\left\{-\lambda_0\!\left(2\pi r^2 - 2r^2\arccos\!\tfrac{t}{2r} + \tfrac{t}{2}\sqrt{4r^2-t^2}\right)\right\} dt$$

Esta integral no tiene forma cerrada elemental y se evalúa numéricamente.


**Tramo 3: $h > 2r$.**

Para $t > 2r$, los dos discos no se solapan y $V_r(t) = 2\pi r^2$, por lo que:

$$K(h) = 2\pi\, e^{2\lambda_0\pi r^2} \left[\int_r^{2r} t\,e^{-\lambda_0 V_r(t)}\,dt + \int_{2r}^h t\,e^{-\lambda_0 \cdot 2\pi r^2}\,dt\right]$$

$$= 2\pi\, e^{2\lambda_0\pi r^2} \int_r^{2r} t\,e^{-\lambda_0 V_r(t)}\,dt + 2\pi \int_{2r}^h t\,dt$$
$$= 2\pi\, e^{2\lambda_0\pi r^2} \int_r^{2r} t\,e^{-\lambda_0 V_r(t)}\,dt + \pi(h^2 - 4r^2)$$
 

Definiendo la constante $C_r = 2\pi e^{2\lambda_0\pi r^2}\int_r^{2r} t\,e^{-\lambda_0 V_r(t)}\,dt > 0$, la expresión compacta para $h > 2r$ es:

$$\boxed{K(h) = \pi h^2 - 4\pi r^2 + C_r, \qquad h > 2r}$$

**Verificación asintótica:** Para $h \to \infty$, el término dominante es $\pi h^2$, que coincide con $K_{\text{PPP}}(h)$. Esto es consistente con que a grandes escalas, la regularidad de corto alcance del modelo Matérn-I se vuelve irrelevante y el proceso se asemeja a un PPP. $\blacksquare$

### Parte c

¿Cómo se ve cualitativamente $L(h) - h$ para $X_1$? Describa el comportamiento para $h < r$, $h \approx r$ y $h \gg r$, y explique por qué $X_1$ es más regular que un PPP.

**Solución**: La función $L$ se define como la transformación estabilizadora de $K$:  

$$L(h) = \sqrt{\frac{K(h)}{\pi}}$$

de modo que bajo CSR, $L(h) = h$ y $L(h) - h = 0$ para todo $h$. Analizamos $X_1$ por tramos.

  **Tramo 1: $h < r$ (zona hard-core).**

$$K(h) = 0 \implies L(h) = 0 \implies L(h) - h = -h < 0$$

La curva cae linealmente desde el origen hasta el valor $-r$ en $h = r$. La ausencia total de pares a distancia $< r$ produce la mayor desviación negativa posible respecto al PPP.

**Tramo 2: $h \approx r$ (inicio de la zona de transición).**

Para $h$ ligeramente mayor que $r$, $K(h)$ comienza a crecer desde cero conforme se integran los primeros pares. Sin embargo, $K(h) < \pi h^2$ (hay menos pares que bajo CSR), por lo que:

$$L(h) = \sqrt{K(h)/\pi} < h \implies L(h) - h < 0$$
La curva tiene su mínimo absoluto en $h = r$ (donde vale $-r$) y comienza a ascender para $h > r$, pero permanece por debajo de cero.

**Tramo 3: $h \gg r$ (régimen asintótico).**

Para $h > 2r$, $K(h) = \pi h^2 - 4\pi r^2 + C_r$, por lo que:

$$L(h) = \sqrt{h^2 - 4r^2 + C_r/\pi} \approx h\sqrt{1 - \frac{4r^2 - C_r/\pi}{h^2}} \approx h - \frac{4r^2 - C_r/\pi}{2h} + O(h^{-3})$$
Luego $L(h) - h \to 0^-$ cuando $h \to \infty$ (asintóticamente cero por abajo, salvo que $C_r > 4\pi r^2$, lo que no ocurre en general).

**Perfil cualitativo:**

  

$$L(h) - h \begin{cases} = -h < 0, & h < r \quad \text{(descenso lineal)} \\ = -r \quad \text{(mínimo global)}, & h = r \\ < 0 \text{ y creciente}, & r < h \leq 2r \\ \to 0^-, & h \to \infty \end{cases}$$

La curva es **siempre negativa**, con mínimo en $h = r$ y recuperación asintótica hacia cero.

**¿Por qué $X_1$ es más regular que un PPP?**

En un PPP no hay interacción entre puntos: la presencia de un punto en $\xi$ no altera la distribución del resto del proceso (propiedad de Slivnyak). En $X_1$, en cambio, existe una **inhibición de corto alcance**: ningún par de puntos puede estar a distancia $< r$. Esta restricción obliga a los puntos a estar más separados de lo que estarían bajo CSR, produciendo:

1. **Ausencia de pares cercanos** ($K(h) = 0$ para $h < r$): no hay cúmulos a escala pequeña.

2. **Curva $L(h)-h$ negativa en todo su rango**: siempre hay menos pares acumulados que bajo CSR.

3. **Patrón visualmente más uniforme**: los puntos "llenan" el espacio de manera más homogénea que bajo CSR, sin zonas densas ni vacíos extremos.

En definitiva, $L(h) - h < 0$ para todo $h > 0$ es la firma característica de la **regularidad** o **inhibición**, y el modelo Matérn-I es el ejemplo paradigmático de un proceso hard-core con esta propiedad. $\blacksquare$

  



  