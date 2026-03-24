# Cómo acceder y mantenerse actualizado con el repositorio del curso

## 1. Instalar Git

Descarga e instala Git desde **https://git-scm.com/downloads**, elige tu sistema operativo y sigue el instalador con las opciones por defecto.

Para verificar que quedó instalado, abre una terminal (CMD en Windows, Terminal en Mac/Linux) y ejecuta:

```bash
git --version
```

Debería aparecer algo como `git version 2.x.x`.

---

## 2. Clonar el repositorio (solo la primera vez)

Clonar significa descargar una copia completa del repositorio en tu computador. Abre la terminal, navega a la carpeta donde quieras guardar los materiales y ejecuta:

```bash
git clone https://github.com/jipino/EYP3417-Estadistica-Espacial.git
```

Esto creará una carpeta llamada `EYP3417-Estadistica-Espacial` con todo el contenido del curso.

---

## 3. Actualizar cada semana (el comando más importante)

Cada vez que se suba material nuevo — script, enunciado, README — solo debes entrar a la carpeta del repositorio y ejecutar:

```bash
cd EYP3417-Estadistica-Espacial
git pull
```

Eso descarga automáticamente todo lo que se haya agregado desde la última vez que actualizaste. No necesitas volver a clonar.

---

## 4. Ver qué cambió

Si quieres saber qué archivos se agregaron o modificaron en la última actualización:

```bash
git log --oneline -5
```

Muestra los últimos 5 commits con su descripción.

```bash
git diff HEAD~1 --name-only
```

Muestra exactamente qué archivos cambiaron respecto a la versión anterior.

---

## Flujo semanal resumido

```
Semana 1:             git clone https://github.com/jipino/EYP3417-Estadistica-Espacial.git
Semana 2 en adelante: git pull
```
