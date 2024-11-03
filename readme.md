# OpenMP Distribuido - Envío de mensajes

## Descripción

El presente proyecto tiene como objetivo extender la funcionalidad de la librería OpenMP para el lenguaje de programación C, la cual tradicionalmente es utilizada para la paralelización de aplicaciones en sistemas de memoria compartida. Esta extensión permite a esta librería aumentar su alcance de funcionamiento, permitiendo a su vez gestionar sistemas de memoria distribuida. Este proyecto surge de la necesidad de simplificar la programación en entornos distribuidos, dado que actualmente la curva de aprendizaje de esta tecnología es excesivamente desafiante. Mediante esta nueva funcionalidad, los desarrolladores podrán utilizar OpenMP para coordinar tareas en clústeres de computadoras de manera más accesible que con la herramienta MPI.
Se han añadido nuevas directivas a OpenMP con el objetivo de reconocer y gestionar la distribución de tareas y la comunicación entre nodos. Para lograr este objetivo, se ha integrado una gramática modificada de la librería OpenMP con la gramática del estándar C99 y se ha diseñado una Tabla de Símbolos para manejar la traducción de directivas OpenMP a MPI. Para probar esta funcionalidad, se han implementado y probado diversos ejemplos prácticos que demuestran la funcionalidad y eficiencia de esta ampliación.
Los resultados indican que la ampliación de OpenMP mejora significativamente la complejidad del desarrollo de aplicaciones paralelas, manteniendo la accesibilidad característica de OpenMP, facilitando la adopción de técnicas de programación paralela y distribuida

### Prerrequisitos 📋

- Sistema Operativo: cualquier distribución linux.

## Ejecutando las Pruebas ⚙️

```bash 
./exec.sh [-h] [-c]  [-i argumento] [-l argumento] [-e argumento] [-o argumento]
```

Donde las opciones corresponden a:
- *-h:* Muestra un mensaje de ayuda.
- *-c:* Compilar el código fuente antes de ejecutar.
- *-i:* Especifica un fichero de entrada distinto al de por defecto (input.c).
- *-l:* Especifica un fichero de logs distinto al de por defecto (log.txt).
- *-e:* Especifica un fichero de volcado de errores distinto al de por defecto (error.txt).
- *-o:* Especifica un fichero de salida distinto al de por defecto (output.c).


Breve explicación de los archivos:
- *Archivo input:* corresponde al archivo el cual se quiere traducir.
- *Arvhivo de logs:* corresponde al fichero donde se podrá ver los logs de la ejecución.
- *Archivo de volcado de errores:* corresponde al fichero donde se podrán ver los errores de ejecución en el caso de que se hayan producido.
- *Archivo output:* donde se podrá observar el resultado de la traducción.
