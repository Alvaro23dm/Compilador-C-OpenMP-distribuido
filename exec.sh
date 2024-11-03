#!/bin/bash

mostrar_ayuda() {
    echo "Uso: $0 [-h]  [-i argumento] [-l argumento] [-e argumento] [-o argumento]"
    echo
    echo "Opciones:"
    echo "  -h                Muestra este mensaje de ayuda"
	echo "  -c                Compilar el código fuente antes de ejecutar"
    echo "  -i argumento      Especifica un fichero de entrada distinto al de por defecto (input.c)"
    echo "  -l argumento      Especifica un fichero de logs distinto al de por defecto (log.txt)"
	echo "  -e argumento      Especifica un fichero de volcado de errores distinto al de por defecto (error.txt)"
	echo "  -o argumento      Especifica un fichero de salida distinto al de por defecto (output.c)"
}

ENTRADA=input.c
SALIDA=output.c
LOGS=log.txt
ERRORES=error.txt

while getopts ":hci:l:e:o:" opcion; do
    case $opcion in
        h)  # Opción -h (ayuda)
            mostrar_ayuda
            exit 0
            ;;
		c)  # Opción -c (compilación)
            echo COMPILACIÓN:
			make clean
			make
			echo
			echo
            ;;
        i)  # Opción -i (con argumento)
            if [ -z "$OPTARG" ]; then
    			echo "La opción -i no se proporcionó o no tiene argumento."
				exit 1
			fi
			ENTRADA=$OPTARG
			;;
        l)  # Opción -l (con argumento)
            if [ -z "$OPTARG" ]; then
    			echo "La opción -l no se proporcionó o no tiene argumento."
				exit 1
			fi
			LOGS=$OPTARG
			;;
        e)  # Opción -e (con argumento)
            if [ -z "$OPTARG" ]; then
    			echo "La opción -e no se proporcionó o no tiene argumento."
				exit 1
			fi
			ERRORES=$OPTARG
			;;
        o)  # Opción -o (con argumento)
            if [ -z "$OPTARG" ]; then
    			echo "La opción -o no se proporcionó o no tiene argumento."
				exit 1
			fi
			SALIDA=$OPTARG
			;;
        \?) # Opción no reconocida
            echo "Opción inválida: -$OPTARG" >&2
            mostrar_ayuda
            exit 1
            ;;
        :)  # Opción sin argumento requerido
            echo "La opción -$OPTARG requiere un argumento" >&2
            mostrar_ayuda
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))

echo EJECUCIÓN:
echo ./fparse $ENTRADA $LOGS $ERRORES $SALIDA
./fparse $ENTRADA $LOGS $ERRORES $SALIDA > /dev/null 2>&1

cat $ERRORES
