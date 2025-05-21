import serial
import matplotlib.pyplot as plt

try:
    ser = serial.Serial('COM6', 115200, timeout=1)
except Exception as e:
    print("No se pudo abrir el puerto serie:", e)
    raise SystemExit(1)

archivo = open("datos_STM32.txt", "w")

x_values = []
y_values = []
# z_values = []  # Si deseas también z

try:
    print("Iniciando la recepción de datos...")
    while True:
        line = ser.readline().decode('utf-8').strip()
        if line:
            print("Recibido:", line)
            archivo.write(line + "\n")  # Guardar la línea cruda en el archivo

            # Buscamos el formato "FracGL -> x=..., y=..., z=..."
            if line.startswith("FracGL ->"):
                contenido = line.replace("FracGL ->", "").strip()
                parts = contenido.split(",")
                if len(parts) >= 2:
                    try:
                        x_str = parts[0].split("=")[1].strip()
                        y_str = parts[1].split("=")[1].strip()
                        x_val = float(x_str)
                        y_val = float(y_str)
                        
                        x_values.append(x_val)
                        y_values.append(y_val)
                        
                        # Cada vez que tengamos 1000 muestras, graficamos xy
                        if len(x_values) >= 5000:
                            plt.figure()
                            plt.plot(x_values, y_values, linestyle='-')
                            plt.title("Gráfica XY")
                            plt.xlabel("x")
                            plt.ylabel("y")
                            plt.grid(True)
                            plt.tight_layout()
                            plt.show()

                            x_values.clear()
                            y_values.clear()
                    except Exception as parse_err:
                        print("Error al parsear la línea:", line, " =>", parse_err)
                else:
                    print("Formato no reconocido (se esperan x=, y=, z=):", line)
            else:
                print("Línea con formato inesperado, se omite.")
            
except KeyboardInterrupt:
    print("Interrupción por teclado (Ctrl+C).")
except Exception as e:
    print("Error inesperado:", e)
finally:
    print("Finalizando recepción de datos...")
    archivo.close()
    ser.close()
